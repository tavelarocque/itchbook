#!/usr/bin/env python3
"""Phase 3c-2 driver: rho* with the k-tilde normalisation (NASDAQ cells).

study    N_PANELS panels per (cell, generator) for G0, GT, GF (and the 3b
         tail-lag sign model at phase3c.SENS_CELLS), sized like the cell's
         confirmation set with exploration day lengths, R2 clip on the moments'
         sizes (phase3c.setup and panel3c.make_panel, as in 3c). Per panel:
         pooled.intervals2 (norm="tilde", B = 200 percentile bootstrap and the
         Student t day jackknife, the wider, 90 %) with a T/N band mask per bin;
         the stock-cluster jackknife diagnostic on the first STOCK_JK_PANELS
         panels; at COMPARE_CELLS also the 3c path (norm="old") on the same
         panels. Chunks of 10 panels under PHASE3C2/study/, resumable.
analyse  bias, 90 % coverage and its Monte Carlo SE by bin against the 3c
         truth (recomputed on the 3c-2 ranges and a grid); rho*, binding
         constraints, old vs new, reach.

  uv run --group analysis tools/est/phase3c2.py study|analyse [--panels N] [--workers 4]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import explore as E
import panel3c as PC
import phase3c as P3
import pooled as PO
import truth3 as TR

PHASE3C2 = E.data_dir() / "phase3c2"
OUT = P3.OUT
N_PANELS = 400
SEED = 20261106
EDGES = (0.0, 1 / 160, 1 / 80, 1 / 40, 1 / 20, 1 / 10)
# rho* is evaluated on the ceiling estimator (stock s at T iff T <= rho N_s, the rule Phase 4 applies), one
# mask per bin edge; the 3c-style bands (each increment restricted to lo < T / N_s <= hi) are a diagnostic
# on the first BAND_PANELS panels (CHOICE, see the report: in a band the matched set changes from one
# increment to the next, so the small-sample bias of each log of a matched mean does not telescope)
CEILS = [("ceil", hi) for hi in EDGES[1:]]
BANDS = [("band", lo, hi) for lo, hi in zip(EDGES[:-1], EDGES[1:], strict=True)]
BINS = {i: EDGES[i + 1] for i in range(5)}  # mask index -> bin edge (same order for CEILS and BANDS)
BAND_PANELS = 100
CRITERION = ("zeta_V", "zeta_K", "pf", "Rratio")
STOCK_JK_PANELS = 50
COMPARE_CELLS = ("Q-T1S2", "Q-T3S1")
BIAS_MAX, COV_MIN, MIN_EVAL = 0.02, 0.85, 0.5
GENS = P3.GENS
RANGE_NAMES = [f"[{lo},{hi}]" for lo, hi in PO.RANGES2]


# ---- truth on the 3c-2 ranges ---------------------------------------------------------


def truth2(tag: str, cell: str, gen: str) -> dict:
    """The 3c truth run's day-mean curves on RANGES2 and A_GRID2: estimand -> (value, jackknife SE)."""
    z = np.load(P3.truth_path(tag, cell, gen))
    a = TR.A_GRID
    sel = np.isin(a, PO.A_GRID2)
    T = TR.T_GRID
    lt2 = np.log2(T.astype(float))
    mids = [int(np.flatnonzero(lt2 == np.floor((np.log2(lo) + np.log2(hi)) / 2))[0]) for lo, hi in PO.RANGES2]
    w_pf = PO._ols_weights(PO.A_GRID2[PO.A_GRID2 <= PO.A_MAX2])
    h = PO.A_GRID2 <= PO.A_RATIO
    w_h = PO._ols_weights(PO.A_GRID2[h])

    def fn(c):
        zV = np.stack([TR.ols_zeta(T, c["v"][sel], lo, hi) for lo, hi in PO.RANGES2], axis=1)  # (na, nr)
        zK = np.stack([TR.ols_zeta(T, c["k"][sel], lo, hi) for lo, hi in PO.RANGES2], axis=1)
        R = c["R"][sel]
        ia0, ih = int(np.flatnonzero(PO.A_GRID2 == 0)[0]), int(np.flatnonzero(PO.A_GRID2 == 0.5)[0])
        return {
            "zeta_V": zV,
            "zeta_K": zK,
            "pf": w_pf @ zK[PO.A_GRID2 <= PO.A_MAX2],
            "vk_diff": w_h @ zV[h] - 2 * (w_h @ zK[h]),
            "Rratio": np.array([R[ih, t] / R[ia0, t] for t in mids]),
        }

    est = fn(TR.curves(z["Vd"].mean(0), z["Kd"].mean(0), z["Sd"].mean(0)))
    se = TR.jackknife(z["Vd"], z["Kd"], z["Sd"], fn=fn)
    meta = json.loads(str(z["meta"]))
    return {k: (np.asarray(v, float), np.asarray(se[k], float)) for k, v in est.items()} | {
        "_meta": {"n_days": len(z["Vd"]), "se_max": meta.get("se_max"), "se_met": meta.get("se_met")}
    }


# ---- study ------------------------------------------------------------------------------


def chunk_path(tag, cell, gen, start):
    return PHASE3C2 / "study" / f"{tag}_{cell}_{gen}_{start:03d}.npz"


def study_job(args) -> str:
    tag, cell, gen, start, count = args
    path = chunk_path(tag, cell, gen, start)
    if path.exists():
        return str(path)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    cl, params, kind, sm = P3.setup(cell, gen, tag)
    nc = P3.n_confirm(cell)
    compare = tag == "main" and cell in COMPARE_CELLS
    rec: dict[str, list] = {}

    def put(k, v):
        rec.setdefault(k, []).append(np.asarray(v, dtype=np.float32))

    t0 = time.time()
    for i in range(start, start + count):
        rng = np.random.default_rng([SEED, P3.q_cells().index(cell), GENS.index(gen), 0 if tag == "main" else 1, i])
        Nsd = PC.day_lengths(cell, nc, 8, rng)
        p = PC.make_panel(gen, cl, Nsd, rng, sign_model=kind, sm=sm, params=dict(params))
        pm = PO.panel_moments(
            p.eps, p.q_mom, p.lm, np.ones(len(p.eps), bool), p.starts, p.lengths, T_grid=PO.T_GRID2, a_grid=PO.A_GRID2
        )
        raw = PO.raw_panel(pm, 8)
        iv = PO.intervals2(raw, CEILS, rng, stock_jk=i < STOCK_JK_PANELS)
        if i < BAND_PANELS:
            bv = PO.intervals2(raw, BANDS, np.random.default_rng([SEED, 98, i]))
            for key in PO.ESTIMANDS2:
                e, lo, hi, _ = bv[key]
                put(f"band_{key}_est", e)
                put(f"band_{key}_lo", lo)
                put(f"band_{key}_hi", hi)
        for key in PO.ESTIMANDS2:
            e, lo, hi, _ = iv[key]
            put(f"{key}_est", e)
            put(f"{key}_lo", lo)
            put(f"{key}_hi", hi)
            sj = iv.get("stock_jk", {}).get(key)
            put(f"{key}_slo", sj[0] if sj is not None else np.full(e.shape, np.nan))
            put(f"{key}_shi", sj[1] if sj is not None else np.full(e.shape, np.nan))
        put("dropped_K", iv["dropped_K"])
        put("Ns", np.percentile(Nsd.min(1), [10, 50, 90]))
        if compare:
            ov = PO.intervals2(raw, CEILS, np.random.default_rng([SEED, 99, i]), norm="old")
            for key in PO.ESTIMANDS2:
                e, lo, hi, _ = ov[key]
                put(f"old_{key}_est", e)
                put(f"old_{key}_lo", lo)
                put(f"old_{key}_hi", hi)
            put("old_excl", np.mean(~(raw.K1_t > PO.K1_Z), axis=(0, 1)))
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, seconds=time.time() - t0, **{k: np.array(v) for k, v in rec.items()})
    return str(path)


def jobs(n_panels: int) -> list[tuple]:
    out = []
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            for g in GENS:
                out += [(tag, c, g, s, min(10, n_panels - s)) for s in range(0, n_panels, 10)]
    return sorted(out, key=lambda j: -P3.n_confirm(j[1]))


def run_study(n_panels: int, workers: int, log=print) -> None:
    js = jobs(n_panels)
    todo = [j for j in js if not chunk_path(*j[:4]).exists()]
    log(f"study: {len(js)} chunks of up to 10 panels, {len(todo)} to run on {workers} workers", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for k, _ in enumerate(ex.map(study_job, todo), 1):
            if k % 40 == 0 or k == len(todo):
                el = time.time() - t0
                log(f"  {k}/{len(todo)}, {el / 60:.1f} min", flush=True)


def load(tag, cell, gen, n_panels) -> dict:
    parts = [np.load(chunk_path(tag, cell, gen, s)) for s in range(0, n_panels, 10)]
    keys = {k for p in parts for k in p.files if k != "seconds"}
    out = {k: np.concatenate([p[k] for p in parts if k in p.files]) for k in keys}
    out["seconds"] = sum(float(p["seconds"]) for p in parts)
    return out


# ---- analysis ----------------------------------------------------------------------------


def stats(st, truth, key, m, prefix=""):
    """Per element under mask m: bias, coverage, MC SE of coverage, share of panels with an estimate."""
    e = st[f"{prefix}{key}_est"][..., m].astype(float)
    lo = st[f"{prefix}{key}_lo"][..., m].astype(float)
    hi = st[f"{prefix}{key}_hi"][..., m].astype(float)
    t = truth[key][0]
    fin = np.isfinite(e)
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bias = np.nanmean(np.where(fin, e - t, np.nan), axis=0)
    cov = np.mean((lo <= t) & (t <= hi), axis=0)
    return bias, cov, np.sqrt(cov * (1 - cov) / len(e)), fin.mean(0)


def element_names(key) -> list[str]:
    if key in ("zeta_V", "zeta_K"):
        return [f"{key}({a:g}){r}" for a in PO.A_GRID2 for r in RANGE_NAMES]
    return [f"{key}{r}" for r in RANGE_NAMES]


def evaluate(st, truth, m, prefix=""):
    """{estimator: dict(evaluable, pass, worst bias, worst coverage, failing elements)} at mask m."""
    out = {}
    for key in CRITERION:
        b, cv, se, fin = (x.ravel() for x in stats(st, truth, key, m, prefix))
        names = element_names(key)
        ev = fin >= MIN_EVAL
        ok = (np.abs(b) < BIAS_MAX) & (cv >= COV_MIN)
        fails = [(names[i], float(b[i]), float(cv[i]), float(se[i])) for i in np.flatnonzero(ev & ~ok)]
        out[key] = {
            "evaluable": bool(ev.any()),
            "n_eval": int(ev.sum()),
            "pass": bool(ev.any() and ok[ev].all()),
            "max_abs_bias": float(np.max(np.abs(b[ev]))) if ev.any() else np.nan,
            "min_cov": float(np.min(cv[ev])) if ev.any() else np.nan,
            "min_cov_se": float(se[ev][np.argmin(cv[ev])]) if ev.any() else np.nan,
            "fails": fails,
        }
    return out


def edge_str(e) -> str:
    return "none" if not np.isfinite(e) else f"1/{round(1 / e)}"


def f(x, d=3) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "-" if not np.isfinite(v) else f"{v:.{d}f}"


def analyse(n_panels: int, log=print) -> None:
    rep = PHASE3C2 / "report"
    rep.mkdir(parents=True, exist_ok=True)
    rows, ev_all, timing, truths, dropped = [], {}, [], {}, []
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            for g in GENS:
                tr = truth2(tag, c, g)
                truths[(tag, c, g)] = tr
                st = load(tag, c, g, n_panels)
                timing.append((tag, c, g, len(st["zeta_K_est"]), st["seconds"]))
                dk = st["dropped_K"]  # (P, na, nT, nm)
                dropped.append((tag, c, g, int(dk[..., -1].any(axis=(1, 2)).sum()), int(dk.any(axis=(1, 2, 3)).sum())))
                for m, edge in BINS.items():
                    ev = evaluate(st, tr, m)
                    ev_all[(tag, c, g, edge)] = ev
                    ev_all[("band", tag, c, g, edge)] = evaluate(st, tr, m, "band_")
                    for key, d in ev.items():
                        rows.append(
                            {
                                "tag": tag,
                                "cell": c,
                                "gen": g,
                                "estimator": key,
                                "bin": edge,
                                **{k: v for k, v in d.items() if k != "fails"},
                                "n_fail": len(d["fails"]),
                            }
                        )
                if tag == "main" and c in COMPARE_CELLS:
                    for m, edge in BINS.items():
                        ev_all[("old", c, g, edge)] = evaluate(st, tr, m, "old_")
                        ev_all[("stockjk", c, g, edge)] = stockjk_eval(st, tr, m)
    P3.write_tsv(
        OUT / "p3c2_bins.tsv",
        rows,
        f"Phase 3c-2 bias and 90 % coverage by T/N bin, {n_panels} panels, k-tilde normalisation",
    )
    # rho* per (tag, cell, gen, estimator) and per cell
    rho_rows, cell_rho = [], {}
    L = [
        "| sign model | cell | " + " | ".join(f"{g} {k}" for g in GENS for k in CRITERION) + " | rho* (cell) |",
        "|---|---|" + "---|" * (len(GENS) * len(CRITERION) + 1),
    ]
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            vals = []
            for g in GENS:
                for key in CRITERION:
                    passing = [e for e in BINS.values() if ev_all[(tag, c, g, e)][key]["pass"]]
                    r = max(passing) if passing else np.nan
                    vals.append(r)
                    rho_rows.append({"tag": tag, "cell": c, "gen": g, "estimator": key, "rho_star": r})
            ok_bins = [
                e for e in BINS.values() if all(ev_all[(tag, c, g, e)][k]["pass"] for g in GENS for k in CRITERION)
            ]
            cell_rho[(tag, c)] = max(ok_bins) if ok_bins else np.nan
            L.append(
                f"| {tag} | {c} | " + " | ".join(edge_str(v) for v in vals) + f" | {edge_str(cell_rho[(tag, c)])} |"
            )
    (rep / "rho.md").write_text("\n".join(L) + "\n")
    P3.write_tsv(
        OUT / "p3c2_rho.tsv",
        rho_rows
        + [{"tag": t, "cell": c, "gen": "all", "estimator": "all", "rho_star": v} for (t, c), v in cell_rho.items()],
        "Phase 3c-2 rho*",
    )
    # band diagnostic (first BAND_PANELS panels): does every estimator pass in each band?
    L = [
        "| sign model | cell | "
        + " | ".join(f"band <= {edge_str(e)}: failing / evaluable, worst bias" for e in BINS.values())
        + " |",
        "|---|---|" + "---|" * len(BINS),
    ]
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            cells_ = []
            for e in BINS.values():
                fails, n_ev = [], 0
                for g in GENS:
                    for key in CRITERION:
                        d = ev_all[("band", tag, c, g, e)][key]
                        n_ev += d["n_eval"]
                        fails += [(g,) + x for x in d["fails"]]
                wb = max(fails, key=lambda x: abs(x[2])) if fails else None
                cells_.append(f"{len(fails)} / {n_ev}" + (f", {wb[0]} {wb[1]} {wb[2]:+.3f}" if wb else ""))
            L.append(f"| {tag} | {c} | " + " | ".join(cells_) + " |")
    (rep / "bands.md").write_text("\n".join(L) + "\n")
    # binding constraints per cell and bin
    L = [
        "| sign model | cell | T/N bin | passes | not evaluable | failing elements (of evaluable) | worst bias | "
        "lowest coverage (MC SE) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            for e in BINS.values():
                fails, ne, n_ev = [], [], 0
                for g in GENS:
                    for key in CRITERION:
                        d = ev_all[(tag, c, g, e)][key]
                        if not d["evaluable"]:
                            ne.append(f"{g} {key}")
                        n_ev += d["n_eval"]
                        fails += [(g,) + x for x in d["fails"]]
                wb = max(fails, key=lambda x: abs(x[2])) if fails else None
                wc = min(fails, key=lambda x: x[3]) if fails else None
                L.append(
                    f"| {tag} | {c} | {edge_str(e)} | {'yes' if not fails and not ne else 'no'} | "
                    f"{', '.join(ne) if ne else '-'} | {len(fails)} of {n_ev} | "
                    f"{f'{wb[0]} {wb[1]} {wb[2]:+.3f}' if wb else '-'} | "
                    f"{f'{wc[0]} {wc[1]} {wc[3]:.3f} ({wc[4]:.3f})' if wc else '-'} |"
                )
    (rep / "binding.md").write_text("\n".join(L) + "\n")
    # old vs new, and the stock jackknife, at the compare cells
    L = [
        "| cell | gen | T/N bin | estimator | new: max abs bias / min coverage | old: max abs bias / min coverage | "
        "stock jackknife (new, first 50 panels): min coverage |",
        "|---|---|---|---|---|---|---|",
    ]
    for c in COMPARE_CELLS:
        for g in GENS:
            for e in BINS.values():
                for key in CRITERION:
                    n, o, s = (ev_all[(x, c, g, e)][key] for x in ("main", "old", "stockjk"))
                    L.append(
                        f"| {c} | {g} | {edge_str(e)} | {key} | {f(n['max_abs_bias'])} / {f(n['min_cov'], 2)} | "
                        f"{f(o['max_abs_bias'])} / {f(o['min_cov'], 2)} | {f(s['min_cov'], 2)} |"
                    )
    (rep / "oldnew.md").write_text("\n".join(L) + "\n")
    # truth table
    L = [
        "| sign model | cell | gen | zeta_K(a) [8,64], a = 0..1 | zeta_V(a) [8,64] | "
        "P-F slope [4,32] / [8,64] / [16,256] | "
        "R1/2/R0 [8,64] | vk_diff [8,64] (descriptive) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for (tag, c, g), tr in truths.items():
        L.append(
            f"| {tag} | {c} | {g} | {' '.join(f(x) for x in tr['zeta_K'][0][:, 1])} | "
            f"{' '.join(f(x) for x in tr['zeta_V'][0][:, 1])} | {' / '.join(f(x) for x in tr['pf'][0])} | "
            f"{f(tr['Rratio'][0][1])} | {f(tr['vk_diff'][0][1])} ({f(tr['vk_diff'][1][1])}) |"
        )
    (rep / "truth.md").write_text("\n".join(L) + "\n")
    # dropped horizons
    L = [
        "| sign model | cell | gen | panels with a dropped T (ceiling 1/10) | panels with a dropped T at any ceiling |",
        "|---|---|---|---|---|",
    ]
    L += [f"| {t} | {c} | {g} | {a} | {b} |" for t, c, g, a, b in dropped]
    (rep / "dropped.md").write_text("\n".join(L) + "\n")
    # reach
    reach(cell_rho, rep)
    L = ["| sign model | cell | gen | panels | worker seconds |", "|---|---|---|---|---|"]
    L += [f"| {t} | {c} | {g} | {n} | {s:.0f} |" for t, c, g, n, s in timing]
    (rep / "time.md").write_text("\n".join(L) + "\n")
    log(f"rho*: {json.dumps({f'{k[0]} {k[1]}': edge_str(v) for k, v in cell_rho.items()})}")


def stockjk_eval(st, tr, m):
    """Coverage of the stock-cluster jackknife intervals (first STOCK_JK_PANELS panels), per estimator."""
    out = {}
    for key in CRITERION:
        lo = st[f"{key}_slo"][:STOCK_JK_PANELS, ..., m].astype(float)
        hi = st[f"{key}_shi"][:STOCK_JK_PANELS, ..., m].astype(float)
        e = st[f"{key}_est"][:STOCK_JK_PANELS, ..., m].astype(float)
        t = tr[key][0]
        ev = np.isfinite(e).mean(0) >= MIN_EVAL
        cov = np.mean((lo <= t) & (t <= hi), axis=0)
        out[key] = {"min_cov": float(cov[ev].min()) if ev.any() else np.nan}
    return out


def reach(cell_rho, rep) -> None:
    head, nmin = None, {}
    for line in open(E.UNIVERSE):
        if line.startswith("#"):
            continue
        row = line.rstrip("\n").split("\t")
        if head is None:
            head = {k: i for i, k in enumerate(row)}
            continue
        if row[head["split_v2"]] == "confirm" and row[head["cell"]].startswith("Q-"):
            nmin.setdefault(row[head["cell"]], []).append(int(row[head["n_min"]]))
    rows = []
    L = [
        "| cell | rho used | confirm stocks | "
        + " | ".join(f"{r}: contributing / N_floor 100, 200, 400, 800" for r in RANGE_NAMES)
        + " |",
        "|---|---|---|---|---|---|",
    ]
    for c in P3.q_cells():
        rho = cell_rho[("main", c)]
        used = rho if np.isfinite(rho) else 1 / 40
        nm = np.array(nmin[c])
        cells = []
        for (_, hi), rn in zip(PO.RANGES2, RANGE_NAMES, strict=True):
            contrib = nm >= hi / used
            fl = [int(np.sum(contrib & (nm >= x))) for x in (100, 200, 400, 800)]
            cells.append(f"{int(contrib.sum())} / {', '.join(map(str, fl))}")
            rows.append(
                {
                    "cell": c,
                    "rho_used": used,
                    "rho_found": rho,
                    "range": rn,
                    "contributing": int(contrib.sum()),
                    **{f"floor_{x}": v for x, v in zip((100, 200, 400, 800), fl, strict=True)},
                }
            )
        L.append(
            f"| {c} | {edge_str(used)}{'' if np.isfinite(rho) else ' (none found)'} | {len(nm)} | "
            + " | ".join(cells)
            + " |"
        )
    (rep / "reach.md").write_text("\n".join(L) + "\n")
    P3.write_tsv(
        OUT / "p3c2_reach.tsv",
        rows,
        "Phase 3c-2 reach: confirmation stocks with N_min >= T2 / rho (N_min from universe.tsv)",
    )


# ---- finite-N truth (diagnostic) ------------------------------------------------------------

FINITE_N = (300, 1000, 3000, 10_000)
FINITE_RUNS = 2_000_000


def finite_truth(workers: int, log=print) -> None:
    """G0 and GT population curves at day lengths FINITE_N (about FINITE_RUNS runs each), next to the
    10^6-run truth: the day-start transient (no impact before the first run of a day)."""
    out_dir = PHASE3C2 / "finite"
    rows = []
    for c in P3.q_cells():
        cl = P3.G.load_cell(c)
        cap = float(np.exp(cl.mu + P3.Z999 * cl.sigma_l))
        for g in ("G0", "GT"):
            params = ({"beta": P3.G.g0_beta(cl)} if g == "G0" else PC.gt_params(cl)) | {"q_cap": cap}
            for N in FINITE_N + (10**6,):
                if N == 10**6:
                    z = np.load(P3.truth_path("main", c, g))
                    Vd, Kd, Sd = z["Vd"], z["Kd"], z["Sd"]
                else:
                    days = max(100, FINITE_RUNS // N)
                    r = TR.run_truth(
                        g,
                        c,
                        params,
                        name=f"{g}_{c}_N{N}",
                        N=N,
                        min_days=days,
                        max_days=days,
                        chunk=days,
                        workers=workers,
                        out_dir=out_dir,
                        n_collapse=0,
                        log=lambda *_: None,
                    )
                    Vd, Kd, Sd = r.days
                cu = TR.curves(Vd.mean(0), Kd.mean(0), Sd.mean(0))
                sel = np.isin(TR.A_GRID, PO.A_GRID2)
                row = {"cell": c, "gen": g, "N": N, "days": len(Vd)}
                for lo, hi in PO.RANGES2:
                    zk = TR.ols_zeta(TR.T_GRID, cu["k"][sel], lo, hi)
                    zv = TR.ols_zeta(TR.T_GRID, cu["v"][sel], lo, hi)
                    row[f"zK0_{lo}_{hi}"], row[f"zK1_{lo}_{hi}"] = zk[0], zk[-1]
                    row[f"zV0_{lo}_{hi}"] = zv[0]
                rows.append(row)
                log(f"{c} {g} N={N}: zK(0) [4,32] {row['zK0_4_32']:.3f} [8,64] {row['zK0_8_64']:.3f}", flush=True)
    P3.write_tsv(
        OUT / "p3c2_finite_truth.tsv",
        rows,
        "Phase 3c-2 population curves at finite day length N (G0, GT; R2 cap fixed) next to the 10^6-run truth",
    )
    L = [
        "| cell | gen | "
        + " | ".join(f"N = {n}: zK(0) [4,32] / [8,64] / [16,256]" for n in FINITE_N + (10**6,))
        + " |",
        "|---|---|" + "---|" * (len(FINITE_N) + 1),
    ]
    for c in P3.q_cells():
        for g in ("G0", "GT"):
            cells = [
                " / ".join(f(r[f"zK0_{lo}_{hi}"]) for lo, hi in PO.RANGES2)
                for r in rows
                if r["cell"] == c and r["gen"] == g
            ]
            L.append(f"| {c} | {g} | " + " | ".join(cells) + " |")
    (PHASE3C2 / "report").mkdir(parents=True, exist_ok=True)
    (PHASE3C2 / "report" / "finite.md").write_text("\n".join(L) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["study", "analyse", "finite"])
    ap.add_argument("--panels", type=int, default=N_PANELS)
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    a = ap.parse_args()
    t0 = time.time()
    if a.step == "study":
        run_study(a.panels, a.workers)
    elif a.step == "finite":
        finite_truth(a.workers)
    else:
        analyse(a.panels)
    print(f"{a.step}: {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
