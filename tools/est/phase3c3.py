#!/usr/bin/env python3
"""Phase 3c-3 driver: stationary days and the final rho* (NASDAQ cells).

brule    the specified burn-in rule B = max(20000, 20 x the longest lag at which the kernel or the sign
         ACF is still above 1 % of its lag-1 value), evaluated per cell and generator
truth    G0, GT, GF truth with burn-in BURN (independent 10^6-run days, each the last 10^6 runs of a
         series of BURN + 10^6 runs), main sign model at every NASDAQ cell and the 3b tail-lag model at
         phase3c.SENS_CELLS
study    N_PANELS stationary-day panels per (cell, generator) as in 3c-2 (confirmation-sized, exploration
         day lengths, R2 clip, ceiling estimator, B = 200 bootstrap and the t day jackknife, the wider),
         ceilings 1/160 .. 1/5; COMPARE_PANELS fresh-start panels at COMPARE_CELLS for the comparison
analyse  bias, 90 % coverage with Monte Carlo SE, rho*, binding constraints, comparison, reach

  uv run --group analysis tools/est/phase3c3.py brule|truth|study|analyse [--panels N] [--workers 4]
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
import gen3 as G
import panel3c as PC
import phase3c as P3
import pooled as PO
import truth3 as TR

PHASE3C3 = E.data_dir() / "phase3c3"
OUT = P3.OUT
BURN = 20_000
N_PANELS = 400
COMPARE_PANELS = 200
COMPARE_CELLS = ("Q-T1S2", "Q-T3S1")
SEED = 20261107
EDGES = (1 / 160, 1 / 80, 1 / 40, 1 / 20, 1 / 10, 1 / 5)
CEILS = [("ceil", e) for e in EDGES]
GENS = P3.GENS
BIAS_MAX, COV_MIN, MIN_EVAL = 0.02, 0.85, 0.5
RANGE_NAMES = [f"[{lo},{hi}]" for lo, hi in PO.RANGES2]


def criterion(gen: str) -> tuple[str, ...]:
    """Estimators counted toward rho*: R_{1/2}/R_0 under G0 and GT only (reported, not used, under GF)."""
    return ("zeta_V", "zeta_K", "pf") + (() if gen == "GF" else ("Rratio",))


# ---- the burn-in rule ----------------------------------------------------------------------------


def one_percent_lag(f, l_max: float = 1e12) -> float:
    """Largest lag on a log grid up to l_max at which f(l) >= 0.01 f(1); inf when f never falls below."""
    l = np.unique(np.round(np.logspace(0, np.log10(l_max), 2000)))
    y = f(l)
    above = np.flatnonzero(y >= 0.01 * y[0])
    return np.inf if above[-1] == len(l) - 1 else float(l[above[-1]])


def b_rule(cell: str, gen: str) -> dict:
    cl = G.load_cell(cell)
    sm = G.gf_sign_model(cell) if gen == "GF" else cl.sign
    sign_lag = one_percent_lag(lambda l: G.sign_model_C(l, sm["w1"], sm["H"], sm["w2"], sm["phi"]))
    if gen == "G0":
        beta = G.g0_beta(cl)
        kern_lag = one_percent_lag(lambda l: (1 + l) ** -beta)
    elif gen == "GT":
        p = PC.gt_params(cl)
        betas = G.tent_betas(G.size_bins(cl)[1], p["beta_max"], p["c"], p["l_star"])
        kern_lag = np.inf if betas.min() <= 0 else one_percent_lag(lambda l: (1 + l) ** -betas.min())
    else:
        kern_lag = np.inf  # permanent impact of the informed runs
    lag = max(sign_lag, kern_lag)
    return {
        "cell": cell,
        "gen": gen,
        "sign_1pct_lag": sign_lag,
        "kernel_1pct_lag": kern_lag,
        "B_rule": max(20_000, 20 * lag),
        "B_used": BURN,
    }


# ---- truth ---------------------------------------------------------------------------------------


def truth_jobs() -> list[tuple]:
    out = []
    for c in P3.q_cells():
        cl = G.load_cell(c)
        cap = float(np.exp(cl.mu + P3.Z999 * cl.sigma_l))
        out.append(("G0", c, {"beta": G.g0_beta(cl), "q_cap": cap, "burn": BURN}, f"main_{c}_G0"))
        out.append(("GT", c, PC.gt_params(cl) | {"q_cap": cap, "burn": BURN}, f"main_{c}_GT"))
        out.append(("GF", c, {"q_cap": cap, "burn": BURN}, f"main_{c}_GF"))
    for c in P3.SENS_CELLS:
        cl, _, _, sm = P3.setup(c, "G0", "tail")
        cap = float(np.exp(cl.mu + P3.Z999 * cl.sigma_l))
        out.append(("G0", cl, {"beta": G.g0_beta(G.load_cell(c)), "q_cap": cap, "burn": BURN}, f"tail_{c}_G0"))
        out.append(("GT", cl, PC.gt_params(cl) | {"q_cap": cap, "burn": BURN}, f"tail_{c}_GT"))
        out.append(("GF", cl, {"sign_model": sm, "q_cap": cap, "burn": BURN}, f"tail_{c}_GF"))
    return out


def run_truths(workers: int, log=print) -> None:
    tdir = PHASE3C3 / "truth"
    jobs = truth_jobs()
    log(f"truth: {len(jobs)} runs of up to 64 days of 10^6 runs (+ {BURN} burn-in each), about 14 min", flush=True)
    for gen, cell, params, name in jobs:
        path = tdir / f"{name}.npz"
        if path.exists() and "chi" in np.load(path).files:
            continue
        r = TR.run_truth(
            gen,
            cell,
            params,
            name=name,
            min_days=16,
            chunk=16,
            max_days=64,
            workers=workers,
            out_dir=tdir,
            n_collapse=0,
            log=lambda *_: None,
        )
        log(f"{name}: {r.n_days} days, domain SE {r.se_max:.4f}, {r.seconds:.0f}s", flush=True)


def truth_path(tag: str, cell: str, gen: str, fresh: bool = False):
    return P3.truth_path(tag, cell, gen) if fresh else PHASE3C3 / "truth" / f"{tag}_{cell}_{gen}.npz"


def truth_values(path) -> dict:
    """Estimands on RANGES2 and A_GRID2 from a truth run's day means: key -> (value, day-jackknife SE)."""
    z = np.load(path)
    sel = np.isin(TR.A_GRID, PO.A_GRID2)
    T = TR.T_GRID
    lt2 = np.log2(T.astype(float))
    mids = [int(np.flatnonzero(lt2 == np.floor((np.log2(lo) + np.log2(hi)) / 2))[0]) for lo, hi in PO.RANGES2]
    w_pf = PO._ols_weights(PO.A_GRID2[PO.A_GRID2 <= PO.A_MAX2])
    h = PO.A_GRID2 <= PO.A_RATIO
    w_h = PO._ols_weights(PO.A_GRID2[h])
    ia0, ih = int(np.flatnonzero(PO.A_GRID2 == 0)[0]), int(np.flatnonzero(PO.A_GRID2 == 0.5)[0])

    def fn(c):
        zV = np.stack([TR.ols_zeta(T, c["v"][sel], lo, hi) for lo, hi in PO.RANGES2], axis=1)
        zK = np.stack([TR.ols_zeta(T, c["k"][sel], lo, hi) for lo, hi in PO.RANGES2], axis=1)
        R = c["R"][sel]
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
        "_meta": {"n_days": len(z["Vd"]), "se_max": meta.get("se_max")}
    }


# ---- study ---------------------------------------------------------------------------------------


def chunk_path(kind, tag, cell, gen, start):
    return PHASE3C3 / "study" / f"{kind}_{tag}_{cell}_{gen}_{start:03d}.npz"


def study_job(args) -> str:
    kind, tag, cell, gen, start, count = args
    path = chunk_path(*args[:5])
    if path.exists():
        return str(path)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    cl, params, model, sm = P3.setup(cell, gen, tag)
    nc = P3.n_confirm(cell)
    rec: dict[str, list] = {}
    t0 = time.time()
    for i in range(start, start + count):
        rng = np.random.default_rng([SEED, P3.q_cells().index(cell), GENS.index(gen), 0 if tag == "main" else 1, i])
        Nsd = PC.day_lengths(cell, nc, 8, rng)
        p = PC.make_panel(
            gen, cl, Nsd, rng, sign_model=model, sm=sm, params=dict(params), burn=BURN if kind == "stat" else None
        )
        pm = PO.panel_moments(
            p.eps, p.q_mom, p.lm, np.ones(len(p.eps), bool), p.starts, p.lengths, T_grid=PO.T_GRID2, a_grid=PO.A_GRID2
        )
        iv = PO.intervals2(PO.raw_panel(pm, 8), CEILS, rng, **PO.LEGACY_RULE)
        for key in PO.ESTIMANDS2:
            e, lo, hi, _ = iv[key]
            for nm, arr in (("est", e), ("lo", lo), ("hi", hi)):
                rec.setdefault(f"{key}_{nm}", []).append(arr.astype(np.float32))
        rec.setdefault("dropped_K", []).append(iv["dropped_K"])
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, seconds=time.time() - t0, **{k: np.array(v) for k, v in rec.items()})
    return str(path)


def jobs(n_panels: int, n_compare: int) -> list[tuple]:
    out = []
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            for g in GENS:
                out += [("stat", tag, c, g, s, min(10, n_panels - s)) for s in range(0, n_panels, 10)]
    for c in COMPARE_CELLS:
        for g in GENS:
            out += [("fresh", "main", c, g, s, min(10, n_compare - s)) for s in range(0, n_compare, 10)]
    return sorted(out, key=lambda j: -P3.n_confirm(j[2]))


def run_study(n_panels: int, n_compare: int, workers: int, log=print) -> None:
    js = jobs(n_panels, n_compare)
    todo = [j for j in js if not chunk_path(*j[:5]).exists()]
    log(f"study: {len(js)} chunks, {len(todo)} to run on {workers} workers", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for k, _ in enumerate(ex.map(study_job, todo), 1):
            if k % 50 == 0 or k == len(todo):
                log(f"  {k}/{len(todo)}, {(time.time() - t0) / 60:.1f} min", flush=True)


def load(kind, tag, cell, gen, n) -> dict:
    parts = [np.load(chunk_path(kind, tag, cell, gen, s)) for s in range(0, n, 10)]
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0].files if k != "seconds"}
    out["seconds"] = sum(float(p["seconds"]) for p in parts)
    return out


# ---- analysis ------------------------------------------------------------------------------------


def element_names(key) -> list[str]:
    if key in ("zeta_V", "zeta_K"):
        return [f"{key}({a:g}){r}" for a in PO.A_GRID2 for r in RANGE_NAMES]
    return [f"{key}{r}" for r in RANGE_NAMES]


def evaluate(st, tr, m, keys) -> dict:
    """{estimator: summary} under ceiling mask m; coverage counts an undefined interval as not covering."""
    out = {}
    for key in keys:
        e = st[f"{key}_est"][..., m].astype(float)
        lo, hi = st[f"{key}_lo"][..., m].astype(float), st[f"{key}_hi"][..., m].astype(float)
        t = tr[key][0]
        fin = np.isfinite(e)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            b = np.nanmean(np.where(fin, e - t, np.nan), axis=0).ravel()
        cv = np.mean((lo <= t) & (t <= hi), axis=0).ravel()
        se = np.sqrt(cv * (1 - cv) / len(e))
        ev = fin.mean(0).ravel() >= MIN_EVAL
        ok = (np.abs(b) < BIAS_MAX) & (cv >= COV_MIN)
        names = element_names(key)
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


ALL_KEYS = ("zeta_V", "zeta_K", "pf", "Rratio")


def analyse(n_panels: int, n_compare: int, log=print) -> None:
    rep = PHASE3C3 / "report"
    rep.mkdir(parents=True, exist_ok=True)
    ev, rows, timing, dropped, truths = {}, [], [], [], {}
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            for g in GENS:
                tr = truth_values(truth_path(tag, c, g))
                truths[(tag, c, g)] = tr
                st = load("stat", tag, c, g, n_panels)
                timing.append(("stationary", tag, c, g, len(st["zeta_K_est"]), st["seconds"]))
                dk = st["dropped_K"]
                dropped.append((tag, c, g, int(dk.any(axis=(1, 2, 3)).sum()), len(dk)))
                for m, e in enumerate(EDGES):
                    ev[(tag, c, g, e)] = evaluate(st, tr, m, ALL_KEYS)
                    for key, d in ev[(tag, c, g, e)].items():
                        rows.append(
                            {
                                "tag": tag,
                                "cell": c,
                                "gen": g,
                                "estimator": key,
                                "ceiling": e,
                                "in_criterion": int(key in criterion(g)),
                                **{k: v for k, v in d.items() if k != "fails"},
                                "n_fail": len(d["fails"]),
                            }
                        )
    for c in COMPARE_CELLS:
        for g in GENS:
            st = load("fresh", "main", c, g, n_compare)
            timing.append(("fresh", "main", c, g, len(st["zeta_K_est"]), st["seconds"]))
            trf = truth_values(truth_path("main", c, g, fresh=True))
            for m, e in enumerate(EDGES):
                ev[("fresh_vs_fresh", c, g, e)] = evaluate(st, trf, m, ALL_KEYS)
                ev[("fresh_vs_stat", c, g, e)] = evaluate(st, truths[("main", c, g)], m, ALL_KEYS)
    P3.write_tsv(
        OUT / "p3c3_ceilings.tsv",
        rows,
        f"Phase 3c-3 bias and 90 % coverage by ceiling rho, {n_panels} "
        "stationary-day panels; in_criterion marks the estimators counted toward rho*",
    )
    # rho*
    cell_rho, rho_rows = {}, []
    L = [
        "| sign model | cell | " + " | ".join(f"{g} {k}" for g in GENS for k in ALL_KEYS) + " | rho* (cell) |",
        "|---|---|" + "---|" * (len(GENS) * len(ALL_KEYS) + 1),
    ]
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            vals = []
            for g in GENS:
                for key in ALL_KEYS:
                    passing = [e for e in EDGES if ev[(tag, c, g, e)][key]["pass"]]
                    r = max(passing) if passing else np.nan
                    rho_rows.append(
                        {
                            "tag": tag,
                            "cell": c,
                            "gen": g,
                            "estimator": key,
                            "rho_star": r,
                            "in_criterion": int(key in criterion(g)),
                        }
                    )
                    vals.append(edge_str(r) + ("" if key in criterion(g) else " (not counted)"))
            ok = [e for e in EDGES if all(ev[(tag, c, g, e)][k]["pass"] for g in GENS for k in criterion(g))]
            cell_rho[(tag, c)] = max(ok) if ok else np.nan
            L.append(f"| {tag} | {c} | " + " | ".join(vals) + f" | {edge_str(cell_rho[(tag, c)])} |")
    (rep / "rho.md").write_text("\n".join(L) + "\n")
    P3.write_tsv(
        OUT / "p3c3_rho.tsv",
        rho_rows
        + [
            {"tag": t, "cell": c, "gen": "all", "estimator": "criterion", "rho_star": v, "in_criterion": 1}
            for (t, c), v in cell_rho.items()
        ],
        "Phase 3c-3 rho* (largest ceiling passing |bias| < 0.02 and coverage >= 0.85)",
    )
    # binding constraints
    L = [
        "| sign model | cell | ceiling | passes | not evaluable | failing elements | worst bias | "
        "lowest coverage (MC SE) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for tag, cells in (("main", P3.q_cells()), ("tail", P3.SENS_CELLS)):
        for c in cells:
            for e in EDGES:
                fails, ne, n_ev = [], [], 0
                for g in GENS:
                    for key in criterion(g):
                        d = ev[(tag, c, g, e)][key]
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
    # fresh-start vs stationary
    L = [
        "| cell | gen | ceiling | estimator | stationary panels vs stationary truth: max abs bias / min cov | "
        "fresh panels vs fresh truth | fresh panels vs stationary truth |",
        "|---|---|---|---|---|---|---|",
    ]
    for c in COMPARE_CELLS:
        for g in GENS:
            for e in EDGES:
                for key in ALL_KEYS:
                    a = ev[("main", c, g, e)][key]
                    b = ev[("fresh_vs_fresh", c, g, e)][key]
                    d = ev[("fresh_vs_stat", c, g, e)][key]
                    L.append(
                        f"| {c} | {g} | {edge_str(e)} | {key} | {f(a['max_abs_bias'])} / {f(a['min_cov'], 2)} | "
                        f"{f(b['max_abs_bias'])} / {f(b['min_cov'], 2)} | "
                        f"{f(d['max_abs_bias'])} / {f(d['min_cov'], 2)} |"
                    )
    (rep / "compare.md").write_text("\n".join(L) + "\n")
    # truth table, stationary next to fresh
    L = [
        "| sign model | cell | gen | days | stationary zeta_K(a) [8,64], a = 0..1 | fresh (3c) zeta_K(a) [8,64] | "
        "stationary zeta_V(0) [8,64] (SE) | P-F slope [8,64] | R1/2/R0 [8,64] |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for (tag, c, g), tr in truths.items():
        fr = truth_values(truth_path(tag, c, g, fresh=True))
        L.append(
            f"| {tag} | {c} | {g} | {tr['_meta']['n_days']} | {' '.join(f(x) for x in tr['zeta_K'][0][:, 1])} | "
            f"{' '.join(f(x) for x in fr['zeta_K'][0][:, 1])} | "
            f"{f(tr['zeta_V'][0][0, 1])} ({f(tr['zeta_V'][1][0, 1])}) | "
            f"{f(tr['pf'][0][1])} | {f(tr['Rratio'][0][1])} |"
        )
    (rep / "truth.md").write_text("\n".join(L) + "\n")
    L = ["| sign model | cell | gen | panels with a dropped T (any ceiling) | panels |", "|---|---|---|---|---|"]
    L += [f"| {t} | {c} | {g} | {a} | {n} |" for t, c, g, a, n in dropped]
    (rep / "dropped.md").write_text("\n".join(L) + "\n")
    reach(cell_rho, rep)
    L = ["| panels | sign model | cell | gen | panels | worker seconds |", "|---|---|---|---|---|---|"]
    L += [f"| {k} | {t} | {c} | {g} | {n} | {s:.0f} |" for k, t, c, g, n, s in timing]
    (rep / "time.md").write_text("\n".join(L) + "\n")
    log("rho*: " + json.dumps({f"{t} {c}": edge_str(v) for (t, c), v in cell_rho.items()}))


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
    found = [v for (t, c), v in cell_rho.items() if t == "main" and np.isfinite(v)]
    rhos = sorted(set(found)) if found else [1 / 40, 1 / 20]
    rows = []
    L = [
        "| cell | rho | confirm stocks | " + " | ".join(f"{r} contributing" for r in RANGE_NAMES) + " |",
        "|---|---|---|---|---|---|",
    ]
    for rho in rhos:
        for c in P3.q_cells():
            nm = np.array(nmin[c])
            n = [int(np.sum(nm >= hi / rho)) for _, hi in PO.RANGES2]
            rows.append({"group": c, "rho": rho, "confirm": len(nm), **dict(zip(RANGE_NAMES, n, strict=True))})
            L.append(f"| {c} | {edge_str(rho)} | {len(nm)} | " + " | ".join(map(str, n)) + " |")
        for tb in ("T1", "T2", "T3"):  # pooled across the sigma_l cells of a tick bin
            nm = np.concatenate([nmin[c] for c in P3.q_cells() if c[2:4] == tb])
            n = [int(np.sum(nm >= hi / rho)) for _, hi in PO.RANGES2]
            rows.append(
                {
                    "group": f"Q-{tb} (S1+S2+S3)",
                    "rho": rho,
                    "confirm": len(nm),
                    **dict(zip(RANGE_NAMES, n, strict=True)),
                }
            )
            L.append(f"| Q-{tb} (S1+S2+S3) | {edge_str(rho)} | {len(nm)} | " + " | ".join(map(str, n)) + " |")
    (rep / "reach.md").write_text("\n".join(L) + "\n")
    P3.write_tsv(OUT / "p3c3_reach.tsv", rows, "Phase 3c-3 reach: confirmation stocks with N_min >= T2 / rho")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["brule", "truth", "study", "analyse"])
    ap.add_argument("--panels", type=int, default=N_PANELS)
    ap.add_argument("--compare", type=int, default=COMPARE_PANELS)
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    a = ap.parse_args()
    t0 = time.time()
    if a.step == "brule":
        rows = [b_rule(c, g) for c in P3.q_cells() for g in GENS]
        P3.write_tsv(
            OUT / "p3c3_burnin.tsv", rows, "Phase 3c-3 burn-in rule per cell and generator (inf: never below 1 %)"
        )
        for r in rows:
            print(r)
    elif a.step == "truth":
        run_truths(a.workers)
    elif a.step == "study":
        run_study(a.panels, a.compare, a.workers)
    else:
        analyse(a.panels, a.compare)
    print(f"{a.step}: {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
