#!/usr/bin/env python3
"""Phase 3c-4 driver: independent stationary days, the widest of three intervals, final rho* (NASDAQ).

check    stationarity check on the new path: G0 zeta_K(0), zeta_V(0) on [8, 64] at Q-T1S3 with day
         lengths N = 300 and 3000, every day its own series after a BURN-run burn-in, against the 10^6-run
         truth (phase3c3 truth, built the same way)
study    N_PANELS panels per NASDAQ cell for G0, GT and GF (GF used only for the P-F slope), independent
         stationary days (panel3c.make_panel(burn=BURN)), pooled.intervals3 (widest of the day-and-stock
         bootstrap, the t day jackknife, the whole-stock bootstrap), ceilings CEILS, ranges RANGES4; and the
         continuous-series stress test (the 3c-3 path) at STRESS_CELLS
topup    panels N_PANELS .. N_TOPUP - 1 for every (cell, generator) with an element inside the Monte Carlo
         band at some ceiling (coverage within 2 MC SE of 0.85)
analyse  bias, coverage with MC SE, pass rule, contiguous rho*, binding constraints, interval method
         shares, reach

  uv run --group analysis tools/est/phase3c4.py check|study|topup|analyse [--panels N] [--workers 3]
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
import phase3c3 as C3
import pooled as PO

PHASE3C4 = E.data_dir() / "phase3c4"
OUT = P3.OUT
BURN = 20_000
N_PANELS, N_TOPUP = 200, 400
SEED = 20261109
EDGES = (1 / 80, 1 / 40, 1 / 20, 1 / 10, 1 / 5)
CEILS = [("ceil", e) for e in EDGES]
T_GRID4 = 2 ** np.arange(7)  # 1 .. 64
RANGES4 = ((4, 32), (8, 64))
RANGE_NAMES = [f"[{lo},{hi}]" for lo, hi in RANGES4]
GENS = ("G0", "GT", "GF")
STRESS_CELLS = ("Q-T3S1",)
STRESS_GENS = ("G0", "GT")
STRESS_PANELS = 100
BIAS_MAX, COV_MIN, MIN_EVAL = 0.02, 0.85, 0.5
KW = dict(T_grid=T_GRID4, ranges=RANGES4)


def criterion(gen: str) -> tuple[str, ...]:
    """Estimators counted toward rho*: GF is used only for the P-F slope; R_{1/2}/R_0 under G0 and GT."""
    return ("pf",) if gen == "GF" else ("zeta_V", "zeta_K", "pf", "Rratio")


def passes(bias, cov, mc_se) -> np.ndarray:
    """|bias| < 0.02, and coverage fails only when coverage + 2 MC SE < 0.85."""
    return (np.abs(bias) < BIAS_MAX) & ~(np.asarray(cov) + 2 * np.asarray(mc_se) < COV_MIN)


def in_mc_band(cov, mc_se) -> np.ndarray:
    return np.abs(np.asarray(cov) - COV_MIN) <= 2 * np.asarray(mc_se)


def rho_star(status: dict) -> float:
    """status: ceiling -> True (passes), False (fails), None (not evaluable). The largest evaluable ceiling such
    that EVERY evaluable ceiling up to and including it passes; NaN when the first evaluable one fails."""
    best = np.nan
    for e in sorted(status):
        if status[e] is None:
            continue
        if not status[e]:
            break
        best = e
    return best


# ---- step 2 check ----------------------------------------------------------------------------------


def run_check(log=print) -> dict:
    cl = G.load_cell("Q-T1S3")
    beta = G.g0_beta(cl)
    T = 2 ** np.arange(7)
    out = {}
    for N, D in ((300, 2000), (3000, 400)):
        stocks = 8
        Nsd = np.full((stocks, D // stocks), N)
        V, K = [], []
        for s in range(stocks):  # one "stock" at a time keeps memory small; stocks are independent replicates
            p = PC.make_panel(
                "G0",
                cl,
                Nsd[s : s + 1],
                np.random.default_rng([SEED, N, s]),
                params={"beta": beta},
                burn=BURN,
                clip=None,
            )
            pm = PO.panel_moments(
                p.eps, p.q, p.lm, np.ones(len(p.eps), bool), p.starts, p.lengths, T_grid=T, a_grid=np.array([0.0])
            )
            V.append(pm.V[:, 0].mean(0))
            K.append(pm.K[:, 0].mean(0))
        V, K = np.array(V), np.array(K)

        def zeta(v, k):
            sel = (T >= 8) & (T <= 64)
            x = np.log(T[sel]) - np.log(T[sel]).mean()
            return (1 + np.log(v[sel] / (T[sel] * v[0])) @ x / (x @ x), 1 + np.log(k[sel] / T[sel]) @ x / (x @ x))

        est = zeta(V.mean(0), K.mean(0))
        n = len(V)
        loo = np.array([zeta(np.delete(V, i, 0).mean(0), np.delete(K, i, 0).mean(0)) for i in range(n)])
        se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean(0)) ** 2, axis=0))
        out[N] = (est, se)
        log(f"N={N}: zeta_V(0) {est[0]:.4f} ({se[0]:.4f})  zeta_K(0) {est[1]:.4f} ({se[1]:.4f})", flush=True)
    tr = C3.truth_values(C3.truth_path("main", "Q-T1S3", "G0"))
    out["truth"] = (tr["zeta_V"][0][0, 1], tr["zeta_V"][1][0, 1], tr["zeta_K"][0][0, 1], tr["zeta_K"][1][0, 1])
    log(
        f"truth (10^6-run days, same burn-in): zeta_V(0) {out['truth'][0]:.4f} ({out['truth'][1]:.4f})  "
        f"zeta_K(0) {out['truth'][2]:.4f} ({out['truth'][3]:.4f})",
        flush=True,
    )
    (PHASE3C4 / "report").mkdir(parents=True, exist_ok=True)
    (PHASE3C4 / "report" / "check.json").write_text(
        json.dumps(
            {
                str(k): [list(map(float, v[0])), list(map(float, v[1]))] if k != "truth" else list(map(float, v))
                for k, v in out.items()
            }
        )
    )
    return out


# ---- study -----------------------------------------------------------------------------------------


def chunk_path(kind, cell, gen, start):
    return PHASE3C4 / "study" / f"{kind}_{cell}_{gen}_{start:03d}.npz"


def study_job(args) -> str:
    kind, cell, gen, start, count = args
    path = chunk_path(*args[:4])
    if path.exists():
        return str(path)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    cl, params, model, sm = P3.setup(cell, gen, "main")
    nc = P3.n_confirm(cell)
    keys = ("pf",) if gen == "GF" else PO.ESTIMANDS2
    rec: dict[str, list] = {}
    t0 = time.time()
    for i in range(start, start + count):
        rng = np.random.default_rng([SEED, P3.q_cells().index(cell), GENS.index(gen), i])
        Nsd = PC.day_lengths(cell, nc, 8, rng)
        p = PC.make_panel(
            gen, cl, Nsd, rng, sign_model=model, sm=sm, params=dict(params), burn=BURN, continuous=kind == "cont"
        )
        pm = PO.panel_moments(
            p.eps, p.q_mom, p.lm, np.ones(len(p.eps), bool), p.starts, p.lengths, T_grid=T_GRID4, a_grid=PO.A_GRID2
        )
        iv = PO.intervals3(PO.raw_panel(pm, 8), CEILS, rng, keys=keys, **KW)
        for key in keys:
            e, lo, hi, meth = iv[key]
            for nm, arr in (("est", e), ("lo", lo), ("hi", hi)):
                rec.setdefault(f"{key}_{nm}", []).append(arr.astype(np.float32))
            rec.setdefault(f"{key}_method", []).append(meth.astype(np.int8))
        rec.setdefault("dropped_K", []).append(iv["dropped_K"])
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, seconds=time.time() - t0, **{k: np.array(v) for k, v in rec.items()})
    return str(path)


def study_jobs(n_from: int, n_to: int, combos) -> list[tuple]:
    js = [(k, c, g, s, min(10, n_to - s)) for k, c, g in combos for s in range(n_from, n_to, 10)]
    return sorted(js, key=lambda j: -P3.n_confirm(j[1]) * (3 if j[2] == "GT" else 1))


def run_jobs(js, workers: int, log=print) -> None:
    todo = [j for j in js if not chunk_path(*j[:4]).exists()]
    log(f"{len(js)} chunks, {len(todo)} to run on {workers} workers", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for k, _ in enumerate(ex.map(study_job, todo), 1):
            if k % 25 == 0 or k == len(todo):
                log(f"  {k}/{len(todo)}, {(time.time() - t0) / 60:.1f} min", flush=True)


def main_combos():
    return [("ind", c, g) for c in P3.q_cells() for g in GENS]


def stress_combos():
    return [("cont", c, g) for c in STRESS_CELLS for g in STRESS_GENS]


def load(kind, cell, gen) -> dict:
    files = sorted((PHASE3C4 / "study").glob(f"{kind}_{cell}_{gen}_*.npz"))
    parts = [np.load(f) for f in files]
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0].files if k != "seconds"}
    out["seconds"] = sum(float(p["seconds"]) for p in parts)
    return out


# ---- analysis -------------------------------------------------------------------------------------


def element_names(key) -> list[str]:
    if key in ("zeta_V", "zeta_K"):
        return [f"{key}({a:g}){r}" for a in PO.A_GRID2 for r in RANGE_NAMES]
    return [f"{key}{r}" for r in RANGE_NAMES]


def truth4(cell, gen) -> dict:
    tr = C3.truth_values(C3.truth_path("main", cell, gen))  # 3c-3 truth: independent days with the same burn-in
    out = {}
    for key, (v, se) in tr.items():
        if key == "_meta":
            continue
        out[key] = (v[..., :2], se[..., :2])  # ranges [4, 32] and [8, 64]
    return out


def evaluate(st, tr, m, keys) -> dict:
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
        ok = passes(b, cv, se)
        band = in_mc_band(cv, se) & ev
        names = element_names(key)
        out[key] = {
            "evaluable": bool(ev.any()),
            "n_eval": int(ev.sum()),
            "pass": bool(ev.any() and ok[ev].all()),
            "max_abs_bias": float(np.max(np.abs(b[ev]))) if ev.any() else np.nan,
            "min_cov": float(np.min(cv[ev])) if ev.any() else np.nan,
            "in_band": int(band.sum()),
            "panels": len(e),
            "fails": [(names[i], float(b[i]), float(cv[i]), float(se[i])) for i in np.flatnonzero(ev & ~ok)],
        }
    return out


def ceiling_status(ev_cg: dict, gens, cell_keys) -> dict:
    """ceiling -> True / False / None for the cell: None when some criterion estimator has no evaluable element."""
    st = {}
    for e in EDGES:
        ds = [ev_cg[(g, e)][k] for g in gens for k in cell_keys(g)]
        st[e] = None if not all(d["evaluable"] for d in ds) else all(d["pass"] for d in ds)
    return st


def topup_combos() -> list[tuple]:
    out = []
    for c in P3.q_cells():
        for g in GENS:
            st = load("ind", c, g)
            if len(st[f"{criterion(g)[0]}_est"]) >= N_TOPUP:
                continue
            tr = truth4(c, g)
            if any(evaluate(st, tr, m, criterion(g))[k]["in_band"] for m in range(len(EDGES)) for k in criterion(g)):
                out.append(("ind", c, g))
    return out


def edge_str(e) -> str:
    return "none" if not np.isfinite(e) else f"1/{round(1 / e)}"


def f(x, d=3) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "-" if not np.isfinite(v) else f"{v:.{d}f}"


def analyse(log=print) -> None:
    rep = PHASE3C4 / "report"
    rep.mkdir(parents=True, exist_ok=True)
    rows, ev, timing, meth = [], {}, [], {}
    for kind, cells, gens in (("ind", P3.q_cells(), GENS), ("cont", STRESS_CELLS, STRESS_GENS)):
        for c in cells:
            for g in gens:
                st = load(kind, c, g)
                tr = truth4(c, g)
                timing.append((kind, c, g, len(st[f"{criterion(g)[0]}_est"]), st["seconds"]))
                for m, e in enumerate(EDGES):
                    ev[(kind, c, g, e)] = evaluate(st, tr, m, criterion(g))
                    for key, d in ev[(kind, c, g, e)].items():
                        rows.append(
                            {
                                "paths": kind,
                                "cell": c,
                                "gen": g,
                                "estimator": key,
                                "ceiling": e,
                                **{k: v for k, v in d.items() if k != "fails"},
                                "n_fail": len(d["fails"]),
                            }
                        )
                for key in criterion(g):
                    mm = st[f"{key}_method"].ravel()
                    meth.setdefault((kind, key), np.zeros(4, int))
                    meth[(kind, key)] += np.bincount(mm[mm >= 0], minlength=4)[:4]
    P3.write_tsv(
        OUT / "p3c4_ceilings.tsv",
        rows,
        "Phase 3c-4 bias, 90 % coverage (widest of three intervals) by "
        "ceiling; ind = independent stationary days, cont = continuous-series stress test",
    )
    cell_rho, L = (
        {},
        [
            "| days | cell | panels (G0 / GT / GF) | " + " | ".join(edge_str(e) for e in EDGES) + " | rho* |",
            "|---|---|---|" + "---|" * len(EDGES) + "---|",
        ],
    )
    for kind, cells, gens in (("ind", P3.q_cells(), GENS), ("cont", STRESS_CELLS, STRESS_GENS)):
        for c in cells:
            evc = {(g, e): ev[(kind, c, g, e)] for g in gens for e in EDGES}
            status = ceiling_status(evc, gens, criterion)
            r = rho_star(status)
            cell_rho[(kind, c)] = r
            pan = " / ".join(str(ev[(kind, c, g, EDGES[0])][criterion(g)[0]]["panels"]) for g in gens)
            cells_ = []
            for e in EDGES:
                if status[e] is None:
                    cells_.append("n/e")
                    continue
                fails = [(g,) + x for g in gens for k in criterion(g) for x in ev[(kind, c, g, e)][k]["fails"]]
                if not fails:
                    cells_.append("pass")
                else:
                    wb = max(fails, key=lambda x: abs(x[2]) / BIAS_MAX + max(0.0, COV_MIN - x[3] - 2 * x[4]) / 0.05)
                    cells_.append(f"fail {len(fails)}: {wb[0]} {wb[1]} b {wb[2]:+.3f} c {wb[3]:.2f}+-{wb[4]:.2f}")
            L.append(f"| {kind} | {c} | {pan} | " + " | ".join(cells_) + f" | {edge_str(r)} |")
    (rep / "rho.md").write_text("\n".join(L) + "\n")
    # per estimator and generator: contiguous rho* of each alone
    L = [
        "| days | cell | " + " | ".join(f"{g} {k}" for g in GENS for k in criterion(g)) + " |",
        "|---|---|" + "---|" * sum(len(criterion(g)) for g in GENS),
    ]
    est_rows = []
    for kind, cells, gens in (("ind", P3.q_cells(), GENS), ("cont", STRESS_CELLS, STRESS_GENS)):
        for c in cells:
            vals = []
            for g in gens:
                for k in criterion(g):
                    stt = {
                        e: (ev[(kind, c, g, e)][k]["pass"] if ev[(kind, c, g, e)][k]["evaluable"] else None)
                        for e in EDGES
                    }
                    r = rho_star(stt)
                    vals.append(edge_str(r))
                    est_rows.append({"paths": kind, "cell": c, "gen": g, "estimator": k, "rho_star": r})
            L.append(f"| {kind} | {c} | " + " | ".join(vals) + " |")
    (rep / "rho_est.md").write_text("\n".join(L) + "\n")
    P3.write_tsv(
        OUT / "p3c4_rho.tsv",
        est_rows
        + [
            {"paths": k, "cell": c, "gen": "all", "estimator": "criterion", "rho_star": v}
            for (k, c), v in cell_rho.items()
        ],
        "Phase 3c-4 rho*: largest ceiling with every evaluable ceiling up to it passing",
    )
    L = ["| days | estimator | elements x panels | " + " | ".join(PO.METHODS3) + " |", "|---|---|---|---|---|---|"]
    for (kind, key), cnt in sorted(meth.items()):
        tot = cnt[1:].sum()
        L.append(f"| {kind} | {key} | {tot} | " + " | ".join(f"{cnt[i] / tot:.2f}" for i in (1, 2, 3)) + " |")
    (rep / "methods.md").write_text("\n".join(L) + "\n")
    L = ["| days | cell | gen | panels | worker seconds |", "|---|---|---|---|---|"]
    L += [f"| {k} | {c} | {g} | {n} | {s:.0f} |" for k, c, g, n, s in timing]
    (rep / "time.md").write_text("\n".join(L) + "\n")
    reach(cell_rho, rep)
    log("rho*: " + json.dumps({f"{k} {c}": edge_str(v) for (k, c), v in cell_rho.items()}))


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
    L = ["| group | rho used | confirm stocks | " + " | ".join(RANGE_NAMES) + " |", "|---|---|---|---|---|"]
    for c in P3.q_cells():
        r = cell_rho[("ind", c)]
        used = r if np.isfinite(r) else 1 / 20
        nm = np.array(nmin[c])
        n = [int(np.sum(nm >= hi / used)) for _, hi in RANGES4]
        rows.append(
            {"group": c, "rho_used": used, "rho_found": r, "confirm": len(nm), **dict(zip(RANGE_NAMES, n, strict=True))}
        )
        L.append(
            f"| {c} | {edge_str(used)}{'' if np.isfinite(r) else ' (none found)'} | {len(nm)} | "
            + " | ".join(map(str, n))
            + " |"
        )
    for tb in ("T1", "T2", "T3"):
        cs = [c for c in P3.q_cells() if c[2:4] == tb]
        tot = [0, 0]
        for c in cs:
            r = cell_rho[("ind", c)]
            used = r if np.isfinite(r) else 1 / 20
            nm = np.array(nmin[c])
            for i, (_, hi) in enumerate(RANGES4):
                tot[i] += int(np.sum(nm >= hi / used))
        conf = sum(len(nmin[c]) for c in cs)
        rows.append(
            {
                "group": f"Q-{tb}",
                "rho_used": "per cell",
                "rho_found": "",
                "confirm": conf,
                **dict(zip(RANGE_NAMES, tot, strict=True)),
            }
        )
        L.append(
            f"| Q-{tb} (S1+S2+S3, each at its cell's rho) | per cell | {conf} | " + " | ".join(map(str, tot)) + " |"
        )
    (rep / "reach.md").write_text("\n".join(L) + "\n")
    P3.write_tsv(OUT / "p3c4_reach.tsv", rows, "Phase 3c-4 reach: confirmation stocks with N_min >= T2 / rho")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["check", "study", "topup", "analyse"])
    ap.add_argument("--panels", type=int, default=N_PANELS)
    ap.add_argument("--workers", type=int, default=min(3, os.cpu_count() or 1))
    a = ap.parse_args()
    t0 = time.time()
    if a.step == "check":
        run_check()
    elif a.step == "study":
        run_jobs(study_jobs(0, a.panels, main_combos()) + study_jobs(0, STRESS_PANELS, stress_combos()), a.workers)
    elif a.step == "topup":
        combos = topup_combos()
        print(f"top-up to {N_TOPUP} panels: {combos}", flush=True)
        run_jobs(study_jobs(N_PANELS, N_TOPUP, combos), a.workers)
    else:
        analyse()
    print(f"{a.step}: {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
