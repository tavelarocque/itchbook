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
analyse  bias, coverage with MC SE, pass rule, contiguous rho* per cell and per tick bin, binding
         estimators, interval method shares, reach. Applies the rules decided 31 Oct 2026 from
         docs/phase3c4_diagnosis.md to the saved panels: range reach (pooled.range_rule, at least 5 stocks at
         T2, full range, evaluable if defined in 90 % of panels), the P-F slope by coverage only, and tick-bin
         rho*. The saved chunks were computed with pooled.LEGACY_RULE; study keeps it so chunks resume exactly.

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
BIAS_MAX, COV_MIN = 0.02, 0.85
MIN_EVAL = 0.9  # rule 1 (31 Oct 2026, docs/phase3c4_diagnosis.md): defined in >= 90 % of panels; was 0.5
KW = dict(T_grid=T_GRID4, ranges=RANGES4)


def criterion(gen: str) -> tuple[str, ...]:
    """Estimators counted toward rho*: GF is used only for the P-F slope; R_{1/2}/R_0 under G0 and GT."""
    return ("pf",) if gen == "GF" else ("zeta_V", "zeta_K", "pf", "Rratio")


def passes(bias, cov, mc_se, coverage_only: bool = False) -> np.ndarray:
    """|bias| < 0.02, and coverage fails only when coverage + 2 MC SE < 0.85. coverage_only (rule 2, 31 Oct 2026,
    docs/phase3c4_diagnosis.md: the P-F slope) drops the bias half."""
    cov_ok = ~(np.asarray(cov) + 2 * np.asarray(mc_se) < COV_MIN)
    return cov_ok if coverage_only else (np.abs(bias) < BIAS_MAX) & cov_ok


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
        iv = PO.intervals3(PO.raw_panel(pm, 8), CEILS, rng, keys=keys, **KW, **PO.LEGACY_RULE)  # as saved
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
# Rule changes decided 31 Oct 2026 from docs/phase3c4_diagnosis.md, applied to the saved panels (no new
# simulation; bias bound, coverage bound, Monte Carlo band and panels unchanged):
#   1. range reach (pooled.MIN_REACH, pooled.range_rule): an element is defined for a panel only if at least 5
#      stocks reach T2 under the ceiling, and zeta_K only over the full range; it is evaluable at a ceiling only
#      if defined in at least MIN_EVAL = 90 % of panels. Bias and coverage are over the defined panels.
#   2. the P-F slope is judged by coverage only (G0, GT and GF); its bias is still reported.
#   3. rho* per cell and per tick bin (the minimum of its cells' rho*), with the binding estimator named.
COVERAGE_ONLY = ("pf",)


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


def panel_ns(cell, gen, n) -> np.ndarray:
    """(n, S) each stock's fewest runs on any study day, per panel: the day lengths are re-drawn from the panel's
    seeded generator (day_lengths is study_job's first draw, before any simulation)."""
    nc = P3.n_confirm(cell)
    out = []
    for i in range(n):
        rng = np.random.default_rng([SEED, P3.q_cells().index(cell), GENS.index(gen), i])
        out.append(PC.day_lengths(cell, nc, 8, rng).min(axis=1))
    return np.array(out)


def defined_masks(st, cell, gen) -> tuple[dict, np.ndarray]:
    """Rule 1 on saved panels: key -> (P, ..., n_ceilings) defined, and the reach (P, n_ranges, n_ceilings).
    The saved point estimates' dropped_K (non-positive pooled k~ at a horizon) marks a shortened zeta_K range."""
    Ns = panel_ns(cell, gen, len(st["dropped_K"]))
    reach = np.stack([PO.range_reach(ns, CEILS, RANGES4) for ns in Ns])
    return PO.range_rule(reach, kept_full_from_dropped(st["dropped_K"]), PO.A_GRID2, PO.A_MAX2), reach


def kept_full_from_dropped(dk, T_grid=T_GRID4, ranges=RANGES4) -> np.ndarray:
    """dropped_K (..., na, nT, nm) -> (..., na, nr, nm): no horizon of the range dropped. With at least one stock
    at T2 every horizon has eligible stocks, so this is core2's 'every horizon kept'."""
    return np.stack([~dk[..., (T_grid >= lo) & (T_grid <= hi), :].any(axis=-2) for lo, hi in ranges], axis=-2)


def evaluate(st, tr, m, keys, dm) -> dict:
    out = {}
    for key in keys:
        d = dm[key][..., m]
        e = np.where(d, st[f"{key}_est"][..., m].astype(float), np.nan)
        lo, hi = st[f"{key}_lo"][..., m].astype(float), st[f"{key}_hi"][..., m].astype(float)
        t = tr[key][0]
        fin = np.isfinite(e)
        nd = fin.sum(0).ravel()
        with warnings.catch_warnings(), np.errstate(invalid="ignore", divide="ignore"):
            warnings.simplefilter("ignore", RuntimeWarning)
            b = np.nanmean(e - t, axis=0).ravel()
            cv = (np.sum(fin & (lo <= t) & (t <= hi), axis=0).ravel() / nd)  # over defined panels
            se = np.sqrt(cv * (1 - cv) / nd)
        dfrac = fin.mean(0).ravel()
        ev = dfrac >= MIN_EVAL
        ok = passes(b, cv, se, coverage_only=key in COVERAGE_ONLY)
        band = in_mc_band(cv, se) & ev
        names = element_names(key)
        out[key] = {
            "evaluable": bool(ev.any()),
            "n_eval": int(ev.sum()),
            "pass": bool(ev.any() and ok[ev].all()),
            "max_abs_bias": float(np.max(np.abs(b[ev]))) if ev.any() else np.nan,
            "min_cov": float(np.min(cv[ev])) if ev.any() else np.nan,
            "min_defined": float(dfrac.min()),
            "min_defined_eval": float(dfrac[ev].min()) if ev.any() else np.nan,
            "in_band": int(band.sum()),
            "panels": len(e),
            "fails": [(names[i], float(b[i]), float(cv[i]), float(se[i]), float(dfrac[i]))
                      for i in np.flatnonzero(ev & ~ok)],
        }
    return out


def ceiling_status(ev_cg: dict, gens, cell_keys) -> dict:
    """ceiling -> True / False / None for the cell: None when some criterion estimator has no evaluable element."""
    st = {}
    for e in EDGES:
        ds = [ev_cg[(g, e)][k] for g in gens for k in cell_keys(g)]
        st[e] = None if not all(d["evaluable"] for d in ds) else all(d["pass"] for d in ds)
    return st


def tick_bin_rho(cell_rho: dict) -> float:
    """The minimum of the cells' rho*; none (NaN) if any cell has none."""
    v = np.array(list(cell_rho.values()), dtype=float)
    return float(np.min(v)) if np.all(np.isfinite(v)) else np.nan


def fail_score(f, key) -> float:
    """The worst-element score of rho.md: |bias| / 0.02 (not for a coverage-only estimator) + the coverage
    shortfall below 0.85 - 2 MC SE in units of 0.05."""
    b = 0.0 if key in COVERAGE_ONLY else abs(f[1]) / BIAS_MAX
    return b + max(0.0, COV_MIN - f[2] - 2 * f[3]) / 0.05


def topup_combos() -> list[tuple]:
    out = []
    for c in P3.q_cells():
        for g in GENS:
            st = load("ind", c, g)
            if len(st[f"{criterion(g)[0]}_est"]) >= N_TOPUP:
                continue
            tr = truth4(c, g)
            dm, _ = defined_masks(st, c, g)
            evs = [evaluate(st, tr, m, criterion(g), dm) for m in range(len(EDGES))]
            if any(x[k]["in_band"] for x in evs for k in criterion(g)):
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


ABBR = {"zeta_V": "V", "zeta_K": "K", "pf": "PF", "Rratio": "R"}


def analyse(log=print) -> None:
    rep = PHASE3C4 / "report"
    rep.mkdir(parents=True, exist_ok=True)
    rows, ev, timing, meth, reach_med, reach_def = [], {}, [], {}, {}, {}
    for kind, cells, gens in (("ind", P3.q_cells(), GENS), ("cont", STRESS_CELLS, STRESS_GENS)):
        for c in cells:
            for g in gens:
                st = load(kind, c, g)
                tr = truth4(c, g)
                dm, rch = defined_masks(st, c, g)
                reach_med[(kind, c, g)] = np.median(rch, axis=0)
                reach_def[(kind, c, g)] = (rch >= PO.MIN_REACH).mean(axis=0)
                timing.append((kind, c, g, len(st[f"{criterion(g)[0]}_est"]), st["seconds"]))
                for m, e in enumerate(EDGES):
                    ev[(kind, c, g, e)] = evaluate(st, tr, m, criterion(g), dm)
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
                    mm = np.where(dm[key], st[f"{key}_method"], -1).ravel()  # defined elements only
                    meth.setdefault((kind, key), np.zeros(4, int))
                    meth[(kind, key)] += np.bincount(mm[mm >= 0], minlength=4)[:4]
    P3.write_tsv(
        OUT / "p3c4_ceilings.tsv",
        rows,
        "Phase 3c-4 bias, 90 % coverage (widest of three intervals) by ceiling, under the 31 Oct 2026 rules "
        "(range reach >= 5 stocks at T2, evaluable if defined in >= 90 % of panels, P-F slope by coverage only); "
        "ind = independent stationary days, cont = continuous-series stress test",
    )
    combos = [("ind", c, GENS) for c in P3.q_cells()] + [("cont", c, STRESS_GENS) for c in STRESS_CELLS]
    cell_rho, binding = {}, {}
    S = [
        "| days | cell | " + " | ".join(edge_str(e) for e in EDGES) + " | ceilings that pass | rho* |",
        "|---|---|" + "---|" * len(EDGES) + "---|---|",
    ]
    D = [
        "| days | cell | panels (G0 / GT / GF) | " + " | ".join(edge_str(e) for e in EDGES) + " | rho* |",
        "|---|---|---|" + "---|" * len(EDGES) + "---|",
    ]
    for kind, c, gens in combos:
        evc = {(g, e): ev[(kind, c, g, e)] for g in gens for e in EDGES}
        status = ceiling_status(evc, gens, criterion)
        r = rho_star(status)
        cell_rho[(kind, c)] = r
        pan = " / ".join(str(ev[(kind, c, g, EDGES[0])][criterion(g)[0]]["panels"]) for g in gens)
        cs, ds, ok_e = [], [], []
        for e in EDGES:
            fails = [(g, k) + x for g in gens for k in criterion(g) for x in ev[(kind, c, g, e)][k]["fails"]]
            fam = " ".join(f"{ABBR[k]} {sum(x[1] == k for x in fails)}" for k in ABBR if any(x[1] == k for x in fails))
            if status[e] is None:
                cs.append("n/e" + (f" ({fam})" if fails else ""))
                ds.append("n/e")
                continue
            if not fails:
                cs.append("pass")
                ds.append("pass")
                ok_e.append(edge_str(e))
                continue
            wb = max(fails, key=lambda x: fail_score(x[2:], x[1]))
            binding.setdefault((kind, c), {})[e] = wb
            cs.append(f"fail: {fam}")
            ds.append(f"fail {len(fails)}: binding {wb[0]} {wb[1]}, worst {wb[2]} b {wb[3]:+.3f} "
                      f"c {wb[4]:.2f}+-{wb[5]:.2f} def {wb[6]:.2f}")
        S.append(f"| {kind} | {c} | " + " | ".join(cs) + f" | {', '.join(ok_e) or 'none'} | {edge_str(r)} |")
        D.append(f"| {kind} | {c} | {pan} | " + " | ".join(ds) + f" | {edge_str(r)} |")
    (rep / "status.md").write_text("\n".join(S) + "\n")
    (rep / "rho.md").write_text("\n".join(D) + "\n")
    # tick bins: the minimum of the cells' rho*, and the cell and estimator that bind it
    bins, B = {}, ["| tick bin | cells' rho* | rho* | binding |", "|---|---|---|---|"]
    for tb in ("T1", "T2", "T3"):
        cs = [c for c in P3.q_cells() if c[2:4] == tb]
        rb = tick_bin_rho({c: cell_rho[("ind", c)] for c in cs})
        bins[tb] = rb
        why = []
        for c in cs:
            rc = cell_rho[("ind", c)]
            if not (rc == rb or not (np.isfinite(rc) or np.isfinite(rb))):
                continue  # only the cells at the tick bin's rho* bind it
            nxt = [e for e in EDGES if (not np.isfinite(rc) or e > rc) and e in binding.get(("ind", c), {})]
            if nxt:
                wb = binding[("ind", c)][nxt[0]]
                why.append(f"{c} fails at {edge_str(nxt[0])}: {wb[0]} {wb[1]} ({wb[2]})")
        B.append(f"| Q-{tb} | " + ", ".join(f"{c} {edge_str(cell_rho[('ind', c)])}" for c in cs)
                 + f" | {edge_str(rb)} | {'; '.join(why)} |")
    (rep / "rho_bins.md").write_text("\n".join(B) + "\n")
    # per estimator and generator: contiguous rho* of each alone
    L = [
        "| days | cell | " + " | ".join(f"{g} {k}" for g in GENS for k in criterion(g)) + " |",
        "|---|---|" + "---|" * sum(len(criterion(g)) for g in GENS),
    ]
    est_rows = []
    for kind, c, gens in combos:
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
        ]
        + [{"paths": "ind", "cell": f"Q-{tb}", "gen": "all", "estimator": "tick bin", "rho_star": v}
           for tb, v in bins.items()],
        "Phase 3c-4 rho* under the 31 Oct 2026 rules: largest ceiling with every evaluable ceiling up to it "
        "passing; tick bin = the minimum of its cells' rho*",
    )
    L = ["| days | estimator | elements x panels | " + " | ".join(PO.METHODS3) + " |", "|---|---|---|---|---|---|"]
    for (kind, key), cnt in sorted(meth.items()):
        tot = cnt[1:].sum()
        L.append(f"| {kind} | {key} | {tot} | " + " | ".join(f"{cnt[i] / tot:.2f}" for i in (1, 2, 3)) + " |")
    (rep / "methods.md").write_text("\n".join(L) + "\n")
    L = ["| days | cell | gen | panels | worker seconds |", "|---|---|---|---|---|"]
    L += [f"| {k} | {c} | {g} | {n} | {s:.0f} |" for k, c, g, n, s in timing]
    (rep / "time.md").write_text("\n".join(L) + "\n")
    L = ["| days | cell | gen | " + " | ".join(f"{r} {edge_str(e)}" for r in RANGE_NAMES for e in EDGES) + " |",
         "|---|---|---|" + "---|" * (len(RANGE_NAMES) * len(EDGES))]
    for (kind, c, g), md in reach_med.items():
        L.append(f"| {kind} | {c} | {g} | " + " | ".join(f"{md[r, m]:g}" for r in range(len(RANGES4))
                                                         for m in range(len(EDGES))) + " |")
    (rep / "panel_reach.md").write_text("\n".join(L) + "\n")
    L = ["| days | cell | gen | range | " + " | ".join(edge_str(e) for e in EDGES) + " |",
         "|---|---|---|---|" + "---|" * len(EDGES)]
    for (kind, c, g), fr in reach_def.items():
        for r, rn in enumerate(RANGE_NAMES):
            L.append(f"| {kind} | {c} | {g} | {rn} | " + " | ".join(f"{x:.2f}" for x in fr[r]) + " |")
    (rep / "defined.md").write_text("\n".join(L) + "\n")
    reach(bins, rep)
    log("rho*: " + json.dumps({f"{k} {c}": edge_str(v) for (k, c), v in cell_rho.items()}))
    log("tick bins: " + json.dumps({tb: edge_str(v) for tb, v in bins.items()}))


def reach(bins, rep) -> None:
    """Confirmation stocks with N_min >= T2 / rho at their tick bin's rho* (none found: the 1/20 fallback,
    marked)."""
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
    for tb in ("T1", "T2", "T3"):
        r = bins[tb]
        used = r if np.isfinite(r) else 1 / 20
        lab = edge_str(used) + ("" if np.isfinite(r) else " (none found)")
        cs = [c for c in P3.q_cells() if c[2:4] == tb]
        tot, conf = [0, 0], 0
        for c in cs:
            nm = np.array(nmin[c])
            n = [int(np.sum(nm >= hi / used)) for _, hi in RANGES4]
            tot = [a + b for a, b in zip(tot, n, strict=True)]
            conf += len(nm)
            rows.append({"group": c, "rho_used": used, "rho_found": r, "confirm": len(nm),
                         **dict(zip(RANGE_NAMES, n, strict=True))})
            L.append(f"| {c} | {lab} | {len(nm)} | " + " | ".join(map(str, n)) + " |")
        rows.append({"group": f"Q-{tb}", "rho_used": used, "rho_found": r, "confirm": conf,
                     **dict(zip(RANGE_NAMES, tot, strict=True))})
        L.append(f"| **Q-{tb}** (pooled) | {lab} | {conf} | " + " | ".join(map(str, tot)) + " |")
    (rep / "reach.md").write_text("\n".join(L) + "\n")
    P3.write_tsv(OUT / "p3c4_reach.tsv", rows,
                 "Phase 3c-4 reach: confirmation stocks with N_min >= T2 / rho at the tick bin's rho*")


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
