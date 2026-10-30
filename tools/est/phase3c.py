#!/usr/bin/env python3
"""Phase 3c driver: operating characteristics of the pooled estimator (NASDAQ cells).

Steps (each checkpointed under PHASE3C = $ITCH_DATA/phase3c/, resumable):

  truth   infinite-N truth with the R2 clip as a fixed size quantile: G0 and GF
          per cell (cap exp(mu + 3.0902 sigma_l)), GT from the 3b runs with the
          same cap, seed and settings, GM (PAPER2) and GM_floor_G01 once (cap =
          the 99.9 % quantile of trade sizes of a pilot day), and the 3b
          tail-lag sign model at SENS_CELLS (G0, GT, GF).
  study   N_PANELS panels per (cell, generator, sign model), sized like the
          cell's confirmation set, each stock's 8 day lengths copied from a
          random exploration stock; pooled.intervals with every band and
          ceiling mask (B = 200 and the day jackknife, the wider, 90 %).
  analyse bias and 90 % coverage against the truth by T/N band, rho*, the
          coverage-vs-ESS threshold and a_max, the N_floor table.
  tmin    T_min from exploration mids alone; l_max and beta-hat bias on G0.
  semi    the semi-synthetic generator on real exploration flow (smoke check).

    uv run --group analysis tools/est/phase3c.py truth|study|analyse|tmin|semi [--panels N] [--workers 4]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path

import numpy as np

import explore as E
import gen3 as G
import panel3c as PC
import pooled as PO
import propagator as PR
import truth3 as TR

PHASE3C = E.data_dir() / "phase3c"
P3B = E.data_dir() / "phase3" / "p3b"
OUT = Path(__file__).resolve().parents[2] / "data" / "phase3"
Z999 = 3.090232306167813
GENS = ("G0", "GT", "GF")
N_PANELS = 200
SEED = 20261104
EDGES = (0.0, 1 / 80, 1 / 40, 1 / 20, 1 / 10, 1 / 5)
SPECS = ([("all",)] + [("band", lo, hi) for lo, hi in zip(EDGES[:-1], EDGES[1:], strict=True)]
         + [("ceil", hi) for hi in EDGES[1:]])
BAND = {i + 1: EDGES[i + 1] for i in range(5)}  # mask index -> band upper edge
CEIL = {EDGES[i + 1]: 6 + i for i in range(5)}  # ceiling rho -> mask index
SENS_CELLS = ("Q-T1S2", "Q-T3S1")
BIAS_MAX, COV_MIN, MIN_EVAL = 0.02, 0.85, 0.5
NA = len(PO.A_GRID)


def q_cells() -> list[str]:
    return [c for c in E.cells() if c.startswith("Q-")]


def n_confirm(cell: str) -> int:
    _, u = E.read_universe()
    return sum(1 for c, sp in u.values() if c == cell and sp == "confirm")


def tail_model(cell: str) -> dict:
    import csv

    rows = [l for l in open(OUT / "p3b_signs.tsv") if not l.startswith("#")]
    r = next(x for x in csv.DictReader(rows, delimiter="\t") if x["cell"] == cell)
    return {"H": float(r["tail_H"]), "tau": float(r["tail_tau"])}


def setup(cell: str, gen: str, tag: str):
    """(Cell for generation, params, sign_model kind, explicit sign model or None)."""
    cl = G.load_cell(cell)
    params = {"beta": G.g0_beta(cl)} if gen == "G0" else (PC.gt_params(cl) if gen == "GT" else {})
    if tag == "main":
        return cl, params, "fgn_ar1", None
    return replace(cl, sign=tail_model(cell)), params, "two_component", tail_model(cell)


# ---- step 1: truth ------------------------------------------------------------------


def truth_quantities(Vd, Kd, Sd) -> dict:
    """The estimands on the truth's day-mean curves (a <= 1.5), and their day-jackknife SE."""
    a = TR.A_GRID
    sel = a <= 1.5
    lt2 = np.log2(TR.T_GRID.astype(float))
    mids = [int(np.flatnonzero(lt2 == np.floor((np.log2(lo) + np.log2(hi)) / 2))[0]) for lo, hi in TR.RANGES]

    h = a <= PO.A_RATIO
    w_pf, w_h = PO._ols_weights(a[sel]), PO._ols_weights(a[h])

    def fn(c):
        out = {}
        for r, (lo, hi) in enumerate(TR.RANGES):
            zV, zK = c[f"zV_{lo}_{hi}"][sel], c[f"zK_{lo}_{hi}"][sel]
            out[f"zV_{r}"], out[f"zK_{r}"] = zV, zK
            out[f"pf_{r}"] = zK @ w_pf
            out[f"slopeV_{r}"] = c[f"zV_{lo}_{hi}"][h] @ w_h
            out[f"slopeK_{r}"] = c[f"zK_{lo}_{hi}"][h] @ w_h
            out[f"ratio_{r}"] = out[f"slopeV_{r}"] / out[f"slopeK_{r}"]
            out[f"vk_diff_{r}"] = out[f"slopeV_{r}"] - 2 * out[f"slopeK_{r}"]
            ih, i0 = int(np.flatnonzero(a == 0.5)[0]), 0
            out[f"Rratio_{r}"] = c["R"][ih, mids[r]] / c["R"][i0, mids[r]]
        return out

    est = fn(TR.curves(Vd.mean(0), Kd.mean(0), Sd.mean(0)))
    se = TR.jackknife(Vd, Kd, Sd, fn=fn)
    return {k: (np.asarray(v, float), np.asarray(se[k], float)) for k, v in est.items()}


def gm_cap(name: str) -> float:
    d = PC.M.simulate_day(G.GM_PARAMS[name], 400_000, np.random.default_rng(SEED))
    return float(np.quantile(d.q, 0.999))


def run_truths(workers: int, log=print) -> None:
    tdir = PHASE3C / "truth"
    tdir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for c in q_cells():
        cl = G.load_cell(c)
        cap = float(np.exp(cl.mu + Z999 * cl.sigma_l))
        jobs.append(("G0", c, {"beta": G.g0_beta(cl), "q_cap": cap}, f"main_{c}_G0"))
        jobs.append(("GF", c, {"q_cap": cap}, f"main_{c}_GF"))
    for c in SENS_CELLS:
        cl, _, _, sm = setup(c, "G0", "tail")
        cap = float(np.exp(cl.mu + Z999 * cl.sigma_l))
        jobs.append(("G0", cl, {"beta": G.g0_beta(G.load_cell(c)), "q_cap": cap}, f"tail_{c}_G0"))
        jobs.append(("GT", cl, PC.gt_params(cl) | {"q_cap": cap}, f"tail_{c}_GT"))
        jobs.append(("GF", cl, {"sign_model": sm, "q_cap": cap}, f"tail_{c}_GF"))
    for name in ("GM", "GM_floor_G01"):
        jobs.append((name, "Q-T3S1", {"q_cap": gm_cap(name)}, f"main_cellfree_{name}"))
    log(f"truth: {len(jobs)} runs (+9 GT reused from 3b), about {len(jobs) * 0.6:.0f} min")
    for gen, cell, params, name in jobs:
        if (tdir / f"{name}.npz").exists() and "chi" in np.load(tdir / f"{name}.npz").files:
            continue
        r = TR.run_truth(gen, cell, params, name=name, min_days=16, chunk=16,
                         max_days=48 if gen.startswith("GM") else 64, workers=workers, out_dir=tdir, n_collapse=0,
                         log=lambda *_: None)
        log(f"{name}: {r.n_days} days, domain SE {r.se_max:.4f} (met {r.se_met}), {r.seconds:.0f}s")


def truth_path(tag: str, cell: str, gen: str) -> Path:
    if tag == "main" and gen == "GT":
        return P3B / "truth" / f"GT_R2_{cell}.npz"
    return PHASE3C / "truth" / f"{tag}_{cell}_{gen}.npz"


def load_truth(tag, cell, gen) -> tuple[dict, dict]:
    z = np.load(truth_path(tag, cell, gen))
    tq = truth_quantities(z["Vd"], z["Kd"], z["Sd"])
    meta = json.loads(str(z["meta"]))
    return tq, {"n_days": len(z["Vd"]), "se_max": meta.get("se_max", np.nan), "se_met": meta.get("se_met")}


# ---- step 2: the study ------------------------------------------------------------------


def ess_medians(p: PC.Panel) -> tuple[np.ndarray, np.ndarray]:
    """Median over stock-days of the Kish ESS of q_mom^(2a) and q_mom^a, per a."""
    sd = np.repeat(np.arange(len(p.lengths)), p.lengths)
    outV, outK = np.empty(NA), np.empty(NA)
    for i, a in enumerate(PO.A_GRID):
        for out, w in ((outV, p.q_mom ** (2 * a)), (outK, p.q_mom**a)):
            s1, s2 = np.bincount(sd, w), np.bincount(sd, w * w)
            out[i] = np.median(s1 * s1 / s2)
    return outV, outK


def study_job(args) -> str:
    tag, cell, gen, start, count = args
    path = PHASE3C / "study" / f"{tag}_{cell}_{gen}_{start:03d}.npz"
    if path.exists():
        return str(path)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    cl, params, kind, sm = setup(cell, gen, tag)
    nc = n_confirm(cell)
    rec = {}
    t0 = time.time()
    for i in range(start, start + count):
        rng = np.random.default_rng([SEED, q_cells().index(cell), GENS.index(gen), 0 if tag == "main" else 1, i])
        Nsd = PC.day_lengths(cell, nc, 8, rng)
        p = PC.make_panel(gen, cl, Nsd, rng, sign_model=kind, sm=sm, params=dict(params))
        pm = PO.panel_moments(p.eps, p.q_mom, p.lm, np.ones(len(p.eps), bool), p.starts, p.lengths)
        v, k, z, excl = PO.normalize(pm)
        sh = lambda x: PO.to_panel(x, 8)  # noqa: E731
        iv = PO.intervals(sh(v), sh(k), sh(z), p.Nsd, SPECS, rng)
        for key in ("zeta_V", "zeta_K", "pf", "ratio", "vk_diff", "Rratio"):
            e, lo, hi, _ = iv[key]
            for nm, arr in (("est", e), ("lo", lo), ("hi", hi)):
                rec.setdefault(f"{key}_{nm}", []).append(arr.astype(np.float32))
        eV, eK = ess_medians(p)
        rec.setdefault("essV", []).append(eV)
        rec.setdefault("essK", []).append(eK)
        rec.setdefault("excl", []).append(excl.mean(0))
        n_inc = PO.estimate(sh(v), sh(k), sh(z), p.Nsd, SPECS)["n_inc_V"][0]  # (nJ, nm) stocks per increment
        rec.setdefault("n_inc", []).append(n_inc.astype(np.float32))
        rec.setdefault("Ns", []).append(np.percentile(p.Nsd.min(1), [10, 50, 90]))
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, seconds=time.time() - t0, **{k: np.array(v) for k, v in rec.items()})
    return str(path)


def study_jobs(n_panels: int, chunk: int = 10) -> list[tuple]:
    jobs = []
    for tag, cells in (("main", q_cells()), ("tail", SENS_CELLS)):
        for c in cells:
            for g in GENS:
                jobs += [(tag, c, g, s, min(chunk, n_panels - s)) for s in range(0, n_panels, chunk)]
    # biggest panels first, so the pool ends evenly
    return sorted(jobs, key=lambda j: -n_confirm(j[1]))


def run_study(n_panels: int, workers: int, log=print) -> None:
    jobs = study_jobs(n_panels)
    todo = [j for j in jobs if not (PHASE3C / "study" / f"{j[0]}_{j[1]}_{j[2]}_{j[3]:03d}.npz").exists()]
    log(f"study: {len(jobs)} jobs of up to 10 panels, {len(todo)} to run on {workers} workers")
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for k, _ in enumerate(ex.map(study_job, todo), 1):
            if k % 20 == 0 or k == len(todo):
                el = time.time() - t0
                log(f"  {k}/{len(todo)} jobs, {el / 60:.1f} min, about {el / k * (len(todo) - k) / 60:.1f} min left")


def load_study(tag, cell, gen, n_panels) -> dict:
    parts = [np.load(PHASE3C / "study" / f"{tag}_{cell}_{gen}_{s:03d}.npz") for s in range(0, n_panels, 10)]
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0].files if k != "seconds"} | {
        "seconds": sum(float(p["seconds"]) for p in parts)}


# ---- step 2/3/4 analysis ----------------------------------------------------------------


ESTIMATORS = ("zeta_V", "zeta_K", "pf", "ratio", "vk_diff", "Rratio")
UNDEFINED = {("GF", "ratio")}  # truth slope of zeta_K is 0: the ratio has no finite truth


def truth_array(tq: dict, key: str) -> np.ndarray:
    """Truth shaped like the study's estimate without the mask axis: zeta (na, nr), others (nr,)."""
    if key in ("zeta_V", "zeta_K"):
        return np.stack([tq[f"z{key[-1]}_{r}"][0] for r in range(3)], axis=1)
    return np.array([float(tq[f"{key}_{r}"][0]) for r in range(3)])


def bias_coverage(st: dict, tq: dict, key: str, m: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(bias, coverage, share of panels with an estimate) per element under mask m."""
    e, lo, hi = (st[f"{key}_{x}"][..., m].astype(float) for x in ("est", "lo", "hi"))
    t = truth_array(tq, key)
    fin = np.isfinite(e)
    with np.errstate(invalid="ignore"):
        bias = np.nanmean(np.where(fin, e - t, np.nan), axis=0)
        cov = np.mean((lo <= t) & (t <= hi), axis=0)  # an undefined interval does not cover
    return bias, cov, fin.mean(0)


def rho_star(rows: list[dict]) -> float:
    """Largest band edge such that every evaluable band up to it passes; NaN when the first evaluable band
    fails."""
    best = np.nan
    for r in sorted(rows, key=lambda x: x["edge"]):
        if not r["evaluable"]:
            continue
        if not r["pass"]:
            break
        best = r["edge"]
    return best


def analyse(n_panels: int, log=print) -> dict:

    res = {"bands": [], "rho": [], "truth": [], "ess": [], "time": []}
    for tag, cells in (("main", q_cells()), ("tail", SENS_CELLS)):
        for c in cells:
            for g in GENS:
                tq, tm = load_truth(tag, c, g)
                st = load_study(tag, c, g, n_panels)
                res["time"].append({"tag": tag, "cell": c, "gen": g, "seconds": st["seconds"],
                                    "panels": len(st["zeta_V_est"])})
                if tag == "main":
                    row = {"cell": c, "gen": g} | tm
                    for r, (lo, hi) in enumerate(PO.RANGES):
                        rk = f"{lo}_{hi}"
                        row |= {f"zK_{rk}": " ".join(f"{x:.3f}" for x in tq[f"zK_{r}"][0]),
                                f"zV_{rk}": " ".join(f"{x:.3f}" for x in tq[f"zV_{r}"][0])}
                        for q in ("pf", "ratio", "vk_diff", "slopeV", "slopeK", "Rratio"):
                            row[f"{q}_{rk}"], row[f"{q}_{rk}_se"] = float(tq[f"{q}_{r}"][0]), float(tq[f"{q}_{r}"][1])
                    res["truth"].append(row)
                rows_by_est = {k: [] for k in ESTIMATORS}
                for m, edge in BAND.items():
                    for key in ESTIMATORS:
                        b, cv, fin = bias_coverage(st, tq, key, m)
                        if key in ("zeta_V", "zeta_K"):
                            b, cv, fin = b[PO.A_GRID <= 1.5], cv[PO.A_GRID <= 1.5], fin[PO.A_GRID <= 1.5]
                        ev = fin >= MIN_EVAL
                        evaluable = bool(ev.any())
                        ok = (np.abs(b) < BIAS_MAX) & (cv >= COV_MIN)
                        r = {"tag": tag, "cell": c, "gen": g, "estimator": key, "edge": edge,
                             "evaluable": evaluable, "n_elements": int(ev.sum()),
                             "pass": bool(np.all(ok[ev])) if evaluable else False,
                             "max_abs_bias": float(np.nanmax(np.abs(b[ev]))) if evaluable else np.nan,
                             "min_cov": float(np.nanmin(cv[ev])) if evaluable else np.nan,
                             "share_failing": float(np.mean(~ok[ev])) if evaluable else np.nan}
                        res["bands"].append(r)
                        rows_by_est[key].append(r)
                for key, rows in rows_by_est.items():
                    rs = np.inf if (g, key) in UNDEFINED else rho_star(rows)  # inf: not part of the minimum
                    res["rho"].append({"tag": tag, "cell": c, "gen": g, "estimator": key, "rho_star": rs})
                res.setdefault("study", {})[(tag, c, g)] = (st, tq)
    for name in ("GM", "GM_floor_G01"):  # truth only (cell-free)
        z = np.load(PHASE3C / "truth" / f"main_cellfree_{name}.npz")
        tq = truth_quantities(z["Vd"], z["Kd"], z["Sd"])
        meta = json.loads(str(z["meta"]))
        row = {"cell": "(all)", "gen": name, "n_days": len(z["Vd"]), "se_max": meta["se_max"], "se_met": meta["se_met"]}
        for r, (lo, hi) in enumerate(PO.RANGES):
            rk = f"{lo}_{hi}"
            row |= {f"zK_{rk}": " ".join(f"{x:.3f}" for x in tq[f"zK_{r}"][0]),
                    f"zV_{rk}": " ".join(f"{x:.3f}" for x in tq[f"zV_{r}"][0])}
            for q in ("pf", "ratio", "vk_diff", "slopeV", "slopeK", "Rratio"):
                row[f"{q}_{rk}"], row[f"{q}_{rk}_se"] = float(tq[f"{q}_{r}"][0]), float(tq[f"{q}_{r}"][1])
        res["truth"].append(row)
    for c in SENS_CELLS:  # tail-model truth rows
        for g in GENS:
            tq, tm = load_truth("tail", c, g)
            row = {"cell": f"{c} (tail)", "gen": g} | tm
            for r, (lo, hi) in enumerate(PO.RANGES):
                rk = f"{lo}_{hi}"
                row |= {f"zK_{rk}": " ".join(f"{x:.3f}" for x in tq[f"zK_{r}"][0]),
                        f"zV_{rk}": " ".join(f"{x:.3f}" for x in tq[f"zV_{r}"][0])}
                for q in ("pf", "ratio", "vk_diff", "slopeV", "slopeK", "Rratio"):
                    row[f"{q}_{rk}"], row[f"{q}_{rk}_se"] = float(tq[f"{q}_{r}"][0]), float(tq[f"{q}_{r}"][1])
            res["truth"].append(row)
    return res


# ---- step 3 ESS ---------------------------------------------------------------------------


def ess_threshold(points: list[tuple[float, float]], cov_min: float = COV_MIN) -> float:
    """Smallest ESS such that every point at or above it has coverage >= cov_min (NaN if none)."""
    pts = sorted(points)
    for i, (e, _) in enumerate(pts):
        if all(cv >= cov_min for _, cv in pts[i:]):
            return e
    return np.nan


def real_ess(cell: str) -> dict:
    import csv

    rows = [l for l in open(OUT / "p3b_ess.tsv") if not l.startswith("#")]
    out = {}
    for r in csv.DictReader(rows, delimiter="\t"):
        if r["cell"] == cell and r["rule"] == "R2":
            out[r["weights"]] = {float(k[4:]): float(v) for k, v in r.items() if k.startswith("med_")}
    return out


# ---- step 4 sigma_l at small N ------------------------------------------------------------


def sigma_bias(cell: str, Ns=(50, 100, 200, 400, 800), seed: int = SEED) -> dict:
    """Mean and sd over stocks of sigma_l-hat (N contiguous runs from each day with >= N runs, at a random
    offset) minus the stock's full-sample sigma_l. Sizes only."""
    z = np.load(P3B / f"raw_{cell}.npz")
    qs = np.split(np.log(z["q"].astype(float)), np.cumsum(z["N"])[:-1])
    sym = z["sym"]
    rng = np.random.default_rng(seed)
    out = {}
    for N in Ns:
        d = []
        for s in np.unique(sym):
            days = [x for x, t in zip(qs, sym, strict=True) if t == s]
            full = np.concatenate(days).std()
            sub = [x[o : o + N] for x in days if len(x) >= N for o in [rng.integers(0, len(x) - N + 1)]]
            if sub:
                d.append(np.concatenate(sub).std() - full)
        out[N] = (float(np.mean(d)), float(np.std(d)), len(d)) if d else (np.nan, np.nan, 0)
    return out


# ---- step 5 -------------------------------------------------------------------------------


TMIN_T = 2 ** np.arange(7)


def tmin_cell(cell: str, T_grid=TMIN_T) -> np.ndarray:
    """P(Delta_T != 0) per T, pooled over windows (stride T // 4, both mids valid) of the cell's
    exploration stock-days. Mids alone: no sign or size is read."""
    syms = E.explore_symbols(cell)
    hit, tot = np.zeros(len(T_grid)), np.zeros(len(T_grid))
    days, _ = E.read_universe()
    for d in days:
        runs = E.load_explore_runs(d, cell=cell)
        for s in syms:
            if s not in runs:
                continue
            mid2, ok = runs[s]["mid2"].astype(np.int64), runs[s]["mid_valid"].astype(bool)
            N = len(mid2)
            for j, T in enumerate(T_grid):
                if N <= T:
                    continue
                n = np.arange(0, N - T, max(1, T // 4))
                n = n[ok[n] & ok[n + T]]
                hit[j] += np.sum(mid2[n + T] != mid2[n])
                tot[j] += len(n)
    return hit / tot


def lmax_cell(cell: str, n_floor: int = 100) -> int:
    import csv

    rows = [l for l in open(OUT / "calibration_stocks.tsv") if not l.startswith("#")]
    nmin = np.array([int(r["N_min"]) for r in csv.DictReader(rows, delimiter="\t") if r["cell"] == cell])
    q25 = np.percentile(nmin[nmin > n_floor], 25)
    l = 1
    while 8 * (2 * l) <= q25:
        l *= 2
    return l


def beta_hat_bias(cell: str, l_max: int, n_panels: int = 20, n_floor: int = 100) -> dict:
    """G0 panels (confirmation-sized, day lengths from exploration stocks with N_min > n_floor): pooled C(k)
    and response increments S(l) to lag l_max, propagator.calibrate_G, fit_power on [1, l_max]."""
    cl = G.load_cell(cell)
    beta = G.g0_beta(cl)
    z = np.load(G.phase3_dir() / f"calibration_{cell}.npz")
    rows = np.array([z["N"][z["sym"] == s] for s in np.unique(z["sym"])])
    rows = rows[rows.min(1) > n_floor]
    nc = n_confirm(cell)
    est = []
    for i in range(n_panels):
        rng = np.random.default_rng([SEED, 5, q_cells().index(cell), i])
        Nsd = rows[rng.integers(0, len(rows), nc)]
        p = PC.make_panel("G0", cl, Nsd, rng, params={"beta": beta}, clip=None)
        sd = np.repeat(np.arange(len(p.lengths)), p.lengths)
        n = len(p.eps)
        Cn, Cd = np.zeros(l_max + 1), np.zeros(l_max + 1)
        Sn, Sd = np.zeros(l_max), np.zeros(l_max)
        r = np.diff(p.lm)  # r_n = m_{n+1} - m_n, valid within a day
        for l in range(l_max + 1):
            same = sd[: n - l] == sd[l:]
            Cn[l] += np.sum((p.eps[: n - l] * p.eps[l:])[same])
            Cd[l] += same.sum()
            if l < l_max:
                k = n - 1 - l
                ok = (sd[:k] == sd[l : l + k]) & (sd[l : l + k] == sd[l + 1 : l + k + 1])
                Sn[l] += np.sum((p.eps[:k] * r[l : l + k])[ok])
                Sd[l] += ok.sum()
        Gk = PR.calibrate_G(Cn / Cd, Sn / Sd)
        G0_, l0, b, _ = PR.fit_power(Gk, 1, l_max, l0s=np.array([1.0]))
        est.append(b)
    est = np.array(est)
    return {"beta": beta, "beta_hat_mean": float(est.mean()), "bias": float(est.mean() - beta),
            "sd": float(est.std(ddof=1)), "n_panels": n_panels}


# ---- report tables -------------------------------------------------------------------------


def _fmt(v) -> str:
    if isinstance(v, (bool, np.bool_)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return "" if not np.isfinite(v) else f"{v:.5g}"
    return str(v)


def write_tsv(path: Path, rows: list[dict], header: str) -> None:
    cols = list(rows[0])
    with open(path, "w") as f:
        f.write(f"# {header}\n" + "\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(r.get(c, "")) for c in cols) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["truth", "study", "analyse", "tmin", "semi"])
    ap.add_argument("--panels", type=int, default=N_PANELS)
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    a = ap.parse_args()
    t0 = time.time()
    if a.step == "truth":
        run_truths(a.workers)
    elif a.step == "study":
        run_study(a.panels, a.workers)
    elif a.step == "analyse":
        import phase3c_report as RP

        RP.write_all(analyse(a.panels), a.panels)
    elif a.step == "tmin":
        import phase3c_report as RP

        RP.write_tmin()
    elif a.step == "semi":
        import phase3c_report as RP

        RP.semi_smoke()
    print(f"{a.step}: {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
