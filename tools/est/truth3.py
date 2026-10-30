#!/usr/bin/env python3
"""Phase 3a ground truth: population curves of each generator by long simulation.

Days of N = DAY_N runs, independent (day j from default_rng([seed, j])), each
reduced by moments.day_moments to V_a(T), K_a(T), S(T) on T_GRID (powers of 2,
1..4096) and A_GRID. Population curves are the means over days:

  v_a(T) = V(T) / (T V(1)), k_a(T) = K(T) / (T K(1)), local slopes between
  adjacent T (moments.local_slopes), zeta_V(a), zeta_K(a) = 1 + OLS slope of
  log v (log k) on log T over each range in RANGES, R_a(T) = K / sqrt(S V),
  a*(T) = argmax over A_GRID of R_a(T) refined by a parabola through the
  maximum and its neighbours (endpoint maxima are not refined), chi_a = collapse
  exponent of I^a_T (collapse.chi_hat) over COLLAPSE_T on the first
  COLLAPSE_DAYS days.

Standard errors: delete-one-day jackknife of every reported slope. The run adds
days in chunks until every local slope of v and k in SE_DOMAIN (a <= 1, upper
T <= 128) has SE < target (default 0.005) or max_days is reached; se_met
reports which, and the SE of every slope is saved. CHOICE: the specified domain
(every a up to 3, T up to 4096) is out of reach. An 8-day pilot of G0 at Q-T3S1
gives SE 0.075 for the a = 0 slope 2048 -> 4096 (about 1800 days for 0.005) and
up to 1.19 at a = 3, where q^(2a) under log-normal sizes is carried by a few
runs per day. Each
chunk is checkpointed to PHASE3/truth/<name>.npz (per-day V, K, S, seeds) and a
run resumes from it. The first chunk is a pilot: the days it implies are printed
with a time estimate before any longer job.

    uv run --group analysis tools/est/truth3.py run GEN CELL [--param k=v ...]   # one generator, one cell
    uv run --group analysis tools/est/truth3.py scan [CELL]                      # GT scan (default Q-T3S1)
    uv run --group analysis tools/est/truth3.py all                              # every generator x cell
    uv run --group analysis tools/est/truth3.py summary                          # data/phase3/truth_summary.tsv
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import collapse as CO
import gen3 as G
import moments as MO

DAY_N = 10**6
T_GRID = 2 ** np.arange(13)
A_GRID = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0])
RANGES = ((8, 64), (16, 256), (64, 1024))
COLLAPSE_T = 2 ** np.arange(4, 11)  # 16 .. 1024
COLLAPSE_P = np.round(np.arange(0.05, 0.951, 0.05), 3)
COLLAPSE_DAYS = 4
TARGET_SE = 0.005
TRUTH_DIR = G.phase3_dir() / "truth"


def _day(args):
    gen, cell, N, seed, j, params, keep_series = args
    t0 = time.time()
    cell = cell if isinstance(cell, G.Cell) else G.load_cell(cell)
    params = dict(params)
    cap = params.pop("q_cap", None)  # clip the sizes the moments see; the generator's impact uses the raw q
    d = G.simulate_day(gen, cell, N, np.random.default_rng([seed, j]), **params)
    q = d.q if cap is None else np.minimum(d.q, cap)
    m = MO.day_moments(d.eps, q, d.lm, d.ok, T_GRID, A_GRID)
    series = (d.eps.astype(np.int8), d.q.astype(np.float32)) if keep_series else None
    return m.V, m.K, m.S, series, time.time() - t0


def simulate_days(gen, cell, N, seed, first, n, params, workers, n_series=0):
    jobs = [(gen, cell, N, seed, j, params, j < n_series) for j in range(first, first + n)]
    if workers <= 1:
        return [_day(a) for a in jobs]
    with ProcessPoolExecutor(workers) as ex:
        return list(ex.map(_day, jobs))


# ---- curves -----------------------------------------------------------------------


def ols_zeta(T: np.ndarray, y: np.ndarray, lo: int, hi: int) -> np.ndarray:
    """1 + OLS slope of log y on log T over lo <= T <= hi, along the last axis."""
    sel = (T >= lo) & (T <= hi)
    x = np.log(T[sel])
    x = x - x.mean()
    with np.errstate(invalid="ignore", divide="ignore"):
        ly = np.log(y[..., sel])
    return 1 + (ly - ly.mean(axis=-1, keepdims=True)) @ x / (x @ x)


def a_star(R: np.ndarray, a: np.ndarray) -> np.ndarray:
    """argmax_a R[:, j] per T, parabola-refined at interior maxima (a grid may be uneven)."""
    out = np.full(R.shape[1], np.nan)
    for j in range(R.shape[1]):
        r = R[:, j]
        if not np.any(np.isfinite(r)):
            continue
        i = int(np.nanargmax(r))
        if 0 < i < len(a) - 1 and np.all(np.isfinite(r[i - 1 : i + 2])):
            c = np.polyfit(a[i - 1 : i + 2], r[i - 1 : i + 2], 2)
            out[j] = -c[1] / (2 * c[0]) if c[0] < 0 else a[i]
        else:
            out[j] = a[i]
    return out


def curves(V: np.ndarray, K: np.ndarray, S: np.ndarray) -> dict:
    """Population curves from day-mean moments; V, K (na, nT), S (nT)."""
    T = T_GRID.astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = V / (T[None, :] * V[:, :1])
        k = K / (T[None, :] * K[:, :1])
        R = K / np.sqrt(S[None, :] * V)
    out = {"v": v, "k": k, "R": R, "slope_v": MO.local_slopes(T, v), "slope_k": MO.local_slopes(T, k),
           "a_star": a_star(R, A_GRID)}
    for lo, hi in RANGES:
        out[f"zV_{lo}_{hi}"] = ols_zeta(T, v, lo, hi)
        out[f"zK_{lo}_{hi}"] = ols_zeta(T, k, lo, hi)
    return out


SE_KEYS = ["slope_v", "slope_k"] + [f"z{x}_{lo}_{hi}" for lo, hi in RANGES for x in "VK"]


def jackknife(Vd: np.ndarray, Kd: np.ndarray, Sd: np.ndarray, fn=None) -> dict:
    """Delete-one-day jackknife SE of every curve in SE_KEYS, or of each value of fn(curves) (a dict);
    NaN where a replicate is undefined."""
    fn = (lambda c: {k: c[k] for k in SE_KEYS}) if fn is None else fn
    n = len(Vd)
    sV, sK, sS = Vd.sum(0), Kd.sum(0), Sd.sum(0)
    reps = {}
    for i in range(n):
        c = fn(curves((sV - Vd[i]) / (n - 1), (sK - Kd[i]) / (n - 1), (sS - Sd[i]) / (n - 1)))
        for k, v in c.items():
            reps.setdefault(k, []).append(v)
    return {k: np.sqrt((n - 1) / n * np.sum((np.array(r) - np.mean(r, axis=0)) ** 2, axis=0))
            for k, r in reps.items()}


def v_depth(zK: np.ndarray) -> float:
    """min over a of zeta_K(a) minus the smaller of zeta_K(0) and zeta_K(3): negative iff the minimum is
    interior and below both ends (a V); 0 when the minimum sits at an end."""
    return float(np.min(zK) - min(zK[0], zK[-1]))


def depths(c: dict) -> dict:
    return {f"depth_{lo}_{hi}": v_depth(c[f"zK_{lo}_{hi}"]) for lo, hi in RANGES}


SE_DOMAIN = (1.0, 128)  # (a_max, T_max): the stopping rule covers local slopes with a <= a_max, upper T <= T_max


def max_slope_se(se: dict, domain: tuple[float, int] | None = SE_DOMAIN) -> float:
    """Largest local-slope SE of v and k over the domain (None: all a and T); NaN (a slope undefined in
    some replicate) counts as infinite."""
    ia = slice(None) if domain is None else A_GRID <= domain[0]
    it = slice(None) if domain is None else T_GRID[1:] <= domain[1]
    s = np.concatenate([se["slope_v"][ia][:, it].ravel(), se["slope_k"][ia][:, it].ravel()])
    return float(np.inf if np.any(~np.isfinite(s)) else s.max())


CHI_GRID = np.round(np.arange(0.25, 2.5001, 0.005), 4)


def chi_by_a(series: list[tuple[np.ndarray, np.ndarray]]) -> np.ndarray:
    out = np.full(len(A_GRID), np.nan)
    for i, a in enumerate(A_GRID):
        xs = [e.astype(float) * q.astype(float) ** a for e, q in series]
        Q, dQ, iqr = CO.quantile_curves(xs, COLLAPSE_T, COLLAPSE_P)
        out[i] = CO.chi_hat(Q, dQ, iqr, COLLAPSE_T, COLLAPSE_P, CHI_GRID)[0]
    return out


# ---- the run ----------------------------------------------------------------------


def _cell_name(cell) -> str:
    return cell.name if isinstance(cell, G.Cell) else cell


@dataclass
class TruthResult:
    name: str
    n_days: int
    curves: dict
    se: dict
    se_max: float
    se_met: bool
    chi: np.ndarray
    seconds: float
    params: dict = field(default_factory=dict)
    days: tuple = ()  # (Vd, Kd, Sd), per-day moments


def run_truth(gen: str, cell: str, params: dict | None = None, name: str | None = None, seed: int = 20261102,
              N: int = DAY_N, target_se: float = TARGET_SE, min_days: int = 8, max_days: int = 96, chunk: int = 8,
              workers: int = 4, out_dir: Path | None = TRUTH_DIR, n_collapse: int = COLLAPSE_DAYS,
              budget_s: float | None = None, log=print) -> TruthResult:
    params = {} if params is None else dict(params)
    name = name or f"{gen}_{_cell_name(cell)}"
    path = None if out_dir is None else Path(out_dir) / f"{name}.npz"
    Vd, Kd, Sd, secs, series = [], [], [], 0.0, []
    if path is not None and path.exists():
        z = np.load(path, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        same = (meta["gen"], meta["cell"], meta["params"], meta["seed"], meta["N"])
        if same == (gen, _cell_name(cell), params, seed, N):
            Vd, Kd, Sd = list(z["Vd"]), list(z["Kd"]), list(z["Sd"])
            secs = float(meta["seconds"])
            log(f"{name}: resuming from {len(Vd)} days")
    t_start = time.time()
    n_done = len(Vd)
    need_series = n_collapse if n_done == 0 else 0
    while True:
        n = len(Vd)
        if n >= min_days:
            se = jackknife(np.array(Vd), np.array(Kd), np.array(Sd))
            se_max = max_slope_se(se)
            if se_max < target_se or n >= max_days:
                break
            if budget_s is not None and secs + time.time() - t_start > budget_s:
                break
        add = min(chunk, max_days - n)
        out = simulate_days(gen, cell, N, seed, n, add, params, workers, need_series - n)
        for V, K, S, s, _ in out:
            Vd.append(V)
            Kd.append(K)
            Sd.append(S)
            if s is not None:
                series.append(s)
        if n == n_done and len(Vd) >= min_days and n_done == 0:
            se0 = max_slope_se(jackknife(np.array(Vd), np.array(Kd), np.array(Sd)))
            per_day = (time.time() - t_start) / len(Vd) * min(workers, chunk)
            want = int(np.ceil(len(Vd) * (se0 / target_se) ** 2)) if np.isfinite(se0) else max_days
            log(f"{name}: pilot {len(Vd)} days, max slope SE {se0:.4f}; target {target_se} needs ~{want} days "
                f"(cap {max_days}), ~{min(want, max_days) * per_day / min(workers, chunk) / 60:.1f} min")
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            meta = {"gen": gen, "cell": _cell_name(cell), "params": params, "seed": seed, "N": N,
                    "seconds": secs + time.time() - t_start}
            np.savez(path, Vd=np.array(Vd), Kd=np.array(Kd), Sd=np.array(Sd), meta=json.dumps(meta))
    secs += time.time() - t_start
    Vd, Kd, Sd = np.array(Vd), np.array(Kd), np.array(Sd)
    cur = curves(Vd.mean(0), Kd.mean(0), Sd.mean(0))
    if not series and n_collapse > 0:  # resumed: regenerate the collapse days
        series = [s for *_, s, _ in simulate_days(gen, cell, N, seed, 0, n_collapse, params, workers, n_collapse)]
    chi = chi_by_a(series) if series else np.full(len(A_GRID), np.nan)
    res = TruthResult(name, len(Vd), cur, se, se_max, bool(se_max < target_se), chi, secs, params, (Vd, Kd, Sd))
    if path is not None:
        save_result(res, Vd, Kd, Sd, path, gen, cell, seed, N)
    return res


def save_result(res: TruthResult, Vd, Kd, Sd, path: Path, gen, cell, seed, N) -> None:
    meta = {"gen": gen, "cell": _cell_name(cell), "params": res.params, "seed": seed, "N": N, "seconds": res.seconds,
            "n_days": res.n_days, "se_max": res.se_max, "se_met": res.se_met, "target_se": TARGET_SE}
    arrays = {f"curve_{k}": v for k, v in res.curves.items()} | {f"se_{k}": v for k, v in res.se.items()}
    np.savez(path, Vd=Vd, Kd=Kd, Sd=Sd, T=T_GRID, a=A_GRID, chi=res.chi, meta=json.dumps(meta), **arrays)


# ---- GT scan ----------------------------------------------------------------------

SCAN_GRID = {"beta_max": (0.3, 0.5, 0.7), "c": (0.1, 0.2, 0.4), "l_off": (1.0, 2.0, 3.0), "psi": (0.25, 0.5)}
SCAN_DIR = G.phase3_dir() / "gt_scan"


def gt_scan(cell: str = "Q-T3S1", n_days: int = 16, workers: int = 4, log=print) -> list[dict]:
    """GT on SCAN_GRID at one cell, l_star = mu + l_off sigma_l; n_days days each (checkpointed, resumable).
    Each row: the setting, zeta_K(a) on every range, the V depth and its jackknife SE."""
    import itertools

    cl = G.load_cell(cell)
    rows = []
    keys = list(SCAN_GRID)
    for vals in itertools.product(*SCAN_GRID.values()):
        s = dict(zip(keys, vals, strict=True))
        params = {"beta_max": s["beta_max"], "c": s["c"], "l_star": cl.mu + s["l_off"] * cl.sigma_l, "psi": s["psi"]}
        name = f"GT_{cell}_b{s['beta_max']}_c{s['c']}_l{s['l_off']}_p{s['psi']}"
        r = run_truth("GT", cell, params, name=name, min_days=n_days, max_days=n_days, chunk=n_days,
                      workers=workers, out_dir=SCAN_DIR, n_collapse=0, log=lambda *_: None)
        d = depths(r.curves)
        se = jackknife(*r.days, fn=depths)
        row = s | {"l_star": params["l_star"], "n_days": r.n_days, "seconds": r.seconds}
        for lo, hi in RANGES:
            k = f"{lo}_{hi}"
            row[f"zK_{k}"] = r.curves[f"zK_{k}"]
            row[f"zK_se_{k}"] = r.se[f"zK_{k}"]
            row[f"depth_{k}"], row[f"depth_se_{k}"] = d[f"depth_{k}"], float(se[f"depth_{k}"])
            row[f"argmin_{k}"] = float(A_GRID[int(np.argmin(r.curves[f"zK_{k}"]))])
        rows.append(row)
        log(f"{name}: depth " + " ".join(f"{row[f'depth_{lo}_{hi}']:+.3f}({row[f'depth_se_{lo}_{hi}']:.3f})"
                                         for lo, hi in RANGES) + f"  {r.seconds:.0f}s")
    return rows


# ---- all generators x cells ------------------------------------------------------

GT_TRUTH = {"beta_max": 0.7, "c": 0.1, "l_off": 1.0, "psi": 0.5}  # chosen from the Q-T3S1 scan
CELL_FREE = ("GF_mb", "GM", "GM_floor_G0", "GM_floor_G01")  # own signs and sizes: one run serves every cell
REPO = Path(__file__).resolve().parents[2]
SUMMARY = REPO / "data" / "phase3" / "truth_summary.tsv"


def plan(cells: list[str]) -> list[tuple[str, str, dict, str, int]]:
    """(gen, cell, params, name, max_days), cheapest first."""
    out = [("G0", c, {}, f"G0_{c}", 64) for c in cells]
    out += [("GF", c, {}, f"GF_{c}", 64) for c in cells]
    for c in cells:
        cl = G.load_cell(c)
        p = {"beta_max": GT_TRUTH["beta_max"], "c": GT_TRUTH["c"], "l_star": cl.mu + GT_TRUTH["l_off"] * cl.sigma_l,
             "psi": GT_TRUTH["psi"]}
        out.append(("GT", c, p, f"GT_{c}", 64))
    out += [(g, "Q-T3S1", {}, f"{g}_cellfree", 48) for g in CELL_FREE]
    return out


def run_all(workers: int = 4, log=print) -> None:
    import explore

    jobs = plan(explore.cells())
    est = sum(md * {"G0": 0.37, "GF": 0.3, "GT": 0.6}.get(g, 2.0) + 12 for g, _, _, _, md in jobs)
    log(f"{len(jobs)} runs; at most {est / 60:.0f} min if every run reaches its day cap")
    for gen, cell, params, name, max_days in jobs:
        path = TRUTH_DIR / f"{name}.npz"
        if path.exists() and "chi" in np.load(path).files:
            log(f"{name}: done")
            continue
        t0 = time.time()
        r = run_truth(gen, cell, params, name=name, min_days=16, chunk=16, max_days=max_days, workers=workers,
                      log=log)
        log(f"{name}: {r.n_days} days, max domain slope SE {r.se_max:.4f} (met {r.se_met}), "
            f"zK[8,64] {np.round(r.curves['zK_8_64'], 3)}, {time.time() - t0:.0f}s")


def summarize(path: Path = SUMMARY) -> None:
    rows = []
    for f in sorted(TRUTH_DIR.glob("*.npz")):
        z = np.load(f)
        if "chi" not in z.files:
            continue
        meta = json.loads(str(z["meta"]))
        gen, cell = meta["gen"], ("(all)" if meta["gen"] in CELL_FREE else meta["cell"])
        se_all = max_slope_se({"slope_v": z["se_slope_v"], "slope_k": z["se_slope_k"]}, None)
        T = z["T"]
        for i, a in enumerate(z["a"]):
            r = {"generator": gen, "cell": cell, "a": a, "n_days": meta["n_days"], "se_met": int(meta["se_met"]),
                 "se_max_domain": meta["se_max"], "se_max_all": se_all}
            for lo, hi in RANGES:
                for x in "VK":
                    k = f"z{x}_{lo}_{hi}"
                    r[k], r[k + "_se"] = z[f"curve_{k}"][i], z[f"se_{k}"][i]
            for t in (1, 16, 256, 1024):
                r[f"R_T{t}"] = z["curve_R"][i, T == t][0]
            for t in (1, 16, 256, 1024):
                r[f"a_star_T{t}"] = z["curve_a_star"][T == t][0]
            r["chi"] = z["chi"][i]
            r["params"] = json.dumps(meta["params"], sort_keys=True).replace("\t", " ")
            rows.append(r)
    cols = list(rows[0])
    with open(path, "w") as f:
        f.write("# Phase 3a GROUND TRUTH by long simulation (tools/est/truth3.py): days of 10^6 runs, T = 1..4096 "
                "(powers of 2); zeta = 1 + OLS slope of log v (log k) on log T over [lo, hi]; *_se delete-one-day "
                "jackknife; se_met: every local slope of v and k with a <= 1 and T <= 128 has SE < 0.005 "
                "(se_max_domain), se_max_all over every a and T; a_star_T* repeats per a; cell (all) = generator "
                "with its own signs and sizes\n")
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(r[c]) for c in cols) + "\n")


def _fmt(v) -> str:
    if isinstance(v, (float, np.floating)):
        return "" if not np.isfinite(v) else f"{v:.5g}"
    return str(v)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["run", "scan", "all", "summary"])
    ap.add_argument("gen", nargs="?")
    ap.add_argument("cell", nargs="?", default="Q-T3S1")
    ap.add_argument("--param", action="append", default=[])
    ap.add_argument("--max-days", type=int, default=96)
    ap.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    a = ap.parse_args()
    if a.step == "all":
        run_all(a.workers)
        summarize()
        return 0
    if a.step == "summary":
        summarize()
        return 0
    if a.step == "scan":
        rows = gt_scan(a.cell, workers=a.workers)
        np.save(SCAN_DIR / f"scan_{a.cell}.npy", np.array(rows, dtype=object), allow_pickle=True)
        return 0
    params = {k: float(v) for k, v in (p.split("=") for p in a.param)}
    r = run_truth(a.gen, a.cell, params, max_days=a.max_days, workers=a.workers)
    print(f"{r.name}: {r.n_days} days, max slope SE {r.se_max:.4f} (met {r.se_met}), {r.seconds:.0f}s")
    for lo, hi in RANGES:
        print(f"  zK[{lo},{hi}]", np.round(r.curves[f"zK_{lo}_{hi}"], 3))
    return 0


if __name__ == "__main__":
    sys.exit(main())
