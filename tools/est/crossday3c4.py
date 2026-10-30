#!/usr/bin/env python3
"""Phase 3c-4 step 1: cross-day dependence of exploration stocks' daily flow statistics (signs and sizes
only, read through explore.load_explore_runs; NASDAQ cells).

Per stock-day d of stock s: x = (a) mean sign, (b) mean log run size, (c) sign ACF at lag 10 centred on
the day's own mean sign (stock-days with at least MIN_N_ACF runs; centred so that (c) does not repeat (a)).

Cross-day correlation (CHOICE of estimator): for each pair of study days (i, j), the correlation across
the cell's stocks of x_si and x_sj after removing each day's cross-sectional mean (common day effects),
averaged over pairs; and the same averaged within calendar-distance groups of the pairs. A persistent
stock effect (a stock that is always larger, always bought) counts as cross-day dependence here: with 8
days per stock it cannot be told apart from slow time variation. 95 % interval: bootstrap over stocks
(B = 999), percentile.
"""

from __future__ import annotations

import sys
from datetime import date

import numpy as np

import explore as E

MIN_N_ACF = 50
LAG = 10
N_BOOT = 999
DIST_GROUPS = ((0, 45), (46, 200), (201, 400))  # calendar days between the two study days


def day_stats(runs) -> np.ndarray:
    eps = runs["sign"].astype(float)
    lq = np.log(runs["q"].astype(float))
    out = [eps.mean(), lq.mean(), np.nan]
    if len(eps) >= MIN_N_ACF:
        x = eps - eps.mean()
        out[2] = float(x[:-LAG] @ x[LAG:] / (len(x) - LAG) / (x @ x / len(x)))
    return np.array(out)


def collect(cells):
    days, _ = E.read_universe()
    X = {c: {} for c in cells}
    for k, d in enumerate(days):
        for c in cells:
            for s, r in E.load_explore_runs(d, cell=c).items():
                X[c].setdefault(s, np.full((len(days), 3), np.nan))[k] = day_stats(r)
    return days, {c: np.array([v for _, v in sorted(m.items())]) for c, m in X.items()}  # (S, D, 3)


def pair_corr(x: np.ndarray, pairs) -> float:
    """Mean over day pairs of the correlation across stocks, after removing day means."""
    x = x - np.nanmean(x, axis=0, keepdims=True)
    r = []
    for i, j in pairs:
        ok = np.isfinite(x[:, i]) & np.isfinite(x[:, j])
        if ok.sum() >= 5 and x[ok, i].std() > 0 and x[ok, j].std() > 0:
            r.append(np.corrcoef(x[ok, i], x[ok, j])[0, 1])
    return float(np.mean(r)) if r else np.nan


def main() -> int:
    cells = [c for c in E.cells() if c.startswith("Q-")]
    days, X = collect(cells)
    dt = [date(int(d[4:]), int(d[:2]), int(d[2:4])) for d in days]
    allp = [(i, j) for i in range(len(days)) for j in range(i + 1, len(days))]
    groups = {g: [(i, j) for i, j in allp if g[0] <= abs((dt[j] - dt[i]).days) <= g[1]] for g in DIST_GROUPS}
    rng = np.random.default_rng(20261108)
    names = ("mean sign", "mean log size", "centred sign ACF(10)")
    rows = []
    for c in cells:
        x = X[c]
        S = len(x)
        idx = [rng.integers(0, S, S) for _ in range(N_BOOT)]
        for q, nm in enumerate(names):
            est = pair_corr(x[:, :, q], allp)
            reps = [pair_corr(x[b, :, q], allp) for b in idx]
            lo, hi = np.nanpercentile(reps, [2.5, 97.5])
            row = {
                "cell": c,
                "stat": nm,
                "stocks": S,
                "corr": est,
                "lo": lo,
                "hi": hi,
                "excludes_0": int(lo > 0 or hi < 0),
            }
            for g, ps in groups.items():
                row[f"corr_{g[0]}_{g[1]}d"] = pair_corr(x[:, :, q], ps)
                row[f"pairs_{g[0]}_{g[1]}d"] = len(ps)
            rows.append(row)
            print(
                f"{c} {nm:22s} corr {est:+.3f} [{lo:+.3f}, {hi:+.3f}]  by distance "
                + " ".join(f"{row[f'corr_{g[0]}_{g[1]}d']:+.3f}({len(groups[g])})" for g in DIST_GROUPS),
                flush=True,
            )
    import phase3c as P3

    P3.write_tsv(
        P3.OUT / "p3c4_crossday.tsv",
        rows,
        "Phase 3c-4 cross-day correlation of daily flow statistics, "
        "exploration stocks (signs and sizes only), day means removed; 95 % stock bootstrap",
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
