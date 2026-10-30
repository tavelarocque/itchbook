"""Phase 3b clipping rules and the Kish effective sample size of size weights. Sizes alone.

Rules, applied to one cell's stock-days (q in shares, one array per stock-day):
  none  no clipping
  R1    q clipped at R1_FRAC = 1 % of that stock-day's own total run volume
        (moments.clip_volume(q, q.sum(), 0.01)): the current rule
  R2    q clipped at the stock's 99.9th percentile of run size over its 8 days
  R3    the same at the 99th percentile
Percentiles are numpy's default (linear interpolation) over all runs of the
stock, every study day pooled; one cap per stock, used on each of its days.

Removed shares, pooled over every run of the cell: runs clipped (q > cap),
sum q removed, and sum q^(2a) removed, sum (q^(2a) - min(q, cap)^(2a)) / sum q^(2a).

Kish ESS of weights w on one stock-day: (sum w)^2 / sum w^2, which is n for
equal weights and tends to 1 when one weight dominates. V moments weight runs by
q^(2a), K moments by q^a.
"""
from __future__ import annotations

import numpy as np

import moments as MO

R1_FRAC = 0.01
RULES = ("none", "R1", "R2", "R3")
QUANTILE = {"R2": 0.999, "R3": 0.99}
A_GRID = np.round(np.arange(0.0, 2.001, 0.25), 2)


def kish(w: np.ndarray) -> float:
    w = np.asarray(w, dtype=float)
    s2 = float(np.sum(w * w))
    return float(np.sum(w)) ** 2 / s2 if s2 > 0 else np.nan


def clip_r1(qs: list[np.ndarray]) -> list[np.ndarray]:
    return [MO.clip_volume(q, float(np.sum(q)), R1_FRAC) for q in qs]


def stock_caps(qs: list[np.ndarray], stocks: np.ndarray, p: float) -> dict[str, float]:
    """The p-quantile of run size per stock over all its stock-days."""
    stocks = np.asarray(stocks)
    return {s: float(np.quantile(np.concatenate([q for q, t in zip(qs, stocks, strict=True) if t == s]), p))
            for s in np.unique(stocks)}


def apply_rule(rule: str, qs: list[np.ndarray], stocks: np.ndarray) -> list[np.ndarray]:
    qs = [np.asarray(q, dtype=float) for q in qs]
    if rule == "none":
        return qs
    if rule == "R1":
        return clip_r1(qs)
    if rule in QUANTILE:
        caps = stock_caps(qs, stocks, QUANTILE[rule])
        return [np.minimum(q, caps[s]) for q, s in zip(qs, stocks, strict=True)]
    raise ValueError(f"unknown rule {rule!r}")


def removed_shares(qs, qcs, a_grid=A_GRID) -> dict:
    q, qc = np.concatenate(qs), np.concatenate(qcs)
    out = {"runs_clipped": float(np.mean(qc < q)), "q_removed": float(1 - qc.sum() / q.sum())}
    for a in a_grid:
        out[f"q2a_removed_{a:g}"] = float(1 - np.sum(qc ** (2 * a)) / np.sum(q ** (2 * a)))
    return out


def ess_table(qcs, a_grid=A_GRID) -> tuple[np.ndarray, np.ndarray]:
    """(ESS_V, ESS_K), each (n_stockdays, len(a_grid)): weights q^(2a) and q^a."""
    V = np.array([[kish(q ** (2 * a)) for a in a_grid] for q in qcs])
    K = np.array([[kish(q ** a) for a in a_grid] for q in qcs])
    return V, K


def a_max(median_ess: np.ndarray, threshold: float, a_grid=A_GRID) -> float:
    """Largest a on the grid whose median ESS meets the threshold (NaN when none does)."""
    ok = np.flatnonzero(np.asarray(median_ess) >= threshold)
    return float(a_grid[ok.max()]) if len(ok) else np.nan


def share_clipped_slope(shares: np.ndarray, Ns: np.ndarray, stocks: np.ndarray, n_boot: int,
                        rng: np.random.Generator) -> tuple[float, float, float]:
    """OLS slope of the per-stock-day share of runs clipped on log N, with a percentile 95 % interval from a
    cluster bootstrap over stocks (each draw keeps all of a stock's stock-days)."""
    x, y, stocks = np.log(np.asarray(Ns, float)), np.asarray(shares, float), np.asarray(stocks)
    slope = float(np.polyfit(x, y, 1)[0])
    ids = np.unique(stocks)
    rows = [np.flatnonzero(stocks == s) for s in ids]
    reps = []
    for _ in range(n_boot):
        idx = np.concatenate([rows[i] for i in rng.integers(0, len(ids), len(ids))])
        if np.ptp(x[idx]) > 0:
            reps.append(np.polyfit(x[idx], y[idx], 1)[0])
    lo, hi = np.percentile(reps, [2.5, 97.5])
    return slope, float(lo), float(hi)
