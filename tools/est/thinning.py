"""Observe a subset of simulated flow: the single-venue problem.

A venue sees only some of the market's child executions, while the price is
formed by all of them. Thinning therefore selects trades and never touches the
mid: an observed trade carries the full-market mid read just before it (the
`mid` and `impact` of mbmodel.MBDay at that trade), and the metaorder table is
kept whole.

  random_thinning(day, p): keep each child execution independently with
    probability p.
  router_thinning(day, p, h, concentration): each metaorder routes to the venue
    with its own probability pi_i ~ Beta(kappa p, kappa (1 - p)) (mean p,
    variance p (1 - p) / (kappa + 1)); each of its executions is kept with
    probability pi_i; then a further fraction h of the kept executions is
    hidden (dropped independently with probability h), default h = 0.17. The
    observed share is p (1 - h) on average, and keep decisions of two children
    of one metaorder are correlated,

        ICC = Var(pi) (1 - h)^2 / (p (1 - h) (1 - p (1 - h))),

    against 0 for random thinning. kappa -> infinity recovers random thinning
    with probability p (1 - h); p = 1 makes pi = 1 for every metaorder.

Both return (observed_day, keep), keep a boolean mask over the input trades;
router_thinning also returns pi and the mask of hidden executions.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

TRADE_FIELDS = ("t", "eps", "q", "meta", "rank")


def _observe(day, keep: np.ndarray):
    idx = np.flatnonzero(keep)
    fields = {f: getattr(day, f)[idx] for f in TRADE_FIELDS}
    return replace(day, mid=day.mid[idx], impact=day.impact[idx], **fields)


def random_thinning(day, p: float, rng: np.random.Generator):
    if not 0 <= p <= 1:
        raise ValueError("p must be in [0, 1]")
    keep = rng.random(len(day.t)) < p
    return _observe(day, keep), keep


def router_thinning(day, p: float, rng: np.random.Generator, h: float = 0.17, concentration: float = 4.0):
    """Returns (observed_day, keep, pi, hidden)."""
    if not 0 <= p <= 1 or not 0 <= h < 1 or concentration <= 0:
        raise ValueError("need 0 <= p <= 1, 0 <= h < 1, concentration > 0")
    n_meta = len(day.meta_t0)
    if p >= 1 or p <= 0:
        pi = np.full(n_meta, float(p))
    else:
        pi = rng.beta(concentration * p, concentration * (1 - p), n_meta)
    n = len(day.t)
    routed = rng.random(n) < pi[day.meta]
    hidden = routed & (rng.random(n) < h)
    keep = routed & ~hidden
    return _observe(day, keep), keep, pi, hidden


def keep_icc(keep: np.ndarray, meta: np.ndarray) -> float:
    """Correlation of the keep decisions of consecutive children of the same metaorder."""
    order = np.lexsort((np.arange(len(meta)), meta))
    k, m = keep[order].astype(float), meta[order]
    same = m[1:] == m[:-1]
    return float(np.corrcoef(k[:-1][same], k[1:][same])[0, 1])
