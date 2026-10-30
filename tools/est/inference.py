"""Interval estimates for statistics of per-stock-day arrays.

Data are one or more arrays of shape (n_stocks, n_days, ...), NaN allowed, and
a statistic stat(*arrays) returning a scalar or a 1-D array; stat must itself
tolerate NaN (nanmean and the like). Stocks and days with no finite value in
any array carry no information and are dropped first, so n_days below counts
informative days only.

cluster_bootstrap (two-stage, crossed):
  1. draw n_days days with replacement, ONE draw shared by every stock, so a
     market-wide day shock moves all stocks together inside a replicate;
  2. within each stratum draw its stocks with replacement (stratum sizes fixed);
  3. recompute stat on X[s*][:, d*] and take percentile intervals
     (alpha/2, 1 - alpha/2). Replicates with a non-finite stat are dropped
     and counted.

jackknife_days (delete one day): theta_(i) is stat without day i, i = 1..n,

    se^2 = (n - 1)/n sum_i (theta_(i) - theta_(.))^2,   theta_hat +- t_{n-1}(1 - alpha/2) se,

centred on the full-sample theta_hat. It treats the stock set as fixed: a
stock effect common to all days never enters se.

preregistered(...) applies the pre-registered rule: at each level (0.90 and
0.95) report the WIDER of the two intervals (larger hi - lo), elementwise for
a vector statistic, with the method that supplied it.

Student t quantiles come from the regularized incomplete beta function
(continued fraction, Numerical Recipes 6.4), numpy and math only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

LEVELS = (0.90, 0.95)


# ---- Student t --------------------------------------------------------------------


def _betacf(a: float, b: float, x: float) -> float:
    tiny, eps = 1e-300, 3e-16
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            return h
    raise RuntimeError("incomplete beta continued fraction did not converge")


def betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lbt) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbt) * _betacf(b, a, 1.0 - x) / b


def t_cdf(t: float, df: float) -> float:
    tail = 0.5 * betainc(df / 2.0, 0.5, df / (df + t * t))
    return 1.0 - tail if t > 0 else tail


def t_ppf(p: float, df: float) -> float:
    """Quantile of Student t with df degrees of freedom (bisection on t_cdf)."""
    if not 0.0 < p < 1.0:
        raise ValueError("p must be in (0, 1)")
    if p < 0.5:
        return -t_ppf(1.0 - p, df)
    lo, hi = 0.0, 1.0
    while t_cdf(hi, df) < p:
        lo, hi = hi, 2.0 * hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, df) < p:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-13 * max(1.0, hi):
            break
    return 0.5 * (lo + hi)


# ---- data -------------------------------------------------------------------------


def _as_arrays(arrays) -> tuple[np.ndarray, ...]:
    if isinstance(arrays, np.ndarray):
        arrays = (arrays,)
    arrays = tuple(np.asarray(a, dtype=float) for a in arrays)
    shape = arrays[0].shape[:2]
    for a in arrays:
        if a.ndim < 2 or a.shape[:2] != shape:
            raise ValueError("every array needs the same leading (n_stocks, n_days) shape")
    return arrays


def informative(arrays) -> tuple[np.ndarray, np.ndarray]:
    """(stock_ok, day_ok): stocks and days with at least one finite value in some array."""
    arrays = _as_arrays(arrays)
    fin = np.zeros(arrays[0].shape[:2], dtype=bool)
    for a in arrays:
        f = np.isfinite(a)
        fin |= f.reshape(f.shape[0], f.shape[1], -1).any(axis=2)
    return fin.any(axis=1), fin.any(axis=0)


def _prepare(arrays, strata):
    arrays = _as_arrays(arrays)
    stock_ok, day_ok = informative(arrays)
    sk = stock_ok
    dk = day_ok
    arrays = tuple(a[sk][:, dk] for a in arrays)
    if strata is None:
        strata = np.zeros(int(sk.sum()), dtype=int)
    else:
        strata = np.asarray(strata)
        if strata.shape != stock_ok.shape:
            raise ValueError("strata needs one label per stock")
        strata = strata[sk]
    if arrays[0].shape[1] < 2:
        raise ValueError("need at least two informative days")
    return arrays, strata


# ---- results ----------------------------------------------------------------------


@dataclass
class Interval:
    method: str | np.ndarray  # 'bootstrap' / 'jackknife' (per component for preregistered)
    level: float
    lo: np.ndarray
    hi: np.ndarray

    @property
    def width(self) -> np.ndarray:
        return self.hi - self.lo


@dataclass
class BootstrapResult:
    estimate: np.ndarray
    replicates: np.ndarray  # (n_valid, k)
    n_failed: int
    intervals: dict[float, Interval] = field(default_factory=dict)


@dataclass
class JackknifeResult:
    estimate: np.ndarray
    leave_one_out: np.ndarray  # (n_days, k)
    se: np.ndarray
    df: int
    intervals: dict[float, Interval] = field(default_factory=dict)


# ---- estimators -------------------------------------------------------------------


def cluster_bootstrap(stat, arrays, n_boot: int = 999, strata=None, levels=LEVELS,
                      rng: np.random.Generator | None = None) -> BootstrapResult:
    rng = np.random.default_rng() if rng is None else rng
    arrays, strata = _prepare(arrays, strata)
    n_s, n_d = arrays[0].shape[:2]
    groups = [np.flatnonzero(strata == k) for k in np.unique(strata)]
    est = np.atleast_1d(np.asarray(stat(*arrays), dtype=float))
    reps = np.empty((n_boot, est.size))
    for b in range(n_boot):
        s = np.concatenate([g[rng.integers(0, len(g), len(g))] for g in groups])
        rows = s[:, None]
        cols = rng.integers(0, n_d, n_d)[None, :]
        reps[b] = np.atleast_1d(stat(*(a[rows, cols] for a in arrays)))
    ok = np.all(np.isfinite(reps), axis=1)
    reps = reps[ok]
    out = BootstrapResult(est, reps, int((~ok).sum()))
    for lev in levels:
        if len(reps) == 0:
            raise ValueError("no finite bootstrap replicate")
        lo, hi = np.quantile(reps, [(1 - lev) / 2, (1 + lev) / 2], axis=0)
        out.intervals[lev] = Interval("bootstrap", lev, lo, hi)
    return out


def jackknife_days(stat, arrays, levels=LEVELS) -> JackknifeResult:
    arrays, _ = _prepare(arrays, None)
    n = arrays[0].shape[1]
    est = np.atleast_1d(np.asarray(stat(*arrays), dtype=float))
    keep = ~np.eye(n, dtype=bool)
    loo = np.array([np.atleast_1d(stat(*(a[:, keep[i]] for a in arrays))) for i in range(n)], dtype=float)
    if not np.all(np.isfinite(loo)):
        raise ValueError("statistic is not finite with some day deleted")
    se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean(axis=0)) ** 2, axis=0))
    out = JackknifeResult(est, loo, se, n - 1)
    for lev in levels:
        c = t_ppf((1 + lev) / 2, n - 1)
        out.intervals[lev] = Interval("jackknife", lev, est - c * se, est + c * se)
    return out


def wider(a: Interval, b: Interval) -> Interval:
    """Elementwise the wider of two intervals at the same level (ties go to a)."""
    if a.level != b.level:
        raise ValueError("intervals at different levels")
    take_b = b.width > a.width
    return Interval(np.where(take_b, b.method, a.method), a.level,
                    np.where(take_b, b.lo, a.lo), np.where(take_b, b.hi, a.hi))


def preregistered(stat, arrays, n_boot: int = 999, strata=None, levels=LEVELS,
                  rng: np.random.Generator | None = None) -> dict[float, Interval]:
    """The pre-registered interval at each level: the wider of bootstrap and jackknife."""
    boot = cluster_bootstrap(stat, arrays, n_boot, strata, levels, rng)
    jack = jackknife_days(stat, arrays, levels)
    return {lev: wider(boot.intervals[lev], jack.intervals[lev]) for lev in levels}
