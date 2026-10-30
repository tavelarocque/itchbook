"""Tests for inference.py: exact checks first (t quantiles, jackknife by hand,
percentile endpoints, strata, NaN), then coverage on synthetic panels with
known truth and day shocks, then planted bugs.

    uv run --group dev --group analysis pytest tools/est/test_inference.py
"""
from __future__ import annotations

import warnings
from fractions import Fraction
from typing import NamedTuple

import numpy as np
import pytest

import inference as I

# ---- checks (each is also run against planted bugs below) -------------------------


def check_t_quantiles():
    # reference values: R qt(), 9 significant digits
    for p, df, want in [(0.975, 7, 2.364624252), (0.95, 7, 1.894578605), (0.975, 1, 12.70620474),
                        (0.95, 2, 2.919985580), (0.975, 30, 2.042272456), (0.975, 1000, 1.962339081)]:
        assert abs(I.t_ppf(p, df) - want) < 1e-6 * want, (p, df)
    assert abs(I.t_ppf(0.025, 7) + 2.364624252) < 1e-6
    assert abs(I.t_cdf(2.364624252, 7) - 0.975) < 1e-9


def check_jackknife_hand_mean():
    # one stock, days 1, 2, 3, 10: leave-one-out means 5, 14/3, 13/3, 2, mean 4
    X = np.array([[1.0, 2.0, 3.0, 10.0]])
    j = I.jackknife_days(np.nanmean, X)
    assert np.allclose(j.leave_one_out[:, 0], [5, 14 / 3, 13 / 3, 2])
    se = np.sqrt(3 / 4 * (1 + 4 / 9 + 1 / 9 + 4))  # = s / sqrt(n) for a mean
    assert np.isclose(j.se[0], se) and np.isclose(se, np.std([1, 2, 3, 10], ddof=1) / 2)
    assert j.df == 3
    c = 3.182446305  # t_{3, 0.975}
    iv = j.intervals[0.95]
    assert np.isclose(iv.lo[0], 4 - c * se, atol=1e-7) and np.isclose(iv.hi[0], 4 + c * se, atol=1e-7)


def check_jackknife_hand_ratio():
    # nonlinear: ratio of sums K / V over days, worked with exact fractions
    K = np.array([[1.0, 2.0, 4.0]])
    V = np.array([[1.0, 1.0, 2.0]])
    j = I.jackknife_days(lambda k, v: np.nansum(k) / np.nansum(v), (K, V))
    loo = [Fraction(6, 3), Fraction(5, 3), Fraction(3, 2)]
    bar = sum(loo) / 3
    var = Fraction(2, 3) * sum((x - bar) ** 2 for x in loo)
    se = float(var) ** 0.5
    assert np.allclose(j.leave_one_out[:, 0], [float(x) for x in loo])
    assert np.isclose(j.se[0], se)
    c = 2.919985580  # t_{2, 0.95}
    iv = j.intervals[0.90]
    assert np.isclose(iv.lo[0], 1.75 - c * se, atol=1e-7) and np.isclose(iv.hi[0], 1.75 + c * se, atol=1e-7)


def check_percentile_endpoints():
    rng = np.random.default_rng(1)
    X = rng.standard_normal((10, 6)) + rng.standard_normal(6)
    b = I.cluster_bootstrap(np.nanmean, X, 499, rng=rng)
    for lev in (0.90, 0.95):
        lo, hi = np.quantile(b.replicates[:, 0], [(1 - lev) / 2, (1 + lev) / 2])
        assert np.isclose(b.intervals[lev].lo[0], lo) and np.isclose(b.intervals[lev].hi[0], hi)
    assert b.intervals[0.95].lo[0] <= b.intervals[0.90].lo[0] <= b.intervals[0.90].hi[0] <= b.intervals[0.95].hi[0]


def check_bootstrap_keeps_day_shocks():
    """A pure day effect, identical across 40 stocks: the bootstrap SD of the mean
    must be the day-level SD, sqrt((n-1)/n) s_d / sqrt(n), not 1/sqrt(40) of it.
    2000 replicates give the SD to about 2%; the tolerance is 10%."""
    rng = np.random.default_rng(2)
    days = rng.standard_normal(8)
    X = np.repeat(days[None, :], 40, axis=0)
    b = I.cluster_bootstrap(np.nanmean, X, 2000, rng=rng)
    want = np.sqrt(7 / 8) * np.std(days, ddof=1) / np.sqrt(8)
    assert abs(np.std(b.replicates[:, 0]) / want - 1) < 0.10


def check_strata_sizes_fixed():
    # stat = share of stratum-A stocks in the replicate; fixed at 12/30 when strata are respected
    strata = np.array([0] * 12 + [1] * 18)
    X = np.repeat((strata == 0).astype(float)[:, None], 5, axis=1)
    b = I.cluster_bootstrap(np.nanmean, X, 200, strata=strata, rng=np.random.default_rng(3))
    assert np.all(b.replicates[:, 0] == 12 / 30)


def first_stock_mean(a):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # mean of an empty slice
        return np.nanmean(a[0])


def check_nan_handling():
    rng = np.random.default_rng(4)
    X = rng.standard_normal((6, 7))
    X[rng.random(X.shape) < 0.15] = np.nan
    Xd = np.insert(X, 3, np.nan, axis=1)  # an all-NaN day
    Xs = np.insert(Xd, 2, np.nan, axis=0)  # and an all-NaN stock
    j, jd = I.jackknife_days(np.nanmean, X), I.jackknife_days(np.nanmean, Xs)
    assert jd.df == j.df == 6
    assert np.allclose(jd.leave_one_out, j.leave_one_out) and np.allclose(jd.se, j.se)
    b = I.cluster_bootstrap(np.nanmean, X, 300, rng=np.random.default_rng(5))
    bd = I.cluster_bootstrap(np.nanmean, Xs, 300, rng=np.random.default_rng(5))
    assert np.array_equal(b.replicates, bd.replicates) and b.n_failed == 0
    assert I.informative(Xs)[1].sum() == 7 and I.informative(Xs)[0].sum() == 6
    # a statistic that is NaN when day 0 is missing: failed replicates are dropped and counted
    Y = np.array([[5.0, np.nan, np.nan, np.nan], [1.0, 2.0, 3.0, 4.0]])
    b = I.cluster_bootstrap(first_stock_mean, Y, 2000, strata=[0, 1], rng=np.random.default_rng(6))
    assert len(b.replicates) + b.n_failed == 2000 and np.all(b.replicates == 5.0)
    p = (3 / 4) ** 4  # day 0 never drawn
    assert abs(b.n_failed / 2000 - p) < 4 * np.sqrt(p * (1 - p) / 2000)
    with pytest.raises(ValueError):
        I.jackknife_days(first_stock_mean, Y)


def check_wider_and_vector_stat():
    a = I.Interval("bootstrap", 0.95, np.array([0.0, 0.0]), np.array([1.0, 3.0]))
    b = I.Interval("jackknife", 0.95, np.array([-1.0, 0.5]), np.array([1.0, 1.5]))
    w = I.wider(a, b)
    assert list(w.method) == ["jackknife", "bootstrap"]
    assert np.array_equal(w.lo, [-1.0, 0.0]) and np.array_equal(w.hi, [1.0, 3.0])
    rng = np.random.default_rng(7)
    X = rng.standard_normal((8, 6, 2))
    r = I.preregistered(lambda x: np.nanmean(x, axis=(0, 1)), X, 199, rng=rng)
    assert r[0.95].lo.shape == (2,) and np.all(r[0.95].width >= r[0.90].width)


def panel(rng, R=40, C=8):
    """X = 0.5 + day shock (sd 1) + stock effect (sd 0.3) + noise (sd 1), 10% NaN."""
    X = 0.5 + rng.standard_normal(C)[None, :] + 0.3 * rng.standard_normal(R)[:, None] + rng.standard_normal((R, C))
    X[rng.random((R, C)) < 0.1] = np.nan
    return X


def coverage(reps=400, B=399, seed=8):
    rng = np.random.default_rng(seed)
    hits = {}
    for _ in range(reps):
        X = panel(rng)
        b = I.cluster_bootstrap(np.nanmean, X, B, rng=rng)
        j = I.jackknife_days(np.nanmean, X)
        for lev in I.LEVELS:
            for name, iv in [("boot", b.intervals[lev]), ("jack", j.intervals[lev]),
                             ("wider", I.wider(b.intervals[lev], j.intervals[lev]))]:
                hits[name, lev] = hits.get((name, lev), 0) + int(iv.lo[0] <= 0.5 <= iv.hi[0])
    return {k: v / reps for k, v in hits.items()}


def check_coverage_at_eight_days():
    """Coverage of the mean over 400 panels of 40 stocks x 8 days.

    Measured (4 seeds x 400 panels): jackknife 0.896 / 0.941 at 90 / 95%; percentile
    bootstrap 0.842 / 0.893; the wider interval was the jackknife in every panel.
    The bootstrap under-covers as small-sample theory says it must: with day shocks
    dominant its variance is (n-1)/n of the truth and it uses normal, not t_7,
    quantiles, so it covers P(|t_7| < z sqrt(7/8)) = 0.832 / 0.891. Tolerances are
    3.5 binomial SEs at 400 panels (0.04 to 0.06)."""
    cov = coverage()
    for lev in I.LEVELS:
        se = np.sqrt(lev * (1 - lev) / 400)
        assert abs(cov["jack", lev] - lev) < 3.5 * se, (lev, cov)
        assert abs(cov["wider", lev] - lev) < 3.5 * se, (lev, cov)
        z = {0.90: 1.644853627, 0.95: 1.959963985}[lev]
        theory = 2 * I.t_cdf(z * np.sqrt(7 / 8), 7) - 1
        assert abs(cov["boot", lev] - theory) < 3.5 * np.sqrt(theory * (1 - theory) / 400), (lev, cov)


CHECKS = {f.__name__: f for f in [
    check_t_quantiles, check_jackknife_hand_mean, check_jackknife_hand_ratio, check_percentile_endpoints,
    check_bootstrap_keeps_day_shocks, check_strata_sizes_fixed, check_nan_handling, check_wider_and_vector_stat,
    check_coverage_at_eight_days,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_inference(name):
    CHECKS[name]()


# ---- planted bugs -----------------------------------------------------------------


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("days drawn separately for each stock", "cluster_bootstrap",
        "cols = rng.integers(0, n_d, n_d)[None, :]", "cols = rng.integers(0, n_d, (len(s), n_d))",
        ("check_bootstrap_keeps_day_shocks", "check_coverage_at_eight_days")),
    Bug("stocks drawn ignoring strata", "cluster_bootstrap",
        "groups = [np.flatnonzero(strata == k) for k in np.unique(strata)]", "groups = [np.arange(n_s)]",
        ("check_strata_sizes_fixed",)),
    Bug("percentile at alpha instead of alpha/2", "cluster_bootstrap",
        "[(1 - lev) / 2, (1 + lev) / 2]", "[1 - lev, lev]",
        ("check_percentile_endpoints",)),
    Bug("jackknife variance with 1/(n(n-1)) instead of (n-1)/n", "jackknife_days",
        "(n - 1) / n *", "1 / (n * (n - 1)) *",
        ("check_jackknife_hand_mean", "check_jackknife_hand_ratio")),
    Bug("normal instead of t critical value", "jackknife_days",
        "c = t_ppf((1 + lev) / 2, n - 1)", "c = {0.9: 1.644853627, 0.95: 1.959963985}[lev]",
        ("check_jackknife_hand_mean", "check_jackknife_hand_ratio", "check_coverage_at_eight_days")),
    Bug("jackknife interval centred on the leave-one-out mean", "jackknife_days",
        "est - c * se, est + c * se", "loo.mean(axis=0) - c * se, loo.mean(axis=0) + c * se",
        ("check_jackknife_hand_ratio",)),
    Bug("all-NaN days kept in the day count", "_prepare",
        "dk = day_ok", "dk = np.ones_like(day_ok)",
        ("check_nan_handling",)),
    Bug("t cdf with t instead of t^2", "t_cdf",
        "df / (df + t * t)", "df / (df + t)",
        ("check_t_quantiles",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(I, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
