"""Tests for surrogates.py: exact structure (multisets, bins, blocks, untouched
fields), then the statistics each surrogate must keep or destroy on synthetic
runs, then planted bugs.

Synthetic stock-day: N runs at uniform random times in 10:00-15:30 ET, log q an
AR(1) with coefficient phi, random signs, and a mid built by the propagator
from eps * q^0.5, so the size of a run predicts its own signed return.

    uv run --group dev --group analysis pytest tools/est/test_surrogates.py
"""
from __future__ import annotations

from collections import Counter
from typing import NamedTuple

import numpy as np
import pytest

import runs as R
import surrogates as SU
import synth as Y


def stock_day(seed: int, n: int, phi: float):
    rng = np.random.default_rng(seed)
    z = np.empty(n)
    e = rng.standard_normal(n)
    z[0] = e[0]
    for i in range(1, n):
        z[i] = phi * z[i - 1] + np.sqrt(1 - phi**2) * e[i]
    q = np.exp(4 + z)
    ts = np.sort(rng.integers(R.SESSION_LO, R.SESSION_HI, n))
    eps = rng.choice([-1.0, 1.0], n)
    mid = Y.propagator_mid(eps * q**0.5, Y.power_kernel(n, 1e-3, 1.0, 0.3))
    return q, ts, eps, mid, rng


# ---- checks -----------------------------------------------------------------------


def check_multisets_and_untouched_fields():
    q, ts, eps, mid, rng = stock_day(1, 3000, 0.8)
    for idx in (SU.s1_index(len(q), rng), SU.s2_index(ts, 7, rng)):
        assert np.array_equal(np.sort(idx), np.arange(len(q)))
        assert not np.array_equal(idx, np.arange(len(q)))
    assert np.array_equal(np.sort(SU.s1(q, rng)), np.sort(q))
    assert np.array_equal(np.sort(SU.s2(q, ts, 7, rng)), np.sort(q))
    runs = np.zeros(len(q), dtype=R.RUN_DTYPE)
    runs["locate"], runs["ts"], runs["sign"] = 7, ts, eps
    runs["q"] = np.round(q).astype(np.uint64) + 1
    runs["mid2"] = 20_000 + np.arange(len(q))
    for kind in ("S1", "S2"):
        out = SU.surrogate_runs(runs, kind, rng)
        for f in R.RUN_DTYPE.names:
            if f != "q":
                assert np.array_equal(out[f], runs[f]), (kind, f)
        assert np.array_equal(np.sort(out["q"]), np.sort(runs["q"])) and not np.array_equal(out["q"], runs["q"])


def check_s1_is_uniform():
    # every permutation of 4 runs equally likely: 24 cells, 12000 draws, chi^2 on 23 df (p = 1e-5 at 60)
    rng = np.random.default_rng(2)
    c = Counter(tuple(SU.s1_index(4, rng)) for _ in range(12000))
    assert len(c) == 24
    assert sum((v - 500) ** 2 / 500 for v in c.values()) < 60


def check_s1_destroys_size_autocorrelation():
    """Original lag-1 ACF of log q is 0.8. Measured over 20 seeds at N = 20000: S1's
    ACF at lags 1..12 has sd 0.0072 (= 1/sqrt(N)) and max |r| 0.022; tolerance 5/sqrt(N)."""
    q, ts, eps, mid, rng = stock_day(3, 20000, 0.8)
    assert SU.acf(np.log(q), 1)[1] > 0.78
    r = SU.acf(np.log(SU.s1(q, rng)), 12)[1:]
    assert np.all(np.abs(r) < 5 / np.sqrt(len(q))), r


def check_s2_keeps_bins():
    q, ts, eps, mid, rng = stock_day(4, 20000, 0.8)
    idx = SU.s2_index(ts, 8, rng)
    b = SU.clock_bins(ts)
    assert np.array_equal(b[idx], b)  # every run's volume comes from its own 30-minute bin
    qs = q[idx]
    for k in np.unique(b):
        assert np.array_equal(np.sort(qs[b == k]), np.sort(q[b == k]))
    assert len(np.unique(b)) == 11  # 10:00-15:30


def check_s2_keeps_blocks():
    """Exact: within each bin the surrogate is a reordering of the original blocks
    (cut every L runs from the bin's first run), each block intact and in order."""
    q, ts, eps, mid, rng = stock_day(5, 5000, 0.8)
    L = 6
    idx = SU.s2_index(ts, L, rng)
    b = SU.clock_bins(ts)
    for k in np.unique(b):
        lo, hi = np.flatnonzero(b == k)[[0, -1]]
        hi += 1
        blocks = {s: min(s + L, hi) for s in range(lo, hi, L)}
        p, seen = lo, set()
        while p < hi:
            s = int(idx[p])
            assert s in blocks and s not in seen, (k, p, s)
            e = blocks[s]
            assert np.array_equal(idx[p : p + e - s], np.arange(s, e))
            seen.add(s)
            p += e - s
        assert seen == set(blocks)
        assert not np.array_equal(idx[lo:hi], np.arange(lo, hi))


def check_s2_keeps_short_lag_autocorrelation():
    """S2 keeps the within-block share (1 - k/L) of the lag-k ACF and nothing beyond L.
    Measured over 20 seeds (N = 20000, phi = 0.8, L = 8): ratio to the original ACF
    0.874, 0.749, 0.625, 0.499 at lags 1..4 (sd 0.003 to 0.016), tolerance 0.06; at
    lags 8..12 the raw ACF is within the S1 noise band 5/sqrt(N)."""
    q, ts, eps, mid, rng = stock_day(6, 20000, 0.8)
    L = 8
    r0 = SU.acf(np.log(q), 12)
    r2 = SU.acf(np.log(SU.s2(q, ts, L, rng)), 12)
    for k in range(1, 5):
        assert abs(r2[k] / r0[k] - (1 - k / L)) < 0.06, (k, r2[k] / r0[k])
    assert np.all(np.abs(r2[L:]) < 5 / np.sqrt(len(q))), r2[L:]


def check_sign_size_link_destroyed():
    """corr(q_n^0.5, eps_n r_n), r_n = m_{n+1} - m_n, is 0.94 by construction.
    Measured over 20 seeds: S1 -0.004 (sd 0.007), S2 0.009 (sd 0.011, max 0.033).
    S2's small positive residual is real: each bin keeps its own size level, and
    the bin's mean impact moves with it. Tolerances: 0.035 (S1), 0.05 (S2)."""
    q, ts, eps, mid, rng = stock_day(7, 20000, 0.8)
    er = eps[:-1] * np.diff(mid)
    assert np.corrcoef(q[:-1] ** 0.5, er)[0, 1] > 0.9
    c1 = np.corrcoef(SU.s1(q, rng)[:-1] ** 0.5, er)[0, 1]
    c2 = np.corrcoef(SU.s2(q, ts, SU.choose_L(q, ts), rng)[:-1] ** 0.5, er)[0, 1]
    assert abs(c1) < 0.035 and abs(c2) < 0.05, (c1, c2)


def check_choose_L():
    """AR(1) log q with phi = 0.5: r_k = 0.5^k first falls below 0.05 at k = 5 (0.0625,
    then 0.031); at N = 100000 the ACF SE is 0.004 and 10 of 10 seeds gave 5. iid
    sizes give 1. phi = 0.99 on 2000 runs never falls below 0.05 before the cap,
    median runs per bin // 4."""
    q, ts, *_ = stock_day(8, 100_000, 0.5)
    assert SU.choose_L(q, ts) == 5
    q, ts, *_ = stock_day(9, 20000, 0.0)
    assert SU.choose_L(q, ts) == 1
    q, ts, *_ = stock_day(10, 2000, 0.99)
    cap = int(np.median(np.unique(SU.clock_bins(ts), return_counts=True)[1])) // 4
    assert 40 <= cap <= 50 and SU.choose_L(q, ts) == cap


def check_input_validation():
    rng = np.random.default_rng(11)
    with pytest.raises(ValueError):
        SU.s2_index(np.array([5, 3, 9]), 2, rng)
    with pytest.raises(ValueError):
        SU.s2_index(np.array([1, 2, 3]), 0, rng)
    with pytest.raises(ValueError):
        SU.choose_L(np.array([1.0, 0.0, 2.0]), np.array([1, 2, 3]))
    assert SU.clock_bins(np.array([R.SESSION_LO, R.SESSION_LO + SU.BIN_NS - 1, R.SESSION_LO + SU.BIN_NS]))\
        .tolist() == [20, 20, 21]


CHECKS = {f.__name__: f for f in [
    check_multisets_and_untouched_fields, check_s1_is_uniform, check_s1_destroys_size_autocorrelation,
    check_s2_keeps_bins, check_s2_keeps_blocks, check_s2_keeps_short_lag_autocorrelation,
    check_sign_size_link_destroyed, check_choose_L, check_input_validation,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_surrogates(name):
    CHECKS[name]()


# ---- planted bugs -----------------------------------------------------------------


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("S1 as a random circular shift", "s1_index",
        "return rng.permutation(n)", "return np.roll(np.arange(n), rng.integers(1, n))",
        ("check_s1_is_uniform", "check_s1_destroys_size_autocorrelation")),
    Bug("S2 blocks permuted across the whole day", "s2_index",
        "np.lexsort((np.arange(n), u[g], bins))", "np.lexsort((np.arange(n), u[g]))",
        ("check_s2_keeps_bins", "check_s2_keeps_blocks")),
    Bug("S2 shuffles runs, not blocks (L ignored)", "s2_index",
        "blk = pos // L", "blk = pos",
        ("check_s2_keeps_blocks", "check_s2_keeps_short_lag_autocorrelation")),
    Bug("S2 blocks cut from the day start, not the bin start", "s2_index",
        "blk = pos // L", "blk = np.arange(n) // L",
        ("check_s2_keeps_blocks",)),
    Bug("ACF not centred", "acf",
        "x = x - x.mean()", "x = x",
        ("check_choose_L",)),
    Bug("choose_L returns the last lag above the threshold", "choose_L",
        "int(below[0] + 1)", "int(below[0])",
        ("check_choose_L",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(SU, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
