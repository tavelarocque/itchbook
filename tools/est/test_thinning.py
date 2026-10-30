"""Tests for thinning.py on one simulated day at Paper II's preset (26372 trades,
9985 metaorders), then planted bugs. Tolerances are 4 sd of 50-seed spreads,
quoted in the docstrings.

    uv run --group dev --group analysis pytest tools/est/test_thinning.py
"""
from __future__ import annotations

from functools import cache
from typing import NamedTuple

import numpy as np
import pytest

import mbmodel as M
import thinning as TH

P, H, KAPPA = 0.3, 0.17, 4.0


@cache
def day():
    return M.simulate_day(M.PAPER2, 10000, np.random.default_rng(1))


# ---- checks -----------------------------------------------------------------------


def check_identity_at_p1_h0():
    d = day()
    for obs, keep in (TH.random_thinning(d, 1.0, np.random.default_rng(2)),
                      TH.router_thinning(d, 1.0, np.random.default_rng(3), h=0.0)[:2]):
        assert keep.all()
        for f in ("t", "eps", "q", "meta", "rank", "mid", "impact", "meta_t0", "meta_q", "meta_eps"):
            assert np.array_equal(getattr(obs, f), getattr(d, f)), f


def check_mid_is_the_full_market_mid():
    d = day()
    for obs, keep in (TH.random_thinning(d, P, np.random.default_rng(4)),
                      TH.router_thinning(d, P, np.random.default_rng(5), H, KAPPA)[:2]):
        idx = np.flatnonzero(keep)
        assert 0 < len(idx) < len(d.t)
        for f in ("t", "eps", "q", "meta", "rank", "mid", "impact"):
            assert np.array_equal(getattr(obs, f), getattr(d, f)[idx]), f
        assert np.array_equal(obs.meta_t0, d.meta_t0) and np.array_equal(obs.meta_dF, d.meta_dF)


def check_random_kept_fraction():
    """50 seeds: kept share 0.2998, sd 0.0021 (binomial 0.0028); tolerance 0.011."""
    _, keep = TH.random_thinning(day(), P, np.random.default_rng(6))
    assert abs(keep.mean() - P) < 0.011


def check_router_fractions():
    """50 seeds: observed share 0.2481 (target p (1 - h) = 0.2490), sd 0.0069; hidden share
    of routed executions 0.1705, sd 0.0046; mean pi 0.3000, sd 0.0020; Var(pi) 0.04199
    (theory p (1 - p) / (kappa + 1) = 0.042), sd 0.00054. Tolerances 4 sd."""
    obs, keep, pi, hidden = TH.router_thinning(day(), P, np.random.default_rng(7), H, KAPPA)
    assert abs(keep.mean() - P * (1 - H)) < 0.028
    assert abs(hidden.sum() / (hidden.sum() + keep.sum()) - H) < 0.018
    assert not np.any(keep & hidden)
    assert abs(pi.mean() - P) < 0.008 and abs(pi.var() - P * (1 - P) / (KAPPA + 1)) < 0.0022


def check_router_correlates_keep_decisions():
    """Correlation of keep decisions of consecutive children of one metaorder. 50 seeds:
    random -0.0010 sd 0.0054; router 0.1557 sd 0.0129, theory
    Var(pi) (1 - h)^2 / (p (1 - h) (1 - p (1 - h))) = 0.1547. Tolerances 4 sd."""
    d = day()
    _, keep = TH.random_thinning(d, P, np.random.default_rng(8))
    icc_random = TH.keep_icc(keep, d.meta)
    _, keep, _, _ = TH.router_thinning(d, P, np.random.default_rng(9), H, KAPPA)
    icc_router = TH.keep_icc(keep, d.meta)
    seen = P * (1 - H)
    theory = P * (1 - P) / (KAPPA + 1) * (1 - H) ** 2 / (seen * (1 - seen))
    assert abs(icc_random) < 0.022 and abs(icc_router - theory) < 0.052, (icc_random, icc_router)


def check_input_validation():
    d = day()
    for bad in (-0.1, 1.1):
        with pytest.raises(ValueError):
            TH.random_thinning(d, bad, np.random.default_rng(0))
    for kw in ({"p": 0.5, "h": 1.0}, {"p": 0.5, "h": -0.1}, {"p": 1.5}, {"p": 0.5, "concentration": 0.0}):
        with pytest.raises(ValueError):
            TH.router_thinning(d, rng=np.random.default_rng(0), **kw)
    obs, keep, pi, _ = TH.router_thinning(d, 0.0, np.random.default_rng(0))
    assert not keep.any() and np.all(pi == 0) and len(obs.t) == 0


CHECKS = {f.__name__: f for f in [
    check_identity_at_p1_h0, check_mid_is_the_full_market_mid, check_random_kept_fraction,
    check_router_fractions, check_router_correlates_keep_decisions, check_input_validation,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_thinning(name):
    CHECKS[name]()


# ---- planted bugs -----------------------------------------------------------------


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("router draws a keep probability per execution", "router_thinning",
        "routed = rng.random(n) < pi[day.meta]",
        "routed = rng.random(n) < rng.beta(concentration * p, concentration * (1 - p), n)",
        ("check_router_correlates_keep_decisions",)),
    Bug("hidden step keeps h instead of dropping h", "router_thinning",
        "hidden = routed & (rng.random(n) < h)", "hidden = routed & (rng.random(n) >= h)",
        ("check_router_fractions",)),
    Bug("Beta with mean p / (1 + p)", "router_thinning",
        "pi = rng.beta(concentration * p, concentration * (1 - p), n_meta)",
        "pi = rng.beta(concentration * p, concentration, n_meta)",
        ("check_router_fractions",)),
    Bug("random thinning keeps with probability 1 - p", "random_thinning",
        "keep = rng.random(len(day.t)) < p", "keep = rng.random(len(day.t)) >= p",
        ("check_random_kept_fraction", "check_identity_at_p1_h0")),
    Bug("observed mid taken by position, not by trade", "_observe",
        "mid=day.mid[idx]", "mid=day.mid[: len(idx)]",
        ("check_mid_is_the_full_market_mid",)),
    Bug("no special case at p = 1 (Beta(kappa, 0))", "router_thinning",
        "if p >= 1 or p <= 0:", "if p <= 0:",
        ("check_identity_at_p1_h0",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(TH, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
