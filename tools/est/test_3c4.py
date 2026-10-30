"""Tests for Phase 3c-4: independent stationary days (panel3c), the whole-stock bootstrap and the widest of
three intervals (pooled), the pass rule with the Monte Carlo band and the all-ceilings-pass rho* rule
(phase3c4). No data files are read.

    uv run --group dev --group analysis pytest tools/est/test_3c3.py tools/est/test_3c4.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import gen3 as G
import panel3c as PC
import phase3c4 as C4
import pooled as PO

STRONG = {"w1": 0.6, "H": 0.99, "w2": 0.0, "phi": 0.5}
CELL = G.Cell("TEST", np.arange(100, 200), 3.0, 1.0, (0.3, 0.7), STRONG)
T = PO.T_GRID2.astype(float)
NA, NT = len(PO.A_GRID2), len(T)


def check_days_independent():
    """New path: a stock's days are separate series, so with a strong slowly varying sign component the mean
    signs of its day 1 and day 2 are uncorrelated across 300 stocks (|r| < 0.25); the continuous path
    correlates them (> 0.5)."""
    Nsd = np.full((300, 2), 400)
    for cont, lo, hi in ((False, -0.25, 0.25), (True, 0.5, 1.0)):
        p = PC.make_panel("G0", CELL, Nsd, np.random.default_rng(2), params={"beta": 0.4}, burn=2000, clip=None,
                          continuous=cont)
        assert len(p.eps) == Nsd.sum()
        m = p.eps.reshape(300, 2, 400).mean(axis=2)
        r = np.corrcoef(m[:, 0], m[:, 1])[0, 1]
        assert lo < r < hi, (cont, r)


def check_burn_per_day():
    """Every stock-day gets its own burn-in: the series lengths are burn + N for every day of every stock."""
    Nsd = np.array([[50, 70, 20], [30, 90, 10]])
    s = PC.independent_series(Nsd, 500)
    assert s.shape == (6, 1) and np.array_equal(s.ravel(), 500 + Nsd.ravel())


def raw_from(slopes):
    """Raw moments with K = T^(1 + slope_sd) per stock-day (slopes (S, D)), V = T, S = T."""
    slopes = np.asarray(slopes, float)
    S_, D_ = slopes.shape
    K = T[None, None, :] ** (1 + slopes[:, :, None])
    rep = lambda x: np.repeat(x[:, :, None, :], NA, axis=2)  # noqa: E731
    V = np.broadcast_to(T, (S_, D_, NT)).copy()
    return PO.Raw(rep(V), rep(K), V.copy(), np.full((S_, D_, NA), np.nan), np.full((S_, D_), 10**6))


def check_whole_stock_bootstrap():
    """Whole-stock bootstrap: stock multiplicities sum to S and every day keeps weight 1; its replicates vary
    with stock differences only, and are constant with day differences only."""
    Wd, Ws = PO.whole_stock_weights(6, 8, 50, np.random.default_rng(3))
    assert np.all(Wd == 1) and np.allclose(Ws.sum(axis=1), 6) and Ws.std() > 0
    rng = np.random.default_rng(4)
    stock_only = raw_from(np.repeat(0.2 + 0.1 * rng.standard_normal((6, 1)), 8, axis=1))
    day_only = raw_from(np.repeat(0.2 + 0.1 * rng.standard_normal((1, 8)), 6, axis=0))
    for raw, varies in ((stock_only, True), (day_only, False)):
        Wd, Ws = PO.whole_stock_weights(6, 8, 40, np.random.default_rng(5))
        z = PO.core2(raw, Wd, Ws, [("all",)])["zeta_K"][:, 0, 0, 0]
        assert (z.std() > 1e-3) == varies, (varies, z.std())


def check_widest_of_three():
    """intervals3 reports, per element, the widest of its three candidates."""
    rng = np.random.default_rng(6)
    raw = raw_from(0.2 + 0.05 * rng.standard_normal((6, 8)))
    iv = PO.intervals3(raw, [("all",)], np.random.default_rng(7), n_boot=40, keys=("zeta_K",))
    e, lo, hi, meth = iv["zeta_K"]
    assert np.all(np.isin(meth, (1, 2, 3))) and np.all(hi > lo)
    se, n = PO.jackknife2(raw, [("all",)], "day")
    jlo, jhi = PO.jk_interval(e, se["zeta_K"], n, PO.LEVEL)
    assert np.all(hi - lo >= jhi - jlo - 1e-12)


def check_pass_rule_mc_band():
    """Coverage fails only below 0.85 by more than 2 MC SE; a bin within 2 MC SE of 0.85 is in the band."""
    assert C4.passes(0.0, 0.83, 0.02) and C4.passes(0.0, 0.90, 0.02)
    assert not C4.passes(0.0, 0.80, 0.02) and not C4.passes(0.03, 0.95, 0.01)
    assert C4.in_mc_band(0.83, 0.02) and C4.in_mc_band(0.88, 0.02) and not C4.in_mc_band(0.95, 0.02)


def check_rho_star_all_pass():
    """rho* = the largest ceiling with every evaluable ceiling up to it passing; not-evaluable ones are skipped."""
    assert C4.rho_star({1 / 80: True, 1 / 40: False, 1 / 20: True, 1 / 10: True}) == 1 / 80
    assert C4.rho_star({1 / 80: None, 1 / 40: True, 1 / 20: True, 1 / 10: False}) == 1 / 20
    assert np.isnan(C4.rho_star({1 / 80: False, 1 / 40: True}))


CHECKS = {f.__name__: f for f in [
    check_days_independent, check_burn_per_day, check_whole_stock_bootstrap, check_widest_of_three,
    check_pass_rule_mc_band, check_rho_star_all_pass,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_3c4(name):
    CHECKS[name]()


class Bug(NamedTuple):
    id: str
    module: object
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("days of a stock drawn from one series on the new path", PC, "make_panel",
        "return _independent_panel(gen, cell, Nsd, rng, sign_model, sm, clip, params, int(burn))",
        "return _stationary_panel(gen, cell, Nsd, rng, sign_model, sm, clip, params, int(burn))",
        ("check_days_independent",)),
    Bug("burn-in shared across a stock's days", PC, "independent_series",
        "return (burn + Nsd.reshape(-1, 1)).astype(np.int64)",
        "return (Nsd.reshape(-1, 1) + burn * (np.arange(Nsd.size) % Nsd.shape[1] == 0)[:, None]).astype(np.int64)",
        ("check_burn_per_day",)),
    Bug("whole-stock bootstrap resampling days instead of stocks", PO, "whole_stock_weights",
        "return np.ones((B, n_d)), Ws",
        "return np.stack([np.bincount(rng.integers(0, n_d, n_d), minlength=n_d) for _ in range(B)]).astype(float), "
        "np.ones((B, n_s))",
        ("check_whole_stock_bootstrap",)),
    Bug("rho* taken as the largest passing ceiling", C4, "rho_star",
        "        if not status[e]:\n            break", "        if not status[e]:\n            continue",
        ("check_rho_star_all_pass",)),
    Bug("Monte Carlo band ignored", C4, "passes",
        "~(np.asarray(cov) + 2 * np.asarray(mc_se) < COV_MIN)", "~(np.asarray(cov) < COV_MIN)",
        ("check_pass_rule_mc_band",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(bug.module, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
