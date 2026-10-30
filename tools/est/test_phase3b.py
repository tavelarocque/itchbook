"""Tests for the Phase 3b pieces (signs3b.py, clip3b.py) on hand examples and synthetic signs; no
data files are read.

    uv run --group dev --group analysis pytest tools/est/test_phase3b.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import calib3 as C3
import clip3b as CL
import signs3b as SG


def check_kish_ess():
    """Equal weights give n; one dominant weight gives about 1; scale free."""
    assert np.isclose(CL.kish(np.ones(37)), 37)
    assert np.isclose(CL.kish(np.full(10, 4.2)), 10)
    w = np.concatenate([[1e6], np.ones(99)])
    assert 1 <= CL.kish(w) < 1.001
    assert np.isclose(CL.kish(np.array([1.0, 1.0, 2.0])), 16 / 6)
    V, K = CL.ess_table([np.array([1.0, 1.0, 2.0])], a_grid=np.array([0.0, 0.5]))
    assert np.allclose(V, [[3, 16 / 6]]) and np.allclose(K, [[3, (2 + np.sqrt(2)) ** 2 / 4]])


def check_r1_own_volume():
    """R1 caps each stock-day at 1 % of ITS OWN total run volume."""
    q1 = np.concatenate([np.ones(99), [50.0]])  # volume 149: cap 1.49
    q2 = np.full(10, 1000.0)  # volume 10000: cap 100
    c1, c2 = CL.clip_r1([q1, q2])
    assert np.allclose(c1, np.concatenate([np.ones(99), [1.49]]))
    assert np.allclose(c2, np.full(10, 100.0))
    assert np.allclose(CL.apply_rule("R1", [q1, q2], np.array(["A", "A"]))[1], c2)


def check_r2_per_stock():
    """R2 and R3 cap each stock at its own quantile over all of its stock-days, not the cell's."""
    rng = np.random.default_rng(1)
    qa = [rng.integers(1, 1000, 500).astype(float) for _ in range(3)]
    qb = [rng.integers(10**5, 10**6, 500).astype(float) for _ in range(2)]
    qs, st = qa + qb, np.array(["A", "A", "A", "B", "B"])
    for rule, p in (("R2", 0.999), ("R3", 0.99)):
        caps = CL.stock_caps(qs, st, p)
        assert np.isclose(caps["A"], np.quantile(np.concatenate(qa), p))
        assert np.isclose(caps["B"], np.quantile(np.concatenate(qb), p))
        out = CL.apply_rule(rule, qs, st)
        for q, qc, s in zip(qs, out, st, strict=True):
            assert np.allclose(qc, np.minimum(q, caps[s]))
        assert all(np.mean(qc < q) < 0.02 for q, qc in zip(qs, out, strict=True))


def check_centred_acf_within_day():
    """Each stock-day loses its own mean: two constant days of opposite sign have centred C = 0 at every
    lag (a global mean, 0 here, would leave C = 1); random days match the direct sum."""
    lags = np.array([1, 3, 10])
    C, _, _ = SG.centred_acf([np.ones(40), -np.ones(40)], lags)
    assert np.allclose(C, 0.0)
    rng = np.random.default_rng(2)
    xs = [rng.choice([-1.0, 1.0], n) + 0.3 * k for k, n in enumerate((20, 50, 80))]
    C, _, _ = SG.centred_acf(xs, lags)
    ys = [x - x.mean() for x in xs]
    want = [sum(float(y[:-l] @ y[l:]) for y in ys) / sum(len(y) - l for y in ys) for l in lags]
    assert np.allclose(C, want)


def check_two_component_recovery():
    """300 stock-days of 2000 runs from (H, tau) = (0.7, 0.2): the fit recovers both within 0.03
    (measured 0.000 and 0.0025 to 0.0075 on two seeds); the exact C(l) formula matches the pooled
    sample C within 0.01 at every lag; on the exact C of (0.7, 0.45) the grid fit is exact to its step."""
    lags = C3.sign_lags(200)
    rng = np.random.default_rng(3)
    C, _, _ = C3.pooled_acf(SG.simulate_two_component([2000] * 300, 0.7, 0.2, rng), lags)
    f = SG.fit_two_component(lags, C)
    assert abs(f["H"] - 0.7) < 0.03 and abs(f["tau"] - 0.2) < 0.03, f
    assert not f["edge_first"]
    assert np.max(np.abs(SG.two_component_C(lags, 0.7, 0.2) - C)) < 0.01
    f = SG.fit_two_component(lags, SG.two_component_C(lags, 0.7, 0.45))  # exact target, large day component
    assert abs(f["H"] - 0.7) < 0.006 and abs(f["tau"] - 0.45) < 0.006, f


def check_edge_flag():
    """An anti-persistent source (H = 0.3) sits on the lower edge of H_GRID: the flag fires, the wide
    grid recovers H, and the final fit is off the edge; a persistent source raises no flag."""
    lags = C3.sign_lags(200)
    f = SG.fit_two_component(lags, SG.two_component_C(lags, 0.3, 0.1))
    assert f["edge_first"] and not f["edge_final"] and abs(f["H"] - 0.3) < 0.005, f
    f = SG.fit_two_component(lags, SG.two_component_C(lags, 0.7, 0.1))
    assert not f["edge_first"] and abs(f["H"] - 0.7) < 0.005
    f = SG.fit_two_component(lags, SG.two_component_C(lags, 0.7, 0.1), H_grid=np.arange(0.505, 0.65, 0.005))
    assert f["edge_first"]  # the upper edge too


CHECKS = {f.__name__: f for f in [
    check_kish_ess, check_r1_own_volume, check_r2_per_stock, check_centred_acf_within_day,
    check_two_component_recovery, check_edge_flag,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_phase3b(name):
    CHECKS[name]()


class Bug(NamedTuple):
    id: str
    module: object
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("Kish ESS without the square of the sum", CL, "kish",
        "float(np.sum(w)) ** 2 / s2", "float(np.sum(w)) / s2", ("check_kish_ess",)),
    Bug("R1 on the whole sample's volume", CL, "clip_r1",
        "MO.clip_volume(q, float(np.sum(q)), R1_FRAC)",
        "MO.clip_volume(q, float(sum(np.sum(x) for x in qs)), R1_FRAC)", ("check_r1_own_volume",)),
    Bug("R2 quantile per cell, not per stock", CL, "stock_caps",
        "if t == s]", "if True]", ("check_r2_per_stock",)),
    Bug("centred C(l) on the global mean", SG, "centred_acf",
        "[x - x.mean() for x in series]", "[x - np.concatenate(series).mean() for x in series]",
        ("check_centred_acf_within_day",)),
    Bug("day mean drawn per run, not per stock-day", SG, "simulate_two_component",
        "tau * rng.standard_normal()", "tau * rng.standard_normal(int(N))", ("check_two_component_recovery",)),
    Bug("fit model without the 1 + tau^2 normalisation", SG, "_grid_fit",
        "/ (1 + t2), -1, 1))", "/ (1 + 0 * t2), -1, 1))", ("check_two_component_recovery",)),
    Bug("edge flag tests only the upper edge", SG, "fit_two_component",
        "edge_first = bool(H <= H_grid[0] or H >= H_grid[-1])", "edge_first = bool(H >= H_grid[-1])",
        ("check_edge_flag",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(bug.module, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
