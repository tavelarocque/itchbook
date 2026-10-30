"""Tests for panel3c.py: batched generation against gen3, per-stock R2 caps, day lengths, and the
semi-synthetic generator (flow unchanged, GT with c = 0 is G0, response on raw sizes).

    uv run --group dev --group analysis pytest tools/est/test_panel3c.py
"""
from __future__ import annotations

import json
from typing import NamedTuple

import numpy as np
import pytest

import calib3 as C3
import gen3 as G
import mbmodel as M
import panel3c as PC

SM = {"w1": 0.1, "H": 0.85, "w2": 0.6, "phi": 0.5}
CELL = G.Cell("TEST", np.arange(100, 200), 3.0, 1.2, (0.3, 0.7), SM)


def flow(lengths, seed=1):
    rng = np.random.default_rng(seed)
    n = int(np.sum(lengths))
    return rng.choice([-1.0, 1.0], n), np.exp(3 + 1.2 * rng.standard_normal(n))


def check_batched_conv_equals_gen3():
    rng = np.random.default_rng(2)
    lengths = np.array([1, 7, 64, 65, 300, 1000])
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    x = rng.standard_normal(lengths.sum())
    for beta in (0.0, 0.4, 0.9):
        y = PC.batched_causal_conv(x, starts, lengths, beta)
        for s, L in zip(starts, lengths, strict=True):
            assert np.allclose(y[s:s + L], G.causal_conv(x[s:s + L], M.trade_kernel(L, beta)), atol=1e-10)


def check_batched_signs_acf():
    """fgn_ar1 signs pooled over 400 stock-days of 2000 runs: C(l) within 0.012 of the arcsine-law
    target at lags 1..100; two_component signs within 0.012 of its exact C(l)."""
    import signs3b as SG

    rng = np.random.default_rng(3)
    lengths = np.full(400, 2000)
    lags = C3.sign_lags(100)
    eps = PC.signs(lengths, "fgn_ar1", SM, rng)
    C, _, _ = C3.pooled_acf(np.split(eps, np.cumsum(lengths)[:-1]), lags)
    assert np.max(np.abs(C - G.sign_model_C(lags, **SM))) < 0.012
    eps = PC.signs(lengths, "two_component", {"H": 0.7, "tau": 0.2}, rng)
    C, _, _ = C3.pooled_acf(np.split(eps, np.cumsum(lengths)[:-1]), lags)
    assert np.max(np.abs(C - SG.two_component_C(lags, 0.7, 0.2))) < 0.012


def check_r2_caps_per_stock():
    lengths = np.array([300, 200, 50, 800])  # two stocks of two days
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    _, q = flow(lengths)
    q[300:500] *= 100  # stock 0, day 2 much larger
    caps = PC.stock_caps(q, starts, lengths, 2)
    a, b = np.quantile(q[:500], 0.999), np.quantile(q[500:], 0.999)
    assert np.allclose(caps[:500], a) and np.allclose(caps[500:], b) and a != b


def check_day_lengths_copied(tmp):
    """Each synthetic stock's days are one exploration stock's days, rows drawn with replacement."""
    rows = {"AAA": [5, 6, 7], "BBB": [100, 200, 300], "CCC": [9, 9, 9]}
    sym = np.repeat(list(rows), 3)
    N = np.concatenate(list(rows.values()))
    np.savez(tmp / "calibration_X.npz", N=N, sym=sym, row=json.dumps({}))
    saved = G.phase3_dir
    G.phase3_dir = lambda: tmp
    try:
        out = PC.day_lengths("X", 40, 3, np.random.default_rng(4))
    finally:
        G.phase3_dir = saved
    assert out.shape == (40, 3)
    assert all(list(r) in rows.values() for r in out.tolist())
    assert len({tuple(r) for r in out.tolist()}) == 3


def check_flow_unchanged():
    """semi_synthetic returns the real flow value for value; only the mid is new."""
    Nsd = np.array([[50, 80], [120, 30]])
    eps, q = flow(Nsd.ravel())
    for gen in ("G0", "GT", "GF"):
        p = PC.semi_synthetic(eps.copy(), q.copy(), Nsd, gen, CELL, np.random.default_rng(5))
        assert np.array_equal(p.eps, eps) and np.array_equal(p.q, q)
        assert np.all(p.q_mom <= p.q) and np.any(p.q_mom < p.q)


def check_gt_c0_is_g0():
    Nsd = np.array([[50, 80], [120, 30]])
    eps, q = flow(Nsd.ravel())
    for fresh in (True, False):
        a = PC.semi_synthetic(eps, q, Nsd, "G0", CELL, np.random.default_rng(6), params={"beta": 0.45}, fresh=fresh)
        b = PC.semi_synthetic(eps, q, Nsd, "GT", CELL, np.random.default_rng(6),
                              params={"beta_max": 0.45, "c": 0.0, "l_star": 4.0, "psi": 0.5}, fresh=fresh)
        assert np.array_equal(a.lm, b.lm)


def check_response_uses_raw_sizes():
    """The mid equals the G0 response computed from the raw sizes, while the moments' sizes are clipped."""
    Nsd = np.array([[200, 300], [250, 150]])
    eps, q = flow(Nsd.ravel())
    p = PC.semi_synthetic(eps, q, Nsd, "G0", CELL, np.random.default_rng(7), params={"beta": 0.4, "noise": 0.0},
                          fresh=True)
    want = PC.response("G0", CELL, eps, q, p.starts, p.lengths, None, {"beta": 0.4, "noise": 0.0})
    assert np.allclose(p.lm, want)
    clipped = PC.response("G0", CELL, eps, p.q_mom, p.starts, p.lengths, None, {"beta": 0.4, "noise": 0.0})
    assert not np.allclose(p.lm, clipped)


CHECKS = {f.__name__: f for f in [
    check_batched_conv_equals_gen3, check_batched_signs_acf, check_r2_caps_per_stock, check_day_lengths_copied,
    check_flow_unchanged, check_gt_c0_is_g0, check_response_uses_raw_sizes,
]}
NEEDS_TMP = {"check_day_lengths_copied"}


def run_check(name, tmp):
    return CHECKS[name](tmp) if name in NEEDS_TMP else CHECKS[name]()


@pytest.mark.parametrize("name", list(CHECKS))
def test_panel3c(name, tmp_path):
    run_check(name, tmp_path)


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("batched kernel off by one lag", "batched_causal_conv",
        "g[:h] = M.trade_kernel(h, beta)", "g[1:h] = M.trade_kernel(h, beta)[:-1]",
        ("check_batched_conv_equals_gen3",)),
    Bug("semi-synthetic returns the clipped sizes as the flow", "semi_synthetic",
        "return Panel(eps, q, q_mom, lm,", "return Panel(eps, q_mom, q_mom, lm,", ("check_flow_unchanged",)),
    Bug("GT bins sharing a beta convolved once each", "response",
        "for beta in np.unique(betas):", "for beta in betas:", ("check_gt_c0_is_g0",)),
    Bug("response on the clipped sizes", "semi_synthetic",
        "        lm = response(gen, cell, eps, q, starts, lengths, rng, p)",
        "        lm = response(gen, cell, eps, np.minimum(q, stock_caps(q, starts, lengths, Nsd.shape[1])), starts, "
        "lengths, rng, p)", ("check_response_uses_raw_sizes",)),
    Bug("R2 cap over the whole panel", "stock_caps",
        "stocks = np.repeat(np.arange(len(lengths) // n_days), n_days)", "stocks = np.zeros(len(lengths), int)",
        ("check_r2_caps_per_stock",)),
    Bug("day lengths drawn day by day across stocks", "day_lengths",
        "return rows[rng.integers(0, len(rows), n_stocks)]",
        "return rows[rng.integers(0, len(rows), (n_stocks, n_days)), np.arange(n_days)]",
        ("check_day_lengths_copied",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate, tmp_path):
    mutate(PC, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            run_check(name, tmp_path)
