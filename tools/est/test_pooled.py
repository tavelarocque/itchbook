"""Tests for pooled.py, the primary pooled cell estimator, on constructed per-stock-day curves and
synthetic runs; no data files are read.

    uv run --group dev --group analysis pytest tools/est/test_pooled.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import inference as INF
import moments as MO
import pooled as PO

TG = np.array([1, 2, 4, 8, 16])
AG = np.array([0.0, 0.25, 0.5])
RG = ((2, 16),)
KW = dict(T_grid=TG, a_grid=AG, ranges=RG, a_max=0.5)


def curves(slopes, n_days=2, scale=None):
    """v = k = 2^(slope j) at T = 2^j, z = 0.1, per stock (rows) and day; slopes (S, D) or (S,)."""
    slopes = np.broadcast_to(np.asarray(slopes, float).reshape(len(slopes), -1), (len(slopes), n_days))
    j = np.arange(len(TG))
    v = 2.0 ** (slopes[:, :, None, None] * j[None, None, None, :]) * np.ones((1, 1, len(AG), 1))
    if scale is not None:
        v = v * np.asarray(scale, float)[:, :, None, None]
    return v, v.copy(), np.full(v.shape, 0.1)


def check_moments_match_day_moments():
    """Uncentred panel_moments equal moments.day_moments day by day; centred ones equal the direct
    centred means."""
    rng = np.random.default_rng(1)
    lengths = np.array([5, 40, 1, 300, 17])
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    n = lengths.sum()
    eps, q, lm = rng.choice([-1.0, 1.0], n), rng.lognormal(2, 1, n), np.cumsum(rng.standard_normal(n))
    ok = rng.random(n) > 0.05
    pm = PO.panel_moments(eps, q, lm, ok, starts, lengths)
    pc = PO.panel_moments(eps, q, lm, ok, starts, lengths, centered=True)
    for i, (s, L) in enumerate(zip(starts, lengths, strict=True)):
        m = MO.day_moments(eps[s:s + L], q[s:s + L], lm[s:s + L], ok[s:s + L], PO.T_GRID, PO.A_GRID)
        assert np.allclose(m.V, pm.V[i], equal_nan=True) and np.allclose(m.K, pm.K[i], equal_nan=True)
        assert np.allclose(m.S, pm.S[i], equal_nan=True) and np.array_equal(m.n_windows, pm.nw[i])
    s, L, j = starts[3], lengths[3], 3  # T = 8
    T = int(PO.T_GRID[j])
    w = np.arange(0, L - T, MO.stride(T))
    w = w[ok[s + w] & ok[s + w + T]]
    c = np.concatenate([[0.0], np.cumsum(eps[s:s + L] * q[s:s + L] ** 0.5)])
    I, D = c[w + T] - c[w], lm[s + w + T] - lm[s + w]
    assert np.isclose(pc.V[3, 2, j], np.var(I)) and np.isclose(pc.K[3, 2, j], np.mean((D - D.mean()) * (I - I.mean())))


def check_k1_exclusion():
    """A stock-day whose K_a(1) is negative is excluded from k (not from v), and counted."""
    rng = np.random.default_rng(2)
    lengths = np.array([400, 400])
    starts = np.array([0, 400])
    eps = rng.choice([-1.0, 1.0], 800)
    q = np.ones(800)
    lm = np.concatenate([[0.0], np.cumsum(eps[:-1])])  # mid moves with the previous sign: K(1) = 1
    lm[400:] = -lm[400:]  # second day: K(1) = -1
    lm[400:] -= lm[400]
    v, k, z, excl = PO.normalize(PO.panel_moments(eps, q, lm, np.ones(800, bool), starts, lengths))
    assert not excl[0].any() and excl[1].all()
    assert np.all(np.isnan(k[1])) and np.all(np.isfinite(k[0, :, :8])) and np.all(np.isfinite(v[1, :, :8]))


def check_mean_of_logs():
    """Two stocks with v = 2^(0.2 j) and 2^(0.6 j): the pooled slope is the mean of the stocks' log
    slopes, 0.4, so zeta = 1.4 exactly (the log of the cross-stock mean would not be a power law)."""
    v, k, z = curves([0.2, 0.6])
    e = PO.estimate(v, k, z, np.full((2, 2), 10**6), [("all",)], **KW)
    assert np.isclose(e["zeta_V"][0, 0, 0], 1.4) and np.isclose(e["zeta_K"][2, 0, 0], 1.4)


def check_equal_weights():
    """Days are averaged with equal weight within a stock, stocks with equal weight, whatever their run
    counts: stock 0 has days with slopes 0.2 and 0.6, stock 1 slope 1.0 on both days."""
    v, k, z = curves(np.array([[0.2, 0.6], [1.0, 1.0]]))
    Nsd = np.array([[100, 10**5], [10**4, 10**4]])
    j = np.arange(len(TG))
    s0 = np.log((2.0 ** (0.2 * j) + 2.0 ** (0.6 * j)) / 2)  # stock 0's log curve
    want = 1 + np.polyfit(np.log(TG[1:]), (s0 + np.log(2.0) * j)[1:] / 2, 1)[0]
    e = PO.estimate(v, k, z, Nsd, [("all",)], **KW)
    assert np.isclose(e["zeta_V"][0, 0, 0], want), (e["zeta_V"][0, 0, 0], want)


def check_ceiling_and_matching():
    """Stock 0 (N_s = 32) has slope 0.2 and a level 3 times stock 1's (N_s = 10^4, slope 0.6). With
    rho = 1/4 stock 0 is eligible up to T = 8: increments 1..8 average both stocks (0.4 per log2 step)
    and 8 -> 16 is stock 1's alone (0.6), exactly; zeta follows from those increments."""
    v, k, z = curves([0.2, 0.6], scale=np.array([[1, 1], [1, 1]]))
    v[0, :, :, 1:] *= 3.0  # a level shift of stock 0 above T = 1 (cancels in its own increments)
    Nsd = np.array([[32, 40], [10**4, 10**4]])
    e = PO.estimate(v, v.copy(), z, Nsd, [("ceil", 0.25)], **KW)
    inc = np.array([np.nan, 0.4, 0.4, 0.6]) * np.log(2)  # increments 2->4, 4->8, 8->16 (1->2 holds the shift)
    L = np.concatenate([[0.0], np.cumsum(inc[1:])])
    want = 1 + np.polyfit(np.log(TG[1:]), L, 1)[0]
    assert np.isclose(e["zeta_V"][0, 0, 0], want), (e["zeta_V"][0, 0, 0], want)
    m = PO.increment_masks(np.array([32, 10**4]), TG, [("ceil", 0.25)])[..., 0]
    assert m[0].tolist() == [True, True, True, False] and m[1].all()


def check_range_weights():
    """zeta - 1 from the increment weights equals the OLS slope of the chained curve over the range."""
    rng = np.random.default_rng(4)
    d = rng.standard_normal(len(PO.T_GRID) - 1)
    L = np.concatenate([[0.0], np.cumsum(d)])
    W = PO.range_weights()
    for r, (lo, hi) in enumerate(PO.RANGES):
        sel = (PO.T_GRID >= lo) & (PO.T_GRID <= hi)
        assert np.isclose(W[r] @ d, np.polyfit(np.log(PO.T_GRID[sel]), L[sel], 1)[0])


def check_jackknife_deletes_days():
    """Day shocks common to every stock and no stock differences: the day jackknife has SE > 0; stock
    differences and no day shocks: SE = 0 (deleting stocks would give the reverse)."""
    rng = np.random.default_rng(5)
    days = 0.3 + 0.1 * rng.standard_normal(8)
    v, k, z = curves(np.tile(days, (6, 1)), n_days=8)
    se, n = PO.jackknife_se(v, k, z, np.full((6, 8), 10**6), [("all",)], **KW)
    assert n == 8 and se["zeta_V"][0, 0, 0] > 0.01
    stocks = 0.3 + 0.1 * rng.standard_normal(6)
    v, k, z = curves(np.tile(stocks[:, None], (1, 8)), n_days=8)
    se, _ = PO.jackknife_se(v, k, z, np.full((6, 8), 10**6), [("all",)], **KW)
    assert se["zeta_V"][0, 0, 0] < 1e-12


def check_bootstrap_matches_inference():
    """The weighted bootstrap replicates equal inference.cluster_bootstrap's (same generator, same draw
    order) for zeta_V(0) over the range."""
    rng = np.random.default_rng(6)
    v, k, z = curves(0.3 + 0.1 * rng.standard_normal((5, 4)), n_days=4)
    Nsd = np.full((5, 4), 10**6)

    def stat(v_, k_, z_):
        return PO.estimate(v_, k_, z_, Nsd[: len(v_), : v_.shape[1]], [("all",)], **KW)["zeta_V"][0, 0, 0]

    ref = INF.cluster_bootstrap(stat, (v, k, z), n_boot=30, rng=np.random.default_rng(7))
    Wd, Ws = PO.bootstrap_weights(5, 4, 30, np.random.default_rng(7))
    mine = PO.core(v, k, z, Nsd, Wd, Ws, [("all",)], **KW)["zeta_V"][:, 0, 0, 0]
    assert np.allclose(mine, ref.replicates[:, 0])


CHECKS = {f.__name__: f for f in [
    check_moments_match_day_moments, check_k1_exclusion, check_mean_of_logs, check_equal_weights,
    check_ceiling_and_matching, check_range_weights, check_jackknife_deletes_days, check_bootstrap_matches_inference,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_pooled(name):
    CHECKS[name]()


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


UNMATCHED = '''    fin_l = np.isfinite(lv)
    Et = np.concatenate([E[:, :1], E], axis=1)  # eligibility at each horizon
    num_t = np.einsum("bs,stm,bsat->batm", ws, Et, np.where(fin_l, lv, 0.0))
    den_t = np.einsum("bs,stm,bsat->batm", ws, Et, fin_l.astype(float))
    mean_t = num_t / den_t
    return mean_t[:, :, 1:] - mean_t[:, :, :-1], den_t[:, :, 1:]
    inc = lv[..., 1:] - lv[..., :-1]'''

BUGS = [
    Bug("logs of raw cross-stock averages", "core",
        "lx = np.where(xs > 0, np.log(np.where(xs > 0, xs, 1.0)), np.nan)",
        "lx = np.broadcast_to(np.log(np.nanmean(xs, axis=1, keepdims=True)), xs.shape)",
        ("check_mean_of_logs",)),
    Bug("days and stocks weighted by run count", "core",
        "wd, ws = Wd[:, None, :] * np.ones(Nsd.shape)[None], Ws * 1.0",
        "wd, ws = Wd[:, None, :] * Nsd[None], Ws * Nsd.sum(axis=1)[None]",
        ("check_equal_weights",)),
    Bug("stock used above its ceiling", "increment_masks",
        "return masks(Ns, T_grid[1:], specs)", "return masks(Ns, T_grid[:-1], specs)",
        ("check_ceiling_and_matching",)),
    Bug("non-matched slopes (per-horizon means differenced)", "_increments",
        "    inc = lv[..., 1:] - lv[..., :-1]", UNMATCHED,
        ("check_ceiling_and_matching",)),
    Bug("jackknife deleting stocks instead of days", "jackknife_weights",
        "return 1.0 - np.eye(n_d), np.ones((n_d, n_s))", "return np.ones((n_s, n_d)), 1.0 - np.eye(n_s)",
        ("check_jackknife_deletes_days",)),
    Bug("K_a(1) exclusion ignored", "normalize",
        "    k[excl] = np.nan\n", "\n", ("check_k1_exclusion",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(PO, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
