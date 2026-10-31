"""Tests for the Phase 3c-2 path of pooled.py: the k-tilde normalisation, horizon dropping, the Student t
day jackknife and the stock-cluster jackknife diagnostic. Hand-built moments; no data files.

    uv run --group dev --group analysis pytest tools/est/test_pooled2.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import inference as INF
import pooled as PO

T = PO.T_GRID2.astype(float)
NA, NT = len(PO.A_GRID2), len(PO.T_GRID2)
ALL = [("all",)]


def raw(K, V=None, Sm=None, n_days=None, N=10**6):
    """Raw from K (S, D, nT) (same for every a); V and Sm default to T (so v = 1 and S(1) = 1)."""
    K = np.asarray(K, float)
    S_, D_ = K.shape[:2]
    V = np.broadcast_to(T, (S_, D_, NT)) if V is None else np.asarray(V, float)
    Sm = np.broadcast_to(T, (S_, D_, NT)) if Sm is None else np.asarray(Sm, float)
    rep = lambda x: np.repeat(x[:, :, None, :], NA, axis=2)  # noqa: E731
    return PO.Raw(rep(V).copy(), rep(K).copy(), Sm.copy(), np.full((S_, D_, NA), np.nan), np.full((S_, D_), N))


def check_day_averaged_denominators():
    """k~ = mean_d K_d(T) / (T sqrt(mean_d S_d(1) mean_d V_d(1))): with V(1) = 1, 3 and S(1) = 2, 4 on the two
    days the denominator is sqrt(3 x 2), not the per-day sqrt(2) and sqrt(12)."""
    K = np.stack([2.0 * T**1.2, 1.0 * T**1.2])[None]
    V = np.stack([1.0 * T, 3.0 * T])[None]
    Sm = np.stack([2.0 * T, 4.0 * T])[None]
    r = raw(K, V, Sm)
    v, k = PO.stock_curves_tilde(r.V, r.K, r.Sm, np.ones((1, 1, 2)), PO.T_GRID2)
    assert np.allclose(k[0, 0, 0], 1.5 * T**1.2 / (T * np.sqrt(3.0 * 2.0)))
    assert np.allclose(v[0, 0, 0], 1.0)


def check_no_k1_exclusion():
    """A stock-day with negative K at every T (so K(1) < 0) stays in the day average."""
    K = np.stack([3.0 * T**1.2, -1.0 * T**1.2])[None]
    r = raw(K)
    _, k = PO.stock_curves_tilde(r.V, r.K, r.Sm, np.ones((1, 1, 2)), PO.T_GRID2)
    assert np.allclose(k[0, 0, 0], 1.0 * T**0.2)


def check_log_slope_equals_K_over_T():
    """With S and V varying with T but constant single-run moments, zeta_K equals 1 + the log slope of
    K(T) / T: K = c T^1.3 for two stocks with different c gives 1.3 on every range."""
    K = np.stack([np.stack([2.0 * T**1.3] * 2), np.stack([5.0 * T**1.3] * 2)])
    V = np.broadcast_to(T**1.6, (2, 2, NT))
    Sm = np.broadcast_to(3.0 * T**1.4, (2, 2, NT))
    e = PO.estimate2(raw(K, V, Sm), ALL, **PO.LEGACY_RULE)  # 2 stocks: below the range-reach rule's 5
    assert np.allclose(e["zeta_K"][:, :, 0], 1.3) and np.allclose(e["zeta_V"][:, :, 0], 1.6)
    assert not e["dropped_K"].any()


def check_nonpositive_pooled_drops_T():
    """Stock 0 has k~ = -10 at T = 8, so the pooled k~ there is negative: T = 8 is dropped (flagged), the
    [4, 32] curve joins T = 4, 16, 32 with both stocks, and no stock-day is excluded elsewhere. This is the
    estimator before 31 Oct 2026 (LEGACY_RULE); under the range-reach rule the shortened range is undefined
    (test_3c4.check_full_range_only)."""
    base = T**1.25
    K = np.stack([np.stack([base] * 2), np.stack([3 * base] * 2)])
    K[0, :, 3] = -10 * 8.0  # T = 8: k~ = K / T = -10
    e = PO.estimate2(raw(K), ALL, **PO.LEGACY_RULE)
    assert e["dropped_K"][0, 3, 0] and e["dropped_K"][:, [0, 1, 2, 4, 5, 6, 7, 8], 0].sum() == 0
    pts = [2, 4, 5]  # T = 4, 16, 32
    m = np.array([(K[0, 0, t] + K[1, 0, t]) / 2 / T[t] for t in pts])
    want = 1 + np.polyfit(np.log(T[pts]), np.log(m), 1)[0]
    assert np.isclose(e["zeta_K"][0, 0, 0], want), (e["zeta_K"][0, 0, 0], want)


def check_jackknife_t_quantile():
    """The day jackknife interval is est -+ t_{G-1} se, Student t with G - 1 df (G = 8: 1.8946)."""
    lo, hi = PO.jk_interval(np.array([0.0]), np.array([1.0]), 8, 0.90)
    assert np.isclose(hi[0], INF.t_ppf(0.95, 7)) and np.isclose(hi[0], 1.894579, atol=1e-5)
    assert np.isclose(lo[0], -hi[0])


def check_stock_jackknife():
    """Delete-one-stock: SE > 0 with stock differences only, 0 with day differences only."""
    rng = np.random.default_rng(1)
    g = 1.2 + 0.1 * rng.standard_normal(6)
    K = np.stack([np.stack([T**x] * 8) for x in g])
    se, n = PO.jackknife2(raw(K), ALL, unit="stock")
    assert n == 6 and se["zeta_K"][0, 0, 0] > 0.01
    d = 1.2 + 0.1 * rng.standard_normal(8)
    K = np.stack([np.stack([T**x for x in d])] * 6)
    se, _ = PO.jackknife2(raw(K), ALL, unit="stock")
    assert se["zeta_K"][0, 0, 0] < 1e-12
    se, n = PO.jackknife2(raw(K), ALL, unit="day")
    assert n == 8 and se["zeta_K"][0, 0, 0] > 0.01


def check_old_path_is_3c():
    """norm="old" reproduces the Phase 3c core on the same moments."""
    rng = np.random.default_rng(2)
    K = np.abs(1 + 0.2 * rng.standard_normal((4, 3, NT))) * T**1.1
    r = raw(K)
    r.K1_t[:] = 5.0
    e = PO.estimate2(r, ALL, norm="old")
    with np.errstate(invalid="ignore", divide="ignore"):
        v, k = r.V / (T * r.V[..., :1]), r.K / (T * r.K[..., :1])
    o = PO.estimate(v, k, r.z, r.Nsd, ALL, T_grid=PO.T_GRID2, a_grid=PO.A_GRID2, ranges=PO.RANGES2, a_max=PO.A_MAX2)
    assert np.allclose(e["zeta_K"], o["zeta_K"], equal_nan=True)


CHECKS = {f.__name__: f for f in [
    check_day_averaged_denominators, check_no_k1_exclusion, check_log_slope_equals_K_over_T,
    check_nonpositive_pooled_drops_T, check_jackknife_t_quantile, check_stock_jackknife, check_old_path_is_3c,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_pooled2(name):
    CHECKS[name]()


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("per-stock-day denominators instead of day-averaged", "stock_curves_tilde",
        "k = Kb / (T * np.sqrt(Sb[:, :, None, :1] * Vb[..., :1]))",
        "k = _stock_mean(K / (T * np.sqrt(Sm[:, :, None, :1] * V[..., :1])), wd)",
        ("check_day_averaged_denominators",)),
    Bug("stock-days excluded on the sign of K_a(1)", "stock_curves_tilde",
        "Vb, Kb, Sb = _stock_mean(V, wd), _stock_mean(K, wd), _stock_mean(Sm, wd)",
        "Vb, Kb, Sb = _stock_mean(V, wd), _stock_mean(np.where(K[..., :1] > 0, K, np.nan), wd), _stock_mean(Sm, wd)",
        ("check_no_k1_exclusion",)),
    Bug("denominator at T instead of T = 1 (log slope of k~ differs from K/T)", "stock_curves_tilde",
        "np.sqrt(Sb[:, :, None, :1] * Vb[..., :1])", "np.sqrt(Sb[:, :, None, :] * Vb[..., :1])",
        ("check_log_slope_equals_K_over_T",)),
    Bug("jackknife with the normal quantile", "jk_interval",
        "c = INF.t_ppf((1 + level) / 2, n - 1)", "c = 1.6448536269514722",
        ("check_jackknife_t_quantile",)),
    Bug("stock jackknife deleting days", "jackknife2",
        "Wd, Ws = np.ones((S_, D_)), 1.0 - np.eye(S_)", "Wd, Ws = jackknife_weights(S_, D_)",
        ("check_stock_jackknife",)),
    Bug("non-positive pooled k~ kept (logs of a negative mean)", "chain_zeta",
        "keep = pooled > 0", "keep = np.isfinite(pooled)",
        ("check_nonpositive_pooled_drops_T",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(PO, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
