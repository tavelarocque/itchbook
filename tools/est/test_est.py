"""Tests for the estimator library: exact properties first, then recovery of
known parameters on synthetic data (fixed seeds, tolerances set from measured
spread, see the module docstrings).

    uv run pytest tools/est
"""

from __future__ import annotations

import numpy as np
import pytest

import longmem as LM
import moments as MO
import propagator as P
import synth as Y

# ---- generators ------------------------------------------------------------------


@pytest.mark.parametrize("H", [0.6, 0.75, 0.9])
def test_fgn_has_the_fgn_autocovariance(H):
    rng = np.random.default_rng(10)
    x = np.stack([Y.fgn(1 << 15, H, rng) for _ in range(8)])
    var = np.mean(x * x)
    lag1 = np.mean(x[:, 1:] * x[:, :-1])
    assert abs(var - 1) < 0.08
    # under strong memory the sample variance and lag-1 covariance share one slowly
    # converging factor, so test their ratio (the sample correlation)
    assert abs(lag1 / var - (2 ** (2 * H - 1) - 1)) < 0.02


def test_fgn_signs_follow_the_arcsine_law():
    rng = np.random.default_rng(11)
    C = np.mean([P.sign_acf(Y.fgn_signs(1 << 17, 0.75, rng), 20) for _ in range(4)], axis=0)
    for k in (1, 5, 20):
        assert abs(C[k] - Y.fgn_sign_acf(np.array([k]), 0.75)[0]) < 0.012


def test_lmf_signs_are_balanced_and_persistent():
    rng = np.random.default_rng(12)
    e = Y.lmf_signs(1 << 15, 1.5, 5, rng)
    assert set(np.unique(e)) == {-1, 1}
    assert abs(e.mean()) < 0.2
    assert P.sign_acf(e, 1)[1] > 0.05


def test_propagator_mid_matches_direct_sum():
    rng = np.random.default_rng(13)
    x = rng.standard_normal(300)
    G = Y.power_kernel(300, 1.3, 4.0, 0.3)
    m = Y.propagator_mid(x, G)
    direct = np.array([sum(G[n - k] * x[k] for k in range(n)) for n in range(300)])
    assert np.allclose(m, direct, atol=1e-9)
    assert abs(m[0]) < 1e-12
    G2 = G.copy()
    G2[0] = 5.0  # G(0) must never enter: an event does not move the mid before it
    assert np.allclose(Y.propagator_mid(x, G2), m, atol=1e-9)


# ---- moments ----------------------------------------------------------------------


def brute_moments(eps, q, lm, ok, T, a):
    N = len(eps)
    Is, Ds = [], []
    for n in range(0, N - T, MO.stride(T)):
        if not (ok[n] and ok[n + T]):
            continue
        Is.append(sum(eps[k] * q[k] ** a for k in range(n, n + T)))
        Ds.append(lm[n + T] - lm[n])
    Is, Ds = np.array(Is), np.array(Ds)
    return np.mean(Is * Is), np.mean(Ds * Is), np.mean(Ds * Ds), len(Is)


def test_day_moments_match_brute_force():
    rng = np.random.default_rng(14)
    N = 400
    eps = rng.choice([-1, 1], N)
    q = rng.lognormal(3, 1, N)
    lm = np.cumsum(rng.standard_normal(N)) * 1e-4
    ok = rng.random(N) > 0.05
    T_grid, a_grid = [1, 2, 3, 7, 16, 50], [0.0, 0.5, 1.0, 1.7]
    m = MO.day_moments(eps, q, lm, ok, T_grid, a_grid)
    for j, T in enumerate(T_grid):
        for i, a in enumerate(a_grid):
            V, K, S, nw = brute_moments(eps, q, lm, ok, T, a)
            assert np.isclose(m.V[i, j], V) and np.isclose(m.K[i, j], K)
            assert np.isclose(m.S[j], S) and m.n_windows[j] == nw


def test_moments_are_not_centered():
    # a day with net buying: centering would remove the mean imbalance
    eps = np.ones(100)
    m = MO.day_moments(eps, np.ones(100), np.zeros(100), np.ones(100, bool), [10], [0.0])
    assert m.V[0, 0] == 100.0


def test_normalized_moments_equal_one_under_independence():
    rng = np.random.default_rng(15)
    N = 200_000
    eps = rng.choice([-1.0, 1.0], N)
    q = rng.lognormal(4, 1, N)
    G = Y.power_kernel(N, 1e-4, 1.0, 0.0)  # constant G: permanent impact
    lm = Y.propagator_mid(eps, G)
    T = MO.log_grid(256)
    v, k = MO.normalized(MO.day_moments(eps, q, lm, np.ones(N, bool), T, [0.0, 0.5, 1.0]))
    sel = (T >= 4) & (T <= 256)
    assert np.all(np.abs(v[:, sel] - 1) < 0.1)
    assert np.all(np.abs(k[:, sel] - 1) < 0.1)


def test_sign_imbalance_variance_scales_with_memory():
    rng = np.random.default_rng(16)
    H = 0.75  # gamma = 0.5: Var(I^0_T) ~ T^(2 - gamma), so v_0 slope ~ 1 - gamma = 0.5
    slopes = []
    for _ in range(3):
        eps = Y.fgn_signs(1 << 17, H, rng)
        T = np.array([16, 32, 64, 128, 256])
        v, _ = MO.normalized(MO.day_moments(eps, np.ones(len(eps)), np.zeros(len(eps)), np.ones(len(eps), bool),
                                            np.concatenate([[1], T]), [0.0]))
        slopes.append(np.polyfit(np.log(T), np.log(v[0, 1:]), 1)[0])
    assert abs(np.mean(slopes) - 0.5) < 0.08


def test_stride_and_window_count():
    assert [MO.stride(T) for T in (1, 3, 4, 8, 100)] == [1, 1, 1, 2, 25]
    N = 101
    m = MO.day_moments(np.ones(N), np.ones(N), np.zeros(N), np.ones(N, bool), [1, 8, 100], [0.0])
    # starts n = 0, s, 2s, ... with n + T <= N - 1
    assert list(m.n_windows) == [100, 47, 1]


def test_log_grid_and_clip():
    g = MO.log_grid(1000)
    assert g[0] == 1 and g[-1] == 1000 and np.all(np.diff(g) > 0)
    assert np.array_equal(MO.clip_volume([5, 50, 500], 10_000, 0.01), [5, 50, 100])
    assert np.array_equal(MO.clip_volume([5, 500], 10_000, None), [5, 500])


# ---- long memory ------------------------------------------------------------------


@pytest.mark.parametrize("H", [0.6, 0.75, 0.9])
def test_local_whittle_recovers_d_on_fgn(H):
    rng = np.random.default_rng(17)
    d = np.mean([LM.local_whittle(Y.fgn(1 << 16, H, rng))[0] for _ in range(10)])
    assert abs(d - (H - 0.5)) < 0.015


def test_local_whittle_ignores_the_mean():
    rng = np.random.default_rng(18)
    x = Y.fgn(1 << 14, 0.7, rng)
    assert np.isclose(LM.local_whittle(x)[0], LM.local_whittle(x + 3.0)[0], atol=1e-6)


def test_local_whittle_standard_error():
    d, se, m = LM.local_whittle(np.random.default_rng(19).standard_normal(10_000), power=0.65)
    assert m == int(10_000**0.65) and np.isclose(se, 1 / (2 * np.sqrt(m)))
    assert abs(d) < 4 * se


def test_gamma_is_one_minus_twice_d():
    x = Y.fgn(1 << 15, 0.75, np.random.default_rng(21))
    d, se, m = LM.local_whittle(x)
    g, gse, gm = LM.gamma_lw(x)
    assert np.isclose(g, 1 - 2 * d) and np.isclose(gse, 2 * se) and gm == m
    assert abs(g - 0.5) < 0.06


def test_pool_days():
    mean, se, Q, df = LM.pool_days([0.4, 0.6], [0.1, 0.1])
    assert np.isclose(mean, 0.5) and np.isclose(se, 0.1 / np.sqrt(2)) and np.isclose(Q, 2.0) and df == 1


# ---- propagator -------------------------------------------------------------------


def test_sign_acf_matches_brute_force_at_long_lags():
    rng = np.random.default_rng(22)
    e = rng.choice([-1.0, 1.0], 300) + 0.3  # nonzero mean: must not be removed
    C = P.sign_acf(e, 250)
    for k in (0, 1, 50, 249, 250):
        assert np.isclose(C[k], np.mean(e[: len(e) - k] * e[k:]))


def test_fit_power_is_exact_on_the_model():
    G = Y.power_kernel(2000, 0.7, 5.0, 0.25)
    G0, l0, beta, rss = P.fit_power(G, 2, 128)
    assert abs(beta - 0.25) < 0.006 and rss < 1e-6


@pytest.mark.parametrize("beta", [0.2, 0.4])
def test_propagator_recovered_from_signs_and_mids(beta):
    rng = np.random.default_rng(20)
    e = Y.fgn_signs(1 << 17, 0.75, rng).astype(float)
    G = Y.power_kernel(1 << 17, 1.0, 5.0, beta)
    m = Y.propagator_mid(e, G)
    Gh = P.calibrate_G(P.sign_acf(e, 1024), P.response_increments(e, m, 1024))
    assert np.allclose(Gh[1:5], G[1:5], atol=0.01)
    assert abs(P.fit_power(Gh, 2, 128)[2] - beta) < 0.02


def test_response_increments_skip_invalid_mids():
    e = np.array([1.0, 1, 1, 1])
    m = np.array([0.0, 1, 100, 50])
    ok = np.array([True, True, False, True])
    S = P.response_increments(e, m, 1, ok)
    assert S[0] == 1.0  # only the (0, 1) step is valid
