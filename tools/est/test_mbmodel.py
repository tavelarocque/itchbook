"""Tests for mbmodel.py.

1. Exact and distributional checks of each ingredient against the papers'
   formulas (Paper I = arXiv:2506.07711, Paper II = arXiv:2509.05065), on small
   days: parameters, metaorder table, child counts and schedule, the mid against
   direct sums of II Algorithm 1 line 16 in both clocks, the cross-metaorder
   sign correlation, the informed component.
2. Paper II's qualitative results at its settings, on pooled synthetic days,
   with tolerances from the seed spreads quoted in each docstring.
3. Planted bugs (on the cheap checks).

    uv run --group dev --group analysis pytest tools/est/test_mbmodel.py
"""
from __future__ import annotations

from dataclasses import replace
from typing import NamedTuple

import numpy as np
import pytest

import collapse as C
import mbmodel as M
import moments as MO
import propagator as P

P2 = M.PAPER2


def pooled_moments(days, T_grid, a_grid):
    V = K = S = 0
    for d in days:
        m = MO.day_moments(d.eps.astype(float), d.q, d.mid, np.ones(len(d.t), bool), T_grid, a_grid)
        V, K, S = V + m.V, K + m.K, S + m.S
    return V / len(days), K / len(days), S / len(days)


T_FIT = np.unique(np.concatenate([[1], np.round(np.logspace(1, 3, 9)).astype(int)]))  # 1, 10 .. 1000


def exponents_vs_a(days, a_grid):
    """Log-log slopes over T = 10..1000 of V_a(T) / V_a(1) and of K_a(T)."""
    V, K, _ = pooled_moments(days, T_FIT, a_grid)
    sel = T_FIT >= 10
    lt = np.log(T_FIT[sel])
    zV = np.array([np.polyfit(lt, np.log(V[i, sel] / V[i, 0]), 1)[0] for i in range(len(a_grid))])
    zK = np.array([np.polyfit(lt, np.log(K[i, sel]), 1)[0] for i in range(len(a_grid))])
    return zV, zK


# ---- ingredients ------------------------------------------------------------------


def check_parameters_and_configurations():
    assert np.isclose(P2.lam * P2.sigma_l**2, 1 / 8) and np.isclose(P2.lam_p, 2 * P2.lam)
    assert (P2.mu_m, P2.beta_m, P2.sigma_l, P2.n0, P2.m, P2.nu, P2.phi) == (1.5, 0.25, 1.0, 3.0, 3.0, 1.5e-3, 2e-3)
    assert np.isclose(P2.mu_1, 1.125) and np.isclose(P2.beta_1, 1.0)  # mu_m = mu_1 + lam m, beta_m = beta_1 - lam' m
    nc = M.configuration("NC-NVD-NVF")
    assert (nc.Gamma, nc.lam, nc.lam_p, nc.sigma_l, nc.m) == (0.0, 0.0, 0.0, 0.0, 0.0)
    assert M.configuration("C-VD-VF") == P2
    d = M.simulate_day(nc, 300, np.random.default_rng(1))
    assert np.all(d.q == 1.0) and np.all(d.meta_mu == 1.5) and np.all(d.meta_beta == 0.25)
    for bad in ("C-VD", "X-VD-VF", "C-VD-NV"):
        with pytest.raises(ValueError):
            M.configuration(bad)


def check_metaorder_table():
    d = M.simulate_day(P2, 20000, np.random.default_rng(2))
    n = len(d.meta_t0)
    assert np.all(np.diff(d.meta_t0) > 0) and d.meta_t0[0] >= 0 and d.meta_t0[-1] < d.t_end
    assert abs(n - 20000) < 4 * np.sqrt(20000)  # Poisson count
    ell = np.log(d.meta_q)
    assert np.all(d.meta_q >= 1)
    assert np.allclose(d.meta_mu, 1.125 + ell / 8)  # I (21): mu_q = mu_1 + lam log q
    assert np.allclose(d.meta_beta, np.clip(1.0 - ell / 4, 0.0, 0.99))  # I (62): beta_q = beta_1 - lam' log q
    assert abs(np.median(ell) - 3.0) < 0.05  # mode m of log q; truncation at 0 is 3 sigma away
    # truncation below q = 1 matters at m = 0.5: log q ~ N(0.5, 1) | log q >= 0 has mean 0.5 + phi(0.5) / Phi(0.5)
    d = M.simulate_day(replace(P2, m=0.5, lam=0.0, lam_p=0.0), 20000, np.random.default_rng(3))
    ell = np.log(d.meta_q)
    want = 0.5 + np.exp(-0.125) / np.sqrt(2 * np.pi) / 0.6914625
    assert ell.min() >= 0 and abs(ell.mean() - want) < 4 * 0.75 / np.sqrt(len(ell))


def check_child_counts():
    """P(s >= l | q) = l^-mu_q (discrete I (14), s0 = 1): compare the share of metaorders
    with s >= l to the mean of l^-mu_i, within 4 binomial SE."""
    d = M.simulate_day(P2, 40000, np.random.default_rng(4))
    for s in (2, 5, 20):
        want = np.mean(float(s) ** -d.meta_mu)
        got = np.mean(d.meta_n >= s)
        assert abs(got - want) < 4 * np.sqrt(want * (1 - want) / len(d.meta_n)), (s, got, want)
    assert d.meta_n.min() == 1


def check_child_schedule():
    d = M.simulate_day(P2, 5000, np.random.default_rng(5))
    assert np.all(np.diff(d.t) >= 0) and d.t[0] >= 0 and d.t[-1] < d.t_end
    assert np.array_equal(d.eps, d.meta_eps[d.meta]) and np.array_equal(d.q, d.meta_q[d.meta])
    counts = np.bincount(d.meta, minlength=len(d.meta_t0))
    assert np.all(counts >= 1) and np.all(counts <= d.meta_n)
    order = np.lexsort((d.rank, d.meta))
    m, r, t = d.meta[order], d.rank[order], d.t[order]
    first = np.r_[True, m[1:] != m[:-1]]
    assert np.all(r[first] == 1) and np.all(np.diff(r)[~first[1:]] == 1)
    assert np.allclose(t[first], d.meta_t0[m[first]])  # the first child executes at the start
    gaps = np.diff(t)[~first[1:]]
    assert abs(gaps.mean() * P2.phi - 1) < 4 / np.sqrt(len(gaps))  # Exp(phi): mean 1/phi
    assert abs(gaps.std() / gaps.mean() - 1) < 0.05


def check_mid_direct_sum_trade_clock():
    """mid_k = sum_{j<k} eps_j theta sqrt(q_j) (phi (t_j - t0_j) + n0)^(beta_j - 1/2) (1/(1+k-j))^beta_j,
    written from II Algorithm 1 line 16 with t_k - t_j = (k - j) tau0, from the day's own fields."""
    d = M.simulate_day(P2, 600, np.random.default_rng(6))
    b = d.meta_beta[d.meta]
    amp = d.eps * P2.theta * np.sqrt(d.q) * (P2.phi * (d.t - d.meta_t0[d.meta]) + P2.n0) ** (b - 0.5)
    want = np.array([np.sum(amp[:k] * (1.0 / (1.0 + k - np.arange(k))) ** b[:k]) for k in range(len(d.t))])
    assert len(d.t) > 1000 and want[0] == 0.0
    assert np.max(np.abs(d.impact - want)) < 1e-9 * np.max(np.abs(want))
    assert np.array_equal(d.mid, d.impact)  # sigma_F = 0


def check_mid_direct_sum_physical_clock():
    p = P2
    d = M.simulate_day(p, 120, np.random.default_rng(7), clock="physical")
    n = len(d.t)
    tau0 = (d.t[-1] - d.t[0]) / (n - 1)
    b = d.meta_beta[d.meta]
    want = np.zeros(n)
    for k in range(n):
        for j in range(k):
            a = d.eps[j] * np.sqrt(d.q[j]) * (p.phi * (d.t[j] - d.meta_t0[d.meta[j]]) + p.n0) ** (b[j] - 0.5)
            want[k] += a * (tau0 / (d.t[k] - d.t[j] + tau0)) ** b[j]
    assert n > 200 and np.isclose(d.tau0, tau0)
    assert np.max(np.abs(d.impact - want)) < 1e-12 * np.max(np.abs(want))


def check_clocks_agree():
    """Trade clock against II's physical clock on the same days (same seeds, same flow).
    Measured on 3 seeds (about 4000 trades): correlation of 10-trade mid changes 0.978 to
    0.986, of 1-trade changes 0.948 to 0.954 (the trade clock's 1-trade variance is 14 to
    17% lower); variance ratio of 100-trade changes 0.93 to 1.01. Tolerances: 0.96, 0.8 to 1.2."""
    dp = M.simulate_day(P2, 1500, np.random.default_rng(8), clock="physical")
    dt = M.simulate_day(P2, 1500, np.random.default_rng(8), clock="trade")
    assert np.array_equal(dp.t, dt.t) and np.array_equal(dp.eps, dt.eps)
    a, c = np.diff(dp.mid[::10]), np.diff(dt.mid[::10])
    assert np.corrcoef(a, c)[0, 1] > 0.96
    a, c = np.diff(dp.mid[::100]), np.diff(dt.mid[::100])
    assert 0.8 < np.var(c) / np.var(a) < 1.2


def check_metaorder_sign_correlation():
    """Measured over 20 seeds of 2^16 metaorders (Gamma = 0.1, gamma_x = 0.6): lag-1 ACF
    0.1161 (exact 0.1161, sd 0.0038), mean over lags 10..50 0.0136 (exact 0.0143, sd
    0.0018). Tolerances 4 sd. The exact form tends to Gamma k^-gamma_x (II (1))."""
    e = M.metaorder_signs(1 << 16, 0.1, 0.6, np.random.default_rng(9)).astype(float)
    C_ = P.sign_acf(e, 50)
    assert abs(C_[1] - M.sign_correlation(np.array([1]), 0.1, 0.6)[0]) < 0.015
    assert abs(C_[10:51].mean() - M.sign_correlation(np.arange(10, 51), 0.1, 0.6).mean()) < 0.0072
    k = np.array([100.0, 1000.0])
    assert np.allclose(M.sign_correlation(k, 0.1, 0.6), 0.1 * k**-0.6, rtol=0.01)
    e0 = M.metaorder_signs(1 << 16, 0.0, 0.6, np.random.default_rng(10)).astype(float)
    assert np.all(np.abs(P.sign_acf(e0, 20)[1:]) < 4 / np.sqrt(len(e0)))
    with pytest.raises(ValueError):
        M.metaorder_signs(100, 0.2, 0.6, np.random.default_rng(0))  # Gamma above (2/pi) H (2H - 1) = 0.178


def check_informed_component():
    """I Sec. 6.2: with theta = 0 the mid is the fundamental part alone. Exact: the move of
    metaorder i lands after its first child. Statistical: per metaorder corr(eps_i, dF_i) = rho,
    and over windows of T trades corr(sum eps over metaorders starting there, change of F) =
    rho at every T, i.e. E[eps xi] = rho / sqrt(nu T). Gamma = 0, rho = 0.3; tolerances 4 SE."""
    p = replace(P2, theta=0.0, Gamma=0.0, sigma_F=1.0, rho=0.3)
    d = M.simulate_day(p, 20000, np.random.default_rng(11))
    first = np.flatnonzero(d.rank == 1)  # increasing
    cum = np.concatenate([[0.0], np.cumsum(d.meta_dF[d.meta[first]])])
    F = cum[np.searchsorted(first, np.arange(len(d.t)), side="left")]  # metaorders whose first child precedes k
    assert np.allclose(d.mid, F) and np.all(d.impact == 0)
    assert np.isclose(np.std(d.meta_dF), 1 / np.sqrt(p.nu), rtol=0.03)
    r = np.corrcoef(d.meta_eps, d.meta_dF)[0, 1]
    assert abs(r - 0.3) < 4 * (1 - 0.09) / np.sqrt(len(d.meta_dF))
    days = M.simulate(p, 6, 20000, 12)
    for T in (100, 1000):
        se_, dF_ = [], []
        for dd in days:
            f = np.flatnonzero(dd.rank == 1)
            w = f // T
            nw = len(dd.t) // T
            ok = w < nw
            se_.append(np.bincount(w[ok], dd.meta_eps[dd.meta[f[ok]]], minlength=nw))
            dF_.append(np.bincount(w[ok], dd.meta_dF[dd.meta[f[ok]]], minlength=nw))
        se_, dF_ = np.concatenate(se_), np.concatenate(dF_)
        r = np.corrcoef(se_, dF_)[0, 1]
        assert abs(r - 0.3) < 4 * (1 - 0.09) / np.sqrt(len(se_)), (T, r)


def check_input_validation():
    with pytest.raises(ValueError):
        M.simulate_day(replace(P2, mu_m=1.1, lam=0.1), 200, np.random.default_rng(0))  # mu_q <= 1 below log q = 2
    with pytest.raises(ValueError):
        M.simulate_day(P2, 200, np.random.default_rng(0), clock="wall")


# ---- Paper II's qualitative results -------------------------------------------------


def check_correlation_peaks_near_half():
    """R_a(T) = K / sqrt(S V) at T = 10, a = 0..1.5 step 0.05, 30 days of 10^4 metaorders.
    With lam = lam' = 0 the peak is at a = 1/2 exactly in expectation: impact ~ sqrt(q)
    and R ~ E[q^(a+1/2)] / sqrt(E[q^2a]), largest at a = 1/2 for lognormal q; measured 0.50
    on 3 of 3 seeds. At C-VD-VF the q-dependent exponents move it to 0.65 (3 of 3 seeds;
    Paper II Fig. 7: a = 0.5 to 1). Required: 0.45..0.55 and 0.5..0.8."""
    a = np.round(np.arange(0, 1.51, 0.05), 2)
    for p, lo, hi, seed in ((M.configuration("NC-NVD-VF"), 0.45, 0.55, 13), (P2, 0.5, 0.8, 14)):
        V, K, S = pooled_moments(M.simulate(p, 30, 10000, seed), [1, 10], a)
        R = K[:, 1] / np.sqrt(S[1] * V[:, 1])
        assert lo <= a[np.argmax(R)] <= hi, (p, a[np.argmax(R)])


def check_sign_imbalance_collapse():
    """chi_hat of I^0_T (collapse.py) on T = 512..8192, 16 days of 4 x 10^4 metaorders
    (about 1.7M trades). Measured on 4 seeds: C-NVD-VF (Paper II Fig. 5's configuration)
    0.6925 sd 0.009, NC-NVD-VF 0.661 sd 0.013, against 1/mu = 0.667. Required: within 0.04.
    Short horizons sit higher, 0.745 to 0.78 on T = 32..1024, which is Paper II's 0.75;
    at the full C-VD-VF preset chi is 0.737 sd 0.027 on T = 512..8192, as small-q
    metaorders have heavier duration tails (I (29) at a = 0: mu_eff = mu_m - lam^2 sigma^2 log T / 2)."""
    p_grid = np.round(np.arange(0.05, 0.951, 0.05), 3)
    T = 2 ** np.arange(9, 14)
    for name, seed in (("C-NVD-VF", 15), ("NC-NVD-VF", 16)):
        series = [d.eps.astype(float) for d in M.simulate(M.configuration(name), 16, 40000, seed)]
        Q, dQ, iqr = C.quantile_curves(series, T, p_grid)
        chi = C.chi_hat(Q, dQ, iqr, T, p_grid)[0]
        assert abs(chi - 1 / 1.5) < 0.04, (name, chi)


def check_covariance_exponent_at_paper2_settings():
    """Paper II Fig. 6 reports a V-shaped exponent of K_a(T) = E[Delta_T I^a_T] in a for
    C-VD-VF. It is NOT reproduced at the stated settings: fitted on T = 10..1000 over
    40 days of 10^4 metaorders, the exponent falls at every step of 0.5 in a over 0..2 on 6
    of 6 seeds (1.33 to 1.27; smallest step 0.0046) and keeps falling to 1.19 at a = 3.
    I (66)'s rising branch is T^(1 - beta_hat(a)) with beta_hat(a) = beta_m - (a + 1/2) lam'
    sigma_l^2, which needs beta_q < 0; II's 0 < beta_q < 1 clips it (16% of metaorders sit at
    beta = 0), capping that branch near T^1, below the off-diagonal T^(2 - gamma_x - beta_hat(0))
    ~ T^1.3. Unclipped beta (impact growing after execution) gives a rising curve from a ~ 0.25
    but a K_a(T) that is unstable across seeds (sd up to 0.19). This check pins the measured
    shape so a change to it is noticed."""
    a = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
    zV, zK = exponents_vs_a(M.simulate(P2, 40, 10000, 17), a)
    assert np.all(np.diff(zK) < 0) and zK[0] - zK[-1] > 0.03, zK
    assert 1.25 < zK[0] < 1.4  # off-diagonal dominated: 2 - gamma_x - beta_hat(0) = 1.275
    assert np.all(np.diff(zV) < 0)  # II Fig. 2: the imbalance exponent falls with a


def check_independent_volume_null():
    """With lam = lam' = 0 and independent metaorder signs (NC-NVD-VF), q is constant within
    a metaorder and independent of everything else, so E[V_a(T)] = E[q^2a] E[V_0(T)] and
    E[K_a(T)] = E[q^(a+1/2)] E[K_0(T)] / E[q^1/2]: the a-dependence of both exponents vanishes.
    Measured over 4 seeds of 80 days x 10^4 metaorders: zV(0.5) - zV(0) = 0.004 sd 0.017,
    zK(0.5) - zK(0) = 0.002 sd 0.005, zK(1) - zK(0) = 0.006 sd 0.014; with lam != 0
    (NC-VD-VF) zV(0.5) - zV(0) = -0.074 sd 0.009. Required: |.| < 0.05, 0.02, 0.05 under
    the null, < -0.04 with lam. At larger a, rare large-q metaorders dominate q^2a and the
    sample exponents scatter (sd 0.11 at a = 1.5). With Gamma > 0 the null does not hold
    exactly: pairs from different metaorders carry E[q^a]^2 / E[q^2a] = exp(-sigma^2 a^2)."""
    a = np.array([0.0, 0.5, 1.0])
    zV, zK = exponents_vs_a(M.simulate(M.configuration("NC-NVD-VF"), 80, 10000, 18), a)
    assert abs(zV[1] - zV[0]) < 0.05 and abs(zK[1] - zK[0]) < 0.02 and abs(zK[2] - zK[0]) < 0.05, (zV, zK)
    zV, _ = exponents_vs_a(M.simulate(M.configuration("NC-VD-VF"), 80, 10000, 19), a)
    assert zV[1] - zV[0] < -0.04, zV


CHECKS = {f.__name__: f for f in [
    check_parameters_and_configurations, check_metaorder_table, check_child_counts, check_child_schedule,
    check_mid_direct_sum_trade_clock, check_mid_direct_sum_physical_clock, check_clocks_agree,
    check_metaorder_sign_correlation, check_informed_component, check_input_validation,
    check_correlation_peaks_near_half, check_sign_imbalance_collapse, check_covariance_exponent_at_paper2_settings,
    check_independent_volume_null,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_mbmodel(name):
    CHECKS[name]()


# ---- planted bugs -----------------------------------------------------------------


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("beta_q = beta_m + lam' (log q - m)", "exponents",
        "p.beta_m - p.lam_p * (ell - p.m)", "p.beta_m + p.lam_p * (ell - p.m)",
        ("check_metaorder_table",)),
    Bug("mu_q = mu_m - lam (log q - m)", "exponents",
        "mu = p.mu_m + p.lam * (ell - p.m)", "mu = p.mu_m - p.lam * (ell - p.m)",
        ("check_metaorder_table",)),
    Bug("q not truncated below 1", "child_sizes",
        "while np.any(bad):", "while False:",
        ("check_metaorder_table",)),
    Bug("child count rounded up", "simulate_day",
        "n_draw = np.floor(rng.random(n_m) ** (-1.0 / mu))", "n_draw = np.ceil(rng.random(n_m) ** (-1.0 / mu))",
        ("check_child_counts",)),
    Bug("child gaps with mean phi instead of 1/phi", "simulate_day",
        "gaps = rng.exponential(1 / p.phi, int(n_gen.sum()))", "gaps = rng.exponential(p.phi, int(n_gen.sum()))",
        ("check_child_schedule",)),
    Bug("mid includes the trade's own impact", "trade_kernel",
        "g[0] = 0.0", "g[0] = 1.0",
        ("check_mid_direct_sum_trade_clock",)),
    Bug("kernel without the tau0 offset", "trade_kernel",
        "g = (1.0 + np.arange(n, dtype=float)) ** -beta", "g = np.maximum(np.arange(n, dtype=float), 1.0) ** -beta",
        ("check_mid_direct_sum_trade_clock",)),
    Bug("amplitude exponent 1/2 - beta", "simulate_day",
        "(p.phi * tau + p.n0) ** (beta[owner] - 0.5)", "(p.phi * tau + p.n0) ** (0.5 - beta[owner])",
        ("check_mid_direct_sum_trade_clock",)),
    Bug("two Chebyshev nodes", "impact_trade_clock",
        "n_nodes: int = 20", "n_nodes: int = 2",
        ("check_mid_direct_sum_trade_clock",)),
    Bug("sign mixing probability Gamma/Gmax, not its square root", "metaorder_signs",
        "p = np.sqrt(Gamma / gmax)", "p = Gamma / gmax",
        ("check_metaorder_sign_correlation",)),
    Bug("informed sign flipped", "simulate_day",
        "p.rho * sgn", "-p.rho * sgn",
        ("check_informed_component",)),
    Bug("fundamental move before the first child", "simulate_day",
        "np.add.at(F, first + 1, dF[owner[first]])", "np.add.at(F, first, dF[owner[first]])",
        ("check_informed_component",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(M, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
