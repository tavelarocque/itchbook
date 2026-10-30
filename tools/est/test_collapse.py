"""Tests for collapse.py: the local fit and S against slow references written from
the Houdayer-Hartmann text, S on curves with known errors, the quantile errors,
then recovery of chi on three processes with known exponents, then planted bugs.

    uv run --group dev --group analysis pytest tools/est/test_collapse.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import collapse as C
import propagator as P
import synth as Y

P_GRID = np.round(np.arange(0.05, 0.951, 0.05), 3)


def stable(alpha: float, n: int, rng: np.random.Generator) -> np.ndarray:
    """Symmetric alpha-stable variables (Chambers-Mallows-Stuck): sums of T of them are
    exactly T^(1/alpha) times one of them, so I_T collapses exactly at chi = 1/alpha."""
    V = rng.uniform(-np.pi / 2, np.pi / 2, n)
    W = rng.exponential(1.0, n)
    return np.sin(alpha * V) / np.cos(V) ** (1 / alpha) * (np.cos((1 - alpha) * V) / W) ** ((1 - alpha) / alpha)


def lmf_one(n: int, alpha: float, rng: np.random.Generator) -> np.ndarray:
    """Lillo-Mike-Farmer signs with one live metaorder: back-to-back runs of one random
    sign with P(L >= l) = l^-alpha. synth.lmf_signs(n, alpha, 1) is this same process
    one step at a time; repeating whole runs makes 32M steps affordable."""
    L = np.floor(rng.random(n // 2 + 10) ** (-1 / alpha)).astype(np.int64)
    k = np.searchsorted(np.cumsum(L), n) + 1
    return np.repeat(rng.choice(np.array([-1, 1], dtype=np.int8), k), L[:k])[:n]


def chi_of(series, T, p=P_GRID) -> float:
    Q, dQ, iqr = C.quantile_curves(series, T, p)
    return C.chi_hat(Q, dQ, iqr, T, p)[0]


def slow_quality(xs, ys, dys):
    """S from the paper's recipe with plain loops and a matrix least-squares fit."""
    terms = []
    for i in range(len(xs)):
        for j in range(len(xs[i])):
            x0, pts = xs[i][j], []
            for k in range(len(xs)):
                xk = xs[k]
                if k == i or not xk[0] <= x0 <= xk[-1]:
                    continue
                jj = len(xk) - 2 if x0 == xk[-1] else max(m for m in range(len(xk) - 1) if xk[m] <= x0)
                pts += [(xk[jj], ys[k][jj], dys[k][jj]), (xk[jj + 1], ys[k][jj + 1], dys[k][jj + 1])]
            if not pts:
                continue
            X = np.array([[1.0, px] for px, _, _ in pts])
            W = np.diag([1 / d**2 for _, _, d in pts])
            M = X.T @ W @ X
            coef = np.linalg.solve(M, X.T @ W @ np.array([py for _, py, _ in pts]))
            u = np.array([1.0, x0])
            Yfit, dY2 = coef @ u, u @ np.linalg.inv(M) @ u
            terms.append((ys[i][j] - Yfit) ** 2 / (dys[i][j] ** 2 + dY2))
    return float(np.mean(terms)), len(terms)


def known_error_curves(seed: int, chi: float = 0.7):
    """Three curves T^chi f(x) + noise with known dy on one grid of 201 x: at the right
    chi they agree within their errors, so E[S] = 1 exactly (the local line through
    the other curves' nodes at x_j and x_j+1 is unbiased at x_j)."""
    rng = np.random.default_rng(seed)
    T, x = np.array([64.0, 256.0, 1024.0]), np.linspace(0.05, 0.95, 201)
    dy = (0.02 + 0.03 * x)[None, :] * T[:, None] ** chi
    y = T[:, None] ** chi * np.tan(np.pi * (x - 0.5))[None, :] + dy * rng.standard_normal((3, 201))
    return T, x, y, dy


# ---- checks -----------------------------------------------------------------------


def check_local_fit_matches_least_squares():
    rng = np.random.default_rng(1)
    x, dy = rng.uniform(0, 3, 9), rng.uniform(0.1, 0.5, 9)
    y = 1.0 - 2.0 * x + dy * rng.standard_normal(9)
    X, W = np.column_stack([np.ones(9), x]), np.diag(1 / dy**2)
    M = X.T @ W @ X
    coef = np.linalg.solve(M, X.T @ W @ y)
    for x0 in (0.0, 1.3, 4.0):
        Yfit, dY = C.local_fit(x, y, dy, x0)
        u = np.array([1.0, x0])
        assert np.isclose(Yfit, coef @ u) and np.isclose(dY**2, u @ np.linalg.inv(M) @ u)


def check_quality_matches_slow_reference():
    rng = np.random.default_rng(2)
    for trial in range(6):
        n_c = rng.integers(2, 5)
        xs = [np.sort(rng.uniform(rng.uniform(0, 1), rng.uniform(1.5, 3), rng.integers(3, 12))) for _ in range(n_c)]
        if trial == 0:
            xs = [np.linspace(0, 1, 7)] * n_c  # coincident nodes: the tie rule matters
        ys = [np.sin(x) + 0.1 * rng.standard_normal(len(x)) for x in xs]
        dys = [rng.uniform(0.05, 0.2, len(x)) for x in xs]
        S, N = C.quality(xs, ys, dys)
        S_ref, N_ref = slow_quality(xs, ys, dys)
        assert N == N_ref and np.isclose(S, S_ref, rtol=1e-10), (trial, S, S_ref)


def check_quality_invariances():
    T, x, y, dy = known_error_curves(3)
    S = C.quality([x] * 3, list(y), list(dy))[0]
    assert np.isclose(C.quality([x] * 3, list(7.5 * y), list(7.5 * dy))[0], S)  # common scale of y and dy
    assert np.isclose(C.quality([3 * x - 1] * 3, list(y), list(dy))[0], S)  # affine map of x


def check_S_near_one_for_a_perfect_collapse():
    """Measured over 200 seeds (N = 603 points): S = 0.992, sd 0.070, range 0.85 to
    1.19; tolerance 0.25. A wrong chi by 0.1 gave S >= 36 on 300 seeds of a smaller
    design; required here: S > 10."""
    T, x, y, dy = known_error_curves(4)
    S = lambda chi: C.quality([x] * 3, list(y / T[:, None] ** chi), list(dy / T[:, None] ** chi))[0]  # noqa: E731
    assert abs(S(0.7) - 1) < 0.25, S(0.7)
    assert S(0.6) > 10 and S(0.8) > 10


def check_quantile_errors_are_calibrated():
    """dQ against the spread of sample quantiles over 400 iid normal samples of 4000:
    measured ratio 0.966 to 1.084 over p = 0.05..0.95; tolerance 15%."""
    rng = np.random.default_rng(5)
    Qs, dQs = [], []
    for _ in range(400):
        Q, dQ, _ = C.quantile_curves([rng.standard_normal(4000)], [1], P_GRID)
        Qs.append(Q[0])
        dQs.append(dQ[0])
    ratio = np.mean(dQs, axis=0) / np.std(Qs, axis=0)
    assert np.all(np.abs(ratio - 1) < 0.15), ratio


def check_window_sums_and_lattice_errors():
    x = np.arange(10.0)
    assert np.array_equal(C.window_sums(x, 4), [0 + 1 + 2 + 3, 1 + 2 + 3 + 4, 2 + 3 + 4 + 5, 3 + 4 + 5 + 6,
                                                4 + 5 + 6 + 7, 5 + 6 + 7 + 8, 6 + 7 + 8 + 9])
    assert np.array_equal(C.window_sums(x, 8), [28, 44])  # stride 2: starts 0 and 2
    e = np.random.default_rng(6).choice([-1.0, 1.0], 4096)  # sums of 16 signs live on even integers
    Q, dQ, iqr = C.quantile_curves([e], [16], P_GRID)
    assert np.all(np.isfinite(dQ)) and np.all(dQ > 0)


def check_chi_on_stable_sums():
    """Exact self-similarity, chi = 1/alpha. Measured on 3 seeds (2^20 variables, T = 4..512):
    within 0.004 of 1/alpha for alpha = 1.5 and 1.2; required: within 0.03."""
    T = 2 ** np.arange(2, 10)
    for alpha, seed in ((1.5, 7), (1.2, 8)):
        chi = chi_of([stable(alpha, 1 << 20, np.random.default_rng(seed))], T)
        assert abs(chi - 1 / alpha) < 0.03, (alpha, chi)


def check_chi_on_fgn_signs():
    """sign(fGn), H = 0.75: Gaussian-like aggregation, chi = H. The variance is
    T + c T^(2H) at finite T, so the exact effective exponent on T = 256..8192 is
    0.7456. Measured on 4 seeds (4 series of 2^20): 0.754, sd 0.005. Shorter series or
    T from 64 drift low (0.729 +- 0.008 on 8 x 2^18, T = 64..2048), from the short-range
    term and the slow convergence of long-memory samples. Required: within 0.03 of H."""
    rng = np.random.default_rng(9)
    chi = chi_of([Y.fgn_signs(1 << 20, 0.75, rng).astype(float) for _ in range(4)], 2 ** np.arange(8, 14))
    assert abs(chi - 0.75) < 0.03, chi


def check_lmf_generator():
    """lmf_one has the exact stationary renewal autocorrelation
    C(k) = sum_{j > k} j^-alpha / zeta(alpha) = 0.6172, 0.2362 at k = 1, 10, and so does
    synth.lmf_signs with one live order. Measured over 30 seeds: lmf_one (2^20 steps)
    0.6162, 0.2338 with sd 0.0061, 0.0112; synth (2^18) 0.6151, 0.2334 with sd 0.0103,
    0.0199. Tolerances 4 sd."""
    j = np.arange(1, 2_000_001, dtype=float) ** -1.5
    zeta = j.sum() + 2 / np.sqrt(2_000_000.5)  # tail of the sum by the integral
    want = np.array([(zeta - j[:k].sum()) / zeta for k in (1, 10)])
    a = P.sign_acf(lmf_one(1 << 20, 1.5, np.random.default_rng(10)).astype(float), 10)[[1, 10]]
    b = P.sign_acf(Y.lmf_signs(1 << 18, 1.5, 1, np.random.default_rng(11)).astype(float), 10)[[1, 10]]
    assert np.all(np.abs(a - want) < [0.025, 0.045]) and np.all(np.abs(b - want) < [0.045, 0.08]), (a, b, want)


def check_chi_on_lmf_signs():
    """LMF, alpha = 1.5: Levy-like aggregation, chi = 1/alpha = 0.667, not the standard
    deviation exponent (3 - alpha)/2 = 0.75. Convergence is slow and from below.
    Measured on 4 seeds (32 series of 2^20, one live order): 0.651 (sd 0.005) on
    T = 1024..16384, 0.656 (sd 0.006) on 2048..32768, 0.662 (sd 0.010) on 4096..65536.
    With 5 live orders the |I| <= T/N truncation pulls chi up (0.72 to 0.76 on
    T = 256..8192). Required on T = 2048..32768: within 0.03 of 1/alpha."""
    rng = np.random.default_rng(12)
    chi = chi_of([lmf_one(1 << 20, 1.5, rng) for _ in range(32)], 2 ** np.arange(11, 16))
    assert abs(chi - 1 / 1.5) < 0.03, chi


CHECKS = {f.__name__: f for f in [
    check_local_fit_matches_least_squares, check_quality_matches_slow_reference, check_quality_invariances,
    check_S_near_one_for_a_perfect_collapse, check_quantile_errors_are_calibrated,
    check_window_sums_and_lattice_errors, check_chi_on_stable_sums, check_chi_on_fgn_signs, check_lmf_generator,
    check_chi_on_lmf_signs,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_collapse(name):
    CHECKS[name]()


# ---- planted bugs -----------------------------------------------------------------


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("local fit uses the curve's own points", "quality",
        "if k == i:\n                continue", "if False:\n                continue",
        ("check_quality_matches_slow_reference", "check_S_near_one_for_a_perfect_collapse")),
    Bug("dY^2 left out of the denominator", "quality",
        "(1.0 / ws[i][sel] + dY2)", "(1.0 / ws[i][sel])",
        ("check_quality_matches_slow_reference", "check_S_near_one_for_a_perfect_collapse")),
    Bug("dY^2 with +2 x Kx", "quality",
        "dY2 = (Kxx - 2 * xi * Kx + xi * xi * K) / D", "dY2 = (Kxx + 2 * xi * Kx + xi * xi * K) / D",
        ("check_quality_matches_slow_reference",)),
    Bug("bracketing pair shifted right", "quality",
        'np.searchsorted(xk, x0, side="right") - 1', 'np.searchsorted(xk, x0, side="right")',
        ("check_quality_matches_slow_reference",)),
    Bug("each curve divided by its own IQR (chi cancels)", "rescale",
        "c = np.median(np.asarray(iqr)[:, None] / s)", "c = np.asarray(iqr)[:, None] / s",
        ("check_chi_on_stable_sums",)),
    Bug("rescaled by T^-chi instead of T^chi", "rescale",
        "s = np.asarray(T_grid, dtype=float)[:, None] ** chi", "s = np.asarray(T_grid, dtype=float)[:, None] ** -chi",
        ("check_chi_on_stable_sums",)),
    Bug("quantile error without the square root", "quantile_curves",
        "np.sqrt(p * (1 - p) / n_eff)", "(p * (1 - p) / n_eff)",
        ("check_quantile_errors_are_calibrated",)),
    Bug("windows of T - 1 terms", "window_sums",
        "return c[n + T] - c[n]", "return c[n + T - 1] - c[n]",
        ("check_window_sums_and_lattice_errors",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(C, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
