"""Synthetic order flow and prices with known exponents, for estimator validation.

Signs
  fgn(n, H)              fractional Gaussian noise (Davies-Harte, exact)
  fgn_signs(n, H)        sign(fGn): C(k) = (2/pi) arcsin rho_H(k) ~ k^-gamma, gamma = 2 - 2H;
                         Gaussian-like aggregation, so the sign-imbalance width scales as T^H
  lmf_signs(n, alpha, N) Lillo-Mike-Farmer splitting: N live metaorders with Pareto(alpha)
                         sizes, one picked uniformly per step; C(k) ~ k^-gamma, gamma = alpha - 1;
                         Levy-like aggregation, width ~ T^(1/alpha) (the paper's chi = 1/mu)
Sizes
  lognormal_sizes(n, m, sigma)   q = exp(m + sigma Z), independent of signs
Prices
  propagator_mid(x, G)   m_n = sum_{k<n} G(n-k) x_k: the mid just before event n,
                         x_k = eps_k q_k^psi, G(l) = G0 l0^beta / (l0 + l)^beta
"""
from __future__ import annotations

import numpy as np


def fgn_acov(n: int, H: float) -> np.ndarray:
    k = np.arange(n, dtype=float)
    return 0.5 * ((k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))


def fgn(n: int, H: float, rng: np.random.Generator) -> np.ndarray:
    """Exact fGn with unit variance by circulant embedding (Davies-Harte)."""
    if not 0 < H < 1:
        raise ValueError("H must be in (0, 1)")
    g = fgn_acov(n + 1, H)
    c = np.concatenate([g, g[-2:0:-1]])  # first row of the 2n circulant
    m = len(c)
    lam = np.fft.fft(c).real
    if lam.min() < -1e-9 * lam.max():
        raise RuntimeError("circulant embedding not non-negative definite")
    lam = np.clip(lam, 0.0, None)
    z = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    return np.fft.fft(np.sqrt(lam / m) * z).real[:n]


def fgn_signs(n: int, H: float, rng: np.random.Generator) -> np.ndarray:
    return np.where(fgn(n, H, rng) >= 0, 1, -1).astype(np.int8)


def fgn_sign_acf(k: np.ndarray, H: float) -> np.ndarray:
    """Exact autocorrelation of sign(fGn) at lags k >= 0 (arcsine law)."""
    k = np.asarray(k, dtype=float)
    rho = 0.5 * ((k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))
    return (2 / np.pi) * np.arcsin(np.clip(rho, -1, 1))


def lmf_signs(n: int, alpha: float, n_orders: int, rng: np.random.Generator, burn: int | None = None) -> np.ndarray:
    """Lillo-Mike-Farmer: sizes P(L >= l) = l^-alpha (l >= 1, integer), signs +-1 at random."""
    if not 1 < alpha < 2:
        raise ValueError("alpha must be in (1, 2) for long memory with gamma = alpha - 1")
    burn = 20 * n_orders if burn is None else burn
    total = n + burn

    def new_size() -> int:
        return int(np.floor(rng.random() ** (-1.0 / alpha)))

    remaining = np.array([new_size() for _ in range(n_orders)], dtype=np.int64)
    sign = rng.choice(np.array([-1, 1], dtype=np.int8), size=n_orders)
    picks = rng.integers(0, n_orders, size=total)
    flips = rng.random(total)
    u = rng.random(total)
    out = np.empty(total, dtype=np.int8)
    for t in range(total):
        i = picks[t]
        out[t] = sign[i]
        remaining[i] -= 1
        if remaining[i] <= 0:
            remaining[i] = int(np.floor(u[t] ** (-1.0 / alpha)))
            sign[i] = 1 if flips[t] < 0.5 else -1
    return out[burn:]


def lognormal_sizes(n: int, m: float, sigma: float, rng: np.random.Generator) -> np.ndarray:
    return np.exp(m + sigma * rng.standard_normal(n))


def power_kernel(lmax: int, G0: float, l0: float, beta: float) -> np.ndarray:
    """G(l) for l = 0..lmax, with G(0) = 0 (an event does not move the mid before it)."""
    l = np.arange(lmax + 1, dtype=float)
    g = G0 * l0**beta / (l0 + l) ** beta
    g[0] = 0.0
    return g


def propagator_mid(x: np.ndarray, G: np.ndarray, noise: float = 0.0, rng: np.random.Generator | None = None):
    """m_n = sum_{k<n} G(n-k) x_k (+ an optional random walk of step `noise`), n = 0..len(x)-1.

    G is given for lags 0..len(x)-1 or longer; G(0) is ignored."""
    n = len(x)
    g = np.zeros(n)
    g[1:] = G[1:n]
    size = 1 << int(np.ceil(np.log2(2 * n)))
    m = np.fft.irfft(np.fft.rfft(x, size) * np.fft.rfft(g, size), size)[:n]
    if noise > 0:
        if rng is None:
            raise ValueError("noise needs rng")
        m = m + noise * np.concatenate([[0.0], np.cumsum(rng.standard_normal(n - 1))])
    return m
