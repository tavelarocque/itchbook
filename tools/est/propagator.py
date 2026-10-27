"""Trade-level propagator kernel G from signs and mid returns.

Model (transient impact): r_n = m_{n+1} - m_n = sum_{j>=0} K(j) eps_{n-j} + noise,
with K(j) = G(j+1) - G(j), G(0) = 0. Then, for l = 0..L-1,

    S(l) = E[eps_n r_{n+l}] = sum_{j=0}^{L-1} K(j) C(|l - j|)        (truncated at L)

with C(k) = E[eps_n eps_{n+k}] (uncentered, within day). Solving this Toeplitz
system gives K, and G(l) = sum_{j<l} K(j) for l = 1..L. fit_power fits
G(l) = G0 l0^beta / (l0 + l)^beta on l in [l_min, l_max] by weighted least
squares over a (beta, l0) grid with G0 in closed form; log-spaced lags get
equal weight.
"""
from __future__ import annotations

import numpy as np


def sign_acf(eps: np.ndarray, L: int) -> np.ndarray:
    """Uncentered C(k), k = 0..L."""
    e = np.asarray(eps, dtype=float)
    n = len(e)
    size = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(e, size)
    ac = np.fft.irfft(f * np.conj(f), size)[: L + 1]
    return ac / (n - np.arange(L + 1))


def response_increments(eps: np.ndarray, m: np.ndarray, L: int, ok: np.ndarray | None = None) -> np.ndarray:
    """S(l) = mean_n eps_n (m_{n+l+1} - m_{n+l}), l = 0..L-1; pairs with an invalid mid are skipped."""
    e = np.asarray(eps, dtype=float)
    m = np.asarray(m, dtype=float)
    n = len(e)
    r = np.diff(m)  # r_n = m_{n+1} - m_n, n = 0..N-2
    rok = np.ones(n - 1, bool) if ok is None else (ok[:-1] & ok[1:])
    out = np.empty(L)
    for l in range(L):
        k = n - 1 - l
        sel = rok[l : l + k]
        prod = e[:k] * r[l : l + k]
        out[l] = prod[sel].mean()
    return out


def calibrate_G(C: np.ndarray, S: np.ndarray) -> np.ndarray:
    """G(l), l = 0..L, from C(0..L-1) and S(0..L-1)."""
    L = len(S)
    idx = np.abs(np.arange(L)[:, None] - np.arange(L)[None, :])
    K = np.linalg.solve(C[idx], S)
    return np.concatenate([[0.0], np.cumsum(K)])


def fit_power(G: np.ndarray, l_min: int, l_max: int, betas=None, l0s=None):
    """Least-squares fit of G(l) = G0 l0^beta/(l0+l)^beta on log-spaced lags; returns (G0, l0, beta, rss)."""
    lags = np.unique(np.round(np.logspace(np.log10(l_min), np.log10(l_max), 40)).astype(int))
    y = G[lags]
    betas = np.linspace(0.0, 1.2, 241) if betas is None else betas
    l0s = np.logspace(-1, 2.5, 71) if l0s is None else l0s
    best = (np.nan, np.nan, np.nan, np.inf)
    for b in betas:
        for l0 in l0s:
            f = (l0 / (l0 + lags)) ** b
            G0 = float(f @ y / (f @ f))
            rss = float(np.sum((y - G0 * f) ** 2))
            if rss < best[3]:
                best = (G0, l0, b, rss)
    return best
