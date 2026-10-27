"""Long-memory exponent of order signs.

local_whittle(x): Robinson (1995) local Whittle estimate of the memory
parameter d from the first m Fourier frequencies,

    R(d) = log( (1/m) sum_j lambda_j^{2d} I(lambda_j) ) - (2d/m) sum_j log lambda_j,
    lambda_j = 2 pi j / n,  I(l) = |sum_t x_t e^{i t l}|^2 / (2 pi n),  j = 1..m,

minimised over d; sqrt(m)(d_hat - d) -> N(0, 1/4). For signs with C(k) ~ k^-gamma,
d = (1 - gamma) / 2, so gamma_hat = 1 - 2 d_hat with SE = 1/sqrt(m). Only j >= 1
enters, so the day's mean (and any constant) does not matter.

pool_days(gamma, se): inverse-variance pool with Cochran's Q.
"""
from __future__ import annotations

import numpy as np


def _golden(f, lo: float, hi: float, tol: float = 1e-6) -> float:
    g = (np.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    while b - a > tol:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return (a + b) / 2


def local_whittle(x: np.ndarray, power: float = 0.65, m: int | None = None, bounds=(-0.49, 0.99)):
    """Returns (d_hat, se_d, m)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    m = int(n**power) if m is None else m
    if m < 10 or m >= n // 2:
        raise ValueError(f"bandwidth m={m} unusable for n={n}")
    w = np.fft.rfft(x)[1 : m + 1]
    I = (w.real**2 + w.imag**2) / (2 * np.pi * n)
    lam = 2 * np.pi * np.arange(1, m + 1) / n
    loglam = np.log(lam)
    mean_loglam = loglam.mean()

    def R(d: float) -> float:
        return np.log(np.mean(np.exp(2 * d * loglam) * I)) - 2 * d * mean_loglam

    d = _golden(R, *bounds)
    return d, 1 / (2 * np.sqrt(m)), m


def gamma_lw(eps: np.ndarray, power: float = 0.65):
    """gamma_hat = 1 - 2 d_hat and its SE, from one day's sign series."""
    d, se, m = local_whittle(eps, power)
    return 1 - 2 * d, 2 * se, m


def pool_days(g: np.ndarray, se: np.ndarray):
    """Inverse-variance mean, its SE, and Cochran's Q with its degrees of freedom."""
    g, se = np.asarray(g, float), np.asarray(se, float)
    w = 1 / se**2
    mean = np.sum(w * g) / np.sum(w)
    Q = float(np.sum(w * (g - mean) ** 2))
    return mean, 1 / np.sqrt(np.sum(w)), Q, len(g) - 1
