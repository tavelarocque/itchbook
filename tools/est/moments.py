"""Uncentered within-day window moments of the generalized imbalance and the mid.

For one stock-day with runs n = 0..N-1 (sign eps_n, clipped volume q_n, log mid
lm_n read just before run n, validity flag ok_n), a window of T runs starting
at n covers runs n..n+T-1:

    I^a_T(n) = sum_{k=n}^{n+T-1} eps_k q_k^a,      Delta_T(n) = lm_{n+T} - lm_n

Windows start every stride(T) = max(1, T // 4) runs and need ok at n and n+T.
Moments are plain means over windows, never centered (centering under long
memory biases large-T moments by a term of order (T/N)^gamma):

    V_a(T) = mean (I^a_T)^2,   K_a(T) = mean Delta_T I^a_T,   S(T) = mean Delta_T^2
    v_a(T) = V_a(T) / (T V_a(1)),   k_a(T) = K_a(T) / (T K_a(1))

v and k are scale-free and equal 1 at every T under independent runs; their log
slopes are zeta_V(a) - 1 and zeta_K(a) - 1.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def stride(T: int) -> int:
    return max(1, T // 4)


def clip_volume(q: np.ndarray, day_volume: float, frac: float | None) -> np.ndarray:
    """Clip run volumes at frac * day_volume (frac None: no clipping)."""
    q = np.asarray(q, dtype=float)
    return q if frac is None else np.minimum(q, frac * day_volume)


@dataclass
class DayMoments:
    T: np.ndarray  # (nT,)
    a: np.ndarray  # (na,)
    V: np.ndarray  # (na, nT)
    K: np.ndarray  # (na, nT)
    S: np.ndarray  # (nT,)
    n_windows: np.ndarray  # (nT,)

    def R(self) -> np.ndarray:
        with np.errstate(invalid="ignore", divide="ignore"):
            return self.K / np.sqrt(self.S[None, :] * self.V)


def day_moments(eps, q, lm, ok, T_grid, a_grid) -> DayMoments:
    eps = np.asarray(eps, dtype=float)
    q = np.asarray(q, dtype=float)
    lm = np.asarray(lm, dtype=float)
    ok = np.asarray(ok, dtype=bool)
    N = len(eps)
    T_grid = np.asarray(T_grid, dtype=int)
    a_grid = np.asarray(a_grid, dtype=float)
    V = np.full((len(a_grid), len(T_grid)), np.nan)
    K = np.full_like(V, np.nan)
    S = np.full(len(T_grid), np.nan)
    nw = np.zeros(len(T_grid), dtype=int)
    cums = []
    for a in a_grid:
        x = eps * q**a
        cums.append(np.concatenate([[0.0], np.cumsum(x)]))
    for j, T in enumerate(T_grid):
        if T < 1 or T > N - 1:
            continue
        n = np.arange(0, N - T, stride(int(T)))
        n = n[ok[n] & ok[n + T]]
        nw[j] = len(n)
        if len(n) == 0:
            continue
        D = lm[n + T] - lm[n]
        S[j] = np.mean(D * D)
        for i, c in enumerate(cums):
            I = c[n + T] - c[n]
            V[i, j] = np.mean(I * I)
            K[i, j] = np.mean(D * I)
    return DayMoments(T_grid, a_grid, V, K, S, nw)


def normalized(m: DayMoments) -> tuple[np.ndarray, np.ndarray]:
    """v_a(T), k_a(T); needs T = 1 in the grid."""
    j1 = np.flatnonzero(m.T == 1)
    if len(j1) != 1:
        raise ValueError("T grid must contain 1")
    j1 = j1[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        v = m.V / (m.T[None, :] * m.V[:, [j1]])
        k = m.K / (m.T[None, :] * m.K[:, [j1]])
    return v, k


def log_grid(t_max: int, per_decade: int = 8) -> np.ndarray:
    """Integer horizons 1..t_max, roughly log-spaced, always containing 1 and t_max."""
    g = np.unique(np.round(np.logspace(0, np.log10(t_max), int(np.log10(t_max) * per_decade) + 1)).astype(int))
    return np.unique(np.concatenate([[1], g[g <= t_max], [t_max]]))


def local_slopes(T: np.ndarray, y: np.ndarray) -> np.ndarray:
    """d log y / d log T between adjacent horizons (y > 0); shape (..., nT - 1)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.diff(np.log(y), axis=-1) / np.diff(np.log(T))
