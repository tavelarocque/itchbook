"""Data-collapse quality (Houdayer & Hartmann, PRB 70, 014418, 2004) and the
collapse exponent chi of an imbalance I_T.

Quality of a set of curves i with points (x_ij, y_ij, dy_ij): for every point,
take from each OTHER curve i' the two points bracketing x_ij,
x_i'j' <= x_ij < x_i'(j'+1) (the last pair when x_ij equals the last x of i';
curves that do not reach x_ij contribute nothing), and fit a straight line through all of them by
weighted least squares (weights w = 1/dy^2). With K = sum w, Kx = sum w x,
Ky = sum w y, Kxx = sum w x^2, Kxy = sum w x y and D = K Kxx - Kx^2,

    Y_ij   = [(Kxx Ky - Kx Kxy) + x_ij (K Kxy - Kx Ky)] / D,
    dY_ij^2 = (Kxx - 2 x_ij Kx + x_ij^2 K) / D,
    S = (1/N) sum_ij (y_ij - Y_ij)^2 / (dy_ij^2 + dY_ij^2),

summed over the N points that at least one other curve brackets (the pyfssa
documentation states the same formulas). S is about 1 when the curves agree
within their errors and grows as the collapse fails. S is unchanged when every
y and dy is multiplied by one constant, and the local fit is unchanged by an
affine map of x.

Collapse of I_T, assumed self-similar, I_T =d T^chi Z: the curves are the
empirical quantile functions Q_T(p) on a common grid of p (x = p), rescaled to
y_T(p) = Q_T(p) / (T^chi c), with c(chi) = median over T of IQR(I_T) / T^chi.
The scale c is common to all T: dividing each curve by ITS OWN scale would
cancel T^chi and leave S flat in chi. It is an IQR, never a standard
deviation, because the limit is Levy-like: for a truncated Levy imbalance the
standard deviation grows faster than T^chi. chi_hat = argmin S over a chi grid.

Quantile errors (Siddiqui's sparsity estimate): with n_eff independent windows,

    dQ(p) = [Q(p + h) - Q(p - h)] / (2 h) * sqrt(p (1 - p) / n_eff),  h = n_eff^(-1/3) / 2,

h clipped to keep p +- h inside (0, 1) and widened on a lattice (sums of signs
take every other integer) until the two quantiles differ. Overlapping windows (stride T // 4, as
in moments.py) give the quantiles; n_eff counts disjoint windows. Errors at
different p of one curve are correlated; S ignores that, as in the original.
"""
from __future__ import annotations

import numpy as np

import moments


def local_fit(x: np.ndarray, y: np.ndarray, dy: np.ndarray, x0: float) -> tuple[float, float]:
    """Weighted straight-line fit through (x, y, dy), evaluated at x0: (Y, dY)."""
    w = 1.0 / np.asarray(dy, dtype=float) ** 2
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    K, Kx, Ky, Kxx, Kxy = w.sum(), w @ x, w @ y, w @ (x * x), w @ (x * y)
    D = K * Kxx - Kx * Kx
    Y = ((Kxx * Ky - Kx * Kxy) + x0 * (K * Kxy - Kx * Ky)) / D
    return float(Y), float(np.sqrt((Kxx - 2 * x0 * Kx + x0 * x0 * K) / D))


def quality(xs, ys, dys) -> tuple[float, int]:
    """(S, N) for curves given as lists of 1-D arrays (each x sorted ascending)."""
    xs = [np.asarray(v, dtype=float) for v in xs]
    ys = [np.asarray(v, dtype=float) for v in ys]
    ws = [1.0 / np.asarray(v, dtype=float) ** 2 for v in dys]
    total, count = 0.0, 0
    for i, x0 in enumerate(xs):
        sums = np.zeros((5, len(x0)))  # K, Kx, Ky, Kxx, Kxy at each point of curve i
        for k in range(len(xs)):
            if k == i:
                continue
            xk, yk, wk = xs[k], ys[k], ws[k]
            if len(xk) < 2:
                continue
            inside = (x0 >= xk[0]) & (x0 <= xk[-1])
            j = np.clip(np.searchsorted(xk, x0, side="right") - 1, 0, len(xk) - 2)
            for jj in (j, j + 1):
                w = np.where(inside, wk[jj], 0.0)
                sums += np.array([w, w * xk[jj], w * yk[jj], w * xk[jj] ** 2, w * xk[jj] * yk[jj]])
        K, Kx, Ky, Kxx, Kxy = sums
        sel = K > 0
        if not np.any(sel):
            continue
        K, Kx, Ky, Kxx, Kxy, xi = K[sel], Kx[sel], Ky[sel], Kxx[sel], Kxy[sel], x0[sel]
        D = K * Kxx - Kx * Kx
        Y = ((Kxx * Ky - Kx * Kxy) + xi * (K * Kxy - Kx * Ky)) / D
        dY2 = (Kxx - 2 * xi * Kx + xi * xi * K) / D
        total += float(np.sum((ys[i][sel] - Y) ** 2 / (1.0 / ws[i][sel] + dY2)))
        count += int(sel.sum())
    if count == 0:
        raise ValueError("no point of any curve lies inside another curve's range")
    return total / count, count


def window_sums(x: np.ndarray, T: int) -> np.ndarray:
    """I_T(n) = x_n + ... + x_{n+T-1} for n = 0, stride(T), 2 stride(T), ..."""
    c = np.concatenate([[0.0], np.cumsum(np.asarray(x, dtype=float))])
    n = np.arange(0, len(x) - T + 1, moments.stride(T))
    return c[n + T] - c[n]


def sparsity(s: np.ndarray, p: np.ndarray, h0: float) -> np.ndarray:
    """[Q(p + h) - Q(p - h)] / (2 h), h = h0 clipped inside (0, 1); on a lattice (sums of
    signs at small T) the two quantiles can coincide, and h is then doubled until they differ."""
    out = np.empty(len(p))
    for i, pi in enumerate(p):
        h = min(h0, pi / 2, (1 - pi) / 2)
        while True:
            lo, hi = np.quantile(s, [pi - h, pi + h])
            if hi > lo or h >= min(pi, 1 - pi) * 0.999:
                break
            h = min(2 * h, min(pi, 1 - pi) * 0.999)
        if not hi > lo:
            raise ValueError(f"degenerate sample: no spread around p = {pi}")
        out[i] = (hi - lo) / (2 * h)
    return out


def quantile_curves(series, T_grid, p) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Quantiles Q[iT, ip] of I_T pooled over the given series (days), their errors
    dQ, and IQR[iT]. n_eff = number of disjoint T-windows summed over the series."""
    p = np.asarray(p, dtype=float)
    T_grid = np.asarray(T_grid, dtype=int)
    Q = np.empty((len(T_grid), len(p)))
    dQ = np.empty_like(Q)
    iqr = np.empty(len(T_grid))
    for a, T in enumerate(T_grid):
        s = np.concatenate([window_sums(x, int(T)) for x in series])
        n_eff = sum(len(x) // int(T) for x in series)
        Q[a] = np.quantile(s, p)
        dQ[a] = sparsity(s, p, 0.5 * n_eff ** (-1 / 3)) * np.sqrt(p * (1 - p) / n_eff)
        lo, hi = np.quantile(s, [0.25, 0.75])
        iqr[a] = hi - lo
    return Q, dQ, iqr


def rescale(Q, dQ, iqr, T_grid, chi: float) -> tuple[np.ndarray, np.ndarray]:
    """y_T(p) = Q_T(p) / (T^chi c), same for dQ, c = median_T IQR_T / T^chi (common to all T)."""
    s = np.asarray(T_grid, dtype=float)[:, None] ** chi
    c = np.median(np.asarray(iqr)[:, None] / s)
    return Q / (s * c), dQ / (s * c)


def collapse_S(Q, dQ, iqr, T_grid, p, chi: float) -> float:
    y, dy = rescale(Q, dQ, iqr, T_grid, chi)
    return quality([p] * len(y), list(y), list(dy))[0]


def chi_hat(Q, dQ, iqr, T_grid, p, chi_grid=None) -> tuple[float, np.ndarray, np.ndarray]:
    """argmin S over chi_grid (default 0.25..1.25 step 0.0025): (chi_hat, chi_grid, S)."""
    chi_grid = np.arange(0.25, 1.2501, 0.0025) if chi_grid is None else np.asarray(chi_grid, dtype=float)
    S = np.array([collapse_S(Q, dQ, iqr, T_grid, p, c) for c in chi_grid])
    return float(chi_grid[np.argmin(S)]), chi_grid, S
