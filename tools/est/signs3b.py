"""Phase 3b sign model: centred and uncentred sign ACF, the day-mean component, a
two-component sign generator, and gamma estimators run on it. Signs alone.

Day-mean decomposition. For stock-day d with N runs and mean sign m_d, the
uncentred pooled C(l) mixes memory within the day with the spread of m_d across
stock-days. centred_acf demeans every stock-day by its OWN mean before the lag
products (pairs never cross a stock-day, as in calib3.pooled_acf).
mean_sign_scaling bins stock-days by N and fits log var(m_d) on log N: iid signs
give slope -1, a fixed per-day component gives slope 0.

Two-component generator, one stock-day of N runs:

    eps_t = sign(X_t + m_d),   X = fGn(H) with unit variance,   m_d ~ N(0, tau^2) per stock-day.

Over stock-days (X_t + m, X_{t+l} + m) is bivariate normal with variance
1 + tau^2 and covariance rho_H(l) + tau^2, so by the arcsine law the
uncentred C(l) is exactly

    C(l) = (2/pi) arcsin( (rho_H(l) + tau^2) / (1 + tau^2) ).

It tends to (2/pi) arcsin(tau^2 / (1 + tau^2)) at long lags, and C(l) minus that
plateau decays as l^-(2 - 2H): the generator's long-lag exponent is
gamma_true = 2 - 2H. fit_two_component fits (H, tau) by least squares on C over
a log-lag grid (equal weight per lag of LOG_LAGS present), first on H_GRID; if
the best H sits on an edge it refits on H_GRID_WIDE and reports both flags.

gamma estimators (estimator_study), run on simulated panels with the cell's
actual stock-day N values:
  lw     longmem.gamma_lw on stock-days with N >= calib3.LW_MIN_N, pooled by
         inverse variance (calib3's rule)
  power  calib3.sign_acf_fit (power-law OLS of the uncentred C on [10, l_hi]),
         bootstrap B = POWER_BOOT for the l_hi cut
  v0     2 - zeta_V(0), zeta_V(0) = 1 + OLS slope of log v_0(T) on log T over
         V0_RANGE, v_0 from the window-pooled V_0(T) of stock-days with N >=
         V0_MIN_N (a fixed set at every T)
"""
from __future__ import annotations

import numpy as np

import calib3 as C3
import gen3 as G
import longmem as LM
import moments as MO
import synth

H_GRID = np.round(np.arange(0.505, 0.9951, 0.005), 4)
H_GRID_WIDE = np.round(np.arange(0.01, 0.9991, 0.001), 4)
TAU_GRID = np.round(np.arange(0.0, 0.6001, 0.0025), 5)
LOG_LAGS = np.unique(np.round(10.0 ** (np.arange(0, 41) / 12)).astype(int))
POWER_BOOT = 199
V0_RANGE = (8, 64)
V0_MIN_N = 256
V0_T = 2 ** np.arange(0, 7)


# ---- day-mean decomposition --------------------------------------------------------


def centred_acf(series, lags: np.ndarray):
    """pooled_acf of each stock-day demeaned by its own mean sign: (C, num, cnt)."""
    return C3.pooled_acf([x - x.mean() for x in series], lags)


def mean_sign_scaling(means: np.ndarray, Ns: np.ndarray, n_bins: int = 8, min_per_bin: int = 10) -> dict:
    """var(m_d) in quantile bins of log N (bins with >= min_per_bin stock-days), and the OLS of
    log var on log (median N of the bin)."""
    means, Ns = np.asarray(means, float), np.asarray(Ns, float)
    edges = np.unique(np.quantile(np.log(Ns), np.linspace(0, 1, n_bins + 1)))
    b = np.clip(np.searchsorted(edges, np.log(Ns), side="right") - 1, 0, len(edges) - 2)
    rows = []
    for k in range(len(edges) - 1):
        sel = b == k
        if sel.sum() >= min_per_bin:
            rows.append((float(np.median(Ns[sel])), float(np.var(means[sel])), int(sel.sum())))
    out = {"var_all": float(np.var(means)), "bins": rows, "slope": np.nan, "intercept": np.nan}
    if len(rows) >= 3 and all(r[1] > 0 for r in rows):
        x, y = np.log([r[0] for r in rows]), np.log([r[1] for r in rows])
        out["slope"], out["intercept"] = (float(v) for v in np.polyfit(x, y, 1))
    return out


# ---- two-component generator ------------------------------------------------------


def two_component_C(l, H: float, tau: float) -> np.ndarray:
    tau2 = tau * tau
    return (2 / np.pi) * np.arcsin(np.clip((G.fgn_rho(l, H) + tau2) / (1 + tau2), -1, 1))


def _grid_fit(l: np.ndarray, C: np.ndarray, H_grid: np.ndarray, tau_grid: np.ndarray):
    best = (np.inf, np.nan, np.nan)
    t2 = tau_grid[:, None] ** 2
    for H in H_grid:
        model = (2 / np.pi) * np.arcsin(np.clip((G.fgn_rho(l, H)[None, :] + t2) / (1 + t2), -1, 1))
        rss = np.sum((model - C[None, :]) ** 2, axis=1)
        i = int(np.argmin(rss))
        if rss[i] < best[0]:
            best = (float(rss[i]), float(H), float(tau_grid[i]))
    return best


def fit_two_component(lags: np.ndarray, C: np.ndarray, H_grid=H_GRID, H_grid_wide=H_GRID_WIDE,
                      tau_grid=TAU_GRID) -> dict:
    """(H, tau) by least squares on C over the LOG_LAGS present in lags; refit on the wide H grid when
    the first fit sits on an edge of H_grid. edge_first / edge_final flag a fit on an H grid edge."""
    lags, C = np.asarray(lags), np.asarray(C, float)
    use = np.isin(lags, LOG_LAGS) & np.isfinite(C)
    l, c = lags[use].astype(float), C[use]
    rss, H, tau = _grid_fit(l, c, H_grid, tau_grid)
    edge_first = bool(H <= H_grid[0] or H >= H_grid[-1])
    grid = H_grid
    if edge_first:
        rss, H, tau = _grid_fit(l, c, H_grid_wide, tau_grid)
        grid = H_grid_wide
    return {"H": H, "tau": tau, "gamma_true": 2 - 2 * H, "rmse": float(np.sqrt(rss / len(l))), "n_lags": int(len(l)),
            "edge_first": edge_first, "edge_final": bool(H <= grid[0] or H >= grid[-1]),
            "tau_edge": bool(tau >= tau_grid[-1])}


def simulate_two_component(Ns, H: float, tau: float, rng: np.random.Generator) -> list[np.ndarray]:
    out = []
    for N in Ns:
        x = synth.fgn(int(N), H, rng) + tau * rng.standard_normal()
        out.append(np.where(x >= 0, 1.0, -1.0))
    return out


# ---- gamma estimators -------------------------------------------------------------


def gamma_estimators(series, rng: np.random.Generator) -> dict:
    """The three estimators on one panel of sign series (one per stock-day)."""
    Ns = np.array([len(x) for x in series])
    out = {}
    lw = [LM.gamma_lw(x)[:2] for x in series if len(x) >= C3.LW_MIN_N]
    out["lw"] = LM.pool_days(*np.array(lw).T)[0] if len(lw) >= 2 else np.nan
    lags = C3.sign_lags()
    _, num, cnt = C3.pooled_acf(series, lags)
    out["power"] = C3.sign_acf_fit(lags, num, cnt, Ns, rng, POWER_BOOT)["gamma_fit"]
    long = [x for x in series if len(x) >= V0_MIN_N]
    if len(long) >= 2:
        V = np.zeros(len(V0_T))
        W = np.zeros(len(V0_T))
        for x in long:
            m = MO.day_moments(x, np.ones(len(x)), np.zeros(len(x)), np.ones(len(x), bool), V0_T, [0.0])
            V += m.V[0] * m.n_windows
            W += m.n_windows
        v = (V / W) / (V0_T * (V / W)[0])
        sel = (V0_T >= V0_RANGE[0]) & (V0_T <= V0_RANGE[1])
        out["v0"] = float(1 - np.polyfit(np.log(V0_T[sel]), np.log(v[sel]), 1)[0])
    else:
        out["v0"] = np.nan
    return out


def estimator_study(Ns, H: float, tau: float, n_rep: int, seed: int) -> np.ndarray:
    """(n_rep, 3) estimates (lw, power, v0) on panels of the fitted generator with stock-day lengths Ns."""
    out = np.empty((n_rep, 3))
    for r in range(n_rep):
        rng = np.random.default_rng([seed, r])
        e = gamma_estimators(simulate_two_component(Ns, H, tau, rng), rng)
        out[r] = (e["lw"], e["power"], e["v0"])
    return out
