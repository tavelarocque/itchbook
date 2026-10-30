"""The primary pooled cell estimator (pre-registered); Phase 4 uses this code.

Input: one cell's stock-days, every stock with the same study days, given as
concatenated run arrays (eps, q, lm, ok) and per stock-day (start, length).

1. Moments per stock-day (panel_moments), windows of T runs starting every
   stride(T) = max(1, T // 4) runs and never crossing a day, exactly as
   moments.day_moments:
       V_a(T) = mean I^2, K_a(T) = mean Delta I, S(T) = mean Delta^2   (primary: uncentred)
   centered=True subtracts the stock-day's window means first (robustness arm).
2. Normalisation (normalize): v_a(T) = V_a(T) / (T V_a(1)), k_a(T) = K_a(T) /
   (T K_a(1)), and z_a(T) = atanh(K / sqrt(S V)) (the Fisher z of R). A stock-day
   whose K_a(1) is not positive at the 95 % level (one-sided t of the T = 1
   window products below 1.6449) is excluded from k at that a; the count is
   reported. v and z keep every stock-day.
3. Pooling (core), for every eligibility mask m:
   - stock level: each stock's days are averaged first (equal weight per day
     with a finite value), v_s(T), k_s(T), z_s(T);
   - matched increments: between adjacent horizons T_j < T_j+1, the increment
     is the equal-weight mean over stocks eligible at T_j+1 (hence at T_j) of
     the stock's OWN log v_s(T_j+1) - log v_s(T_j); logs are never taken of a
     cross-stock average. The pooled log curve is the chain of these
     increments (log v(1) = 0);
   - zeta over [T1, T2] = 1 + OLS slope of the chained curve on log T over the
     grid points in the range (only the increments inside the range enter);
   - R_a(T) = tanh of the mean over eligible stocks of z_s(T); R_{1/2}/R_0 at
     the range's middle horizon 2^floor((log2 T1 + log2 T2) / 2); a*(T) = vertex
     of the least-squares quadratic in a of R_a(T) on a in [0, min(1.5, a_max)]
     (NaN unless the quadratic is concave);
   - P-F slope = OLS slope of zeta_K(a) on a over [0, a_max]; candidate ratio =
     slope of zeta_V(a) / slope of zeta_K(a), both OLS on a in [0, 0.5], and the
     difference vk_diff = slope zeta_V - 2 slope zeta_K (defined when the zeta_K
     slope is 0, where the ratio is not).
   Masks: ("all",) every stock; ("ceil", rho) stock s at T iff T <= rho N_s,
   N_s = its fewest runs on any study day; ("band", lo, hi) iff lo < T / N_s <= hi.
   An increment's eligibility is decided at its upper horizon T_j+1.
4. Inference (intervals): the day-and-stock cluster bootstrap of
   inference.cluster_bootstrap (one draw of days shared by all stocks, stocks
   drawn with replacement, the same random draws in the same order), B = 200,
   percentile interval; and the delete-one-day jackknife with its t interval
   (n_days - 1 df). The reported interval is the wider of the two, per element.
   Unlike inference.py, a replicate that is undefined for one element is
   dropped for that element only (outputs here are vectors with undefined
   bands), and an element whose jackknife is undefined has no interval.

Every pooling step is a weighted sum with day multiplicities Wd (B, D) and stock
multiplicities Ws (B, S), so the point estimate (all weights 1), the bootstrap
(multiplicities of the draws) and the jackknife (one day weight 0) run through
the same core.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

import inference as INF
import moments as MO

T_GRID = 2 ** np.arange(13)
A_GRID = np.round(np.arange(0.0, 1.501, 0.25), 2)
RANGES = ((8, 64), (16, 256), (64, 1024))
K1_Z = 1.6448536269514722  # one-sided 95 % normal quantile
A_PF = 1.5
A_RATIO = 0.5
N_BOOT = 200
LEVEL = 0.90


# ---- 1. moments ---------------------------------------------------------------------


@dataclass
class PanelMoments:
    V: np.ndarray  # (n_sd, na, nT)
    K: np.ndarray  # (n_sd, na, nT)
    S: np.ndarray  # (n_sd, nT)
    nw: np.ndarray  # (n_sd, nT) windows
    K1_t: np.ndarray  # (n_sd, na) t statistic of K_a(1)
    N: np.ndarray  # (n_sd,) runs


def _windows(starts: np.ndarray, lengths: np.ndarray, T: int):
    """Window starts n = start + k stride(T), 0 <= k stride < N - T, and the stock-day of each."""
    s = MO.stride(int(T))
    cnt = np.where(lengths > T, (lengths - T + s - 1) // s, 0)
    sd = np.repeat(np.arange(len(lengths)), cnt)
    within = np.arange(int(cnt.sum())) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    return starts[sd] + s * within, sd


def panel_moments(eps, q, lm, ok, starts, lengths, T_grid=T_GRID, a_grid=A_GRID, centered=False) -> PanelMoments:
    eps, q, lm = (np.asarray(x, dtype=float) for x in (eps, q, lm))
    ok = np.asarray(ok, dtype=bool)
    starts, lengths = np.asarray(starts, dtype=np.int64), np.asarray(lengths, dtype=np.int64)
    n_sd, na, nT = len(lengths), len(a_grid), len(T_grid)
    if T_grid[0] != 1:
        raise ValueError("T grid must start at 1")
    cums = [np.concatenate([[0.0], np.cumsum(eps * q**a)]) for a in a_grid]
    V = np.full((n_sd, na, nT), np.nan)
    K = np.full_like(V, np.nan)
    S = np.full((n_sd, nT), np.nan)
    nw = np.zeros((n_sd, nT), dtype=np.int64)
    K1_t = np.full((n_sd, na), np.nan)

    def mean(x, sd, c):
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.bincount(sd, x, minlength=n_sd) / c

    for j, T in enumerate(T_grid):
        n, sd = _windows(starts, lengths, int(T))
        keep = ok[n] & ok[n + T]
        n, sd = n[keep], sd[keep]
        c = np.bincount(sd, minlength=n_sd).astype(float)
        c[c == 0] = np.nan
        nw[:, j] = np.nan_to_num(c).astype(np.int64)
        D = lm[n + T] - lm[n]
        mD = mean(D, sd, c)
        S[:, j] = mean(D * D, sd, c) - (mD * mD if centered else 0.0)
        for i, cum in enumerate(cums):
            I = cum[n + T] - cum[n]
            mI = mean(I, sd, c)
            V[:, i, j] = mean(I * I, sd, c) - (mI * mI if centered else 0.0)
            p = D * I
            mp = mean(p, sd, c)
            K[:, i, j] = mp - (mD * mI if centered else 0.0)
            if T == 1:
                var = (mean(p * p, sd, c) - mp * mp) * c / (c - 1)
                with np.errstate(invalid="ignore", divide="ignore"):
                    K1_t[:, i] = mp / np.sqrt(var / c)
    return PanelMoments(V, K, S, nw, K1_t, lengths)


# ---- 2. normalisation ---------------------------------------------------------------


def normalize(pm: PanelMoments, T_grid=T_GRID) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(v, k, z, excluded): per stock-day arrays (n_sd, na, nT); excluded (n_sd, na) marks K_a(1) not
    positive at the 95 % level (k set to NaN there)."""
    T = np.asarray(T_grid, dtype=float)
    excl = ~(pm.K1_t > K1_Z)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = pm.V / (T[None, None, :] * pm.V[:, :, :1])
        k = pm.K / (T[None, None, :] * pm.K[:, :, :1])
        R = pm.K / np.sqrt(pm.S[:, None, :] * pm.V)
    k[excl] = np.nan
    z = np.arctanh(np.clip(R, -1 + 1e-12, 1 - 1e-12))
    return v, k, z, excl


# ---- 3. pooling -----------------------------------------------------------------------


def masks(Ns: np.ndarray, T: np.ndarray, specs) -> np.ndarray:
    """(S, len(T), n_masks) eligibility of stock s at horizon T under each spec."""
    r = np.asarray(T, dtype=float)[None, :] / np.asarray(Ns, dtype=float)[:, None]
    out = []
    for sp in specs:
        if sp[0] == "all":
            out.append(np.ones_like(r, dtype=bool))
        elif sp[0] == "ceil":
            out.append(r <= sp[1])
        elif sp[0] == "band":
            out.append((r > sp[1]) & (r <= sp[2]))
        else:
            raise ValueError(f"unknown mask {sp!r}")
    return np.stack(out, axis=-1)


def increment_masks(Ns: np.ndarray, T_grid: np.ndarray, specs) -> np.ndarray:
    """Eligibility of each increment T_j -> T_j+1, decided at the upper horizon T_j+1."""
    return masks(Ns, T_grid[1:], specs)


def _stock_mean(x, wd):
    """Weighted mean over each stock's days: x (S, D, ...) with weights wd (B, S, D) -> (B, S, ...)."""
    fin = np.isfinite(x)
    num = np.einsum("bsd,sd...->bs...", wd, np.where(fin, x, 0.0))
    den = np.einsum("bsd,sd...->bs...", wd, fin.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den


def _increments(lv, ws, E):
    """Matched increments (B, na, nJ, nm): mean over eligible stocks of each stock's own log change."""
    inc = lv[..., 1:] - lv[..., :-1]  # one stock's own change between adjacent horizons
    fin = np.isfinite(inc)
    num = np.einsum("bs,sjm,bsaj->bajm", ws, E, np.where(fin, inc, 0.0))
    den = np.einsum("bs,sjm,bsaj->bajm", ws, E, fin.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        return num / den, den


def _ols_weights(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    d = x - x.mean()
    return d / (d @ d)


def range_weights(T_grid=T_GRID, ranges=RANGES) -> np.ndarray:
    """(n_ranges, nJ): zeta - 1 = sum_j w_j Delta_j, the OLS slope of the chained curve over the range."""
    lt = np.log(np.asarray(T_grid, dtype=float))
    W = np.zeros((len(ranges), len(T_grid) - 1))
    for r, (lo, hi) in enumerate(ranges):
        sel = np.flatnonzero((T_grid >= lo) & (T_grid <= hi))
        c = _ols_weights(lt[sel])
        for j in range(sel[0], sel[-1]):  # L(T_t) = sum_{j < t} Delta_j; only in-range increments matter
            W[r, j] = c[sel > j].sum()
    return W


def core(v, k, z, Nsd, Wd, Ws, specs, T_grid=T_GRID, a_grid=A_GRID, ranges=RANGES, a_max=A_PF) -> dict:
    """Pooled quantities with a leading replicate axis B (one per row of Wd, Ws)."""
    S_, D_, na, nT = v.shape
    Ns = Nsd.min(axis=1)
    wd, ws = Wd[:, None, :] * np.ones(Nsd.shape)[None], Ws * 1.0
    E = increment_masks(Ns, T_grid, specs)
    M = masks(Ns, T_grid, specs)
    Wr = range_weights(T_grid, ranges)
    out = {}
    for name, x in (("V", v), ("K", k)):
        xs = _stock_mean(x, wd)
        with np.errstate(invalid="ignore", divide="ignore"):
            lx = np.where(xs > 0, np.log(np.where(xs > 0, xs, 1.0)), np.nan)
        dlt, den = _increments(lx, ws, E)
        out[f"n_inc_{name}"] = den
        # zeta (B, na, nr, nm): NaN when any in-range increment is undefined
        zeta = np.full((len(Ws), na, len(ranges), len(specs)), np.nan)
        for r in range(len(ranges)):
            js = np.flatnonzero(Wr[r])
            zeta[:, :, r, :] = 1 + np.einsum("bajm,j->bam", dlt[:, :, js, :], Wr[r, js])
        out[f"zeta_{name}"] = zeta
    zs = _stock_mean(z, wd)  # (B, S, na, nT)
    fin = np.isfinite(zs)
    num = np.einsum("bs,stm,bsat->batm", ws, M, np.where(fin, zs, 0.0))
    den = np.einsum("bs,stm,bsat->batm", ws, M, fin.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        R = np.tanh(num / den)
    out["R"] = R
    lt2 = np.log2(np.asarray(T_grid, dtype=float))
    mids = [int(np.flatnonzero(lt2 == np.floor((np.log2(lo) + np.log2(hi)) / 2))[0]) for lo, hi in ranges]
    ia0, ih = int(np.flatnonzero(a_grid == 0)[0]), int(np.flatnonzero(a_grid == 0.5)[0])
    with np.errstate(invalid="ignore", divide="ignore"):
        out["Rratio"] = np.stack([R[:, ih, t, :] / R[:, ia0, t, :] for t in mids], axis=1)  # (B, nr, nm)
    qa = a_grid <= min(1.5, a_max)
    X = np.stack([np.ones(qa.sum()), a_grid[qa], a_grid[qa] ** 2], axis=1)
    coef = np.einsum("ca,batm->bctm", np.linalg.pinv(X), R[:, qa])
    with np.errstate(invalid="ignore", divide="ignore"):
        out["a_star"] = np.where(coef[:, 2] < 0, -coef[:, 1] / (2 * coef[:, 2]), np.nan)  # (B, nT, nm)
    wpf = _ols_weights(a_grid[a_grid <= a_max])
    out["pf"] = np.einsum("barm,a->brm", out["zeta_K"][:, a_grid <= a_max], wpf)
    wra = _ols_weights(a_grid[a_grid <= A_RATIO])
    sV = np.einsum("barm,a->brm", out["zeta_V"][:, a_grid <= A_RATIO], wra)
    sK = np.einsum("barm,a->brm", out["zeta_K"][:, a_grid <= A_RATIO], wra)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["ratio"] = sV / sK
    out["vk_diff"] = sV - 2 * sK  # the candidate's equality, defined even when the zeta_K slope is 0
    return out


ESTIMANDS = ("zeta_V", "zeta_K", "pf", "ratio", "vk_diff", "Rratio", "R", "a_star")


def to_panel(x: np.ndarray, n_days: int) -> np.ndarray:
    """(n_sd, ...) in stock-major order -> (S, D, ...)."""
    return x.reshape(len(x) // n_days, n_days, *x.shape[1:])


def estimate(v, k, z, Nsd, specs, **kw) -> dict:
    S_, D_ = Nsd.shape
    out = core(v, k, z, Nsd, np.ones((1, D_)), np.ones((1, S_)), specs, **kw)
    return {key: val[0] for key, val in out.items()}


def bootstrap_weights(n_s: int, n_d: int, B: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Multiplicities of inference.cluster_bootstrap's draws (one stratum): per replicate, stocks first,
    then days shared by every stock, from the same generator calls."""
    Ws, Wd = np.zeros((B, n_s)), np.zeros((B, n_d))
    for b in range(B):
        Ws[b] = np.bincount(rng.integers(0, n_s, n_s), minlength=n_s)
        Wd[b] = np.bincount(rng.integers(0, n_d, n_d), minlength=n_d)
    return Wd, Ws


def jackknife_weights(n_s: int, n_d: int) -> tuple[np.ndarray, np.ndarray]:
    """Delete one day at a time: (Wd, Ws) with n_d rows."""
    return 1.0 - np.eye(n_d), np.ones((n_d, n_s))


def jackknife_se(v, k, z, Nsd, specs, **kw) -> tuple[dict, int]:
    """Delete-one-day jackknife SE of every estimand (NaN where a leave-one-out value is undefined), and
    the number of deleted units."""
    S_, D_ = Nsd.shape
    jWd, jWs = jackknife_weights(S_, D_)
    jack = core(v, k, z, Nsd, jWd, jWs, specs, **kw)
    n = jWd.shape[0]
    out = {}
    for key in ESTIMANDS:
        loo = jack[key]
        se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean(axis=0)) ** 2, axis=0))
        se[~np.all(np.isfinite(loo), axis=0)] = np.nan
        out[key] = se
    return out, n


def intervals(v, k, z, Nsd, specs, rng: np.random.Generator, n_boot: int = N_BOOT, level: float = LEVEL,
              chunk: int = 50, **kw) -> dict:
    """{estimand: (estimate, lo, hi, method)} with the wider of bootstrap and day jackknife per element;
    method 1 = bootstrap, 2 = jackknife, 0 = no interval."""
    S_, D_ = Nsd.shape
    est = estimate(v, k, z, Nsd, specs, **kw)
    Wd, Ws = bootstrap_weights(S_, D_, n_boot, rng)
    boot = {key: [] for key in ESTIMANDS}
    for i in range(0, n_boot, chunk):
        o = core(v, k, z, Nsd, Wd[i : i + chunk], Ws[i : i + chunk], specs, **kw)
        for key in ESTIMANDS:
            boot[key].append(o[key])
    se, n = jackknife_se(v, k, z, Nsd, specs, **kw)
    c = INF.t_ppf((1 + level) / 2, n - 1)
    out = {}
    for key in ESTIMANDS:
        b = np.concatenate(boot[key])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            blo = np.nanpercentile(b, 100 * (1 - level) / 2, axis=0)
            bhi = np.nanpercentile(b, 100 * (1 + level) / 2, axis=0)
        jok = np.isfinite(se[key])
        jlo, jhi = est[key] - c * se[key], est[key] + c * se[key]
        use_j = ~((bhi - blo) >= (jhi - jlo))  # the wider; ties to the bootstrap
        lo = np.where(use_j, jlo, blo)
        hi = np.where(use_j, jhi, bhi)
        method = np.where(use_j, 2, 1)
        undefined = ~np.isfinite(est[key]) | ~jok | ~np.isfinite(blo)
        lo[undefined], hi[undefined], method[undefined] = np.nan, np.nan, 0
        out[key] = (est[key], lo, hi, method)
    return out
