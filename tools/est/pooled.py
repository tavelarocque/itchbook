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
   Range reach (core2, from 31 Oct 2026; see MIN_REACH): an estimate over [T1, T2] is defined only if at
   least 5 stocks are eligible at T2, and zeta_K only if every horizon of the range is kept.
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
    out = {}
    for key in ESTIMANDS:
        b = np.concatenate(boot[key])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            blo = np.nanpercentile(b, 100 * (1 - level) / 2, axis=0)
            bhi = np.nanpercentile(b, 100 * (1 + level) / 2, axis=0)
        jok = np.isfinite(se[key])
        jlo, jhi = jk_interval(est[key], se[key], n, level)
        use_j = ~((bhi - blo) >= (jhi - jlo))  # the wider; ties to the bootstrap
        lo = np.where(use_j, jlo, blo)
        hi = np.where(use_j, jhi, bhi)
        method = np.where(use_j, 2, 1)
        undefined = ~np.isfinite(est[key]) | ~jok | ~np.isfinite(blo)
        lo[undefined], hi[undefined], method[undefined] = np.nan, np.nan, 0
        out[key] = (est[key], lo, hi, method)
    return out


# ---- Phase 3c-2: the k-tilde normalisation (primary) ---------------------------------------
#
# norm="tilde" (primary since 3c-2; pooling of v revised in 3c-3): per stock s, with day-averaged single-run
# moments Sbar(1), Vbar_a(1),
#     v_a(T) = Vbar_a(T) / (T Vbar_a(1)),   k~_a(T) = Kbar_a(T) / (T sqrt(Sbar(1) Vbar_a(1))),
# where Xbar is the stock's day average (equal weight per day with a value) under the replicate's day
# weights. No stock-day is excluded for the sign of K_a(1). Pooling of v (3c-3): the equal-weight mean
# over stocks of each stock's own log increment, matched sets, chained (zeta_mean_of_logs). Pooling of
# k~ (it can be negative per stock): equal weight per stock, matched sets, chained logs of matched means:
# between consecutive kept horizons T_p < T_t of a range, the increment is
#     log( mean_s x_s(T_t) ) - log( mean_s x_s(T_p) ),
# both means over the SAME stocks (eligible at T_t, with values at both). The logs are of means of curves
# already normalised by each stock's own constants, never of raw cross-stock averages. A horizon whose
# pooled k~ (mean over the stocks eligible there) is not positive is dropped for that estimate (counted
# in "dropped"); the chain then joins the neighbouring kept horizons, again with a matched set. zeta =
# 1 + OLS slope of the chained curve on log T over the kept horizons of the range (at least MIN_KEPT).
# norm="old" is the Phase 3c path (per-stock-day normalisation, K_a(1) exclusion, mean of stock logs).

T_GRID2 = 2 ** np.arange(9)  # 1 .. 256
A_GRID2 = np.round(np.arange(0.0, 1.001, 0.25), 2)
RANGES2 = ((4, 32), (8, 64), (16, 256))
A_MAX2 = 1.0
MIN_KEPT = 3
ESTIMANDS2 = ("zeta_V", "zeta_K", "pf", "vk_diff", "Rratio")
# Range-reach rule (decided 31 Oct 2026 from docs/phase3c4_diagnosis.md, which found zeta_K "estimated" from
# the lower horizons of a range that no stock reached): an estimate over [T1, T2] is defined for a panel only
# if at least MIN_REACH stocks are eligible at T2 under the mask, and zeta_K only if every horizon of the range
# is kept (never a slope over a shortened range). Undefined elements are NaN.
MIN_REACH = 5
LEGACY_RULE = dict(min_reach=0, full_range=False)  # the estimator before 31 Oct 2026 (3c-2, 3c-3, 3c-4 chunks)


def range_reach(Ns: np.ndarray, specs, ranges=RANGES2) -> np.ndarray:
    """(n_ranges, n_masks): the number of stocks eligible at each range's top horizon T2."""
    return masks(Ns, np.array([hi for _, hi in ranges]), specs).sum(axis=0)


def range_rule(reach: np.ndarray, kept_full: np.ndarray, a_grid=A_GRID2, a_max=A_MAX2, min_reach: int = MIN_REACH,
               full_range: bool = True) -> dict:
    """Where each estimand is defined under the range-reach rule. reach (..., nr, nm) from range_reach;
    kept_full (..., na, nr, nm): zeta_K kept every horizon of the range. Returns key -> boolean mask with the
    estimand's shape: zeta_V, zeta_K (..., na, nr, nm); pf, vk_diff, Rratio (..., nr, nm). The P-F slope and
    vk_diff need zeta_K over a <= a_max (a <= 0.5) to be defined; R_{1/2}/R_0 over [T1, T2] needs the reach."""
    ok = np.asarray(reach) >= min_reach
    zV = np.broadcast_to(ok[..., None, :, :], np.shape(kept_full))
    zK = zV & np.asarray(kept_full) if full_range else zV
    a_grid = np.asarray(a_grid)
    return {"zeta_V": zV, "zeta_K": zK, "pf": zK[..., a_grid <= a_max, :, :].all(axis=-3),
            "vk_diff": zK[..., a_grid <= A_RATIO, :, :].all(axis=-3), "Rratio": ok}


@dataclass
class Raw:
    """One cell's per-stock-day moments in panel shape (S stocks, D days)."""

    V: np.ndarray  # (S, D, na, nT)
    K: np.ndarray  # (S, D, na, nT)
    Sm: np.ndarray  # (S, D, nT)
    K1_t: np.ndarray  # (S, D, na)
    Nsd: np.ndarray  # (S, D)

    @property
    def z(self) -> np.ndarray:
        with np.errstate(invalid="ignore", divide="ignore"):
            R = self.K / np.sqrt(self.Sm[:, :, None, :] * self.V)
        return np.arctanh(np.clip(R, -1 + 1e-12, 1 - 1e-12))


def raw_panel(pm: PanelMoments, n_days: int) -> Raw:
    return Raw(to_panel(pm.V, n_days), to_panel(pm.K, n_days), to_panel(pm.S, n_days), to_panel(pm.K1_t, n_days),
               to_panel(pm.N, n_days))


def jk_interval(est, se, n: int, level: float):
    """est -+ t_{n-1}((1 + level) / 2) se: Student t with n - 1 df for n deleted units."""
    c = INF.t_ppf((1 + level) / 2, n - 1)
    return est - c * se, est + c * se


def stock_curves_tilde(V, K, Sm, wd, T_grid=T_GRID2):
    """Stock-level v and k~ (B, S, na, nT) from day-averaged raw moments (day weights wd (B, S, D))."""
    T = np.asarray(T_grid, dtype=float)
    Vb, Kb, Sb = _stock_mean(V, wd), _stock_mean(K, wd), _stock_mean(Sm, wd)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = Vb / (T * Vb[..., :1])
        k = Kb / (T * np.sqrt(Sb[:, :, None, :1] * Vb[..., :1]))
    return v, k


def _pair_means(x, ws, Mt, t1: int, t2: int):
    """Means over the stocks eligible at t2 with values at t1 and t2, of x at t1 and at t2: (B, na, nm)."""
    fin = (np.isfinite(x[..., t1]) & np.isfinite(x[..., t2])).astype(float)
    den = np.einsum("bs,sm,bsa->bam", ws, Mt[:, t2, :], fin)
    with np.errstate(invalid="ignore", divide="ignore"):
        m1 = np.einsum("bs,sm,bsa->bam", ws, Mt[:, t2, :], np.where(fin > 0, x[..., t1], 0.0)) / den
        m2 = np.einsum("bs,sm,bsa->bam", ws, Mt[:, t2, :], np.where(fin > 0, x[..., t2], 0.0)) / den
    return m1, m2


def chain_zeta(x, ws, Mt, T_grid=T_GRID2, ranges=RANGES2):
    """zeta (B, na, nr, nm) from logs of matched means, dropping horizons whose pooled value is not
    positive; also the pooled value per horizon (B, na, nT, nm) and the number of kept horizons."""
    fin = np.isfinite(x)
    num = np.einsum("bs,stm,bsat->batm", ws, Mt, np.where(fin, x, 0.0))
    den = np.einsum("bs,stm,bsat->batm", ws, Mt, fin.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = num / den
    keep = pooled > 0
    lt = np.log(np.asarray(T_grid, dtype=float))
    B, na, _, nm = pooled.shape
    zeta = np.full((B, na, len(ranges), nm), np.nan)
    n_kept = np.zeros((B, na, len(ranges), nm), dtype=int)
    for r, (lo, hi) in enumerate(ranges):
        idx = np.flatnonzero((T_grid >= lo) & (T_grid <= hi))
        prev = np.full((B, na, nm), -1)
        L = np.zeros((B, na, nm))
        bad = np.zeros((B, na, nm), dtype=bool)
        acc = {k_: np.zeros((B, na, nm)) for k_ in ("n", "x", "xx", "L", "xL")}
        pm = {}
        for p, t in enumerate(idx):
            kt = keep[:, :, t, :]
            inc = np.full((B, na, nm), np.nan)
            for pv in idx[:p]:
                if (pv, t) not in pm:
                    pm[(pv, t)] = _pair_means(x, ws, Mt, pv, t)
                m1, m2 = pm[(pv, t)]
                with np.errstate(invalid="ignore", divide="ignore"):
                    pos = (m1 > 0) & (m2 > 0)
                    d = np.where(pos, np.log(np.where(pos, m2, 1.0)) - np.log(np.where(pos, m1, 1.0)), np.nan)
                inc = np.where(prev == pv, d, inc)
            start = kt & (prev < 0)
            step = kt & (prev >= 0)
            bad |= step & ~np.isfinite(inc)
            ok = step & np.isfinite(inc)
            L = np.where(ok, L + np.nan_to_num(inc), np.where(start, 0.0, L))
            add = start | ok
            for k_, val in (("n", 1.0), ("x", lt[t]), ("xx", lt[t] ** 2)):
                acc[k_] += np.where(add, val, 0.0)
            acc["L"] += np.where(add, L, 0.0)
            acc["xL"] += np.where(add, lt[t] * L, 0.0)
            prev = np.where(add, t, prev)
        n = acc["n"]
        with np.errstate(invalid="ignore", divide="ignore"):
            slope = (n * acc["xL"] - acc["x"] * acc["L"]) / (n * acc["xx"] - acc["x"] ** 2)
        zeta[:, :, r, :] = np.where((n >= MIN_KEPT) & ~bad, 1 + slope, np.nan)
        n_kept[:, :, r, :] = n
    return zeta, pooled, n_kept


def zeta_mean_of_logs(v, ws, Ns, specs, T_grid=T_GRID2, ranges=RANGES2):
    """zeta from the equal-weight mean over stocks of each stock's own log increment (days averaged within
    the stock first), on matched sets (eligibility at the upper horizon), chained (Phase 3c-3 rule for v,
    which is positive for every stock)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        lv = np.where(v > 0, np.log(np.where(v > 0, v, 1.0)), np.nan)
    dlt, _ = _increments(lv, ws, increment_masks(Ns, T_grid, specs))
    Wr = range_weights(T_grid, ranges)
    zeta = np.full((len(ws), v.shape[2], len(ranges), len(specs)), np.nan)
    for r in range(len(ranges)):
        js = np.flatnonzero(Wr[r])
        zeta[:, :, r, :] = 1 + np.einsum("bajm,j->bam", dlt[:, :, js, :], Wr[r, js])
    return zeta


def core2(raw: Raw, Wd, Ws, specs, norm: str = "tilde", T_grid=T_GRID2, a_grid=A_GRID2, ranges=RANGES2,
          a_max=A_MAX2, min_reach: int = MIN_REACH, full_range: bool = True) -> dict:
    """The pooled estimands of ESTIMANDS2 (and R, a_star, diagnostics) with a leading replicate axis. The
    range-reach rule (range_rule) applies on the tilde path; reach is the panel's, the same in every replicate.
    **LEGACY_RULE reproduces the estimator before 31 Oct 2026."""
    if norm == "old":
        T = np.asarray(T_grid, dtype=float)
        excl = ~(raw.K1_t > K1_Z)
        with np.errstate(invalid="ignore", divide="ignore"):
            v = raw.V / (T * raw.V[..., :1])
            k = raw.K / (T * raw.K[..., :1])
        k = np.where(excl[..., None], np.nan, k)
        o = core(v, k, raw.z, raw.Nsd, Wd, Ws, specs, T_grid=T_grid, a_grid=a_grid, ranges=ranges, a_max=a_max)
        return {key: o[key] for key in ESTIMANDS2 + ("R", "a_star")}
    if norm != "tilde":
        raise ValueError("norm is 'tilde' or 'old'")
    Ns = raw.Nsd.min(axis=1)
    wd, ws = Wd[:, None, :] * np.ones(raw.Nsd.shape)[None], Ws * 1.0
    Mt = masks(Ns, T_grid, specs)
    v, k = stock_curves_tilde(raw.V, raw.K, raw.Sm, wd, T_grid)
    out = {}
    out["zeta_V"] = zeta_mean_of_logs(v, ws, Ns, specs, T_grid, ranges)
    out["zeta_K"], pooled_k, out["kept_K"] = chain_zeta(k, ws, Mt, T_grid, ranges)
    out["dropped_K"] = ~(pooled_k > 0) & np.isfinite(pooled_k)  # (B, na, nT, nm): pooled k~ not positive
    zs = _stock_mean(raw.z, wd)
    fin = np.isfinite(zs)
    num = np.einsum("bs,stm,bsat->batm", ws, Mt, np.where(fin, zs, 0.0))
    den = np.einsum("bs,stm,bsat->batm", ws, Mt, fin.astype(float))
    with np.errstate(invalid="ignore", divide="ignore"):
        R = np.tanh(num / den)
    out["R"] = R
    lt2 = np.log2(np.asarray(T_grid, dtype=float))
    mids = [int(np.flatnonzero(lt2 == np.floor((np.log2(lo) + np.log2(hi)) / 2))[0]) for lo, hi in ranges]
    ia0, ih = int(np.flatnonzero(a_grid == 0)[0]), int(np.flatnonzero(a_grid == 0.5)[0])
    with np.errstate(invalid="ignore", divide="ignore"):
        out["Rratio"] = np.stack([R[:, ih, t, :] / R[:, ia0, t, :] for t in mids], axis=1)
    qa = a_grid <= min(1.5, a_max)
    X = np.stack([np.ones(qa.sum()), a_grid[qa], a_grid[qa] ** 2], axis=1)
    coef = np.einsum("ca,batm->bctm", np.linalg.pinv(X), R[:, qa])
    with np.errstate(invalid="ignore", divide="ignore"):
        out["a_star"] = np.where(coef[:, 2] < 0, -coef[:, 1] / (2 * coef[:, 2]), np.nan)
    sel = a_grid <= a_max
    out["pf"] = np.einsum("barm,a->brm", out["zeta_K"][:, sel], _ols_weights(a_grid[sel]))
    h = a_grid <= A_RATIO
    wh = _ols_weights(a_grid[h])
    out["vk_diff"] = (np.einsum("barm,a->brm", out["zeta_V"][:, h], wh)
                      - 2 * np.einsum("barm,a->brm", out["zeta_K"][:, h], wh))  # descriptive only
    n_range = np.array([np.sum((T_grid >= lo) & (T_grid <= hi)) for lo, hi in ranges])
    out["reach"] = range_reach(Ns, specs, ranges)
    rule = range_rule(out["reach"][None], out["kept_K"] == n_range[None, None, :, None], a_grid, a_max,
                      min_reach, full_range)
    for key, ok in rule.items():
        out[key] = np.where(ok, out[key], np.nan)
    return out


def estimate2(raw: Raw, specs, norm: str = "tilde", **kw) -> dict:
    S_, D_ = raw.Nsd.shape
    o = core2(raw, np.ones((1, D_)), np.ones((1, S_)), specs, norm, **kw)
    return {key: val[0] for key, val in o.items()}


def jackknife2(raw: Raw, specs, unit: str = "day", norm: str = "tilde", **kw) -> tuple[dict, int]:
    """Delete-one-day (primary) or delete-one-stock (diagnostic) jackknife SE of ESTIMANDS2."""
    S_, D_ = raw.Nsd.shape
    if unit == "day":
        Wd, Ws = jackknife_weights(S_, D_)
    elif unit == "stock":
        Wd, Ws = np.ones((S_, D_)), 1.0 - np.eye(S_)
    else:
        raise ValueError("unit is 'day' or 'stock'")
    n = Wd.shape[0]
    o = core2(raw, Wd, Ws, specs, norm, **kw)
    out = {}
    for key in ESTIMANDS2:
        loo = o[key]
        se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean(axis=0)) ** 2, axis=0))
        se[~np.all(np.isfinite(loo), axis=0)] = np.nan
        out[key] = se
    return out, n


def intervals2(raw: Raw, specs, rng: np.random.Generator, norm: str = "tilde", n_boot: int = N_BOOT,
               level: float = LEVEL, chunk: int = 50, stock_jk: bool = False, **kw) -> dict:
    """{estimand: (estimate, lo, hi, method)}: the wider of the percentile bootstrap (B = n_boot) and the
    day jackknife with Student t (G - 1 df); method 1 = bootstrap, 2 = jackknife, 0 = none. With stock_jk,
    also "stock_jk": {estimand: (lo, hi)}, the delete-one-stock jackknife t interval (S - 1 df), a
    diagnostic only. "dropped_K" carries the point estimate's dropped horizons."""
    S_, D_ = raw.Nsd.shape
    e1 = core2(raw, np.ones((1, D_)), np.ones((1, S_)), specs, norm, **kw)
    est = {key: val[0] for key, val in e1.items()}
    Wd, Ws = bootstrap_weights(S_, D_, n_boot, rng)
    boot = {key: [] for key in ESTIMANDS2}
    for i in range(0, n_boot, chunk):
        o = core2(raw, Wd[i : i + chunk], Ws[i : i + chunk], specs, norm, **kw)
        for key in ESTIMANDS2:
            boot[key].append(o[key])
    se, n = jackknife2(raw, specs, "day", norm, **kw)
    out = {}
    for key in ESTIMANDS2:
        b = np.concatenate(boot[key])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            blo = np.nanpercentile(b, 100 * (1 - level) / 2, axis=0)
            bhi = np.nanpercentile(b, 100 * (1 + level) / 2, axis=0)
        jlo, jhi = jk_interval(est[key], se[key], n, level)
        use_j = ~((bhi - blo) >= (jhi - jlo))
        lo, hi = np.where(use_j, jlo, blo), np.where(use_j, jhi, bhi)
        method = np.where(use_j, 2, 1)
        undefined = ~np.isfinite(est[key]) | ~np.isfinite(se[key]) | ~np.isfinite(blo)
        lo[undefined], hi[undefined], method[undefined] = np.nan, np.nan, 0
        out[key] = (est[key], lo, hi, method)
    if stock_jk:
        sse, ns = jackknife2(raw, specs, "stock", norm, **kw)
        out["stock_jk"] = {key: jk_interval(est[key], sse[key], ns, level) for key in ESTIMANDS2}
    if "dropped_K" in est:
        out["dropped_K"] = est["dropped_K"]
    return out


# ---- Phase 3c-4: the widest of three intervals ----------------------------------------------------

METHODS3 = ("day_stock_bootstrap", "day_jackknife_t", "whole_stock_bootstrap")


def whole_stock_weights(n_s: int, n_d: int, B: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Whole-stock bootstrap: stocks drawn with replacement, each with ALL its days (day weights 1)."""
    Ws = np.zeros((B, n_s))
    for b in range(B):
        Ws[b] = np.bincount(rng.integers(0, n_s, n_s), minlength=n_s)
    return np.ones((B, n_d)), Ws


def _percentile(reps: np.ndarray, level: float):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lo = np.nanpercentile(reps, 100 * (1 - level) / 2, axis=0)
        return lo, np.nanpercentile(reps, 100 * (1 + level) / 2, axis=0)


def intervals3(raw: Raw, specs, rng: np.random.Generator, n_boot: int = N_BOOT, level: float = LEVEL,
               chunk: int = 50, keys=ESTIMANDS2, **kw) -> dict:
    """{estimand: (estimate, lo, hi, method)}: the WIDEST of (i) the day-and-stock cluster bootstrap
    (inference.cluster_bootstrap's draws), (ii) the delete-one-day jackknife with Student t (G - 1 df),
    (iii) the whole-stock bootstrap (stocks with all their days), B = n_boot each, percentile intervals for
    the bootstraps. method = 1, 2, 3 for (i), (ii), (iii), 0 when the estimate is undefined. A method that is
    undefined for an element is skipped for it; no method defined: no interval."""
    S_, D_ = raw.Nsd.shape
    e1 = core2(raw, np.ones((1, D_)), np.ones((1, S_)), specs, "tilde", **kw)
    est = {key: e1[key][0] for key in keys}
    reps = {}
    draws = (("ds", bootstrap_weights(S_, D_, n_boot, rng)), ("ws", whole_stock_weights(S_, D_, n_boot, rng)))
    for name, (Wd, Ws) in draws:
        acc = {key: [] for key in keys}
        for i in range(0, n_boot, chunk):
            o = core2(raw, Wd[i : i + chunk], Ws[i : i + chunk], specs, "tilde", **kw)
            for key in keys:
                acc[key].append(o[key])
        reps[name] = {key: np.concatenate(v) for key, v in acc.items()}
    se, n = jackknife2(raw, specs, "day", "tilde", **kw)
    out = {}
    for key in keys:
        cands = [_percentile(reps["ds"][key], level), jk_interval(est[key], se[key], n, level),
                 _percentile(reps["ws"][key], level)]
        widths = np.stack([np.where(np.isfinite(hi - lo), hi - lo, -np.inf) for lo, hi in cands])
        pick = np.argmax(widths, axis=0)  # ties to the first
        lo = np.choose(pick, [c[0] for c in cands])
        hi = np.choose(pick, [c[1] for c in cands])
        method = pick + 1
        undefined = ~np.isfinite(est[key]) | ~np.isfinite(widths.max(axis=0))
        lo, hi = np.where(undefined, np.nan, lo), np.where(undefined, np.nan, hi)
        out[key] = (est[key], lo, hi, np.where(undefined, 0, method))
    out["dropped_K"] = e1["dropped_K"][0]
    return out
