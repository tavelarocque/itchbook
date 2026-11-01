"""Phase 3d-1: the four primary tests P-V, P-R, P-M, P-F and Holm-Bonferroni, as frozen functions.

Constants fixed by the pre-registration live in prereg3d. Every decision taken here is marked CHOICE and
listed in docs/phase3d_tests.md. No function here reads data files: callers pass Flow objects.

Pooling unit (tick_bin_unit). The confirmatory cells of one tick bin are pooled. Each stock enters horizon T
only while T <= rho_cell N_s, through the mask ("cellceil", rho per stock), so through the matched sets of
pooled.core2 every stock contributes only up to its own cell's ceiling. Descriptive cells are dropped and
logged. Per-cell units (cell_unit) are secondary.

Range rule (select_range). [8, 64] if at least pooled.MIN_REACH stocks reach T2 = 64 under the ceilings,
else [4, 32] if at least MIN_REACH reach 32, else not evaluable. It is decided from day lengths and
ceilings alone, before any estimate. Every test reports both ranges and the one it used.

Inference. One set of draws per unit, shared by every test: the day-and-stock bootstrap
(pooled.bootstrap_weights), the delete-one-day jackknife with Student t (D - 1 df), and the whole-stock
bootstrap (pooled.whole_stock_weights), in that order from one generator, as pooled.intervals3. The interval
is the widest of the three at prereg3d.LEVEL (0.90). CHOICE: a p-value is the LARGEST of the three
methods' p-values, the conservative counterpart of the widest interval. Percentile bootstrap: the share of
replicates on the null side (two-sided: twice the smaller share). Jackknife: Student t at D - 1 df. An
undefined estimate has p = 1.

The tests (each returns a TestResult; Holm then sets rejected and passed):
  P-V  (i) leave-one-day-out model comparison on log v_a(T): M2 must beat M1 and M3. (ii) the slope in a,
       over a below M2's fitted a_c, of zeta_V(a) real minus surrogate is negative, for S1 and for S2.
       p = max(p_S1, p_S2) when (i) holds, else 1 (CHOICE: (i) is a gate without a p-value).
  P-R  equivalence: a*(T1) and a*(T2) of the range used both within 1/2 +- 0.15 (TOST at each end;
       intersection-union, p = the larger end's TOST p). Secondary: R_{1/2}/R_0 against e^{sigma_l^2/8},
       psi-hat, and da*/dlog T with a two-sided interval.
  P-M  size-binned trade-level propagator: the slope of beta-hat on log q over bins is negative (one-sided)
       AND the top bin's beta-hat has its upper bound below 0. p = max of the two one-sided p-values.
  P-F  the slope of zeta_K(a) on a over [0, 1]: two-sided against 0. zeta_K levels reported, not tested.
A test passes when Holm rejects its p-value AND every sub-condition holds.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

import clip3b as CL
import inference as INF
import pooled as PO
import prereg3d as C
import propagator as PR
import surrogates as SU

KW = dict(T_grid=C.T_GRID, a_grid=C.A_GRID, ranges=C.RANGES, a_max=C.A_MAX)
RANGE_NAMES = tuple(f"[{lo},{hi}]" for lo, hi in C.RANGES)
CHUNK = 50
N_SUR = 4  # CHOICE: surrogate draws per surrogate kind, zeta_V averaged over them
PSI_BINS = 10  # CHOICE: deciles of run size for psi-hat
PSI_GRID = np.round(np.arange(0.0, 1.5001, 0.01), 2)
PM_L_MIN = 1  # CHOICE: the propagator fit starts at lag 1
PM_BETAS = np.round(np.arange(-1.0, 1.5001, 0.01), 2)  # CHOICE: propagator.fit_power's form, beta < 0 allowed
PM_L0S = np.logspace(-1, 2.5, 71)  # propagator.fit_power's grid
M2_KBELOW = (2, 3, 4, 5)  # CHOICE: a_c leaves 2..5 grid points below it (a_c in (0.25, 0.5], ..., or > 1)
SNAP = 1e-12  # CHOICE: a slope below this in size is 0 (an exactly flat input gives +-3e-17 from the OLS weights)


# ---- input ------------------------------------------------------------------------------------------------


@dataclass
class Flow:
    """Stock-days of some stocks, stock-major, concatenated runs, every stock with n_days study days.
    q is the RAW run size; the R2 clip is applied in make_unit."""

    eps: np.ndarray
    q: np.ndarray
    lm: np.ndarray
    ok: np.ndarray
    ts: np.ndarray  # ns since midnight ET, for the S2 surrogate's 30-minute bins
    starts: np.ndarray
    lengths: np.ndarray
    n_days: int

    @property
    def Nsd(self) -> np.ndarray:
        return self.lengths.reshape(-1, self.n_days)

    @property
    def n_stocks(self) -> int:
        return len(self.lengths) // self.n_days

    def day(self, i: int) -> slice:
        return slice(int(self.starts[i]), int(self.starts[i] + self.lengths[i]))


def concat_flows(flows: list[Flow]) -> Flow:
    if len({f.n_days for f in flows}) != 1:
        raise ValueError("every stock needs the same number of study days")
    lengths = np.concatenate([f.lengths for f in flows]).astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int64)
    cat = lambda k: np.concatenate([np.concatenate([getattr(f, k)[f.day(i)] for i in range(len(f.lengths))])  # noqa: E731
                                    for f in flows])
    return Flow(cat("eps").astype(float), cat("q").astype(float), cat("lm").astype(float), cat("ok").astype(bool),
                cat("ts").astype(np.int64), starts, lengths, flows[0].n_days)


def r2_clip(q: np.ndarray, lengths: np.ndarray, n_days: int) -> np.ndarray:
    """Each stock's sizes clipped at its own 99.9th percentile over its study days (clip3b R2)."""
    qs = np.split(np.asarray(q, dtype=float), np.cumsum(lengths)[:-1])
    stocks = np.repeat(np.arange(len(lengths) // n_days), n_days)
    caps = CL.stock_caps(qs, stocks, C.CLIP_QUANTILE)
    return np.minimum(q, np.repeat(np.array([caps[s] for s in stocks]), lengths))


@dataclass
class Unit:
    """One pooling unit: R2-clipped flow, per-stock ceilings and cells, the panel moments."""

    name: str
    flow: Flow  # q clipped
    rho: np.ndarray  # (S,)
    cells: np.ndarray  # (S,)
    excluded_cells: tuple
    raw: PO.Raw

    @property
    def specs(self):
        return [("cellceil", self.rho)]

    @property
    def Ns(self) -> np.ndarray:
        return self.flow.Nsd.min(axis=1)


def moments_raw(flow: Flow, q: np.ndarray | None = None) -> PO.Raw:
    pm = PO.panel_moments(flow.eps, flow.q if q is None else q, flow.lm, flow.ok, flow.starts, flow.lengths,
                          C.T_GRID, C.A_GRID)
    return PO.raw_panel(pm, flow.n_days)


def make_unit(name: str, flows: dict, cells: list) -> Unit:
    """Pool the given confirmatory cells; each stock gets its cell's ceiling from prereg3d."""
    for c in cells:
        if c not in C.CONFIRMATORY:
            raise ValueError(f"{c} is not a confirmatory cell")
    fl = concat_flows([flows[c] for c in cells])
    fl = Flow(fl.eps, r2_clip(fl.q, fl.lengths, fl.n_days), fl.lm, fl.ok, fl.ts, fl.starts, fl.lengths, fl.n_days)
    rho = np.concatenate([np.full(flows[c].n_stocks, C.CONFIRMATORY[c]) for c in cells])
    cl = np.concatenate([np.full(flows[c].n_stocks, c, dtype=object) for c in cells])
    excluded = tuple(sorted(c for c in flows if c in C.DESCRIPTIVE_CELLS and C.tick_bin(c) == C.tick_bin(cells[0])))
    return Unit(name, fl, rho, cl, excluded, moments_raw(fl))


def tick_bin_unit(flows: dict, tb: str) -> Unit:
    """The primary pooling unit: the tick bin's confirmatory cells; its descriptive cells are dropped."""
    cells = sorted(c for c in flows if C.tick_bin(c) == tb and c in C.CONFIRMATORY)
    if not cells:
        raise ValueError(f"no confirmatory cell of {tb} in the input")
    return make_unit(tb, flows, cells)


def cell_unit(flows: dict, cell: str) -> Unit:
    """Secondary: one confirmatory cell alone."""
    return make_unit(cell, flows, [cell])


# ---- range rule -------------------------------------------------------------------------------------------


def select_range(Ns: np.ndarray, rho: np.ndarray) -> tuple[tuple | None, dict]:
    """[8, 64] if at least MIN_REACH stocks reach 64 under their ceilings, else [4, 32] if MIN_REACH reach
    32, else None. Uses day lengths and ceilings only, never an estimate."""
    reach = PO.range_reach(np.asarray(Ns), [("cellceil", np.asarray(rho))], C.RANGES)[:, 0]
    by = {r: int(n) for r, n in zip(C.RANGES, reach, strict=True)}
    if by[C.PRIMARY_RANGE] >= PO.MIN_REACH:
        return C.PRIMARY_RANGE, by
    if by[C.FALLBACK_RANGE] >= PO.MIN_REACH:
        return C.FALLBACK_RANGE, by
    return None, by


def range_label(rng_used) -> str:
    if rng_used is None:
        return "not evaluable (fewer than MIN_REACH stocks at T2 of either range)"
    return "primary [8,64]" if rng_used == C.PRIMARY_RANGE else "fallback [4,32] ([8,64] below the reach)"


def stocks_at_top(Ns, rho) -> dict:
    """Indices of the stocks eligible at each range's T2."""
    M = PO.masks(np.asarray(Ns), np.array([hi for _, hi in C.RANGES]), [("cellceil", np.asarray(rho))])[:, :, 0]
    return {r: np.flatnonzero(M[:, i]).tolist() for i, r in enumerate(C.RANGES)}


# ---- replicates, intervals, p-values ------------------------------------------------------------------------


@dataclass
class Draws:
    full: tuple
    ds: tuple
    jk: tuple
    ws: tuple


def draws(n_s: int, n_d: int, rng: np.random.Generator, n_boot: int = C.N_BOOT) -> Draws:
    ds = PO.bootstrap_weights(n_s, n_d, n_boot, rng)
    ws = PO.whole_stock_weights(n_s, n_d, n_boot, rng)
    return Draws((np.ones((1, n_d)), np.ones((1, n_s))), ds, PO.jackknife_weights(n_s, n_d), ws)


@dataclass
class Reps:
    est: np.ndarray  # (...)
    ds: np.ndarray  # (B, ...)
    jk: np.ndarray  # (D, ...) leave-one-day-out values
    ws: np.ndarray  # (B, ...)

    def __getitem__(self, idx) -> Reps:
        sl = idx if isinstance(idx, tuple) else (idx,)
        return Reps(np.asarray(self.est)[sl], self.ds[(slice(None), *sl)], self.jk[(slice(None), *sl)],
                    self.ws[(slice(None), *sl)])


def run_stats(stat, dr: Draws, chunk: int = CHUNK) -> dict:
    """stat(Wd, Ws) -> {name: (B, ...)}; returns {name: Reps}."""
    def go(W):
        Wd, Ws = W
        parts = [stat(Wd[i : i + chunk], Ws[i : i + chunk]) for i in range(0, len(Wd), chunk)]
        return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}

    full, ds, jk, ws = go(dr.full), go(dr.ds), go(dr.jk), go(dr.ws)
    return {k: Reps(full[k][0], ds[k], jk[k], ws[k]) for k in full}


def _jk_se(loo: np.ndarray) -> tuple[np.ndarray, int]:
    n = loo.shape[0]
    se = np.sqrt((n - 1) / n * np.sum((loo - loo.mean(axis=0)) ** 2, axis=0))
    return np.where(np.all(np.isfinite(loo), axis=0), se, np.nan), n


def _pct(reps: np.ndarray, level: float):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return (np.nanpercentile(reps, 100 * (1 - level) / 2, axis=0),
                np.nanpercentile(reps, 100 * (1 + level) / 2, axis=0))


def interval(r: Reps, level: float = C.LEVEL) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(lo, hi, method): the widest of day-and-stock bootstrap (1), day jackknife t (2), whole-stock
    bootstrap (3); method 0 and NaN bounds where the estimate or every method is undefined."""
    est = np.asarray(r.est, dtype=float)
    se, n = _jk_se(r.jk)
    c = INF.t_ppf((1 + level) / 2, n - 1)
    cands = [_pct(r.ds, level), (est - c * se, est + c * se), _pct(r.ws, level)]
    widths = np.stack([np.where(np.isfinite(hi - lo), hi - lo, -np.inf) for lo, hi in cands])
    pick = np.argmax(widths, axis=0)
    lo = np.choose(pick, [x[0] for x in cands])
    hi = np.choose(pick, [x[1] for x in cands])
    bad = ~np.isfinite(est) | ~np.isfinite(widths.max(axis=0))
    return np.where(bad, np.nan, lo), np.where(bad, np.nan, hi), np.where(bad, 0, pick + 1)


def _boot_p(reps: np.ndarray, null, alternative: str) -> np.ndarray:
    fin = np.isfinite(reps)
    n = fin.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ge = np.sum(fin & (reps >= null), axis=0) / n
        le = np.sum(fin & (reps <= null), axis=0) / n
    p = {"less": ge, "greater": le, "two-sided": np.minimum(1.0, 2 * np.minimum(ge, le))}[alternative]
    return np.where(n > 0, p, np.nan)


def _jk_p(est, loo, null, alternative: str) -> np.ndarray:
    se, n = _jk_se(loo)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = (np.asarray(est, dtype=float) - null) / se
    t = np.where((se == 0) & (np.asarray(est) == null), np.nan, t)
    cdf = np.vectorize(lambda x: INF.t_cdf(x, n - 1) if np.isfinite(x) else (1.0 if x > 0 else 0.0))
    tt = np.where(np.isnan(t), 0.0, t)
    F = cdf(tt).astype(float)
    if alternative == "less":
        p = F
    elif alternative == "greater":
        p = 1 - F
    else:
        p = np.minimum(1.0, 2 * np.minimum(F, 1 - F))
    p = np.where((se == 0) & (np.asarray(est) == null), 1.0, p)
    return np.where(np.isfinite(se), p, np.nan)


def p_value(r: Reps, null, alternative: str) -> np.ndarray:
    """The largest of the three methods' p-values (methods undefined for the element skipped); 1 where the
    estimate is undefined or no method is defined. alternative: 'less', 'greater' or 'two-sided'."""
    ps = np.stack([_boot_p(r.ds, null, alternative), _jk_p(r.est, r.jk, null, alternative),
                   _boot_p(r.ws, null, alternative)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        p = np.nanmax(ps, axis=0)
    return np.where(np.isfinite(r.est) & np.isfinite(p), p, 1.0)


def summary(r: Reps, level: float = C.LEVEL) -> dict:
    lo, hi, m = interval(r, level)
    return {"estimate": np.asarray(r.est).tolist(), "lo": lo.tolist(), "hi": hi.tolist(),
            "method": [PO.METHODS3[i - 1] if i else None for i in np.ravel(m)]}


# ---- results and Holm -------------------------------------------------------------------------------------


@dataclass
class TestResult:
    name: str
    p: float
    alternative: str
    conditions: dict  # sub-condition name -> bool
    statistic: dict  # the primary range's statistics and intervals
    range_used: tuple | None
    range_rule: str
    reach: dict
    stocks_used: dict
    by_range: dict = field(default_factory=dict)  # both ranges, reported
    secondary: dict = field(default_factory=dict)  # reported only, never in Holm
    rejected: bool | None = None  # set by apply_holm
    passed: bool | None = None

    @property
    def all_conditions(self) -> bool:
        return bool(self.conditions) and all(self.conditions.values())


def holm(pvals: dict, alpha: float = C.ALPHA_FW) -> dict:
    """Holm-Bonferroni over EXACTLY the four primary tests. Ascending p; the i-th smallest (i = 0..3) is
    compared with alpha / (4 - i); stop at the first that fails. A NaN p counts as 1."""
    if len(pvals) != len(C.PRIMARY_TESTS) or set(pvals) != set(C.PRIMARY_TESTS):
        raise ValueError(f"Holm takes exactly {C.PRIMARY_TESTS}, got {sorted(pvals)}")
    names = list(C.PRIMARY_TESTS)
    p = np.array([1.0 if not np.isfinite(pvals[k]) else float(pvals[k]) for k in names])
    m = len(p)
    order = np.argsort(p, kind="stable")
    out = {k: {"p": float(p[i]), "rank": None, "threshold": None, "rejected": False, "p_adj": None}
           for i, k in enumerate(names)}
    running, stop = 0.0, False
    for i, j in enumerate(order):
        thr = alpha / (m - i)
        running = max(running, min(1.0, (m - i) * p[j]))
        rec = out[names[j]]
        rec.update(rank=i + 1, threshold=thr, p_adj=running)
        if not stop and p[j] <= thr:
            rec["rejected"] = True
        else:
            stop = True
    return out


def apply_holm(results: dict, alpha: float = C.ALPHA_FW) -> dict:
    h = holm({k: results[k].p for k in C.PRIMARY_TESTS}, alpha)
    for k in C.PRIMARY_TESTS:
        results[k].rejected = h[k]["rejected"]
        results[k].passed = bool(h[k]["rejected"] and results[k].all_conditions)
    return h


# ---- decision layers (pure functions of replicates) ----------------------------------------------------------


def excludes_zero_below(r: Reps) -> bool:
    return bool(interval(r)[1] < 0)


def decide_pf(pf: Reps) -> dict:
    lo, hi, _ = interval(pf)
    return {"p": float(p_value(pf, 0.0, "two-sided")), "alternative": "two-sided",
            "conditions": {"interval_excludes_0": bool((lo > 0) | (hi < 0))}}


def decide_pr(ends: Reps) -> dict:
    """ends: a*(T1), a*(T2). TOST at each end against 1/2 -+ 0.15; p = the larger end's TOST p."""
    lo_m, hi_m = C.PR_CENTRE - C.PR_MARGIN, C.PR_CENTRE + C.PR_MARGIN
    p_lo = p_value(ends, lo_m, "greater")  # H0: a* <= 0.35
    p_hi = p_value(ends, hi_m, "less")  # H0: a* >= 0.65
    p_end = np.maximum(p_lo, p_hi)
    lo, hi, _ = interval(ends)
    inside = (lo > lo_m) & (hi < hi_m)
    return {"p": float(np.max(p_end)), "alternative": "equivalence (TOST, both ends)",
            "p_end": p_end.tolist(), "conditions": {"T1_within_margin": bool(inside[0]),
                                                     "T2_within_margin": bool(inside[1])}}


def decide_pm(slope: Reps, top: Reps) -> dict:
    p_s, p_t = float(p_value(slope, 0.0, "less")), float(p_value(top, 0.0, "less"))
    return {"p": max(p_s, p_t), "alternative": "less (both)", "p_slope": p_s, "p_top": p_t,
            "conditions": {"slope_declines_with_log_q": excludes_zero_below(slope),
                           "top_bin_upper_below_0": excludes_zero_below(top)}}


def decide_pv(errors: dict, s1: Reps, s2: Reps) -> dict:
    """errors: {'M1', 'M2', 'M3'} leave-one-day-out totals (NaN when undefined)."""
    e1, e2, e3 = (errors.get(k, np.nan) for k in ("M1", "M2", "M3"))
    beats1, beats3 = bool(e2 < e1), bool(e2 < e3)
    p1, p2 = float(p_value(s1, 0.0, "less")), float(p_value(s2, 0.0, "less"))
    gate = beats1 and beats3
    return {"p": max(p1, p2) if gate else 1.0, "alternative": "less (S1 and S2), gated by (i)",
            "p_S1": p1, "p_S2": p2,
            "conditions": {"i_M2_beats_M1": beats1, "i_M2_beats_M3": beats3,
                           "ii_S1_negative": excludes_zero_below(s1), "ii_S2_negative": excludes_zero_below(s2)}}


# ---- shared statistics -------------------------------------------------------------------------------------


def _ols(x) -> np.ndarray:
    return PO._ols_weights(np.asarray(x, dtype=float))


def snap(x: np.ndarray) -> np.ndarray:
    """Slopes within SNAP of 0 are 0, so rounding never decides a test."""
    return np.where(np.abs(x) < SNAP, 0.0, x)


def _range_idx(r) -> np.ndarray:
    lo, hi = r
    return np.flatnonzero((C.T_GRID >= lo) & (C.T_GRID <= hi))


def vertex(R: np.ndarray, a_grid=C.A_GRID) -> np.ndarray:
    """a* = vertex of the least-squares quadratic in a of R (..., na); NaN unless concave (as pooled.core2)."""
    a = np.asarray(a_grid, dtype=float)
    X = np.stack([np.ones_like(a), a, a * a], axis=1)
    coef = np.einsum("ca,...a->...c", np.linalg.pinv(X), R)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(coef[..., 2] < 0, -coef[..., 1] / (2 * coef[..., 2]), np.nan)


def _pooled(x_sd: np.ndarray, Wd, Ws, n_days: int) -> np.ndarray:
    """Equal weight per stock of each stock's day average (finite days only), under replicate weights:
    x_sd (n_sd, ...) -> (B, ...). Pooling coefficients are built once per stock-day."""
    x = x_sd.reshape(-1, n_days, *x_sd.shape[1:])
    fin = np.isfinite(x).reshape(*x.shape[:2], -1).all(axis=2)  # a stock-day is used whole or not at all
    wd = Wd[:, None, :] * fin[None]  # (B, S, D)
    nd = wd.sum(axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        sw = np.where(nd > 0, Ws / np.where(nd > 0, nd, 1.0), 0.0)
        coef = wd * sw[:, :, None]
        coef = coef / (Ws * (nd > 0)).sum(axis=1)[:, None, None]
    return np.einsum("bsd,sd...->b...", coef, np.where(fin.reshape(*fin.shape, *[1] * (x.ndim - 2)), x, 0.0))


def sig2_per_sd(flow: Flow) -> np.ndarray:
    """Variance of log run size (R2-clipped) per stock-day. CHOICE: the sizes the moments see."""
    return np.array([np.var(np.log(flow.q[flow.day(i)])) for i in range(len(flow.lengths))])


def psi_moments(flow: Flow, n_bins: int = PSI_BINS, grid=PSI_GRID) -> tuple[np.ndarray, np.ndarray]:
    """Per stock-day, over deciles of run size: y = mean eps_n (lm_{n+1} - lm_n) and X_g = mean (q / gm)^g,
    gm the stock-day's geometric mean size; only steps with ok at both ends."""
    y = np.full((len(flow.lengths), n_bins), np.nan)
    X = np.full((len(flow.lengths), n_bins, len(grid)), np.nan)
    for i in range(len(flow.lengths)):
        s = flow.day(i)
        e, q, lm, ok = flow.eps[s][:-1], flow.q[s][:-1], flow.lm[s], flow.ok[s]
        v = ok[:-1] & ok[1:]
        if v.sum() < 2 * n_bins:
            continue
        r = np.diff(lm)
        e, q, r = e[v], q[v], r[v]
        b = np.searchsorted(np.quantile(q, np.linspace(0, 1, n_bins + 1)[1:-1]), q, side="left")
        gm = np.exp(np.mean(np.log(q)))
        cnt = np.bincount(b, minlength=n_bins).astype(float)
        if np.any(cnt == 0):
            continue
        y[i] = np.bincount(b, e * r, minlength=n_bins) / cnt
        P = (q / gm)[:, None] ** np.asarray(grid)[None, :]
        X[i] = np.stack([np.bincount(b, P[:, g], minlength=n_bins) for g in range(len(grid))], axis=1) / cnt[:, None]
    return y, X


def psi_hat(y_sd, X_sd, Wd, Ws, n_days: int, grid=PSI_GRID) -> np.ndarray:
    """CHOICE (psi-hat): common psi over stocks minimising the sum over stocks of 1 - R^2 of the stock's
    day-averaged binned response on c0 + c1 X_psi (stock-specific c0, c1). (B,)"""
    S_ = len(y_sd) // n_days
    yb = PO._stock_mean(y_sd.reshape(S_, n_days, -1), Wd[:, None, :] * np.ones((S_, n_days))[None])
    Xb = PO._stock_mean(X_sd.reshape(S_, n_days, *X_sd.shape[1:]), Wd[:, None, :] * np.ones((S_, n_days))[None])
    yc = yb - yb.mean(axis=2, keepdims=True)  # (B, S, nb)
    Xc = Xb - Xb.mean(axis=2, keepdims=True)  # (B, S, nb, ng)
    with np.errstate(invalid="ignore", divide="ignore"):
        sxy = np.einsum("bsk,bskg->bsg", yc, Xc)
        sxx = np.einsum("bskg,bskg->bsg", Xc, Xc)
        tss = np.einsum("bsk,bsk->bs", yc, yc)
        loss = np.where(sxx > 0, (tss[..., None] - sxy**2 / sxx) / tss[..., None], 1.0)  # 1 - R^2; psi = 0: flat
    fin = np.isfinite(loss).all(axis=2) & (tss > 0)
    w = Ws * fin
    score = np.einsum("bs,bsg->bg", w, np.where(fin[..., None], loss, 0.0)) / w.sum(axis=1, keepdims=True)
    out = np.asarray(grid)[np.argmin(np.where(np.isfinite(score), score, np.inf), axis=1)]
    return np.where(np.isfinite(score).any(axis=1), out, np.nan)


# ---- P-V (i): the three models ---------------------------------------------------------------------------------


def m3_V(flow: Flow) -> np.ndarray:
    """M3, the independent-volume prediction per stock-day, (n_sd, na, nT):
    V_a(T) = T E[q^2a] + 2 E[q^a]^2 sum_{l=1}^{T-1} (T - l) C(l), C the stock-day's uncentred sign ACF,
    the moments over its (clipped) sizes. It is the expectation of V_a(T) under within-day size shuffling."""
    T = C.T_GRID
    out = np.empty((len(flow.lengths), len(C.A_GRID), len(T)))
    for i in range(len(flow.lengths)):
        s = flow.day(i)
        Cl = PR.sign_acf(flow.eps[s], int(T[-1]) - 1)
        tri = np.array([np.sum((t - np.arange(1, t)) * Cl[1:t]) for t in T])
        for j, a in enumerate(C.A_GRID):
            qa = flow.q[s] ** a
            out[i, j] = T * np.mean(qa * qa) + 2 * np.mean(qa) ** 2 * tri
    return out


def logv_rel(V: np.ndarray, Ns, specs, Wd, Ws, r) -> np.ndarray:
    """Pooled log v_a(T) - log v_a(T1) over the range's horizons (B, na, nT_r): the chain of matched
    increments (each stock's own log change, mean over stocks eligible at the upper horizon), as
    pooled.zeta_mean_of_logs, from V (S, D, na, nT)."""
    wd = Wd[:, None, :] * np.ones(V.shape[:2])[None]
    Vb = PO._stock_mean(V, wd)
    T = C.T_GRID.astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        v = Vb / (T * Vb[..., :1])
        lv = np.where(v > 0, np.log(np.where(v > 0, v, 1.0)), np.nan)
    dlt, _ = PO._increments(lv, Ws * 1.0, PO.increment_masks(Ns, C.T_GRID, specs))  # (B, na, nJ, 1)
    idx = _range_idx(r)
    inc = dlt[:, :, idx[0] : idx[-1], 0]
    return np.concatenate([np.zeros(inc.shape[:2] + (1,)), np.cumsum(inc, axis=2)], axis=2)


def fit_m1(s: np.ndarray) -> dict:
    """M1: zeta_V(a) - 1 = c for every a; s (na,) per-a slopes. Least squares: c = mean s."""
    c = float(np.mean(s))
    return {"zeta_minus_1": np.full(len(s), c), "zeta": 1 + c}


def fit_m2(s: np.ndarray, a=C.A_GRID) -> dict:
    """M2: zeta_V(a) = 3 - mu_m - 2 a lambda sigma_l^2 below a_c, 1 above; free mu_m, lambda sigma_l^2 >= 0
    (CHOICE: the model's sign), a_c with 2..5 grid points below it (CHOICE). Least squares on the per-a
    slopes s; ties go to the smaller a_c."""
    a = np.asarray(a, dtype=float)
    best = None
    for k in M2_KBELOW:
        x, y = a[:k], s[:k]
        b, c0 = np.polyfit(x, y, 1)
        if b > 0:  # lambda sigma_l^2 < 0 is not M2
            b, c0 = 0.0, float(np.mean(y))
        z = np.concatenate([c0 + b * x, np.zeros(len(a) - k)])
        rss = float(np.sum((s - z) ** 2))
        if best is None or rss < best["rss"] - 1e-15:
            best = {"rss": rss, "k_below": k, "zeta_minus_1": z, "mu_m": 2 - c0, "lambda_sigma2": -b / 2,
                    "a_c": (float(a[k - 1]), float(a[k]) if k < len(a) else np.inf)}
    return best


def _slopes_a(D: np.ndarray, x: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted LS slope through the origin of D (na, nT_r) on x (nT_r,) per a."""
    return (D * (w * x)[None]).sum(axis=1) / np.sum(w * x * x)


def range_wx(r) -> tuple[np.ndarray, np.ndarray]:
    """x = log T - log T1 and the weights w = 1/T (CHOICE: about the inverse variance of log v_a(T),
    which grows like T / N), over the range's horizons."""
    T = C.T_GRID[_range_idx(r)].astype(float)
    return np.log(T / T[0]), 1.0 / T


def model_errors(Dtr: np.ndarray, Dte: np.ndarray, D3: np.ndarray, r) -> dict:
    """Leave-one-day-out weighted squared errors of log v_a(T) (relative to T1) over the range.
    Dtr, Dte, D3: (n_days, na, nT_r): training curve without day d, day d's curve, day d's M3 curve.
    CHOICE: a day with any undefined value in any of the three is skipped for every model."""
    x, w = range_wx(r)
    E = {"M1": 0.0, "M2": 0.0, "M3": 0.0}
    used = 0
    for d in range(len(Dtr)):
        if not (np.all(np.isfinite(Dtr[d])) and np.all(np.isfinite(Dte[d])) and np.all(np.isfinite(D3[d]))):
            continue
        s = _slopes_a(Dtr[d], x, w)
        preds = {"M1": fit_m1(s)["zeta_minus_1"][:, None] * x[None],
                 "M2": fit_m2(s)["zeta_minus_1"][:, None] * x[None], "M3": D3[d]}
        for k, pr in preds.items():
            E[k] += float(np.sum(w[None] * (Dte[d] - pr) ** 2))
        used += 1
    if used == 0:
        return {"M1": np.nan, "M2": np.nan, "M3": np.nan, "days_used": 0}
    return {**E, "days_used": used}


# ---- P-M: size-binned propagator -------------------------------------------------------------------------------


@dataclass(frozen=True)
class PMParams:
    """Chosen in Phase 3d-2; validated here."""

    n_bins: int
    l_max: int
    L: int
    min_runs_per_bin: int

    def __post_init__(self):
        if self.n_bins not in C.PM_BINS_ALLOWED:
            raise ValueError(f"n_bins must be one of {C.PM_BINS_ALLOWED}")
        if self.l_max not in C.PM_LMAX_ALLOWED:
            raise ValueError(f"l_max must be one of {C.PM_LMAX_ALLOWED}")
        if self.L < C.PM_L_FACTOR * self.l_max:
            raise ValueError("L must be at least 8 l_max")
        if self.min_runs_per_bin < 1:
            raise ValueError("min_runs_per_bin must be positive")


def pm_moments(flow: Flow, rho: np.ndarray, pm: PMParams) -> dict:
    """Per stock-day, with x_b(n) = eps_n 1[q_n in size bin b] (bins: the stock-day's own size quantiles,
    ties to the lower bin) and r_n = lm_{n+1} - lm_n over valid steps, divided by the stock-day's rms r
    (CHOICE: comparable scale across stocks):
        Cp[b, b', k] = mean_n x_b(n) x_b'(n + k), k = 0..L-1;   Sr[b, l] = mean_n x_b(n) r_{n+l}, l = 0..L-1;
        xlog[b] = mean log q in bin b minus the stock-day's mean log q.
    A stock-day is used only if every bin holds >= min_runs_per_bin runs, N > L + 1, and its stock reaches
    L under its cell's ceiling (CHOICE: L <= rho N_s); unused stock-days are NaN."""
    nb, L = pm.n_bins, pm.L
    n_sd = len(flow.lengths)
    Ns = np.repeat(flow.Nsd.min(axis=1), flow.n_days)
    rho_sd = np.repeat(np.asarray(rho, dtype=float), flow.n_days)
    Cp = np.full((n_sd, nb, nb, L), np.nan)
    Sr = np.full((n_sd, nb, L), np.nan)
    xlog = np.full((n_sd, nb), np.nan)
    used = np.zeros(n_sd, dtype=bool)
    for i in range(n_sd):
        N = int(flow.lengths[i])
        if N <= L + 1 or L > rho_sd[i] * Ns[i]:
            continue
        s = flow.day(i)
        q, e, lm, ok = flow.q[s], flow.eps[s], flow.lm[s], flow.ok[s]
        b = np.searchsorted(np.quantile(q, np.linspace(0, 1, nb + 1)[1:-1]), q, side="left")
        cnt = np.bincount(b, minlength=nb)
        if cnt.min() < pm.min_runs_per_bin:
            continue
        v = ok[:-1] & ok[1:]
        r = np.where(v, np.diff(lm), 0.0)
        rms = np.sqrt(np.mean(r[v] ** 2)) if v.any() else 0.0
        if not rms > 0:
            continue
        x = np.zeros((nb, N))
        x[b, np.arange(N)] = e
        size = 1 << int(np.ceil(np.log2(2 * N)))
        F = np.fft.rfft(x, size, axis=1)
        cc = np.fft.irfft(np.conj(F)[:, None, :] * F[None, :, :], size, axis=-1)[..., :L]
        Cp[i] = cc / (N - np.arange(L))
        rr = np.zeros(N)
        rr[: N - 1] = r / rms
        sr = np.fft.irfft(np.conj(F) * np.fft.rfft(rr, size)[None], size, axis=-1)[:, :L]
        vv = np.zeros(N)
        vv[: N - 1] = v
        cnt_l = np.cumsum(vv[::-1])[::-1][:L]  # valid steps at m >= l
        Sr[i] = sr / cnt_l
        lq = np.log(q)
        xlog[i] = np.bincount(b, lq, minlength=nb) / cnt - lq.mean()
        used[i] = True
    return {"Cp": Cp, "Sr": Sr, "xlog": xlog, "used": used}


def pm_solve(Cp: np.ndarray, Sr: np.ndarray) -> np.ndarray:
    """Kernels of r_n = sum_b' sum_{j<L} K_b'(j) x_b'(n - j): solve S_b(l) = sum_b' sum_j K_b'(j)
    C_bb'(l - j), C_bb'(-k) = C_b'b(k). Returns G_b(l) = sum_{j<l} K_b(j), (B, nb, L + 1)."""
    B, nb, _, L = Cp.shape
    neg = np.swapaxes(Cp, 1, 2)[..., 1:][..., ::-1]
    Cf = np.concatenate([neg, Cp], axis=-1)  # index k + L - 1
    idx = np.arange(L)[:, None] - np.arange(L)[None, :] + L - 1
    M = Cf[:, :, :, idx].transpose(0, 1, 3, 2, 4).reshape(B, nb * L, nb * L)
    rhs = Sr.reshape(B, nb * L)
    K = np.full((B, nb * L), np.nan)
    for i in range(B):
        if np.all(np.isfinite(M[i])) and np.all(np.isfinite(rhs[i])):
            try:
                K[i] = np.linalg.solve(M[i], rhs[i])
            except np.linalg.LinAlgError:
                pass
    K = K.reshape(B, nb, L)
    return np.concatenate([np.zeros((B, nb, 1)), np.cumsum(K, axis=2)], axis=2)


def fit_beta(G: np.ndarray, l_min: int, l_max: int, betas=PM_BETAS, l0s=PM_L0S) -> np.ndarray:
    """beta of G(l) = G0 l0^beta / (l0 + l)^beta, least squares on propagator.fit_power's log-spaced lags,
    over the (beta, l0) grid with G0 in closed form, for every leading index of G (..., L + 1)."""
    lags = np.unique(np.round(np.logspace(np.log10(l_min), np.log10(l_max), 40)).astype(int))
    y = G[..., lags].reshape(-1, len(lags))
    f = (l0s[None, :, None] / (l0s[None, :, None] + lags[None, None, :])) ** np.asarray(betas)[:, None, None]
    ff = np.sum(f * f, axis=2)
    out = np.full(len(y), np.nan)
    ok = np.all(np.isfinite(y), axis=1)
    for i0 in range(0, len(y), 64):
        sel = np.flatnonzero(ok[i0 : i0 + 64]) + i0
        if len(sel) == 0:
            continue
        fy = np.einsum("gkn,rn->rgk", f, y[sel])
        rss = np.sum(y[sel] ** 2, axis=1)[:, None, None] - fy**2 / ff[None]
        out[sel] = np.asarray(betas)[np.argmin(rss.reshape(len(sel), -1), axis=1) // len(l0s)]
    return out.reshape(G.shape[:-1])


def pm_stats(mom: dict, pm: PMParams, n_days: int, Wd, Ws) -> dict:
    """beta per bin, relative log size per bin, slope of beta on log q over bins, top-bin beta: (B, ...)."""
    Cp = _pooled(mom["Cp"], Wd, Ws, n_days)
    Sr = _pooled(mom["Sr"], Wd, Ws, n_days)
    x = _pooled(mom["xlog"], Wd, Ws, n_days)
    beta = fit_beta(pm_solve(Cp, Sr), PM_L_MIN, pm.l_max)
    xc = x - x.mean(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        slope = np.sum(xc * (beta - beta.mean(axis=1, keepdims=True)), axis=1) / np.sum(xc * xc, axis=1)
    return {"beta": beta, "xlog": x, "slope": snap(slope), "top": beta[:, -1]}


# ---- the unit run --------------------------------------------------------------------------------------------


def surrogate_raws(u: Unit, kind: str, n: int, rng: np.random.Generator) -> list[PO.Raw]:
    """Panel moments of n surrogate draws of the unit's flow: q permuted within each stock-day (S1) or in
    blocks of L within 30-minute bins (S2, L from surrogates.choose_L per stock-day); signs, mids, ok fixed."""
    f = u.flow
    out = []
    for _ in range(n):
        q = f.q.copy()
        for i in range(len(f.lengths)):
            s = f.day(i)
            if kind == "S1":
                q[s] = SU.s1(f.q[s], rng)
            elif kind == "S2":
                q[s] = SU.s2(f.q[s], f.ts[s], SU.choose_L(f.q[s], f.ts[s]), rng)
            else:
                raise ValueError("kind is 'S1' or 'S2'")
        out.append(moments_raw(f, q))
    return out


def _empty(name: str, alternative: str, conditions: tuple, range_used, by) -> TestResult:
    return TestResult(name, 1.0, alternative, {c: False for c in conditions}, {}, range_used,
                      range_label(range_used), by, {})


def run_unit(u: Unit, pm: PMParams, rng: np.random.Generator, n_boot: int = C.N_BOOT, n_sur: int = N_SUR,
             chunk: int = CHUNK) -> dict:
    """The four primary tests on one pooling unit, with Holm over exactly those four. Returns
    {'tests': {name: TestResult}, 'holm': ..., 'descriptive': {...}, 'unit': {...}}."""
    S_, D_ = u.flow.Nsd.shape
    r_used, by = select_range(u.Ns, u.rho)
    top = stocks_at_top(u.Ns, u.rho)
    ir = {r: i for i, r in enumerate(C.RANGES)}
    dr = draws(S_, D_, rng, n_boot)

    # P-V (i): leave-one-day-out curves, and M2's a_c on all days for (ii)
    V3 = m3_V(u.flow).reshape(S_, D_, len(C.A_GRID), len(C.T_GRID))
    eye, ones_s = np.eye(D_), np.ones((D_, S_))
    pv_i, kbelow = {}, {}
    for r in C.RANGES:
        Dtr = logv_rel(u.raw.V, u.Ns, u.specs, 1.0 - eye, ones_s, r)
        Dte = logv_rel(u.raw.V, u.Ns, u.specs, eye, ones_s, r)
        D3 = logv_rel(V3, u.Ns, u.specs, eye, ones_s, r)
        pv_i[r] = model_errors(Dtr, Dte, D3, r) if by[r] >= PO.MIN_REACH else {"M1": np.nan, "M2": np.nan,
                                                                              "M3": np.nan, "days_used": 0}
        Dall = logv_rel(u.raw.V, u.Ns, u.specs, *dr.full, r)[0]
        x, w = range_wx(r)
        if np.all(np.isfinite(Dall)) and by[r] >= PO.MIN_REACH:
            m2 = fit_m2(_slopes_a(Dall, x, w))
            kbelow[r] = m2["k_below"]
            pv_i[r]["M2_fit_all_days"] = {k: m2[k] for k in ("mu_m", "lambda_sigma2", "a_c", "k_below")}
        else:
            kbelow[r] = None
    sur = {k: surrogate_raws(u, k, n_sur, rng) for k in ("S1", "S2")}
    sig2_sd = sig2_per_sd(u.flow)
    psi_y, psi_X = psi_moments(u.flow)
    i1 = {r: _range_idx(r) for r in C.RANGES}

    def stat(Wd, Ws):
        o = PO.core2(u.raw, Wd, Ws, u.specs, "tilde", **KW)
        reach_ok = o["reach"][:, 0] >= PO.MIN_REACH
        B = len(Wd)
        out = {"pf": snap(o["pf"][..., 0]), "zeta_K": o["zeta_K"][..., 0], "zeta_V": o["zeta_V"][..., 0],
               "vk_diff": o["vk_diff"][..., 0], "Rratio": o["Rratio"][..., 0]}
        ast = o["a_star"][..., 0]  # (B, nT)
        ends = np.full((B, len(C.RANGES), 2), np.nan)
        slope_a = np.full((B, len(C.RANGES)), np.nan)
        for r, j in ir.items():
            if reach_ok[j]:
                ends[:, j] = ast[:, [i1[r][0], i1[r][-1]]]
                slope_a[:, j] = ast[:, i1[r]] @ _ols(np.log(C.T_GRID[i1[r]]))
        out["a_star_ends"], out["dastar_dlogT"] = ends, snap(slope_a)
        s2 = _pooled(sig2_sd, Wd, Ws, D_)
        out["sig2"] = s2
        out["Rratio_minus_target"] = out["Rratio"] - np.exp(s2 / 8)[:, None]
        out["psi"] = psi_hat(psi_y, psi_X, Wd, Ws, D_)
        for kind, raws in sur.items():
            zs = np.mean([PO.core2(rw, Wd, Ws, u.specs, "tilde", **KW)["zeta_V"][..., 0] for rw in raws], axis=0)
            d = out["zeta_V"] - zs  # (B, na, nr)
            sl = np.full((B, len(C.RANGES)), np.nan)
            for r, j in ir.items():
                k = kbelow[r]
                if k is not None:
                    sl[:, j] = d[:, :k, j] @ _ols(C.A_GRID[:k])
            out[f"pv_slope_{kind}"] = snap(sl)
        return out

    reps = run_stats(stat, dr, chunk)
    mom = pm_moments(u.flow, u.rho, pm)
    pm_reps = run_stats(lambda Wd, Ws: pm_stats(mom, pm, D_, Wd, Ws), dr, max(1, chunk // 5))
    pm_stocks = np.flatnonzero(mom["used"].reshape(S_, D_).any(axis=1))

    tests = {}
    # P-F
    by_r = {}
    for r, j in ir.items():
        d = decide_pf(reps["pf"][j])
        by_r[r] = {**d, "pf": summary(reps["pf"][j]), "zeta_K_levels": summary(reps["zeta_K"][:, j])}
    tests["P-F"] = _assemble("P-F", by_r, r_used, by, top, ("interval_excludes_0",))
    # P-R
    by_r = {}
    for r, j in ir.items():
        d = decide_pr(reps["a_star_ends"][j])
        by_r[r] = {**d, "a_star_T1_T2": summary(reps["a_star_ends"][j])}
    tests["P-R"] = _assemble("P-R", by_r, r_used, by, top, ("T1_within_margin", "T2_within_margin"))
    tests["P-R"].secondary = {
        "Rratio": {RANGE_NAMES[j]: summary(reps["Rratio"][j]) for j in range(len(C.RANGES))},
        "Rratio_minus_exp_sig2_over_8": {RANGE_NAMES[j]: {**summary(reps["Rratio_minus_target"][j]),
                                                          "p_two_sided": float(p_value(reps["Rratio_minus_target"][j],
                                                                                       0.0, "two-sided"))}
                                         for j in range(len(C.RANGES))},
        "sigma_l2": summary(reps["sig2"]),
        "psi_hat": summary(reps["psi"]),
        "dastar_dlogT": {RANGE_NAMES[j]: {**summary(reps["dastar_dlogT"][j]),
                                          "p_two_sided": float(p_value(reps["dastar_dlogT"][j], 0.0, "two-sided"))}
                         for j in range(len(C.RANGES))},
    }
    # P-V
    by_r = {}
    for r, j in ir.items():
        d = decide_pv(pv_i[r], reps["pv_slope_S1"][j], reps["pv_slope_S2"][j])
        by_r[r] = {**d, "loo_errors": pv_i[r], "slope_S1": summary(reps["pv_slope_S1"][j]),
                   "slope_S2": summary(reps["pv_slope_S2"][j]), "zeta_V": summary(reps["zeta_V"][:, j])}
    tests["P-V"] = _assemble("P-V", by_r, r_used, by, top,
                             ("i_M2_beats_M1", "i_M2_beats_M3", "ii_S1_negative", "ii_S2_negative"))
    # P-M (no T range: lags 1..l_max, kernels to L)
    if len(pm_stocks) >= PO.MIN_REACH:
        d = decide_pm(pm_reps["slope"], pm_reps["top"])
        tests["P-M"] = TestResult("P-M", d["p"], d["alternative"], d["conditions"],
                                  {**d, "beta_by_bin": summary(pm_reps["beta"]),
                                   "xlog_by_bin": summary(pm_reps["xlog"]), "slope": summary(pm_reps["slope"]),
                                   "top": summary(pm_reps["top"])},
                                  None, f"lags [{PM_L_MIN},{pm.l_max}], L = {pm.L}", by,
                                  {"stocks": pm_stocks.tolist(), "stock_days": int(mom["used"].sum())})
    else:
        tests["P-M"] = _empty("P-M", "less (both)", ("slope_declines_with_log_q", "top_bin_upper_below_0"), None, by)
        tests["P-M"].range_rule = f"not evaluable: {len(pm_stocks)} stocks reach L = {pm.L} under their ceilings"
        tests["P-M"].stocks_used = {"stocks": pm_stocks.tolist(), "stock_days": int(mom["used"].sum())}
    h = apply_holm(tests)
    descriptive = {"P-K": {RANGE_NAMES[j]: summary(reps["vk_diff"][j]) for j in range(len(C.RANGES))}}
    return {"tests": tests, "holm": h, "descriptive": descriptive,
            "unit": {"name": u.name, "cells": sorted(set(u.cells)), "excluded_cells": list(u.excluded_cells),
                     "n_stocks": S_, "n_days": D_, "pm_params": pm}}


def _assemble(name, by_r: dict, r_used, by, top, conds) -> TestResult:
    if r_used is None:
        t = _empty(name, by_r[C.RANGES[0]]["alternative"], conds, None, by)
        t.by_range = {RANGE_NAMES[i]: v for i, v in enumerate(by_r.values())}
        return t
    d = by_r[r_used]
    return TestResult(name, d["p"], d["alternative"], d["conditions"], d, r_used, range_label(r_used), by,
                      {r_used: top[r_used]}, {RANGE_NAMES[i]: v for i, v in enumerate(by_r.values())})


def run_tick_bin(flows: dict, tb: str, pm: PMParams, rng: np.random.Generator, **kw) -> dict:
    """Primary: the tick bin's confirmatory cells pooled; Holm over its four p-values."""
    return run_unit(tick_bin_unit(flows, tb), pm, rng, **kw)


def run_cell(flows: dict, cell: str, pm: PMParams, rng: np.random.Generator, **kw) -> dict:
    """Secondary: one confirmatory cell alone (its Holm is reported, not a primary decision)."""
    return run_unit(cell_unit(flows, cell), pm, rng, **kw)
