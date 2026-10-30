#!/usr/bin/env python3
"""Descriptive universe table for the eight study days: data/universe.tsv.

DESCRIPTIVE ONLY. Nothing here combines trade signs with price changes:
no sign autocorrelation, propagator, impact decay, imbalance moments, R_a
or collapse statistic touches real data in this step (those are
pre-registered estimands). The tick-size eta below uses mid changes alone.

Inputs per study day d (in ITCH_DATA, default ~/CODE/itch-data):
runs/{d}.runs.npz (runs.build_runs output, trimmed to 10:00-15:30 ET),
tape/{d}.symbols.tsv, tape/{d}.tape.gz (streamed in chunks, never written
out). Locates differ by day, so everything is joined by symbol. Prices are
ITCH units (1e-4 $).

Columns of the table (one row per symbol in any study-day symbol table):

  listing          market_category on the last day it appears
  eligible, reason 1 if every filter passes; else the failed filters:
                   listing (NASDAQ Q/G/S or NYSE N on every day), issue (C on
                   every day), etp (N on every day), presence (runs on all 8
                   days), price (VWAP of displayed executions, sum first_px q /
                   sum q over runs, >= $5 on every day)
  n_min, n_median  min and median over days of runs per day
  vwap_min         lowest daily VWAP, $
  sigma_l          sd of log q over all runs of all days
  q_d1..q_d9       deciles of q (pooled runs; observed values, inverted CDF)
  top10_min/share  top size decile = runs with q >= q_d9: the fewest such runs
                   on any day, and their pooled share (above 0.1 when q ties
                   at q_d9); top20 likewise with q >= q_d8 (top quintile)
  frac_roundlot    share of runs whose q is a multiple of 100 shares (the lattice
                   that a discrete power law on every integer cannot follow)
  csn_*            Clauset-Shalizi-Newman (arXiv:0706.1062) discrete power law
                   p(x) = x^-alpha / zeta(alpha, x_min) on q >= x_min: alpha
                   (density exponent) by exact MLE, x_min by smallest KS
                   distance, n_tail, ks, the KS p-value from 200
                   semiparametric bootstrap refits (blank below 2000 runs; every csn column blank below 50),
                   and the Vuong statistic against a discrete log-normal on
                   the same tail (positive favours the power law) with its
                   two-sided p-value
  s1               PRIMARY TICK-SIZE VARIABLE: share of all runs (8 days) whose
                   pre-run spread ask - bid is exactly one tick (valid quote;
                   tick = 100 units for prices >= $1, 1 unit below)
  eta, eta_nc/na   descriptor only. Dayri-Rosenbaum (arXiv:1207.6325) eta = N_c /
                   (2 N_a), with the mid m = mid2 / 2 read at consecutive runs of
                   a day (valid mids only): among the mid moves of exactly one
                   tick, a move continues (N_c) or reverses (N_a) the previous
                   one-tick move of that stock-day; moves of other sizes are
                   skipped. Pooled over days; undefined without an alternation.
                   It no longer defines the strata: its low bin (eta < 0.35)
                   filled with illiquid wide-spread stocks whose few one-tick
                   moves alternate, while large-tick stocks (s1 >= 0.9) have
                   eta near 0.4 (between consecutive runs their mid moves by a
                   full tick about 90% of the time, the spread having closed
                   again). s1 measures the tick constraint directly.
  eta_half(_nc/na) the same on half-tick mid moves (one tick of mid2)
  frac_onetick     share of nonzero mid moves that are one tick
  frac_halftick    share of nonzero mid moves that are half a tick
  frac_spread1     share of runs WITH A VALID QUOTE with ask - bid = one tick
  spread_ticks     mean of (ask - bid) / tick over runs with a valid quote
  vol_{day}        NASDAQ displayed volume (sum of q over runs) in the TRIMMED
                   10:00-15:30 session, shares
  hid_n            regular-session P executions in 10:00-15:30 with trading state
                   T (from the tape), all days
  hid_mid          share of them printing exactly at their own pre-trade mid
                   (2 price = bid + ask, 0 < bid < ask)
  hid_unsigned     share left unsigned by runs.sign_hidden
  tick_bin,        strata v2, for NASDAQ-listed and (separately) NYSE-listed
  sigma_t          eligible stocks: tick bin 1 = small tick (s1 < 0.4), 2 =
                   intermediate (0.4 <= s1 < 0.8), 3 = large tick (s1 >= 0.8); sigma_l
                   terciles computed within the listing group
  cell, split_v2   the v2 stratum after merges ('Q-T3S2' NASDAQ, 'N-T3S2' NYSE) and
                   'explore' / 'confirm': NASDAQ with default_rng(SPLIT_SEED_V2),
                   NYSE with its own default_rng(SPLIT_SEED_NYSE)
  eta_bin_v1,      the superseded v1 strata (NASDAQ only, eta < 0.35, 0.35 to 0.5,
  cell_v1,         >= 0.5 by the same sigma_l terciles, default_rng(SPLIT_SEED_V1);
  split_v1         undefined eta left out), kept for the record

Split: cells are taken in sorted label order, symbols sorted inside a cell,
one generator per listing group draws a permutation per cell, and the first
round(0.3 n) are 'explore'. A cell whose confirmation count n - round(0.3 n)
is below 10 is merged first into its sigma_l neighbour in the same tick bin
(tercile 1 or 3 into 2; tercile 2 into the smaller neighbour, ties to 1),
repeatedly.

CSN details: x_min is scanned over observed values of q with at least 50
runs at or above them, at most 40 of them, log-spaced (CSN scan every
value; this caps the cost). The discrete power law is sampled exactly
(continuous-approximation start, then corrected with the exact survival
function); the bootstrap below x_min resamples the empirical body. The
log-normal is truncated at x_min with mass Q(x - 1/2) - Q(x + 1/2) on integer
x; its parameters maximise the midpoint-density likelihood, and the Vuong
ratio uses the exact interval masses.

    python3 tools/est/describe.py [--data DIR] [--out data/universe.tsv] [--workers 4]
"""
from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tape"))
import runs as R  # noqa: E402
import tape  # noqa: E402

STUDY_DAYS = ("01302019", "03272019", "07302019", "08302019", "10182019", "10302019", "12302019", "01302020")
NASDAQ = ("Q", "G", "S")
NYSE = ("N",)
PRICE_MIN = 5 * 10_000  # $5 in ITCH units
SPLIT_SEED_V1 = 20261029  # superseded strata (eta bins)
SPLIT_SEED_V2 = 20261030  # NASDAQ-listed, s1 tick bins
SPLIT_SEED_NYSE = 20261031  # NYSE-listed, s1 tick bins
CSN_SEED = 20261029
EXPLORE_FRAC = 0.3
MIN_CONFIRM = 10
ETA_EDGES = (0.35, 0.5)
TICK_EDGES = (0.4, 0.8)
N_BOOT = 200
MIN_P_RUNS = 2000


# ---- Hurwitz zeta and the discrete power law -------------------------------------

_BERNOULLI = (1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66)  # B2 .. B10


def hurwitz(s, a) -> np.ndarray:
    """zeta(s, a) = sum_{k>=0} (a + k)^-s for s > 1, a > 0: direct terms until a + m >= 10,
    then Euler-Maclaurin with Bernoulli numbers to B10 (error below 1e-13 relative)."""
    s, a = np.broadcast_arrays(np.asarray(s, dtype=float), np.asarray(a, dtype=float))
    out = np.zeros(s.shape)
    m = np.clip(np.ceil(10.0 - a), 0, None)
    for k in range(int(m.max()) if m.size else 0):
        sel = k < m
        out[sel] += (a[sel] + k) ** -s[sel]
    b = a + m
    out += b ** (1 - s) / (s - 1) + 0.5 * b**-s
    rising, fact, pw = s.copy(), 2.0, b ** (-s - 1)
    for j, B in enumerate(_BERNOULLI, start=1):
        out += B / fact * rising * pw
        rising = rising * (s + 2 * j - 1) * (s + 2 * j)
        fact *= (2 * j + 1) * (2 * j + 2)
        pw = pw / (b * b)
    return out


def alpha_mle(xmin, n_tail, sum_log, lo: float = 1.0001, hi: float = 10.0, iters: int = 45) -> np.ndarray:
    """Exact discrete MLE: argmax over alpha of -n log zeta(alpha, xmin) - alpha sum log x
    (concave), by golden section on [lo, hi] (width below 3e-9 after 45 steps), vectorized over
    candidate x_min."""
    xmin, n_tail, sum_log = (np.asarray(v, dtype=float) for v in (xmin, n_tail, sum_log))
    g = (np.sqrt(5) - 1) / 2

    def f(al):  # negative log-likelihood
        return n_tail * np.log(hurwitz(al, xmin)) + al * sum_log

    a, b = np.full(xmin.shape, lo), np.full(xmin.shape, hi)
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        left = fc < fd  # keep [a, d]: old c becomes the new d; else keep [c, b]: old d becomes the new c
        a, b = np.where(left, a, c), np.where(left, d, b)
        x = np.where(left, b - g * (b - a), a + g * (b - a))
        fx = f(x)
        c, d, fc, fd = np.where(left, x, d), np.where(left, c, x), np.where(left, fx, fd), np.where(left, fc, fx)
    return (a + b) / 2


def ks_distance(u: np.ndarray, c: np.ndarray, i0: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    """KS distance between the tail u[i0:] (counts c) and the discrete power law with x_min = u[i0],
    for each candidate. The model CDF also rises at integers between observed values, so it is
    compared with the empirical CDF both at each u_j and just below it, at u_j - 1."""
    cum = np.cumsum(c)
    before = np.where(i0 > 0, cum[np.maximum(i0 - 1, 0)], 0)[:, None]
    n_t = cum[-1] - before
    j = np.arange(len(u))[None, :]
    inside = j >= i0[:, None]
    xmin = u[i0][:, None]
    al = alpha[:, None]
    zmin = hurwitz(al, xmin)
    F_emp = (cum[None, :] - before) / n_t
    F_at = 1 - hurwitz(al, u[None, :] + 1.0) / zmin  # P(X <= u_j)
    F_below = 1 - hurwitz(al, u[None, :]) / zmin  # P(X <= u_j - 1)
    F_emp_prev = np.concatenate([np.zeros((len(i0), 1)), F_emp[:, :-1]], axis=1)
    F_emp_prev = np.where(j == i0[:, None], 0.0, F_emp_prev)
    d = np.maximum(np.abs(F_emp - F_at), np.abs(F_emp_prev - F_below))
    return np.max(np.where(inside, d, 0.0), axis=1)


def xmin_candidates(u: np.ndarray, c: np.ndarray, min_tail: int = 50, max_cand: int = 40) -> np.ndarray:
    """Indices into u of candidate x_min: values with >= min_tail runs at or above, at most
    max_cand of them, log-spaced in value."""
    tail = np.cumsum(c[::-1])[::-1]
    ok = np.flatnonzero(tail >= min_tail)
    if len(ok) <= max_cand:
        return ok
    grid = np.exp(np.linspace(np.log(u[ok[0]]), np.log(u[ok[-1]]), max_cand))
    return np.unique(ok[np.clip(np.searchsorted(u[ok], grid), 0, len(ok) - 1)])


def fit_powerlaw(u: np.ndarray, c: np.ndarray) -> tuple[float, float, int, float]:
    """(alpha, x_min, n_tail, ks) for a histogram (sorted unique values u, counts c)."""
    i0 = xmin_candidates(u, c)
    if len(i0) == 0:
        raise ValueError("too few observations for a tail fit")
    lc = c * np.log(u)
    n_t = np.cumsum(c[::-1])[::-1][i0]
    s_l = np.cumsum(lc[::-1])[::-1][i0]
    al = alpha_mle(u[i0], n_t, s_l)
    D = ks_distance(u, c, i0, al)
    k = int(np.argmin(D))
    return float(al[k]), float(u[i0[k]]), int(n_t[k]), float(D[k])


def sample_powerlaw(n: int, alpha: float, xmin: float, rng: np.random.Generator) -> np.ndarray:
    """Exact draws from p(x) = x^-alpha / zeta(alpha, xmin), x = xmin, xmin + 1, ...:
    X = max{x : S(x) >= v}, S(x) = zeta(alpha, x) / zeta(alpha, xmin), v ~ U(0, 1]."""
    v = 1.0 - rng.random(n)
    x = np.maximum(xmin, np.floor((xmin - 0.5) * v ** (-1 / (alpha - 1)) + 0.5))
    zmin = float(hurwitz(alpha, xmin))
    active = np.arange(n)
    for _ in range(200):
        if len(active) == 0:
            break
        xa, va = x[active], v[active]
        up = hurwitz(alpha, xa + 1) / zmin >= va
        down = (hurwitz(alpha, xa) / zmin < va) & (xa > xmin)
        x[active[up]] += 1
        x[active[down]] -= 1
        active = active[up | down]
    return x


def _merge_hist(u1, c1, u2, c2):
    u = np.concatenate([u1, u2])
    c = np.concatenate([c1, c2])
    uu, inv = np.unique(u, return_inverse=True)
    return uu, np.bincount(inv, weights=c).astype(np.int64)


def ks_pvalue(u, c, alpha, xmin, ks, n_boot: int, rng: np.random.Generator) -> float:
    """CSN semiparametric bootstrap: body below x_min resampled from the data, tail drawn from
    the fitted law with probability n_tail / n, every synthetic set refitted (x_min scan too)."""
    n = int(c.sum())
    body = u < xmin
    ub, cb = u[body], c[body]
    p_tail = c[~body].sum() / n
    worse = 0
    for _ in range(n_boot):
        nt = rng.binomial(n, p_tail)
        cb_s = rng.multinomial(n - nt, cb / cb.sum()) if len(ub) else np.zeros(0, dtype=np.int64)
        ut, ct = np.unique(sample_powerlaw(nt, alpha, xmin, rng), return_counts=True)
        us, cs = _merge_hist(ub[cb_s > 0], cb_s[cb_s > 0], ut, ct)
        worse += fit_powerlaw(us, cs)[3] >= ks
    return worse / n_boot


def _log_sf(z: np.ndarray) -> np.ndarray:
    """log of the standard normal upper tail Q(z), accurate in the far tail."""
    z = np.asarray(z, dtype=float)
    out = np.empty(z.shape)
    small = z < 30
    erfc = np.frompyfunc(math.erfc, 1, 1)
    out[small] = np.log(0.5 * erfc(z[small] / math.sqrt(2)).astype(float))
    zb = z[~small]
    out[~small] = -0.5 * zb**2 - np.log(zb * math.sqrt(2 * math.pi)) + np.log1p(-1 / zb**2 + 3 / zb**4)
    return out


def _nelder_mead(f, x0, step, iters: int = 400, tol: float = 1e-10):
    pts = [np.asarray(x0, dtype=float)] + [np.asarray(x0, dtype=float) + np.eye(len(x0))[i] * step[i]
                                           for i in range(len(x0))]
    vals = [f(p) for p in pts]
    for _ in range(iters):
        order = np.argsort(vals)
        pts, vals = [pts[i] for i in order], [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol * (1 + abs(vals[0])):
            break
        cen = np.mean(pts[:-1], axis=0)
        xr = cen + (cen - pts[-1])
        fr = f(xr)
        if fr < vals[0]:
            xe = cen + 2 * (cen - pts[-1])
            fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = cen + 0.5 * (pts[-1] - cen)
            fc = f(xc)
            if fc < vals[-1]:
                pts[-1], vals[-1] = xc, fc
            else:
                pts = [pts[0]] + [pts[0] + 0.5 * (p - pts[0]) for p in pts[1:]]
                vals = [vals[0]] + [f(p) for p in pts[1:]]
    return pts[int(np.argmin(vals))]


def lognormal_tail_fit(u: np.ndarray, c: np.ndarray, xmin: float) -> tuple[float, float]:
    """(mu, sigma) of a log-normal truncated at x_min - 1/2, by the midpoint-density likelihood."""
    tail = u >= xmin
    lg, w = np.log(u[tail]), c[tail].astype(float)
    n, s1, s2 = w.sum(), w @ lg, w @ (lg * lg)
    lx = math.log(xmin - 0.5)
    m0 = s1 / n
    sd0 = math.sqrt(max(s2 / n - m0 * m0, 1e-6))

    def nll(th):
        mu, ls = th
        sg = math.exp(ls)
        ll = -s1 - n * ls - (s2 - 2 * mu * s1 + n * mu * mu) / (2 * sg * sg)
        return -(ll - n * float(_log_sf(np.array([(lx - mu) / sg]))[0]))

    mu, ls = _nelder_mead(nll, [m0, math.log(sd0)], [0.5, 0.3])
    return float(mu), float(math.exp(ls))


def vuong(u: np.ndarray, c: np.ndarray, alpha: float, xmin: float) -> tuple[float, float]:
    """Normalized log-likelihood ratio (power law minus log-normal) on q >= x_min and its two-sided p."""
    mu, sg = lognormal_tail_fit(u, c, xmin)
    tail = u >= xmin
    x, w = u[tail], c[tail].astype(float)
    l_pl = -alpha * np.log(x) - np.log(hurwitz(alpha, xmin))
    za, zb = (np.log(x - 0.5) - mu) / sg, (np.log(x + 0.5) - mu) / sg
    la, lb = _log_sf(za), _log_sf(zb)
    l_ln = la + np.log(-np.expm1(lb - la)) - _log_sf(np.array([(math.log(xmin - 0.5) - mu) / sg]))[0]
    d = l_pl - l_ln
    n = w.sum()
    mean = (w @ d) / n
    sd = math.sqrt(max((w @ (d - mean) ** 2) / n, 1e-300))
    stat = float(w @ d / (sd * math.sqrt(n)))
    return stat, float(math.erfc(abs(stat) / math.sqrt(2)))


CSN_FIELDS = ("csn_alpha", "csn_xmin", "csn_ntail", "csn_ks", "csn_p", "vuong", "vuong_p")


def csn(q: np.ndarray, rng: np.random.Generator, n_boot: int = N_BOOT, min_p_runs: int = MIN_P_RUNS) -> dict:
    u, c = np.unique(np.asarray(q, dtype=float), return_counts=True)
    if len(xmin_candidates(u, c)) == 0:  # fewer than 50 runs: no tail to fit
        return dict.fromkeys(CSN_FIELDS, float("nan"))
    alpha, xmin, n_tail, ks = fit_powerlaw(u, c)
    p = ks_pvalue(u, c, alpha, xmin, ks, n_boot, rng) if c.sum() >= min_p_runs and n_boot > 0 else float("nan")
    v, vp = vuong(u, c, alpha, xmin)
    return dict(csn_alpha=alpha, csn_xmin=xmin, csn_ntail=n_tail, csn_ks=ks, csn_p=p, vuong=v, vuong_p=vp)


# ---- per-day statistics ------------------------------------------------------------


def read_symbols(path) -> dict[int, tuple[str, str, str, str]]:
    """locate -> (symbol, market_category, issue_classification, etp_flag)."""
    out = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head != ["locate", "symbol", "market_category", "issue_classification", "etp_flag"]:
            raise ValueError(f"unexpected symbols header {head}")
        for line in f:
            loc, sym, mc, ic, etp = line.rstrip("\n").split("\t")
            out[int(loc)] = (sym, mc, ic, etp)
    return out


def tick_of(mid2: np.ndarray) -> np.ndarray:
    """Tick in ITCH units: 100 ($0.01) at prices >= $1 (mid2 >= 20000), else 1."""
    return np.where(np.asarray(mid2) >= 20_000, 100, 1).astype(np.int64)


def eta_counts(locate: np.ndarray, mid2: np.ndarray, n_loc: int, size: int = 2):
    """Per locate (index): continuations, alternations and number of the mid2 moves of exactly
    `size` ticks (2: one-tick mid moves, 1: half-tick mid moves), and all nonzero moves, from
    consecutive runs (all with valid mids, in time order within each locate)."""
    loc = np.asarray(locate, dtype=np.int64)
    m2 = np.asarray(mid2, dtype=np.int64)
    dm = m2[1:] - m2[:-1]
    same = loc[1:] == loc[:-1]
    one = same & (np.abs(dm) == size * tick_of(m2[:-1]))  # mid2 = 2 mid: a one-tick mid move is 2 ticks of mid2
    moves = same & (dm != 0)
    p = np.flatnonzero(one)
    s, lp = np.sign(dm[p]), loc[1:][p]
    pair = lp[1:] == lp[:-1]
    cont = pair & (s[1:] == s[:-1])
    alt = pair & (s[1:] != s[:-1])
    nc = np.bincount(lp[1:][cont], minlength=n_loc)
    na = np.bincount(lp[1:][alt], minlength=n_loc)
    n1 = np.bincount(lp, minlength=n_loc)
    nm = np.bincount(loc[1:][moves], minlength=n_loc)
    return nc, na, n1, nm


def eta(nc, na):
    """N_c / (2 N_a); undefined (NaN) without any alternation."""
    nc, na = np.asarray(nc, dtype=float), np.asarray(na, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(na > 0, nc / (2 * na), np.nan)


def day_stats(runs: np.ndarray) -> dict[str, np.ndarray]:
    """Per-locate statistics of one day's runs (indexed by locate)."""
    n_loc = int(runs["locate"].max()) + 1 if len(runs) else 1
    loc = runs["locate"].astype(np.int64)
    q = runs["q"].astype(float)
    out = {
        "n": np.bincount(loc, minlength=n_loc),
        "sum_q": np.bincount(loc, q, minlength=n_loc),
        "sum_pq": np.bincount(loc, runs["first_px"].astype(float) * q, minlength=n_loc),
    }
    v = runs[runs["mid_valid"] == 1]
    lv = v["locate"].astype(np.int64)
    out["nc"], out["na"], out["n1"], out["nm"] = eta_counts(lv, v["mid2"], n_loc, size=2)
    out["nch"], out["nah"], out["nh"], _ = eta_counts(lv, v["mid2"], n_loc, size=1)
    spread = v["ask_px"].astype(np.int64) - v["bid_px"].astype(np.int64)
    tk = tick_of(v["mid2"])
    out["n_quote"] = np.bincount(lv, minlength=n_loc)
    out["n_spread1"] = np.bincount(lv, (spread == tk).astype(float), minlength=n_loc)
    out["sum_spread_ticks"] = np.bincount(lv, spread / tk, minlength=n_loc)
    return out


def tick_stats(a: dict) -> dict:
    """s1 over all runs, frac_spread1 and spread_ticks over runs with a valid quote, from summed day_stats."""
    return {"s1": a["n_spread1"] / a["n_runs"], "frac_spread1": a["n_spread1"] / a["n_quote"],
            "spread_ticks": a["sum_spread_ticks"] / a["n_quote"]}


def spearman(x, y) -> float:
    """Spearman rank correlation with average ranks for ties, over pairs where both are finite."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)

    def rank(v):
        order = np.argsort(v, kind="stable")
        r = np.empty(len(v))
        r[order] = np.arange(1, len(v) + 1)
        u, inv = np.unique(v, return_inverse=True)
        return (np.bincount(inv, r) / np.bincount(inv))[inv]  # average rank of each tie group

    rx, ry = rank(x[ok]), rank(y[ok])
    return float(np.corrcoef(rx, ry)[0, 1])


def iter_tape_chunks(path, records_per_chunk: int = 1 << 20):
    """Stream a tape in numpy chunks (no full decompression in memory or on disk)."""
    dt = tape.numpy_dtype()
    with gzip.open(path, "rb") as f:
        tape._check_header(f.read(tape.HEADER.size))
        while True:
            b = f.read(tape.RECORD.size * records_per_chunk)
            if not b:
                return
            if len(b) % tape.RECORD.size:
                raise ValueError("truncated tape")
            yield np.frombuffer(b, dtype=dt)


def hidden_prints(chunks, session=(R.SESSION_LO, R.SESSION_HI)) -> np.ndarray:
    """Regular-session P executions in the session window with trading state T."""
    keep = []
    for t in chunks:
        sel = (t["type"] == b"P") & (t["phase"] == tape.REGULAR) & (t["tstate"] == b"T")
        sel &= (t["ts"] >= session[0]) & (t["ts"] < session[1])
        keep.append(t[sel].copy())
    return np.concatenate(keep) if keep else np.zeros(0, dtype=tape.numpy_dtype())


def hidden_stats(p: np.ndarray, runs: np.ndarray, n_loc: int) -> dict[str, np.ndarray]:
    loc = p["locate"].astype(np.int64)
    bid, ask = p["bid_px"].astype(np.int64), p["ask_px"].astype(np.int64)
    at_mid = (bid > 0) & (bid < ask) & (2 * p["price"].astype(np.int64) == bid + ask)
    unsigned = R.sign_hidden(p, runs) == 0
    return {"hid_n": np.bincount(loc, minlength=n_loc),
            "hid_mid": np.bincount(loc, at_mid.astype(float), minlength=n_loc),
            "hid_uns": np.bincount(loc, unsigned.astype(float), minlength=n_loc)}


# ---- eligibility, strata, split ----------------------------------------------------

FILTERS = ("listing", "issue", "etp", "presence", "price")


def eligibility(days: list[dict]) -> dict[str, list[str]]:
    """Failed filters per symbol. days: one dict per study day, symbol -> dict(mc, ic, etp, n, vwap)
    (n and vwap absent or 0 when the symbol has no runs that day)."""
    syms = sorted(set().union(*[d.keys() for d in days]))
    out = {}
    for s in syms:
        rows = [d.get(s) for d in days]
        fail = []
        if not all(r is not None and r["mc"] in NASDAQ + NYSE for r in rows):
            fail.append("listing")
        if not all(r is not None and r["ic"] == "C" for r in rows):
            fail.append("issue")
        if not all(r is not None and r["etp"] == "N" for r in rows):
            fail.append("etp")
        if not all(r is not None and r.get("n", 0) > 0 for r in rows):
            fail.append("presence")
        if not all(r is not None and r.get("n", 0) > 0 and r["vwap"] >= PRICE_MIN for r in rows):
            fail.append("price")
        out[s] = fail
    return out


def eta_bin(e: np.ndarray) -> np.ndarray:
    return 1 + (np.asarray(e) >= ETA_EDGES[0]).astype(int) + (np.asarray(e) >= ETA_EDGES[1]).astype(int)


def tick_bin(s1: np.ndarray) -> np.ndarray:
    """1 = small tick (s1 < 0.4), 2 = intermediate, 3 = large tick (s1 >= 0.8)."""
    s1 = np.asarray(s1)
    return 1 + (s1 >= TICK_EDGES[0]).astype(int) + (s1 >= TICK_EDGES[1]).astype(int)


def terciles(x: np.ndarray, edges=None) -> np.ndarray:
    """1, 2, 3 by the tercile edges of x (or the given edges)."""
    x = np.asarray(x, dtype=float)
    lo, hi = np.quantile(x, [1 / 3, 2 / 3]) if edges is None else edges
    return 1 + (x >= lo).astype(int) + (x >= hi).astype(int)


def n_explore(n: int) -> int:
    return int(math.floor(EXPLORE_FRAC * n + 0.5))


def merge_cells(eb: np.ndarray, st: np.ndarray, prefix: str = "E") -> tuple[np.ndarray, list[str]]:
    """Cell labels '{prefix}{bin}S{terciles}' after merges. Within each bin, while some cell has
    fewer than MIN_CONFIRM confirmation stocks and more than one cell is left, the smallest such cell
    (ties: lower tercile) joins its smaller adjacent cell (ties: lower)."""
    eb, st = np.asarray(eb, dtype=int), np.asarray(st, dtype=int)
    label, merges = {}, []
    for e in sorted(set(eb.tolist())):
        groups = [(t,) for t in sorted(set(st[eb == e].tolist()))]
        size = {g: int(np.sum((eb == e) & (st == g[0]))) for g in groups}
        while len(groups) > 1:
            small = [g for g in groups if size[g] - n_explore(size[g]) < MIN_CONFIRM]
            if not small:
                break
            g = min(small, key=lambda h: (size[h], h))
            i = groups.index(g)
            target = min((groups[j] for j in (i - 1, i + 1) if 0 <= j < len(groups)), key=lambda h: (size[h], h))
            new = tuple(sorted(g + target))
            gl, tl = "+".join(map(str, g)), "+".join(map(str, target))
            merges.append(f"{prefix}{e}: S{gl} ({size[g]} stocks) into S{tl}")
            size[new] = size.pop(g) + size.pop(target)
            groups = sorted([h for h in groups if h not in (g, target)] + [new])
        for g in groups:
            for t in g:
                label[(e, t)] = f"{prefix}{e}S" + "+".join(map(str, g))
    return np.array([label[(e, t)] for e, t in zip(eb.tolist(), st.tolist(), strict=True)]), merges


def split(symbols, cells, seed: int) -> dict[str, str]:
    """30/70 explore/confirm within each cell: cells in sorted order, symbols sorted, one generator."""
    rng = np.random.default_rng(seed)
    out = {}
    symbols, cells = np.asarray(symbols), np.asarray(cells)
    for cell in sorted(set(cells.tolist())):
        members = np.sort(symbols[cells == cell])
        perm = rng.permutation(len(members))
        k = n_explore(len(members))
        for rank, i in enumerate(perm):
            out[str(members[i])] = "explore" if rank < k else "confirm"
    return out


# ---- driver -------------------------------------------------------------------------


def _csn_task(args):
    sym, q, n_boot = args
    rng = np.random.default_rng([CSN_SEED, zlib.crc32(sym.encode())])
    return sym, csn(q, rng, n_boot)


ACC_DAY = ("nc", "na", "n1", "nm", "nch", "nah", "nh", "n_quote", "n_spread1", "sum_spread_ticks")


def collect(data: Path, log=print):
    """One pass over the study days, one day in memory at a time: per-day symbol rows, pooled
    run sizes and summed quote, eta and hidden-print counts of every candidate symbol."""
    per_day, qs, acc = [], {}, {}
    for d in STUDY_DAYS:
        syms = read_symbols(data / "tape" / f"{d}.symbols.tsv")
        runs = np.load(data / "runs" / f"{d}.runs.npz")["runs"]
        st = day_stats(runs)
        n_loc = len(st["n"])
        p = hidden_prints(iter_tape_chunks(data / "tape" / f"{d}.tape.gz"))
        hs = hidden_stats(p, runs, max(n_loc, int(p["locate"].max()) + 1 if len(p) else 1))
        day = {}
        loc_sorted = runs["locate"]  # runs are sorted by (locate, ts)
        for loc, (sym, mc, ic, etp) in syms.items():
            row = {"mc": mc, "ic": ic, "etp": etp, "n": 0, "vwap": 0.0}
            if loc < n_loc and st["n"][loc] > 0:
                row.update(n=int(st["n"][loc]), vwap=st["sum_pq"][loc] / st["sum_q"][loc], vol=st["sum_q"][loc])
                if mc in NASDAQ + NYSE and ic == "C" and etp == "N":
                    lo, hi = np.searchsorted(loc_sorted, [loc, loc + 1])
                    qs.setdefault(sym, []).append(runs["q"][lo:hi].astype(np.uint32))
                    a = acc.setdefault(sym, dict.fromkeys(ACC_DAY + ("n_runs", "hid_n", "hid_mid", "hid_uns"), 0.0))
                    for k in ACC_DAY:
                        a[k] += float(st[k][loc])
                    a["n_runs"] += float(st["n"][loc])
                    for k in ("hid_n", "hid_mid", "hid_uns"):
                        a[k] += float(hs[k][loc]) if loc < len(hs[k]) else 0.0
            day[sym] = row
        per_day.append(day)
        log(f"{d}: {len(runs):,} runs, {len(syms):,} symbols, {len(p):,} hidden prints")
        del runs, p, st, hs
    return per_day, qs, acc


def analyse(per_day, qs, acc, workers: int, n_boot: int = N_BOOT, log=print) -> tuple[list[dict], dict]:
    fails = eligibility(per_day)
    elig = sorted(s for s, f in fails.items() if not f)
    log(f"eligible: {len(elig)}; running CSN with {n_boot} bootstraps on {workers} workers")
    fits = {}
    tasks = [(s, np.concatenate(qs[s]), n_boot) for s in elig]
    if workers > 1:
        with ProcessPoolExecutor(workers) as ex:
            for i, (s, r) in enumerate(ex.map(_csn_task, tasks, chunksize=8)):
                fits[s] = r
                if (i + 1) % 250 == 0:
                    log(f"  CSN {i + 1}/{len(tasks)}")
    else:
        fits = dict(map(_csn_task, tasks))
    rows = []
    for s in sorted(fails):
        days = [d.get(s) for d in per_day]
        last = [d for d in days if d is not None][-1]
        row = {"symbol": s, "listing": last["mc"], "eligible": int(not fails[s]), "reason": ",".join(fails[s])}
        ns = [d["n"] if d else 0 for d in days]
        row["n_min"], row["n_median"] = int(min(ns)), float(np.median(ns))
        vw = [d["vwap"] for d in days if d and d["n"] > 0]
        row["vwap_min"] = min(vw) / 10_000 if vw else float("nan")
        for dname, dd in zip(STUDY_DAYS, days, strict=True):
            row[f"vol_{dname}"] = int(dd.get("vol", 0)) if dd else 0
        if not fails[s]:
            q = np.concatenate(qs[s]).astype(float)
            row["sigma_l"] = float(np.std(np.log(q)))
            dec = np.quantile(q, np.arange(1, 10) / 10, method="inverted_cdf")
            for i, v in enumerate(dec, start=1):
                row[f"q_d{i}"] = int(v)
            row["top10_min"] = int(min(int(np.sum(x >= dec[8])) for x in qs[s]))
            row["top20_min"] = int(min(int(np.sum(x >= dec[7])) for x in qs[s]))
            row["top10_share"] = float(np.mean(q >= dec[8]))
            row["top20_share"] = float(np.mean(q >= dec[7]))
            row["frac_roundlot"] = float(np.mean(q % 100 == 0))
            row.update(fits[s])
            a = acc[s]
            row["eta_nc"], row["eta_na"] = int(a["nc"]), int(a["na"])
            row["eta"] = float(eta(a["nc"], a["na"]))
            row["eta_half_nc"], row["eta_half_na"] = int(a["nch"]), int(a["nah"])
            row["eta_half"] = float(eta(a["nch"], a["nah"]))
            row["frac_onetick"] = a["n1"] / a["nm"] if a["nm"] else float("nan")
            row["frac_halftick"] = a["nh"] / a["nm"] if a["nm"] else float("nan")
            row.update(tick_stats(a))
            row["hid_n"] = int(a["hid_n"])
            row["hid_mid"] = a["hid_mid"] / a["hid_n"] if a["hid_n"] else float("nan")
            row["hid_unsigned"] = a["hid_uns"] / a["hid_n"] if a["hid_n"] else float("nan")
        rows.append(row)
    meta = {"csn_seed": CSN_SEED, "n_boot": n_boot, "numpy": np.__version__}
    meta["v1"] = stratify_v1(rows)
    meta["v2"] = stratify_v2(rows)
    el = [r for r in rows if r["eligible"]]
    meta["spearman"] = {g: {k: spearman([r["s1"] for r in m], [r[k] for r in m]) for k in ("eta", "eta_half")}
                        for g, m in (("NASDAQ", [r for r in el if r["listing"] in NASDAQ]),
                                     ("NYSE", [r for r in el if r["listing"] in NYSE]), ("all", el))}
    return rows, meta


def stratify_group(rows: list[dict], listings, value: str, bin_fn, prefix: str, seed: int, cols) -> dict:
    """Strata and split of the eligible rows of one listing group, in place: bins of `value`
    (bin_fn) by sigma_l terciles computed on the whole group, merges, then the split with its own
    generator. A stock whose `value` is undefined gets no stratum. cols names the bin, tercile
    (or None), cell and split fields."""
    group = [r for r in rows if r["eligible"] and r["listing"] in listings]
    edges = np.quantile([r["sigma_l"] for r in group], [1 / 3, 2 / 3])
    strat = [r for r in group if np.isfinite(r[value])]
    b = bin_fn(np.array([r[value] for r in strat]))
    t = terciles([r["sigma_l"] for r in strat], edges)
    cells, merges = merge_cells(b, t, prefix)
    assign = split([r["symbol"] for r in strat], cells, seed)
    bin_col, terc_col, cell_col, split_col = cols
    for r, bb, tt, c in zip(strat, b, t, cells, strict=True):
        r[bin_col], r[cell_col], r[split_col] = int(bb), str(c), assign[r["symbol"]]
        if terc_col:
            r[terc_col] = int(tt)
    return {"seed": seed, "merges": merges, "sigma_terciles": edges.tolist(), "n": len(strat),
            "unstratified": sorted(r["symbol"] for r in group if not np.isfinite(r[value]))}


def stratify_v1(rows: list[dict]) -> dict:
    """The superseded v1 strata: NASDAQ-listed, eta bins, seed SPLIT_SEED_V1 (reproduces commit 73fd12f)."""
    return stratify_group(rows, NASDAQ, "eta", eta_bin, "E", SPLIT_SEED_V1, ("eta_bin_v1", None, "cell_v1", "split_v1"))


def stratify_v2(rows: list[dict]) -> dict:
    """v2 strata on s1 tick bins: NASDAQ-listed (SPLIT_SEED_V2) and NYSE-listed (SPLIT_SEED_NYSE), each
    with its own sigma_l terciles and its own generator, so neither split depends on the other group."""
    cols = ("tick_bin", "sigma_t", "cell", "split_v2")
    return {"NASDAQ": stratify_group(rows, NASDAQ, "s1", tick_bin, "Q-T", SPLIT_SEED_V2, cols),
            "NYSE": stratify_group(rows, NYSE, "s1", tick_bin, "N-T", SPLIT_SEED_NYSE, cols)}


COLUMNS = (["symbol", "listing", "eligible", "reason", "n_min", "n_median", "vwap_min", "sigma_l"]
           + [f"q_d{i}" for i in range(1, 10)]
           + ["top10_min", "top20_min", "top10_share", "top20_share", "frac_roundlot"] + list(CSN_FIELDS)
           + ["s1", "eta", "eta_nc", "eta_na", "eta_half", "eta_half_nc", "eta_half_na", "frac_onetick",
              "frac_halftick", "frac_spread1", "spread_ticks"]
           + [f"vol_{d}" for d in STUDY_DAYS]
           + ["hid_n", "hid_mid", "hid_unsigned", "tick_bin", "sigma_t", "cell", "split_v2",
              "eta_bin_v1", "cell_v1", "split_v1"])


def _fmt(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def write_tsv(rows, meta, path: Path) -> None:
    with open(path, "w") as f:
        f.write("# itchbook universe: descriptive statistics of the 8 study days (tools/est/describe.py)\n")
        f.write(f"# study_days: {' '.join(STUDY_DAYS)}\n")
        v1, v2 = meta["v1"], meta["v2"]
        f.write(f"# numpy {meta['numpy']}; explore fraction {EXPLORE_FRAC} within cells; merge below {MIN_CONFIRM} "
                "confirmation stocks into the sigma_l neighbour in the same bin\n")
        f.write("# split_v2 (CURRENT): strata changed in step 3b because eta mismeasures tick size (its low bin "
                "held illiquid wide-spread stocks, while large-tick stocks have eta near 0.4); tick bins now "
                f"use s1: < {TICK_EDGES[0]}, {TICK_EDGES[0]} to {TICK_EDGES[1]}, >= {TICK_EDGES[1]}; new seeds\n")
        for g in ("NASDAQ", "NYSE"):
            m = v2[g]
            lo, hi = m["sigma_terciles"]
            f.write(f"#   {g}-listed: default_rng({m['seed']}), {m['n']} stocks, sigma_l terciles {lo:.6g} {hi:.6g}, "
                    f"merges: {'; '.join(m['merges']) or 'none'}\n")
        lo, hi = v1["sigma_terciles"]
        f.write(f"# split_v1 (SUPERSEDED, commit 73fd12f): NASDAQ-listed only, default_rng({v1['seed']}), eta bins < "
                f"{ETA_EDGES[0]}, {ETA_EDGES[0]} to {ETA_EDGES[1]}, >= {ETA_EDGES[1]}, sigma_l terciles {lo:.6g} "
                f"{hi:.6g}, merges: {'; '.join(v1['merges']) or 'none'}; {len(v1['unstratified'])} with undefined "
                f"eta left out: {' '.join(v1['unstratified'])}\n")
        sp = meta["spearman"]
        f.write("# spearman(s1, eta) / (s1, eta_half): " + "; ".join(
            f"{g} {sp[g]['eta']:.3f} / {sp[g]['eta_half']:.3f}" for g in ("NASDAQ", "NYSE", "all")) + "\n")
        f.write(f"# csn: {meta['n_boot']} bootstraps, rng [{meta['csn_seed']}, crc32(symbol)]; "
                f"p-value blank below {MIN_P_RUNS} runs\n")
        f.write("# vol_*: NASDAQ displayed volume in the trimmed 10:00-15:30 session (shares)\n")
        f.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(r.get(c)) for c in COLUMNS) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", type=Path, default=Path(os.environ.get("ITCH_DATA", "~/CODE/itch-data")).expanduser())
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[2] / "data" / "universe.tsv")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=N_BOOT)
    a = ap.parse_args()
    rows, meta = analyse(*collect(a.data), a.workers, a.n_boot)
    write_tsv(rows, meta, a.out)
    print(f"wrote {a.out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
