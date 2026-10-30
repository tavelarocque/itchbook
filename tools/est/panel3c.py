"""Phase 3c panels: batched generators for many short stock-days, and the
semi-synthetic generator.

A panel is n_stocks x n_days stock-days stored stock-major as concatenated runs
(eps, q, q_mom, lm) with per-stock-day (start, length). The generators follow
gen3 (G0, GT, GF, same parameters and conventions) but generate every
stock-day of the panel at once: stock-days are grouped by the circulant size
M = 2^ceil(log2(2 N)) and each group's sign series and kernel convolutions are
one batched FFT. A Gaussian series of length N <= M / 2 from a circulant
embedding of size M of its ACF is exact whenever the embedding is
non-negative definite (checked), as in gen3.gaussian_acf with M = 2N.

Signs: "fgn_ar1" is the 3a model gen3.sign_model_r (w1, H, w2, phi);
"two_component" is the 3b model sign(fGn(H) + m_d), m_d ~ N(0, tau^2) per
stock-day. Sizes: log-normal with the cell's (mu, sigma_l).

Clip (the moments' sizes only, never the impact): q_mom = min(q, cap_s) with
cap_s the stock's own 99.9th percentile of run size over its days (R2, as on
real data, clip3b.stock_caps), or a fixed cap.

Day lengths: each synthetic stock copies the n_days day lengths of an
exploration stock of the cell drawn at random (with replacement), so N_s
follows the real distribution.

Semi-synthetic (semi_synthetic): the real exploration signs and sizes of each
stock-day pass through unchanged; only the mid is synthetic, from the G0, GT
or GF response to the RAW sizes. On real flow GF's informed runs are a random
subset of real runs whose signs are not independent of the rest of the flow, so
GF's exact zeta_K = 1 does not carry over (a Phase 3d caveat).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import clip3b as CL
import gen3 as G
import mbmodel as M


@dataclass
class Panel:
    eps: np.ndarray  # float +-1
    q: np.ndarray  # raw sizes (the impact's)
    q_mom: np.ndarray  # sizes the moments see (clipped)
    lm: np.ndarray
    starts: np.ndarray
    lengths: np.ndarray
    n_days: int

    @property
    def Nsd(self) -> np.ndarray:
        return self.lengths.reshape(-1, self.n_days)


_SMOOTH = np.array(sorted({a * 2**k for a in (1, 3, 5, 9, 15) for k in range(1, 40) if a * 2**k < 2**40}))


def embed_size(N: np.ndarray) -> np.ndarray:
    """Circulant / FFT size for series of length N: the smallest even 5-smooth number of the form
    a 2^k (a in 1, 3, 5, 9, 15) that is at least 2 N (3c-4; was the next power of 2, up to 33 % larger)."""
    need = np.maximum(2 * np.asarray(N, dtype=np.int64), 2)
    return _SMOOTH[np.searchsorted(_SMOOTH, need)].astype(np.int64)


def kernel_spectrum(m: int, beta: float) -> np.ndarray:
    """rfft of the trade kernel (1 + l)^-beta on lags 0 .. m/2 - 1 (0 at lag 0), zero-padded to m."""
    g = np.zeros(m)
    g[: m // 2] = M.trade_kernel(m // 2, beta)
    return np.fft.rfft(g)


def _groups(lengths: np.ndarray):
    Ms = embed_size(lengths)
    for m in np.unique(Ms):
        yield int(m), np.flatnonzero(Ms == m)


def gaussian_batch(lengths: np.ndarray, acf, rng: np.random.Generator) -> list[np.ndarray]:
    """One unit-variance Gaussian series per length, ACF acf(l) for l >= 1 (vectorised callable)."""
    out = [None] * len(lengths)
    for m, idx in _groups(lengths):
        h = m // 2
        r = np.empty(h + 1)
        r[0] = 1.0
        r[1:] = acf(np.arange(1, h + 1))
        lam = np.fft.rfft(np.concatenate([r, r[h - 1 : 0 : -1]])).real
        if lam.min() < -1e-9 * lam.max():
            raise RuntimeError("circulant embedding not non-negative definite")
        lam = np.sqrt(np.clip(lam, 0.0, None))
        z = rng.standard_normal((len(idx), h + 1)) + 1j * rng.standard_normal((len(idx), h + 1))
        z[:, 0] = z[:, 0].real * np.sqrt(2)
        z[:, -1] = z[:, -1].real * np.sqrt(2)
        x = np.fft.irfft(lam[None, :] * z, m, axis=1) * np.sqrt(m / 2)
        for row, i in enumerate(idx):
            out[i] = x[row, : lengths[i]]
    return out


def signs(lengths: np.ndarray, model: str, sm: dict, rng: np.random.Generator) -> np.ndarray:
    if model == "fgn_ar1":
        g = gaussian_batch(lengths, lambda l: G.sign_model_r(l, sm["w1"], sm["H"], sm["w2"], sm["phi"]), rng)
        x = np.concatenate(g)
    elif model == "two_component":
        g = gaussian_batch(lengths, lambda l: G.fgn_rho(l, sm["H"]), rng)
        m = sm["tau"] * rng.standard_normal(len(lengths))
        x = np.concatenate(g) + np.repeat(m, lengths)
    else:
        raise ValueError(f"unknown sign model {model!r}")
    return np.where(x >= 0, 1.0, -1.0)


def batched_causal_conv(x: np.ndarray, starts: np.ndarray, lengths: np.ndarray, beta: float) -> np.ndarray:
    """Per stock-day y_n = sum_{k<n} (1 + n - k)^-beta x_k (gen3.causal_conv with mbmodel.trade_kernel)."""
    y = np.zeros(len(x))
    for m, idx in _groups(lengths):
        X = np.zeros((len(idx), m))
        for row, i in enumerate(idx):
            X[row, : lengths[i]] = x[starts[i] : starts[i] + lengths[i]]
        Y = np.fft.irfft(np.fft.rfft(X, axis=1) * kernel_spectrum(m, beta)[None, :], m, axis=1)
        for row, i in enumerate(idx):
            y[starts[i] : starts[i] + lengths[i]] = Y[row, : lengths[i]]
    return y


def batched_multi_conv(xs: list[np.ndarray], betas: list[float], starts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """sum_i batched_causal_conv(xs[i], ..., betas[i]) with one inverse FFT per group: the forward transforms
    of every (x_i, kernel_i) pair are multiplied and summed in frequency, then inverted once."""
    y = np.zeros(len(xs[0]))
    for m, idx in _groups(lengths):
        acc = None
        for x, beta in zip(xs, betas, strict=True):
            X = np.zeros((len(idx), m))
            for row, i in enumerate(idx):
                X[row, : lengths[i]] = x[starts[i] : starts[i] + lengths[i]]
            F = np.fft.rfft(X, axis=1) * kernel_spectrum(m, beta)[None, :]
            acc = F if acc is None else acc + F
        Y = np.fft.irfft(acc, m, axis=1)
        for row, i in enumerate(idx):
            y[starts[i] : starts[i] + lengths[i]] = Y[row, : lengths[i]]
    return y


def day_cumsum(x: np.ndarray, starts: np.ndarray, lengths: np.ndarray, exclusive: bool = False) -> np.ndarray:
    """Cumulative sum restarting at every stock-day; exclusive = sum of the strictly earlier runs."""
    c = np.cumsum(x)
    base = np.repeat(c[starts] - x[starts], lengths)
    out = c - base
    return out - x if exclusive else out


def random_walk(starts, lengths, noise: float, rng) -> np.ndarray:
    """Per stock-day W_0 = 0, W_n = sum of n unit normals, times noise."""
    e = rng.standard_normal(int(lengths.sum()))
    return noise * day_cumsum(e, starts, lengths, exclusive=True)


def response(gen: str, cell: G.Cell, eps, q, starts, lengths, rng, params: dict) -> np.ndarray:
    """The synthetic mid of G0, GT or GF given the flow (eps, raw q)."""
    p = dict(params)
    psi = p.get("psi", 0.5)
    noise = p.get("noise", 1.0)
    x = eps * q**psi / G.mean_q_psi(cell, psi, "lognormal")
    if gen == "G0":
        m = batched_causal_conv(x, starts, lengths, p.get("beta", G.g0_beta(cell)))
    elif gen == "GT":
        inner, centres = G.size_bins(cell)
        betas = G.tent_betas(centres, p["beta_max"], p["c"], p["l_star"])
        b = np.searchsorted(inner, np.log(q), side="right")
        xs, bs = [], []
        for beta in np.unique(betas):
            sel = np.isin(b, np.flatnonzero(betas == beta))
            xs.append(np.where(sel, x, 0.0) if not sel.all() else x)
            bs.append(float(beta))
        if len(xs) == 1:
            m = batched_causal_conv(xs[0], starts, lengths, bs[0])
        else:
            m = batched_multi_conv(xs, bs, starts, lengths)
    elif gen == "GF":
        informed = p["informed"]
        m = day_cumsum(np.where(informed, x, 0.0) * p.get("kappa", 1.0), starts, lengths, exclusive=True)
    else:
        raise ValueError(f"unknown generator {gen!r}")
    return m + (random_walk(starts, lengths, noise, rng) if noise > 0 else 0.0)


def gt_params(cell: G.Cell) -> dict:
    import truth3 as TR

    g = TR.GT_TRUTH
    return {"beta_max": g["beta_max"], "c": g["c"], "l_star": cell.mu + g["l_off"] * cell.sigma_l, "psi": g["psi"]}


def stock_caps(q: np.ndarray, starts, lengths, n_days: int, p: float = 0.999) -> np.ndarray:
    """Per run, the cap of its stock: the stock's p-quantile over all its days (R2)."""
    qs = np.split(q, np.cumsum(lengths)[:-1])
    stocks = np.repeat(np.arange(len(lengths) // n_days), n_days)
    caps = CL.stock_caps(qs, stocks, p)
    return np.repeat(np.array([caps[s] for s in stocks]), lengths)


def day_lengths(cell_name: str, n_stocks: int, n_days: int, rng) -> np.ndarray:
    """(n_stocks, n_days): each row the day lengths of a random exploration stock of the cell."""
    z = np.load(G.phase3_dir() / f"calibration_{cell_name}.npz")
    N, sym = z["N"], z["sym"]
    rows = np.array([N[sym == s] for s in np.unique(sym)])
    if rows.shape[1] != n_days:
        raise ValueError("every exploration stock must have n_days study days")
    return rows[rng.integers(0, len(rows), n_stocks)]


def make_panel(gen: str, cell: G.Cell, Nsd: np.ndarray, rng: np.random.Generator, sign_model: str = "fgn_ar1",
               sm: dict | None = None, clip: str | float | None = "R2", params: dict | None = None,
               burn: int | None = None, continuous: bool = False) -> Panel:
    """burn=None: fresh start (every stock-day its own series, no history; the 3c path).
    burn=B (3c-4 primary): independent stationary days, every stock-day its own series of B + N runs with
    the first B discarded, independent of the stock's other days (stock-level parameters are the cell's).
    burn=B, continuous=True (the 3c-3 path, kept as a stress test): each stock one series of B + the sum of
    its day lengths, cut into its days."""
    if burn is not None:
        if continuous:
            return _stationary_panel(gen, cell, Nsd, rng, sign_model, sm, clip, params, int(burn))
        return _independent_panel(gen, cell, Nsd, rng, sign_model, sm, clip, params, int(burn))
    lengths = Nsd.ravel().astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    n_days = Nsd.shape[1]
    p = dict(params or {})
    if gen == "GT" and not p:
        p = gt_params(cell)
    if gen == "GF":
        f_inf = p.pop("f_inf", G.GF_INFORMED)
        sm_h = sm if sm is not None else (G.gf_sign_model(cell.name, f_inf) if sign_model == "fgn_ar1" else None)
        herd = signs(lengths, sign_model, sm_h, rng)
        informed = rng.random(len(herd)) < f_inf
        u = np.where(rng.random(len(herd)) < 0.5, 1.0, -1.0)
        eps = np.where(informed, u, herd)
        p["informed"] = informed
    else:
        eps = signs(lengths, sign_model, cell.sign if sm is None else sm, rng)
    q = np.exp(cell.mu + cell.sigma_l * rng.standard_normal(len(eps)))
    lm = response(gen, cell, eps, q, starts, lengths, rng, p)
    if clip == "R2":
        q_mom = np.minimum(q, stock_caps(q, starts, lengths, n_days))
    elif clip is None:
        q_mom = q
    else:
        q_mom = np.minimum(q, float(clip))
    return Panel(eps, q, q_mom, lm, starts, lengths, n_days)


def semi_synthetic(eps, q, Nsd: np.ndarray, gen: str, cell: G.Cell, rng: np.random.Generator,
                   clip: str | None = "R2", params: dict | None = None, fresh: bool = False) -> Panel:
    """Real flow (eps, q in stock-major order, day lengths Nsd) with a synthetic G0, GT or GF mid; the flow
    is returned unchanged (the same arrays' values), the clip touches only q_mom.

    CHOICE (3c-3): the response of a stock is one convolution over its study days concatenated in calendar
    order as a single trade-time sequence, after a burn-in made by prepending a copy of the same days
    (discarded), so the impact of a run reaches later days. fresh=True is the 3c path (each day alone)."""
    eps = np.asarray(eps, dtype=float)
    q = np.asarray(q, dtype=float)
    lengths = Nsd.ravel().astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    p = dict(params or {})
    if gen == "GT" and not p:
        p = gt_params(cell)
    if fresh:
        if gen == "GF":
            p["informed"] = rng.random(len(eps)) < p.pop("f_inf", G.GF_INFORMED)
        lm = response(gen, cell, eps, q, starts, lengths, rng, p)
    else:
        per = Nsd.sum(axis=1).astype(np.int64)  # runs per stock
        s_st = np.concatenate([[0], np.cumsum(per)[:-1]])
        ee = np.concatenate([np.concatenate([eps[a : a + n], eps[a : a + n]]) for a, n in zip(s_st, per, strict=True)])
        qq = np.concatenate([np.concatenate([q[a : a + n], q[a : a + n]]) for a, n in zip(s_st, per, strict=True)])
        series = 2 * per
        sst = np.concatenate([[0], np.cumsum(series)[:-1]])
        if gen == "GF":
            p["informed"] = rng.random(len(ee)) < p.pop("f_inf", G.GF_INFORMED)
        full = response(gen, cell, ee, qq, sst, series, rng, p)
        lm = np.concatenate([full[b + n : b + 2 * n] for b, n in zip(sst, per, strict=True)])
    q_mom = np.minimum(q, stock_caps(q, starts, lengths, Nsd.shape[1])) if clip == "R2" else q.copy()
    return Panel(eps, q, q_mom, lm, starts, lengths, Nsd.shape[1])


# ---- Phase 3c-3: stationary days ------------------------------------------------------------------
#
# A synthetic stock's days are cut from ONE continuous trade-time series (signs, sizes, impact): length
# burn + sum of its day lengths, the first `burn` runs discarded, the rest cut into consecutive days. The
# impact of every run reaches all later runs of the stock, across day boundaries; windows still never
# cross a day boundary (they are formed per day by pooled.panel_moments). make_panel(..., burn=None) is
# the old fresh-start path (each day its own series).


def stationary_stock(gen: str, cell: G.Cell, lengths: np.ndarray, burn: int, rng: np.random.Generator,
                     sign_model: str = "fgn_ar1", sm: dict | None = None, params: dict | None = None):
    """(eps, q, lm) of one stock's days, cut from a single series after a burn-in of `burn` runs."""
    lengths = np.asarray(lengths, dtype=np.int64)
    total = int(burn + lengths.sum())
    L = np.array([total])
    S0 = np.array([0])
    p = dict(params or {})
    if gen == "GT" and not p:
        p = gt_params(cell)
    if gen == "GF":
        f_inf = p.pop("f_inf", G.GF_INFORMED)
        sm_h = sm if sm is not None else (G.gf_sign_model(cell.name, f_inf) if sign_model == "fgn_ar1" else None)
        herd = signs(L, sign_model, sm_h, rng)
        informed = rng.random(total) < f_inf
        eps = np.where(informed, np.where(rng.random(total) < 0.5, 1.0, -1.0), herd)
        p["informed"] = informed
    else:
        eps = signs(L, sign_model, cell.sign if sm is None else sm, rng)
    q = np.exp(cell.mu + cell.sigma_l * rng.standard_normal(total))
    lm = response(gen, cell, eps, q, S0, L, rng, p)
    return eps[burn:], q[burn:], lm[burn:]


def cut_days(x: np.ndarray, series_lengths: np.ndarray, burn: int) -> np.ndarray:
    """Concatenated per-stock series (each burn + its days) -> concatenated days with every burn-in dropped."""
    starts = np.concatenate([[0], np.cumsum(series_lengths)[:-1]])
    keep = np.ones(len(x), dtype=bool)
    for s0 in starts:
        keep[s0 : s0 + burn] = False
    return x[keep]


STOCK_GROUP = 32  # stocks generated together (bounds the memory of the burn-in series)


def _stationary_panel(gen, cell, Nsd, rng, sign_model, sm, clip, params, burn) -> Panel:
    n_days = Nsd.shape[1]
    lengths = Nsd.ravel().astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    parts = [_stationary_group(gen, cell, Nsd[i : i + STOCK_GROUP], rng, sign_model, sm, params, burn)
             for i in range(0, len(Nsd), STOCK_GROUP)]
    eps, q, lm = (np.concatenate([x[k] for x in parts]) for k in range(3))
    if clip == "R2":
        q_mom = np.minimum(q, stock_caps(q, starts, lengths, n_days))
    elif clip is None:
        q_mom = q
    else:
        q_mom = np.minimum(q, float(clip))
    return Panel(eps, q, q_mom, lm, starts, lengths, n_days)


def _stationary_group(gen, cell, Nsd, rng, sign_model, sm, params, burn):
    """(eps, q, lm) of a group of stocks, each one series of burn + its days, burn-in dropped."""
    series = burn + Nsd.sum(axis=1).astype(np.int64)  # one series per stock
    s_starts = np.concatenate([[0], np.cumsum(series)[:-1]])
    p = dict(params or {})
    if gen == "GT" and not p:
        p = gt_params(cell)
    total = int(series.sum())
    if gen == "GF":
        f_inf = p.pop("f_inf", G.GF_INFORMED)
        sm_h = sm if sm is not None else (G.gf_sign_model(cell.name, f_inf) if sign_model == "fgn_ar1" else None)
        herd = signs(series, sign_model, sm_h, rng)
        informed = rng.random(total) < f_inf
        eps = np.where(informed, np.where(rng.random(total) < 0.5, 1.0, -1.0), herd)
        p["informed"] = informed
    else:
        eps = signs(series, sign_model, cell.sign if sm is None else sm, rng)
    q = np.exp(cell.mu + cell.sigma_l * rng.standard_normal(total))
    lm = response(gen, cell, eps, q, s_starts, series, rng, p)
    return tuple(cut_days(x, series, burn) for x in (eps, q, lm))


DAY_GROUP = 64  # stock-days generated together on the independent-days path


def independent_series(Nsd: np.ndarray, burn: int) -> np.ndarray:
    """(S * D, 1) series lengths of the independent-days path: every stock-day gets its own burn-in."""
    return (burn + Nsd.reshape(-1, 1)).astype(np.int64)


def _independent_panel(gen, cell, Nsd, rng, sign_model, sm, clip, params, burn) -> Panel:
    n_days = Nsd.shape[1]
    lengths = Nsd.ravel().astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    ser = independent_series(Nsd, burn) - burn  # one row per stock-day: its own series
    parts = [_stationary_group(gen, cell, ser[i : i + DAY_GROUP], rng, sign_model, sm, params, burn)
             for i in range(0, len(ser), DAY_GROUP)]
    eps, q, lm = (np.concatenate([x[k] for x in parts]) for k in range(3))
    if clip == "R2":
        q_mom = np.minimum(q, stock_caps(q, starts, lengths, n_days))
    elif clip is None:
        q_mom = q
    else:
        q_mom = np.minimum(q, float(clip))
    return Panel(eps, q, q_mom, lm, starts, lengths, n_days)
