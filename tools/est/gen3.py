"""Phase 3a generators: synthetic panels calibrated to a v2 cell.

A panel is n_stocks x n_days stock-days; each day d has its own length N drawn
from the cell's empirical distribution of runs per stock-day (resampled from
the calibration's stock-day N values), and gives the arrays that
moments.day_moments takes: eps (sign), q (size), lm (mid read just before each
run), ok (all True).

Ingredients shared by G0 and GT
  signs   eps_t = sign(g_t), g a stationary unit-variance Gaussian series with
          ACF r(l) = w1 rho_H(l) + w2 phi^l (l >= 1), rho_H the fGn
          autocorrelation; by the arcsine law C(l) = (2/pi) arcsin r(l). The
          four parameters are fitted to the cell's measured C(l) on its lags
          1..L_max (fit_sign_model, least squares on the relative error of r);
          w1 + w2 <= 1, the rest of the variance is white. Asymptotically C(l) ~
          (2/pi) w1 H (2H - 1) l^-(2 - 2H): gamma_gen = 2 - 2H. One circulant
          embedding generates g (gaussian_acf).
  sizes   log q = mu + sigma_l z with the cell's (mu, sigma_l), z iid N(0, 1)
          ("lognormal"), or q resampled from the cell's pooled empirical sizes
          ("empirical"). copula=True replaces the iid z by a Gaussian series with
          ACF w rho_H'(l) (l >= 1), (w, H') fitted to the cell's log q ACF
          (fit_size_copula); with "empirical" the day's resampled sizes are
          sorted and assigned by the rank of z (a Gaussian copula with the
          empirical margin; its log q ACF is close to, not exactly, w rho_H').
  mid     m_n = theta sum_{k<n} G(n - k) x_k + noise W_n, x_k = eps_k q_k^psi,
          G(l) = (1 + l)^-beta for l >= 1 and G(0) = 0 (mbmodel.trade_kernel),
          theta = 1 / E[q^psi] under the size model, W a Gaussian random walk
          with unit steps times `noise` (W_0 = 0).

Families
  G0  "monotone null": one size-independent kernel, default beta = (1 -
      gamma_gen) / 2, the exponent that keeps the price diffusive under the
      generator's own sign memory (Bouchaud et al. 2004); psi = 0.5 default.
  GT  "tent oracle", a CHOICE (power device, not a mechanism): log q is cut into
      N_BINS bins of equal width on [mu - 3 sigma_l, mu + 5 sigma_l] (outer bins
      open; the top edge is about the largest log q of a 10^6-run day, where the
      q^(a + psi) weight of a = 3 sits: a range ending at mu + 3 sigma_l, the first
      version, left that weight in one open bin and never tested the tent there),
      and bin b with centre l_b decays with
          beta_b = max(BETA_FLOOR, beta_max - c |l_b - l_star|);
      m_n = theta sum_b sum_{k<n} G_b(n - k) x_k 1[q_k in b]. Bins with equal
      beta share one FFT convolution, so c = 0 gives G0's arrays exactly.
  GF  "informed flow" (P-F world), run level: each run is informed with
      probability f_inf = GF_INFORMED, independently of everything else. An
      informed run's sign is a fair coin u_t independent of all other signs and
      sizes; an uninformed run's sign is sign(g_t), g the cell's sign model
      refitted to C(l) / (1 - f_inf)^2, so that the mixture's C(l) =
      (1 - f_inf)^2 (2/pi) arcsin r_g(l) is the cell's. Only informed runs move
      the price, permanently: m_n = kappa sum_{k<n} 1_k eps_k q_k^psi / E[q^psi]
      + noise W_n. Then E[Delta_T(n) I^a_T(n)] = T f_inf kappa E[q^(psi+a)] /
      E[q^psi] exactly, since a pair of different runs always involves a u
      independent of the other factor: K_a(T) is linear in T and zeta_K(a) = 1
      for every a at every T, while the flow keeps the cell's sign memory.
  GF_mb  mbmodel's existing informed component (Paper I Sec. 6.2): PAPER2
      metaorders with no mechanical impact (theta = 0), independent metaorder
      signs (Gamma = 0) and a permanent fundamental move per metaorder,
      dF_i = (sigma_F / sqrt(nu)) (rho eps_i + sqrt(1 - rho^2) eta_i), landing after
      its first child; GF_MB_PARAMS lists the values. Its own signs and sizes, not
      the cell's; only day lengths come from the cell. MEASURED: it is not a P-F
      world at finite T. A metaorder's children straddle window edges, so
      K(T) = T c_tot - offset(T) with offset / T ~ T^(1 - mu) (mu = 1.125 to 1.9
      here): zeta_K ~ 1.18 on [8, 64] and 1.08 on [64, 1024], nearly flat in a.
  GM  "Paper II model": mbmodel.PAPER2 ("GM"), and the beta-floor variants
      beta_clip = (-0.6, 0.99) with Gamma = 0 ("GM_floor_G0") and Gamma = 0.1
      ("GM_floor_G01"). Days are the first N trades of an mbmodel day.

simulate_day(gen, cell, N, rng, **params) draws one stock-day of length N;
panel(gen, cell, n_stocks, n_days, seed, **params) draws the panel with
stock-day (i, d) from default_rng([seed, i, d]).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, replace
from functools import cache
from pathlib import Path

import numpy as np

import mbmodel as M

N_BINS = 12
BIN_RANGE = (-3.0, 5.0)  # bins span mu + BIN_RANGE sigma_l
BETA_FLOOR = 0.0
GENERATORS = ("G0", "GT", "GF", "GF_mb", "GM", "GM_floor_G0", "GM_floor_G01")
GF_INFORMED = 0.1
GF_MB_PARAMS = replace(M.PAPER2, theta=0.0, Gamma=0.0, sigma_F=1.0, rho=0.5)
GM_PARAMS = {"GM": M.PAPER2,
             "GM_floor_G0": replace(M.PAPER2, beta_clip=(-0.6, 0.99), Gamma=0.0),
             "GM_floor_G01": replace(M.PAPER2, beta_clip=(-0.6, 0.99), Gamma=0.1)}


# ---- Gaussian series --------------------------------------------------------------


def fgn_rho(l, H: float) -> np.ndarray:
    """fGn autocorrelation at lags l (any sign)."""
    l = np.abs(np.asarray(l, dtype=float))
    return 0.5 * ((l + 1) ** (2 * H) - 2 * l ** (2 * H) + np.abs(l - 1) ** (2 * H))


def gaussian_acf(n: int, r: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Stationary unit-variance Gaussian series of length n with ACF r(0..n), r(0) = 1, by circulant
    embedding (Davies-Harte); exact when the embedding is non-negative definite (checked)."""
    c = np.concatenate([r[: n + 1], r[n - 1 : 0 : -1]])
    lam = np.fft.rfft(c).real
    if lam.min() < -1e-9 * lam.max():
        raise RuntimeError("circulant embedding not non-negative definite")
    lam = np.clip(lam, 0.0, None)
    m = len(c)
    z = rng.standard_normal(m // 2 + 1) + 1j * rng.standard_normal(m // 2 + 1)
    z[0] = z[0].real * np.sqrt(2)
    if m % 2 == 0:
        z[-1] = z[-1].real * np.sqrt(2)
    return np.fft.irfft(np.sqrt(lam) * z, m)[:n] * np.sqrt(m / 2)


def sign_model_r(l, w1: float, H: float, w2: float, phi: float) -> np.ndarray:
    """Gaussian ACF of the sign model at lags l >= 1."""
    l = np.asarray(l, dtype=float)
    return w1 * fgn_rho(l, H) + w2 * phi**l


def sign_model_C(l, w1: float, H: float, w2: float, phi: float) -> np.ndarray:
    """C(l) of eps = sign(g) (arcsine law)."""
    return (2 / np.pi) * np.arcsin(np.clip(sign_model_r(l, w1, H, w2, phi), -1, 1))


def _nn2(A: np.ndarray, y: np.ndarray, wt: np.ndarray) -> np.ndarray:
    """Weighted least squares for two coefficients with both >= 0 and their sum <= 1."""
    Aw, yw = A * wt[:, None], y * wt
    best, cand = None, []
    try:
        cand.append(np.linalg.lstsq(Aw, yw, rcond=None)[0])
    except np.linalg.LinAlgError:
        pass
    for j in range(2):
        a = Aw[:, j]
        cand.append(np.eye(2)[j] * max(0.0, float(a @ yw / (a @ a))))
    d = Aw[:, 0] - Aw[:, 1]
    t = float(np.clip(d @ (yw - Aw[:, 1]) / (d @ d), 0, 1)) if d @ d > 0 else 0.0
    cand.append(np.array([t, 1 - t]))
    for x in cand:
        x = np.clip(x, 0, None)
        if x.sum() > 1:
            x = x / x.sum()
        rss = float(np.sum((Aw @ x - yw) ** 2))
        if best is None or rss < best[0]:
            best = (rss, x)
    return best[1]


def fit_sign_model(lags: np.ndarray, C: np.ndarray, H_grid=None, phi_grid=None) -> dict:
    """(w1, H, w2, phi) minimising sum over lags with C > 0 of ((r_model - r) / r)^2,
    r = sin(pi C / 2) (the Gaussian ACF that sign(g) needs); H and phi on grids,
    (w1, w2) by constrained least squares at each grid point."""
    lags = np.asarray(lags, dtype=float)
    C = np.asarray(C, dtype=float)
    use = np.isfinite(C) & (C > 0) & (lags >= 1)
    l, r = lags[use], np.sin(np.pi * C[use] / 2)
    H_grid = np.arange(0.51, 0.996, 0.005) if H_grid is None else H_grid
    phi_grid = np.arange(0.02, 0.981, 0.02) if phi_grid is None else phi_grid
    wt = 1 / r
    best = (np.inf, None)
    for H in H_grid:
        f1 = fgn_rho(l, H)
        for phi in phi_grid:
            A = np.stack([f1, phi**l], axis=1)
            x = _nn2(A, r, wt)
            rss = float(np.sum(((A @ x - r) * wt) ** 2))
            if rss < best[0]:
                best = (rss, (float(x[0]), float(H), float(x[1]), float(phi)))
    w1, H, w2, phi = best[1]
    return {"w1": w1, "H": H, "w2": w2, "phi": phi, "rmse": float(np.sqrt(best[0] / len(l)))}


def fit_size_copula(acf: np.ndarray, lags: np.ndarray, H_grid=None) -> tuple[float, float]:
    """(w, H) of rho(l) = w rho_H(l), least squares on the log-spaced subset of lags with acf > 0."""
    lags = np.asarray(lags, dtype=float)
    acf = np.asarray(acf, dtype=float)
    pick = np.unique(np.round(np.logspace(0, np.log10(lags.max()), 25)))
    use = np.isin(lags, pick) & (acf > 0)
    l, y = lags[use], acf[use]
    H_grid = np.arange(0.51, 0.996, 0.005) if H_grid is None else H_grid
    best = (np.inf, 0.0, 0.5)
    for H in H_grid:
        f = fgn_rho(l, H) / y
        w = float(np.clip(f @ np.ones_like(y) / (f @ f), 0, 1))
        rss = float(np.sum((w * f - 1) ** 2))
        if rss < best[0]:
            best = (rss, w, float(H))
    return best[1], best[2]


def signs(n: int, sm: dict, rng: np.random.Generator) -> np.ndarray:
    """sign(g), g with the fGn + AR(1) ACF of sm; or, when sm has "tau" (the 3b two-component model),
    sign(X + m), X = fGn(H), m ~ N(0, tau^2) drawn once for the day."""
    r = np.empty(n + 1)
    r[0] = 1.0
    if "tau" in sm:
        r[1:] = fgn_rho(np.arange(1, n + 1), sm["H"])
        x = gaussian_acf(n, r, rng) + sm["tau"] * rng.standard_normal()
        return np.where(x >= 0, 1, -1).astype(np.int8)
    r[1:] = sign_model_r(np.arange(1, n + 1), sm["w1"], sm["H"], sm["w2"], sm["phi"])
    return np.where(gaussian_acf(n, r, rng) >= 0, 1, -1).astype(np.int8)


# ---- cells ------------------------------------------------------------------------


def phase3_dir() -> Path:
    import explore  # the data directory lives with the guard

    return explore.data_dir() / "phase3"


@dataclass(frozen=True)
class Cell:
    name: str
    Ns: np.ndarray  # runs per stock-day, one per calibration stock-day
    mu: float
    sigma_l: float
    copula: tuple[float, float]  # (w, H) of the log q ACF
    sign: dict  # w1, H, w2, phi

    @property
    def gamma_gen(self) -> float:
        return 2 - 2 * self.sign["H"]

    def pool(self) -> np.ndarray:
        return _pool(self.name)


@cache
def _pool(name: str) -> np.ndarray:
    return np.load(phase3_dir() / f"sizes_{name}.npy").astype(float)


@cache
def load_cell(name: str) -> Cell:
    z = np.load(phase3_dir() / f"calibration_{name}.npz")
    row = json.loads(str(z["row"]))
    return Cell(name, z["N"].astype(int), float(row["mu_logq"]), float(row["sigma_l"]),
                (float(row["copula_w"]), float(row["copula_H"])),
                {k: float(row[f"sgn_{k}"]) for k in ("w1", "H", "w2", "phi")})


def draw_N(cell: Cell, rng: np.random.Generator) -> int:
    return int(cell.Ns[rng.integers(0, len(cell.Ns))])


def sizes(n: int, cell: Cell, rng: np.random.Generator, kind: str = "lognormal", copula: bool = False) -> np.ndarray:
    if copula:
        w, H = cell.copula
        r = np.empty(n + 1)
        r[0] = 1.0
        r[1:] = w * fgn_rho(np.arange(1, n + 1), H)
        z = gaussian_acf(n, r, rng)
    else:
        z = rng.standard_normal(n) if kind == "lognormal" else None
    if kind == "lognormal":
        return np.exp(cell.mu + cell.sigma_l * z)
    if kind == "empirical":
        pool = cell.pool()
        q = pool[rng.integers(0, len(pool), n)]
        if copula:
            q = np.sort(q)[np.argsort(np.argsort(z))]
        return q
    raise ValueError("kind is 'lognormal' or 'empirical'")


def mean_q_psi(cell: Cell, psi: float, kind: str) -> float:
    if kind == "lognormal":
        return float(np.exp(psi * cell.mu + 0.5 * (psi * cell.sigma_l) ** 2))
    return float(np.mean(cell.pool() ** psi))


# ---- impact -----------------------------------------------------------------------


def causal_conv(x: np.ndarray, g: np.ndarray) -> np.ndarray:
    """y_n = sum_{k<=n} g(n - k) x_k by FFT (g given at lags 0..n-1)."""
    return M._causal_conv(x, g)


def size_bins(cell: Cell, n_bins: int = N_BINS) -> tuple[np.ndarray, np.ndarray]:
    """Inner edges (n_bins - 1) and centres (n_bins) of the log q bins."""
    lo, hi = cell.mu + BIN_RANGE[0] * cell.sigma_l, cell.mu + BIN_RANGE[1] * cell.sigma_l
    edges = np.linspace(lo, hi, n_bins + 1)
    return edges[1:-1], 0.5 * (edges[:-1] + edges[1:])


def tent_betas(centres: np.ndarray, beta_max: float, c: float, l_star: float) -> np.ndarray:
    return np.maximum(BETA_FLOOR, beta_max - c * np.abs(centres - l_star))


def binned_impact(x: np.ndarray, b: np.ndarray, betas: np.ndarray) -> np.ndarray:
    """m_n = sum_{k<n} (1 + n - k)^-beta_{b_k} x_k: one FFT convolution per distinct beta."""
    n = len(x)
    out = np.zeros(n)
    for beta in np.unique(betas):
        sel = np.isin(b, np.flatnonzero(betas == beta))
        out += causal_conv(np.where(sel, x, 0.0) if not sel.all() else x, M.trade_kernel(n, beta))
    return out


def direct_binned_impact(x: np.ndarray, b: np.ndarray, betas: np.ndarray) -> np.ndarray:
    """The same sum, term by term (O(n^2), for tests)."""
    n = len(x)
    out = np.zeros(n)
    for t in range(n):
        for k in range(t):
            out[t] += (1.0 + t - k) ** -betas[b[k]] * x[k]
    return out


@dataclass
class Day:
    eps: np.ndarray
    q: np.ndarray
    lm: np.ndarray
    ok: np.ndarray


def _propagator_day(cell: Cell, N: int, rng: np.random.Generator, betas_fn, psi: float, noise: float,
                    size_kind: str, copula: bool) -> Day:
    eps = signs(N, cell.sign, rng)
    q = sizes(N, cell, rng, size_kind, copula)
    x = eps * q**psi / mean_q_psi(cell, psi, size_kind)
    m = betas_fn(x, q)
    if noise > 0:
        m = m + noise * np.concatenate([[0.0], np.cumsum(rng.standard_normal(N - 1))])
    return Day(eps, q, m, np.ones(N, bool))


def g0_beta(cell: Cell) -> float:
    return (1 - cell.gamma_gen) / 2


def simulate_g0(cell: Cell, N: int, rng: np.random.Generator, beta: float | None = None, psi: float = 0.5,
                noise: float = 1.0, size_kind: str = "lognormal", copula: bool = False) -> Day:
    beta = g0_beta(cell) if beta is None else beta
    return _propagator_day(cell, N, rng, lambda x, q: causal_conv(x, M.trade_kernel(len(x), beta)),
                           psi, noise, size_kind, copula)


def simulate_gt(cell: Cell, N: int, rng: np.random.Generator, beta_max: float | None = None, c: float = 0.2,
                l_star: float | None = None, psi: float = 0.5, noise: float = 1.0, size_kind: str = "lognormal",
                copula: bool = False) -> Day:
    beta_max = g0_beta(cell) if beta_max is None else beta_max
    l_star = cell.mu + 2 * cell.sigma_l if l_star is None else l_star
    inner, centres = size_bins(cell)
    betas = tent_betas(centres, beta_max, c, l_star)

    def impact(x, q):
        return binned_impact(x, np.searchsorted(inner, np.log(q), side="right"), betas)

    return _propagator_day(cell, N, rng, impact, psi, noise, size_kind, copula)


@cache
def gf_sign_model(name: str, f_inf: float = GF_INFORMED) -> dict:
    """The cell's sign model refitted to C(l) / (1 - f_inf)^2 (the uninformed runs' share of the memory)."""
    z = np.load(phase3_dir() / f"calibration_{name}.npz")
    sm = fit_sign_model(z["sign_lags"], np.minimum(z["C"] / (1 - f_inf) ** 2, 0.999))
    return {k: sm[k] for k in ("w1", "H", "w2", "phi")}


def simulate_gf(cell: Cell, N: int, rng: np.random.Generator, f_inf: float = GF_INFORMED, kappa: float = 1.0,
                psi: float = 0.5, noise: float = 1.0, size_kind: str = "lognormal", copula: bool = False,
                sign_model: dict | None = None) -> Day:
    sm = gf_sign_model(cell.name, f_inf) if sign_model is None else sign_model
    herd = signs(N, sm, rng)
    informed = rng.random(N) < f_inf
    u = np.where(rng.random(N) < 0.5, 1, -1).astype(np.int8)
    eps = np.where(informed, u, herd).astype(np.int8)
    q = sizes(N, cell, rng, size_kind, copula)
    x = np.where(informed, eps * q**psi, 0.0) * (kappa / mean_q_psi(cell, psi, size_kind))
    m = np.concatenate([[0.0], np.cumsum(x[:-1])])
    if noise > 0:
        m = m + noise * np.concatenate([[0.0], np.cumsum(rng.standard_normal(N - 1))])
    return Day(eps, q, m, np.ones(N, bool))


@cache
def _trades_per_meta(p: M.MBParams) -> float:
    d = M.simulate_day(p, 4000, np.random.default_rng(0))
    return len(d.t) / 4000


def mb_cut(p: M.MBParams, N: int, rng: np.random.Generator) -> Day:
    """The first N trades of an mbmodel day long enough to hold them."""
    n_meta = int(1.25 * N / _trades_per_meta(p)) + 50
    while True:
        d = M.simulate_day(p, n_meta, rng)
        if len(d.t) >= N:
            break
        n_meta *= 2
    return Day(d.eps[:N].astype(np.int8), d.q[:N], d.mid[:N], np.ones(N, bool))


def simulate_day(gen: str, cell: Cell, N: int, rng: np.random.Generator, **params) -> Day:
    if gen == "G0":
        return simulate_g0(cell, N, rng, **params)
    if gen == "GT":
        return simulate_gt(cell, N, rng, **params)
    if gen == "GF":
        return simulate_gf(cell, N, rng, **params)
    if gen == "GF_mb":
        return mb_cut(replace(GF_MB_PARAMS, **params), N, rng)
    if gen in GM_PARAMS:
        return mb_cut(replace(GM_PARAMS[gen], **params), N, rng)
    raise ValueError(f"unknown generator {gen!r}")


def panel(gen: str, cell: Cell | str, n_stocks: int, n_days: int, seed: int, **params) -> list[list[Day]]:
    cell = load_cell(cell) if isinstance(cell, str) else cell
    out = []
    for i in range(n_stocks):
        row = []
        for d in range(n_days):
            rng = np.random.default_rng([seed, i, d])
            row.append(simulate_day(gen, cell, draw_N(cell, rng), rng, **params))
        out.append(row)
    return out
