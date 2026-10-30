"""Generator for the Maitrier-Bouchaud model: correlated metaorders with
square-root, slowly decaying impact.

References: G. Maitrier & J.-P. Bouchaud, arXiv:2506.07711 (Paper I, Sections
3 to 6 and Appendix A) and G. Maitrier, G. Loeper & J.-P. Bouchaud,
arXiv:2509.05065 (Paper II, Section 3 and Algorithm 1). "I (n)" is Paper I's
equation n, "II (n)" Paper II's.

One day, following Paper II, Section 3:

  Metaorders start as a Poisson process of rate nu on [-warmup, n_meta / nu)
    (Paper I, Sec. 3.1; II, Alg. 1 step 2). Start times are distinct and the
    metaorder id is its rank in start order.
  Child size q, constant within a metaorder: log q ~ N(m, sigma_l^2), so m is
    the mode (and mean) of log q (I (27)), truncated below q = 1 by redrawing
    (II, Sec. 3).
  Number of children s ~ Psi_q, P(s >= l) = l^(-mu_q) for integer l >= 1, the
    discrete form of I (14) with s0 = 1, and mu_q = mu_1 + lambda log q (I (21)),
    written here as mu_q = mu_m + lambda (log q - m) with mu_m = mu_1 + lambda m.
    Every mu_q must exceed 1 (finite mean duration, I Sec. 3.2).
  Children execute at rate phi: the first at the start, then Exp(phi) gaps
    (I Sec. 3.1; II, Alg. 1 steps 8 and 9). Children after the end of the day
    are dropped.
  Sign: +-1 with probability 1/2, or with cross-metaorder correlation
    E[eps_i eps_{i+k}] ~ Gamma k^(-gamma_x) in metaorder rank (II (1), I Sec. 4.3),
    built as eps_i = sign(fGn_i) with probability p and an independent random
    sign otherwise, H = 1 - gamma_x / 2. Then E[eps_i eps_{i+k}] =
    p^2 (2/pi) arcsin rho_H(k) -> Gamma k^(-gamma_x) with
    Gamma = p^2 (2/pi) H (2H - 1), which bounds Gamma (0.178 at gamma_x = 0.6).
  Impact of child j of a metaorder that started at t0_j, executed at t_j, on
    the mid at t > t_j, the generalized propagator I (49) = II (2):

        G = theta sqrt(q_j) (phi (t_j - t0_j) + n0)^(beta_j - 1/2) (tau0 / (t - t_j + tau0))^beta_j,
        beta_q = beta_1 - lambda' log q (I (62)) = beta_m - lambda' (log q - m),

    clipped into beta_clip (II, Sec. 3: "we ensure that 0 < beta_q < 1"; I (49)
    asks beta < 1/2, II allows up to 1). The mid read just before trade k is
    m_k = sum_{j < k} G_j(t_k) (II, Alg. 1 line 16), so m_0 = 0.
  Informed component (I Sec. 6.2 and (A8)): a fundamental move
    dF_i = (sigma_F / sqrt(nu)) (rho eps_i + sqrt(1 - rho^2) eta_i) per metaorder,
    eta_i ~ N(0, 1), added to the mid after the metaorder's first child. Over a
    window holding N = nu T metaorders the normalized move xi = sum dF /
    (sigma_F sqrt(T)) has E[eps_i xi] = rho / sqrt(nu T), as in Sec. 6.2, and
    sigma_F^2 is the fundamental variance per unit time.

Clock of the decay kernel. clock="physical" is II Algorithm 1 literally, with
tau0 the realized mean time between trades; it costs O(N^2). clock="trade"
(default) puts t_k - t_j = (k - j) tau0, so the kernel is (1 + k - j)^(-beta_j)
and m is a sum of convolutions: (1 + l)^(-beta) is interpolated in beta on 20
Chebyshev nodes spanning the day's betas (barycentric Lagrange weights,
relative error below 1e-9 for lags up to 10^6), and each node's convolution
is one FFT. The per-trade amplitude keeps the physical t_j - t0_j.

PAPER2 is Paper II's realistic configuration C-VD-VF (II Sec. 3 and Table 2):
mu_m = 1.5, beta_m = 0.25, lambda sigma_l^2 = 1/8, lambda' = 2 lambda,
sigma_l = 1, n0 = 3, m = 3, nu = 1.5e-3, phi = 2e-3. Paper II does not list
gamma_x or Gamma for its Figs. 2 and 6 to 8; Gamma = 0.1 is its Fig. 3 value
and gamma_x = 0.6 is Paper I's choice (Appendix A). configuration(name) gives
the other rows of II Table 2 (NC: Gamma = 0; NVD: lambda = lambda' = 0;
NVF: q = 1).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

import synth


@dataclass(frozen=True)
class MBParams:
    nu: float = 1.5e-3  # metaorder starts per unit time
    phi: float = 2e-3  # child executions per unit time within a metaorder
    m: float = 3.0  # mode of log q
    sigma_l: float = 1.0  # width of log q
    mu_m: float = 1.5  # duration tail exponent at log q = m
    lam: float = 1 / 8  # mu_q = mu_m + lam (log q - m)
    beta_m: float = 0.25  # impact decay exponent at log q = m
    lam_p: float = 1 / 4  # beta_q = beta_m - lam_p (log q - m)
    beta_clip: tuple[float, float] = (0.0, 0.99)
    n0: float = 3.0  # detection lag in child orders
    theta: float = 1.0  # theta(q) = theta sqrt(q)
    Gamma: float = 0.1  # cross-metaorder sign correlation amplitude
    gamma_x: float = 0.6  # and its decay exponent in metaorder rank
    sigma_F: float = 0.0  # fundamental volatility per unit time
    rho: float = 0.0  # information of a metaorder sign about its fundamental move

    @property
    def mu_1(self) -> float:
        return self.mu_m - self.lam * self.m

    @property
    def beta_1(self) -> float:
        return self.beta_m + self.lam_p * self.m


PAPER2 = MBParams()


def configuration(name: str, base: MBParams = PAPER2) -> MBParams:
    """A row of Paper II Table 2: [N]C-[N]VD-[N]VF."""
    parts = name.split("-")
    if len(parts) != 3 or parts[0] not in ("C", "NC") or parts[1] not in ("VD", "NVD") or parts[2] not in ("VF", "NVF"):
        raise ValueError(f"unknown configuration {name!r}")
    p = base
    if parts[0] == "NC":
        p = replace(p, Gamma=0.0)
    if parts[1] == "NVD":
        p = replace(p, lam=0.0, lam_p=0.0)
    if parts[2] == "NVF":
        p = replace(p, m=0.0, sigma_l=0.0)
    return p


@dataclass
class MBDay:
    """One simulated day. Trade arrays in trade order; meta_* indexed by metaorder id."""

    t: np.ndarray  # execution time
    eps: np.ndarray  # sign (int8)
    q: np.ndarray  # child size
    meta: np.ndarray  # metaorder id
    rank: np.ndarray  # child rank within its metaorder, 1 = first
    mid: np.ndarray  # mid just before the trade: impact + fundamental
    impact: np.ndarray  # impact part of mid
    meta_t0: np.ndarray
    meta_eps: np.ndarray
    meta_q: np.ndarray
    meta_n: np.ndarray  # children drawn (before the end-of-day cut)
    meta_mu: np.ndarray
    meta_beta: np.ndarray
    meta_dF: np.ndarray
    tau0: float  # realized mean time between trades
    t_end: float


def metaorder_signs(n: int, Gamma: float, gamma_x: float, rng: np.random.Generator) -> np.ndarray:
    """Signs with E[eps_i eps_{i+k}] = p^2 (2/pi) arcsin rho_H(k), H = 1 - gamma_x/2 (see module doc)."""
    if Gamma == 0:
        return rng.choice(np.array([-1, 1], dtype=np.int8), n)
    H = 1 - gamma_x / 2
    gmax = (2 / np.pi) * H * (2 * H - 1)
    if not 0 < gamma_x < 1 or not 0 < Gamma <= gmax:
        raise ValueError(f"need 0 < gamma_x < 1 and 0 < Gamma <= {gmax:.4f}")
    p = np.sqrt(Gamma / gmax)
    s = synth.fgn_signs(n, H, rng)
    r = rng.choice(np.array([-1, 1], dtype=np.int8), n)
    return np.where(rng.random(n) < p, s, r).astype(np.int8)


def sign_correlation(k: np.ndarray, Gamma: float, gamma_x: float) -> np.ndarray:
    """Exact E[eps_i eps_{i+k}] of metaorder_signs, k >= 1."""
    if Gamma == 0:
        return np.zeros(len(np.atleast_1d(k)))
    H = 1 - gamma_x / 2
    return Gamma / ((2 / np.pi) * H * (2 * H - 1)) * synth.fgn_sign_acf(k, H)


def child_sizes(n: int, p: MBParams, rng: np.random.Generator) -> np.ndarray:
    """q with log q ~ N(m, sigma_l^2) truncated at q >= 1."""
    ell = p.m + p.sigma_l * rng.standard_normal(n)
    bad = ell < 0
    while np.any(bad):
        ell[bad] = p.m + p.sigma_l * rng.standard_normal(int(bad.sum()))
        bad = ell < 0
    return np.exp(ell)


def exponents(q: np.ndarray, p: MBParams) -> tuple[np.ndarray, np.ndarray]:
    """mu_q (I (21)) and beta_q (I (62), clipped)."""
    ell = np.log(q)
    mu = p.mu_m + p.lam * (ell - p.m)
    beta = np.clip(p.beta_m - p.lam_p * (ell - p.m), *p.beta_clip)
    if np.any(mu <= 1):
        raise ValueError("mu_q <= 1 for some q: infinite mean duration")
    return mu, beta


def _cheb_weights(beta: np.ndarray, lo: float, hi: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Chebyshev nodes on [lo, hi] and the Lagrange basis l_c(beta_j), shape (n, len(beta))."""
    k = np.arange(n)
    nodes = (lo + hi) / 2 + (hi - lo) / 2 * np.cos((2 * k + 1) * np.pi / (2 * n))
    w = (-1.0) ** k * np.sin((2 * k + 1) * np.pi / (2 * n))
    d = beta[None, :] - nodes[:, None]
    exact = d == 0
    d[exact] = 1.0
    terms = w[:, None] / d
    L = terms / terms.sum(axis=0)
    hit = exact.any(axis=0)
    L[:, hit] = exact[:, hit].astype(float)
    return nodes, L


def _causal_conv(x: np.ndarray, g: np.ndarray) -> np.ndarray:
    n = len(x)
    size = 1 << int(np.ceil(np.log2(2 * n)))
    return np.fft.irfft(np.fft.rfft(x, size) * np.fft.rfft(g, size), size)[:n]


def trade_kernel(n: int, beta: float) -> np.ndarray:
    """(tau0 / (l tau0 + tau0))^beta = (1 + l)^-beta at lags l = 0..n-1, with 0 at l = 0:
    the mid read before trade k excludes trade k itself."""
    g = (1.0 + np.arange(n, dtype=float)) ** -beta
    g[0] = 0.0
    return g


def impact_trade_clock(amp: np.ndarray, beta: np.ndarray, n_nodes: int = 20) -> np.ndarray:
    """m_k = sum_{j<k} amp_j (1 + k - j)^(-beta_j)."""
    n = len(amp)
    lo, hi = float(beta.min()), float(beta.max())
    if hi - lo < 1e-12:
        return _causal_conv(amp, trade_kernel(n, lo))
    nodes, L = _cheb_weights(beta, lo, hi, n_nodes)
    out = np.zeros(n)
    for c, b in enumerate(nodes):
        out += _causal_conv(amp * L[c], trade_kernel(n, b))
    return out


def impact_physical_clock(amp: np.ndarray, beta: np.ndarray, t: np.ndarray, tau0: float,
                          chunk: int = 64) -> np.ndarray:
    """m_k = sum_{j<k} amp_j (tau0 / (t_k - t_j + tau0))^beta_j, II Algorithm 1 line 16, O(N^2)."""
    n = len(amp)
    out = np.zeros(n)
    for a in range(1, n, chunk):
        b = min(n, a + chunk)
        dt = t[a:b, None] - t[None, :b - 1]
        k = (tau0 / (np.maximum(dt, 0.0) + tau0)) ** beta[None, :b - 1]
        k[np.arange(a, b)[:, None] <= np.arange(b - 1)[None, :]] = 0.0  # only j < k
        out[a:b] = k @ amp[: b - 1]
    return out


def simulate_day(p: MBParams, n_meta: int, rng: np.random.Generator, warmup: float = 0.0,
                 clock: str = "trade") -> MBDay:
    """One day of n_meta / nu time units; trades before time 0 (warmup) are dropped from
    the output, but their impact stays in the mid."""
    t_end = n_meta / p.nu
    # metaorders
    t0 = -warmup + np.cumsum(rng.exponential(1 / p.nu, int(1.2 * p.nu * (t_end + warmup)) + 100))
    while t0[-1] < t_end:
        t0 = np.concatenate([t0, t0[-1] + np.cumsum(rng.exponential(1 / p.nu, 1000))])
    t0 = t0[t0 < t_end]
    if len(np.unique(t0)) != len(t0):
        raise RuntimeError("tied metaorder start times")
    n_m = len(t0)
    q = child_sizes(n_m, p, rng)
    mu, beta = exponents(q, p)
    n_draw = np.floor(rng.random(n_m) ** (-1.0 / mu))
    n_draw = np.minimum(n_draw, 1e15).astype(np.int64)
    sgn = metaorder_signs(n_m, p.Gamma, p.gamma_x, rng)
    eta = rng.standard_normal(n_m)
    dF = p.sigma_F / np.sqrt(p.nu) * (p.rho * sgn + np.sqrt(1 - p.rho**2) * eta)
    # children: enough Exp(phi) gaps to pass the end of the day with overwhelming probability
    room = p.phi * (t_end - t0)
    n_gen = np.minimum(n_draw, np.ceil(room + 10 * np.sqrt(room) + 20).astype(np.int64))
    start = np.concatenate([[0], np.cumsum(n_gen)[:-1]])
    owner = np.repeat(np.arange(n_m), n_gen)
    gaps = rng.exponential(1 / p.phi, int(n_gen.sum()))
    gaps[start] = 0.0
    cs = np.cumsum(gaps)
    tau = cs - cs[start][owner]
    t = t0[owner] + tau
    rank = np.arange(len(t)) - start[owner] + 1
    keep = t < t_end
    owner, t, tau, rank = owner[keep], t[keep], tau[keep], rank[keep]
    order = np.argsort(t, kind="stable")
    owner, t, tau, rank = owner[order], t[order], tau[order], rank[order]
    n = len(t)
    tau0 = (t[-1] - t[0]) / (n - 1)
    amp = p.theta * sgn[owner] * np.sqrt(q[owner]) * (p.phi * tau + p.n0) ** (beta[owner] - 0.5)
    if clock == "trade":
        impact = impact_trade_clock(amp, beta[owner])
    elif clock == "physical":
        impact = impact_physical_clock(amp, beta[owner], t, tau0)
    else:
        raise ValueError("clock is 'trade' or 'physical'")
    # fundamental moves land after each metaorder's first child
    F = np.zeros(n + 1)
    first = np.flatnonzero(rank == 1)
    np.add.at(F, first + 1, dF[owner[first]])
    F = np.cumsum(F)[:n]
    s = t >= 0
    return MBDay(t[s], sgn[owner][s], q[owner][s], owner[s], rank[s], (impact + F)[s], impact[s],
                 t0, sgn, q, n_draw, mu, beta, dF, tau0, t_end)


def simulate(p: MBParams, n_days: int, n_meta: int, seed: int, **kw) -> list[MBDay]:
    """Independent days from one seed."""
    return [simulate_day(p, n_meta, np.random.default_rng([seed, d]), **kw) for d in range(n_days)]
