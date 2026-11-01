"""Phase 3d-2 worlds: the panels on which the four primary tests' operating characteristics are measured.

Every panel holds the six confirmatory cells at their confirmation stock counts. Real data enter only as the
exploration stocks' signs, sizes and timestamps, read through explore.load_explore_runs (flow_cache keeps
exactly those three fields, never a price). Every price response here is synthetic.

Semi-synthetic panels (W-null, W0, W0-psi3, WB). Per confirmatory cell, the cell's exploration stocks are
resampled with replacement up to the cell's confirmation count; each drawn stock keeps its real 8 days of
signs, raw sizes and timestamps (a stock drawn twice appears twice with independent responses). The mid is
synthetic, on independent stationary days (3c-4): every stock-day is the end of its own series of
BURN + N runs, the first BURN discarded.
  CHOICE: on real flow the BURN runs before a stock-day are a contiguous cyclic stretch of the same stock's
  OTHER seven days (concatenated in calendar order, random start), so the history has the stock's own size
  and sign statistics and is independent of the day (the cross-day mean-sign correlation was about 0 in
  3c-4). The burn-in only matters for decaying kernels: under GF (permanent impact) it adds a constant to
  the day's mid, which no statistic sees, so GF skips it (exact).

Worlds (flow, response):
  W-null   real flow, sizes permuted within each stock-day (surrogates.s1), GF response (informed runs a
           random 10 %, permanent impact eps q^psi, psi = 0.5). Phase 3d-2b: the informed runs' signs are
           independent fair coins (as GF in 3a), so the price is a random walk independent of the real signs
           and zeta_K(a) is flat; every other run keeps its real sign
  W0       real flow, G0 response (one size-independent kernel (1 + l)^-beta, beta = gen3.g0_beta(cell),
           impact eps q^psi, psi = 0.5)
  W0-psi3  as W0 with psi = 0.3
  WB       real flow, size-dependent decay: beta(q) falls linearly in log q from 0.6 at the median log size
           of the stock's bottom quintile to -0.3 at the median of its top quintile, clipped to
           [-0.3, 0.6] outside (CHOICE: the two anchors are the quintiles' medians, per stock over its 8
           days, raw sizes); impact eps q^psi, psi = 0.5
  WM       the Paper II model, mbmodel.PAPER2, its own flow; day lengths from exploration stocks of each
           cell (panel3c.day_lengths); every stock-day is the first N trades after BURN trades of an
           mbmodel day (the burn-in trades' impact stays in the mid)
Diagnostic worlds (Phase 3d-2b, not in WORLDS): W-null-real is the pilot's W-null (informed runs keep their
real signs); W0-shuf is W0 on the real flow with every stock-day's signs permuted within the day (no sign
memory), burn-in included.
CORRECT gives, per world, the outcome that is correct for each primary test (None: no label imposed).

Thinning (WM only, descriptive): the observed flow of a venue that sees share s of the market. CHOICE: the
observed day keeps the cell's day length N (the venue's run count is what real data show), so the market
day behind it is about N / s trades; all settings of one stock-day thin the SAME market day (paired), and
share 1.0 is that market unthinned. Random rule: keep each trade with probability s. Router rule: route
with metaorder probability pi ~ Beta(kappa p, kappa (1 - p)), then hide h; the observed share is p (1 - h),
so p = s / (1 - h) (router_p). The observed trades carry the full-market mid (thinning.py).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

import explore as E
import gen3 as G
import mbmodel as M
import panel3c as PC
import phase3c as P3
import prereg3d as C
import primary3d as P
import surrogates as SU
import thinning as TH

BURN = 20_000  # 3c-4
N_DAYS = 8
PSI = 0.5
GF_F_INF = G.GF_INFORMED
WB_BETA = (0.6, -0.3)  # bottom-quintile and top-quintile anchors
WM_PARAMS = M.PAPER2
DAY_NS = (34_200 * 10**9, 57_600 * 10**9)  # 9:30 to 16:00 ET, for WM timestamps
WORLDS = ("W-null", "W0", "W0-psi3", "WB", "WM")
DIAG_WORLDS = ("W-null-real", "W0-shuf")
SEMI = ("W-null", "W0", "W0-psi3", "WB") + DIAG_WORLDS

# Correct outcome per (world, test): "reject" (Holm rejects and every condition holds), "no reject", "pass"
# (P-R: the equivalence established, which is its rejection with both conditions), "fail" (P-R not
# established), None (reported as observed).
CORRECT = {
    "W-null": {"P-V": "no reject", "P-R": None, "P-M": "no reject", "P-F": "no reject"},
    "W0": {"P-V": None, "P-R": "pass", "P-M": "no reject", "P-F": None},
    "W0-psi3": {"P-V": None, "P-R": "fail", "P-M": None, "P-F": None},
    "WB": {"P-V": None, "P-R": None, "P-M": "reject", "P-F": None},
    "WM": {"P-V": None, "P-R": None, "P-M": None, "P-F": None},
}
SIZE_WORLDS = ("W0", "W-null")  # P-M design rule: size constraint
POWER_WORLD = "WB"  # P-M design rule: objective
PV_SIZE_WORLD, PV_POWER_WORLD = "W-null", "WM"


def world_spec(world: str) -> dict:
    """flow ('real', 'real_s1' or 'mb'), signs ('real', 'coin_informed': GF's informed runs get fair coins,
    'day_shuffled': permuted within each stock-day, 'model'), response generator and psi."""
    spec = {
        "W-null": {"flow": "real_s1", "signs": "coin_informed", "gen": "GF", "psi": PSI},
        "W0": {"flow": "real", "signs": "real", "gen": "G0", "psi": PSI},
        "W0-psi3": {"flow": "real", "signs": "real", "gen": "G0", "psi": 0.3},
        "WB": {"flow": "real", "signs": "real", "gen": "WB", "psi": PSI},
        "WM": {"flow": "mb", "signs": "model", "gen": "MB", "psi": None},
        "W-null-real": {"flow": "real_s1", "signs": "real", "gen": "GF", "psi": PSI},
        "W0-shuf": {"flow": "real", "signs": "day_shuffled", "gen": "G0", "psi": PSI},
    }
    if world not in spec:
        raise ValueError(f"unknown world {world!r}")
    return spec[world]


# ---- real flow (exploration stocks, through the guard) ---------------------------------------------------


def flow_dir():
    return E.data_dir() / "phase3d2" / "flow"


def build_flow_cache(cell: str) -> dict:
    """The cell's exploration stocks with all 8 study days: sign, raw size and timestamp only, stock-major."""
    if cell not in C.CONFIRMATORY:
        raise ValueError(f"{cell} is not a confirmatory cell")
    days, _ = E.read_universe()
    per = {d: E.load_explore_runs(d, cell=cell) for d in days}
    syms = sorted(set.intersection(*(set(r) for r in per.values())))
    eps, q, ts, Nsd = [], [], [], []
    for s in syms:
        row = []
        for d in days:
            r = per[d][s]
            eps.append(r["sign"].astype(np.int8))
            q.append(r["q"].astype(np.float64))
            ts.append(r["ts"].astype(np.int64))
            row.append(len(r))
        Nsd.append(row)
    del per
    out = {"eps": np.concatenate(eps), "q": np.concatenate(q), "ts": np.concatenate(ts),
           "Nsd": np.array(Nsd, dtype=np.int64), "syms": np.array(syms)}
    path = flow_dir() / f"{cell}.npz"
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **out)
    return out


_CACHE: dict = {}


def flow_cache(cell: str) -> dict:
    if cell not in _CACHE:
        path = flow_dir() / f"{cell}.npz"
        if not path.exists():
            build_flow_cache(cell)
        z = np.load(path)
        _CACHE[cell] = {k: z[k] for k in z.files}
    return _CACHE[cell]


def resample_stocks(fc: dict, n: int, rng: np.random.Generator) -> dict:
    """n exploration stocks drawn with replacement, each with its real 8 days (eps, raw q, ts), stock-major."""
    Nsd = fc["Nsd"]
    per = Nsd.sum(axis=1)
    s0 = np.concatenate([[0], np.cumsum(per)[:-1]])
    pick = rng.integers(0, len(Nsd), n)
    idx = np.concatenate([np.arange(s0[s], s0[s] + per[s]) for s in pick])
    return {"eps": fc["eps"][idx].astype(float), "q": fc["q"][idx], "ts": fc["ts"][idx], "Nsd": Nsd[pick],
            "stocks": pick}


def shuffle_sizes(q: np.ndarray, lengths: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """S1 within every stock-day: q permuted over the day (surrogates.s1)."""
    out = np.empty_like(q)
    a = 0
    for n in lengths:
        out[a : a + n] = SU.s1(q[a : a + n], rng)
        a += n
    return out


def coin_informed(eps: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """GF's informed runs (a random GF_F_INF) with independent fair-coin signs; the rest keep theirs."""
    informed = rng.random(len(eps)) < GF_F_INF
    coin = np.where(rng.random(len(eps)) < 0.5, 1.0, -1.0)
    return np.where(informed, coin, eps), informed


# ---- responses on stationary days ------------------------------------------------------------------------


def burn_index(stock_start: int, Nd: np.ndarray, d: int, burn: int, rng: np.random.Generator) -> np.ndarray:
    """Indices (into the flow) of the burn-in of day d: a cyclic stretch of the stock's other days."""
    s0 = stock_start + np.concatenate([[0], np.cumsum(Nd)[:-1]])
    other = np.concatenate([np.arange(s0[j], s0[j] + Nd[j]) for j in range(len(Nd)) if j != d])
    start = rng.integers(0, len(other))
    return other[(start + np.arange(burn)) % len(other)]


def wb_anchors(q: np.ndarray) -> tuple[float, float]:
    """Median log size of the bottom and of the top quintile (raw sizes, ties inside the quintile)."""
    lq = np.log(q)
    lo, hi = np.quantile(lq, [0.2, 0.8])
    return float(np.median(lq[lq <= lo])), float(np.median(lq[lq >= hi]))


def wb_beta(q: np.ndarray, m1: float, m5: float, betas=WB_BETA) -> np.ndarray:
    """beta linear in log q: betas[0] at m1, betas[1] at m5, clipped between them."""
    b0, b1 = betas
    if not m5 > m1:
        return np.full(len(q), 0.5 * (b0 + b1))
    t = (np.log(q) - m1) / (m5 - m1)
    return np.clip(b0 + (b1 - b0) * t, min(b0, b1), max(b0, b1))


def wb_impact(x: np.ndarray, beta: np.ndarray, starts, lengths, n_nodes: int = 20) -> np.ndarray:
    """m_n = sum_{k<n} x_k (1 + n - k)^-beta_k per series: (1 + l)^-beta interpolated in beta on Chebyshev
    nodes (mbmodel._cheb_weights), one batched convolution per node, summed in frequency."""
    lo, hi = float(beta.min()), float(beta.max())
    if hi - lo < 1e-12:
        return PC.batched_causal_conv(x, starts, lengths, lo)
    nodes, L = M._cheb_weights(beta, lo, hi, n_nodes)
    return PC.batched_multi_conv([x * L[c] for c in range(len(nodes))], list(map(float, nodes)), starts, lengths)


def semi_response(world: str, cell: G.Cell, eps, q, Nsd, rng, burn: int = BURN, group: int | None = None,
                  informed: np.ndarray | None = None):
    """The synthetic mid of a semi-synthetic world on the given flow (raw q), independent stationary days.
    Stock-days are generated in groups of `group` (default panel3c.DAY_GROUP; 8 for WB, whose Chebyshev
    weights and 20 node series each hold a copy of the group's series: 64 peaked near 1 GB, measured).
    GF: `informed` is the informed-run mask (default: a random GF_F_INF of the runs, drawn here)."""
    sp = world_spec(world)
    group = group or (8 if sp["gen"] == "WB" else PC.DAY_GROUP)
    lengths = Nsd.ravel().astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    psi = sp["psi"]
    if sp["gen"] == "GF":
        if informed is None:
            informed = rng.random(len(eps)) < GF_F_INF
        return PC.response("GF", cell, eps, q, starts, lengths, rng, {"informed": informed, "psi": psi})
    S_, D_ = Nsd.shape
    per = Nsd.sum(axis=1)
    s_start = np.concatenate([[0], np.cumsum(per)[:-1]])
    theta = 1.0 / G.mean_q_psi(cell, psi, "lognormal")
    if sp["gen"] == "WB":
        anchors = [wb_anchors(q[s_start[s] : s_start[s] + per[s]]) for s in range(S_)]
    lm = np.empty(len(eps))
    sd_all = np.arange(S_ * D_)
    for g0 in range(0, len(sd_all), group):
        sds = sd_all[g0 : g0 + group]
        idx, ser, beta = [], [], []
        for i in sds:
            s, d = divmod(int(i), D_)
            b = burn_index(int(s_start[s]), Nsd[s], d, burn, rng)
            day = np.arange(starts[i], starts[i] + lengths[i])
            ii = np.concatenate([b, day])
            idx.append(ii)
            ser.append(len(ii))
            if sp["gen"] == "WB":
                beta.append(wb_beta(q[ii], *anchors[s]))
        idx = np.concatenate(idx)
        ser = np.array(ser, dtype=np.int64)
        sst = np.concatenate([[0], np.cumsum(ser)[:-1]])
        x = eps[idx] * q[idx] ** psi * theta
        if sp["gen"] == "G0":
            m = PC.batched_causal_conv(x, sst, ser, G.g0_beta(cell))
        else:
            m = wb_impact(x, np.concatenate(beta), sst, ser)
        m = m + PC.random_walk(sst, ser, 1.0, rng)
        for k, i in enumerate(sds):
            lm[starts[i] : starts[i] + lengths[i]] = m[sst[k] + burn : sst[k] + ser[k]]
    return lm


def semi_cell_flow(world: str, cell: str, rng: np.random.Generator, n_stocks: int | None = None) -> P.Flow:
    """One confirmatory cell of a semi-synthetic panel: resampled exploration stocks, synthetic mid."""
    if world not in SEMI:
        raise ValueError(f"{world} is not a semi-synthetic world")
    fc = flow_cache(cell)
    n = P3.n_confirm(cell) if n_stocks is None else n_stocks
    f = resample_stocks(fc, n, rng)
    lengths = f["Nsd"].ravel().astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int64)
    sp = world_spec(world)
    q = shuffle_sizes(f["q"], lengths, rng) if sp["flow"] == "real_s1" else f["q"]
    eps, informed = f["eps"], None
    if sp["signs"] == "coin_informed":
        eps, informed = coin_informed(eps, rng)
    elif sp["signs"] == "day_shuffled":
        eps = shuffle_sizes(eps, lengths, rng)  # S1's within-day permutation, applied to the signs
    lm = semi_response(world, G.load_cell(cell), eps, q, f["Nsd"], rng, informed=informed)
    return P.Flow(eps, q, lm, np.ones(len(lm), bool), f["ts"], starts, lengths, N_DAYS)


# ---- WM: the Paper II model ------------------------------------------------------------------------------


def mb_market_day(n_trades: int, rng: np.random.Generator, p: M.MBParams = WM_PARAMS, burn: int = BURN):
    """An mbmodel day with at least burn + n_trades trades, its first `burn` trades cut (their impact stays
    in the mid). Returns an MBDay whose trade arrays start after the burn-in."""
    n_meta = int(1.25 * (burn + n_trades) / G._trades_per_meta(p)) + 50
    while True:
        d = M.simulate_day(p, n_meta, rng)
        if len(d.t) >= burn + n_trades:
            break
        n_meta *= 2
    sl = slice(burn, None)
    return replace(d, t=d.t[sl], eps=d.eps[sl], q=d.q[sl], meta=d.meta[sl], rank=d.rank[sl], mid=d.mid[sl],
                   impact=d.impact[sl])


def clock_ts(t: np.ndarray) -> np.ndarray:
    """Model time of a day's trades spread linearly over 9:30 to 16:00 ET in ns (CHOICE: for S2's bins)."""
    t = np.asarray(t, dtype=float)
    span = t[-1] - t[0] if len(t) > 1 and t[-1] > t[0] else 1.0
    return (DAY_NS[0] + (t - t[0]) / span * (DAY_NS[1] - DAY_NS[0] - 1)).astype(np.int64)


def router_p(share: float, h: float) -> float:
    """Routing probability whose observed share p (1 - h) equals share."""
    p = share / (1.0 - h)
    if not 0 < p <= 1:
        raise ValueError(f"observed share {share} is not reachable with hidden share {h}")
    return p


def thin_day(day, setting: tuple, rng: np.random.Generator):
    """setting: ('full',), ('random', share) or ('router', share, h, kappa). Returns the observed MBDay."""
    if setting[0] == "full":
        return day
    if setting[0] == "random":
        return TH.random_thinning(day, setting[1], rng)[0]
    if setting[0] == "router":
        _, share, h, kappa = setting
        return TH.router_thinning(day, router_p(share, h), rng, h, kappa)[0]
    raise ValueError(f"unknown thinning rule {setting[0]!r}")


def min_share(setting: tuple) -> float:
    return 1.0 if setting[0] == "full" else float(setting[1])


def wm_cell_flows(cell: str, rng: np.random.Generator, settings=(("full",),), n_stocks: int | None = None,
                  thin_rngs: list | None = None, keep=None) -> list:
    """WM flows of one cell, one per kept thinning setting, all from the same market days (paired). rng draws
    the day lengths and the markets; thin_rngs (one per setting, default rng) the thinning. Every setting
    in `settings` is thinned and must reach N observed trades (else the market day is redrawn longer), but
    only those in `keep` (indices, default all) are stored: a job that keeps a subset still makes the same
    redraw decisions, so it sees the same markets. Each observed day is the first N observed trades after
    the burn-in, N the cell's day length. Returns compact dicts (to_flow makes the Flow; q as float32)."""
    n = P3.n_confirm(cell) if n_stocks is None else n_stocks
    thin_rngs = [rng] * len(settings) if thin_rngs is None else thin_rngs
    keep = list(range(len(settings))) if keep is None else list(keep)
    Nsd = PC.day_lengths(cell, n, N_DAYS, rng)
    lengths = Nsd.ravel().astype(np.int64)
    smin = min(min_share(s) for s in settings)
    parts = [{k: [] for k in ("eps", "q", "lm", "ts")} for _ in keep]
    for N in lengths:
        need = int(np.ceil(1.3 * N / smin)) + 50
        while True:
            day = mb_market_day(need, rng)
            obs = [thin_day(day, s, r) for s, r in zip(settings, thin_rngs, strict=True)]
            if all(len(o.t) >= N for o in obs):
                break
            need *= 2
        for k, o in enumerate(obs[j] for j in keep):
            parts[k]["eps"].append(o.eps[:N].astype(np.int8))
            parts[k]["q"].append(o.q[:N].astype(np.float32))
            parts[k]["lm"].append(o.mid[:N])
            parts[k]["ts"].append(clock_ts(o.t[:N]))
        del day, obs
    return [{**{k: np.concatenate(v) for k, v in pt.items()}, "lengths": lengths} for pt in parts]


def to_flow(d: dict) -> P.Flow:
    lengths = d["lengths"]
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int64)
    return P.Flow(d["eps"].astype(float), d["q"].astype(float), d["lm"], np.ones(len(d["lm"]), bool), d["ts"],
                  starts, lengths, N_DAYS)


# ---- panels -----------------------------------------------------------------------------------------------


def cells_of(tb: str) -> list[str]:
    return sorted(c for c in C.CONFIRMATORY if C.tick_bin(c) == tb)


def tick_bin_flows(world: str, tb: str, rng: np.random.Generator) -> dict:
    """{cell: Flow} of the tick bin's confirmatory cells in one panel of the world."""
    if world == "WM":
        return {c: to_flow(wm_cell_flows(c, rng)[0]) for c in cells_of(tb)}
    return {c: semi_cell_flow(world, c, rng) for c in cells_of(tb)}


# ---- P-M design rule (stated before the study's results were seen) ----------------------------------------

PM_DESIGNS = tuple(P.PMParams(nb, lm, C.PM_L_FACTOR * lm, mr) for nb in C.PM_BINS_ALLOWED
                   for lm in C.PM_LMAX_ALLOWED for mr in (20, 50, 100))
SIZE_LEVEL = 0.05


def design_key(pm: P.PMParams) -> str:
    return f"b{pm.n_bins}_l{pm.l_max}_L{pm.L}_m{pm.min_runs_per_bin}"


def choose_design(rates: dict, n_panels: dict) -> dict:
    """The P-M design rule, fixed before the study (CHOICE where the instruction left room):
    rates[(world, design key, tick bin)] = share of that world's panels in which P-M rejects at the
    UNADJUSTED level with both conditions (p <= 0.05, slope and top-bin upper bounds below 0); a pair that
    is not evaluable counts as not rejected. n_panels[world] = panels.
      admissible: in every tick bin, for W0 and for W-null, rate <= 0.05 + 2 MC SE, MC SE = sqrt(0.05 x 0.95
        / panels) (the SE at the nominal level, so the bound does not loosen with the observed rate);
      objective: the WB rate averaged over the three tick bins;
      ties (equal objective): fewer bins, then smaller l_max, then smaller min runs per bin.
    If no design is admissible, none is chosen and the design with the smallest worst null rate is reported
    as the fallback (same ties). Returns {'chosen', 'admissible', 'objective', 'worst_null', 'fallback'}."""
    keys = [design_key(d) for d in PM_DESIGNS]
    order = {design_key(d): (d.n_bins, d.l_max, d.min_runs_per_bin) for d in PM_DESIGNS}
    obj, worst, adm = {}, {}, {}
    for k in keys:
        obj[k] = float(np.mean([rates[(POWER_WORLD, k, tb)] for tb in C.TICK_BINS]))
        ok, w = True, 0.0
        for wd in SIZE_WORLDS:
            bound = SIZE_LEVEL + 2 * np.sqrt(SIZE_LEVEL * (1 - SIZE_LEVEL) / n_panels[wd])
            for tb in C.TICK_BINS:
                r = rates[(wd, k, tb)]
                w = max(w, r)
                ok &= r <= bound
        adm[k], worst[k] = bool(ok), w
    cand = [k for k in keys if adm[k]]
    chosen = min(cand, key=lambda k: (-obj[k], order[k])) if cand else None
    fallback = None if cand else min(keys, key=lambda k: (worst[k], order[k]))
    return {"chosen": chosen, "admissible": adm, "objective": obj, "worst_null": worst, "fallback": fallback}


def design_by_key(key: str) -> P.PMParams:
    return next(d for d in PM_DESIGNS if design_key(d) == key)
