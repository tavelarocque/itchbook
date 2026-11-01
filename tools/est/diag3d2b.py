"""Phase 3d-2b: four short diagnostics before the 3d-2 study (docs/phase3d2b_report.md).

  occupancy  tie-aware size bins (primary3d.tie_bins) on the exploration stocks' sizes only: realised bins
             per stock-day (P-M) and per stock (psi-hat), share of stock-days usable, stocks with a psi-hat
             slope defined by the bins, against the old quantile bins (ties to the lower bin)
  pnull      the clean P-F null: W-null (informed runs with fair-coin signs) against the pilot's W-null-real,
             P-F alone, N_PNULL panels per tick bin
  wm         the Paper II model at the cells' day lengths (as in the pilot), N_WM panels per tick bin
  long       the Paper II model on long days: LONG["n_stocks"] stocks x 8 days x N = 100000 runs, fitting
             ranges [64, 512] and [128, 1024], the same estimators (primary3d with its T grid and ranges
             replaced), N_LONG panels
  abias      a*(T) at every T in W0 against W0-shuf (signs permuted within the day), N_ABIAS panels
  vertex     sizes only: the expected T = 1 curve R(a) of W0 without sign memory, its argmax against the
             quadratic vertex that a* is

Every job writes one JSON under OUT/<kind>/ and is skipped when its file exists (resumable). At most 3
workers. Real data enter only through worlds3d2's flow cache (signs, sizes, timestamps; no prices).

    uv run --group analysis tools/est/diag3d2b.py occupancy|vertex|pnull|wm|long|abias|all|summarise [--workers 3]
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import explore as E
import mbmodel as M
import pooled as PO
import prereg3d as C
import primary3d as P
import worlds3d2 as W

OUT = Path(os.environ.get("PHASE3D2B_OUT", E.data_dir() / "phase3d2" / "diag3d2b"))
SEED = 20261101
N_PNULL, N_WM, N_LONG, N_ABIAS = 20, 30, 30, 30
PM_DIAG = (P.PMParams(3, 16, 128, 50), P.PMParams(3, 32, 256, 50))  # CHOICE: the pilot's design and l_max 32
OCC_BINS, OCC_MIN_RUNS = C.PM_BINS_ALLOWED, (20, 50, 100)
LONG = {"n_stocks": 50, "N": 100_000, "rho": 0.2,  # CHOICE: ceiling 1/5 (any rho >= 0.0103 reaches 1024)
        "T_grid": 2 ** np.arange(11), "ranges": ((64, 512), (128, 1024))}
KINDS = ("pnull", "wm", "long", "abias")


def _f(x):
    x = np.asarray(x, dtype=float)
    if x.ndim == 0:
        return float(x) if np.isfinite(x) else None
    return np.where(np.isfinite(x), x, np.nan).tolist()


def _summ(s: dict) -> dict:
    return {k: _f(s[k]) for k in ("estimate", "lo", "hi")}


# ---- 1. occupancy on exploration sizes ----------------------------------------------------------------------


def old_bins_full(q: np.ndarray, nb: int, min_runs: int) -> bool:
    """The pilot's rule: quantile edges, ties to the lower bin; usable when every bin holds min_runs."""
    b = np.searchsorted(np.quantile(q, np.linspace(0, 1, nb + 1)[1:-1]), q, side="left")
    return bool(np.bincount(b, minlength=nb).min() >= min_runs)


def occupancy() -> dict:
    """Sizes only (R2-clipped as make_unit does), every exploration stock of each confirmatory cell once."""
    out = {}
    for tb in C.TICK_BINS:
        rec = {"pm": {}, "psi": {}}
        qs_sd, qs_st = [], []
        for c in W.cells_of(tb):
            fc = W.flow_cache(c)
            lengths = fc["Nsd"].ravel()
            q = P.r2_clip(fc["q"], lengths, W.N_DAYS)
            days = np.split(q, np.cumsum(lengths)[:-1])
            qs_sd += days
            # psi-hat bins live on the steps n with n + 1 in the day: every run but the day's last
            qs_st += [np.concatenate([d[:-1] for d in days[s * W.N_DAYS : (s + 1) * W.N_DAYS]])
                      for s in range(len(fc["Nsd"]))]
        rec["n_stocks"], rec["n_stock_days"] = len(qs_st), len(qs_sd)
        rec["round_lot_share"] = float(np.mean(np.concatenate(qs_sd) % 100 == 0))
        for nb in OCC_BINS:
            for mr in OCC_MIN_RUNS:
                k = np.array([len(P.tie_bins(d, nb, mr)[1]) + 1 if len(d) else 0 for d in qs_sd])
                old = np.array([old_bins_full(d, nb, mr) if len(d) >= nb else False for d in qs_sd])
                rec["pm"][f"b{nb}_m{mr}"] = {
                    "usable_share": float(np.mean(k == nb)), "old_usable_share": float(np.mean(old)),
                    "realised": {int(v): int(np.sum(k == v)) for v in np.unique(k)}}
        kp = np.array([len(P.tie_bins(d, P.PSI_BINS, P.PSI_MIN_RUNS)[1]) + 1 if len(d) else 0 for d in qs_st])
        old_full = np.array([len(d) >= P.PSI_BINS and np.bincount(
            np.searchsorted(np.quantile(d, np.linspace(0, 1, P.PSI_BINS + 1)[1:-1]), d, side="left"),
            minlength=P.PSI_BINS).min() > 0 for d in qs_st])
        rec["psi"] = {"realised": {int(v): int(np.sum(kp == v)) for v in np.unique(kp)},
                      "slope_defined_by_bins": int(np.sum(kp >= P.PSI_MIN_BINS)),
                      "old_all_deciles_nonempty": int(np.sum(old_full)), "n_stocks": len(qs_st)}
        out[tb] = rec
    path = OUT / "occupancy.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    return out


def vertex(psi: float = 0.5, k: float = 0.3) -> dict:
    """Sizes only: the T = 1 correlation R(a) that W0 without sign memory has in expectation given the sizes,
    per stock-day c(a) = mean(qc^a q^psi) / sqrt(mean(qc^2a) mean(q^2psi)) (response on raw q, regressor
    R2-clipped, as in W0), pooled as core2 pools z = arctanh(k c) (k: signal scale, the vertex barely moves
    with it). Its argmax is psi on every day (Cauchy-Schwarz); the record holds the quadratic vertex over
    the estimator's a range [0, 1.5], over [0, 1], the argmax on a fine grid and the shape R(a) / R(psi)."""
    g = PO.A_GRID[PO.A_GRID <= PO.A_PF]
    fine = np.linspace(0.0, float(g[-1]), 151)

    def vx(a, r):
        c = np.linalg.pinv(np.stack([np.ones_like(a), a, a * a], axis=1)) @ r
        return float(-c[1] / (2 * c[2]))

    out = {}
    for tb in C.TICK_BINS:
        zs, zf = [], []
        for c in W.cells_of(tb):
            fc = W.flow_cache(c)
            lengths = fc["Nsd"].ravel()
            qc = np.split(P.r2_clip(fc["q"], lengths, W.N_DAYS), np.cumsum(lengths)[:-1])
            qr = np.split(fc["q"], np.cumsum(lengths)[:-1])
            for s in range(len(fc["Nsd"])):
                sl = slice(s * W.N_DAYS, (s + 1) * W.N_DAYS)
                rows = [(r, x) for r, x in zip(qr[sl], qc[sl], strict=True) if len(r) >= 2]
                if not rows:
                    continue
                for grid, acc in ((g, zs), (fine, zf)):
                    cv = [[np.mean(x**a * r**psi) / np.sqrt(np.mean(x ** (2 * a)) * np.mean(r ** (2 * psi)))
                           for a in grid] for r, x in rows]
                    acc.append(np.arctanh(k * np.array(cv)).mean(axis=0))
        R, Rf = np.tanh(np.mean(zs, axis=0)), np.tanh(np.mean(zf, axis=0))
        out[tb] = {"n_stocks": len(zs), "a_grid": g.tolist(), "shape": (R / R[g == psi][0]).tolist(),
                   "vertex_0_1.5": vx(g, R), "vertex_0_1": vx(g[g <= 1], R[g <= 1]),
                   "argmax_fine": float(fine[np.argmax(Rf)])}
    path = OUT / "vertex.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    return out


# ---- shared: one unit's records ------------------------------------------------------------------------------


def a_star_all_T(u: P.Unit) -> list:
    """a*(T) at every T of the grid, all days and stocks (pooled.core2 as run_unit calls it)."""
    S_, D_ = u.flow.Nsd.shape
    return _f(PO.core2(u.raw, np.ones((1, D_)), np.ones((1, S_)), u.specs, "tilde", **P.KW)["a_star"][0, :, 0])


def logq_sd(fl: P.Flow, used: np.ndarray, d: P.PMParams) -> np.ndarray:
    """Mean log size per tie-aware bin of every used stock-day (n_sd, n_bins), NaN elsewhere."""
    x = np.full((len(fl.lengths), d.n_bins), np.nan)
    for i in np.flatnonzero(used):
        q = fl.q[fl.day(i)]
        b, _ = P.tie_bins(q, d.n_bins, d.min_runs_per_bin)
        x[i] = np.bincount(b, np.log(q), minlength=d.n_bins) / np.bincount(b, minlength=d.n_bins)
    return x


def pooled_logq(x: np.ndarray, S_: int, D_: int) -> list:
    """Absolute mean log size per bin, equal weight per stock of its days' values (as _pooled)."""
    return _f(P._pooled(x, np.ones((1, D_)), np.ones((1, S_)), D_)[0])


def unit_record(res: dict, u: P.Unit | None = None, pm_logq: dict | None = None) -> dict:
    """p, conditions and the statistics of every computed test, both ranges."""
    rec = {}
    for name, t in res["tests"].items():
        r = {"p": float(t.p), "conditions": {k: bool(v) for k, v in t.conditions.items()},
             "range_used": list(t.range_used) if t.range_used else None}
        if name == "P-V":
            r["by_range"] = {rn: {"p": v["p"], "conditions": v["conditions"],
                                  "errors": {k: _f(v["loo_errors"][k]) for k in ("M1", "M2", "M3")},
                                  "days_used": v["loo_errors"]["days_used"],
                                  "m2_all_days": {k: _f(v["loo_errors"]["M2_fit_all_days"][k])
                                                  for k in ("mu_m", "lambda_sigma2", "a_c", "k_below")}
                                  if "M2_fit_all_days" in v["loo_errors"] else None,
                                  "slope_S1": _summ(v["slope_S1"]), "slope_S2": _summ(v["slope_S2"]),
                                  "zeta_V": _summ(v["zeta_V"])} for rn, v in t.by_range.items()}
        elif name == "P-R":
            r["by_range"] = {rn: {"p": v["p"], "conditions": v["conditions"], "a_star_T1_T2": _summ(v["a_star_T1_T2"])}
                             for rn, v in t.by_range.items()}
            r["psi_hat"] = _summ(t.secondary["psi_hat"])
            r["psi_n_stocks"] = t.secondary["psi_n_stocks"]
        elif name == "P-F":
            r["by_range"] = {rn: {"p": v["p"], "conditions": v["conditions"], "pf": _summ(v["pf"]),
                                  "zeta_K": _summ(v["zeta_K_levels"])} for rn, v in t.by_range.items()}
        elif name == "P-M":
            st = t.statistic
            r.update(stocks=len(t.stocks_used.get("stocks", [])), stock_days=t.stocks_used.get("stock_days", 0))
            if st:
                r.update(beta=_summ(st["beta_by_bin"]), xlog=_summ(st["xlog_by_bin"]), slope=_summ(st["slope"]),
                         top=_summ(st["top"]), p_slope=st["p_slope"], p_top=st["p_top"])
        rec[name] = r
    if u is not None:
        rec["a_star_all_T"] = a_star_all_T(u)
        rec["T_grid"] = C.T_GRID.tolist()
    if pm_logq is not None:
        rec["pm_logq"] = pm_logq
    return rec


def pm_records(u: P.Unit, saved: np.random.Generator, by, logq) -> dict:
    """P-M at every PM_DIAG design on run_unit's draws (P.draws from a copy of the generator as it was before
    run_unit); logq(design) gives the per-stock-day bin log sizes."""
    S_, D_ = u.flow.Nsd.shape
    out = {}
    for d in PM_DIAG:
        ev = P.evaluability(u.Ns, u.rho, d)["P-M"]
        if not ev["evaluable"]:
            out[W.design_key(d)] = {"evaluable": False, "reason": ev["reason"]}
            continue
        tests: dict = {}
        P._pm_test(u, d, P.draws(S_, D_, copy.deepcopy(saved), C.N_BOOT), P.CHUNK, by, tests)
        rec = unit_record({"tests": tests}, pm_logq=pooled_logq(logq(d), S_, D_))
        out[W.design_key(d)] = {"evaluable": True, **rec["P-M"], "logq": rec["pm_logq"]}
    return out


def run_semi_unit(flows: dict, tb: str, rng: np.random.Generator, tests: tuple, pm: bool) -> dict:
    """run_unit on the tick bin's unit for the evaluable range tests, plus every PM_DIAG design if pm."""
    u = P.tick_bin_unit(flows, tb)
    ev = P.evaluability(u.Ns, u.rho, PM_DIAG[0])
    want = tuple(t for t in tests if ev[t]["evaluable"])
    saved = copy.deepcopy(rng)
    res = P.run_unit(u, PM_DIAG[0], rng, tests=want)
    rec = unit_record(res, u)
    rec["evaluable"] = {t: ev[t]["evaluable"] for t in tests}
    rec["n_stocks"] = int(u.flow.n_stocks)
    if pm:
        _, by = P.select_range(u.Ns, u.rho)
        rec["pm"] = pm_records(u, saved, by, lambda d: logq_sd(u.flow, P.pm_moments(u.flow, u.rho, d)["used"], d))
    return rec


# ---- 3b. long days: the unit built one stock at a time ---------------------------------------------------------


@contextmanager
def grid(T_grid, ranges):
    """primary3d with another T grid and fitting ranges (the later range primary), restored after."""
    saved = (C.T_GRID, C.RANGES, C.PRIMARY_RANGE, C.FALLBACK_RANGE, P.KW, P.RANGE_NAMES)
    C.T_GRID, C.RANGES = np.asarray(T_grid), tuple(tuple(r) for r in ranges)
    C.FALLBACK_RANGE, C.PRIMARY_RANGE = C.RANGES[0], C.RANGES[-1]
    P.KW = dict(T_grid=C.T_GRID, a_grid=C.A_GRID, ranges=C.RANGES, a_max=C.A_MAX)
    P.RANGE_NAMES = tuple(f"[{lo},{hi}]" for lo, hi in C.RANGES)
    try:
        yield
    finally:
        C.T_GRID, C.RANGES, C.PRIMARY_RANGE, C.FALLBACK_RANGE, P.KW, P.RANGE_NAMES = saved


PARTS = ("m3_V", "surrogate_raws", "sig2_per_sd", "psi_moments", "pm_moments")


def stock_parts(fl: P.Flow, rho: float, rng: np.random.Generator, designs, n_sur: int = P.N_SUR) -> dict:
    """Every per-stock-day input of run_unit for one stock's (clipped) flow, in run_unit's functions."""
    st = SimpleNamespace(flow=fl)
    return {"raw": P.moments_raw(fl), "m3_V": P.m3_V(fl),
            "sur": {k: P.surrogate_raws(st, k, n_sur, rng) for k in ("S1", "S2")},
            "sig2_per_sd": P.sig2_per_sd(fl), "psi_moments": P.psi_moments(fl),
            "pm": {W.design_key(d): P.pm_moments(fl, np.array([rho]), d) for d in designs}}


def _cat_raw(raws: list) -> PO.Raw:
    return PO.Raw(*(np.concatenate([getattr(r, k) for r in raws]) for k in ("V", "K", "Sm", "K1_t", "Nsd")))


def assemble(parts: list, n_days: int, rho: float) -> tuple[P.Unit, dict]:
    """The unit of all stocks (stub flow: day lengths only) and the concatenated per-stock-day inputs."""
    lengths = np.concatenate([p["raw"].Nsd.ravel() for p in parts]).astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int64)
    z = np.zeros(0)
    stub = P.Flow(z, z, z, z.astype(bool), z.astype(np.int64), starts, lengths, n_days)
    S_ = len(parts)
    u = P.Unit("long", stub, np.full(S_, rho), np.full(S_, "long", dtype=object), (),
               _cat_raw([p["raw"] for p in parts]))
    psi = {k: np.concatenate([p["psi_moments"][k] for p in parts]) for k in ("ysum", "cnt", "lmed", "n_real")}
    pm = {key: {k: np.concatenate([p["pm"][key][k] for p in parts]) for k in ("Cp", "Sr", "xlog", "used")}
          for key in parts[0]["pm"]}
    n_sur = len(parts[0]["sur"]["S1"])
    inputs = {"m3_V": np.concatenate([p["m3_V"] for p in parts]),
              "sur": {k: [_cat_raw([p["sur"][k][j] for p in parts]) for j in range(n_sur)] for k in ("S1", "S2")},
              "sig2_per_sd": np.concatenate([p["sig2_per_sd"] for p in parts]), "psi_moments": psi, "pm": pm}
    return u, inputs


@contextmanager
def patched_inputs(inputs: dict):
    """run_unit's per-stock-day functions return the precomputed inputs (run_unit itself unchanged)."""
    saved = {k: getattr(P, k) for k in PARTS}
    P.m3_V = lambda flow: inputs["m3_V"]
    P.surrogate_raws = lambda u, kind, n, rng: inputs["sur"][kind][:n]
    P.sig2_per_sd = lambda flow: inputs["sig2_per_sd"]
    P.psi_moments = lambda flow, *a, **k: inputs["psi_moments"]
    P.pm_moments = lambda flow, rho, pm: inputs["pm"][W.design_key(pm)]
    try:
        yield
    finally:
        for k, f in saved.items():
            setattr(P, k, f)


def long_stock(rng: np.random.Generator, n_days: int, N: int) -> P.Flow:
    """One stock: n_days Paper II days of N trades after the burn-in, sizes R2-clipped."""
    eps, q, lm, ts = [], [], [], []
    for _ in range(n_days):
        d = W.mb_market_day(N, rng)
        eps.append(d.eps[:N].astype(float))
        q.append(d.q[:N].astype(float))
        lm.append(d.mid[:N])
        ts.append(W.clock_ts(d.t[:N]))
    lengths = np.full(n_days, N, dtype=np.int64)
    starts = np.arange(n_days, dtype=np.int64) * N
    qq = np.concatenate(q)
    return P.Flow(np.concatenate(eps), P.r2_clip(qq, lengths, n_days), np.concatenate(lm), np.ones(len(qq), bool),
                  np.concatenate(ts), starts, lengths, n_days)


def long_panel(rng: np.random.Generator, n_stocks: int = LONG["n_stocks"], N: int = LONG["N"],
               rho: float = LONG["rho"]) -> dict:
    """One long-day panel: stocks generated and reduced one at a time (RAM), then run_unit on the whole."""
    with grid(LONG["T_grid"], LONG["ranges"]):
        parts = []
        for _ in range(n_stocks):
            fl = long_stock(rng, W.N_DAYS, N)
            p = stock_parts(fl, rho, rng, PM_DIAG)
            p["logq"] = {W.design_key(d): logq_sd(fl, p["pm"][W.design_key(d)]["used"], d) for d in PM_DIAG}
            parts.append(p)
            del fl
        u, inputs = assemble(parts, W.N_DAYS, rho)
        logq = {k: np.concatenate([p["logq"][k] for p in parts]) for k in parts[0]["logq"]}
        del parts
        saved = copy.deepcopy(rng)
        with patched_inputs(inputs):
            rec = unit_record(P.run_unit(u, PM_DIAG[0], rng, tests=("P-V", "P-R", "P-F")), u)
            _, by = P.select_range(u.Ns, u.rho)
            rec["pm"] = pm_records(u, saved, by, lambda d: logq[W.design_key(d)])
        rec["n_stocks"] = n_stocks
    return rec


# ---- the paper's predictions for PAPER2 --------------------------------------------------------------------------


def paper_zeta_V(a, p: M.MBParams = W.WM_PARAMS, T: float | None = None) -> np.ndarray:
    """Paper I Eq. 30 (T = None): zeta_V(a) = 3 - mu_m - 2 a lambda sigma_l^2 below a_c = (2 - mu_m) /
    (2 lambda sigma_l^2) (Eq. 31), 1 above. With T, Eq. 29's effective exponent mu_m + lambda sigma_l^2
    (2 a - lambda log T / 2) instead."""
    a = np.asarray(a, dtype=float)
    ls = p.lam * p.sigma_l**2
    mu = p.mu_m + ls * 2 * a if T is None else p.mu_m + ls * (2 * a - p.lam * np.log(T) / 2)
    return np.where(3 - mu > 1, 3 - mu, 1.0)


def paper_a_star_band(T, p: M.MBParams = W.WM_PARAMS) -> tuple[np.ndarray, np.ndarray]:
    """Paper I Eq. 74: a* is a weighted average of 1/2 (the A, diagonal term) and lambda log T (the B term),
    so it lies between them; Eqs. 71 and 72 give B / A growing like T^((mu_m/2 - gamma_x - beta_m) -
    (1 - mu_m - lambda sigma_l^2)/2) (0.21 for PAPER2, reading the B term's size exponent at beta_m), so a*
    moves toward lambda log T as T grows."""
    c = p.lam * np.log(np.asarray(T, dtype=float))
    return np.minimum(0.5, c), np.maximum(0.5, c)


def paper_ba_exponent(p: M.MBParams = W.WM_PARAMS) -> float:
    return (p.mu_m / 2 - p.gamma_x - p.beta_m) - (1 - p.mu_m - p.lam * p.sigma_l**2) / 2


def paper_beta(logq, p: M.MBParams = W.WM_PARAMS) -> np.ndarray:
    """The generator's decay exponent beta_q = beta_m - lambda' (log q - m), clipped into beta_clip."""
    return np.clip(p.beta_m - p.lam_p * (np.asarray(logq, dtype=float) - p.m), *p.beta_clip)


# ---- jobs ------------------------------------------------------------------------------------------------------


def job_path(kind: str, world: str, tb: str, i: int) -> Path:
    return OUT / kind / f"{world}_{tb}_{i:03d}.json"


def job_rng(kind: str, world: str, tb: str, i: int) -> np.random.Generator:
    return np.random.default_rng([SEED, KINDS.index(kind), (W.WORLDS + W.DIAG_WORLDS).index(world),
                                  (C.TICK_BINS + ("long",)).index(tb), i])


def job(args) -> str:
    kind, world, tb, i = args
    path = job_path(kind, world, tb, i)
    if path.exists():
        return str(path)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    rng = job_rng(kind, world, tb, i)
    t0 = time.time()
    if kind == "long":
        rec = long_panel(rng)
    else:
        flows = W.tick_bin_flows(world, tb, rng)
        tg = time.time() - t0
        tests = {"pnull": ("P-F",), "wm": C.PRIMARY_TESTS[:2] + ("P-F",), "abias": ("P-R",)}[kind]
        rec = run_semi_unit(flows, tb, rng, tests, pm=kind == "wm")
        rec["generate_s"] = tg
    rec.update(kind=kind, world=world, tick_bin=tb, panel=i, seconds=time.time() - t0)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec))
    tmp.rename(path)
    return str(path)


def jobs(kind: str) -> list:
    if kind == "pnull":
        return [(kind, w, tb, i) for w in ("W-null", "W-null-real") for tb in C.TICK_BINS for i in range(N_PNULL)]
    if kind == "wm":
        return [(kind, "WM", tb, i) for tb in C.TICK_BINS for i in range(N_WM)]
    if kind == "long":
        return [(kind, "WM", "long", i) for i in range(N_LONG)]
    if kind == "abias":
        return [(kind, w, tb, i) for w in ("W0", "W0-shuf") for tb in C.TICK_BINS for i in range(N_ABIAS)]
    raise ValueError(kind)


def run(kinds, workers: int) -> None:
    js = [j for k in kinds for j in jobs(k)]
    js.sort(key=lambda j: (j[0] != "long", j[2] != "Q-T1"))  # the longest jobs first
    todo = [j for j in js if not job_path(*j).exists()]
    print(f"{len(js)} jobs, {len(todo)} to run on {workers} workers", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for k, _ in enumerate(ex.map(job, todo), 1):
            if k % 10 == 0 or k == len(todo):
                print(f"  {k}/{len(todo)}, {(time.time() - t0) / 60:.1f} min", flush=True)


# ---- summary ---------------------------------------------------------------------------------------------------


def load(kind: str) -> list:
    return [json.loads(f.read_text()) for f in sorted((OUT / kind).glob("*.json"))]


def _ms(x) -> tuple[float, float, int]:
    """Mean and MC SE over panels of the finite values, and their count."""
    x = np.asarray([v for v in x if v is not None and np.isfinite(v)], dtype=float)
    if len(x) == 0:
        return np.nan, np.nan, 0
    return float(x.mean()), float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else np.nan, len(x)


def _rate(flags) -> str:
    f = np.asarray(list(flags), dtype=float)
    if len(f) == 0:
        return "n/a"
    se = np.sqrt(max(f.mean() * (1 - f.mean()), 1e-12) / len(f))
    return f"{int(f.sum())}/{len(f)} = {f.mean():.2f} ({se:.2f})"


def _fmt(m, se=None, d=3) -> str:
    if m is None or not np.isfinite(m):
        return "n/a"
    return f"{m:.{d}f}" if se is None or not np.isfinite(se) else f"{m:.{d}f} ({se:.{d}f})"


def _cols(recs, get) -> list:
    """Per-element mean (MC SE) of a list-valued statistic across panels."""
    vals = [get(r) for r in recs]
    vals = [v for v in vals if v is not None]
    if not vals:
        return []
    n = max(len(v) for v in vals)
    return [_fmt(*_ms([v[i] if i < len(v) else None for v in vals])[:2]) for i in range(n)]


def passed(t: dict) -> bool:
    """Unadjusted decision: p <= 0.05 and every condition (the Holm decision at its most lenient step)."""
    return bool(t["p"] <= C.ALPHA_FW and t["conditions"] and all(t["conditions"].values()))


def _iv(xs: list, get) -> str:
    """Mean estimate [mean lower, mean upper bound] over panels of get(x) -> {'estimate', 'lo', 'hi'}."""
    m = {k: _fmt(_ms([get(x)[k] for x in xs])[0]) for k in ("estimate", "lo", "hi")}
    return f"{m['estimate']} [{m['lo']}, {m['hi']}]"


def _row(*cells) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def _head(*cols) -> list:
    return [_row(*cols), "|" + "---|" * len(cols)]


def _occupancy_tables(occ: dict) -> list:
    L = ["## 1. Tie-aware size bins on exploration sizes", "",
         *_head("tick bin", "stocks", "stock-days", "round-lot share", "psi-hat: stocks with >= 3 bins (new)",
                "psi-hat: all 10 deciles non-empty (old)", "realised psi bins (bins: stocks)")]
    for tb, r in occ.items():
        L.append(_row(tb, r["n_stocks"], r["n_stock_days"], f"{r['round_lot_share']:.2f}",
                      r["psi"]["slope_defined_by_bins"], r["psi"]["old_all_deciles_nonempty"],
                      ", ".join(f"{k}: {v}" for k, v in r["psi"]["realised"].items())))
    L += ["", *_head("tick bin", "bins, min runs", "usable stock-days, new", "usable, old",
                     "realised bins (bins: stock-days)")]
    for tb, r in occ.items():
        for k, v in r["pm"].items():
            L.append(_row(tb, k, f"{v['usable_share']:.2f}", f"{v['old_usable_share']:.2f}",
                          ", ".join(f"{a}: {b}" for a, b in v["realised"].items())))
    return L + [""]


def _rng_name(t: dict) -> str:
    return f"[{t['range_used'][0]},{t['range_used'][1]}]"


def _pnull_tables(pn: list) -> list:
    L = ["## 2. The clean P-F null", "",
         *_head("world", "tick bin", "panels", "P-F passes (unadjusted)", "p <= 0.05", "slope of zeta_K on a",
                "zeta_K(a), a = 0, 0.25, 0.5, 0.75, 1")]
    for w in ("W-null", "W-null-real"):
        for tb in C.TICK_BINS:
            fs = [r["P-F"] for r in pn if r["world"] == w and r["tick_bin"] == tb and "P-F" in r]
            if not fs:
                continue
            br = [f["by_range"][_rng_name(f)] for f in fs]
            L.append(_row(w, tb, len(fs), _rate(passed(f) for f in fs), _rate(f["p"] <= 0.05 for f in fs),
                          _fmt(*_ms([x["pf"]["estimate"] for x in br])[:2]),
                          ", ".join(_cols(br, lambda x: x["zeta_K"]["estimate"]))))
    return L + [""]


def _wm_tables(title: str, rs_all: list, units) -> list:
    L = [f"## {title}", ""]
    groups = [(tb, [r for r in rs_all if r["tick_bin"] == tb]) for tb in units]
    T_grid = rs_all[0]["T_grid"]
    lo, hi = paper_a_star_band(T_grid)
    L += ["a*(T), all days and stocks pooled (the paper: Eq. 74 puts a* between 1/2 and lambda log T):", "",
          *_head("unit", "panels", *(f"T = {t}" for t in T_grid)),
          _row("Eq. 74 band", "", *(f"{a:.2f} to {b:.2f}" for a, b in zip(lo, hi, strict=True)))]
    for tb, rs in groups:
        L.append(_row(tb, len(rs), *_cols(rs, lambda r: r["a_star_all_T"])))
    z30 = ", ".join(f"{v:.2f}" for v in paper_zeta_V(C.A_GRID))
    L += ["", "P-V, P-R, P-F by range:", "",
          *_head("unit", "range", "M1 error", "M2 error", "M3 error", "M2 < M1", "M2 < M3",
                 "real - S1 slope [lo, hi]", "real - S2 slope [lo, hi]", "S1 upper < 0", "S2 upper < 0",
                 "P-V passes", "zeta_V(a)", "Eq. 30", "a*(T1), a*(T2)", "P-R passes", "P-F passes")]
    for tb, rs in groups:
        rs_v = [r for r in rs if "P-V" in r]
        for rn in (rs_v[0]["P-V"]["by_range"] if rs_v else []):
            v = [r["P-V"]["by_range"][rn] for r in rs_v]
            pr = [r["P-R"]["by_range"][rn] for r in rs_v]
            pf = [r["P-F"]["by_range"][rn] for r in rs_v]
            c = lambda k, xs=v: _rate(x["conditions"][k] for x in xs)  # noqa: E731
            L.append(_row(tb, rn, *(_fmt(*_ms([x["errors"][k] for x in v])[:2], d=4) for k in ("M1", "M2", "M3")),
                          c("i_M2_beats_M1"), c("i_M2_beats_M3"), _iv(v, lambda x: x["slope_S1"]),
                          _iv(v, lambda x: x["slope_S2"]), c("ii_S1_negative"), c("ii_S2_negative"),
                          _rate(passed(x) for x in v), ", ".join(_cols(v, lambda x: x["zeta_V"]["estimate"])), z30,
                          ", ".join(_cols(pr, lambda x: x["a_star_T1_T2"]["estimate"])), _rate(passed(x) for x in pr),
                          _rate(passed(x) for x in pf)))
    L += ["", "P-M (tie-aware bins): beta-hat by size bin against the generator's beta_q at the bin's mean log q:", "",
          *_head("unit", "design", "panels evaluable", "stocks", "mean log q by bin", "beta-hat by bin",
                 "beta_q (generator)", "slope [lo, hi]", "top-bin upper < 0", "slope upper < 0", "P-M passes")]
    for tb, rs in groups:
        for key in [W.design_key(d) for d in PM_DIAG]:
            pm = [r["pm"][key] for r in rs if r["pm"][key].get("evaluable") and "beta" in r["pm"][key]]
            if not pm:
                L.append(_row(tb, key, f"0/{len(rs)}", *[""] * 8))
                continue
            lq = [_ms([x["logq"][i] for x in pm])[0] for i in range(len(pm[0]["logq"]))]
            c = lambda k, xs=pm: _rate(x["conditions"][k] for x in xs)  # noqa: E731
            L.append(_row(tb, key, f"{len(pm)}/{len(rs)}", _fmt(_ms([x["stocks"] for x in pm])[0], d=0),
                          ", ".join(f"{v:.2f}" for v in lq), ", ".join(_cols(pm, lambda x: x["beta"]["estimate"])),
                          ", ".join(f"{v:.2f}" for v in paper_beta(lq)), _iv(pm, lambda x: x["slope"]),
                          c("top_bin_upper_below_0"), c("slope_declines_with_log_q"), _rate(passed(x) for x in pm)))
    return L + [""]


def _abias_tables(ab: list) -> list:
    T_grid = ab[0]["T_grid"]
    L = ["## 4. a* bias source (W0, psi = 0.5)", "",
         *_head("world", "tick bin", "panels", *(f"T = {t}" for t in T_grid), "P-R passes", "psi-hat",
                "stocks with psi-hat")]
    for w in ("W0", "W0-shuf"):
        for tb in C.TICK_BINS:
            rs = [r for r in ab if r["world"] == w and r["tick_bin"] == tb]
            pr = [r["P-R"] for r in rs if "P-R" in r]
            if rs:
                L.append(_row(w, tb, len(rs), *_cols(rs, lambda r: r["a_star_all_T"]), _rate(passed(x) for x in pr),
                              _fmt(*_ms([x["psi_hat"]["estimate"] for x in pr])[:2]),
                              _fmt(_ms([x["psi_n_stocks"] for x in pr])[0], d=0)))
    return L + [""]


def _vertex_tables(vx: dict) -> list:
    g = next(iter(vx.values()))["a_grid"]
    L = ["Sizes only, no signs or prices: W0's expected T = 1 curve without sign memory (diag3d2b.vertex).", "",
         *_head("tick bin", "stocks", *(f"R({a:g}) / R(0.5)" for a in g), "argmax (fine grid)",
                "quadratic vertex on [0, 1.5] (a*)", "vertex on [0, 1]")]
    for tb, r in vx.items():
        L.append(_row(tb, r["n_stocks"], *(f"{x:.3f}" for x in r["shape"]), f"{r['argmax_fine']:.2f}",
                      f"{r['vertex_0_1.5']:.3f}", f"{r['vertex_0_1']:.3f}"))
    return L + [""]


def summarise() -> str:
    L = ["# Phase 3d-2b diagnostics: tables", "",
         "Mean over panels with the Monte Carlo SE in brackets; rates are k/n = rate (MC SE). A test passes here "
         "when p <= 0.05 unadjusted and every condition holds (the most lenient Holm step).", ""]
    if (OUT / "occupancy.json").exists():
        L += _occupancy_tables(json.loads((OUT / "occupancy.json").read_text()))
    if pn := load("pnull"):
        L += _pnull_tables(pn)
    if wm := load("wm"):
        L += _wm_tables("3a. Paper II model, the cells' day lengths", wm, C.TICK_BINS)
    if lg := load("long"):
        L += _wm_tables("3b. Paper II model, long days (N = 100000)", lg, ("long",))
    if ab := load("abias"):
        L += _abias_tables(ab)
    if (OUT / "vertex.json").exists():
        L += _vertex_tables(json.loads((OUT / "vertex.json").read_text()))
    text = "\n".join(L)
    (OUT / "tables.md").write_text(text)
    return text


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["occupancy", "vertex", *KINDS, "all", "summarise"])
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    if a.workers > 3:
        raise SystemExit("at most 3 workers")
    t0 = time.time()
    if a.step in ("occupancy", "all"):
        occupancy()
    if a.step in ("vertex", "all"):
        vertex()
    if a.step not in ("occupancy", "vertex", "summarise"):
        run(KINDS if a.step == "all" else (a.step,), a.workers)
    if a.step in ("all", "summarise"):
        print(summarise())
    print(f"{a.step}: {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
