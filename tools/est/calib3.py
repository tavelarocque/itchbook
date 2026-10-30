#!/usr/bin/env python3
"""Phase 3a calibration of the exploration stocks, per v2 cell (both listings).

CALIBRATION VALUES, not estimands: everything here uses run signs and run
sizes alone (no mid, no price change), read only through explore.py.

Per stock-day (one stock, one study day, runs in time order):
  N             number of runs
  eps_t, q_t    run sign and run volume (shares)

Per cell (all exploration stocks of the cell, all 8 days):
  N quantiles   10/25/50/75/90 % of N over stock-days; N_min per stock (min
                over its 8 days, in calibration_stocks.tsv)
  sizes         log-normal fit by maximum likelihood, mu = mean log q and
                sigma_l = sd log q (denominator n), over all runs of the cell;
                the pooled sizes go to PHASE3/sizes_<cell>.npy for resampling
  log q ACF     lags 1..200: log q demeaned within each stock-day, then
                rho(l) = [sum_sd sum_t x_t x_{t+l} / sum_sd (N - l)] /
                         [sum_sd sum_t x_t^2 / sum_sd N],
                pairs never cross a stock-day; fitted by the copula model
                rho(l) = w rho_H(l) (fGn autocorrelation, see gen3.size_copula)
  sign ACF      C(l) = sum_sd sum_t eps_t eps_{t+l} / sum_sd (N - l), uncentred
                (as in moments.py), pairs within a stock-day only, at log-spaced
                lags; stock-day sd contributes the lags l <= min(2000, N_sd / 10),
                and L_max is the largest lag with at least MIN_STOCKDAYS (20)
                contributors (CHOICE: the cell median N / 10 would stop at lag 1
                to 155, as the median N is 18 to 1557); bands from a cluster
                bootstrap over stock-days (B = 999, percentile 2.5 / 97.5 %)
  gamma         (i) OLS of log C on log l over the fit range [FIT_LO, l_hi],
                l_hi = L_max capped at the last lag before C(l) first
                fails to exceed twice its bootstrap SE, unweighted on the
                log-spaced lags, C(l) = c l^-gamma; bootstrap percentile
                interval from the same replicates; (ii) local Whittle
                (longmem.gamma_lw, m = N^0.65) on every stock-day with N >=
                LW_MIN_N, pooled by inverse variance (longmem.pool_days),
                with Cochran's Q
  fGn signs     the Gaussian ACF r(l) = sin(pi C(l) / 2) that sign(g) needs
                (arcsine law), fitted by gen3.fit_sign_model

    uv run --group analysis tools/est/calib3.py collect   # one pass over the data
    uv run --group analysis tools/est/calib3.py fit       # table, summary, curves
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

import explore as E
import gen3 as G
import longmem as LM

REPO = Path(__file__).resolve().parents[2]
PHASE3 = E.data_dir() / "phase3"
OUT_TSV = REPO / "data" / "phase3" / "calibration.tsv"
OUT_STOCKS = REPO / "data" / "phase3" / "calibration_stocks.tsv"
OUT_MD = REPO / "data" / "phase3" / "calibration.md"
SIGN_LAG_MAX = 2000
LOGQ_LAG_MAX = 200
FIT_LO = 10
LW_MIN_N = 500
MIN_STOCKDAYS = 20
N_BOOT = 999
BOOT_SEED = 20261101


def sign_lags(l_max: int = SIGN_LAG_MAX, per_decade: int = 12) -> np.ndarray:
    """1..10, then 10^(k / per_decade) rounded (exact at 10, 100, 1000), and l_max."""
    k = np.arange(int(np.floor(np.log10(l_max) * per_decade)) + 1)
    g = np.round(10.0 ** (k / per_decade)).astype(int)
    return np.unique(np.concatenate([np.arange(1, 11), g[g <= l_max], [l_max]]))


def lag_products(x: np.ndarray, lags: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """sum_t x_t x_{t+l} and the pair count N - l (0 when l >= N), within one series."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    num = np.zeros(len(lags))
    cnt = np.maximum(n - lags, 0).astype(float)
    use = lags < n
    if n > 1 and use.any():
        size = 1 << int(np.ceil(np.log2(2 * n)))
        f = np.fft.rfft(x, size)
        ac = np.fft.irfft(f * np.conj(f), size)
        num[use] = ac[lags[use]]
    return num, cnt


def pooled_acf(series, lags: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """C(l) pooled over series without crossing a series boundary, and the per-series
    (num, cnt) rows the cluster bootstrap resamples."""
    num = np.zeros((len(series), len(lags)))
    cnt = np.zeros_like(num)
    for i, x in enumerate(series):
        num[i], cnt[i] = lag_products(x, lags)
    with np.errstate(invalid="ignore", divide="ignore"):
        return num.sum(0) / cnt.sum(0), num, cnt


def cluster_boot(num: np.ndarray, cnt: np.ndarray, n_boot: int, rng: np.random.Generator) -> np.ndarray:
    """Replicates of sum num / sum cnt with the rows (stock-days) drawn with replacement."""
    n = len(num)
    out = np.empty((n_boot, num.shape[1]))
    for b in range(n_boot):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[b] = (w @ num) / (w @ cnt)
    return out


def power_fit(l: np.ndarray, C: np.ndarray) -> tuple[float, float]:
    """OLS of log C on log l: (c, gamma) with C = c l^-gamma."""
    b, a = np.polyfit(np.log(l), np.log(C), 1)
    return float(np.exp(a)), float(-b)


# ---- pass 1: collect ------------------------------------------------------------


def collect(log=print) -> None:
    PHASE3.mkdir(parents=True, exist_ok=True)
    days, u = E.read_universe()
    cell_of = {s: c for s, (c, sp) in u.items() if sp == "explore"}
    lags = sign_lags()
    qlags = np.arange(0, LOGQ_LAG_MAX + 1)
    acc = {}
    for d in days:
        t0 = time.time()
        runs = E.load_explore_runs(d)
        for c in sorted({cell_of[s] for s in runs}):
            syms = sorted(s for s in runs if cell_of[s] == c)
            a = acc.setdefault(c, {"sym": [], "day": [], "N": [], "q": [], "snum": [], "scnt": [],
                                   "qnum": np.zeros(len(qlags)), "qcnt": np.zeros(len(qlags)),
                                   "lw_g": [], "lw_se": []})
            eps = [runs[s]["sign"].astype(float) for s in syms]
            lq = [np.log(runs[s]["q"].astype(float)) for s in syms]
            _, num, cnt = pooled_acf(eps, lags)  # one row per stock-day
            a["snum"] += list(num)
            a["scnt"] += list(cnt)
            _, num, cnt = pooled_acf([x - x.mean() for x in lq], qlags)
            a["qnum"] += num.sum(0)
            a["qcnt"] += cnt.sum(0)
            for s, e in zip(syms, eps, strict=True):
                a["sym"].append(s)
                a["day"].append(d)
                a["N"].append(len(e))
                a["q"].append(runs[s]["q"].astype(np.uint32))
                if len(e) >= LW_MIN_N:
                    g, se, _ = LM.gamma_lw(e)
                    a["lw_g"].append(g)
                    a["lw_se"].append(se)
        log(f"{d}: {len(runs)} stocks, {sum(len(v) for v in runs.values()):,} runs, {time.time() - t0:.1f}s")
        del runs
    for c, a in sorted(acc.items()):
        np.save(PHASE3 / f"sizes_{c}.npy", np.concatenate(a["q"]))
        np.savez(PHASE3 / f"calib_raw_{c}.npz", sym=np.array(a["sym"]), day=np.array(a["day"]),
                 N=np.array(a["N"]), sign_lags=lags, snum=np.array(a["snum"]), scnt=np.array(a["scnt"]),
                 qlags=qlags, qnum=a["qnum"], qcnt=a["qcnt"], lw_g=np.array(a["lw_g"]),
                 lw_se=np.array(a["lw_se"]))
    log(f"wrote {len(acc)} cells to {PHASE3}")


# ---- pass 2: fit ------------------------------------------------------------------


def sign_acf_fit(lags: np.ndarray, snum: np.ndarray, scnt: np.ndarray, N: np.ndarray, rng: np.random.Generator,
                 n_boot: int = N_BOOT) -> dict:
    """The calibration's sign-ACF rule on per-stock-day lag products (rows = stock-days): each stock-day
    contributes l <= min(SIGN_LAG_MAX, N / 10), L_max = last lag with >= MIN_STOCKDAYS contributors,
    stock-day cluster bootstrap, power-law OLS on [FIT_LO, l_hi] with l_hi cut where C <= 2 SE."""
    use = lags[None, :] <= np.minimum(SIGN_LAG_MAX, N[:, None] // 10)  # stock-day sd contributes l <= N_sd / 10
    n_sd = use.sum(0)
    L_max = int(lags[n_sd >= MIN_STOCKDAYS].max())
    sel = lags <= L_max
    lags, use, n_sd = lags[sel], use[:, sel], n_sd[sel]
    num, cnt = np.where(use, snum[:, sel], 0.0), np.where(use, scnt[:, sel], 0.0)
    C = num.sum(0) / cnt.sum(0)
    reps = cluster_boot(num, cnt, n_boot, rng)
    lo, hi = np.nanpercentile(reps, [2.5, 97.5], axis=0)
    se = np.nanstd(reps, axis=0)
    weak = np.flatnonzero((lags >= FIT_LO) & ~(C > 2 * se))
    l_hi = int(lags[weak[0] - 1]) if len(weak) else L_max
    fit = (lags >= FIT_LO) & (lags <= l_hi)
    out = {"lags": lags, "n_sd": n_sd, "C": C, "lo": lo, "hi": hi, "se": se, "L_max": L_max, "l_hi": l_hi,
           "fit_n_lags": int(fit.sum())}
    if fit.sum() >= 3 and np.all(C[fit] > 0):
        cc, gg = power_fit(lags[fit], C[fit])
        g_reps = [power_fit(lags[fit], r[fit])[1] for r in reps if np.all(r[fit] > 0)]
        out.update(c_fit=cc, gamma_fit=gg, gamma_fit_nboot=len(g_reps))
        out["gamma_fit_lo"], out["gamma_fit_hi"] = (float(v) for v in np.percentile(g_reps, [2.5, 97.5]))
    else:
        out.update(c_fit=np.nan, gamma_fit=np.nan, gamma_fit_lo=np.nan, gamma_fit_hi=np.nan, gamma_fit_nboot=0)
    return out


def fit_cell(c: str, rng: np.random.Generator) -> tuple[dict, dict]:
    z = np.load(PHASE3 / f"calib_raw_{c}.npz")
    q = np.load(PHASE3 / f"sizes_{c}.npy").astype(float)
    N = z["N"]
    lq = np.log(q)
    row = {"cell": c, "listing": "NASDAQ" if c.startswith("Q") else "NYSE",
           "n_stocks": len(np.unique(z["sym"])), "n_stockdays": len(N), "n_runs": int(N.sum())}
    for p in (10, 25, 50, 75, 90):
        row[f"N_q{p}"] = float(np.percentile(N, p))
    per_stock = {s: int(N[z["sym"] == s].min()) for s in np.unique(z["sym"])}
    row["N_min_min"] = min(per_stock.values())
    row["N_min_med"] = float(np.median(list(per_stock.values())))
    row["mu_logq"] = float(lq.mean())
    row["sigma_l"] = float(lq.std())
    # log q ACF
    qacf = (z["qnum"] / z["qcnt"]) / (z["qnum"][0] / z["qcnt"][0])
    for l in (1, 10, 100):
        row[f"logq_acf_{l}"] = float(qacf[l])
    w, H = G.fit_size_copula(qacf[1:], np.arange(1, LOGQ_LAG_MAX + 1))
    row["copula_w"], row["copula_H"] = w, H
    # sign ACF
    sf = sign_acf_fit(z["sign_lags"], z["snum"], z["scnt"], N, rng, N_BOOT)
    lags, n_sd, C, lo, hi, se, L_max = (sf[k] for k in ("lags", "n_sd", "C", "lo", "hi", "se", "L_max"))
    row["L_max"] = L_max
    for l in (1, 10, 100):
        row[f"C_{l}"] = float(C[lags == l][0]) if l in lags else np.nan
    row["fit_lo"], row["fit_hi"], row["fit_n_lags"] = FIT_LO, sf["l_hi"], sf["fit_n_lags"]
    row.update({k: sf[k] for k in ("c_fit", "gamma_fit", "gamma_fit_lo", "gamma_fit_hi", "gamma_fit_nboot")})
    if len(z["lw_g"]) >= 2:
        g, s, Q, df = LM.pool_days(z["lw_g"], z["lw_se"])
        row.update(gamma_lw=g, gamma_lw_se=s, lw_Q=Q, lw_df=df, lw_n=len(z["lw_g"]))
    else:
        row.update(gamma_lw=np.nan, gamma_lw_se=np.nan, lw_Q=np.nan, lw_df=0, lw_n=len(z["lw_g"]))
    # fGn sign model on the measured C(l) over 1..L_max
    sm = G.fit_sign_model(lags, C)
    row.update({f"sgn_{k}": v for k, v in sm.items()})
    curves = {"sign_lags": lags, "sign_n_stockdays": n_sd, "C": C, "C_lo": lo, "C_hi": hi, "C_se": se,
              "logq_lags": np.arange(LOGQ_LAG_MAX + 1),
              "logq_acf": qacf, "N": N, "sym": z["sym"], "day": z["day"]}
    return row, curves | {"per_stock": per_stock}


COLS = ["cell", "listing", "n_stocks", "n_stockdays", "n_runs", "N_q10", "N_q25", "N_q50", "N_q75", "N_q90",
        "N_min_min", "N_min_med", "mu_logq", "sigma_l", "logq_acf_1", "logq_acf_10", "logq_acf_100", "copula_w",
        "copula_H", "L_max", "C_1", "C_10", "C_100", "fit_lo", "fit_hi", "fit_n_lags", "c_fit", "gamma_fit",
        "gamma_fit_lo", "gamma_fit_hi", "gamma_fit_nboot", "gamma_lw", "gamma_lw_se", "lw_Q", "lw_df", "lw_n",
        "sgn_w1", "sgn_H", "sgn_w2", "sgn_phi", "sgn_rmse"]


def _fmt(v) -> str:
    if isinstance(v, (float, np.floating)):
        return "" if not np.isfinite(v) else f"{v:.6g}"
    return str(v)


def fit(log=print) -> None:
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    rows, stocks = [], []
    rng = np.random.default_rng(BOOT_SEED)
    for c in E.cells():
        row, cur = fit_cell(c, rng)
        rows.append(row)
        stocks += [(s, c, n) for s, n in sorted(cur.pop("per_stock").items())]
        np.savez(PHASE3 / f"calibration_{c}.npz", **cur, row=json.dumps({k: _fmt(v) for k, v in row.items()}))
        log(f"{c}: N50 {row['N_q50']:.0f} sigma_l {row['sigma_l']:.3f} gamma fit {row['gamma_fit']:.3f} "
            f"[{row['fit_lo']}, {row['fit_hi']}] LW {row['gamma_lw']:.3f}")
    with open(OUT_TSV, "w") as f:
        f.write("# Phase 3a CALIBRATION values (not estimands): run signs and sizes of exploration stocks only, "
                "8 study days, per v2 cell (tools/est/calib3.py)\n")
        f.write("# sign ACF uncentred, pooled within cell and within stock-day, stock-day lags <= "
                f"min({SIGN_LAG_MAX}, N / 10), "
                f"L_max = last lag with >= {MIN_STOCKDAYS} stock-days; "
                f"power fit OLS on log-spaced lags [{FIT_LO}, fit_hi]; bootstrap B = {N_BOOT} over stock-days, "
                f"default_rng({BOOT_SEED}); local Whittle m = N^0.65 on stock-days with N >= {LW_MIN_N}\n")
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(r[k]) for k in COLS) + "\n")
    with open(OUT_STOCKS, "w") as f:
        f.write("# N_min = fewest runs over the 8 study days, exploration stocks (tools/est/calib3.py)\n")
        f.write("symbol\tcell\tN_min\n")
        for s, c, n in stocks:
            f.write(f"{s}\t{c}\t{n}\n")
    write_md(rows)
    log(f"wrote {OUT_TSV}, {OUT_STOCKS}, {OUT_MD}")


def write_md(rows) -> None:
    def f(v, d=3):
        return "" if not np.isfinite(v) else f"{v:.{d}f}"

    lines = ["# Phase 3a calibration summary", "",
             "CALIBRATION values from exploration stocks only (split_v2 = explore), 8 study days, run signs and",
             "sizes alone. Produced by `tools/est/calib3.py`; full table in `calibration.tsv`, per-stock N_min in",
             "`calibration_stocks.tsv`, curves in `$ITCH_DATA/phase3/calibration_<cell>.npz`.", "",
             "| cell | stocks | stock-days | N10 | N50 | N90 | sigma_l | rho_logq(1) | rho_logq(10) | rho_logq(100) "
             "| C(1) | fit range | gamma fit [95%] | gamma LW (SE) |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(
            f"| {r['cell']} | {r['n_stocks']} | {r['n_stockdays']} | {r['N_q10']:.0f} | {r['N_q50']:.0f} | "
            f"{r['N_q90']:.0f} | {f(r['sigma_l'])} | {f(r['logq_acf_1'])} | {f(r['logq_acf_10'])} | "
            f"{f(r['logq_acf_100'])} | {f(r['C_1'])} | {r['fit_lo']}..{r['fit_hi']} | {f(r['gamma_fit'])} "
            f"[{f(r['gamma_fit_lo'])}, {f(r['gamma_fit_hi'])}] | {f(r['gamma_lw'])} ({f(r['gamma_lw_se'])}) |")
    g = np.array([r["gamma_fit"] for r in rows], float)
    lw = np.array([r["gamma_lw"] for r in rows], float)
    edge = [r["cell"] for r in rows if r["sgn_H"] >= 0.995 - 1e-9]
    lines += ["", "Measurements:", "",
              f"- gamma by the power fit spans {np.nanmin(g):.2f} to {np.nanmax(g):.2f}; by local Whittle "
              f"{np.nanmin(lw):.2f} to {np.nanmax(lw):.2f}. The uncentred pooled C(l) carries the variance of "
              "each stock-day's mean sign (a near-flat tail beyond lag 10), which local Whittle, using "
              "frequencies j >= 1 only, does not see.",
              f"- The sign model's fGn component sits on the H = 0.995 grid edge in {len(edge)} of {len(rows)} "
              f"cells ({', '.join(edge)}): a near-constant per-day sign bias plus short memory."]
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["collect", "fit"])
    a = ap.parse_args()
    collect() if a.step == "collect" else fit()
    return 0


if __name__ == "__main__":
    sys.exit(main())
