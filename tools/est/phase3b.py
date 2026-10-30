#!/usr/bin/env python3
"""Phase 3b driver: sign model, clipping, usable a range, GT under clipping.

Real data (signs and sizes of exploration stocks, 8 study days) are read once,
through explore.load_explore_runs, into PHASE3/p3b/raw_<cell>.npz (per
stock-day eps, q, N, symbol, day). Every step checkpoints per cell under
PHASE3/p3b/ and skips cells already done.

  signs   C(l) uncentred and centred on each stock-day's own mean (lags of
          calib3.sign_lags up to the cell's L_max), var of the daily mean sign
          and its scaling with N, the two-component fit next to the 3a
          one-component fit; fit quality as chi2/df = mean over LOG_LAGS of
          ((C_model - C) / SE)^2 with the calibration's bootstrap SE.
  est     the gamma estimator study on the fitted two-component generator,
          N_REP panels per cell with the cell's own stock-day N values.
  clip    removed shares per rule, the R1 share-clipped slope on log N, the
          Kish ESS tables and a_max.
  gt      GT ground truth with sizes clipped at the cell's log-normal 99.9 %
          quantile (R2 in the simulation), NASDAQ cells, paired with the 3a run
          (same seed), and the 3a values on the same a range.

    uv run --group analysis tools/est/phase3b.py collect|signs|est|clip|gt|all [--workers 4]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import calib3 as C3
import clip3b as CL
import explore as E
import gen3 as G
import signs3b as SG
import truth3 as TR

P3B = E.data_dir() / "phase3" / "p3b"
REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "data" / "phase3"
N_REP = 50
SEED = 20261103
FAIL_CHI2 = 4.0
Z999 = 3.090232306167813  # standard normal 0.999 quantile
GT_A = np.round(np.arange(0.0, 1.501, 0.25), 2)
GT_RANGES = ((8, 64), (16, 256))
TAIL_LO = 10


def _fmt(v) -> str:
    if isinstance(v, (bool, np.bool_)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return "" if not np.isfinite(v) else f"{v:.5g}"
    return str(v)


def write_tsv(path: Path, rows: list[dict], header: str) -> None:
    cols = list(rows[0])
    with open(path, "w") as f:
        f.write(f"# {header}\n" + "\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(r[c]) for c in cols) + "\n")


# ---- data -------------------------------------------------------------------------


def collect(log=print) -> None:
    P3B.mkdir(parents=True, exist_ok=True)
    days, u = E.read_universe()
    cell_of = {s: c for s, (c, sp) in u.items() if sp == "explore"}
    acc: dict[str, dict] = {}
    for d in days:
        runs = E.load_explore_runs(d)
        for s in sorted(runs):
            a = acc.setdefault(cell_of[s], {"eps": [], "q": [], "N": [], "sym": [], "day": []})
            a["eps"].append(runs[s]["sign"].astype(np.int8))
            a["q"].append(runs[s]["q"].astype(np.uint32))
            a["N"].append(len(runs[s]))
            a["sym"].append(s)
            a["day"].append(d)
        log(f"{d}: {len(runs)} stocks")
    for c, a in acc.items():
        np.savez(P3B / f"raw_{c}.npz", eps=np.concatenate(a["eps"]), q=np.concatenate(a["q"]), N=np.array(a["N"]),
                 sym=np.array(a["sym"]), day=np.array(a["day"]))


def load_raw(c: str):
    z = np.load(P3B / f"raw_{c}.npz")
    cut = np.cumsum(z["N"])[:-1]
    return (np.split(z["eps"].astype(float), cut), np.split(z["q"].astype(float), cut), z["N"], z["sym"],
            z["day"])


# ---- 1. sign model ----------------------------------------------------------------


def sign_cell(c: str) -> dict:
    eps, _, N, sym, _ = load_raw(c)
    cal = np.load(G.phase3_dir() / f"calibration_{c}.npz")
    row3a = json.loads(str(cal["row"]))
    lags = cal["sign_lags"]
    se = cal["C_se"]
    allL = C3.sign_lags()
    use = allL[None, :] <= np.minimum(C3.SIGN_LAG_MAX, N[:, None] // 10)
    _, num, cnt = C3.pooled_acf(eps, allL)
    _, numc, cntc = SG.centred_acf(eps, allL)
    sel = np.isin(allL, lags)
    Cu = np.where(use, num, 0).sum(0)[sel] / np.where(use, cnt, 0).sum(0)[sel]
    Cc = np.where(use, numc, 0).sum(0)[sel] / np.where(use, cntc, 0).sum(0)[sel]
    assert np.allclose(Cu, cal["C"]), "uncentred C differs from the calibration"
    means = np.array([x.mean() for x in eps])
    ms = SG.mean_sign_scaling(means, N)
    fit = SG.fit_two_component(lags, Cu)
    sm = {k: float(row3a[f"sgn_{k}"]) for k in ("w1", "H", "w2", "phi")}
    fl = np.isin(lags, SG.LOG_LAGS)
    C2 = SG.two_component_C(lags, fit["H"], fit["tau"])
    C1 = G.sign_model_C(lags, **sm)
    row = {"cell": c, "n_stockdays": len(N), "L_max": int(lags.max()), "var_mean_sign": ms["var_all"],
           "var_slope_logN": ms["slope"], "var_intercept": ms["intercept"], "var_n_bins": len(ms["bins"])}
    for l in (1, 10, 100):
        row[f"Cu_{l}"] = float(Cu[lags == l][0]) if l in lags else np.nan
        row[f"Cc_{l}"] = float(Cc[lags == l][0]) if l in lags else np.nan
    row.update({f"two_{k}": fit[k] for k in ("H", "tau", "gamma_true", "rmse", "edge_first", "edge_final",
                                              "tau_edge")})
    row["two_plateau"] = float((2 / np.pi) * np.arcsin(fit["tau"] ** 2 / (1 + fit["tau"] ** 2)))
    tail = lags >= TAIL_LO  # sensitivity: the same fit on the tail lags only
    ft = SG.fit_two_component(lags[tail], Cu[tail])
    row.update({f"tail_{k}": ft[k] for k in ("H", "tau", "gamma_true", "rmse", "edge_first", "edge_final")})
    row["one_rmse"] = float(np.sqrt(np.mean((C1[fl] - Cu[fl]) ** 2)))
    row["two_chi2df"] = float(np.mean(((C2[fl] - Cu[fl]) / se[fl]) ** 2))
    row["one_chi2df"] = float(np.mean(((C1[fl] - Cu[fl]) / se[fl]) ** 2))
    row["two_fails"] = row["two_chi2df"] > FAIL_CHI2
    row["one_fails"] = row["one_chi2df"] > FAIL_CHI2
    worst = int(np.argmax(np.abs(C2[fl] - Cu[fl])))
    row["two_worst_lag"], row["two_worst_resid"] = int(lags[fl][worst]), float((C2 - Cu)[fl][worst])
    np.savez(P3B / f"signs_{c}.npz", lags=lags, Cu=Cu, Cc=Cc, se=se, C_two=C2, C_one=C1, bins=np.array(ms["bins"]),
             row=json.dumps({k: _fmt(v) for k, v in row.items()}))
    return row


# ---- 1c. estimator study -----------------------------------------------------------


def est_cell(c: str) -> dict:
    srow = json.loads(str(np.load(P3B / f"signs_{c}.npz")["row"]))
    eps, _, N, _, _ = load_raw(c)
    real = SG.gamma_estimators(eps, np.random.default_rng(SEED))
    row = {"cell": c} | {f"{k}_real": real[k] for k in ("lw", "power", "v0")}
    for fit in ("two", "tail"):  # the all-lags fit (the spec's) and the tail-lags sensitivity fit
        H, tau = float(srow[f"{fit}_H"]), float(srow[f"{fit}_tau"])
        path = P3B / f"est_{fit}_{c}.npz"
        if path.exists():
            est = np.load(path)["est"]
        else:
            est = SG.estimator_study(N, H, tau, N_REP, SEED)
            np.savez(path, est=est)
        row |= {f"{fit}_H": H, f"{fit}_tau": tau, f"{fit}_gamma_true": 2 - 2 * H}
        for i, k in enumerate(("lw", "power", "v0")):
            x = est[:, i][np.isfinite(est[:, i])]
            row[f"{fit}_{k}_mean"] = float(x.mean()) if len(x) else np.nan
            row[f"{fit}_{k}_bias"] = float(x.mean() - (2 - 2 * H)) if len(x) else np.nan
            row[f"{fit}_{k}_sd"] = float(x.std(ddof=1)) if len(x) > 1 else np.nan
            row[f"{fit}_{k}_n"] = len(x)
    return row


# ---- 2 and 3. clipping and ESS ------------------------------------------------------


def clip_cell(c: str) -> tuple[list[dict], list[dict], dict]:
    _, qs, N, sym, _ = load_raw(c)
    clip_rows, ess_rows = [], []
    r1_share = None
    for rule in CL.RULES:
        qcs = CL.apply_rule(rule, qs, sym)
        row = {"cell": c, "rule": rule} | CL.removed_shares(qs, qcs)
        clip_rows.append(row)
        if rule == "R1":
            r1_share = np.array([np.mean(qc < q) for q, qc in zip(qs, qcs, strict=True)])
        V, K = CL.ess_table(qcs)
        pooled = {s: np.concatenate([qc for qc, t in zip(qcs, sym, strict=True) if t == s]) for s in np.unique(sym)}
        Vp, Kp = CL.ess_table(list(pooled.values()))
        for kind, X, Xp in (("V", V, Vp), ("K", K, Kp)):
            med, p10, medp = np.median(X, 0), np.percentile(X, 10, axis=0), np.median(Xp, 0)
            er = {"cell": c, "rule": rule, "weights": kind}
            er |= {f"med_{a:g}": m for a, m in zip(CL.A_GRID, med, strict=True)}
            er |= {f"p10_{a:g}": m for a, m in zip(CL.A_GRID, p10, strict=True)}
            for thr in (10, 30, 100):
                er[f"a_max_{thr}"] = CL.a_max(med, thr)
            for thr in (10, 30, 100):
                er[f"a_max_pooled_{thr}"] = CL.a_max(medp, thr)
            ess_rows.append(er)
    slope, lo, hi = CL.share_clipped_slope(r1_share, N, sym, 999, np.random.default_rng(SEED))
    r1 = {"cell": c, "r1_slope_logN": slope, "r1_slope_lo": lo, "r1_slope_hi": hi,
          "r1_share_N_le_100": float(np.mean(r1_share[N <= 100])) if np.any(N <= 100) else np.nan,
          "r1_share_N_gt_1000": float(np.mean(r1_share[N > 1000])) if np.any(N > 1000) else np.nan,
          "r1_stockdays_all_clipped_frac": float(np.mean(r1_share[N <= 100] > 0)) if np.any(N <= 100) else np.nan}
    return clip_rows, ess_rows, r1


# ---- 4. GT under clipping -----------------------------------------------------------


def depth15(zK: np.ndarray) -> float:
    z = zK[TR.A_GRID <= 1.5]
    return float(np.min(z) - min(z[0], z[-1]))


def gt_summary(c, r) -> dict:
    out = {}
    se = TR.jackknife(*r.days, fn=lambda cu: {f"d_{lo}_{hi}": depth15(cu[f"zK_{lo}_{hi}"]) for lo, hi in GT_RANGES})
    for lo, hi in GT_RANGES:
        z = r.curves[f"zK_{lo}_{hi}"][TR.A_GRID <= 1.5]
        out[f"depth_{lo}_{hi}"] = depth15(r.curves[f"zK_{lo}_{hi}"])
        out[f"depth_se_{lo}_{hi}"] = float(se[f"d_{lo}_{hi}"])
        out[f"argmin_{lo}_{hi}"] = float(GT_A[int(np.argmin(z))])
        out[f"zKmin_{lo}_{hi}"] = float(np.min(z))
        out[f"zK_se_max_{lo}_{hi}"] = float(np.max(r.se[f"zK_{lo}_{hi}"][TR.A_GRID <= 1.5]))
    return out


class _Saved:  # a 3a truth run read back from its npz, shaped like TruthResult
    def __init__(self, path):
        z = np.load(path)
        self.days = (z["Vd"], z["Kd"], z["Sd"])
        self.curves = TR.curves(*(a.mean(0) for a in self.days))
        self.se = TR.jackknife(*self.days)
        self.n_days, self.se_max = len(z["Vd"]), TR.max_slope_se(self.se)


def gt_clip(workers: int, log=print) -> list[dict]:
    rows = []
    cells = [c for c in E.cells() if c.startswith("Q-")]
    log(f"GT under R2 clipping: {len(cells)} cells, at most 64 days each, about 1 min per cell")
    for c in cells:
        cl = G.load_cell(c)
        g = TR.GT_TRUTH
        params = {"beta_max": g["beta_max"], "c": g["c"], "l_star": cl.mu + g["l_off"] * cl.sigma_l, "psi": g["psi"]}
        cap = float(np.exp(cl.mu + Z999 * cl.sigma_l))
        r = TR.run_truth("GT", c, params | {"q_cap": cap}, name=f"GT_R2_{c}", min_days=16, chunk=16, max_days=64,
                         workers=workers, out_dir=P3B / "truth", n_collapse=0, log=log)
        base = _Saved(TR.TRUTH_DIR / f"GT_{c}.npz")
        row = {"cell": c, "q_cap": cap, "n_days": r.n_days, "se_max_domain": r.se_max, "se_met": r.se_met}
        row |= {f"clip_{k}": v for k, v in gt_summary(c, r).items()}
        row |= {"base_n_days": base.n_days, "base_se_max_domain": base.se_max}
        row |= {f"base_{k}": v for k, v in gt_summary(c, base).items()}
        for lo, hi in GT_RANGES:
            row[f"base_depth_fullgrid_{lo}_{hi}"] = TR.v_depth(base.curves[f"zK_{lo}_{hi}"])
            row[f"clip_zK_{lo}_{hi}"] = " ".join(f"{v:.3f}" for v in r.curves[f"zK_{lo}_{hi}"][TR.A_GRID <= 1.5])
            row[f"base_zK_{lo}_{hi}"] = " ".join(f"{v:.3f}" for v in base.curves[f"zK_{lo}_{hi}"][TR.A_GRID <= 1.5])
        rows.append(row)
        log(f"{c}: clipped depth [8,64] {row['clip_depth_8_64']:+.3f} ({row['clip_depth_se_8_64']:.3f}) "
            f"vs {row['base_depth_8_64']:+.3f}")
    return rows


# ---- main ---------------------------------------------------------------------------


def _pool(fn, items, workers):
    if workers <= 1:
        return [fn(x) for x in items]
    with ProcessPoolExecutor(workers) as ex:
        return list(ex.map(fn, items))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["collect", "signs", "est", "clip", "gt", "all"])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    steps = ["collect", "signs", "est", "clip", "gt"] if a.step == "all" else [a.step]
    cells = E.cells()
    for st in steps:
        t0 = time.time()
        if st == "collect":
            collect()
        elif st == "signs":
            write_tsv(OUT / "p3b_signs.tsv", _pool(sign_cell, cells, a.workers),
                      "Phase 3b sign model per cell (tools/est/phase3b.py signs): exploration signs only")
        elif st == "est":
            write_tsv(OUT / "p3b_gamma_estimators.tsv", _pool(est_cell, cells, a.workers),
                      f"Phase 3b gamma estimators on the fitted two-component generator, {N_REP} panels per cell "
                      "with the cell's stock-day N values; *_real = the same estimator on the exploration signs")
        elif st == "clip":
            res = _pool(clip_cell, cells, a.workers)
            write_tsv(OUT / "p3b_clip.tsv", [r for cr, _, _ in res for r in cr],
                      "Phase 3b clipping: shares removed, pooled over the cell's runs (exploration sizes only)")
            write_tsv(OUT / "p3b_r1_slope.tsv", [r1 for _, _, r1 in res],
                      "Phase 3b R1: OLS slope of the stock-day share of runs clipped on log N, 95 % cluster "
                      "bootstrap over stocks (B = 999)")
            write_tsv(OUT / "p3b_ess.tsv", [r for _, er, _ in res for r in er],
                      "Phase 3b Kish ESS per stock-day: median and 10th percentile by a; a_max = largest a with "
                      "median ESS >= threshold; pooled = ESS of each stock's 8 days together, median over stocks")
        elif st == "gt":
            write_tsv(OUT / "p3b_gt_clip.tsv", gt_clip(a.workers),
                      "Phase 3b GT ground truth with sizes clipped at exp(mu + 3.0902 sigma_l) (R2 in simulation), "
                      "a <= 1.5; base = the unclipped Phase 3a run (same seed, same days)")
        print(f"{st}: {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
