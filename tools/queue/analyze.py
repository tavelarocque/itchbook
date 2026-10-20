#!/usr/bin/env python3
"""M4 analysis: how wrong are L2 queue-position models against L3 truth?

Reads the itch-queue outputs for each day (PREFIX.orders.tsv.gz,
PREFIX.points.tsv.gz, PREFIX.cancels.tsv, and the out-of-sample calibrated
reruns PREFIX.cal_OTHER.orders.tsv.gz) and writes:

  results/m4_models_DAY.tsv      per model: estimation error, fill prediction, timing
  results/m4_strata_DAY.tsv      the same headline numbers by queue length and spread
  results/m4_population.tsv      who was tracked and what happened to them
  results/m4_cancels.tsv         where cancellations come from (both days)
  docs/figures/m4_*.png          figures for the walkthrough

Definitions (all per tracked order, then aggregated across orders):
  relative error   mean over the order's evaluation points of |estimate - truth|,
                   divided by the shares ahead when it joined (ahead0)
  bias             the same with the sign kept: positive = the model thinks more
                   is ahead than really is (pessimistic)
  predicted fill   the model's queue ahead reached zero and a print went past it
                   while the order was alive
  overfill ratio   shares the model would have filled / shares really filled

  analyze.py --data ~/CODE/itch-data/m4 --days 12302019 01302020 --repo .
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

MODELS = ["pessimistic", "optimistic", "uniform", "power2", "power3", "log", "front3", "calibrated"]
FIXED = MODELS[:-1]
QUEUE_BUCKETS = [(1, 1), (2, 4), (5, 9), (10, 24), (25, 99), (100, 10**9)]


def dtypes(columns: list[str]) -> dict[str, str]:
    """Compact types: a day has millions of tracked orders and 60-odd columns."""
    out = {}
    for c in columns:
        if c.endswith(("_abs", "_err", "_err_fill", "_pred_filled")):
            out[c] = "float32"
        elif c.endswith("_pred_trades") or c in ("n_eval", "n_trades", "fill_trades", "orders_ahead0", "own0"):
            out[c] = "int32"
        elif c in ("side", "end"):
            out[c] = "category"
        else:
            out[c] = "int64"
    return out


def read_tsv(path: Path, usecols: list[str] | None = None) -> pd.DataFrame:
    cols = list(pd.read_csv(path, sep="\t", nrows=0).columns)
    if usecols is not None:
        cols = [c for c in cols if c in usecols]
    return pd.read_csv(path, sep="\t", usecols=cols, dtype=dtypes(cols))


def load_orders(path: Path) -> pd.DataFrame:
    df = read_tsv(path)
    tick = np.where(df["price"] >= 10000, 100, 1)  # $0.01 at or above $1, $0.0001 below
    df["spread_ticks"] = df["spread"] / tick
    df["filled_any"] = df["first_fill_ts"] > 0
    df["hour"] = (df["add_ts"] // 3_600_000_000_000).astype("int8")
    df["spread_key"] = pd.Categorical(np.where(df["spread_ticks"] <= 1, "1 tick", ">1 tick"))
    df["queue_bucket"] = pd.cut(
        df["orders_ahead0"],
        bins=[0] + [hi for _, hi in QUEUE_BUCKETS],
        labels=[f"{lo}" if lo == hi else (f"{lo}-{hi}" if hi < 10**9 else f"{lo}+") for lo, hi in QUEUE_BUCKETS],
    )
    return df


def model_metrics(df: pd.DataFrame, models: list[str]) -> pd.DataFrame:
    ev = df[df["n_eval"] > 0]
    filled = df[df["filled_any"]]
    rows = []
    for m in models:
        rel = ev[f"{m}_abs"] / ev["n_eval"] / ev["ahead0"]
        bias = ev[f"{m}_err"] / ev["n_eval"] / ev["ahead0"]
        pred = df[f"{m}_pred_ts"] > 0
        both = filled[filled[f"{m}_pred_ts"] > 0]
        dt = (both[f"{m}_pred_ts"] - both["first_fill_ts"]) / 1e9
        same_print = (both[f"{m}_pred_trades"] == both["fill_trades"]).mean() if len(both) else np.nan
        rows.append(
            {
                "model": m,
                "median_rel_err": rel.median(),
                "mean_rel_err": rel.mean(),
                "mean_bias": bias.mean(),
                "within_10pct": (rel <= 0.10).mean(),
                "err_at_fill_rel_median": (filled[f"{m}_err_fill"] / filled["ahead0"]).median(),
                "pred_fill_rate": pred.mean(),
                "actual_fill_rate": df["filled_any"].mean(),
                "precision": df.loc[pred, "filled_any"].mean() if pred.any() else np.nan,
                "recall": (filled[f"{m}_pred_ts"] > 0).mean() if len(filled) else np.nan,
                "overfill_ratio": df[f"{m}_pred_filled"].sum() / max(df["filled"].sum(), 1),
                "fill_same_print": same_print,
                "fill_early": (dt < 0).mean() if len(dt) else np.nan,
                "fill_dt_median_s": dt.median() if len(dt) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def strata(df: pd.DataFrame, models: list[str]) -> pd.DataFrame:
    out = []
    for name, col in [("queue_orders", "queue_bucket"), ("spread", "spread_key")]:
        for key, g in df.groupby(col, observed=True):
            if len(g) < 200:
                continue
            mm = model_metrics(g, models)
            mm.insert(0, "stratum", f"{name}={key}")
            mm.insert(1, "orders", len(g))
            out.append(mm)
    return pd.concat(out, ignore_index=True)


def load_cancels(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def fig_phi(cancels: dict[str, pd.DataFrame], out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.4), dpi=150)
    r = np.linspace(0, 1, 201)
    ax.plot(r, r, color="#888", lw=1, ls="--", label="uniform: φ = r")
    ax.plot(r, r**3, color="#c66", lw=1, ls=":", label="hftbacktest front3: φ = r³")
    ax.plot(r, r**2 / (r**2 + (1 - r) ** 2), color="#6a6", lw=1, ls="-.", label="power2")
    for (day, c), colour in zip(cancels.items(), ["#1f4e79", "#e08a00"], strict=False):
        mid = (c["lo"] + c["hi"]) / 2
        ax.plot(mid, c["phi"], marker="o", ms=3.5, lw=1.6, color=colour, label=f"measured, {day}")
    ax.set_xlabel("true position of the order: r = shares ahead / other shares at the price")
    ax.set_ylabel("share of cancelled volume that was ahead (φ)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_error_by_queue(per_day: dict[str, pd.DataFrame], models: list[str], out: Path) -> None:
    day, df = next(reversed(per_day.items()))
    s = strata(df, models)
    s = s[s["stratum"].str.startswith("queue_orders")]
    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
    labels = s["stratum"].str.replace("queue_orders=", "").unique()
    for m in models:
        y = s[s["model"] == m]["median_rel_err"].to_numpy()
        ax.plot(labels, y, marker="o", ms=3.5, lw=1.4, label=m)
    ax.set_xlabel("orders ahead when the order joined")
    ax.set_ylabel("median relative error of shares ahead")
    ax.set_title(f"{day}", fontsize=9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_fill(per_day: dict[str, pd.DataFrame], models: list[str], out: Path) -> None:
    """Backtest fills: precision (simulated fills that really happened) against
    recall (real fills the simulation caught), all orders and long queues."""
    day, df = next(reversed(per_day.items()))
    panels = [("all tracked orders", df), ("joined behind 10 or more orders", df[df["orders_ahead0"] >= 10])]
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.9), dpi=150, sharey=True)
    for ax, (title, g) in zip(axes, panels, strict=True):
        mm = model_metrics(g, models)
        ax.scatter(mm["recall"], mm["precision"], s=22, color="#1f4e79", zorder=3)
        for _, r in mm.iterrows():
            ax.annotate(
                r["model"], (r["recall"], r["precision"]), textcoords="offset points", xytext=(4, 3), fontsize=7
            )
        ax.set_title(f"{title} ({len(g):,})", fontsize=9)
        ax.set_xlabel("recall: real fills the model predicted")
        ax.set_xlim(0.3, 1.02)
        ax.set_ylim(0.55, 1.01)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("precision: predicted fills that were real")
    fig.suptitle(day, fontsize=9)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_trajectory(points_path: Path, df: pd.DataFrame, models: list[str], out: Path) -> None:
    """Mean shares ahead (as a share of the queue at entry) after k trade prints."""
    pts = pd.read_csv(points_path, sep="\t")
    pts = pts.merge(df[["ref", "ahead0"]], on="ref", how="inner")
    pts["k"] = pts.groupby("ref").cumcount()
    pts = pts[pts["k"] < 30]
    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
    g = pts.groupby("k")
    truth = (g["truth"].sum() / g["ahead0"].sum()).to_numpy()
    ax.plot(g.size().index, truth, color="black", lw=2.2, label="truth (L3)")
    for m in models:
        if m in pts.columns:
            ax.plot(g.size().index, (g[m].sum() / g["ahead0"].sum()).to_numpy(), lw=1.2, label=m)
    ax.set_xlabel("trade prints at the order's price since it joined")
    ax.set_ylabel("shares ahead / shares ahead at entry")
    ax.set_ylim(0, 1.02)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def population(df: pd.DataFrame, day: str) -> dict[str, object]:
    ends = df["end"].value_counts(normalize=True)
    return {
        "day": day,
        "tracked_orders": len(df),
        "symbols": df["locate"].nunique(),
        "filled_any": df["filled_any"].mean(),
        **{f"end_{k}": ends.get(k, 0.0) for k in ["filled", "cancelled", "deleted", "replaced", "close"]},
        "median_orders_ahead0": df["orders_ahead0"].median(),
        "median_shares_ahead0": df["ahead0"].median(),
        "one_tick_spread": (df["spread_ticks"] <= 1).mean(),
        "median_life_s": ((df["end_ts"] - df["add_ts"]) / 1e9).median(),
        "with_eval_points": (df["n_eval"] > 0).mean(),
        "filled_at_head": (df.loc[df["filled_any"], "truth_at_fill"] == 0).mean(),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--days", nargs="+", required=True)
    ap.add_argument("--repo", type=Path, default=Path("."))
    a = ap.parse_args()
    res = a.repo / "results"
    figs = a.repo / "docs" / "figures"
    res.mkdir(exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    per_day: dict[str, pd.DataFrame] = {}
    calibrated_days: set[str] = set()
    pops = []
    cancels = {}
    for day in a.days:
        per_day.clear()  # only the last day is kept, for the figures
        df = load_orders(a.data / f"{day}.orders.tsv.gz")
        others = [d for d in a.days if d != day]
        cal_path = a.data / f"{day}.cal_{others[0]}.orders.tsv.gz" if others else None
        models = list(FIXED)
        if cal_path and cal_path.exists():
            head = list(pd.read_csv(cal_path, sep="\t", nrows=0).columns)
            cal = read_tsv(cal_path, ["ref"] + [c for c in head if c.startswith("calibrated_")])
            if len(cal) != len(df) or not (cal["ref"].to_numpy() == df["ref"].to_numpy()).all():
                raise SystemExit(f"{cal_path}: tracked orders differ from the uncalibrated run")
            for c in cal.columns:
                if c != "ref":
                    df[c] = cal[c].to_numpy()
            models.append("calibrated")
            calibrated_days.add(day)
        per_day[day] = df
        pops.append(population(df, day))
        model_metrics(df, models).to_csv(res / f"m4_models_{day}.tsv", sep="\t", index=False, float_format="%.4f")
        strata(df, models).to_csv(res / f"m4_strata_{day}.tsv", sep="\t", index=False, float_format="%.4f")
        c = load_cancels(a.data / f"{day}.cancels.tsv")
        c.insert(0, "day", day)
        cancels[day] = c
    pd.DataFrame(pops).to_csv(res / "m4_population.tsv", sep="\t", index=False, float_format="%.4f")
    pd.concat(cancels.values()).to_csv(res / "m4_cancels.tsv", sep="\t", index=False, float_format="%.6f")
    models_all = FIXED + (["calibrated"] if a.days[-1] in calibrated_days else [])
    fig_phi(cancels, figs / "m4_phi.png")
    fig_error_by_queue(per_day, models_all, figs / "m4_error_by_queue.png")
    fig_fill(per_day, models_all, figs / "m4_fill.png")
    last = a.days[-1]
    # Points come from pass 1, where the calibrated column is still the uniform model: plot fixed models only.
    fig_trajectory(a.data / f"{last}.points.tsv.gz", per_day[last], FIXED, figs / "m4_trajectory.png")
    for day in a.days:
        print(f"== {day}")
        print(pd.read_csv(res / f"m4_models_{day}.tsv", sep="\t").to_string(index=False))
    print(pd.DataFrame(pops).T.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
