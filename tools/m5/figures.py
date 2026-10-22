#!/usr/bin/env python3
"""M5 figure: forecast edge per trade against the cost of trading on it.

figures.py --results results --out docs/figures/m5_edge.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COLOURS = {"ridge": "#1f4e79", "deeplob": "#e08a00", "deepofi": "#2e8b57"}
LABELS = {"ridge": "ridge on OFI", "deeplob": "DeepLOB", "deepofi": "deep model on OFI"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=Path("results"))
    ap.add_argument("--out", type=Path, default=Path("docs/figures/m5_edge.png"))
    a = ap.parse_args()
    t = pd.read_csv(a.results / "m5_trading.tsv", sep="\t")
    fr = t[t["costs"].str.startswith("frictionless")]
    models = [m for m in ("ridge", "deeplob", "deepofi") if m in set(fr["model"])]
    horizons = ["1s", "5s", "30s"]
    x = np.arange(len(horizons))
    w = 0.8 / len(models)
    fig, ax = plt.subplots(figsize=(6.6, 4.0), dpi=150)
    for i, m in enumerate(models):
        y = [fr[(fr["model"] == m) & (fr["horizon"] == h)]["gross_per_trade"].iloc[0] for h in horizons]
        ax.bar(x - 0.4 + w / 2 + i * w, y, w * 0.95, color=COLOURS[m], label=LABELS[m])
    ax.axhline(1.0, color="#999", ls="--", lw=1.2)
    ax.text(len(horizons) - 0.5, 1.03, "1-tick spread, no fee", ha="right", va="bottom", fontsize=7.5, color="#666")
    ax.axhline(1.6, color="#c33", ls="--", lw=1.2)
    ax.text(
        len(horizons) - 0.5,
        1.63,
        "1-tick spread + $0.003/share taker fee each way",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#c33",
    )
    ax.set_xticks(x, horizons)
    ax.set_xlabel("forecast horizon")
    ax.set_ylabel("edge per trade, mid to mid (ticks)")
    ax.set_ylim(0, 1.9)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out)
    print(a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
