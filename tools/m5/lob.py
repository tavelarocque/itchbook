"""M5 features and targets from itch-lob arrays.

One symbol-day is an int32 array of shape (rows, 56) on a fixed 250 ms grid
(layout in include/itchbook/lob/recorder.hpp). Everything here is computed on
the fly from it; nothing is cached, because the raw arrays already take most
of the disk.

Three views of each row t:
  lob(t)     the DeepLOB input: 10 levels x (ask price, ask size, bid price,
             bid size), prices in ticks from the mid, sizes as log1p(round lots)
  ofi(t)     the order-flow input: multi-level OFI and signed volume in the
             interval ending at t, scaled by the symbol's typical best-level depth
  linear(t)  the baseline's features: OFI and signed volume summed over several
             look-back windows, book imbalance, spread
Targets: the mid-price change from t to t+h in ticks, and its sign.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DEPTH = 10
BID_PX, BID_SZ, ASK_PX, ASK_SZ, OFI = 0, 10, 20, 30, 40
BUY_VOL, SELL_VOL, HIDDEN_B, HIDDEN_S, EVENTS, TRADING = 50, 51, 52, 53, 54, 55
GRID_S = 0.25
HORIZONS = {"1s": 4, "5s": 20, "30s": 120}
WINDOWS = (1, 4, 20, 120)  # look-backs for the linear model, in grid steps
WARMUP = 240  # skip the first minute: windows and the deep model's 100-step input need history
LOOKBACK = 100  # DeepLOB input length (25 s)


def tick_of(price: np.ndarray) -> np.ndarray:
    """US tick size in ITCH price units (1e-4 dollars): $0.01 at or above $1."""
    return np.where(price >= 10000, 100, 1)


@dataclass
class SymbolDay:
    symbol: str
    day: str
    raw: np.ndarray  # int32 (rows, 56)

    @classmethod
    def load(cls, data: Path, day: str, symbol: str) -> SymbolDay:
        return cls(symbol, day, np.load(data / f"{day}.{symbol}.npy"))

    @property
    def bid(self) -> np.ndarray:
        return self.raw[:, BID_PX].astype(np.float64)

    @property
    def ask(self) -> np.ndarray:
        return self.raw[:, ASK_PX].astype(np.float64)

    def valid(self) -> np.ndarray:
        r = self.raw
        return (r[:, TRADING] == 1) & (r[:, BID_PX] > 0) & (r[:, ASK_PX] > r[:, BID_PX])

    def mid(self) -> np.ndarray:
        return (self.bid + self.ask) / 2

    def tick(self) -> np.ndarray:
        return tick_of(self.mid()).astype(np.float64)

    def depth_scale(self) -> float:
        """Typical best-level depth (shares): the OFI and volume scale."""
        v = self.valid()
        d = (self.raw[v, BID_SZ] + self.raw[v, ASK_SZ]) / 2
        return float(np.median(d)) if d.size else 1.0

    # ---- targets --------------------------------------------------------------
    def targets(self, h: int) -> tuple[np.ndarray, np.ndarray]:
        """Mid change t -> t+h in ticks, and a mask of rows where it is defined."""
        n = len(self.raw)
        mid, tick, v = self.mid(), self.tick(), self.valid()
        y = np.zeros(n)
        ok = np.zeros(n, dtype=bool)
        y[: n - h] = (mid[h:] - mid[: n - h]) / tick[: n - h]
        ok[: n - h] = v[: n - h] & v[h:]
        ok[:WARMUP] = False
        return y, ok

    # ---- deep inputs ----------------------------------------------------------
    def lob(self) -> np.ndarray:
        """(rows, 40) float32 in DeepLOB order: per level ask px, ask sz, bid px, bid sz."""
        r = self.raw
        mid, tick = self.mid(), self.tick()
        out = np.zeros((len(r), 4 * DEPTH), dtype=np.float32)
        last_a = np.zeros(len(r))
        last_b = np.zeros(len(r))
        for m in range(DEPTH):
            ap = r[:, ASK_PX + m].astype(np.float64)
            bp = r[:, BID_PX + m].astype(np.float64)
            ra = np.where(ap > 0, (ap - mid) / tick, last_a + 1)  # missing level: one tick beyond the last
            rb = np.where(bp > 0, (bp - mid) / tick, last_b - 1)
            out[:, 4 * m + 0] = ra
            out[:, 4 * m + 1] = np.log1p(r[:, ASK_SZ + m] / 100.0)
            out[:, 4 * m + 2] = rb
            out[:, 4 * m + 3] = np.log1p(r[:, BID_SZ + m] / 100.0)
            last_a, last_b = ra, rb
        return out

    def ofi(self, scale: float) -> np.ndarray:
        """(rows, 12) float32: OFI levels 1-10, signed displayed volume, log events."""
        r = self.raw
        out = np.zeros((len(r), DEPTH + 2), dtype=np.float32)
        out[:, :DEPTH] = r[:, OFI : OFI + DEPTH] / scale
        out[:, DEPTH] = (r[:, BUY_VOL] - r[:, SELL_VOL]) / scale
        out[:, DEPTH + 1] = np.log1p(r[:, EVENTS])
        return out

    # ---- baseline features ----------------------------------------------------
    def linear(self, scale: float) -> np.ndarray:
        """(rows, 48) float64: OFI levels 1-10 and signed volume summed over each
        look-back window, then L1 and L1-5 imbalance, spread in ticks, log events."""
        r = self.raw
        base = np.zeros((len(r), DEPTH + 1))
        base[:, :DEPTH] = r[:, OFI : OFI + DEPTH] / scale
        base[:, DEPTH] = (r[:, BUY_VOL] - r[:, SELL_VOL]) / scale
        c = np.vstack([np.zeros((1, base.shape[1])), np.cumsum(base, axis=0)])  # c[i] = sum of base[:i]
        idx = np.arange(len(r))
        cols = []
        for w in WINDOWS:
            s = c[idx + 1] - c[np.maximum(0, idx + 1 - w)]  # sum over the w intervals ending at t
            cols.append(s[:, :DEPTH])
            cols.append(s[:, DEPTH : DEPTH + 1])
        b1, a1 = r[:, BID_SZ].astype(np.float64), r[:, ASK_SZ].astype(np.float64)
        b5 = r[:, BID_SZ : BID_SZ + 5].sum(axis=1).astype(np.float64)
        a5 = r[:, ASK_SZ : ASK_SZ + 5].sum(axis=1).astype(np.float64)
        imb1 = (b1 - a1) / np.maximum(b1 + a1, 1)
        imb5 = (b5 - a5) / np.maximum(b5 + a5, 1)
        spread = (self.ask - self.bid) / self.tick()
        extra = np.column_stack([imb1, imb5, spread, np.log1p(r[:, EVENTS])])
        return np.hstack(cols + [extra])


LINEAR_NAMES = [n for w in WINDOWS for n in [f"ofi{m + 1}_w{w}" for m in range(DEPTH)] + [f"vol_w{w}"]] + [
    "imb1",
    "imb5",
    "spread",
    "log_events",
]


@dataclass
class Study:
    """The day split and universe of an M5 run."""

    data: Path
    universe: list[str]
    train: list[str]
    val: list[str]
    test: list[str]

    @classmethod
    def from_json(cls, path: Path, data: Path) -> Study:
        j = json.loads(path.read_text())
        return cls(data, j["universe"], j["train"], j["val"], j["test"])

    def days(self) -> list[str]:
        return self.train + self.val + self.test

    def load(self, day: str, symbol: str) -> SymbolDay:
        return SymbolDay.load(self.data, day, symbol)

    def scales(self) -> dict[str, float]:
        """Per-symbol depth scale from the training days only (no look-ahead)."""
        out = {}
        for s in self.universe:
            out[s] = float(np.median([self.load(d, s).depth_scale() for d in self.train]))
        return out
