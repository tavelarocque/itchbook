"""Volume surrogates for one stock-day of runs.

A surrogate keeps every run's sign, timestamp and mid and changes only which
volume q goes with which run. Each is an index idx with q_new = q[idx]:

  S1  a uniform random permutation of q over the stock-day. Destroys size
      autocorrelation, intraday size seasonality and every sign/size and
      size/return link.
  S2  a block permutation within 30-minute clock bins. Bins are ts // 30 min
      with ts in ns since midnight ET, so they start on the hour and the half
      hour. Inside a bin the runs are cut into consecutive blocks of L runs
      starting at the bin's first run (the last block may be shorter), and the
      blocks are put in random order, each keeping its internal order. Keeps
      each bin's q multiset (intraday size seasonality) and the size
      autocorrelation at lags 1..L-1 for pairs inside a block; breaks the
      sign/size link. With an AR(1)-like size series the lag-k ACF survives as
      roughly (1 - k/L) of the original.

choose_L(q, ts): the smallest lag k >= 1 at which the (centred) ACF of log q
falls below 0.05, capped at (median runs per occupied 30-minute bin) // 4 so a
typical bin keeps at least four blocks to permute; never below 1.

surrogate_runs applies either surrogate to a runs array (runs.RUN_DTYPE, one
stock-day): only the 'q' field moves. Fields derived from the fills (n_fills,
first_px, last_px, hid_q) stay with their run, so they no longer match q.
"""
from __future__ import annotations

import numpy as np

BIN_NS = 30 * 60 * 10**9


def clock_bins(ts: np.ndarray) -> np.ndarray:
    """30-minute clock bin of each run (ts in ns since midnight)."""
    return np.asarray(ts, dtype=np.int64) // BIN_NS


def s1_index(n: int, rng: np.random.Generator) -> np.ndarray:
    return rng.permutation(n)


def s2_index(ts: np.ndarray, L: int, rng: np.random.Generator) -> np.ndarray:
    ts = np.asarray(ts, dtype=np.int64)
    if L < 1:
        raise ValueError("block length L must be >= 1")
    if np.any(np.diff(ts) < 0):
        raise ValueError("runs must be in time order")
    n = len(ts)
    bins = clock_bins(ts)
    new_bin = np.ones(n, dtype=bool)
    new_bin[1:] = bins[1:] != bins[:-1]
    first = np.maximum.accumulate(np.where(new_bin, np.arange(n), 0))  # first run of each run's bin
    pos = np.arange(n) - first
    blk = pos // L
    new_blk = np.ones(n, dtype=bool)
    new_blk[1:] = new_bin[1:] | (blk[1:] != blk[:-1])
    g = np.cumsum(new_blk) - 1  # global block id
    u = rng.random(int(g[-1]) + 1 if n else 0)
    # within a bin, blocks in random order, runs in original order inside a block
    return np.lexsort((np.arange(n), u[g], bins))


def s1(q: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    q = np.asarray(q)
    return q[s1_index(len(q), rng)]


def s2(q: np.ndarray, ts: np.ndarray, L: int, rng: np.random.Generator) -> np.ndarray:
    q = np.asarray(q)
    if len(q) != len(ts):
        raise ValueError("q and ts differ in length")
    return q[s2_index(ts, L, rng)]


def acf(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Centred sample autocorrelation r_k, k = 0..max_lag (denominator n at every lag)."""
    x = np.asarray(x, dtype=float)
    x = x - x.mean()
    n = len(x)
    size = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(x, size)
    ac = np.fft.irfft(f * np.conj(f), size)[: max_lag + 1]
    return ac / ac[0]


def choose_L(q: np.ndarray, ts: np.ndarray, threshold: float = 0.05) -> int:
    q = np.asarray(q, dtype=float)
    if np.any(q <= 0):
        raise ValueError("volumes must be positive")
    counts = np.unique(clock_bins(ts), return_counts=True)[1]
    cap = max(1, int(np.median(counts)) // 4)
    r = acf(np.log(q), cap)
    below = np.flatnonzero(r[1:] < threshold)
    return int(below[0] + 1) if len(below) else cap


def surrogate_runs(runs: np.ndarray, kind: str, rng: np.random.Generator, L: int | None = None) -> np.ndarray:
    """Copy of one stock-day's runs with q replaced by surrogate 'S1' or 'S2' (L from choose_L if None)."""
    if len(np.unique(runs["locate"])) > 1:
        raise ValueError("one stock-day at a time")
    out = runs.copy()
    if kind == "S1":
        idx = s1_index(len(runs), rng)
    elif kind == "S2":
        L = choose_L(runs["q"], runs["ts"]) if L is None else L
        idx = s2_index(runs["ts"], L, rng)
    else:
        raise ValueError("kind is 'S1' or 'S2'")
    out["q"] = runs["q"][idx]
    return out
