#!/usr/bin/env python3
"""Runs (reconstructed market orders) and the pre-run NASDAQ mid, from a tape.

A run is a maximal set of eligible displayed executions on one symbol sharing
one nanosecond timestamp. Timestamps never decrease in file order, so on one
symbol these fills are contiguous in that symbol's sequence of eligible fills.
A displayed fill (E, or C with printable 'Y') is eligible when ALL hold:

  F1  session phase is regular (between the S 'Q' and 'M' events)
  F2  the symbol's trading state is 'T'
  F3  no NASDAQ operational halt ('h') is in effect for the symbol
  F4  no Q cross print exists at the same (locate, timestamp)
  F5  its match number is not cancelled by any B (broken trade) message
  F6  it is not an E sharing a match number with a non-printable C
  F7  its resting order was known (not an orphan), so it has a sign
  F8  its timestamp lies in the trimmed session [session_lo, session_hi)

Each run gets the sign of its fills (the opposite of the resting side), its
volume (sum of fill shares), the quote of its first fill, read before that
fill, and the mid twice-encoded as mid2 = bid + ask (exact at half ticks),
valid only when 0 < bid < ask. Runs whose fills disagree on sign are kept but
flagged (none were seen on eight days).

Hidden P executions sharing a run's (locate, ts) are counted on the run; their
sensitivity-arm sign is set by sign_hidden.

    python3 runs.py TAPE.tape.gz OUT.npz      # build and save one day
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tape  # noqa: E402

NS_PER_MIN = 60 * 10**9
SESSION_LO = 600 * NS_PER_MIN  # 10:00:00 ET
SESSION_HI = 930 * NS_PER_MIN  # 15:30:00 ET
FULL_LO = 0
FULL_HI = 24 * 60 * NS_PER_MIN

RUN_DTYPE = np.dtype(
    [
        ("locate", "<u2"),
        ("sign", "i1"),
        ("mid_valid", "u1"),
        ("mixed", "u1"),
        ("n_fills", "<u4"),
        ("ts", "<u8"),
        ("seq", "<u4"),  # seq of the first fill: global file order
        ("q", "<u8"),  # run volume in shares
        ("bid_px", "<u4"),
        ("ask_px", "<u4"),
        ("bid_sz", "<u4"),
        ("ask_sz", "<u4"),
        ("mid2", "<u8"),  # bid + ask, 0 when invalid
        ("first_px", "<u4"),
        ("last_px", "<u4"),
        ("hid_n", "<u4"),  # P executions with the same (locate, ts)
        ("hid_q", "<u8"),
    ]
)


def _key(locate: np.ndarray, ts: np.ndarray) -> np.ndarray:
    # ts < 2^47 (a day is 8.64e13 ns), so locate fits above bit 48
    return (locate.astype(np.uint64) << np.uint64(48)) | ts.astype(np.uint64)


def eligible_mask(t: np.ndarray, session: tuple[int, int] = (SESSION_LO, SESSION_HI)) -> dict[str, np.ndarray]:
    """Per-record booleans for every rule, and 'eligible' (all of them, on displayed fills)."""
    typ = t["type"]
    is_e, is_c, is_p, is_q, is_b = (typ == b"E"), (typ == b"C"), (typ == b"P"), (typ == b"Q"), (typ == b"B")
    displayed = is_e | (is_c & (t["printable"] == b"Y"))
    key = _key(t["locate"], t["ts"])
    cross_keys = np.unique(key[is_q])
    broken = np.unique(t["match"][is_b])
    cn_match = np.unique(t["match"][is_c & (t["printable"] == b"N")])
    m = {
        "displayed": displayed,
        "F1_regular": t["phase"] == tape.REGULAR,
        "F2_state_T": t["tstate"] == b"T",
        "F3_no_halt": (t["flags"] & tape.OPERATIONAL_HALT) == 0,
        "F4_not_at_cross": ~np.isin(key, cross_keys),
        "F5_not_broken": ~np.isin(t["match"], broken),
        "F6_not_E_paired_with_CN": ~(is_e & np.isin(t["match"], cn_match)),
        "F7_signed": t["sign"] != 0,
        "F8_in_session": (t["ts"] >= session[0]) & (t["ts"] < session[1]),
        "is_p": is_p,
    }
    ok = displayed.copy()
    for k, v in m.items():
        if k.startswith("F"):
            ok &= v
    m["eligible"] = ok
    return m


def build_runs(t: np.ndarray, session: tuple[int, int] = (SESSION_LO, SESSION_HI), rules=None) -> np.ndarray:
    """Runs of one day, ordered by (locate, ts). `rules` restricts F1..F8 (for cross-checks)."""
    m = eligible_mask(t, session)
    if rules is not None:
        ok = m["displayed"].copy()
        for r in rules:
            ok &= m[r]
    else:
        ok = m["eligible"]
    f = t[ok]
    order = np.argsort(f["locate"], kind="stable")  # file order kept within a symbol
    f = f[order]
    loc, ts = f["locate"], f["ts"]
    same_loc = loc[1:] == loc[:-1]
    if np.any(same_loc & (ts[1:] < ts[:-1])):
        raise ValueError("timestamp regression within a symbol: runs undefined")
    start = np.ones(len(f), dtype=bool)
    if len(f):
        start[1:] = ~same_loc | (ts[1:] != ts[:-1])
    idx = np.flatnonzero(start)
    n = len(idx)
    out = np.zeros(n, dtype=RUN_DTYPE)
    if n == 0:
        return out
    ends = np.append(idx[1:], len(f))
    out["locate"] = loc[idx]
    out["ts"] = ts[idx]
    out["seq"] = f["seq"][idx]
    out["n_fills"] = ends - idx
    out["q"] = np.add.reduceat(f["shares"].astype(np.uint64), idx)
    sgn = f["sign"].astype(np.int64)
    s_min = np.minimum.reduceat(sgn, idx)
    s_max = np.maximum.reduceat(sgn, idx)
    out["sign"] = sgn[idx]
    out["mixed"] = (s_min != s_max).astype(np.uint8)
    for k in ("bid_px", "ask_px", "bid_sz", "ask_sz"):
        out[k] = f[k][idx]
    out["first_px"] = f["price"][idx]
    out["last_px"] = f["price"][ends - 1]
    bid = out["bid_px"].astype(np.uint64)
    ask = out["ask_px"].astype(np.uint64)
    valid = (bid > 0) & (ask > 0) & (bid < ask)
    out["mid_valid"] = valid.astype(np.uint8)
    out["mid2"] = np.where(valid, bid + ask, 0)
    # hidden P at the same (locate, ts) as a run
    p = t[m["is_p"]]
    if len(p):
        pk = _key(p["locate"], p["ts"])
        rk = _key(out["locate"], out["ts"])
        pos = np.searchsorted(rk, pk)  # rk is sorted: locate major, ts minor
        hit = (pos < n) & (rk[np.minimum(pos, n - 1)] == pk)
        np.add.at(out["hid_n"], pos[hit], 1)
        np.add.at(out["hid_q"], pos[hit], p["shares"][hit].astype(np.uint64))
    return out


def sign_hidden(t: np.ndarray, runs: np.ndarray) -> np.ndarray:
    """Sensitivity-arm sign for each P record (in tape order; 0 = unsigned).

    The price side of a displayed mid decides when it is informative: +1 above
    the mid, -1 below. For a P sharing (locate, ts) with a run, the mid is the
    run's pre-run quote and the run's sign breaks a tie at the mid. Any other
    P uses its own pre-trade quote and stays unsigned at the mid or when the
    mid is undefined (locked, crossed, one-sided).

    Measured on 01302020 (P sharing a run's (locate, ts)): prints after the
    run's first fill sit at or beyond the pre-run touch and agree with the run
    100%; prints before it sit inside the spread and agree 97%, 44% of them
    at the mid; 0.7% sit at the touch opposite the run (another incoming
    order in the same nanosecond), which the price side signs correctly."""
    p = t[t["type"] == b"P"]
    out = np.zeros(len(p), dtype=np.int8)
    if len(p) == 0:
        return out
    bid = p["bid_px"].astype(np.int64)
    ask = p["ask_px"].astype(np.int64)
    run_sign = np.zeros(len(p), dtype=np.int8)
    hit = np.zeros(len(p), dtype=bool)
    if len(runs):
        rk = _key(runs["locate"], runs["ts"])
        pk = _key(p["locate"], p["ts"])
        pos = np.searchsorted(rk, pk)
        hit = (pos < len(runs)) & (rk[np.minimum(pos, len(runs) - 1)] == pk)
        hr = runs[pos[hit]]
        bid[hit] = hr["bid_px"]
        ask[hit] = hr["ask_px"]
        run_sign[hit] = hr["sign"]
    valid = (bid > 0) & (ask > 0) & (bid < ask)
    side = np.sign(2 * p["price"].astype(np.int64) - (bid + ask)).astype(np.int8)
    out[valid] = side[valid]
    tie = hit & (out == 0)
    out[tie] = run_sign[tie]
    return out


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    t = tape.read_tape(sys.argv[1])
    m = eligible_mask(t)
    runs = build_runs(t)
    np.savez_compressed(sys.argv[2], runs=runs)
    disp = m["displayed"]
    print(f"records {len(t):,}  displayed fills {disp.sum():,}  eligible {m['eligible'].sum():,}  runs {len(runs):,}")
    for k in [k for k in m if k.startswith("F")]:
        print(f"  {k:26s} removes {int((disp & ~m[k]).sum()):>10,} displayed fills (alone)")
    print(f"  mixed-sign runs {int(runs['mixed'].sum())}, runs without a valid mid {int((runs['mid_valid'] == 0).sum())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
