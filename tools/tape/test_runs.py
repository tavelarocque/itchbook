"""Tests for runs.py.

1. One unit test per filter rule and per run field.
2. A slow reference written as plain loops from the rule text in runs.py's
   docstring (no numpy grouping), compared with build_runs and sign_hidden on
   random tapes with heavy timestamp ties, interleaved symbols and every
   record type.

    uv run pytest tools/tape
"""

from __future__ import annotations

import random

import numpy as np
import pytest

import runs as R
import tape

LO, HI = R.SESSION_LO, R.SESSION_HI
MID = (LO + HI) // 2


def mk(rows: list[dict]) -> np.ndarray:
    """Tape array from dicts; defaults give an eligible E fill at 12:45."""
    a = np.zeros(len(rows), dtype=tape.numpy_dtype())
    for i, r in enumerate(rows):
        d = dict(ts=MID, match=10_000 + i, ref=1, seq=i, shares=100, price=10_0000, bid_px=9_9900, ask_px=10_0100,
                 bid_sz=500, ask_sz=500, locate=1, type=b"E", printable=b"-", side=b"S", sign=1, tstate=b"T",
                 phase=tape.REGULAR, at_head=1, flags=0)
        d.update(r)
        for k, v in d.items():
            a[i][k] = v
    return a


def test_single_fill_is_one_run():
    r = R.build_runs(mk([{}]))
    assert len(r) == 1
    assert (r[0]["sign"], r[0]["q"], r[0]["n_fills"], r[0]["mid2"], r[0]["mid_valid"]) == (1, 100, 1, 20_0000, 1)


def test_sweep_aggregates_fills_and_keeps_first_quote():
    t = mk([
        {"shares": 100, "price": 10_0100, "ask_px": 10_0100, "ask_sz": 100},
        {"shares": 50, "price": 10_0200, "ask_px": 10_0200, "ask_sz": 300},
    ])
    r = R.build_runs(t)
    assert len(r) == 1
    assert (r[0]["q"], r[0]["n_fills"], r[0]["ask_px"], r[0]["ask_sz"]) == (150, 2, 10_0100, 100)
    assert (r[0]["first_px"], r[0]["last_px"]) == (10_0100, 10_0200)


def test_interleaved_symbols_split_but_do_not_break_runs():
    t = mk([{"locate": 1}, {"locate": 2, "sign": -1, "side": b"B"}, {"locate": 1}, {"locate": 1, "ts": MID + 1}])
    r = R.build_runs(t)
    assert list(zip(r["locate"], r["n_fills"], r["ts"])) == [(1, 2, MID), (1, 1, MID + 1), (2, 1, MID)]


@pytest.mark.parametrize(
    "bad",
    [
        {"phase": tape.PRE},  # F1
        {"phase": tape.POST},
        {"tstate": b"H"},  # F2
        {"tstate": b"P"},
        {"flags": tape.OPERATIONAL_HALT},  # F3
        {"sign": 0, "side": b"?", "flags": tape.ORPHAN},  # F7
        {"ts": LO - 1},  # F8
        {"ts": HI},
        {"type": b"C", "printable": b"N"},  # not displayed
        {"type": b"P", "sign": 0},
    ],
)
def test_single_rules_exclude(bad):
    assert len(R.build_runs(mk([bad]))) == 0


def test_session_bounds_are_half_open():
    assert len(R.build_runs(mk([{"ts": LO}]))) == 1
    assert len(R.build_runs(mk([{"ts": HI - 1}]))) == 1


def test_printable_c_counts():
    r = R.build_runs(mk([{"type": b"C", "printable": b"Y", "price": 10_0050}]))
    assert len(r) == 1 and r[0]["first_px"] == 10_0050


def test_cross_timestamp_excludes_fills_of_that_symbol_only():  # F4
    t = mk([
        {"type": b"Q", "sign": 0, "side": b"H", "locate": 1},
        {"locate": 1},
        {"locate": 2},
        {"locate": 1, "ts": MID + 5},
    ])
    r = R.build_runs(t)
    assert list(zip(r["locate"], r["ts"])) == [(1, MID + 5), (2, MID)]


def test_broken_match_excluded_wherever_the_B_appears():  # F5
    t = mk([{"match": 7}, {"match": 8, "ts": MID + 1}, {"type": b"B", "match": 7, "sign": 0, "phase": tape.POST}])
    r = R.build_runs(t)
    assert list(r["ts"]) == [MID + 1]


def test_e_paired_with_nonprintable_c_excluded():  # F6
    t = mk([
        {"match": 5},
        {"type": b"C", "printable": b"N", "match": 5, "ts": MID + 9},
        {"match": 6, "ts": MID + 1},
    ])
    r = R.build_runs(t)
    assert list(r["ts"]) == [MID + 1]


@pytest.mark.parametrize(
    ("bid", "ask", "valid", "mid2"),
    [(9_9900, 10_0100, 1, 20_0000), (10_0000, 10_0100, 1, 20_0100), (10_0000, 10_0000, 0, 0),
     (10_0100, 10_0000, 0, 0), (0, 10_0000, 0, 0), (10_0000, 0, 0, 0)],
)
def test_mid_validity_and_half_tick(bid, ask, valid, mid2):
    r = R.build_runs(mk([{"bid_px": bid, "ask_px": ask}]))
    assert (r[0]["mid_valid"], r[0]["mid2"]) == (valid, mid2)


def test_mixed_sign_run_is_flagged_not_dropped():
    r = R.build_runs(mk([{}, {"sign": -1, "side": b"B"}]))
    assert len(r) == 1 and r[0]["mixed"] == 1


def test_timestamp_regression_raises():
    with pytest.raises(ValueError):
        R.build_runs(mk([{"ts": MID + 5}, {"ts": MID}]))


def test_hidden_attach_and_sign_rule():
    # run at MID, seller-initiated, pre-run quote 9.99 / 10.01 (mid 10.00)
    t = mk([
        {"type": b"P", "sign": 0, "side": b"B", "price": 10_0000, "bid_px": 1, "ask_px": 2},  # before the run, at
        # the pre-run mid: run sign (its own quote is ignored for a same-timestamp P)
        {"locate": 1, "sign": -1, "side": b"B"},
        {"type": b"P", "sign": 0, "side": b"B", "price": 9_9900, "bid_px": 9_9800},  # at pre-run bid: -1
        {"type": b"P", "sign": 0, "side": b"B", "price": 10_0100},  # at pre-run ask, opposite the run: +1
        {"type": b"P", "sign": 0, "side": b"B", "ts": MID + 1, "price": 10_0050},  # no run: above mid, +1
        {"type": b"P", "sign": 0, "side": b"B", "ts": MID + 2, "price": 9_9950},  # below mid: -1
        {"type": b"P", "sign": 0, "side": b"B", "ts": MID + 3, "price": 10_0000},  # at mid, no run: 0
        {"type": b"P", "sign": 0, "side": b"B", "ts": MID + 4, "bid_px": 10_0000, "ask_px": 10_0000},  # locked: 0
        {"type": b"P", "sign": 0, "side": b"B", "locate": 2, "price": 9_9900},  # other symbol, no run: -1
    ])
    r = R.build_runs(t)
    assert (r[0]["hid_n"], r[0]["hid_q"]) == (3, 300)
    assert list(R.sign_hidden(t, r)) == [-1, -1, 1, 1, -1, 0, 0, -1]


# ---- slow reference ------------------------------------------------------------------


def slow_runs(t, lo=LO, hi=HI):
    rows = [dict(zip(tape.FIELDS, (x.item() if hasattr(x, "item") else x for x in rec))) for rec in t]
    crosses = {(r["locate"], r["ts"]) for r in rows if r["type"] == b"Q"}
    broken = {r["match"] for r in rows if r["type"] == b"B"}
    cn = {r["match"] for r in rows if r["type"] == b"C" and r["printable"] == b"N"}
    open_run: dict[int, dict] = {}
    done = []
    for r in rows:
        displayed = r["type"] == b"E" or (r["type"] == b"C" and r["printable"] == b"Y")
        if not displayed:
            continue
        if r["phase"] != tape.REGULAR or r["tstate"] != b"T" or r["flags"] & tape.OPERATIONAL_HALT:
            continue
        if (r["locate"], r["ts"]) in crosses or r["match"] in broken or r["sign"] == 0:
            continue
        if r["type"] == b"E" and r["match"] in cn:
            continue
        if not (lo <= r["ts"] < hi):
            continue
        cur = open_run.get(r["locate"])
        if cur is not None and cur["ts"] == r["ts"]:
            cur["q"] += r["shares"]
            cur["n_fills"] += 1
            cur["last_px"] = r["price"]
            cur["signs"].add(r["sign"])
        else:
            if cur is not None:
                done.append(cur)
            open_run[r["locate"]] = dict(locate=r["locate"], ts=r["ts"], seq=r["seq"], q=r["shares"], n_fills=1,
                                         sign=r["sign"], signs={r["sign"]}, bid_px=r["bid_px"], ask_px=r["ask_px"],
                                         bid_sz=r["bid_sz"], ask_sz=r["ask_sz"], first_px=r["price"],
                                         last_px=r["price"])
    done.extend(open_run.values())
    for d in done:
        d["mixed"] = int(len(d["signs"]) > 1)
        ok = 0 < d["bid_px"] < d["ask_px"] and d["ask_px"] > 0
        d["mid_valid"] = int(ok)
        d["mid2"] = d["bid_px"] + d["ask_px"] if ok else 0
        d["hid_n"] = sum(1 for r in rows if r["type"] == b"P" and (r["locate"], r["ts"]) == (d["locate"], d["ts"]))
        d["hid_q"] = sum(r["shares"] for r in rows
                         if r["type"] == b"P" and (r["locate"], r["ts"]) == (d["locate"], d["ts"]))
    done.sort(key=lambda d: (d["locate"], d["ts"]))
    return done


def slow_sign_hidden(t, runs_list):
    by_key = {(d["locate"], d["ts"]): d for d in runs_list}
    out = []
    for rec in t:
        if rec["type"] != b"P":
            continue
        k = (int(rec["locate"]), int(rec["ts"]))
        p = int(rec["price"])
        if k in by_key:
            run = by_key[k]
            b, a = run["bid_px"], run["ask_px"]
        else:
            run = None
            b, a = int(rec["bid_px"]), int(rec["ask_px"])
        side = ((2 * p > b + a) - (2 * p < b + a)) if 0 < b < a else 0
        if side == 0 and run is not None:
            side = run["sign"]
        out.append(side)
    return out


def random_tape(rng: random.Random, n: int) -> np.ndarray:
    rows, ts, phase, match = [], LO - 3_000, tape.PRE, 1
    states = {loc: "T" for loc in range(1, 5)}
    for i in range(n):
        if rng.random() < 0.45:
            ts += rng.randint(1, 900)
        if i == n // 10:
            phase = tape.REGULAR
        if i == (9 * n) // 10:
            phase = tape.POST
        loc = rng.randint(1, 4)
        if rng.random() < 0.03:
            states[loc] = rng.choice("THPQ")
        match += 1
        k = rng.random()
        bid = rng.choice([0, 9_9900, 10_0000, 10_0000, 10_0000])
        ask = rng.choice([0, 10_0000, 10_0100, 10_0100, 10_0200])
        base = dict(ts=ts, match=match, seq=i, locate=loc, shares=rng.randint(1, 900), price=rng.randint(99_800, 100_300),
                    bid_px=bid, ask_px=ask, tstate=states[loc].encode(), phase=phase,
                    flags=tape.OPERATIONAL_HALT if rng.random() < 0.02 else 0)
        if k < 0.6:
            side = rng.choice("BS")
            orphan = rng.random() < 0.01
            typ, pr = (b"E", b"-") if rng.random() < 0.8 else (b"C", rng.choice([b"Y", b"Y", b"N"]))
            base.update(type=typ, printable=pr, side=b"?" if orphan else side.encode(),
                        sign=0 if orphan else (-1 if side == "B" else 1))
            if orphan:
                base["flags"] |= tape.ORPHAN
            if rng.random() < 0.02:
                base["match"] = rng.randint(1, match)  # reuse: pairs with C(N) or B
        elif k < 0.85:
            base.update(type=b"P", side=b"B", sign=0)
        elif k < 0.93:
            base.update(type=b"Q", side=rng.choice([b"O", b"C", b"H"]), sign=0)
        else:
            base.update(type=b"B", match=rng.randint(1, match), sign=0)
        rows.append(base)
    return mk(rows)


@pytest.mark.parametrize("seed", range(40))
def test_matches_slow_reference(seed):
    rng = random.Random(seed)
    t = random_tape(rng, rng.randint(1, 1500))
    lo, hi = (LO, HI) if seed % 2 else (R.FULL_LO, R.FULL_HI)
    got = R.build_runs(t, session=(lo, hi))
    want = slow_runs(t, lo, hi)
    assert len(got) == len(want)
    for g, w in zip(got, want):
        for k in ("locate", "ts", "seq", "q", "n_fills", "sign", "mixed", "bid_px", "ask_px", "bid_sz", "ask_sz",
                  "first_px", "last_px", "mid_valid", "mid2", "hid_n", "hid_q"):
            assert int(g[k]) == int(w[k]), (seed, k, g, w)
    assert list(R.sign_hidden(t, got)) == slow_sign_hidden(t, want)
