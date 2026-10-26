#!/usr/bin/env python3
"""Independent Python reference for the itch-trades tape.

Decodes with itchref's struct layouts and replays bookref.Book (neither shares
code with the C++ side), and re-implements every field rule of
include/itchbook/tape/record.hpp from its comments: quote read before the
message, sign opposite the resting side, C priced at its execution price, E at
the resting order's price, flags, phase, trading state, NASDAQ 'h' halts.
Yields tuples in tape.FIELDS order.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "tape"))

import bookref  # noqa: E402
import itchref  # noqa: E402
import tape  # noqa: E402

U32 = (1 << 32) - 1


class IndexedBook(bookref.Book):
    """bookref.Book plus a price set per (locate, side) for a fast best quote."""

    def __init__(self) -> None:
        super().__init__()
        self.prices: dict[tuple, set] = defaultdict(set)

    def _insert(self, ref, locate, side, shares, price):
        new_level = (locate, side, price) not in self.levels
        super()._insert(ref, locate, side, shares, price)
        if new_level:
            self.prices[(locate, side)].add(price)

    def _erase(self, ref, o):
        key = (o[0], o[1], o[2])
        super()._erase(ref, o)
        if key not in self.levels:
            self.prices[(o[0], o[1])].discard(o[2])

    def quote(self, loc: int) -> tuple[int, int, int, int, int]:
        """bid_px, ask_px, bid_sz, ask_sz, flags (0 = empty side)."""
        flags = 0
        bid_px = ask_px = bid_sz = ask_sz = 0
        bids = self.prices.get((loc, 0))
        if bids:
            bid_px = max(bids)
            bid_sz = self.levels[(loc, 0, bid_px)][0]
            if bid_sz > U32:
                bid_sz, flags = U32, flags | tape.BID_SZ_SATURATED
        asks = self.prices.get((loc, 1))
        if asks:
            ask_px = min(asks)
            ask_sz = self.levels[(loc, 1, ask_px)][0]
            if ask_sz > U32:
                ask_sz, flags = U32, flags | tape.ASK_SZ_SATURATED
        return bid_px, ask_px, bid_sz, ask_sz, flags


def build(path: str, max_messages: int | None = None) -> Iterator[tuple]:
    book = IndexedBook()
    phase = tape.PRE
    tstate: dict[int, str] = {}
    halted: dict[int, bool] = {}
    S, L = itchref.STRUCTS, itchref.LENGTHS
    seq = 0

    def rec(ts, match, ref, shares, price, loc, typ, printable, side, sign, at_head, flags):
        bid_px, ask_px, bid_sz, ask_sz, qflags = book.quote(loc)
        flags |= qflags
        if halted.get(loc):
            flags |= tape.OPERATIONAL_HALT
        return (ts, match, ref, seq, shares, price, bid_px, ask_px, bid_sz, ask_sz, loc, typ, printable, side,
                sign, tstate.get(loc, "?"), phase, at_head, flags)

    f = itchref.open_any(path)
    buf = b""
    while True:
        data = f.read(1 << 22)
        if not data:
            break
        buf += data
        pos, n = 0, len(buf)
        while pos + 2 <= n:
            ln = (buf[pos] << 8) | buf[pos + 1]
            if pos + 2 + ln > n:
                break
            if ln and max_messages is not None and seq == max_messages:
                return
            s = pos + 2
            pos = s + ln
            if ln == 0:
                continue
            t = chr(buf[s])
            if L.get(t) != ln:
                continue
            v = S[t].unpack_from(buf, s + 1)
            loc, ts = v[0], int.from_bytes(v[2], "big")
            if t == "S":
                if v[3] == b"Q":
                    phase = tape.REGULAR
                elif v[3] == b"M":
                    phase = tape.POST
            elif t == "H":
                tstate[loc] = v[4].decode("latin-1")
            elif t == "h":
                if v[4] == b"Q":
                    halted[loc] = v[5] == b"H"
            elif t in "AF":
                book.add(loc, v[3], v[4], v[5], v[7])
            elif t == "X":
                book.cancel(loc, v[3], v[4])
            elif t == "D":
                book.delete(loc, v[3])
            elif t == "U":
                book.replace(loc, v[3], v[4], v[5], v[6])
            elif t in "EC":
                ref, n_sh, match = v[3], v[4], v[5]
                o = book.orders.get(ref)
                if t == "C":
                    printable, xprice = v[6].decode("latin-1"), v[7]
                else:
                    printable, xprice = "-", None
                if o is not None:
                    rloc, side_i, oprice, oshares = o
                    q = book.levels[(rloc, side_i, oprice)][1]
                    at_head = 1 if next(iter(q)) == ref else 0
                    price = xprice if t == "C" else oprice
                    flags = tape.OVER_EXECUTE if n_sh > oshares else 0
                    r = rec(ts, match, ref, n_sh, price, rloc, t, printable, "BS"[side_i],
                            -1 if side_i == 0 else 1, at_head, flags)
                else:
                    r = rec(ts, match, ref, n_sh, xprice if t == "C" else 0, loc, t, printable, "?", 0, 0,
                            tape.ORPHAN)
                book.execute(loc, ref, n_sh)
                yield r
            elif t == "P":
                yield rec(ts, v[8], v[3], v[5], v[7], loc, "P", "-", v[4].decode("latin-1"), 0, 0, 0)
            elif t == "Q":
                sh, flags = v[3], 0
                if sh > U32:
                    sh, flags = U32, tape.SHARES_SATURATED
                yield rec(ts, v[6], 0, sh, v[5], loc, "Q", "-", v[7].decode("latin-1"), 0, 0, flags)
            elif t == "B":
                yield rec(ts, v[3], 0, 0, 0, loc, "B", "-", "-", 0, 0, 0)
            seq += 1
        buf = buf[pos:]
