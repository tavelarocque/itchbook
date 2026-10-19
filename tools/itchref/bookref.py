#!/usr/bin/env python3
"""Independent reference L3 book: reproduces `itch-book`'s stdout exactly.

Shares nothing with the C++ engine. Orders live in a dict keyed by reference
number; each price level is an insertion-ordered dict of the refs queued
there, so dict order *is* time priority (a new key goes to the tail, deleting
a key keeps the rest in order). Event semantics are written out in
include/itchbook/book/reference_book.hpp and repeated here independently.

    python3 bookref.py [--every 30] [--symbols AAPL,MSFT] [--depth 5] [--max-messages N] FILE

Slow (a few microseconds per message), so a full day takes a while.
"""

from __future__ import annotations

import argparse
import struct
import sys

import itchref

MASK = (1 << 64) - 1
DEFAULT_SYMBOLS = "AAPL,MSFT,SPY,QQQ,TSLA,AMD,ZVZZT"

# Field extraction formats, from byte 0 of the message (x = skipped byte).
ADD = struct.Struct(">xH8xQcI8xI")  # locate, ref, side, shares, price (A and F)
REF_SHARES = struct.Struct(">xH8xQI")  # locate, ref, shares (E, C, X)
REF = struct.Struct(">xH8xQ")  # locate, ref (D)
REPLACE = struct.Struct(">xH8xQQII")  # locate, old, new, shares, price (U)
DIRECTORY = struct.Struct(">xH8x8s")  # locate, stock (R)

STAT_NAMES = [
    "adds", "executes", "cancels", "deletes", "replaces",
    "orphan_execute", "orphan_cancel", "orphan_delete", "orphan_replace",
    "duplicate_add", "duplicate_replace", "bad_side", "over_execute", "over_cancel",
    "locate_mismatch", "replace_to_zero", "executed_at_head", "executed_behind_head",
    "cancel_to_zero", "max_live_orders",
]  # fmt: skip


def splitmix64(x: int) -> int:
    x = (x + 0x9E3779B97F4A7C15) & MASK
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK
    return x ^ (x >> 31)


def mix(*vs: int) -> int:
    h = 0
    for v in vs:
        h = splitmix64(h ^ v)
    return h


def hms(ns: int) -> str:
    s = ns // 1_000_000_000
    return f"{s // 3600:02d}:{s // 60 % 60:02d}:{s % 60:02d}"


class Book:
    def __init__(self) -> None:
        self.orders: dict[int, list] = {}  # ref -> [locate, side, price, shares]
        self.levels: dict[tuple, list] = {}  # (locate, side, price) -> [total_shares, {ref: None}]
        self.st = dict.fromkeys(STAT_NAMES, 0)

    def _insert(self, ref, locate, side, shares, price):
        key = (locate, side, price)
        lvl = self.levels.get(key)
        if lvl is None:
            lvl = self.levels[key] = [0, {}]
        lvl[0] += shares
        lvl[1][ref] = None
        self.orders[ref] = [locate, side, price, shares]
        if len(self.orders) > self.st["max_live_orders"]:
            self.st["max_live_orders"] = len(self.orders)

    def _erase(self, ref, o):
        key = (o[0], o[1], o[2])
        lvl = self.levels[key]
        lvl[0] -= o[3]
        del lvl[1][ref]
        if not lvl[1]:
            del self.levels[key]
        del self.orders[ref]

    def _reduce(self, ref, o, n):
        if n >= o[3]:
            self._erase(ref, o)
        else:
            self.levels[(o[0], o[1], o[2])][0] -= n
            o[3] -= n

    def add(self, locate, ref, side, shares, price):
        if side not in (b"B", b"S"):
            self.st["bad_side"] += 1
            return
        if ref in self.orders:
            self.st["duplicate_add"] += 1
            return
        self.st["adds"] += 1
        self._insert(ref, locate, 0 if side == b"B" else 1, shares, price)

    def execute(self, locate, ref, n):
        o = self.orders.get(ref)
        if o is None:
            self.st["orphan_execute"] += 1
            return
        self.st["executes"] += 1
        if locate != o[0]:
            self.st["locate_mismatch"] += 1
        q = self.levels[(o[0], o[1], o[2])][1]
        if next(iter(q)) == ref:
            self.st["executed_at_head"] += 1
        else:
            self.st["executed_behind_head"] += 1
        if n > o[3]:
            self.st["over_execute"] += 1
        self._reduce(ref, o, n)

    def cancel(self, locate, ref, n):
        o = self.orders.get(ref)
        if o is None:
            self.st["orphan_cancel"] += 1
            return
        self.st["cancels"] += 1
        if locate != o[0]:
            self.st["locate_mismatch"] += 1
        if n > o[3]:
            self.st["over_cancel"] += 1
        if n >= o[3]:
            self.st["cancel_to_zero"] += 1
        self._reduce(ref, o, n)

    def delete(self, locate, ref):
        o = self.orders.get(ref)
        if o is None:
            self.st["orphan_delete"] += 1
            return
        self.st["deletes"] += 1
        if locate != o[0]:
            self.st["locate_mismatch"] += 1
        self._erase(ref, o)

    def replace(self, locate, old, new, shares, price):
        o = self.orders.get(old)
        if o is None:
            self.st["orphan_replace"] += 1
            return
        if new != old and new in self.orders:
            self.st["duplicate_replace"] += 1
            return
        self.st["replaces"] += 1
        if locate != o[0]:
            self.st["locate_mismatch"] += 1
        loc, side = o[0], o[1]
        self._erase(old, o)
        if shares == 0:
            self.st["replace_to_zero"] += 1
            return
        self._insert(new, loc, side, shares, price)

    def checksum(self) -> int:
        total = 0
        for (loc, side, price), (shares, q) in self.levels.items():
            total += mix(loc, side, price, shares, len(q))
            for rank, ref in enumerate(q):
                total += mix(ref, loc, side, price, self.orders[ref][3], rank)
        return total & MASK

    def top(self, locate, side, depth):
        prices = sorted((p for (lc, s, p) in self.levels if lc == locate and s == side), reverse=(side == 0))
        return [(p, *self._level_info(locate, side, p)) for p in prices[:depth]]

    def _level_info(self, locate, side, price):
        shares, q = self.levels[(locate, side, price)]
        return shares, len(q)


def run(path: str, every_min: int, symbols: list[str], depth: int, max_messages: int | None) -> str:
    book = Book()
    minute = 60_000_000_000
    times = list(range(210 * minute, 1230 * minute + 1, every_min * minute))
    want = set(symbols)
    names: dict[int, str] = {}
    top_locates: list[int] = []
    out_snaps: list[str] = []
    next_snap = 0
    messages = 0

    def snapshot(at):
        out_snaps.append(f"snap\t{hms(at)}\t{book.checksum()}\t{len(book.orders)}\t{len(book.levels)}")
        rows = []
        for loc in top_locates:
            for side in (0, 1):
                for rank, (p, sh, n) in enumerate(book.top(loc, side, depth)):
                    rows.append(
                        (names[loc], side, rank, f"top\t{hms(at)}\t{names[loc]}\t{'BS'[side]}\t{rank}\t{p}\t{sh}\t{n}")
                    )
        out_snaps.extend(r[3] for r in sorted(rows))

    lengths = itchref.LENGTHS
    f = itchref.open_any(path)
    buf = b""
    stopped = False
    while not stopped:
        data = f.read(1 << 22)
        if not data:
            break
        buf += data
        pos, n = 0, len(buf)
        while pos + 2 <= n:
            ln = (buf[pos] << 8) | buf[pos + 1]
            if pos + 2 + ln > n:
                break
            if ln and max_messages is not None and messages == max_messages:
                stopped = True
                break
            s = pos + 2
            pos = s + ln
            if ln == 0:
                continue
            t = buf[s]
            if lengths.get(chr(t)) != ln:
                continue
            messages += 1
            if next_snap < len(times):
                ts = int.from_bytes(buf[s + 5 : s + 11], "big")
                while next_snap < len(times) and ts > times[next_snap]:
                    snapshot(times[next_snap])
                    next_snap += 1
            if t == 65 or t == 70:  # A, F
                book.add(*ADD.unpack_from(buf, s))
            elif t == 68:  # D
                book.delete(*REF.unpack_from(buf, s))
            elif t == 85:  # U
                book.replace(*REPLACE.unpack_from(buf, s))
            elif t == 69 or t == 67:  # E, C
                book.execute(*REF_SHARES.unpack_from(buf, s))
            elif t == 88:  # X
                book.cancel(*REF_SHARES.unpack_from(buf, s))
            elif t == 82:  # R
                loc, stock = DIRECTORY.unpack_from(buf, s)
                name = stock.decode("ascii").rstrip(" ")
                names[loc] = name
                if name in want:
                    top_locates.append(loc)
        buf = buf[pos:]
    while next_snap < len(times):
        snapshot(times[next_snap])
        next_snap += 1

    lines = [f"messages\t{messages}"]
    lines += [f"stat\t{k}\t{book.st[k]}" for k in STAT_NAMES]
    lines += [f"end_live_orders\t{len(book.orders)}", f"end_levels\t{len(book.levels)}"]
    lines += out_snaps
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("file")
    ap.add_argument("--every", type=int, default=30)
    ap.add_argument("--symbols", default=DEFAULT_SYMBOLS)
    ap.add_argument("--depth", type=int, default=5)
    ap.add_argument("--max-messages", type=int)
    a = ap.parse_args()
    sys.stdout.write(run(a.file, a.every, a.symbols.split(","), a.depth, a.max_messages))
    return 0


if __name__ == "__main__":
    sys.exit(main())
