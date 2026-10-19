#!/usr/bin/env python3
"""Compare the C++ L3 book (itch-book) with the Python reference (bookref).

Synthetic mode generates a realistic order-event stream (stock directory,
adds, executions, cancels, deletes and replaces across a whole trading day,
plus deliberate anomalies) and requires byte-identical output:

    crosscheck_book.py --itch-book build/release/itch-book --synthetic 50000

File mode runs both over a real day and diffs them:

    crosscheck_book.py --itch-book build/release/itch-book --file 01302020.NASDAQ_ITCH50.gz
"""

from __future__ import annotations

import argparse
import difflib
import gzip
import random
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bookref  # noqa: E402

SYMBOLS = ["AAPL", "MSFT", "SPY", "ZVZZT", "TINY", "ODD"]


def frame(p: bytes) -> bytes:
    return len(p).to_bytes(2, "big") + p


def hdr(t: str, locate: int, ts: int) -> bytes:
    return t.encode() + locate.to_bytes(2, "big") + b"\x00\x00" + ts.to_bytes(6, "big")


def alpha(s: str, n: int) -> bytes:
    return s.encode().ljust(n, b" ")


def synthetic_stream(n: int, seed: int) -> bytes:
    rng = random.Random(seed)
    out = bytearray()
    ts = 3 * 3600 * 10**9
    step = (17 * 3600 * 10**9) // max(n, 1)  # spread events from 03:00 to 20:00
    out += frame(hdr("S", 0, ts) + b"O")
    for i, s in enumerate(SYMBOLS, start=1):
        body = alpha(s, 8) + b"QN" + (100).to_bytes(4, "big") + b"NCZ PNN1N" + bytes(4) + b"N"
        out += frame(hdr("R", i, ts) + body)
    live: dict[int, int] = {}  # ref -> locate
    issued: list[int] = []
    next_ref = 1
    for _ in range(n):
        ts += rng.randint(1, 2 * step)
        k = rng.random()
        loc = rng.randint(1, len(SYMBOLS))
        price = rng.randint(1000, 1012) * 100  # 13 ticks: plenty of shared levels
        shares = rng.randint(1, 500)
        target = rng.choice(issued) if issued and rng.random() < 0.97 else next_ref + 10**6
        tloc = live.get(target, loc)
        if k < 0.40 or not live:
            side = rng.choice(b"BS") if rng.random() > 0.001 else ord("Z")
            t = "A" if rng.random() < 0.9 else "F"
            body = next_ref.to_bytes(8, "big") + bytes([side]) + shares.to_bytes(4, "big")
            body += alpha(SYMBOLS[loc - 1], 8) + price.to_bytes(4, "big")
            if t == "F":
                body += b"NSDQ"
            out += frame(hdr(t, loc, ts) + body)
            if side != ord("Z"):
                live[next_ref] = loc
            issued.append(next_ref)
            next_ref += 1
        elif k < 0.55:
            if rng.random() < 0.5:
                out += frame(hdr("E", tloc, ts) + target.to_bytes(8, "big") + shares.to_bytes(4, "big") + bytes(8))
            else:
                body = (
                    target.to_bytes(8, "big") + shares.to_bytes(4, "big") + bytes(8) + b"Y" + price.to_bytes(4, "big")
                )
                out += frame(hdr("C", tloc, ts) + body)
        elif k < 0.65:
            out += frame(hdr("X", tloc, ts) + target.to_bytes(8, "big") + shares.to_bytes(4, "big"))
        elif k < 0.85:
            out += frame(hdr("D", tloc, ts) + target.to_bytes(8, "big"))
            live.pop(target, None)
        elif k < 0.98:
            new = next_ref if rng.random() > 0.01 else rng.choice(issued)
            body = target.to_bytes(8, "big") + new.to_bytes(8, "big") + rng.randint(0, 500).to_bytes(4, "big")
            out += frame(hdr("U", tloc, ts) + body + price.to_bytes(4, "big"))
            issued.append(next_ref)
            next_ref += 1
        else:  # non-book traffic: hidden trade, NOII
            body = bytes(8) + b"B" + shares.to_bytes(4, "big") + alpha(SYMBOLS[loc - 1], 8)
            out += frame(hdr("P", loc, ts) + body + price.to_bytes(4, "big") + bytes(8))
    out += frame(hdr("S", 0, ts + 1) + b"C")
    return bytes(out)


def run_cpp(binary: str, path: str, extra: list[str]) -> str:
    p = subprocess.run([binary, *extra, path], capture_output=True, text=True, errors="backslashreplace")
    sys.stderr.write(p.stderr)
    if p.returncode not in (0, 2):
        raise SystemExit(f"itch-book failed with exit code {p.returncode}")
    return p.stdout


def same(a: str, b: str) -> bool:
    if a == b:
        return True
    sys.stdout.writelines(difflib.unified_diff(a.splitlines(True), b.splitlines(True), "itch-book", "bookref", n=1))
    return False


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--itch-book", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--synthetic", type=int, metavar="N")
    g.add_argument("--file")
    ap.add_argument("--seed", type=int, default=20261019)
    ap.add_argument("--max-messages", type=int)
    a = ap.parse_args()
    syms = ["--symbols", ",".join(SYMBOLS)] if a.synthetic else []
    lim = ["--max-messages", str(a.max_messages)] if a.max_messages else []

    if a.synthetic:
        raw = synthetic_stream(a.synthetic, a.seed)
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "book.itch.gz"
            p.write_bytes(gzip.compress(raw))
            cpp = run_cpp(a.itch_book, str(p), syms)
            py = bookref.run(str(p), 30, SYMBOLS, 5, None)
    else:
        cpp = run_cpp(a.itch_book, a.file, lim)
        py = bookref.run(a.file, 30, bookref.DEFAULT_SYMBOLS.split(","), 5, a.max_messages)
    ok = same(cpp, py)
    n_snap = sum(1 for line in py.splitlines() if line.startswith("snap\t"))
    n_top = sum(1 for line in py.splitlines() if line.startswith("top\t"))
    print(f"book cross-check: {n_snap} snapshots, {n_top} top-of-book rows: {'MATCH' if ok else 'MISMATCH'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
