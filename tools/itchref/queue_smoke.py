#!/usr/bin/env python3
"""End-to-end check of itch-queue on a synthetic day.

Wraps the synthetic stream from crosscheck_book.py in a market open (S 'Q'),
trading-state messages for every symbol and a market close (S 'M'), runs
itch-queue with every eligible order tracked and the truth verified against
the L3 book at every evaluation point, and checks the outputs are consistent.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from crosscheck_book import SYMBOLS, alpha, frame, hdr, synthetic_stream  # noqa: E402


def with_session(raw: bytes) -> bytes:
    """Insert Q and H messages after the directory, and M before the end."""
    out = bytearray()
    i = 0
    inserted = False
    ts_open = 9 * 3600 * 10**9
    while i < len(raw):
        n = int.from_bytes(raw[i : i + 2], "big")
        payload = raw[i + 2 : i + 2 + n]
        if not inserted and payload[:1] not in (b"S", b"R"):
            out += frame(hdr("S", 0, ts_open) + b"Q")
            for loc, s in enumerate(SYMBOLS, start=1):
                out += frame(hdr("H", loc, ts_open) + alpha(s, 8) + b"T" + b" " + b"    ")
            inserted = True
        if payload[:1] == b"S" and payload[11:12] == b"C":
            ts = int.from_bytes(payload[5:11], "big")
            out += frame(hdr("S", 0, ts) + b"M")
        out += raw[i : i + 2 + n]
        i += 2 + n
    return bytes(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--itch-queue", required=True)
    ap.add_argument("--synthetic", type=int, default=30000)
    ap.add_argument("--seed", type=int, default=20261019)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory() as d:
        stream = Path(d) / "day.itch"
        stream.write_bytes(with_session(synthetic_stream(a.synthetic, a.seed)))
        prefix = str(Path(d) / "q")
        res = subprocess.run(
            [a.itch_queue, "--sample", "1", "--events", "1", "--verify", "1", "--out", prefix, str(stream)],
            capture_output=True,
            text=True,
            check=False,
        )
        print(res.stdout, end="")
        if res.returncode != 0:
            print(res.stderr, file=sys.stderr)
            print(f"FAIL: exit {res.returncode}")
            return 1
        stats = dict(line.split("\t", 1) for line in res.stdout.splitlines())
        verified, failures = stats["verified"].split("\tfailures\t")
        with gzip.open(prefix + ".orders.tsv.gz", "rt") as f:
            orders = list(csv.DictReader(f, delimiter="\t"))
        with gzip.open(prefix + ".points.tsv.gz", "rt") as f:
            points = list(csv.DictReader(f, delimiter="\t"))
        problems = []
        if int(verified) == 0 or int(failures) != 0:
            problems.append(f"verified {verified}, failures {failures}")
        if len(orders) != int(stats["tracked_orders"]) or len(orders) < 50:
            problems.append(f"{len(orders)} order rows vs tracked_orders {stats['tracked_orders']}")
        if len(points) != sum(int(o["n_eval"]) for o in orders):
            problems.append("points file does not hold every evaluation point")
        if sum(int(o["n_eval"]) for o in orders) != int(verified):
            problems.append("verify did not cover every evaluation point")
        for o in orders:
            if int(o["filled"]) > int(o["own0"]) or int(o["n_eval"]) > int(o["n_trades"]):
                problems.append(f"inconsistent row {o['ref']}")
                break
            if int(o["first_fill_ts"]) and o["end"] == "close" and int(o["filled"]) == 0:
                problems.append(f"fill time without filled shares {o['ref']}")
                break
        ends = {o["end"] for o in orders}
        print(f"orders {len(orders)}, points {len(points)}, end reasons {sorted(ends)}")
        if problems:
            print("FAIL: " + "; ".join(problems))
            return 1
        print("OK")
        return 0


if __name__ == "__main__":
    sys.exit(main())
