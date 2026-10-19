#!/usr/bin/env python3
"""Compare the C++ engine (itch-stats) against the Python reference (itchref).

Synthetic mode builds a stream containing every message type (random field
bytes, zero-length frames, an unknown type and a wrong-length frame), writes
it raw and gzipped, and requires all four outputs (C++ and Python, raw and
gzip) to be identical:

    crosscheck.py --itch-stats build/itch-stats --synthetic 20000

File mode runs both parsers over a real NASDAQ file and diffs them. With
--published it also compares message counts against a published table:

    crosscheck.py --itch-stats build/itch-stats --file 01302020.NASDAQ_ITCH50.gz \
        --fingerprint --published 01302020
"""

from __future__ import annotations

import argparse
import difflib
import gzip
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import itchref  # noqa: E402


def frame(payload: bytes) -> bytes:
    return len(payload).to_bytes(2, "big") + payload


def synthetic_stream(n: int, seed: int) -> bytes:
    rng = random.Random(seed)
    types = sorted(itchref.LENGTHS, key=ord)
    out = bytearray()
    order = types + [rng.choice(types) for _ in range(max(0, n - len(types)))]
    rng.shuffle(order)
    for i, t in enumerate(order):
        body = bytes(rng.getrandbits(8) for _ in range(itchref.LENGTHS[t] - 1))
        out += frame(t.encode() + body)
        if i % 997 == 0:
            out += b"\x00\x00"  # end-of-session marker mid-stream: must be skipped
    out += frame(b"~" + bytes(10))  # unknown type
    out += frame(b"A" + bytes(itchref.LENGTHS["A"]))  # one byte too long
    out += b"\x00\x00"
    return bytes(out)


def run_cpp(itch_stats: str, path: str, fingerprint: bool, max_messages: int | None) -> str:
    cmd = [itch_stats] + (["--fingerprint"] if fingerprint else [])
    if max_messages is not None:
        cmd += ["--max-messages", str(max_messages)]
    cmd.append(path)
    p = subprocess.run(cmd, capture_output=True, text=True, errors="backslashreplace")
    if p.returncode not in (0, 2):
        sys.stderr.write(p.stderr)
        raise SystemExit(f"itch-stats failed with exit code {p.returncode}")
    sys.stderr.write(p.stderr)
    return p.stdout


def run_py(path: str, fingerprint: bool, max_messages: int | None) -> str:
    return itchref.render(itchref.parse(path, fingerprint, max_messages), fingerprint)


def same(label_a: str, a: str, label_b: str, b: str) -> bool:
    if a == b:
        return True
    sys.stdout.writelines(difflib.unified_diff(a.splitlines(True), b.splitlines(True), label_a, label_b))
    return False


def check_published(name: str, text: str) -> bool:
    table = json.loads((HERE / "published_counts.json").read_text())[name]
    ours = {}
    for line in text.splitlines():
        parts = line.split("\t")
        if parts[0] == "count":
            ours[parts[1]] = int(parts[2])
    ok = True
    print(f"published counts ({table['source']}):")
    for t, want in sorted(table["counts"].items(), key=lambda kv: -kv[1]):
        got = ours.get(t, 0)
        status = "ok" if got == want else ("disputed" if t in table.get("disputed", {}) else "MISMATCH")
        ok &= status != "MISMATCH"
        print(f"  {t}  published {want:>13,}  ours {got:>13,}  {status}")
    unlisted = sorted(set(ours) - set(table["counts"]), key=ord)
    if unlisted:
        print("  types the published parser does not decode: " + ", ".join(f"{t}={ours[t]:,}" for t in unlisted))
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--itch-stats", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--synthetic", type=int, metavar="N")
    g.add_argument("--file")
    ap.add_argument("--fingerprint", action="store_true", help="file mode: also compare per-field sums")
    ap.add_argument("--max-messages", type=int)
    ap.add_argument("--published", help="file mode: compare counts with published_counts.json[NAME]")
    ap.add_argument("--seed", type=int, default=20261018)
    a = ap.parse_args()

    if a.synthetic is not None:
        raw = synthetic_stream(a.synthetic, a.seed)
        with tempfile.TemporaryDirectory() as d:
            rp, gp = Path(d) / "s.itch", Path(d) / "s.itch.gz"
            rp.write_bytes(raw)
            gp.write_bytes(gzip.compress(raw))
            outs = {
                "cpp-raw": run_cpp(a.itch_stats, str(rp), True, None),
                "cpp-gz": run_cpp(a.itch_stats, str(gp), True, None),
                "py-raw": run_py(str(rp), True, None),
                "py-gz": run_py(str(gp), True, None),
            }
        ok = all(same("cpp-raw", outs["cpp-raw"], k, v) for k, v in outs.items())
        n_types = sum(1 for line in outs["py-raw"].splitlines() if line.startswith("count\t"))
        ok &= n_types == len(itchref.LENGTHS)
        print(f"synthetic: {a.synthetic} messages, {n_types} types: {'MATCH' if ok else 'MISMATCH'}")
        return 0 if ok else 1

    cpp = run_cpp(a.itch_stats, a.file, a.fingerprint, a.max_messages)
    py = run_py(a.file, a.fingerprint, a.max_messages)
    ok = same("itch-stats", cpp, "itchref", py)
    print(f"file {a.file}: {'MATCH' if ok else 'MISMATCH'}")
    if a.published:
        ok &= check_published(a.published, cpp)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
