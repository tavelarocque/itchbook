#!/usr/bin/env python3
"""Compare the itch-trades tape with the independent Python reference (taperef).

Synthetic mode builds a stream exercising every record rule (sweeps sharing a
timestamp, non-printable C, crosses including one above 2^32 shares, broken
trades, halts, trading states, session events, orphans, over-executions,
header-locate mismatches, and NOII, Reg SHO and market-participant messages
that only advance the sequence number), runs itch-trades with BOTH engines, and requires:
  * the fast and reference tapes to be byte-identical after decompression;
  * every record to equal the Python reference, field by field, in order.
File mode does the same on the first N messages of a real day.

    tape_crosscheck.py --itch-trades build/release/itch-trades --synthetic 30000
    tape_crosscheck.py --itch-trades build/release/itch-trades --file DAY.gz --max-messages 3000000
"""
from __future__ import annotations

import argparse
import gzip
import random
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "tape"))

import tape  # noqa: E402
import taperef  # noqa: E402
from crosscheck_book import SYMBOLS, alpha, frame, hdr  # noqa: E402


def synthetic_stream(n: int, seed: int) -> bytes:
    rng = random.Random(seed)
    out = bytearray()
    ts = 4 * 3600 * 10**9
    step = (16 * 3600 * 10**9) // max(n, 1)
    out += frame(hdr("S", 0, ts) + b"O")
    for i, s in enumerate(SYMBOLS, start=1):
        body = alpha(s, 8) + b"QN" + (100).to_bytes(4, "big") + b"NCZ PNN1N" + bytes(4) + b"N"
        out += frame(hdr("R", i, ts) + body)
    live: dict[int, int] = {}
    issued: list[int] = []
    next_ref, match = 1, 1
    phase_events = {n // 4: b"Q", (3 * n) // 4: b"M"}
    for i in range(n):
        if rng.random() < 0.6:  # ties: about 40% of messages share the previous timestamp
            ts += rng.randint(1, 2 * step)
        if i in phase_events:
            out += frame(hdr("S", 0, ts) + phase_events[i])
        k = rng.random()
        loc = rng.randint(1, len(SYMBOLS))
        price = rng.randint(1000, 1008) * 100
        shares = rng.randint(1, 500)
        target = rng.choice(issued) if issued and rng.random() < 0.97 else next_ref + 10**6
        tloc = live.get(target, loc)
        if rng.random() < 0.01:
            tloc = rng.randint(1, len(SYMBOLS))  # header locate mismatch
        match += 1
        if k < 0.38 or not live:
            side = rng.choice(b"BS")
            big = rng.random() < 0.002
            sh = 3_000_000_000 if big else shares
            body = next_ref.to_bytes(8, "big") + bytes([side]) + sh.to_bytes(4, "big")
            body += alpha(SYMBOLS[loc - 1], 8) + price.to_bytes(4, "big")
            out += frame(hdr("A", loc, ts) + body)
            live[next_ref] = loc
            issued.append(next_ref)
            next_ref += 1
        elif k < 0.58:
            # a sweep: several executions with one timestamp on one symbol
            for _ in range(rng.choice([1, 1, 1, 2, 3, 5])):
                tgt = rng.choice(issued)
                tl = live.get(tgt, loc)
                match += 1
                if rng.random() < 0.75:
                    out += frame(hdr("E", tl, ts) + tgt.to_bytes(8, "big") + shares.to_bytes(4, "big")
                                 + match.to_bytes(8, "big"))
                else:
                    pr = b"Y" if rng.random() < 0.7 else b"N"
                    out += frame(hdr("C", tl, ts) + tgt.to_bytes(8, "big") + shares.to_bytes(4, "big")
                                 + match.to_bytes(8, "big") + pr + (price + 37).to_bytes(4, "big"))
        elif k < 0.66:
            out += frame(hdr("X", tloc, ts) + target.to_bytes(8, "big") + shares.to_bytes(4, "big"))
        elif k < 0.80:
            out += frame(hdr("D", tloc, ts) + target.to_bytes(8, "big"))
            live.pop(target, None)
        elif k < 0.90:
            new = next_ref if rng.random() > 0.02 else rng.choice(issued)
            body = target.to_bytes(8, "big") + new.to_bytes(8, "big") + rng.randint(0, 500).to_bytes(4, "big")
            out += frame(hdr("U", tloc, ts) + body + price.to_bytes(4, "big"))
            if target in live and (new == target or new not in live):
                live[new] = live.pop(target)
            issued.append(next_ref)
            next_ref += 1
        elif k < 0.95:
            body = bytes(8) + b"B" + shares.to_bytes(4, "big") + alpha(SYMBOLS[loc - 1], 8)
            out += frame(hdr("P", loc, ts) + body + (price + 50).to_bytes(4, "big") + match.to_bytes(8, "big"))
        elif k < 0.965:
            qsh = rng.choice([0, 12345, 5_000_000_000])
            body = qsh.to_bytes(8, "big") + alpha(SYMBOLS[loc - 1], 8) + price.to_bytes(4, "big")
            out += frame(hdr("Q", loc, ts) + body + match.to_bytes(8, "big") + rng.choice([b"O", b"C", b"H"]))
        elif k < 0.975:
            out += frame(hdr("B", loc, ts) + rng.randint(1, match).to_bytes(8, "big"))
        elif k < 0.982:
            st = rng.choice([b"T", b"H", b"P", b"Q"])
            out += frame(hdr("H", loc, ts) + alpha(SYMBOLS[loc - 1], 8) + st + b" " + b"    ")
        elif k < 0.987:
            out += frame(hdr("h", loc, ts) + alpha(SYMBOLS[loc - 1], 8) + rng.choice([b"Q", b"B", b"X"])
                         + rng.choice([b"H", b"T"]))
        else:  # non-book traffic handled by no specific handler: NOII, Reg SHO, market participant
            sym = alpha(SYMBOLS[loc - 1], 8)
            kind = rng.choice("IYL")
            if kind == "I":
                body = (rng.randint(0, 10**6).to_bytes(8, "big") + rng.randint(0, 10**6).to_bytes(8, "big") + b"B"
                        + sym + price.to_bytes(4, "big") * 3 + b"O" + b" ")
            elif kind == "Y":
                body = sym + b"0"
            else:
                body = b"ABCD" + sym + b"YNA"
            out += frame(hdr(kind, loc, ts) + body)
    out += frame(hdr("S", 0, ts + 1) + b"C")
    return bytes(out)


def run_cpp(binary: str, engine: str, path: str, prefix: str, extra: list[str]) -> Path:
    p = subprocess.run([binary, "--engine", engine, *extra, "--out", prefix, path], capture_output=True, text=True)
    sys.stderr.write(p.stderr)
    if p.returncode != 0:
        raise SystemExit(f"itch-trades ({engine}) failed with exit code {p.returncode}")
    return Path(prefix + ".tape.gz")


def compare(cpp_path: Path, ref_iter) -> tuple[int, bool]:
    n = 0
    cpp = tape.iter_records(str(cpp_path))
    for want in ref_iter:
        got = next(cpp, None)
        if got != want:
            print(f"MISMATCH at record {n}")
            if got is None:
                print("  C++ tape ended early")
            else:
                for name, g, w in zip(tape.FIELDS, got, want):
                    print(f"  {name:10s} cpp={g!r:24} ref={w!r}{'   <--' if g != w else ''}")
            return n, False
        n += 1
    extra = next(cpp, None)
    if extra is not None:
        print(f"MISMATCH: C++ tape has extra records from {n}: {extra}")
        return n, False
    return n, True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--itch-trades", required=True)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--synthetic", type=int, metavar="N")
    g.add_argument("--file")
    ap.add_argument("--seed", type=int, default=20261026)
    ap.add_argument("--max-messages", type=int)
    ap.add_argument("--engines", default="fast,reference")
    a = ap.parse_args()
    engines = a.engines.split(",")
    with tempfile.TemporaryDirectory() as d:
        if a.synthetic:
            src = Path(d) / "synthetic.itch.gz"
            src.write_bytes(gzip.compress(synthetic_stream(a.synthetic, a.seed)))
            path, extra, lim = str(src), [], None
        else:
            path, lim = a.file, a.max_messages
            extra = ["--max-messages", str(lim)] if lim else []
        tapes = {e: run_cpp(a.itch_trades, e, path, str(Path(d) / e), extra) for e in engines}
        ok = True
        if len(tapes) == 2:
            same = gzip.decompress(tapes["fast"].read_bytes()) == gzip.decompress(tapes["reference"].read_bytes())
            print(f"fast vs reference tapes: {'IDENTICAL' if same else 'DIFFERENT'}")
            ok &= same
        n, match = compare(tapes[engines[0]], taperef.build(path, lim))
        types: dict[str, int] = {}
        for r in tape.iter_records(str(tapes[engines[0]])):
            types[r[11]] = types.get(r[11], 0) + 1
        print(f"tape vs Python reference: {n} records ({', '.join(f'{k} {v}' for k, v in sorted(types.items()))}): "
              f"{'MATCH' if match else 'MISMATCH'}")
        ok &= match
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
