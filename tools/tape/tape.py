#!/usr/bin/env python3
"""Reader for the per-execution tape written by itch-trades.

File: gzip( 24-byte header + N x 64-byte records ), little-endian, layout in
include/itchbook/tape/record.hpp. Two readers with the same field order:

    iter_records(path)  standard library only; yields tuples (slow, for tests)
    read_tape(path)     numpy structured array (needs numpy)

    python3 tape.py FILE.tape.gz      # print a summary
"""
from __future__ import annotations

import gzip
import struct
import sys
from collections import Counter
from typing import Iterator

MAGIC = b"ITCHTAPE"
VERSION = 1
HEADER = struct.Struct("<8sIIII")  # magic, version, record_size, endian, reserved
RECORD = struct.Struct("<QQQIIIIIIIHcccbcBBB2x")
FIELDS = (
    "ts", "match", "ref", "seq", "shares", "price", "bid_px", "ask_px", "bid_sz", "ask_sz",
    "locate", "type", "printable", "side", "sign", "tstate", "phase", "at_head", "flags",
)
assert RECORD.size == 64 and HEADER.size == 24

# flags bits (record.hpp)
SHARES_SATURATED, BID_SZ_SATURATED, ASK_SZ_SATURATED = 1, 2, 4
OVER_EXECUTE, ORPHAN, OPERATIONAL_HALT = 8, 16, 32
PRE, REGULAR, POST = 0, 1, 2


def _check_header(b: bytes) -> None:
    magic, version, size, endian, _ = HEADER.unpack(b)
    if magic != MAGIC:
        raise ValueError("not an ITCH tape (bad magic)")
    if version != VERSION or size != RECORD.size:
        raise ValueError(f"unsupported tape: version {version}, record size {size}")
    if endian != 0x01020304:
        raise ValueError("tape written on a big-endian host")


def iter_records(path: str) -> Iterator[tuple]:
    with gzip.open(path, "rb") as f:
        _check_header(f.read(HEADER.size))
        while True:
            chunk = f.read(RECORD.size * 65536)
            if not chunk:
                return
            if len(chunk) % RECORD.size:
                raise ValueError("truncated tape")
            for rec in RECORD.iter_unpack(chunk):
                # chars come back as 1-byte bytes; decode for convenience
                r = list(rec)
                for i in (11, 12, 13, 15):
                    r[i] = r[i].decode("latin-1")
                yield tuple(r)


def numpy_dtype():
    import numpy as np

    return np.dtype(
        {
            "names": list(FIELDS),
            "formats": ["<u8", "<u8", "<u8", "<u4", "<u4", "<u4", "<u4", "<u4", "<u4", "<u4",
                        "<u2", "S1", "S1", "S1", "i1", "S1", "u1", "u1", "u1"],
            "offsets": [0, 8, 16, 24, 28, 32, 36, 40, 44, 48, 52, 54, 55, 56, 57, 58, 59, 60, 61],
            "itemsize": 64,
        }
    )


def read_tape(path: str):
    import numpy as np

    with gzip.open(path, "rb") as f:
        _check_header(f.read(HEADER.size))
        raw = f.read()
    if len(raw) % RECORD.size:
        raise ValueError("truncated tape")
    return np.frombuffer(raw, dtype=numpy_dtype())


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    types: Counter = Counter()
    n = 0
    for r in iter_records(sys.argv[1]):
        types[r[11]] += 1
        n += 1
    print(f"records {n}  " + "  ".join(f"{k} {v}" for k, v in sorted(types.items())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
