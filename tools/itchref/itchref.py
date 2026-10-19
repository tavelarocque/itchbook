#!/usr/bin/env python3
"""Independent reference parser for NASDAQ TotalView-ITCH 5.0.

This is deliberately a separate implementation from the C++ engine: layouts
are written as sequential `struct` format strings (the C++ views use explicit
byte offsets), framing is re-implemented here, and nothing is shared with the
C++ code. It prints exactly the canonical text that `itch-stats` prints, so the
two outputs can be compared with `diff`.

    python3 itchref.py [--fingerprint] [--max-messages N] FILE

Standard library only; slow (about a microsecond per message) but simple.
"""

from __future__ import annotations

import argparse
import gzip
import struct
import sys
from typing import BinaryIO

# Body layout after the 1-byte type: (struct format, field names). The common
# header is H (stock locate), H (tracking number), 6s (timestamp).
HDR = "HH6s"
HDR_FIELDS = ["stock_locate", "tracking_number", "timestamp"]

LAYOUTS: dict[str, tuple[str, list[str]]] = {
    "S": ("c", ["event_code"]),
    "R": (
        "8sccIcc2scccccIc",
        [
            "stock",
            "market_category",
            "financial_status",
            "round_lot_size",
            "round_lots_only",
            "issue_classification",
            "issue_sub_type",
            "authenticity",
            "short_sale_threshold",
            "ipo_flag",
            "luld_tier",
            "etp_flag",
            "etp_leverage_factor",
            "inverse_indicator",
        ],
    ),
    "H": ("8scc4s", ["stock", "trading_state", "reserved", "reason"]),
    "Y": ("8sc", ["stock", "reg_sho_action"]),
    "L": ("4s8sccc", ["mpid", "stock", "primary_market_maker", "market_maker_mode", "market_participant_state"]),
    "V": ("QQQ", ["level1", "level2", "level3"]),
    "W": ("c", ["breached_level"]),
    "K": ("8sIcI", ["stock", "release_time", "release_qualifier", "ipo_price"]),
    "J": ("8sIIII", ["stock", "reference_price", "upper_price", "lower_price", "extension"]),
    "h": ("8scc", ["stock", "market_code", "halt_action"]),
    "A": ("QcI8sI", ["order_ref", "side", "shares", "stock", "price"]),
    "F": ("QcI8sI4s", ["order_ref", "side", "shares", "stock", "price", "attribution"]),
    "E": ("QIQ", ["order_ref", "executed_shares", "match_number"]),
    "C": ("QIQcI", ["order_ref", "executed_shares", "match_number", "printable", "execution_price"]),
    "X": ("QI", ["order_ref", "cancelled_shares"]),
    "D": ("Q", ["order_ref"]),
    "U": ("QQII", ["original_order_ref", "new_order_ref", "shares", "price"]),
    "P": ("QcI8sIQ", ["order_ref", "side", "shares", "stock", "price", "match_number"]),
    "Q": ("Q8sIQc", ["shares", "stock", "cross_price", "match_number", "cross_type"]),
    "B": ("Q", ["match_number"]),
    "I": (
        "QQc8sIIIcc",
        [
            "paired_shares",
            "imbalance_shares",
            "imbalance_direction",
            "stock",
            "far_price",
            "near_price",
            "current_reference_price",
            "cross_type",
            "price_variation_indicator",
        ],
    ),
    "N": ("8sc", ["stock", "interest_flag"]),
    "O": (
        "8scIIIQII",
        [
            "stock",
            "open_eligibility",
            "min_allowable_price",
            "max_allowable_price",
            "near_execution_price",
            "near_execution_time",
            "lower_price_range_collar",
            "upper_price_range_collar",
        ],
    ),
}

# Big-endian, no padding; offset 1 skips the type byte.
STRUCTS = {t: struct.Struct(">" + HDR + fmt) for t, (fmt, _) in LAYOUTS.items()}
FIELDS = {t: HDR_FIELDS + names for t, (_, names) in LAYOUTS.items()}
LENGTHS = {t: 1 + s.size for t, s in STRUCTS.items()}  # total length, type byte included

MASK64 = (1 << 64) - 1


def as_int(v: object) -> int:
    if isinstance(v, bytes):
        return int.from_bytes(v, "big")
    assert isinstance(v, int)
    return v


def open_any(path: str) -> BinaryIO:
    f = sys.stdin.buffer if path == "-" else open(path, "rb")
    head = f.peek(2)[:2] if hasattr(f, "peek") else b""
    if head == b"\x1f\x8b":
        return gzip.GzipFile(fileobj=f)  # type: ignore[return-value]
    return f  # type: ignore[return-value]


def parse(path: str, fingerprint: bool, max_messages: int | None) -> dict:
    counts: dict[str, int] = {}
    sums: dict[str, list[int]] = {t: [0] * len(FIELDS[t]) for t in LAYOUTS}
    st = dict(frames=0, messages=0, end_of_session_frames=0, unknown_type=0, bad_length=0, truncated_tail_bytes=0)
    stopped = False
    f = open_any(path)
    buf = b""
    chunk = 1 << 22
    source_truncated = False
    while not stopped:
        try:
            data = f.read(chunk)
        except EOFError:  # gzip stream cut short (partial download)
            source_truncated = True
            data = b""
        if not data:
            break
        buf += data
        pos = 0
        n = len(buf)
        while pos + 2 <= n:
            ln = (buf[pos] << 8) | buf[pos + 1]
            if pos + 2 + ln > n:
                break
            if ln and max_messages is not None and st["messages"] == max_messages:
                stopped = True
                break
            start = pos + 2
            pos = start + ln
            st["frames"] += 1
            if ln == 0:
                st["end_of_session_frames"] += 1
                continue
            t = chr(buf[start])
            want = LENGTHS.get(t)
            if want is None:
                st["unknown_type"] += 1
                continue
            if ln != want:
                st["bad_length"] += 1
                continue
            st["messages"] += 1
            counts[t] = counts.get(t, 0) + 1
            if fingerprint:
                acc = sums[t]
                for i, v in enumerate(STRUCTS[t].unpack_from(buf, start + 1)):
                    acc[i] += as_int(v)
        buf = buf[pos:]
    if not stopped:
        st["truncated_tail_bytes"] = len(buf)
    return dict(stats=st, counts=counts, sums=sums, stopped=stopped, source_truncated=source_truncated)


def render(result: dict, fingerprint: bool) -> str:
    st = result["stats"]
    out = [
        f"{k}\t{st[k]}"
        for k in ("frames", "messages", "end_of_session_frames", "unknown_type", "bad_length", "truncated_tail_bytes")
    ]
    for t in sorted(result["counts"], key=ord):
        out.append(f"count\t{t}\t{result['counts'][t]}")
    if fingerprint:
        for t in sorted(result["counts"], key=ord):
            for name, s in zip(FIELDS[t], result["sums"][t], strict=True):
                out.append(f"fp\t{t}\t{name}\t{s & MASK64}")
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("file")
    ap.add_argument("--fingerprint", action="store_true")
    ap.add_argument("--max-messages", type=int)
    a = ap.parse_args()
    sys.stdout.write(render(parse(a.file, a.fingerprint, a.max_messages), a.fingerprint))
    return 0


if __name__ == "__main__":
    sys.exit(main())
