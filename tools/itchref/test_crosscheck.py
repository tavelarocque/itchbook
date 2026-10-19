"""Mutation tests for the cross-check itself.

A cross-check is only worth something if it fails when one side is wrong. For
every message type, swap two adjacent fields in the Python layout (a classic
offset bug) and require the fingerprint to stop matching the C++ engine.

    ITCH_STATS=build/itch-stats uv run pytest tools/itchref
"""

from __future__ import annotations

import os
import struct
import tempfile
from pathlib import Path

import pytest

import crosscheck
import itchref

ITCH_STATS = os.environ.get("ITCH_STATS")
needs_binary = pytest.mark.skipif(not ITCH_STATS, reason="set ITCH_STATS to the itch-stats binary")


# Split a struct format into per-field tokens, e.g. "8sccI" -> ["8s", "c", "c", "I"].
def tokens(fmt: str) -> list[str]:
    out, num = [], ""
    for ch in fmt:
        if ch.isdigit():
            num += ch
        else:
            out.append(num + ch)
            num = ""
    return out


def test_lengths_match_spec():
    spec = dict(
        S=12,
        R=39,
        H=25,
        Y=20,
        L=26,
        V=35,
        W=12,
        K=28,
        J=35,
        h=21,
        A=36,
        F=40,
        E=31,
        C=36,
        X=23,
        D=19,
        U=35,
        P=44,
        Q=40,
        B=19,
        I=50,
        N=20,
        O=48,
    )
    assert itchref.LENGTHS == spec


def test_every_format_has_one_name_per_field():
    for t, (fmt, names) in itchref.LAYOUTS.items():
        assert len(tokens(fmt)) == len(names), t


@pytest.fixture(scope="module")
def stream_file():
    raw = crosscheck.synthetic_stream(3000, seed=7)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "s.itch"
        p.write_bytes(raw)
        yield str(p)


@needs_binary
def test_unmutated_reference_matches(stream_file):
    assert crosscheck.run_cpp(ITCH_STATS, stream_file, True, None) == crosscheck.run_py(stream_file, True, None)


@needs_binary
@pytest.mark.parametrize("msg_type", sorted(itchref.LAYOUTS, key=ord))
def test_swapping_adjacent_fields_is_detected(stream_file, msg_type, monkeypatch):
    cpp = crosscheck.run_cpp(ITCH_STATS, stream_file, True, None)
    fmt, names = itchref.LAYOUTS[msg_type]
    toks = tokens(itchref.HDR + fmt)
    all_names = itchref.HDR_FIELDS + names
    assert len(toks) >= 2
    # Swap the formats of the last two fields but keep their names: the bytes
    # are now read at the wrong offsets, as a mis-transcribed spec would.
    i = len(toks) - 2
    if struct.calcsize(">" + toks[i]) == struct.calcsize(">" + toks[i + 1]):
        # Same width: swapping formats would read identical bytes, so swap the
        # names instead, which mislabels the two sums.
        mutated_names = all_names[:i] + [all_names[i + 1], all_names[i]] + all_names[i + 2 :]
        mutated_fmt = toks
    else:
        mutated_names = all_names
        mutated_fmt = toks[:i] + [toks[i + 1], toks[i]] + toks[i + 2 :]
    monkeypatch.setitem(itchref.STRUCTS, msg_type, struct.Struct(">" + "".join(mutated_fmt)))
    monkeypatch.setitem(itchref.FIELDS, msg_type, mutated_names)
    assert crosscheck.run_py(stream_file, True, None) != cpp
