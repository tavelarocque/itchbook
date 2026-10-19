"""Book cross-check tests, including mutation tests of the check itself.

Each mutation plants a plausible book bug in the Python reference; the
cross-check against the C++ book must then fail. If one ever passes, the
synthetic stream or the checksum has a blind spot.

    ITCH_BOOK=build/release/itch-book uv run pytest tools/itchref
"""

from __future__ import annotations

import gzip
import os
import subprocess
import tempfile
from pathlib import Path

import pytest

import bookref
import crosscheck_book

ITCH_BOOK = os.environ.get("ITCH_BOOK")
needs_binary = pytest.mark.skipif(not ITCH_BOOK, reason="set ITCH_BOOK to the itch-book binary")


def test_splitmix_known_values():
    # SplitMix64 reference outputs for seed 0 (first three draws of the generator).
    golden = 0x9E3779B97F4A7C15
    outs = [bookref.splitmix64(k * golden & bookref.MASK) for k in range(3)]
    assert outs == [0xE220A8397B1DCDAF, 0x6E789E6AA1B965F4, 0x06C45D188009454F]


@pytest.fixture(scope="module")
def stream_and_cpp():
    raw = crosscheck_book.synthetic_stream(20000, seed=11)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "s.itch.gz"
        p.write_bytes(gzip.compress(raw))
        out = ""
        if ITCH_BOOK:
            out = subprocess.run(
                [ITCH_BOOK, "--symbols", ",".join(crosscheck_book.SYMBOLS), str(p)], capture_output=True, text=True
            ).stdout
        yield str(p), out


ORIG = {m: getattr(bookref.Book, m) for m in ("replace", "cancel", "execute")}


def py(path):
    return bookref.run(path, 30, crosscheck_book.SYMBOLS, 5, None)


@needs_binary
def test_unmutated_reference_matches(stream_and_cpp):
    path, cpp = stream_and_cpp
    assert py(path) == cpp


def _replace_keeps_priority(self, locate, old, new, shares, price):
    # Bug: rename the order in place instead of re-queueing it at the tail.
    o = self.orders.get(old)
    if o is None or (new != old and new in self.orders):
        return ORIG["replace"](self, locate, old, new, shares, price)
    self.st["replaces"] += 1
    lvl = self.levels[(o[0], o[1], o[2])]
    if o[2] == price:
        q = lvl[1]
        items = [(new if r == old else r) for r in q]
        q.clear()
        q.update(dict.fromkeys(items))
        lvl[0] += shares - o[3]
        del self.orders[old]
        self.orders[new] = [o[0], o[1], price, shares]
    else:
        ORIG["replace"](self, locate, old, new, shares, price)


def _cancel_forgets_level(self, locate, ref, n):
    # Bug: order shrinks but its level total does not.
    o = self.orders.get(ref)
    if o is not None and n < o[3]:
        self.st["cancels"] += 1
        o[3] -= n
        return
    ORIG["cancel"](self, locate, ref, n)


def _execute_ignores_clamp(self, locate, ref, n):
    # Bug: an execution larger than the order is not counted as an anomaly.
    o = self.orders.get(ref)
    if o is not None and n > o[3]:
        n = o[3]
    ORIG["execute"](self, locate, ref, n)


@needs_binary
@pytest.mark.parametrize(
    "method,mutant",
    [("replace", _replace_keeps_priority), ("cancel", _cancel_forgets_level), ("execute", _execute_ignores_clamp)],
)
def test_mutations_are_detected(stream_and_cpp, monkeypatch, method, mutant):
    path, cpp = stream_and_cpp
    monkeypatch.setattr(bookref.Book, method, mutant)
    assert py(path) != cpp
