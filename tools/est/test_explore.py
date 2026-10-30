"""Tests for explore.py, the only gate from Phase 3 code to real data.

1. On a constructed data directory (a small universe table, symbol tables and
   runs files in a temporary directory): the guard refuses confirmation,
   ineligible and unknown symbols, non-study days (even when that day's files
   exist), unknown cells, and returns exactly the exploration stocks' runs,
   also when the runs file is streamed in chunks smaller than the array.
2. On the committed data/universe.tsv (no data files read): its study days
   are describe.STUDY_DAYS, and it refuses a real confirmation stock.
3. Planted bugs.

    uv run --group dev --group analysis pytest tools/est/test_explore.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import describe as D
import explore as E
import runs as R

DAYS = ("01302019", "01302020")
UNIV = [("AAA", "X", "explore"), ("BBB", "X", "confirm"), ("CCC", "Y", "explore"), ("DDD", "", ""),
        ("EEE", "Y", "confirm")]


def refused(f, *args, **kw) -> None:
    """f must raise ValueError (an AssertionError otherwise, so a planted bug's test can catch it)."""
    try:
        f(*args, **kw)
    except ValueError:
        return
    raise AssertionError(f"not refused: {args} {kw}")


def make_data(tmp, days=DAYS + ("12312019",)):
    """A data directory with every symbol present on every day (locates differ by day)."""
    u = tmp / "universe.tsv"
    u.write_text(f"# study_days: {' '.join(DAYS)}\n# a comment\nsymbol\tlisting\tcell\tsplit_v2\n"
                 + "".join(f"{s}\tQ\t{c}\t{sp}\n" for s, c, sp in UNIV))
    (tmp / "tape").mkdir()
    (tmp / "runs").mkdir()
    rows = {}
    for k, d in enumerate(days):
        locs = {s: 10 + 3 * i + k for i, (s, _, _) in enumerate(UNIV)}
        (tmp / "tape" / f"{d}.symbols.tsv").write_text(
            "locate\tsymbol\tmarket_category\tissue_classification\tetp_flag\n"
            + "".join(f"{locs[s]}\t{s}\tQ\tC\tN\n" for s in locs))
        parts = []
        for i, s in enumerate(locs):
            n = 5 + i
            r = np.zeros(n, R.RUN_DTYPE)
            r["locate"] = locs[s]
            r["ts"] = np.arange(n) * 1000 + 7
            r["q"] = 100 * (i + 1) + np.arange(n)
            r["sign"] = np.where(np.arange(n) % 2, 1, -1)
            parts.append(r)
            rows[(d, s)] = r
        np.savez_compressed(tmp / "runs" / f"{d}.runs.npz", runs=np.concatenate(parts))
    return u, rows


def check_guard_refuses_confirmation_stocks(tmp):
    u, _ = make_data(tmp)
    for bad in (["BBB"], ["AAA", "BBB"], ["EEE"], ["DDD"], ["ZZZ"]):
        refused(E.load_explore_runs, DAYS[0], symbols=bad, data=tmp, universe=u)
    refused(E.load_explore_runs, DAYS[0], cell="Y", symbols=["AAA"], data=tmp, universe=u)  # explore, other cell
    refused(E.load_explore_runs, DAYS[0], cell="Z", data=tmp, universe=u)
    assert E.explore_symbols(universe=u) == ["AAA", "CCC"]


def check_guard_refuses_non_study_days(tmp):
    u, _ = make_data(tmp)
    for d in ("12312019", "01302021", "", "1302019"):
        refused(E.load_explore_runs, d, data=tmp, universe=u)


def check_guard_returns_exploration_runs_only(tmp):
    u, rows = make_data(tmp)
    for chunk in (E.CHUNK, 3):
        E.CHUNK = chunk
        try:
            for d in DAYS:
                got = E.load_explore_runs(d, data=tmp, universe=u)
                assert sorted(got) == ["AAA", "CCC"]
                for s, r in got.items():
                    assert np.array_equal(r, rows[(d, s)])
                got = E.load_explore_runs(d, cell="X", data=tmp, universe=u)
                assert sorted(got) == ["AAA"]
                got = E.load_explore_runs(d, symbols=["CCC"], data=tmp, universe=u)
                assert sorted(got) == ["CCC"]
        finally:
            E.CHUNK = 1 << 20


def check_committed_universe(tmp):
    days, u = E.read_universe()
    assert days == D.STUDY_DAYS
    confirm = sorted(s for s, (_, sp) in u.items() if sp == "confirm")
    assert confirm and all(u[s][1] == "explore" for s in E.explore_symbols())
    refused(E.load_explore_runs, days[0], symbols=confirm[:1], data=tmp)  # before any file is opened


CHECKS = {f.__name__: f for f in [
    check_guard_refuses_confirmation_stocks, check_guard_refuses_non_study_days,
    check_guard_returns_exploration_runs_only, check_committed_universe,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_explore(name, tmp_path):
    CHECKS[name](tmp_path)


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("guard lets confirmation stocks through", "_check",
        'if sp == "explore" and', 'if sp in ("explore", "confirm") and',
        ("check_guard_refuses_confirmation_stocks", "check_guard_returns_exploration_runs_only")),
    Bug("guard accepts any day", "_check",
        "if day not in days:", "if False:",
        ("check_guard_refuses_non_study_days",)),
    Bug("chunk filter keeps every row", "_stream_filter",
        'parts.append(a[keep_locate[a["locate"]]].copy())', "parts.append(a.copy())",
        ("check_guard_returns_exploration_runs_only",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate, tmp_path):
    mutate(E, bug.func, bug.old, bug.new)
    for i, name in enumerate(bug.caught_by):
        d = tmp_path / str(i)
        d.mkdir()
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name](d)
