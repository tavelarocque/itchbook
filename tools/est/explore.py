"""The one gate through which Phase 3 code reads real data: exploration stocks only.

data/universe.tsv assigns every eligible stock a v2 cell ('Q-T1S1' .. 'N-T3S3',
column `cell`) and a split (`split_v2`: 'explore' or 'confirm'). The study days
are the table's own header line '# study_days: ...'.

load_explore_runs(day, cell=None, symbols=None) returns {symbol: runs} for the
exploration stocks of one study day (all cells, one cell, or the named
symbols), each value one stock-day of runs.RUN_DTYPE records in time order. It
raises ValueError, before opening any file, for a day outside the study days,
an unknown cell, or any requested symbol that is not an exploration stock
(confirmation stocks, ineligible symbols, unknown symbols).

A runs file is one compressed array holding every stock of the day, so its
bytes must be inflated to be read; it is streamed in chunks of CHUNK records
and each chunk keeps only rows whose locate maps (through that day's
tape/{day}.symbols.tsv) to an exploration symbol, so no confirmation row
outlives the chunk buffer it was inflated into; a final check raises if any
unrequested locate reached the output. Nothing here reads the
replication .gz files.

Data directory: ITCH_DATA, default ~/CODE/itch-data.
"""
from __future__ import annotations

import os
import sys
import zipfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tape"))
import runs as R  # noqa: E402

UNIVERSE = Path(__file__).resolve().parents[2] / "data" / "universe.tsv"
CHUNK = 1 << 20


def data_dir() -> Path:
    return Path(os.environ.get("ITCH_DATA", "~/CODE/itch-data")).expanduser()


def read_universe(path: Path = UNIVERSE) -> tuple[tuple[str, ...], dict[str, tuple[str, str]]]:
    """(study_days, {symbol: (cell, split_v2)}) for every symbol with a v2 cell."""
    days, out, head = None, {}, None
    with open(path) as f:
        for line in f:
            if line.startswith("# study_days:"):
                days = tuple(line.split(":", 1)[1].split())
                continue
            if line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if head is None:
                head = {k: i for i, k in enumerate(row)}
                continue
            cell, split = row[head["cell"]], row[head["split_v2"]]
            if cell:
                out[row[head["symbol"]]] = (cell, split)
    if not days:
        raise ValueError(f"{path}: no '# study_days:' header")
    return days, out


def explore_symbols(cell: str | None = None, universe: Path = UNIVERSE) -> list[str]:
    _, u = read_universe(universe)
    cells = {c for c, _ in u.values()}
    if cell is not None and cell not in cells:
        raise ValueError(f"unknown cell {cell!r}")
    return sorted(s for s, (c, sp) in u.items() if sp == "explore" and (cell is None or c == cell))


def cells(universe: Path = UNIVERSE) -> list[str]:
    _, u = read_universe(universe)
    return sorted({c for c, _ in u.values()})


def _check(day: str, cell: str | None, symbols, universe: Path) -> list[str]:
    days, u = read_universe(universe)
    if day not in days:
        raise ValueError(f"{day!r} is not a study day")
    if cell is not None and cell not in {c for c, _ in u.values()}:
        raise ValueError(f"unknown cell {cell!r}")
    allowed = [s for s, (c, sp) in u.items() if sp == "explore" and (cell is None or c == cell)]
    if symbols is None:
        return sorted(allowed)
    bad = sorted(set(symbols) - set(allowed))
    if bad:
        why = {s: (u[s][1] if s in u else "not in a cell") for s in bad[:5]}
        raise ValueError(f"not exploration stocks{'' if cell is None else ' of ' + cell}: {why}")
    return sorted(set(symbols))


def _read_symbols(path: Path) -> dict[str, int]:
    """symbol -> locate of one day (tape/{day}.symbols.tsv)."""
    out = {}
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        if head[:2] != ["locate", "symbol"]:
            raise ValueError(f"unexpected symbols header {head}")
        for line in f:
            loc, sym = line.split("\t", 2)[:2]
            out[sym] = int(loc)
    return out


def _stream_filter(npz: Path, keep_locate: np.ndarray) -> np.ndarray:
    """Rows of runs.npy inside npz whose locate is flagged in keep_locate (bool, 65536)."""
    with zipfile.ZipFile(npz) as z, z.open("runs.npy") as f:
        version = np.lib.format.read_magic(f)
        read = np.lib.format.read_array_header_1_0 if version == (1, 0) else np.lib.format.read_array_header_2_0
        shape, fortran, dtype = read(f)
        if dtype != R.RUN_DTYPE or fortran or len(shape) != 1:
            raise ValueError(f"{npz}: unexpected array {dtype} {shape}")
        n, parts = shape[0], []
        for lo in range(0, n, CHUNK):
            k = min(CHUNK, n - lo)
            buf = f.read(k * dtype.itemsize)
            a = np.frombuffer(buf, dtype=dtype, count=k)
            parts.append(a[keep_locate[a["locate"]]].copy())
            del a, buf
    return np.concatenate(parts) if parts else np.zeros(0, R.RUN_DTYPE)


def load_explore_runs(day: str, cell: str | None = None, symbols=None, data: Path | None = None,
                      universe: Path = UNIVERSE) -> dict[str, np.ndarray]:
    want = _check(day, cell, symbols, universe)
    data = data_dir() if data is None else Path(data)
    sym_loc = _read_symbols(data / "tape" / f"{day}.symbols.tsv")
    keep = np.zeros(1 << 16, dtype=bool)
    loc_sym = {}
    for s in want:
        if s in sym_loc:
            keep[sym_loc[s]] = True
            loc_sym[sym_loc[s]] = s
    runs = _stream_filter(data / "runs" / f"{day}.runs.npz", keep)
    if not np.all(keep[runs["locate"]]):  # second line of defence: nothing unrequested leaves this function
        raise RuntimeError("runs of a stock that was not requested")
    out = {}
    if len(runs):
        loc = runs["locate"]
        cut = np.flatnonzero(np.diff(loc.astype(np.int64))) + 1
        for part in np.split(runs, cut):
            out[loc_sym[int(part["locate"][0])]] = part
    for s in out.values():
        if np.any(np.diff(s["ts"].astype(np.int64)) < 0):
            raise ValueError("runs of a stock-day out of time order")
    return out
