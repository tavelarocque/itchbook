"""Phase 4: the confirmatory run, executed once.

The plan is frozen at tag prereg-freeze-v1 (commit d126088fba041824950d7ffae4aa75a7a98b2703) and registered on
OSF (osf.io/z4t7x; the registered PDF, Frozen_Preregistration_itchbook_2026-11-02.pdf, SHA-256
79c8984ab832c0a9a0d858ac46959277399ce6394106fb526fe9bd48779ccafc). This driver adds no statistics: it reads the
runs, builds the Flow inputs the frozen functions expect, calls primary3d.run_primary once with the frozen P-M
design and seed, checks the family, and writes what came out.

Data gate (load_runs). It reads one study day's runs of the requested NASDAQ-listed cells for ONE split
('confirm' or 'explore'), streamed and filtered exactly as explore.load_explore_runs (its frozen helpers
_read_symbols and _stream_filter). check_request raises ValueError before any file is opened for a day outside
the study days, a cell that is not one of the nine NASDAQ-listed cells, or a symbol of another split (an
exploration stock in a confirmation request and the reverse), another cell, or no cell; a second check raises
if any returned symbol is not of the requested split. Only runs/{day}.runs.npz and tape/{day}.symbols.tsv of
the eight study days are opened; no replication file is.

Inputs (to_flow), per stock-day in study-day order, stocks sorted within a cell:
  eps = the run's sign; q = its RAW size in shares (the R2 clip is applied by primary3d.make_unit);
  lm = log of the mid read before the run's first fill, mid2 / 2 in dollars, NaN where the book is locked,
  crossed or one-sided; ok = mid_valid; ts = ns since midnight ET; lengths = runs per stock-day.
  A stock without runs on a study day, a zero sign or a non-positive size raises (nothing is dropped).

The run (run_family). primary3d.run_primary(flows, FROZEN_PM, default_rng(SEED)) with its defaults (B = 200,
four surrogate draws). The flows hold all nine NASDAQ-listed cells: tick_bin_unit pools each tick bin's
confirmatory cells and logs the descriptive cells it drops; Q-T2 (Q-T2S3) is computed and never in Holm.
Per-cell results, the descriptive cells' included, are secondary in the frozen plan and are not computed here.
The p-value of each test is the largest of three methods' p-values; observe() records each method's p-value
and interval from the same replicates the frozen decision functions receive (frozen helpers only, no draws
taken), and the record is checked against the frozen p-values. verify_family checks the family and the Holm
decisions before anything is written.

    uv run --group analysis python tools/est/phase4.py dryrun    # exploration stocks; plumbing record only
    uv run --group analysis python tools/est/phase4.py confirm   # the confirmatory run (refuses a second run)
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import subprocess
import sys
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np

import explore as E
import inference as INF
import prereg3d as C
import primary3d as P

REPO = Path(__file__).resolve().parents[2]
FROZEN_TAG = "prereg-freeze-v1"
FROZEN_COMMIT = "d126088fba041824950d7ffae4aa75a7a98b2703"
OSF = "osf.io/z4t7x"
PLAN_PDF_SHA256 = "79c8984ab832c0a9a0d858ac46959277399ce6394106fb526fe9bd48779ccafc"
MANIFEST_SHA256 = "9780239827c00763fd542ad43cd9a474fbfe815a0418847a95a6328804bc5a9d"
UNIVERSE_SHA256 = "5fb7e1b7061ea8433d1916a3ca54c5f97a290e29b6c367f70881b248a3b42d25"
NEW_FILES = ("tools/est/phase4.py", "tools/est/test_phase4.py")  # the only files added after the freeze

SEED = 20261102
FROZEN_PM = P.PMParams(n_bins=3, l_max=16, L=128, min_runs_per_bin=50)  # b3_l16_L128_m50 (plan, page 4)
SPLITS = ("confirm", "explore")
NASDAQ_CELLS = tuple(sorted(set(C.CONFIRMATORY) | set(C.DESCRIPTIVE_CELLS)))
# Confirmation stocks per cell (plan, page 3); the confirmatory run stops if the universe gives others.
PLAN_CONFIRM_COUNTS = {"Q-T1S1": 161, "Q-T1S2": 216, "Q-T1S3": 252, "Q-T2S1": 80, "Q-T2S2": 50, "Q-T2S3": 17,
                       "Q-T3S1": 50, "Q-T3S2": 27, "Q-T3S3": 23}

OUT = REPO / "data" / "phase4"
DOC = REPO / "docs" / "phase4_results.md"
STARTED = "RUN_STARTED"

# Plan, page 6, "Interpretation, fixed in advance", copied exactly.
INTERPRETATION = {
    "pass": ("A P-V, P-R or P-M pass supports that prediction at the stated horizons, which are local exponents "
             "over [4, 32] or [8, 64] runs, not asymptotic ones."),
    "pf_reject": ("A P-F rejection is evidence against the informed-flow regime and nothing more. It does not "
                  "indicate size-dependent impact, because a plain linear propagator also produces it."),
    "no_reject": ("No non-rejection is reported as evidence against the paper. The measured power at reachable "
                  "horizons does not support that reading."),
}


# ---- data gate ---------------------------------------------------------------------------------------------


def check_request(day: str, split: str, cells, symbols=None, universe: Path = E.UNIVERSE) -> dict[str, list[str]]:
    """{cell: sorted symbols of the split} for the request; ValueError before any data file is opened."""
    if split not in SPLITS:
        raise ValueError(f"split must be one of {SPLITS}, got {split!r}")
    days, u = E.read_universe(universe)
    if day not in days:
        raise ValueError(f"{day!r} is not a study day")
    bad_cells = sorted(set(cells) - set(NASDAQ_CELLS))
    if bad_cells:
        raise ValueError(f"not NASDAQ-listed cells of the plan: {bad_cells}")
    want = {c: sorted(s for s, (cc, sp) in u.items() if cc == c and sp == split) for c in cells}
    if symbols is None:
        return want
    allowed = {s for v in want.values() for s in v}
    bad = sorted(set(symbols) - allowed)
    if bad:
        why = {s: (u[s] if s in u else "not in a cell") for s in bad[:5]}
        raise ValueError(f"not {split} stocks of {sorted(cells)}: {why}")
    keep = set(symbols)
    return {c: [s for s in v if s in keep] for c, v in want.items()}


def load_runs(day: str, split: str, cells, symbols=None, data: Path | None = None,
              universe: Path = E.UNIVERSE) -> dict[str, np.ndarray]:
    """{symbol: runs} of one study day for the split's stocks of the cells (streamed as explore does)."""
    want = check_request(day, split, cells, symbols, universe)
    flat = sorted(s for v in want.values() for s in v)
    data = E.data_dir() if data is None else Path(data)
    sym_loc = E._read_symbols(data / "tape" / f"{day}.symbols.tsv")
    keep = np.zeros(1 << 16, dtype=bool)
    loc_sym = {}
    for s in flat:
        if s in sym_loc:
            keep[sym_loc[s]] = True
            loc_sym[sym_loc[s]] = s
    runs = E._stream_filter(data / "runs" / f"{day}.runs.npz", keep)
    if not np.all(keep[runs["locate"]]):
        raise RuntimeError("runs of a stock that was not requested")
    out = {}
    if len(runs):
        loc = runs["locate"]
        cut = np.flatnonzero(np.diff(loc.astype(np.int64))) + 1
        for part in np.split(runs, cut):
            out[loc_sym[int(part["locate"][0])]] = part
    _, u = E.read_universe(universe)
    wrong = sorted(s for s in out if u.get(s, ("", ""))[1] != split)
    if wrong:  # second line of defence
        raise RuntimeError(f"stocks outside the {split} split reached the output: {wrong[:5]}")
    for s in out.values():
        if np.any(np.diff(s["ts"].astype(np.int64)) < 0):
            raise ValueError("runs of a stock-day out of time order")
    return out


def to_flow(stock_days: list[list[np.ndarray]]) -> P.Flow:
    """Flow of stocks x study days of runs (stock-major), with the raw sizes and the pre-run log mid."""
    n_days = len(stock_days[0])
    if any(len(row) != n_days for row in stock_days):
        raise ValueError("every stock needs every study day")
    parts = [r for row in stock_days for r in row]
    lengths = np.array([len(r) for r in parts], dtype=np.int64)
    if np.any(lengths == 0):
        raise ValueError("a stock-day without runs")
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int64)
    cat = np.concatenate(parts)
    eps = cat["sign"].astype(float)
    q = cat["q"].astype(float)
    if not np.all(np.abs(eps) == 1):
        raise ValueError("a run without a sign of +1 or -1")
    if not np.all(q > 0):
        raise ValueError("a run with a non-positive size")
    ok = cat["mid_valid"].astype(bool)
    mid = cat["mid2"].astype(np.float64) / 20_000.0  # mid2 = bid + ask in 1e-4 dollars
    lm = np.where(ok, np.log(np.where(ok, mid, 1.0)), np.nan)
    return P.Flow(eps, q, lm, ok, cat["ts"].astype(np.int64), starts, lengths, n_days)


def build_flows(split: str, cells=NASDAQ_CELLS, data: Path | None = None,
                universe: Path = E.UNIVERSE) -> tuple[dict, dict]:
    """({cell: Flow}, {cell: symbols}) over the eight study days for the split's stocks of the cells."""
    days, _ = E.read_universe(universe)
    want = check_request(days[0], split, cells, universe=universe)
    per = {d: load_runs(d, split, cells, data=data, universe=universe) for d in days}
    flows, syms = {}, {}
    for c, ss in want.items():
        if not ss:
            continue
        missing = [(s, d) for s in ss for d in days if s not in per[d]]
        if missing:
            raise ValueError(f"{c}: stocks without runs on a study day: {missing[:5]}")
        flows[c] = to_flow([[per[d][s] for d in days] for s in ss])
        syms[c] = ss
    return flows, syms


# ---- the run -----------------------------------------------------------------------------------------------


def methods(r: P.Reps, null: float, alternative: str) -> dict:
    """Each of the three methods' p-value and 0.90 interval, from the frozen helpers, and the frozen p-value
    (the largest) for comparison."""
    est = np.asarray(r.est, dtype=float)
    se, n = P._jk_se(r.jk)
    c = INF.t_ppf((1 + C.LEVEL) / 2, n - 1)
    ds, ws = P._pct(r.ds, C.LEVEL), P._pct(r.ws, C.LEVEL)
    return {"estimate": est, "null": null, "alternative": alternative,
            "day_stock_bootstrap": {"p": P._boot_p(r.ds, null, alternative), "lo": ds[0], "hi": ds[1]},
            "day_jackknife_t": {"p": P._jk_p(est, r.jk, null, alternative), "lo": est - c * se, "hi": est + c * se},
            "whole_stock_bootstrap": {"p": P._boot_p(r.ws, null, alternative), "lo": ws[0], "hi": ws[1]},
            "p_frozen": P.p_value(r, null, alternative)}


@contextmanager
def observe():
    """Record, for every call of the frozen decision functions inside run_primary, the unit and each method's
    p-values and intervals. The frozen functions run unchanged and their results are returned as they are."""
    names = ("run_unit", "decide_pv", "decide_pr", "decide_pm", "decide_pf")
    saved = {k: getattr(P, k) for k in names}
    rec, ctx = [], {"unit": None}
    lo_m, hi_m = C.PR_CENTRE - C.PR_MARGIN, C.PR_CENTRE + C.PR_MARGIN

    def note(test, d, m, **extra):
        rec.append({"unit": ctx["unit"], "test": test, "p": d["p"], "methods": m, **extra})
        return d

    def run_unit(u, *a, **kw):
        ctx["unit"] = u.name
        return saved["run_unit"](u, *a, **kw)

    def decide_pv(errors, s1, s2):
        return note("P-V", saved["decide_pv"](errors, s1, s2),
                    {"slope_S1": methods(s1, 0.0, "less"), "slope_S2": methods(s2, 0.0, "less")})

    def decide_pr(ends, cens=None):
        return note("P-R", saved["decide_pr"](ends, cens),
                    {"a_star_vs_0.35": methods(ends, lo_m, "greater"), "a_star_vs_0.65": methods(ends, hi_m, "less")},
                    censored_T1_T2=np.zeros(2) if cens is None else np.asarray(cens))

    def decide_pm(slope, top):
        return note("P-M", saved["decide_pm"](slope, top),
                    {"slope": methods(slope, 0.0, "less"), "top": methods(top, 0.0, "less")})

    def decide_pf(pf):
        return note("P-F", saved["decide_pf"](pf), {"slope": methods(pf, 0.0, "two-sided")})

    for k, f in (("run_unit", run_unit), ("decide_pv", decide_pv), ("decide_pr", decide_pr),
                 ("decide_pm", decide_pm), ("decide_pf", decide_pf)):
        setattr(P, k, f)
    try:
        yield rec
    finally:
        for k, f in saved.items():
            setattr(P, k, f)


def attribute(res: dict, rec: list) -> list:
    """Attach the range to each observed decision (the range tests decide once per range, in prereg3d.RANGES
    order; P-M once) and check its p against the TestResult's. Raises on any mismatch."""
    seen: dict = {}
    for r in rec:
        k = seen[(r["unit"], r["test"])] = seen.get((r["unit"], r["test"]), -1) + 1
        t = res["by_tick_bin"][r["unit"]]["tests"][r["test"]]
        if r["test"] == "P-M":
            if k != 0 or r["p"] != t.p:
                raise RuntimeError(f"P-M decision record does not match {r['unit']}")
            r["range"] = None
        else:
            r["range"] = P.RANGE_NAMES[k]
            if r["p"] != t.by_range[r["range"]]["p"]:
                raise RuntimeError(f"{r['test']} {r['unit']} {r['range']}: decision record does not match")
    return rec


def verify_family(res: dict) -> None:
    """The family is exactly the evaluable (test, tick bin) pairs of Q-T1 and Q-T3; Q-T2 never in it and its
    rejected and passed None; Holm recomputed from the family's p-values gives the same decisions; passed =
    rejected and every condition. Raises RuntimeError otherwise."""
    allowed = {(t, tb) for t in C.PRIMARY_TESTS for tb in C.FAMILY_TICK_BINS}
    fam = set(res["family"])
    if not fam <= allowed or set(res["holm"]) != fam:
        raise RuntimeError(f"family holds pairs outside Q-T1 and Q-T3: {sorted(set(res['holm']) - allowed)}")
    ev = res["evaluability"]
    expect = {(t, tb) for tb in C.FAMILY_TICK_BINS if tb in ev for t in C.PRIMARY_TESTS if ev[tb][t]["evaluable"]}
    if fam != expect:
        raise RuntimeError(f"family {sorted(fam)} is not the evaluable pairs {sorted(expect)}")
    for tb in C.DESCRIPTIVE_TICK_BINS:
        for t in res["by_tick_bin"].get(tb, {}).get("tests", {}).values():
            if t.rejected is not None or t.passed is not None:
                raise RuntimeError(f"{tb} {t.name} carries a Holm decision")
    ps = {k: res["by_tick_bin"][k[1]]["tests"][k[0]].p for k in fam}
    h = P.holm(ps)
    for k in fam:
        t = res["by_tick_bin"][k[1]]["tests"][k[0]]
        if res["holm"][k]["p"] != ps[k] or not (h[k]["rejected"] == res["holm"][k]["rejected"] == t.rejected):
            raise RuntimeError(f"Holm decision of {k} does not reproduce")
        if t.passed != bool(t.rejected and t.all_conditions):
            raise RuntimeError(f"passed of {k} is not rejected and every condition")


def run_family(flows: dict, seed: int = SEED, pm: P.PMParams = FROZEN_PM, **kw) -> dict:
    """The primary analysis: run_primary once with the frozen design and the recorded seed."""
    if pm != FROZEN_PM:
        raise ValueError(f"P-M design {pm} is not the frozen {FROZEN_PM}")
    rng = np.random.default_rng(seed)
    with observe() as rec:
        res = P.run_primary(flows, pm, rng, **kw)
    verify_family(res)
    if res_pm(res) != FROZEN_PM:
        raise RuntimeError(f"units ran with P-M design {res_pm(res)}, not the frozen {FROZEN_PM}")
    return {"seed": seed, "pm": res_pm(res),"n_boot": kw.get("n_boot", C.N_BOOT), "n_sur": kw.get("n_sur", P.N_SUR),
            "result": res, "methods": attribute(res, rec)}


def res_pm(res: dict) -> P.PMParams:
    """The P-M design every unit was actually run with (one design, or raise)."""
    pms = {u["unit"]["pm_params"] for u in res["by_tick_bin"].values()}
    if len(pms) != 1:
        raise RuntimeError(f"units ran with different P-M designs: {pms}")
    return pms.pop()


# ---- integrity and records -----------------------------------------------------------------------------------


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True, text=True).stdout.strip()


def preflight() -> dict:
    """HEAD descends from the tag; no file of the tag under tools/ is changed, deleted or renamed; the only
    files under tools/ not in the tag are NEW_FILES; manifest and universe hashes as registered."""
    head = git("rev-parse", "HEAD")
    if git("rev-parse", f"{FROZEN_TAG}^{{commit}}") != FROZEN_COMMIT:
        raise RuntimeError(f"{FROZEN_TAG} is not {FROZEN_COMMIT}")
    if subprocess.run(["git", "-C", str(REPO), "merge-base", "--is-ancestor", FROZEN_TAG, "HEAD"]).returncode:
        raise RuntimeError("HEAD does not descend from the frozen tag")
    changed = git("diff", "--name-only", "--diff-filter=MDRTCU", FROZEN_TAG, "--", "tools/")
    if changed:
        raise RuntimeError(f"frozen files changed: {changed.split()}")
    added = set(git("diff", "--name-only", "--diff-filter=A", FROZEN_TAG, "--", "tools/").split())
    added |= set(git("ls-files", "--others", "--exclude-standard", "--", "tools/").split())
    if not added <= set(NEW_FILES):
        raise RuntimeError(f"files under tools/ beyond the new driver: {sorted(added - set(NEW_FILES))}")
    shas = {"data/manifest.tsv": sha256(REPO / "data" / "manifest.tsv"),
            "data/universe.tsv": sha256(REPO / "data" / "universe.tsv")}
    if shas != {"data/manifest.tsv": MANIFEST_SHA256, "data/universe.tsv": UNIVERSE_SHA256}:
        raise RuntimeError(f"data hashes differ from the registration: {shas}")
    return {"git_head": head, "frozen_tag": FROZEN_TAG, "frozen_commit": FROZEN_COMMIT, "osf": OSF,
            "plan_pdf_sha256": PLAN_PDF_SHA256, "git_status": git("status", "--porcelain"), "sha256": shas,
            "phase4_py_sha256": sha256(Path(__file__)), "added_files": sorted(added),
            "python": sys.version.split()[0], "numpy": np.__version__}


def claim_once(out: Path) -> Path:
    """The confirmatory run happens once: an exclusive marker; a second run (or a rerun after a failure)
    is refused until the marker is removed by a documented decision."""
    out.mkdir(parents=True, exist_ok=True)
    m = out / STARTED
    try:
        fd = os.open(m, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        raise RuntimeError(f"{m} exists: the confirmatory run was already started; not running again") from None
    with os.fdopen(fd, "w") as f:
        f.write(f"pid {os.getpid()} started {time.strftime('%Y-%m-%d %H:%M:%S %z')}\n")
    return m


def key(k) -> str:
    if isinstance(k, tuple):
        if len(k) == 2 and all(isinstance(x, (int, np.integer)) for x in k):
            return f"[{k[0]},{k[1]}]"
        return "|".join(str(x) for x in k)
    return str(k)


def jsonable(x):
    if dataclasses.is_dataclass(x) and not isinstance(x, type):
        return jsonable(dataclasses.asdict(x))
    if isinstance(x, dict):
        return {key(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        f = float(x)
        return f if np.isfinite(f) else (None if np.isnan(f) else ("inf" if f > 0 else "-inf"))
    return x


def stock_counts(flows: dict, syms: dict) -> dict:
    out = {}
    for c, f in flows.items():
        role = "confirmatory" if c in C.CONFIRMATORY and C.tick_bin(c) in C.FAMILY_TICK_BINS else "descriptive"
        out[c] = {"stocks": f.n_stocks, "role": role, "runs": int(f.lengths.sum()),
                  "stock_days": len(f.lengths), "invalid_mids": int((~f.ok).sum()),
                  "N_s_min": int(f.Nsd.min()), "N_s_median": float(np.median(f.Nsd.min(axis=1))),
                  "symbols": syms[c]}
    return out


def plumbing(rec: dict) -> dict:
    """What the dry run may report: that it ran, the family size, which pairs were evaluable and the range
    each used. No estimate and no p-value."""
    res = rec["result"]
    pairs = []
    for tb, ev in res["evaluability"].items():
        for t, e in ev.items():
            used = None
            if tb in res["by_tick_bin"] and t in res["by_tick_bin"][tb]["tests"]:
                tr = res["by_tick_bin"][tb]["tests"][t]
                used = key(tr.range_used) if tr.range_used else tr.range_rule
            pairs.append({"tick_bin": tb, "test": t, "evaluable": e["evaluable"], "in_family": (t, tb) in res["holm"],
                          "range_used": used, "reason": e["reason"]})
    return {"family_size": len(res["family"]), "pairs": pairs,
            "excluded_cells": {tb: u["unit"]["excluded_cells"] for tb, u in res["by_tick_bin"].items()}}


# ---- tables and the results document ----------------------------------------------------------------------


def f4(x) -> str:
    if x is None:
        return "n/a"
    if isinstance(x, str):
        return x
    x = float(x)
    if not np.isfinite(x):
        return "n/a"
    return f"{x:.2e}" if 0 < abs(x) < 1e-3 else f"{x:.4f}"


def interval_text(s: dict, i: int | None = None) -> str:
    def pick(v):
        if i is not None:
            return v[i]
        return v[0] if isinstance(v, list) and len(v) == 1 else v

    est, lo, hi, m = (pick(s[k]) for k in ("estimate", "lo", "hi", "method"))
    return f"{f4(est)} [{f4(lo)}, {f4(hi)}] ({m or 'undefined'})"


def family_rows(res: dict) -> list[dict]:
    """One row per member of the Holm family (Q-T1 and Q-T3 only), in Holm rank order."""
    rows = []
    for tb in C.FAMILY_TICK_BINS:
        for t in C.PRIMARY_TESTS:
            if (t, tb) not in res["holm"]:
                continue
            h, tr = res["holm"][(t, tb)], res["by_tick_bin"][tb]["tests"][t]
            rows.append({"test": t, "tick_bin": tb, "range": key(tr.range_used) if tr.range_used else tr.range_rule,
                         "p": h["p"], "rank": h["rank"], "threshold": h["threshold"], "p_adj": h["p_adj"],
                         "rejected": h["rejected"], "conditions": tr.conditions, "passed": tr.passed})
    return sorted(rows, key=lambda r: r["rank"])


def interpretation(t: str, tr) -> str:
    if tr.passed:
        return INTERPRETATION["pf_reject"] if t == "P-F" else INTERPRETATION["pass"]
    if not tr.rejected:
        return INTERPRETATION["no_reject"]
    failed = [k for k, v in tr.conditions.items() if not v]
    return (f"Holm rejected but not passed (failed: {', '.join(failed)}). The plan fixes no sentence for this "
            "case; it is not a pass.")


def method_rows(rec: dict, tb: str, t: str, rng_name: str | None) -> list[str]:
    out = []
    for r in rec["methods"]:
        if r["unit"] != tb or r["test"] != t or r["range"] != rng_name:
            continue
        for q, m in r["methods"].items():
            n = np.size(m["estimate"])
            for i in range(n):
                def g(v, i=i):
                    return np.ravel(v)[i]
                lab = q if n == 1 else f"{q} ({['T1', 'T2'][i] if t == 'P-R' else i})"
                cells = [f"{f4(g(m[k]['p']))} / [{f4(g(m[k]['lo']))}, {f4(g(m[k]['hi']))}]"
                         for k in ("day_stock_bootstrap", "day_jackknife_t", "whole_stock_bootstrap")]
                out.append(f"| {tb} | {t} | {lab} | {f4(g(m['estimate']))} | {m['alternative']} {m['null']:g} | "
                           + " | ".join(cells) + f" | {f4(g(m['p_frozen']))} |")
    return out


def statistic_lines(t: str, tr) -> list[str]:
    s = tr.statistic
    if not s:
        return [f"- statistic undefined: {tr.range_rule}"]
    if t == "P-V":
        e = s["loo_errors"]
        return [f"- (i) leave-one-day-out weighted error: M1 {f4(e['M1'])}, M2 {f4(e['M2'])}, M3 {f4(e['M3'])} "
                f"(days used {e['days_used']}); M2 fit on all days: {jsonable(e.get('M2_fit_all_days'))}",
                f"- (ii) slope of zeta_V real minus S1: {interval_text(s['slope_S1'])}; p_S1 {f4(s['p_S1'])}",
                f"- (ii) slope of zeta_V real minus S2: {interval_text(s['slope_S2'])}; p_S2 {f4(s['p_S2'])}"]
    if t == "P-R":
        return [f"- a*(T1): {interval_text(s['a_star_T1_T2'], 0)}; censored: {s['censored'][0]}",
                f"- a*(T2): {interval_text(s['a_star_T1_T2'], 1)}; censored: {s['censored'][1]}",
                f"- TOST p at T1, T2: {f4(s['p_end'][0])}, {f4(s['p_end'][1])} (margin 0.35 to 0.65)"]
    if t == "P-M":
        b = s["beta_by_bin"]
        return ["- beta-hat by size bin (small to large): "
                + "; ".join(interval_text(b, i) for i in range(len(b["estimate"]))),
                f"- slope of beta-hat on log q: {interval_text(s['slope'])}; p_slope {f4(s['p_slope'])}",
                f"- top bin beta-hat: {interval_text(s['top'])}; p_top {f4(s['p_top'])}",
                f"- stocks used {len(tr.stocks_used['stocks'])}, stock-days used {tr.stocks_used['stock_days']}"]
    return [f"- slope of zeta_K(a) on a over [0, 1]: {interval_text(s['pf'])}",
            "- zeta_K(a) levels, a = 0, 0.25, 0.5, 0.75, 1 (reported, not tested): "
            + "; ".join(interval_text(s["zeta_K_levels"], i) for i in range(len(s["zeta_K_levels"]["estimate"])))]


def tables_md(rec: dict) -> str:
    res, L = rec["result"], []
    L += ["# Phase 4 confirmatory run: tables", "",
          f"Seed {rec['seed']}; P-M design {dataclasses.asdict(rec['pm'])}; B = {rec['n_boot']}, surrogate draws "
          f"{rec['n_sur']}; git HEAD {rec['preflight']['git_head']}; wall time {rec['wall_s']:.0f} s.", "",
          "## Holm family (Q-T1 and Q-T3, family-wise 0.05)", "",
          f"Family size m = {len(res['family'])}.", "",
          "| rank | test | tick bin | range used | p | Holm threshold | adjusted p | Holm rejects | sub-conditions | "
          "passed |", "|---|---|---|---|---|---|---|---|---|---|"]
    for r in family_rows(res):
        conds = ", ".join(f"{k} {'yes' if v else 'no'}" for k, v in r["conditions"].items())
        L.append(f"| {r['rank']} | {r['test']} | {r['tick_bin']} | {r['range']} | {f4(r['p'])} | {f4(r['threshold'])} "
                 f"| {f4(r['p_adj'])} | {'yes' if r['rejected'] else 'no'} | {conds} | "
                 f"{'yes' if r['passed'] else 'no'} |")
    L += ["", "## Pairs that left the family (reach rule, before any p-value)", ""]
    L += [f"- {x['test']} {x['tick_bin']}: {x['reason']}" for x in res["left"]] or ["- none"]
    L += ["", "## Each method's p-value and 0.90 interval (range used)", "",
          "Columns: p / [lo, hi] per method; the frozen p is the largest of the three. For P-R the frozen TOST p",
          "at an end is the larger of its two rows, and 1 at a censored end (censoring is not applied in these rows).",
          "",
          "| tick bin | test | quantity | estimate | null | day-and-stock bootstrap | day jackknife t | "
          "whole-stock bootstrap | frozen p |", "|---|---|---|---|---|---|---|---|---|"]
    for tb in C.TICK_BINS:
        for t, tr in res.get("by_tick_bin", {}).get(tb, {}).get("tests", {}).items():
            rn = None if t == "P-M" else (key(tr.range_used) if tr.range_used else None)
            if t == "P-M" or rn:
                L += method_rows(rec, tb, t, rn)
    L += ["", "## Evaluability log", ""]
    for tb, ev in res["evaluability"].items():
        for t, e in ev.items():
            L.append(f"- {tb} {t}: {'evaluable' if e['evaluable'] else 'not evaluable'}; {e['reason']}")
    L += ["", "## Stock counts", "", "| cell | role | stocks | stock-days | runs | invalid mids | min N_s |",
          "|---|---|---|---|---|---|---|"]
    for c, s in rec["stock_counts"].items():
        L.append(f"| {c} | {s['role']} | {s['stocks']} | {s['stock_days']} | {s['runs']} | {s['invalid_mids']} | "
                 f"{s['N_s_min']} |")
    return "\n".join(L) + "\n"


def doc_md(rec: dict) -> str:
    res, pre = rec["result"], rec["preflight"]
    L = ["# Phase 4: the confirmatory run", "",
         f"Run once on {rec['started']} with `tools/est/phase4.py confirm`. The plan is frozen at tag {FROZEN_TAG} "
         f"(commit {FROZEN_COMMIT}) and registered on OSF ({OSF}; registered PDF SHA-256 {PLAN_PDF_SHA256}).", "",
         f"- git HEAD at the run: {pre['git_head']}; files added after the freeze: {', '.join(pre['added_files'])}",
         f"- data/manifest.tsv {pre['sha256']['data/manifest.tsv']}; data/universe.tsv "
         f"{pre['sha256']['data/universe.tsv']} (both as registered)",
         f"- bootstrap seed {rec['seed']} (numpy default_rng, one generator through run_primary); B = {rec['n_boot']}; "
         f"surrogate draws {rec['n_sur']}", f"- P-M design: {dataclasses.asdict(rec['pm'])}",
         f"- wall time {rec['wall_s']:.0f} s (loading {rec['load_s']:.0f} s)",
         "- outputs: data/phase4/results.json (every output), data/phase4/tables.md", "",
         "Labels: **Confirmatory** = a member of the one Holm family; **Descriptive** = computed and reported, "
         "never in Holm. Secondary and robustness analyses are not part of this step.", "",
         "## The family", "",
         f"One Holm-Bonferroni family at family-wise 0.05 over the evaluable (test, tick bin) pairs of Q-T1 and "
         f"Q-T3: m = {len(res['family'])}.", ""]
    L += [f"- left the family: {x['test']} {x['tick_bin']} ({x['reason']})" for x in res["left"]]
    L += ["", "| rank | test | tick bin | range used | p | threshold | adjusted p | Holm rejects | passed |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in family_rows(res):
        L.append(f"| {r['rank']} | {r['test']} | {r['tick_bin']} | {r['range']} | {f4(r['p'])} | {f4(r['threshold'])} "
                 f"| {f4(r['p_adj'])} | {'yes' if r['rejected'] else 'no'} | {'yes' if r['passed'] else 'no'} |")
    L += ["", "## Each test (confirmatory)", ""]
    for r in sorted(family_rows(res), key=lambda r: (r["tick_bin"], C.PRIMARY_TESTS.index(r["test"]))):
        tr = res["by_tick_bin"][r["tick_bin"]]["tests"][r["test"]]
        L += [f"### {r['test']}, {r['tick_bin']}", "",
              f"- range: {tr.range_rule}; stocks at T2 of each range: {jsonable(tr.reach)}",
              f"- p = {f4(tr.p)} ({tr.alternative}); Holm: rank {r['rank']}, threshold {f4(r['threshold'])}, "
              f"{'rejected' if r['rejected'] else 'not rejected'}"]
        L += [f"- sub-condition {k}: {'holds' if v else 'fails'}" for k, v in tr.conditions.items()]
        L += statistic_lines(r["test"], tr)
        L += [f"- passed: {'yes' if tr.passed else 'no'}", "",
              f"Pre-registered interpretation: \"{interpretation(r['test'], tr)}\"", ""]
    L += ["## Descriptive (never in Holm)", "",
          "Q-T2 (its one cell, Q-T2S3, at its ceiling of 1/5) is computed by run_primary and is descriptive "
          "(decision of 2 Nov 2026). p-values are before Holm.", ""]
    q2 = res["by_tick_bin"].get("Q-T2", {}).get("tests", {})
    for t in C.PRIMARY_TESTS:
        if t not in q2:
            e = res["evaluability"].get("Q-T2", {}).get(t, {})
            L += [f"### {t}, Q-T2 (descriptive)", "", f"- not computed: {e.get('reason', 'no unit')}", ""]
            continue
        tr = q2[t]
        L += [f"### {t}, Q-T2 (descriptive)", "", f"- range: {tr.range_rule}",
              f"- p before Holm = {f4(tr.p)} ({tr.alternative})"]
        L += [f"- sub-condition {k}: {'holds' if v else 'fails'}" for k, v in tr.conditions.items()]
        L += statistic_lines(t, tr) + [""]
    L += ["### P-K (descriptive): vk_diff = slope zeta_V - 2 slope zeta_K, by tick bin and range", ""]
    for tb, u in res["by_tick_bin"].items():
        for rn, s in u["descriptive"].get("P-K", {}).items():
            L.append(f"- {tb} {rn}: {interval_text(s)}")
    L += ["", "The descriptive cells Q-T1S3, Q-T2S1 and Q-T2S2 were loaded and dropped from the tick-bin units "
          "(logged by tick_bin_unit). Their per-cell results are a secondary analysis in the frozen plan (page 7) "
          "and are not computed here.", "", "## Stock counts", "",
          "| cell | role | stocks | stock-days | runs | invalid mids |", "|---|---|---|---|---|---|"]
    for c, s in rec["stock_counts"].items():
        L.append(f"| {c} | {s['role']} | {s['stocks']} | {s['stock_days']} | {s['runs']} | {s['invalid_mids']} |")
    L += ["", "Full tables, each method's p-value and interval, and the evaluability log: data/phase4/tables.md.", ""]
    return "\n".join(L)


# ---- main ----------------------------------------------------------------------------------------------------


def runs_hashes(data: Path) -> dict:
    days, _ = E.read_universe()
    return {d: sha256(data / "runs" / f"{d}.runs.npz") for d in days}


def main(argv: list[str]) -> int:
    if len(argv) != 2 or argv[1] not in ("dryrun", "confirm"):
        print(__doc__)
        return 2
    mode = argv[1]
    split = "confirm" if mode == "confirm" else "explore"
    pre = preflight()
    started = time.strftime("%Y-%m-%d %H:%M:%S %z")
    out = OUT if mode == "confirm" else E.data_dir() / "phase4" / "dryrun_explore"
    if mode == "confirm":
        claim_once(OUT)
    data = E.data_dir()
    t0 = time.time()
    rh = runs_hashes(data)
    flows, syms = build_flows(split, data=data)
    counts = stock_counts(flows, syms)
    if mode == "confirm":
        got = {c: s["stocks"] for c, s in counts.items()}
        if got != PLAN_CONFIRM_COUNTS:
            raise RuntimeError(f"confirmation stock counts {got} differ from the plan {PLAN_CONFIRM_COUNTS}")
    t1 = time.time()
    rec = run_family(flows)
    rec.update(preflight=pre, started=started, mode=mode, split=split, runs_npz_sha256=rh, stock_counts=counts,
               load_s=t1 - t0, wall_s=time.time() - t0)
    out.mkdir(parents=True, exist_ok=True)
    if mode == "dryrun":
        p = {"ran": True, "wall_s": rec["wall_s"], "load_s": rec["load_s"], "seed": rec["seed"],
             "pm": dataclasses.asdict(rec["pm"]), "preflight": pre,
             "stock_counts": {c: {k: v for k, v in s.items() if k != "symbols"} for c, s in counts.items()},
             **plumbing(rec)}
        (out / "plumbing.json").write_text(json.dumps(jsonable(p), indent=1))
        print(json.dumps(jsonable({k: v for k, v in p.items() if k != "preflight"}), indent=1))
        return 0
    (out / "results.json").write_text(json.dumps(jsonable(rec), indent=1))
    (out / "tables.md").write_text(tables_md(rec))
    DOC.write_text(doc_md(rec))
    print(tables_md(rec))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
