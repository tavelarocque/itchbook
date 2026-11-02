"""Tests for phase4.py, the Phase 4 driver: the data gate (one split, study days only), the Flow inputs (signs,
raw sizes, pre-run log mid, validity, day lengths), the run through primary3d.run_primary with the frozen
P-M design and the recorded seed (the observer leaves every result unchanged), the family check, the
one-run marker and the written tables. Every input is synthetic; no data file is read.

    uv run --group dev --group analysis pytest tools/est/test_phase4.py
"""
from __future__ import annotations
import __future__

import inspect
import json
import textwrap
from typing import NamedTuple

import numpy as np
import phase4 as F
import pytest
from test_3d1 import g0_flow

import prereg3d as C
import primary3d as P
import runs as R
import worlds3d2 as W

DAYS = ("01302019", "01302020")
UNIV = [("AAA", "Q-T1S1", "explore"), ("BBB", "Q-T1S1", "confirm"), ("CCC", "Q-T3S1", "confirm"),
        ("DDD", "N-T1S1", "confirm"), ("EEE", "", ""), ("FFF", "Q-T2S3", "explore"), ("GGG", "Q-T1S3", "confirm")]
KW = {"n_boot": 20, "n_sur": 1}


def refused(f, *args, **kw) -> None:
    try:
        f(*args, **kw)
    except ValueError:
        return
    raise AssertionError(f"not refused: {args} {kw}")


def make_data(tmp, days=DAYS + ("12312019",)):
    """Every symbol on every day (locates differ by day); run 1 of BBB on the first day has no valid mid."""
    u = tmp / "universe.tsv"
    u.write_text(f"# study_days: {' '.join(DAYS)}\nsymbol\tlisting\tcell\tsplit_v2\n"
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
            n = 6 + i + k
            r = np.zeros(n, R.RUN_DTYPE)
            r["locate"] = locs[s]
            r["ts"] = 36_000 * 10**9 + np.arange(n) * 1000
            r["q"] = 100 * (i + 1) + np.arange(n)
            r["sign"] = np.where(np.arange(n) % 2, 1, -1)
            r["bid_px"] = 200_000 + 100 * np.arange(n)
            r["ask_px"] = r["bid_px"] + 100
            r["mid2"] = r["bid_px"].astype(np.uint64) + r["ask_px"]
            r["mid_valid"] = 1
            if s == "BBB" and k == 0:
                r["mid_valid"][1], r["mid2"][1] = 0, 0
            parts.append(r)
            rows[(d, s)] = r
        np.savez_compressed(tmp / "runs" / f"{d}.runs.npz", runs=np.concatenate(parts))
    return u, rows


# ---- checks --------------------------------------------------------------------------------------------------


def check_gate_refuses_other_split(tmp):
    """A confirmation request refuses exploration stocks (and NYSE, uncelled, unknown symbols and cells), and
    returns exactly the confirmation stocks' runs; the reverse for exploration."""
    u, rows = make_data(tmp)
    d = DAYS[0]
    for bad in (["AAA"], ["BBB", "AAA"], ["FFF"], ["DDD"], ["EEE"], ["ZZZ"]):
        refused(F.check_request, d, "confirm", F.NASDAQ_CELLS, bad, universe=u)
    refused(F.check_request, d, "confirm", ["N-T1S1"], universe=u)
    refused(F.check_request, d, "both", F.NASDAQ_CELLS, universe=u)
    refused(F.check_request, d, "confirm", ["Q-T3S1"], ["BBB"], universe=u)  # confirm, but another cell
    refused(F.check_request, d, "explore", F.NASDAQ_CELLS, ["BBB"], universe=u)
    got = F.load_runs(d, "confirm", F.NASDAQ_CELLS, data=tmp, universe=u)
    assert set(got) == {"BBB", "CCC", "GGG"}
    for s, r in got.items():
        assert r.tobytes() == rows[(d, s)].tobytes()
    assert set(F.load_runs(d, "explore", F.NASDAQ_CELLS, data=tmp, universe=u)) == {"AAA", "FFF"}


def check_gate_refuses_non_study_days(tmp):
    u, _ = make_data(tmp)
    for d in ("12312019", "01302021", "", "1302019"):
        refused(F.check_request, d, "confirm", F.NASDAQ_CELLS, universe=u)
        refused(F.load_runs, d, "confirm", F.NASDAQ_CELLS, data=tmp, universe=u)


def check_flow_inputs(tmp):
    """Signs, RAW sizes, ts, day lengths in study-day order, lm = log(mid2 / 2) in dollars before each run,
    NaN and ok False where the mid is invalid; the R2 clip is make_unit's."""
    u, rows = make_data(tmp)
    flows, syms = F.build_flows("confirm", cells=("Q-T1S1", "Q-T3S1", "Q-T1S3"), data=tmp, universe=u)
    assert syms == {"Q-T1S1": ["BBB"], "Q-T3S1": ["CCC"], "Q-T1S3": ["GGG"]}
    f = flows["Q-T1S1"]
    r = [rows[(d, "BBB")] for d in DAYS]
    assert f.n_days == 2 and f.lengths.tolist() == [len(x) for x in r] and f.starts.tolist() == [0, len(r[0])]
    cat = np.concatenate(r)
    assert np.array_equal(f.eps, cat["sign"].astype(float)) and np.array_equal(f.q, cat["q"].astype(float))
    assert np.array_equal(f.ts, cat["ts"].astype(np.int64))
    ok = cat["mid_valid"] == 1
    assert f.ok.tolist() == ok.tolist() and not f.ok[1] and np.isnan(f.lm[1])
    assert np.allclose(f.lm[ok], np.log((cat["bid_px"][ok] + cat["ask_px"][ok]) / 2 / 1e4))
    big = P.Flow(f.eps, np.where(np.arange(len(f.q)) == 3, 1e9, f.q), f.lm, f.ok, f.ts, f.starts, f.lengths, 2)
    un = P.make_unit("Q-T1S1", {"Q-T1S1": big}, ["Q-T1S1"])
    assert un.flow.q[3] < 1e9 and np.array_equal(un.flow.q, P.r2_clip(big.q, big.lengths, 2))


def check_flow_needs_every_day(tmp):
    u, _ = make_data(tmp)
    stock = [np.zeros(3, R.RUN_DTYPE)]
    with pytest.raises(ValueError):
        F.to_flow([stock * 2, stock])
    with pytest.raises(ValueError):
        F.to_flow([[np.zeros(0, R.RUN_DTYPE), np.zeros(0, R.RUN_DTYPE)]])


def toy_flows() -> dict:
    """Q-T1 (two confirmatory cells and the descriptive Q-T1S3), Q-T2S3 and Q-T3, long enough for every range
    and for L = 128 in Q-T1 and Q-T3."""
    return {"Q-T1S1": g0_flow(3, 4, 1500, 31), "Q-T1S2": g0_flow(3, 4, 3000, 32), "Q-T1S3": g0_flow(2, 4, 1500, 33),
            "Q-T2S3": g0_flow(5, 4, 300, 34), "Q-T3S1": g0_flow(3, 4, 1500, 35), "Q-T3S2": g0_flow(2, 4, 1500, 36)}


_TOY: dict = {}


def plain_run(seed: int) -> dict:
    """run_primary alone, without the driver (cached: the frozen code is not mutated by any planted bug)."""
    if seed not in _TOY:
        _TOY[seed] = P.run_primary(toy_flows(), P.PMParams(3, 16, 128, 50), np.random.default_rng(seed), **KW)
    return _TOY[seed]


def same_results(a: dict, b: dict) -> bool:
    if set(a["by_tick_bin"]) != set(b["by_tick_bin"]):
        return False
    for tb, u in a["by_tick_bin"].items():
        for t, x in u["tests"].items():
            y = b["by_tick_bin"][tb]["tests"].get(t)
            if y is None or not (x.p == y.p and x.conditions == y.conditions and x.rejected == y.rejected
                                 and x.passed == y.passed and x.range_used == y.range_used):
                return False
    return a["holm"] == b["holm"]


def check_run_family():
    """The frozen design, the recorded seed (the run reproduces from it, and equals run_primary alone, so the
    observer changes nothing), Q-T2 and the descriptive cell outside Holm, every method's p recorded with the
    frozen p their largest."""
    rec = F.run_family(toy_flows(), **KW)
    assert rec["seed"] == F.SEED == 20261102
    assert rec["pm"] == P.PMParams(3, 16, 128, 50)
    res = rec["result"]
    assert {u["unit"]["pm_params"] for u in res["by_tick_bin"].values()} == {P.PMParams(3, 16, 128, 50)}
    assert same_results(res, plain_run(rec["seed"]))
    assert not same_results(res, plain_run(rec["seed"] + 1))  # the comparison can fail
    fam = set(res["family"])
    assert fam == set(res["holm"]) and fam == {(t, tb) for t in C.PRIMARY_TESTS for tb in ("Q-T1", "Q-T3")}
    assert all(tb != "Q-T2" for _, tb in res["holm"])
    for t in res["by_tick_bin"]["Q-T2"]["tests"].values():
        assert t.rejected is None and t.passed is None
    assert res["by_tick_bin"]["Q-T1"]["unit"]["excluded_cells"] == ["Q-T1S3"]
    assert "Q-T1S3" not in res["by_tick_bin"]["Q-T1"]["unit"]["cells"]
    seen = {(r["unit"], r["test"], r["range"]) for r in rec["methods"]}
    for tb, u in res["by_tick_bin"].items():
        for t, x in u["tests"].items():
            for rn in ([None] if t == "P-M" else list(x.by_range)):
                assert (tb, t, rn) in seen
    for r in rec["methods"]:
        for m in r["methods"].values():
            ps = np.stack([np.ravel(m[k]["p"]) for k in ("day_stock_bootstrap", "day_jackknife_t",
                                                          "whole_stock_bootstrap")])
            fin = np.isfinite(np.ravel(m["estimate"])) & np.any(np.isfinite(ps), axis=0)
            largest = np.max(np.where(np.isfinite(ps), ps, -1.0), axis=0)
            assert np.array_equal(largest[fin], np.ravel(m["p_frozen"])[fin])
            for k in ("day_stock_bootstrap", "day_jackknife_t", "whole_stock_bootstrap"):
                assert np.all(np.ravel(m[k]["lo"])[fin] <= np.ravel(m[k]["hi"])[fin])


def check_frozen_design():
    assert F.FROZEN_PM == P.PMParams(3, 16, 128, 50) == W.design_by_key("b3_l16_L128_m50")
    with pytest.raises(ValueError):
        F.run_family(toy_flows(), pm=P.PMParams(4, 16, 128, 50), **KW)


def check_records_and_tables(tmp):
    rec = F.run_family(toy_flows(), **KW)
    counts = {c: {"stocks": f.n_stocks, "role": "x", "runs": int(f.lengths.sum()), "stock_days": len(f.lengths),
                  "invalid_mids": 0, "N_s_min": int(f.Nsd.min()), "N_s_median": 0.0, "symbols": []}
              for c, f in toy_flows().items()}
    rec.update(preflight={"git_head": "h", "added_files": list(F.NEW_FILES),
                          "sha256": {"data/manifest.tsv": "m", "data/universe.tsv": "u"}},
               started="now", stock_counts=counts, load_s=1.0, wall_s=2.0)
    js = json.dumps(F.jsonable(rec), allow_nan=False)
    back = json.loads(js)
    assert back["seed"] == F.SEED and back["pm"] == {"n_bins": 3, "l_max": 16, "L": 128, "min_runs_per_bin": 50}
    assert set(back["result"]["holm"]) == {f"{t}|{tb}" for t in C.PRIMARY_TESTS for tb in ("Q-T1", "Q-T3")}
    tab, doc = F.tables_md(rec), F.doc_md(rec)
    for text in (tab, doc, js):
        assert "—" not in text and "–" not in text
    assert "Q-T2" not in tab.split("## Holm family")[1].split("## Pairs")[0]
    assert sum(INTERP in doc for INTERP in F.INTERPRETATION.values()) >= 1
    n_members = sum(1 for line in doc.split("## Descriptive")[0].splitlines() if line.startswith("### P-"))
    assert n_members == len(rec["result"]["family"])
    p = F.plumbing(rec)
    assert p["family_size"] == 8 and len(p["pairs"]) == 12
    assert not keys_of(F.jsonable(p)) & {"p", "p_adj", "estimate", "lo", "hi", "statistic", "conditions", "rejected"}


def keys_of(x) -> set:
    if isinstance(x, dict):
        return set(x) | set().union(*(keys_of(v) for v in x.values()))
    if isinstance(x, list):
        return set().union(*(keys_of(v) for v in x)) if x else set()
    return set()


def check_interpretation_mapping():
    def tr(name, rejected, conds):
        return P.TestResult(name, 0.01, "", conds, {}, None, "", {}, {}, rejected=rejected,
                            passed=bool(rejected and all(conds.values())))
    assert F.interpretation("P-R", tr("P-R", True, {"a": True})) == F.INTERPRETATION["pass"]
    assert F.interpretation("P-F", tr("P-F", True, {"a": True})) == F.INTERPRETATION["pf_reject"]
    assert F.interpretation("P-M", tr("P-M", False, {"a": False})) == F.INTERPRETATION["no_reject"]
    assert F.interpretation("P-F", tr("P-F", False, {"a": True})) == F.INTERPRETATION["no_reject"]
    assert "not a pass" in F.interpretation("P-V", tr("P-V", True, {"a": False}))


def check_claim_once(tmp):
    F.claim_once(tmp / "out")
    try:
        F.claim_once(tmp / "out")
    except RuntimeError:
        return
    raise AssertionError("a second confirmatory run was not refused")


CHECKS = {f.__name__: f for f in [
    check_gate_refuses_other_split, check_gate_refuses_non_study_days, check_flow_inputs,
    check_flow_needs_every_day, check_run_family, check_frozen_design, check_records_and_tables,
    check_interpretation_mapping, check_claim_once,
]}
NEEDS_TMP = {"check_gate_refuses_other_split", "check_gate_refuses_non_study_days", "check_flow_inputs",
             "check_flow_needs_every_day", "check_records_and_tables", "check_claim_once"}


def call(name, tmp_path_factory):
    f = CHECKS[name]
    return f(tmp_path_factory.mktemp(name)) if name in NEEDS_TMP else f()


@pytest.mark.parametrize("name", list(CHECKS))
def test_phase4(name, tmp_path_factory):
    call(name, tmp_path_factory)


class Bug(NamedTuple):
    id: str
    muts: tuple  # (module, function, old, new), applied together
    caught_by: tuple[str, ...]


QT2_FROZEN = (
    (P, "holm", "k[1] in C.FAMILY_TICK_BINS):", "k[1] in C.TICK_BINS):"),
    (P, "combine", "pairs = [(t, tb) for tb in ev if tb in C.FAMILY_TICK_BINS for t in C.PRIMARY_TESTS]",
     "pairs = [(t, tb) for tb in ev for t in C.PRIMARY_TESTS]"),
    (P, "run_primary", "for tb in out if tb in C.FAMILY_TICK_BINS for t in", "for tb in out for t in"),
)
BUGS = [
    Bug("loader lets an exploration stock in (first guard)",
        ((F, "check_request", "if cc == c and sp == split)", "if cc == c and sp in SPLITS)"),),
        ("check_gate_refuses_other_split",)),
    Bug("loader lets an exploration stock in (both guards)",
        ((F, "check_request", "if cc == c and sp == split)", "if cc == c and sp in SPLITS)"),
         (F, "load_runs", "if wrong:  # second line of defence", "if False:")), ("check_gate_refuses_other_split",)),
    Bug("loader lets an NYSE cell in", ((F, "check_request", "bad_cells = sorted(set(cells) - set(NASDAQ_CELLS))",
                                         "bad_cells = []"),), ("check_gate_refuses_other_split",)),
    Bug("loader reads a non-study day", ((F, "check_request", "if day not in days:", "if False:"),),
        ("check_gate_refuses_non_study_days",)),
    Bug("invalid mids treated as valid", ((F, "to_flow", 'ok = cat["mid_valid"].astype(bool)',
                                           "ok = np.ones(len(cat), dtype=bool)"),), ("check_flow_inputs",)),
    Bug("driver clips the sizes itself (R2 twice)", ((F, "to_flow", 'q = cat["q"].astype(float)',
                                                      'q = np.minimum(cat["q"].astype(float), 106.0)'),),
        ("check_flow_inputs",)),
    Bug("Q-T2 enters Holm (frozen layers widened; the driver's family check)", QT2_FROZEN, ("check_run_family",)),
    Bug("Q-T2 enters Holm and the driver's family check is off",
        (*QT2_FROZEN, (F, "verify_family", "if not fam <= allowed or set(res[\"holm\"]) != fam:", "if False:"),
         (F, "verify_family", "if fam != expect:", "if False:"),
         (F, "verify_family", "if t.rejected is not None or t.passed is not None:", "if False:")),
        ("check_run_family",)),
    Bug("non-frozen P-M design used", ((F, "run_family", "res = P.run_primary(flows, pm, rng, **kw)",
                                        "res = P.run_primary(flows, P.PMParams(4, 16, 128, 50), rng, **kw)"),),
        ("check_run_family",)),
    Bug("non-frozen P-M design used, the driver's design checks off",
        ((F, "run_family", "res = P.run_primary(flows, pm, rng, **kw)",
          "res = P.run_primary(flows, P.PMParams(4, 16, 128, 50), rng, **kw)"),
         (F, "run_family", "if res_pm(res) != FROZEN_PM:", "if False:")), ("check_run_family",)),
    Bug("design argument not checked", ((F, "run_family", "if pm != FROZEN_PM:", "if False:"),),
        ("check_frozen_design",)),
    Bug("seed not recorded", ((F, "run_family", '"seed": seed,', '"seed": None,'),), ("check_run_family",)),
    Bug("seed recorded but not used (unseeded generator)",
        ((F, "run_family", "rng = np.random.default_rng(seed)", "rng = np.random.default_rng()"),),
        ("check_run_family",)),
    Bug("seed recorded but another used", ((F, "run_family", "rng = np.random.default_rng(seed)",
                                            "rng = np.random.default_rng(seed + 1)"),), ("check_run_family",)),
    Bug("observer reports the smallest method's p", ((F, "methods", '"p_frozen": P.p_value(r, null, alternative)',
                                                       '"p_frozen": np.fmin(P._boot_p(r.ds, null, alternative), '
                                                       'P._boot_p(r.ws, null, alternative))'),),
        ("check_run_family",)),
    Bug("second confirmatory run allowed", ((F, "claim_once", "os.O_CREAT | os.O_EXCL | os.O_WRONLY",
                                             "os.O_CREAT | os.O_WRONLY"),), ("check_claim_once",)),
    Bug("a non-rejection read as a pass", ((F, "interpretation", "if tr.passed:", "if True:"),),
        ("check_interpretation_mapping",)),
]


def mutate_all(monkeypatch, muts) -> None:
    """conftest.mutate for several replacements, possibly in one function (each target must occur exactly once
    in that function's source), compiled once per function."""
    by: dict = {}
    for module, name, old, new in muts:
        by.setdefault((module, name), []).append((old, new))
    for (module, name), pairs in by.items():
        src = textwrap.dedent(inspect.getsource(getattr(module, name)))
        for old, new in pairs:
            if src.count(old) != 1:
                raise AssertionError(f"mutation target occurs {src.count(old)} times in {name}: {old!r}")
            src = src.replace(old, new)
        code = compile(src, f"<mutant {module.__name__}.{name}>", "exec",
                       flags=__future__.annotations.compiler_flag, dont_inherit=True)
        ns: dict = {}
        exec(code, module.__dict__, ns)
        monkeypatch.setattr(module, name, ns[name])


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, monkeypatch, tmp_path_factory):
    mutate_all(monkeypatch, bug.muts)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            call(name, tmp_path_factory)
