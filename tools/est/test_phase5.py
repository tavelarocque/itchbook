"""Tests for phase5.py (secondary and robustness analyses, descriptive). Every input is synthetic; no data file is
read and nothing is written under data/phase4/ or data/phase4_repro/.

The anchor (the primary setting through phase5's wrappers reproduces run_primary bit for bit), the data gate
(confirmation stocks only), the cells without a registered ceiling, one setting per arm, the seed, the psi-hat
note, the p-value header, the write guard and the document lint, each with planted bugs that must be caught.

    uv run --group dev --group analysis pytest tools/est/test_phase5.py
"""

from __future__ import annotations

import json
from contextlib import contextmanager
from typing import NamedTuple

import numpy as np
import phase4 as F4
import phase5 as F5
import pytest
from test_3d1 import g0_flow
from test_phase4 import make_data, mutate_all, toy_flows

import pooled as PO
import prereg3d as C
import primary3d as P

KW = {"n_boot": 20, "n_sur": 1}
F5_ORIG = dict(vars(F5))
FROZEN_ORIG = {
    "r2_clip": P.r2_clip,
    "draws": P.draws,
    "run_stats": P.run_stats,
    "run_unit": P.run_unit,
    "panel_moments": PO.panel_moments,
}


def tiny_flows() -> dict:
    """One confirmatory cell (Q-T2S3) and one cell without a ceiling (Q-T1S3): quick runs."""
    return {"Q-T1S3": g0_flow(2, 4, 300, 33), "Q-T2S3": g0_flow(5, 4, 300, 34)}


def stored_of(res: dict) -> dict:
    """What results.json would hold for a run_primary result."""
    return json.loads(json.dumps(F4.jsonable({"result": res}), allow_nan=False))


@contextmanager
def refuses(exc):
    """The body must raise exc; otherwise an AssertionError (not pytest's Failed, a BaseException, so that a planted
    bug's test sees it as caught)."""
    try:
        yield
    except exc:
        return
    raise AssertionError(f"not refused: {exc.__name__} expected")


_PLAIN: dict = {}


def plain(flows_fn, seed: int = F5.SEED) -> dict:
    """run_primary alone, no wrapper (the frozen code is never mutated by a planted bug here)."""
    k = (flows_fn.__name__, seed)
    if k not in _PLAIN:
        _PLAIN[k] = P.run_primary(flows_fn(), F5.FROZEN_PM, np.random.default_rng(seed), **KW)
    return _PLAIN[k]


RENDER_ONLY = {
    "doc_a",
    "doc_b",
    "doc_c",
    "render_doc",
    "lint_doc",
    "write_outputs",
    "safe_path",
    "P_HEADER",
    "FORBIDDEN",
    "LABELS",
    "CHOICES",
    "PLANNED_NOT_RUN",
    "method_cells",
    "f4",
    "ivl",
    "yn",
}


def unmutated() -> bool:
    """Nothing compute() uses is mutated (rendering and writing may be)."""
    return all(vars(F5).get(k) is v for k, v in F5_ORIG.items() if k not in RENDER_ONLY)


def rec_flows() -> dict:
    """As toy_flows, with five stocks in Q-T1S1 and Q-T3S1 so that those cells are evaluable alone."""
    return {
        "Q-T1S1": g0_flow(5, 4, 1500, 31),
        "Q-T1S2": g0_flow(3, 4, 3000, 32),
        "Q-T1S3": g0_flow(2, 4, 1500, 33),
        "Q-T2S3": g0_flow(5, 4, 300, 34),
        "Q-T3S1": g0_flow(5, 4, 1500, 35),
        "Q-T3S2": g0_flow(2, 4, 1500, 36),
    }


_REC: dict = {}


def toy_record() -> dict:
    """compute() on the toy flows with the header fields a run adds, as JSON (cached only when nothing is
    mutated, so a planted bug never reaches the cache)."""
    if "rec" in _REC and unmutated():
        return _REC["rec"]
    rec = F5.compute(rec_flows(), stored_of(plain(rec_flows)), **KW)
    rec.update(
        started="2026-11-02 00:00:00 -0700",
        load_s=1.0,
        wall_s=2.0,
        preflight={
            "git_head": "h" * 40,
            "added_files": list(F5.ALLOWED_ADDED),
            "phase5_py_sha256": "s",
            "plan_pdf": "plan.pdf",
            "plan_pdf_sha256": F4.PLAN_PDF_SHA256,
        },
        runs_npz_sha256={"01302019": "a" * 64},
        stock_counts={c: f.n_stocks for c, f in rec_flows().items()},
        protected={"n_files": 5, "unchanged": True, "files": {}},
    )
    rec = json.loads(json.dumps(F4.jsonable(rec), allow_nan=False))
    if unmutated():
        _REC["rec"] = rec
    return rec


# ---- checks --------------------------------------------------------------------------------------------------


def check_anchor():
    """The primary arm through phase5's wrappers reproduces run_primary bit for bit; the comparison fails for
    another arm and another seed; compute stops at a failed anchor."""
    stored = stored_of(plain(toy_flows))
    prim = F5.run_arm(toy_flows(), F5.PRIMARY, **KW)
    a = F5.anchor(prim["result"], stored["result"])
    assert a["identical"] and a["mismatches"] == 0 and a["fields_compared"] > 500
    assert not F5.anchor(F5.run_arm(toy_flows(), F5.Arm("R3", False), **KW)["result"], stored["result"])["identical"]
    assert not F5.anchor(plain(toy_flows, F5.SEED + 1), stored["result"])["identical"]
    with refuses(F5.AnchorFailed):
        F5.compute(toy_flows(), stored_of(plain(toy_flows, F5.SEED + 1)), **KW)
    for k, f in FROZEN_ORIG.items():  # every wrapper restored
        assert (PO.panel_moments if k == "panel_moments" else getattr(P, k)) is f


def check_clip_rules():
    """clip_sizes('R2') is the frozen r2_clip bit for bit; R3 clips more, none not at all, R1 per stock-day."""
    f = g0_flow(3, 4, 500, 7)
    r2 = F5.clip_sizes("R2", f.q, f.lengths, f.n_days)
    assert r2.tobytes() == P.r2_clip(f.q, f.lengths, f.n_days).tobytes()
    r3, none = F5.clip_sizes("R3", f.q, f.lengths, f.n_days), F5.clip_sizes("none", f.q, f.lengths, f.n_days)
    assert np.all(r3 <= r2) and np.any(r3 < r2) and np.array_equal(none, f.q)
    r1 = F5.clip_sizes("R1", f.q, f.lengths, f.n_days)
    for i in range(len(f.lengths)):
        s = f.day(i)
        assert np.array_equal(r1[s], np.minimum(f.q[s], 0.01 * f.q[s].sum()))


def check_loader_confirm_only(tmp):
    """Only confirmation stocks are loaded; an exploration stock is refused by the second guard."""
    u, _ = make_data(tmp)
    flows, syms = F5.load_confirm(data=tmp, universe=u)
    assert syms == {"Q-T1S1": ["BBB"], "Q-T1S3": ["GGG"], "Q-T3S1": ["CCC"]}
    assert set(flows) == set(syms)
    with refuses(RuntimeError):
        F5.check_confirm_only({"Q-T1S1": ["AAA"]}, u)  # an exploration stock of the cell
    with refuses(RuntimeError):
        F5.check_confirm_only({"Q-T1S1": ["CCC"]}, u)  # a confirmation stock of another cell
    F5.check_confirm_only({"Q-T1S1": ["BBB"]}, u)


def check_cells():
    """Cells without a registered ceiling: not estimable, nothing computed; a cell with one: the frozen
    run_cell's results, only its evaluable tests, each method's interval and P-K recorded."""
    for c in F5.NASDAQ_CELLS:
        assert F5.registered_ceiling(c) == F5.PLAN_CEILINGS[c] == C.CONFIRMATORY.get(c)
    assert all(F5.registered_ceiling(c) is None for c in C.DESCRIPTIVE_CELLS)
    out = F5.run_cells(tiny_flows(), **KW)
    assert list(out) == list(F5.NASDAQ_CELLS)
    assert out["Q-T1S3"] == {"estimable": False, "reason": F5.NOT_ESTIMABLE, "n_stocks": 2}
    for c in ("Q-T1S1", "Q-T2S1", "Q-T3S3"):
        assert out[c]["estimable"] is False and "tests" not in out[c]
    r = out["Q-T2S3"]
    assert r["estimable"] and r["ceiling"] == 1 / 5 and r["seed"] == F5.SEED
    want = tuple(t for t in C.PRIMARY_TESTS if r["evaluability"][t]["evaluable"])
    assert tuple(t for t in C.PRIMARY_TESTS if t in r["tests"]) == want and want
    ref = P.run_cell(tiny_flows(), "Q-T2S3", F5.FROZEN_PM, np.random.default_rng(F5.SEED), tests=want, **KW)
    for t in want:
        assert r["tests"][t]["p"] == ref["tests"][t].p
        assert F4.jsonable(r["tests"][t]["statistic"]) == F4.jsonable(F5.tr_dict(ref["tests"][t])["statistic"])
        assert "rejected" not in r["tests"][t] and "passed" not in r["tests"][t]
    assert {(m["test"], m["range"]) for m in r["methods"]} >= {
        (t, rn) for t in want if t != "P-M" for rn in P.RANGE_NAMES
    }
    assert F4.jsonable(r["P-K"]["summary"]) == F4.jsonable(ref["descriptive"]["P-K"])


def check_arms_one_setting():
    """Four arms, each changing exactly one setting; two at once refused; the settings the wrappers applied are
    the arm's; the frozen functions are restored."""
    assert len(F5.ARMS) == 4 and F5.PRIMARY == F5.Arm("R2", False)
    assert all(len(a.changed()) == 1 for a in F5.ARMS)
    assert {(a.clip, a.centered) for a in F5.ARMS} == {("R1", False), ("R3", False), ("none", False), ("R2", True)}
    for bad in (F5.Arm("R1", True), F5.Arm("none", True), F5.Arm("R7", False)):
        with refuses(ValueError):
            F5.check_arm(bad)
        with refuses(ValueError):
            F5.run_arm(tiny_flows(), bad, **KW)
    for arm in (F5.Arm("R1", False), F5.Arm("R2", True)):
        run = F5.run_arm(tiny_flows(), arm, **KW)
        assert run["effective"] == {"clip": [arm.clip], "centered": [arm.centered]}
    for k, f in FROZEN_ORIG.items():
        assert (PO.panel_moments if k == "panel_moments" else getattr(P, k)) is f


def check_seed():
    """Every arm from SEED: the recorded seed, the first unit's draws equal to the primary's (and for the centred
    arm every unit's), and the primary equal to run_primary at SEED alone."""
    prim, arms = F5.run_arms(tiny_flows(), **KW)
    assert F5.anchor(prim["result"], stored_of(plain(tiny_flows))["result"])["identical"]
    first = next(iter(prim["draws"]))
    for a in arms.values():
        assert a["seed"] == F5.SEED
        assert a["draws"][first] == prim["draws"][first]
    assert arms[F5.Arm("R2", True).name]["draws"] == prim["draws"]
    bad = {**arms[F5.Arm("R1", False).name], "seed": F5.SEED + 1}
    with refuses(RuntimeError):
        F5.check_seeds(prim, {"x": bad})
    other = F5.run_arm(tiny_flows(), F5.Arm("R1", False), seed=F5.SEED + 1, **KW)
    with refuses(RuntimeError):
        F5.check_seeds(prim, {"x": {**other, "seed": F5.SEED}})


def check_psi_note():
    """Every line of the document that mentions psi-hat carries primary3d.PSI_NOTE verbatim; a stored record with
    another note is refused; the lint refuses a psi-hat line without the note."""
    text = F5.render_doc(toy_record())
    lines = [x for x in text.splitlines() if "psi-hat" in x.lower()]
    assert len(lines) >= 4  # three tick bins and the Q-T1 paragraph
    assert all(P.PSI_NOTE in x for x in lines)
    stored = stored_of(plain(toy_flows))
    stored["result"]["by_tick_bin"]["Q-T1"]["tests"]["P-R"]["secondary"]["psi_hat_note"] = "corrected"
    with refuses(RuntimeError):
        F5.section_b(stored)
    with refuses(ValueError):
        F5.lint_doc("**Measurement** psi-hat is 0.05.\n")


def check_p_header():
    """Every p column is headed exactly 'p (descriptive, no multiplicity control)'; the lint refuses others."""
    assert F5.P_HEADER == "p (descriptive, no multiplicity control)"
    lines = F5.render_doc(toy_record()).splitlines()
    heads = [x for i, x in enumerate(lines[:-1]) if x.startswith("| ") and lines[i + 1].startswith("|---")]
    assert sum("| p (descriptive, no multiplicity control) |" in h for h in heads) >= 3
    for h in heads:
        cells = [c.strip() for c in h.strip().strip("|").split("|")]
        assert all(c == "p (descriptive, no multiplicity control)" for c in cells if c == "p" or c.startswith("p ("))
    for bad in ("| cell | p |\n", "| cell | p-value |\n", "| cell | p (descriptive) |\n", "| cell | p before x |\n"):
        with refuses(ValueError):
            F5.lint_doc(bad)
    F5.lint_doc("| cell | p (descriptive, no multiplicity control) |\n")


def check_write_guard(tmp):
    """Nothing under data/phase4/ or data/phase4_repro/, and in the repo only data/phase5/ and the document."""
    for p in (
        F5.REPO / "data" / "phase4" / "results.json",
        F5.REPO / "data" / "phase4" / "RUN_STARTED",
        F5.REPO / "data" / "phase4",
        F5.REPO / "data" / "phase4_repro" / "compare.json",
        F5.REPO / "data" / "manifest.tsv",
        F5.REPO / "tools" / "est" / "x.py",
        F5.REPO / "docs" / "x.md",
    ):
        with refuses(PermissionError):
            F5.safe_path(p)
    for p in (F5.OUT / "results.json", F5.DOC, tmp / "x.json"):
        F5.safe_path(p)
    before = F5.protected_snapshot()
    assert before["files"] and all(k.startswith(("data/phase4/", "data/phase4_repro/")) for k in before["files"])
    F5.write_outputs(toy_record(), out=tmp / "out", doc=tmp / "doc.md")
    assert sorted(p.name for p in (tmp / "out").iterdir()) == ["results.json"] and (tmp / "doc.md").is_file()
    assert F5.protected_snapshot() == before


def check_lint():
    """The rendered document passes; dashes, confirmation or refutation wording, family words and unlabelled
    bullets are refused."""
    F5.lint_doc(F5.render_doc(toy_record()))
    for bad in (
        "The estimate \u2014 here.\n",
        "a range 4\u201332\n",
        "This confirms the paper.\n",
        "The paper is supported.\n",
        "It contradicts Eq. 30.\n",
        "P-F rejected.\n",
        "Holm at 0.05\n",
        "P-R passes.\n",
        "- an unlabelled bullet\n",
    ):
        with refuses(ValueError):
            F5.lint_doc(bad)


def check_end_to_end(tmp):
    """compute, write, re-render from the JSON alone: the same document, with sections A, B, C and the not-run
    list, the anchor line, the cells without a ceiling, every planned item."""
    rec = toy_record()
    text = F5.write_outputs(rec, out=tmp / "out", doc=tmp / "doc.md")
    back = json.loads((tmp / "out" / "results.json").read_text())
    assert F5.render_doc(back) == text == (tmp / "doc.md").read_text()
    for h in (
        "## Anchor",
        "## A. Per-cell",
        "## B. P-R secondaries",
        "## C. Robustness arms",
        "## Planned, not run in Phase 5",
        "## Choices made in Phase 5",
    ):
        assert h in text
    assert "bit for bit" in text and "0 mismatches" in text
    assert text.count(F5.NOT_ESTIMABLE) >= 2
    assert all(x in text for x in F5.PLANNED_NOT_RUN)
    assert set(back["arms"]["arms"]) == {a.name for a in F5.ARMS}
    for a in back["arms"]["arms"].values():
        assert set(a["by_tick_bin"]) == set(F5.ARM_TICK_BINS)
        assert "rejected" not in json.dumps(a) and "holm" not in json.dumps(a)


def check_rratio_horizon():
    """RRATIO_T is pooled.core2's middle horizon of each range."""
    lt2 = np.log2(C.T_GRID.astype(float))
    for (lo, hi), rn in zip(C.RANGES, P.RANGE_NAMES, strict=True):
        t = int(np.flatnonzero(lt2 == np.floor((np.log2(lo) + np.log2(hi)) / 2))[0])
        assert F5.RRATIO_T[rn] == int(C.T_GRID[t])


CHECKS = {
    f.__name__: f
    for f in [
        check_anchor,
        check_clip_rules,
        check_loader_confirm_only,
        check_cells,
        check_arms_one_setting,
        check_seed,
        check_psi_note,
        check_p_header,
        check_write_guard,
        check_lint,
        check_end_to_end,
        check_rratio_horizon,
    ]
}
NEEDS_TMP = {"check_loader_confirm_only", "check_write_guard", "check_end_to_end"}


def call(name, tmp_path_factory):
    f = CHECKS[name]
    return f(tmp_path_factory.mktemp(name)) if name in NEEDS_TMP else f()


@pytest.mark.parametrize("name", list(CHECKS))
def test_phase5(name, tmp_path_factory):
    call(name, tmp_path_factory)


class Bug(NamedTuple):
    id: str
    muts: tuple  # (module, function, old, new), applied together
    caught_by: tuple[str, ...]
    sets: tuple = ()  # (module, attribute, value)


BUGS = [
    # (1) an exploration stock let in
    Bug(
        "1a loader asks for the exploration split",
        ((F5, "load_confirm", 'F4.build_flows("confirm",', 'F4.build_flows("explore",'),),
        ("check_loader_confirm_only",),
    ),
    Bug(
        "1b confirmation guard checks the cell only",
        ((F5, "check_confirm_only", "if u.get(s) != (c, SPLIT))", 'if u.get(s, ("", ""))[0] != c)'),),
        ("check_loader_confirm_only",),
    ),
    Bug(
        "1c loader asks for exploration and the guard is weakened",
        (
            (F5, "load_confirm", 'F4.build_flows("confirm",', 'F4.build_flows("explore",'),
            (F5, "check_confirm_only", "if u.get(s) != (c, SPLIT))", 'if u.get(s, ("", ""))[0] != c)'),
        ),
        ("check_loader_confirm_only",),
    ),
    # (2) a descriptive cell given a made-up ceiling
    Bug(
        "2a made-up ceiling for a cell without one",
        ((F5, "registered_ceiling", "rho = C.CONFIRMATORY.get(cell)", "rho = C.CONFIRMATORY.get(cell, 1 / 5)"),),
        ("check_cells",),
    ),
    Bug(
        "2b made-up ceiling and the plan check off",
        (
            (F5, "registered_ceiling", "rho = C.CONFIRMATORY.get(cell)", "rho = C.CONFIRMATORY.get(cell, 1 / 5)"),
            (F5, "registered_ceiling", "if rho != PLAN_CEILINGS[cell]:", "if False:"),
        ),
        ("check_cells",),
    ),
    Bug(
        "2c cell without a ceiling not skipped", ((F5, "run_cells", "if rho is None:", "if False:"),), ("check_cells",)
    ),
    # (3) an arm that changes two settings at once
    Bug(
        "3a an arm of two settings in ARMS",
        (),
        ("check_arms_one_setting",),
        ((F5, "ARMS", (F5.Arm("R1", True), F5.Arm("R3", False), F5.Arm("none", False), F5.Arm("R2", True))),),
    ),
    Bug(
        "3b two settings at once not refused",
        ((F5, "check_arm", "if len(arm.changed()) > 1:", "if False:"),),
        ("check_arms_one_setting",),
    ),
    Bug(
        "3c clip arms also centre the moments",
        ((F5, "arm_setting", "c = arm.centered", 'c = arm.centered or arm.clip != "R2"'),),
        ("check_arms_one_setting",),
    ),
    Bug(
        "3d the centred arm also clips at R1",
        ((F5, "arm_setting", "rule = arm.clip", 'rule = "R1" if arm.centered else arm.clip'),),
        ("check_arms_one_setting",),
    ),
    # (4) psi-hat reported without PSI_NOTE
    Bug(
        "4a psi-hat line without the note",
        ((F5, "doc_b", "f\"stocks with a slope. Note: {s['psi_hat_note']}\"", 'f"stocks with a slope."'),),
        ("check_psi_note",),
    ),
    Bug(
        "4b psi-hat lines without the note and the lint off",
        (
            (F5, "doc_b", "f\"stocks with a slope. Note: {s['psi_hat_note']}\"", 'f"stocks with a slope."'),
            (F5, "lint_doc", "and P.PSI_NOTE not in line:", "and False:"),
        ),
        ("check_psi_note",),
    ),
    Bug(
        "4c the a* paragraph without the note",
        ((F5, "doc_b", "Note: {q['psi_hat_note']}. These values", "These values"),),
        ("check_psi_note",),
    ),
    Bug(
        "4d the stored note not checked and replaced",
        ((F5, "section_b", 'if sec["psi_hat_note"] != P.PSI_NOTE:', "if False:"),),
        ("check_psi_note",),
    ),
    # (5) a p-value column labelled without "descriptive, no multiplicity control"
    Bug(
        "5a per-cell p column headed 'p'",
        (
            (
                F5,
                "doc_a",
                '[8,64] | {P_HEADER} | sub-conditions |",',
                '[8,64] | p | sub-conditions |",',
            ),
        ),
        ("check_p_header",),
    ),
    Bug("5b the header constant shortened", (), ("check_p_header",), ((F5, "P_HEADER", "p (descriptive)"),)),
    Bug(
        "5c section B p column headed 'p' and the lint off",
        (
            (
                F5,
                "doc_b",
                'f"| tick bin | quantity | range | estimate | widest 0.90 interval (method) | {P_HEADER} | '
                'JSON path |"',
                '"| tick bin | quantity | range | estimate | widest 0.90 interval (method) | p | JSON path |"',
            ),
            (F5, "lint_doc", "and c != P_HEADER:", "and False:"),
        ),
        ("check_p_header",),
    ),
    # (6) a robustness arm run with a different seed
    Bug(
        "6a the arms run from another seed",
        (
            (
                F5,
                "run_arms",
                "arms = {a.name: run_arm(flows, a, **kw) for a in ARMS}",
                "arms = {a.name: run_arm(flows, a, seed=SEED + 1, **kw) for a in ARMS}",
            ),
        ),
        ("check_seed",),
    ),
    Bug(
        "6b run_arm ignores its seed",
        ((F5, "run_arm", "rng = np.random.default_rng(seed)", "rng = np.random.default_rng(seed + 1)"),),
        ("check_seed", "check_anchor"),
    ),
    Bug(
        "6c the centred arm from another seed and the seed check off",
        (
            (
                F5,
                "run_arms",
                "arms = {a.name: run_arm(flows, a, **kw) for a in ARMS}",
                "arms = {a.name: run_arm(flows, a, seed=SEED + a.centered, **kw) for a in ARMS}",
            ),
            (F5, "check_seeds", 'if a["seed"] != SEED:', "if False:"),
            (F5, "check_seeds", 'if a["draws"].get(first) != prim["draws"][first]:', "if False:"),
        ),
        ("check_seed",),
    ),
    Bug(
        "6d an arm unseeded",
        ((F5, "run_arm", "rng = np.random.default_rng(seed)", "rng = np.random.default_rng()"),),
        ("check_seed",),
    ),
    # (7) any write under data/phase4/ (the guard stays in place, so nothing is written there)
    Bug(
        "7a the writer targets data/phase4/",
        (
            (
                F5,
                "write_outputs",
                'safe_path(out / "results.json")',
                'safe_path(REPO / "data" / "phase4" / "phase5_results.json")',
            ),
        ),
        ("check_write_guard", "check_end_to_end"),
    ),
    Bug(
        "7b the write guard off",
        (
            (F5, "safe_path", "if p == q or q in p.parents:", "if False:"),
            (F5, "safe_path", "if (p == repo or repo in p.parents) and not", "if False and not"),
        ),
        ("check_write_guard",),
    ),
    # the anchor itself
    Bug(
        "anchor: R2 path clips at R3",
        (
            (
                F5,
                "clip_sizes",
                "CL.apply_rule(rule, qs, stocks)",
                'CL.apply_rule("R3" if rule == "R2" else rule, qs, stocks)',
            ),
        ),
        ("check_anchor", "check_clip_rules"),
    ),
    Bug("anchor: moments always centred", ((F5, "arm_setting", "c = arm.centered", "c = True"),), ("check_anchor",)),
    # the document lint
    Bug(
        "lint: dashes allowed",
        ((F5, "lint_doc", 'if "\\u2014" in text or "\\u2013" in text:', "if False:"),),
        ("check_lint",),
    ),
    Bug(
        "lint: confirmation wording allowed",
        ((F5, "lint_doc", "m = FORBIDDEN.search(", "m = None and FORBIDDEN.search("),),
        ("check_lint",),
    ),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, monkeypatch, tmp_path_factory):
    mutate_all(monkeypatch, bug.muts)
    for module, attr, value in bug.sets:
        monkeypatch.setattr(module, attr, value)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            call(name, tmp_path_factory)
