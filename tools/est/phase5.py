"""Phase 5: the pre-registered secondary and robustness analyses, run once on the confirmation data.

Secondary and descriptive. Nothing here is confirmatory, nothing here can change any Phase 4 result, and nothing
here is a test of the paper. The plan is the PDF registered on OSF (osf.io/z4t7x, SHA-256 in phase4.PLAN_PDF_SHA256),
page 7, "Secondary, descriptive and robustness analyses". Three items only:

A. Per-cell results for the nine NASDAQ-listed cells (run_cells). A cell with a registered ceiling rho* (prereg3d.
   CONFIRMATORY, as on page 3 of the plan) is one pooling unit (primary3d.cell_unit); the frozen reach rule
   (primary3d.evaluability) decides which tests are computed, and primary3d.run_unit computes them with a fresh
   default_rng(SEED). A cell without a registered ceiling (Q-T1S3, Q-T2S1, Q-T2S2) is reported as NOT_ESTIMABLE
   and nothing is computed for it. Each method's interval comes from phase4.observe (P-V, P-R, P-M, P-F) and from
   the recorded replicates of primary3d.run_stats (P-K); no draw is taken by either recorder.
B. P-R secondaries per tick bin, read from data/phase4/results.json (section_b; nothing is recomputed).
C. Robustness arms (run_arm): clip rule R1, R3 or none instead of R2, or moments centred on each stock-day's
   window means instead of uncentred, one setting at a time. Each arm is the frozen primary3d.run_primary over all
   nine cells with FROZEN_PM and default_rng(SEED), inside arm_setting, which (i) replaces primary3d.r2_clip by
   clip_sizes, i.e. clip3b.apply_rule(rule) over the same stock-days, and (ii) calls pooled.panel_moments with
   centered=arm.centered (the frozen flag). Nothing else changes. The primary arm (R2, uncentred) runs through
   the same wrappers and must reproduce data/phase4/results.json["result"] bit for bit (anchor, floats by
   float.hex, phase4_repro.Compare) before anything else is computed or reported.

Guards: the confirmation split only (load_confirm, check_confirm_only); the runs files' SHA-256 equal to the
ones stored in data/phase4/results.json; no write under data/phase4/, data/phase4_repro/ or anywhere in the repo
but data/phase5/ and docs/phase5_secondary.md (safe_path); one run (phase4.claim_once on data/phase5/); every
arm changes at most one setting (check_arm) and the settings the wrappers actually applied are checked; every
arm uses SEED (check_seeds: the first unit's bootstrap draws, taken before any surrogate draw, are identical
across arms); the document is linted before it is written (lint_doc).

    uv run --group analysis python tools/est/phase5.py run      # the one run: compute, write JSON and document
    uv run --group analysis python tools/est/phase5.py render   # re-render the document from the saved JSON only
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import phase4 as F4
from phase4_repro import Compare

import clip3b as CL
import explore as E
import pooled as PO
import prereg3d as C
import primary3d as P

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "data" / "phase5"
DOC = REPO / "docs" / "phase5_secondary.md"
PHASE4_RESULTS = REPO / "data" / "phase4" / "results.json"
PROTECTED = (REPO / "data" / "phase4", REPO / "data" / "phase4_repro")
PLAN_PDF = Path(
    os.environ.get("ITCHBOOK_PLAN_PDF", "~/Downloads/Frozen_Preregistration_itchbook_2026-11-02.pdf")
).expanduser()
ALLOWED_ADDED = (
    "tools/est/phase4.py",
    "tools/est/test_phase4.py",
    "tools/est/phase4_repro.py",
    "tools/est/phase5.py",
    "tools/est/test_phase5.py",
)

SEED = F4.SEED
FROZEN_PM = F4.FROZEN_PM
NASDAQ_CELLS = F4.NASDAQ_CELLS
SPLIT = "confirm"
# The plan, page 3: rho* per NASDAQ-listed cell; None where the plan's table says "none".
PLAN_CEILINGS = {
    "Q-T1S1": 1 / 5,
    "Q-T1S2": 1 / 20,
    "Q-T1S3": None,
    "Q-T2S1": None,
    "Q-T2S2": None,
    "Q-T2S3": 1 / 5,
    "Q-T3S1": 1 / 5,
    "Q-T3S2": 1 / 5,
    "Q-T3S3": 1 / 10,
}
ARM_TICK_BINS = ("Q-T1", "Q-T3")
# pooled.core2: R_{1/2}/R_0 is read at the range's middle horizon on the log2 grid, floor of the mean log2 T.
RRATIO_T = {"[4,32]": 8, "[8,64]": 16}

P_HEADER = "p (descriptive, no multiplicity control)"
NOT_ESTIMABLE = "not estimable: no registered ceiling"
METHODS3 = ("day_stock_bootstrap", "day_jackknife_t", "whole_stock_bootstrap")


# ---- arms ----------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Arm:
    clip: str  # clip3b rule: 'R2' (primary), 'R1', 'R3' or 'none'
    centered: bool  # False (primary): uncentred moments

    @property
    def name(self) -> str:
        if self == PRIMARY:
            return "primary (R2, uncentred)"
        parts = ([] if self.clip == PRIMARY.clip else [{"none": "no clip"}.get(self.clip, f"clip {self.clip}")]) + (
            [] if self.centered == PRIMARY.centered else ["centred moments" if self.centered else "uncentred"]
        )
        return ", ".join(parts)

    def changed(self) -> list[str]:
        """The settings in which the arm differs from the primary."""
        return [f.name for f in dataclasses.fields(self) if getattr(self, f.name) != getattr(PRIMARY, f.name)]


PRIMARY = Arm("R2", False)
ARMS = (Arm("R1", False), Arm("R3", False), Arm("none", False), Arm("R2", True))


def check_arm(arm: Arm) -> None:
    """A known clip rule, a bool centring flag, and at most one setting changed from the primary."""
    if arm.clip not in CL.RULES:
        raise ValueError(f"unknown clip rule {arm.clip!r}")
    if not isinstance(arm.centered, bool):
        raise ValueError("centered must be a bool")
    if len(arm.changed()) > 1:
        raise ValueError(f"arm {arm} changes {arm.changed()} at once; arms change one setting at a time")


def clip_sizes(rule: str, q: np.ndarray, lengths: np.ndarray, n_days: int) -> np.ndarray:
    """New wrapper (the frozen make_unit hardcodes R2): run sizes of stock-major stock-days under clip3b's rule,
    the same arguments and return as primary3d.r2_clip."""
    qs = np.split(np.asarray(q, dtype=float), np.cumsum(lengths)[:-1])
    stocks = np.repeat(np.arange(len(lengths) // n_days), n_days)
    return np.concatenate(CL.apply_rule(rule, qs, stocks))


def draws_hash(d: P.Draws) -> str:
    h = hashlib.sha256()
    for a in (*d.ds, *d.ws):
        h.update(str(a.shape).encode())
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


@contextmanager
def arm_setting(arm: Arm):
    """For the duration of one arm: primary3d.r2_clip -> clip_sizes(arm.clip); pooled.panel_moments called with
    centered=arm.centered; primary3d.draws recorded (hash of its bootstrap weights, no draw added). Yields the
    log of the settings actually applied. The frozen functions are restored afterwards."""
    check_arm(arm)
    saved = {"r2_clip": P.r2_clip, "panel_moments": PO.panel_moments, "draws": P.draws}
    log: dict = {"clip": [], "centered": [], "draws": []}

    def r2_clip(q, lengths, n_days):
        rule = arm.clip
        log["clip"].append(rule)
        return clip_sizes(rule, q, lengths, n_days)

    def panel_moments(*a, **kw):
        if "centered" in kw or len(a) > 8:
            raise RuntimeError("a frozen caller passed centered itself")
        c = arm.centered
        log["centered"].append(c)
        return saved["panel_moments"](*a, centered=c, **kw)

    def draws(*a, **kw):
        d = saved["draws"](*a, **kw)
        log["draws"].append(draws_hash(d))
        return d

    P.r2_clip, PO.panel_moments, P.draws = r2_clip, panel_moments, draws
    try:
        yield log
    finally:
        P.r2_clip, PO.panel_moments, P.draws = saved["r2_clip"], saved["panel_moments"], saved["draws"]


def run_arm(flows: dict, arm: Arm, seed: int = SEED, **kw) -> dict:
    """The frozen run_primary under one arm's setting, from default_rng(seed)."""
    check_arm(arm)
    rng = np.random.default_rng(seed)
    with arm_setting(arm) as log:
        res = P.run_primary(flows, FROZEN_PM, rng, **kw)
    units = list(res["by_tick_bin"])
    if len(log["draws"]) != len(units):
        raise RuntimeError(f"{len(log['draws'])} draw sets for {len(units)} units")
    eff = {"clip": sorted(set(log["clip"])), "centered": sorted(set(log["centered"]))}
    if eff != {"clip": [arm.clip], "centered": [arm.centered]}:
        raise RuntimeError(f"arm {arm} applied {eff}")
    return {
        "arm": arm,
        "seed": seed,
        "result": res,
        "draws": dict(zip(units, log["draws"], strict=True)),
        "effective": eff,
    }


def check_seeds(prim: dict, arms: dict) -> None:
    """Every arm from SEED, and the first unit's bootstrap draws (taken before any surrogate draw) identical to
    the primary's."""
    first = next(iter(prim["draws"]))
    for name, a in {"primary": prim, **arms}.items():
        if a["seed"] != SEED:
            raise RuntimeError(f"arm {name} ran with seed {a['seed']}, not {SEED}")
        if a["draws"].get(first) != prim["draws"][first]:
            raise RuntimeError(f"arm {name}: the first unit's draws differ from the primary's (another seed)")


def run_arms(flows: dict, prim: dict | None = None, **kw) -> tuple[dict, dict]:
    """(primary run, {arm name: run}) for PRIMARY and every arm of ARMS."""
    for a in ARMS:
        if len(a.changed()) != 1:
            raise ValueError(f"arm {a} changes {a.changed()}; each arm changes exactly one setting")
    prim = run_arm(flows, PRIMARY, **kw) if prim is None else prim
    arms = {a.name: run_arm(flows, a, **kw) for a in ARMS}
    check_seeds(prim, arms)
    return prim, arms


def anchor(res: dict, stored_result: dict) -> dict:
    """Every leaf of run_primary's return against results.json["result"] (floats by float.hex)."""
    cmp = Compare()
    cmp(res, stored_result, "result")
    return {
        "identical": not cmp.bad,
        "fields_compared": cmp.n,
        "mismatches": len(cmp.bad),
        "mismatch_list": cmp.bad[:50],
    }


# ---- statistics of a test result -------------------------------------------------------------------------------

# (label, statistic key, index into it or None, p getter, observe() methods key, reference values)
QUANTITIES = {
    "P-V": (
        ("slope in a of zeta_V real minus S1", "slope_S1", None, lambda s: s["p_S1"], "slope_S1", (0.0,)),
        ("slope in a of zeta_V real minus S2", "slope_S2", None, lambda s: s["p_S2"], "slope_S2", (0.0,)),
    ),
    "P-R": (
        ("a*(T1)", "a_star_T1_T2", 0, lambda s: s["p_end"][0], "a_star_vs_0.35", (0.35, 0.65)),
        ("a*(T2)", "a_star_T1_T2", 1, lambda s: s["p_end"][1], "a_star_vs_0.35", (0.35, 0.65)),
    ),
    "P-M": (
        ("slope of beta-hat on log q", "slope", None, lambda s: s["p_slope"], "slope", (0.0,)),
        ("top size bin beta-hat", "top", None, lambda s: s["p_top"], "top", (0.0,)),
    ),
    "P-F": (("slope of zeta_K(a) on a over [0, 1]", "pf", None, lambda s: s["p"], "slope", (0.0,)),),
}


def at(v, i):
    """Element i of a summary field (a list), or the field itself (a scalar, or a list of one)."""
    if i is None:
        return v[0] if isinstance(v, list) and len(v) == 1 else v
    return v[i]


def num(x):
    return None if x is None or isinstance(x, str) else float(x)


def excludes(v: float, lo, hi):
    lo, hi = num(lo), num(hi)
    if lo is None or hi is None or not (np.isfinite(lo) and np.isfinite(hi)):
        return None
    return bool(v < lo or v > hi)


def tr_dict(tr) -> dict:
    """A TestResult as a JSON-ready dict, without the family decisions (rejected, passed)."""
    d = tr if isinstance(tr, dict) else dataclasses.asdict(tr)
    return {k: v for k, v in d.items() if k not in ("rejected", "passed")}


def quantity_rows(test: str, tr: dict) -> list[dict]:
    """One row per statistic of a test at its range used: estimate, widest 0.90 interval and method, censoring."""
    s = tr["statistic"]
    out = []
    for label, k, i, pget, mkey, refs in QUANTITIES[test]:
        if not s:
            out.append({"test": test, "quantity": label, "defined": False, "refs": list(refs), "mkey": mkey, "i": i})
            continue
        sm = s[k]
        row = {
            "test": test,
            "quantity": label,
            "defined": True,
            "estimate": num(at(sm["estimate"], i)),
            "lo": num(at(sm["lo"], i)),
            "hi": num(at(sm["hi"], i)),
            "method": at(sm["method"], i),
            "p": num(pget(s)),
            "refs": list(refs),
            "mkey": mkey,
            "i": i,
            "censored": s["censored"][i] if test == "P-R" else None,
        }
        row["excludes"] = {f"{r:g}": excludes(r, row["lo"], row["hi"]) for r in refs}
        out.append(row)
    return out


# ---- A: per-cell results ---------------------------------------------------------------------------------------


def registered_ceiling(cell: str) -> float | None:
    """The cell's frozen ceiling (prereg3d.CONFIRMATORY), checked against the plan's table; None without one."""
    if cell not in PLAN_CEILINGS:
        raise ValueError(f"{cell} is not a NASDAQ-listed cell of the plan")
    rho = C.CONFIRMATORY.get(cell)
    if rho != PLAN_CEILINGS[cell]:
        raise RuntimeError(f"{cell}: ceiling {rho} differs from the plan's {PLAN_CEILINGS[cell]}")
    return rho


@contextmanager
def capture_reps():
    """Record the replicates primary3d.run_stats returns (the call and its result unchanged; no draw taken)."""
    saved, rec = P.run_stats, []

    def run_stats(*a, **kw):
        out = saved(*a, **kw)
        rec.append(out)
        return out

    P.run_stats = run_stats
    try:
        yield rec
    finally:
        P.run_stats = saved


def pk_methods(reps: P.Reps, stored: dict) -> dict:
    """Each method's interval of vk_diff (P-K) at each range, checked against run_unit's own P-K summary."""
    out = {}
    for j, rn in enumerate(P.RANGE_NAMES):
        r = reps[j]
        if F4.jsonable(P.summary(r)) != F4.jsonable(stored[rn]):
            raise RuntimeError(f"P-K {rn}: recorded replicates do not give run_unit's summary")
        m = F4.methods(r, 0.0, "two-sided")
        out[rn] = {"estimate": m["estimate"], **{k: {"lo": m[k]["lo"], "hi": m[k]["hi"]} for k in METHODS3}}
    return out


def run_cells(flows: dict, **kw) -> dict:
    """{cell: record} for every NASDAQ-listed cell, in NASDAQ_CELLS order."""
    out = {}
    for cell in NASDAQ_CELLS:
        if cell not in flows:
            out[cell] = {"estimable": False, "reason": "no stocks of the cell in the input"}
            continue
        rho = registered_ceiling(cell)
        if rho is None:
            out[cell] = {"estimable": False, "reason": NOT_ESTIMABLE, "n_stocks": flows[cell].n_stocks}
            continue
        u = P.cell_unit(flows, cell)
        ev = P.evaluability(u.Ns, u.rho, FROZEN_PM)
        r_used, reach = P.select_range(u.Ns, u.rho)
        tests = tuple(t for t in C.PRIMARY_TESTS if ev[t]["evaluable"])
        rec = {
            "estimable": True,
            "ceiling": rho,
            "n_stocks": u.flow.n_stocks,
            "evaluability": ev,
            "range_used": None if r_used is None else F4.key(r_used),
            "reach": {F4.key(r): n for r, n in reach.items()},
            "seed": SEED,
            "tests": {},
            "methods": [],
            "P-K": None,
        }
        if tests:
            rng = np.random.default_rng(SEED)
            with F4.observe() as obs, capture_reps() as cap:
                o = P.run_unit(u, FROZEN_PM, rng, tests=tests, **kw)
            F4.attribute({"by_tick_bin": {cell: o}}, obs)
            for t, tr in o["tests"].items():
                if tr.rejected is not None or tr.passed is not None:
                    raise RuntimeError(f"{cell} {t} carries a family decision")
            rec["tests"] = {t: tr_dict(tr) for t, tr in o["tests"].items()}
            rec["methods"] = obs
            rng_reps = [c for c in cap if "vk_diff" in c]
            if r_used is not None:
                if len(rng_reps) != 1:
                    raise RuntimeError(f"{cell}: {len(rng_reps)} range-statistic replicate sets")
                rec["P-K"] = {
                    "summary": o["descriptive"]["P-K"],
                    "methods": pk_methods(rng_reps[0]["vk_diff"], o["descriptive"]["P-K"]),
                }
        out[cell] = rec
    return out


# ---- B: P-R secondaries from Phase 4 ------------------------------------------------------------------------------


def section_b(stored: dict) -> dict:
    """P-R's secondaries per tick bin as stored by Phase 4, with their JSON paths; nothing recomputed."""
    out = {}
    for tb in C.TICK_BINS:
        base = f"result.by_tick_bin.{tb}.tests.P-R"
        t = stored["result"]["by_tick_bin"].get(tb, {}).get("tests", {}).get("P-R")
        if t is None:
            out[tb] = {"stored": False, "path": base}
            continue
        sec = t["secondary"]
        if sec["psi_hat_note"] != P.PSI_NOTE:
            raise RuntimeError(f"{tb}: the stored psi-hat note is not primary3d.PSI_NOTE")
        out[tb] = {
            "stored": True,
            "path": base,
            "range_used": F4.key(tuple(t["range_used"])),
            "a_star": t["statistic"]["a_star_T1_T2"],
            "censored": t["statistic"]["censored"],
            "psi_hat": sec["psi_hat"],
            "psi_n_stocks": sec["psi_n_stocks"],
            "psi_hat_note": sec["psi_hat_note"],
            "Rratio": sec["Rratio"],
            "Rratio_minus_exp_sig2_over_8": sec["Rratio_minus_exp_sig2_over_8"],
            "sigma_l2": sec["sigma_l2"],
            "dastar_dlogT": sec["dastar_dlogT"],
        }
    return out


# ---- C: arm records ----------------------------------------------------------------------------------------------


def arm_record(run: dict) -> dict:
    """The arm's statistics in Q-T1 and Q-T3 (no family decision kept)."""
    res = run["result"]
    tb_out = {}
    for tb in ARM_TICK_BINS:
        u = res["by_tick_bin"].get(tb)
        if u is None:
            continue
        tests = {t: tr_dict(tr) for t, tr in u["tests"].items()}
        pv = tests.get("P-V", {}).get("statistic") or {}
        tb_out[tb] = {
            "rows": [r for t in C.PRIMARY_TESTS if t in tests for r in quantity_rows(t, tests[t])],
            "range_used": {t: tr["range_used"] for t, tr in tests.items()},
            "pv_loo_errors": pv.get("loo_errors"),
            "tests": tests,
        }
    return {
        "arm": dataclasses.asdict(run["arm"]),
        "name": run["arm"].name,
        "seed": run["seed"],
        "draws": run["draws"],
        "effective": run["effective"],
        "by_tick_bin": tb_out,
    }


# ---- data, integrity ------------------------------------------------------------------------------------------------


def check_confirm_only(syms: dict, universe: Path = E.UNIVERSE) -> None:
    """Every loaded symbol is a confirmation stock of its cell; raises otherwise."""
    _, u = E.read_universe(universe)
    bad = sorted((c, s, u.get(s)) for c, ss in syms.items() for s in ss if u.get(s) != (c, SPLIT))
    if bad:
        raise RuntimeError(f"stocks outside the confirmation split reached the analysis: {bad[:5]}")


def load_confirm(data: Path | None = None, universe: Path = E.UNIVERSE) -> tuple[dict, dict]:
    flows, syms = F4.build_flows("confirm", data=data, universe=universe)
    check_confirm_only(syms, universe)
    return flows, syms


def sha256(path: Path) -> str:
    return F4.sha256(path)


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True, text=True).stdout.strip()


def protected_snapshot() -> dict:
    """SHA-256 of every file under data/phase4/ and data/phase4_repro/, and the tags."""
    files = {str(p.relative_to(REPO)): sha256(p) for d in PROTECTED for p in sorted(d.rglob("*")) if p.is_file()}
    return {"files": files, "tags": git("show-ref", "--tags")}


def preflight() -> dict:
    """HEAD descends from the frozen tag; no file of the tag under tools/ changed; files added under tools/ are
    only ALLOWED_ADDED; manifest, universe and the plan PDF as registered."""
    if git("rev-parse", f"{F4.FROZEN_TAG}^{{commit}}") != F4.FROZEN_COMMIT:
        raise RuntimeError(f"{F4.FROZEN_TAG} is not {F4.FROZEN_COMMIT}")
    if subprocess.run(["git", "-C", str(REPO), "merge-base", "--is-ancestor", F4.FROZEN_TAG, "HEAD"]).returncode:
        raise RuntimeError("HEAD does not descend from the frozen tag")
    changed = git("diff", "--name-only", "--diff-filter=MDRTCU", F4.FROZEN_TAG, "--", "tools/")
    if changed:
        raise RuntimeError(f"frozen files changed: {changed.split()}")
    added = set(git("diff", "--name-only", "--diff-filter=A", F4.FROZEN_TAG, "--", "tools/").split())
    added |= set(git("ls-files", "--others", "--exclude-standard", "--", "tools/").split())
    if not added <= set(ALLOWED_ADDED):
        raise RuntimeError(f"files under tools/ beyond the allowed: {sorted(added - set(ALLOWED_ADDED))}")
    shas = {
        "data/manifest.tsv": sha256(REPO / "data" / "manifest.tsv"),
        "data/universe.tsv": sha256(REPO / "data" / "universe.tsv"),
    }
    if shas != {"data/manifest.tsv": F4.MANIFEST_SHA256, "data/universe.tsv": F4.UNIVERSE_SHA256}:
        raise RuntimeError(f"data hashes differ from the registration: {shas}")
    if not PLAN_PDF.is_file() or sha256(PLAN_PDF) != F4.PLAN_PDF_SHA256:
        raise RuntimeError(f"{PLAN_PDF}: missing or not the registered plan")
    return {
        "git_head": git("rev-parse", "HEAD"),
        "git_status": git("status", "--porcelain"),
        "sha256": shas,
        "plan_pdf": str(PLAN_PDF),
        "plan_pdf_sha256": F4.PLAN_PDF_SHA256,
        "added_files": sorted(added),
        "phase5_py_sha256": sha256(Path(__file__)),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    }


def safe_path(p) -> Path:
    """p, if Phase 5 may write it: never under data/phase4/ or data/phase4_repro/, and inside the repo only
    data/phase5/ and docs/phase5_secondary.md. Raises PermissionError otherwise."""
    p = Path(p).resolve()
    for q in PROTECTED:
        q = q.resolve()
        if p == q or q in p.parents:
            raise PermissionError(f"Phase 5 never writes under {q}: {p}")
    repo, out, doc = REPO.resolve(), OUT.resolve(), DOC.resolve()
    if (p == repo or repo in p.parents) and not (p == out or out in p.parents or p == doc):
        raise PermissionError(f"Phase 5 writes only data/phase5/ and docs/phase5_secondary.md in the repo: {p}")
    return p


def write_outputs(rec: dict, out: Path = OUT, doc: Path = DOC) -> str:
    """results.json first (so the document can be re-rendered from it alone), then the linted document."""
    out = safe_path(out)
    out.mkdir(parents=True, exist_ok=True)
    js = json.dumps(F4.jsonable(rec), indent=1, allow_nan=False)
    safe_path(out / "results.json").write_text(js)
    text = render_doc(json.loads(js))
    safe_path(doc).write_text(text)
    return text


# ---- the document ----------------------------------------------------------------------------------------------------

FORBIDDEN = re.compile(
    r"\b(confirms?|confirmed|refutes?|refuted|supports?|supported|supporting|contradicts?|"
    r"contradicted|contradicting|rejects?|rejected|rejection|holm|pass|passes|passed)\b",
    re.I,
)
LABELS = ("- **Measurement", "- **Choice", "- Planned, not run in Phase 5")


def lint_doc(text: str) -> None:
    """No em or en dash; no sentence of confirmation or refutation, no family decision; every bullet labelled;
    every line that mentions psi-hat carries primary3d.PSI_NOTE verbatim; every p column headed P_HEADER."""
    if "\u2014" in text or "\u2013" in text:
        raise ValueError("an em or en dash in the document")
    for n, line in enumerate(text.splitlines(), 1):
        m = FORBIDDEN.search(re.sub(r"`[^`]*`", "", line))  # code spans are names, not sentences
        if m:
            raise ValueError(f"line {n}: forbidden word {m.group(0)!r}: {line[:100]}")
        if line.startswith("- ") and not line.startswith(LABELS):
            raise ValueError(f"line {n}: a bullet without a Measurement or Choice label: {line[:100]}")
        if re.search(r"psi[-_ ]?hat|\u03c8", line, re.I) and P.PSI_NOTE not in line:
            raise ValueError(f"line {n}: psi-hat without its note: {line[:100]}")
        if line.startswith("|") and not set(line) <= set("|-: "):
            for c in (c.strip() for c in line.strip().strip("|").split("|")):
                if (re.fullmatch(r"p(\s.*)?", c) or re.search(r"p-value|p value|\bp_", c, re.I)) and c != P_HEADER:
                    raise ValueError(f"line {n}: p column {c!r} is not headed {P_HEADER!r}")


def f4(x) -> str:
    return F4.f4(x)


def ivl(est, lo, hi, method=None) -> str:
    s = f"{f4(est)} [{f4(lo)}, {f4(hi)}]"
    return s + (f" ({method})" if method else "")


def yn(x) -> str:
    return "n/a" if x is None else ("yes" if x else "no")


CHOICES = (
    (
        "Cells without a registered ceiling",
        "Q-T1S3, Q-T2S1 and Q-T2S2 have rho* 'none' in the plan's table (page "
        '3), and the plan does not say how to estimate them; they are reported as "' + NOT_ESTIMABLE + '" and '
        "nothing is computed for them. No ceiling is invented.",
    ),
    (
        "Per-cell generator",
        "each per-cell unit runs from its own numpy default_rng(20261102), so a cell's numbers "
        "do not depend on which other cells ran or in which order.",
    ),
    (
        "Per-cell tests computed",
        "only the tests the frozen reach rule (primary3d.evaluability, day lengths and "
        "ceilings alone) finds evaluable in the cell are computed, as run_primary does for a tick bin.",
    ),
    (
        "Per-cell P-K",
        "vk_diff is reported at the range the reach rule selects for the cell (both ranges are in "
        "data/phase5/results.json), evaluable when the reach rule gives a range, and with no p-value: the plan fixes "
        "no reference value or threshold for P-K.",
    ),
    (
        "Per-cell p-values",
        "per statistic, the frozen p-value (the largest of the three methods') against the "
        "reference value the primary test uses: P-V slopes one-sided below 0, P-R the TOST p at that end against "
        "0.35 and 0.65 (1 at a censored end), P-M slope and top bin one-sided below 0, P-F two-sided against 0; and "
        "per test, the frozen decision function's p-value. No multiplicity control is applied.",
    ),
    (
        "Q-T2S3 per cell",
        "Q-T2S3 is Q-T2's only cell, so its per-cell unit holds the same stocks as Phase 4's Q-T2 "
        "unit; its numbers here differ from Phase 4's descriptive Q-T2 numbers only through the draws (a fresh "
        "generator here, the generator after Q-T1 in Phase 4).",
    ),
    (
        "How an arm changes a setting",
        "the frozen make_unit hardcodes R2 and the frozen moments_raw does not "
        "expose centring, so two small wrappers in phase5.py are used for the duration of an arm: primary3d.r2_clip "
        "is replaced by clip_sizes, which calls the frozen clip3b.apply_rule(rule) on the same stock-days, and "
        "pooled.panel_moments is called with its frozen centered flag. The primary arm runs through the same "
        "wrappers (R2 by clip3b.apply_rule) and the anchor shows it bit-identical to Phase 4.",
    ),
    (
        "Clip arm scope",
        "the arm's rule replaces R2 wherever the primary applies R2 (every use of the unit's run "
        "sizes: window moments, surrogates, M3, sigma_l^2, the size bins of P-M and of the impact slope), because "
        "the plan caps the sizes entering any moment; nothing else changes.",
    ),
    (
        "R1 definition",
        "R1 is clip3b's frozen R1: each run capped at 1 % of the stock-day's summed run volume in "
        "this sample (NASDAQ displayed runs), the only R1 in the frozen code; consolidated volume is not in the "
        "data.",
    ),
    (
        "Centring scope",
        "centred moments subtract each stock-day's window means in V, K and S (pooled.panel_moments "
        "centered=True) for the real flow and for the surrogates alike; the M3 prediction, sigma_l^2, the impact "
        "slope's inputs and the P-M propagator moments are unchanged, since the frozen code exposes centring only "
        "for the window moments. P-M therefore differs from the primary only if the draws differ.",
    ),
    (
        "Generator in the arms",
        "each arm reruns the frozen run_primary over all nine cells from a fresh "
        "default_rng(20261102), in the frozen order Q-T1, Q-T2, Q-T3. Q-T1's bootstrap draws come before any "
        "surrogate draw and are identical in every arm (checked). A clip rule can change the S2 block lengths "
        "(surrogates.choose_L reads the sizes) and so the draws of later units; whether each unit's draws equal "
        "the primary's is measured and reported. No unit is re-seeded, which would change the primary path.",
    ),
    (
        "Q-T2 in the arms",
        "Q-T2 is computed inside each arm's run_primary (to keep the frozen generator order) and "
        "not reported; the arms' scope is Q-T1 and Q-T3.",
    ),
    (
        "Excludes the reference value",
        "an interval excludes v when v < lo or v > hi; for P-R both 0.35 and 0.65 "
        "are shown; P-M's top-bin statistic uses 0 like its slope. The column is descriptive; nothing is decided "
        "from it.",
    ),
    (
        "Primary column of section C",
        "the primary values are those of the anchor run, bit-identical to data/phase4/results.json.",
    ),
    (
        "Section B ranges",
        "both ranges are shown for the range-dependent secondaries, the range used marked; "
        "values are read from data/phase4/results.json and shown to 4 decimals.",
    ),
)

PLANNED_NOT_RUN = (
    "P-K: the zeta_K(a) curve, unclipped and under each clip rule, with real minus S1 and S2 (plan page 7). "
    "Section A gives only the per-cell vk_diff summary.",
    "P-C: the collapse exponent chi-hat beside 1/(1 + gamma-hat) and 1 - gamma-hat/2, gamma-hat by local "
    "Whittle (bandwidth n^0.65), a power-law fit to the sign autocorrelation and the slope of v_0(T), each with "
    "its measured bias (plan page 7).",
    "Model comparison M0 (linear in T), M1, M2 and M3 for K_a(T) (plan page 7).",
    "The P-R secondaries (R_{1/2}/R_0 against e^{sigma_l^2/8}, the impact slope beside a*(T), da*/dlog T) per "
    "cell; section B gives them per tick bin, from Phase 4.",
    "The clip and centring arms for Q-T2, for the per-cell units and for P-K and P-C; section C covers the four "
    "primary statistics in Q-T1 and Q-T3.",
    "Hidden-order prints signed by the recorded rule (side of the mid; same-timestamp prints at the mid take the "
    "run's sign; others unsigned and dropped, their share reported) (plan page 7).",
    "Fill-level analysis, with sigma_l re-measured on fills (plan page 7).",
    "Hierarchical pooling (DerSimonian-Laird with Hartung-Knapp-Sidik-Jonkman variance) (plan page 7).",
    "The rescaled axis u = a sigma-hat_l^2 for P-V and P-K (plan page 7).",
    "Deseasonalised mid returns (plan page 7).",
    "Sub-windows 10:00 to 12:00, 12:00 to 14:00 and 14:00 to 15:30 (plan page 7).",
    "The NYSE-listed secondary universe, with the same frozen code (plan pages 3 and 7).",
    "The venue-thinning experiment: observed shares {1.0, 0.46, 0.26, 0.17}, hidden share 0.17, random and "
    "router-like thinning, kappa in {1, 4, 16} (plan page 8).",
    "The NASDAQ-share dose-response check, after the consolidated-volume vendor is reconciled on about ten names "
    "(plan page 8).",
    "The nine replication days, with the same frozen code (plan pages 2 and 8).",
)


def method_cells(rec: dict, cell: str, test: str, rng_name, mkey: str, i) -> list[str]:
    """Each method's [lo, hi] of one statistic from the observe() record."""
    for r in rec["methods"]:
        if r["unit"] == cell and r["test"] == test and r["range"] == rng_name:
            m = r["methods"][mkey]
            return [f"[{f4(at(m[k]['lo'], i))}, {f4(at(m[k]['hi'], i))}]" for k in METHODS3]
    return ["n/a"] * 3


def doc_a(rec: dict) -> list[str]:
    cells = rec["cells"]
    L = [
        "## A. Per-cell results, all nine NASDAQ-listed cells",
        "",
        "**Measurement** (computed here; data/phase5/results.json, key `cells`). Each cell with a registered "
        "ceiling is one pooling unit; the reach rule decides which tests are computed. Intervals are 0.90.",
        "",
        "### A.1 Cells and tests",
        "",
        f"| cell | rho* | stocks | test | evaluable (frozen reach rule) | reach rule's reason | range used | "
        f"stocks at T2: [4,32] | stocks at T2: [8,64] | {P_HEADER} | sub-conditions |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for c, r in cells.items():
        if not r["estimable"]:
            L.append(
                f"| {c} | none | {r.get('n_stocks', 'n/a')} | all | n/a | {r['reason']} | n/a | n/a | n/a | n/a | n/a |"
            )
            continue
        rho = f"1/{round(1 / r['ceiling'])}"
        for t in (*C.PRIMARY_TESTS, "P-K"):
            if t == "P-K":
                ev, why = r["range_used"] is not None, r["evaluability"]["P-F"]["reason"]
                used, p, conds = r["range_used"] or "n/a", "n/a (no reference value)", "n/a"
            else:
                ev, why = r["evaluability"][t]["evaluable"], r["evaluability"][t]["reason"]
                tr = r["tests"].get(t)
                used = (
                    "n/a" if tr is None else (F4.key(tuple(tr["range_used"])) if tr["range_used"] else tr["range_rule"])
                )
                p = "n/a (not computed)" if tr is None else f4(tr["p"])
                conds = (
                    "n/a"
                    if tr is None
                    else ", ".join(f"{k} {'holds' if v else 'fails'}" for k, v in tr["conditions"].items())
                )
            L.append(
                f"| {c} | {rho} | {r['n_stocks']} | {t} | {yn(ev)} | {why} | {used} | {r['reach']['[4,32]']} "
                f"| {r['reach']['[8,64]']} | {p} | {conds} |"
            )
    L += [
        "",
        "### A.2 Statistics, each method's interval",
        "",
        "**Measurement** (computed here). Widest-of-three 0.90 interval and the method that gave it; then each "
        "method's 0.90 interval. P-R rows at a censored end show the end's a.",
        "",
        f"| cell | test | statistic | range | estimate | widest 0.90 interval (method) | day-and-stock bootstrap | "
        f"day jackknife t | whole-stock bootstrap | {P_HEADER} |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for c, r in cells.items():
        if not r["estimable"]:
            L.append(f"| {c} | all | {r['reason']} | n/a | n/a | n/a | n/a | n/a | n/a | n/a |")
            continue
        for t in C.PRIMARY_TESTS:
            tr = r["tests"].get(t)
            if tr is None:
                L.append(
                    f"| {c} | {t} | not computed: {r['evaluability'][t]['reason']} | n/a | n/a | n/a | n/a | "
                    "n/a | n/a | n/a |"
                )
                continue
            rn = None if t == "P-M" else (F4.key(tuple(tr["range_used"])) if tr["range_used"] else None)
            rlab = tr["range_rule"] if t == "P-M" else (rn or "n/a")
            for q in quantity_rows(t, tr):
                if not q["defined"]:
                    L.append(
                        f"| {c} | {t} | {q['quantity']} | {rlab} | undefined | n/a | n/a | n/a | n/a | {f4(tr['p'])} |"
                    )
                    continue
                cen = f", censored {q['censored']}" if q["censored"] else ""
                mc = method_cells(r, c, t, rn, q["mkey"], q["i"])
                L.append(
                    f"| {c} | {t} | {q['quantity']}{cen} | {rlab} | {f4(q['estimate'])} | "
                    f"[{f4(q['lo'])}, {f4(q['hi'])}] ({q['method'] or 'undefined'}) | "
                    + " | ".join(mc)
                    + f" | {f4(q['p'])} |"
                )
        pk = r["P-K"]
        if pk is None:
            L.append(
                f"| {c} | P-K | vk_diff | n/a | not computed (no range under the reach rule) | n/a | n/a | n/a | "
                "n/a | n/a |"
            )
        else:
            rn = r["range_used"]
            s, m = pk["summary"][rn], pk["methods"][rn]
            mc = [f"[{f4(at(m[k]['lo'], None))}, {f4(at(m[k]['hi'], None))}]" for k in METHODS3]
            L.append(
                f"| {c} | P-K | vk_diff = slope zeta_V - 2 slope zeta_K | {rn} | {f4(at(s['estimate'], None))} "
                f"| [{f4(at(s['lo'], None))}, {f4(at(s['hi'], None))}] ({at(s['method'], None) or 'undefined'}) | "
                + " | ".join(mc)
                + " | n/a (no reference value) |"
            )
    return L


def doc_b(rec: dict) -> list[str]:
    b = rec["section_b"]
    L = [
        "## B. P-R secondaries per tick bin (from Phase 4)",
        "",
        "**Measurement** (computed in Phase 4 by primary3d.run_unit and stored in data/phase4/results.json; read "
        "here, not recomputed; each value's JSON path is given). Intervals are 0.90, widest of three.",
        "",
    ]
    for tb, s in b.items():
        if not s["stored"]:
            L.append(f"- **Measurement** ({tb}): no P-R record at `{s['path']}`.")
            continue
        ps = s["psi_hat"]
        L.append(
            f"- **Measurement** ({tb}, `{s['path']}.secondary.psi_hat`): psi-hat "
            f"{ivl(ps['estimate'], ps['lo'], ps['hi'], at(ps['method'], None))} over {s['psi_n_stocks']:g} "
            f"stocks with a slope. Note: {s['psi_hat_note']}"
        )
    L += [
        "",
        f"| tick bin | quantity | range | estimate | widest 0.90 interval (method) | {P_HEADER} | JSON path |",
        "|---|---|---|---|---|---|---|",
    ]
    for tb, s in b.items():
        if not s["stored"]:
            continue
        p = s["path"]
        for rn in P.RANGE_NAMES:
            mark = f"{rn}{' (range used)' if rn == s['range_used'] else ''}"
            for lab, k, has_p in (
                (f"R_{{1/2}}/R_0 at T = {RRATIO_T[rn]}", "Rratio", False),
                (f"R_{{1/2}}/R_0 minus e^{{sigma_l^2/8}} at T = {RRATIO_T[rn]}", "Rratio_minus_exp_sig2_over_8", True),
                ("da*/dlog T", "dastar_dlogT", True),
            ):
                v = s[k][rn]
                m = at(v["method"], None) or "undefined"
                pv = f4(v["p_two_sided"]) + " (two-sided)" if has_p else "n/a"
                L.append(
                    f"| {tb} | {lab} | {mark} | {f4(v['estimate'])} | [{f4(v['lo'])}, {f4(v['hi'])}] ({m}) | {pv} | "
                    f"`{p}.secondary.{k}.{rn}` |"
                )
        v = s["sigma_l2"]
        m = at(v["method"], None) or "undefined"
        L.append(
            f"| {tb} | sigma_l^2 (pooled variance of log run size) | n/a | {f4(v['estimate'])} | "
            f"[{f4(v['lo'])}, {f4(v['hi'])}] ({m}) | n/a | `{p}.secondary.sigma_l2` |"
        )
        a = s["a_star"]
        for i, lab in enumerate(("a*(T1)", "a*(T2)")):
            cen = f", censored {s['censored'][i]}" if s["censored"][i] else ""
            L.append(
                f"| {tb} | {lab}{cen} | {s['range_used']} (range used) | {f4(a['estimate'][i])} | "
                f"[{f4(a['lo'][i])}, {f4(a['hi'][i])}] ({a['method'][i] or 'undefined'}) | n/a | "
                f"`{p}.statistic.a_star_T1_T2` |"
            )
    q = b.get("Q-T1", {})
    if q.get("stored"):
        rn = q["range_used"]
        lo_t, hi_t = (int(x) for x in rn.strip("[]").split(","))
        a, cen = q["a_star"], q["censored"]
        rr, rr2 = q["Rratio"][rn], q["Rratio"]["[4,32]" if rn == "[8,64]" else "[8,64]"]
        da, ps = q["dastar_dlogT"][rn], q["psi_hat"]

        def side(v):
            e = num(v["estimate"])
            pos = "below 1" if e < 1 else ("above 1" if e > 1 else "equal to 1")
            return f"{pos}, its interval {'excludes' if excludes(1.0, v['lo'], v['hi']) else 'includes'} 1"

        def ast(i):
            t = (lo_t, hi_t)[i]
            s = f"a*(T = {t}) is {ivl(a['estimate'][i], a['lo'][i], a['hi'][i])}"
            return s + (
                f", censored {cen[i]} (the grid maximum of pooled R_a({t}) is at the "
                f"{'lowest' if cen[i] == 'low' else 'highest'} a of the grid)"
                if cen[i]
                else ""
            )

        L += [
            "",
            "### a* in Q-T1, with the secondaries beside it",
            "",
            f"**Measurement** (Phase 4 values, read from data/phase4/results.json). In Q-T1 over the range used "
            f"{rn}, {ast(0)}, and {ast(1)}. Beside it: R_{{1/2}}/R_0 is {ivl(rr['estimate'], rr['lo'], rr['hi'])} "
            f"at T = {RRATIO_T[rn]} ({side(rr)}) and {ivl(rr2['estimate'], rr2['lo'], rr2['hi'])} at T = "
            f"{RRATIO_T['[4,32]' if rn == '[8,64]' else '[8,64]']} ({side(rr2)}); da*/dlog T over {rn} is "
            f"{ivl(da['estimate'], da['lo'], da['hi'])}, its interval "
            f"{'excludes' if excludes(0.0, da['lo'], da['hi']) else 'includes'} 0; psi-hat is "
            f"{ivl(ps['estimate'], ps['lo'], ps['hi'])}. Note: {q['psi_hat_note']}. These values are placed side "
            "by side as stored; no causal reading is drawn from them.",
        ]
    return L


def doc_c(rec: dict) -> list[str]:
    prim, arms = rec["arms"]["primary"], rec["arms"]["arms"]
    L = [
        "## C. Robustness arms in Q-T1 and Q-T3",
        "",
        "**Measurement** (computed here; data/phase5/results.json, key `arms`). Each arm changes one setting "
        "from the primary (clip R2, uncentred moments), everything else at the primary setting, seed 20261102. "
        "Estimates at the range each test used (the reach rule does not depend on the arm). The excludes "
        "columns are descriptive; nothing is decided from them.",
        "",
        "| tick bin | test | statistic | reference | primary: estimate [0.90 interval] (method) | primary interval "
        "excludes reference | arm | arm: estimate [0.90 interval] (method) | arm interval excludes reference |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    def cell(q):
        if not q["defined"]:
            return "undefined", "n/a"
        cen = f", censored {q['censored']}" if q["censored"] else ""
        ex = (
            "; ".join(f"{k}: {yn(v)}" for k, v in q["excludes"].items())
            if len(q["excludes"]) > 1
            else yn(next(iter(q["excludes"].values())))
        )
        return f"{ivl(q['estimate'], q['lo'], q['hi'], q['method'])}{cen}", ex

    for tb in ARM_TICK_BINS:
        prow = prim["by_tick_bin"].get(tb, {}).get("rows", [])
        for i, pq in enumerate(prow):
            ptxt, pex = cell(pq)
            ref = " and ".join(f"{r:g}" for r in pq["refs"])
            for a in arms.values():
                aq = a["by_tick_bin"][tb]["rows"][i]
                if (aq["test"], aq["quantity"]) != (pq["test"], pq["quantity"]):
                    raise RuntimeError("arm rows out of order")
                atxt, aex = cell(aq)
                L.append(
                    f"| {tb} | {pq['test']} | {pq['quantity']} | {ref} | {ptxt} | {pex} | {a['name']} | {atxt} "
                    f"| {aex} |"
                )
    L += [
        "",
        "### C.1 P-V (i): leave-one-day-out weighted error of log v_a(T), by arm",
        "",
        "**Measurement** (computed here). The P-V statistics above are its part (ii); (i) is a gate without an "
        "interval.",
        "",
        "| tick bin | arm | M1 | M2 | M3 | days used |",
        "|---|---|---|---|---|---|",
    ]
    for tb in ARM_TICK_BINS:
        for a in (prim, *arms.values()):
            e = a["by_tick_bin"].get(tb, {}).get("pv_loo_errors")
            if e:
                L.append(f"| {tb} | {a['name']} | {f4(e['M1'])} | {f4(e['M2'])} | {f4(e['M3'])} | {e['days_used']} |")
    L += [
        "",
        "### C.2 Settings applied and draws",
        "",
        "**Measurement** (computed here). The clip rule and centring flag the wrappers actually applied in each "
        "arm, and whether each unit's bootstrap draws (SHA-256 of the weights) equal the primary's.",
        "",
        "| arm | clip rule applied | centred flag applied | seed | draws equal to the primary's: Q-T1 | Q-T2 | Q-T3 |",
        "|---|---|---|---|---|---|---|",
    ]
    for a in (prim, *arms.values()):
        eq = [yn(a["draws"].get(u) == prim["draws"].get(u)) if u in prim["draws"] else "n/a" for u in C.TICK_BINS]
        L.append(
            f"| {a['name']} | {', '.join(a['effective']['clip'])} | "
            f"{', '.join(str(x) for x in a['effective']['centered'])} | {a['seed']} | " + " | ".join(eq) + " |"
        )
    return L


def render_doc(rec: dict) -> str:
    """The document, from the JSON record alone; linted before it is returned."""
    pre, anc = rec["preflight"], rec["anchor"]
    rh = rec["runs_npz_sha256"]
    L = [
        "# Phase 5: secondary and robustness analyses (descriptive)",
        "",
        "Secondary and descriptive. Nothing here is confirmatory, nothing here changes any Phase 4 result, and "
        "nothing here is a test of the paper. No multiplicity control is applied and no decision is made.",
        "",
        f"- **Measurement.** Run once on {rec['started']} with `tools/est/phase5.py run`, at git HEAD "
        f"{pre['git_head']} (files under tools/ added after the freeze: {', '.join(pre['added_files'])}; "
        f"tools/est/phase5.py SHA-256 {pre['phase5_py_sha256']}), seed {rec['seed']}, B = {rec['n_boot']}, "
        f"surrogate draws {rec['n_sur']}, P-M design {rec['pm']}.",
        f"- **Measurement.** Registered plan read from {pre['plan_pdf']}, SHA-256 {pre['plan_pdf_sha256']} (as "
        f"registered on OSF, {F4.OSF}). data/manifest.tsv and data/universe.tsv as registered.",
        f"- **Measurement.** Runs files: the SHA-256 of all {len(rh)} study days' runs files matched "
        f"data/phase4/results.json (`runs_npz_sha256`) before anything was computed: "
        + "; ".join(f"{d} {h[:16]}" for d, h in rh.items())
        + " (full hashes in data/phase5/results.json).",
        "- **Measurement.** Confirmation stocks only, loaded with `phase4.build_flows('confirm')`; every symbol "
        "checked against the universe's confirmation split; counts per cell as in the plan and as in Phase 4: "
        + ", ".join(f"{c} {n}" for c, n in rec["stock_counts"].items())
        + ".",
        f"- **Measurement.** data/phase4/ and data/phase4_repro/ ({rec['protected']['n_files']} files) and the "
        f"tags were hashed before and after the run: {'unchanged' if rec['protected']['unchanged'] else 'CHANGED'}.",
        f"- **Measurement.** Wall time {rec['wall_s']:.0f} s (loading {rec['load_s']:.0f} s, anchor "
        f"{rec['anchor_s']:.0f} s, arms {rec['arms_s']:.0f} s, cells {rec['cells_s']:.0f} s).",
        "",
        "## Anchor",
        "",
        f"- **Measurement.** phase5's code path with every setting at the primary (clip R2 through "
        f"clip3b.apply_rule, centered=False through pooled.panel_moments, seed 20261102) reproduces "
        f'data/phase4/results.json["result"] '
        f"{'bit for bit' if anc['identical'] else 'NOT bit for bit'}: {anc['fields_compared']} fields compared "
        f"(floats by float.hex), {anc['mismatches']} mismatches.",
        "",
        "## Choices made in Phase 5",
        "",
    ]
    L += [f"- **Choice: {k}.** {v[0].upper()}{v[1:]}" for k, v in CHOICES]
    L += [""] + doc_a(rec) + [""] + doc_b(rec) + [""] + doc_c(rec)
    L += [
        "",
        "## Planned, not run in Phase 5",
        "",
        "Every other secondary, descriptive or robustness item of the registered plan (pages 7 and 8):",
        "",
    ]
    L += [f"- Planned, not run in Phase 5: {x}" for x in PLANNED_NOT_RUN]
    text = "\n".join(L) + "\n"
    lint_doc(text)
    return text


# ---- main ------------------------------------------------------------------------------------------------------------


class AnchorFailed(RuntimeError):
    def __init__(self, anc: dict):
        super().__init__(f"anchor failed: {anc['mismatches']} of {anc['fields_compared']} fields differ")
        self.anchor = anc


def compute(flows: dict, stored: dict, **kw) -> dict:
    """The anchor (raises AnchorFailed before anything else is computed), then the arms, the cells and section B."""
    t1 = time.time()
    prim = run_arm(flows, PRIMARY, **kw)
    anc = anchor(prim["result"], stored["result"])
    if not anc["identical"]:
        raise AnchorFailed(anc)
    t2 = time.time()
    prim, arms = run_arms(flows, prim=prim, **kw)
    t3 = time.time()
    cells = run_cells(flows, **kw)
    t4 = time.time()
    return {
        "seed": SEED,
        "pm": dataclasses.asdict(FROZEN_PM),
        "n_boot": kw.get("n_boot", C.N_BOOT),
        "n_sur": kw.get("n_sur", P.N_SUR),
        "anchor": anc,
        "arms": {"primary": arm_record(prim), "arms": {k: arm_record(v) for k, v in arms.items()}},
        "cells": cells,
        "section_b": section_b(stored),
        "anchor_s": t2 - t1,
        "arms_s": t3 - t2,
        "cells_s": t4 - t3,
    }


def run() -> int:
    t0 = time.time()
    started = time.strftime("%Y-%m-%d %H:%M:%S %z")
    pre = preflight()
    prot0 = protected_snapshot()
    F4.claim_once(safe_path(OUT))
    stored = json.loads(PHASE4_RESULTS.read_text())
    data = E.data_dir()
    rh = F4.runs_hashes(data)
    if rh != stored["runs_npz_sha256"]:
        raise RuntimeError("runs files differ from those of the Phase 4 run")
    flows, syms = load_confirm(data)
    counts = {c: f.n_stocks for c, f in flows.items()}
    if counts != F4.PLAN_CONFIRM_COUNTS:
        raise RuntimeError(f"confirmation counts {counts} differ from the plan")
    if syms != {c: s["symbols"] for c, s in stored["stock_counts"].items()}:
        raise RuntimeError("confirmation symbols differ from Phase 4's")
    t1 = time.time()
    try:
        rec = compute(flows, stored)
    except AnchorFailed as e:
        safe_path(OUT / "anchor_failed.json").write_text(json.dumps(F4.jsonable(e.anchor), indent=1))
        print(json.dumps(F4.jsonable(e.anchor), indent=1))
        raise SystemExit("anchor failed: the primary setting does not reproduce Phase 4; stopping") from None
    rec.update(started=started, preflight=pre, runs_npz_sha256=rh, stock_counts=counts, load_s=t1 - t0)
    prot1 = protected_snapshot()
    if prot1 != prot0:
        raise RuntimeError("data/phase4/, data/phase4_repro/ or the tags changed during the run")
    rec["protected"] = {"n_files": len(prot0["files"]), "unchanged": True, "files": prot0["files"]}
    rec["wall_s"] = time.time() - t0
    text = write_outputs(rec)
    if protected_snapshot() != prot0:
        raise RuntimeError("data/phase4/, data/phase4_repro/ or the tags changed while writing")
    print(text)
    return 0


def render() -> int:
    rec = json.loads((OUT / "results.json").read_text())
    text = render_doc(rec)
    safe_path(DOC).write_text(text)
    print(text)
    return 0


def main(argv: list[str]) -> int:
    if len(argv) != 2 or argv[1] not in ("run", "render"):
        print(__doc__)
        return 2
    return run() if argv[1] == "run" else render()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
