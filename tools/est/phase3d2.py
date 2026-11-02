#!/usr/bin/env python3
"""Phase 3d-2 driver: operating characteristics of the four primary tests (primary3d).

pilot    PILOT_PANELS panels per world of STUDY_WORLDS at the default P-M design, timed per tick bin; prints the projected
         wall time of the study on the given workers
design   stage A: N_STUDY panels of W-null, W0, WB; per panel and tick bin P-V, P-R, P-F (one run_unit, as
         run_primary) and P-M at every design of worlds3d2.PM_DESIGNS (6 since 1 Nov 2026) on the SAME draws
power    applies the P-M design rule (analyse3d2.chosen_design, the rule in worlds3d2.choose_design), then
         stage B: N_STUDY panels of W0-psi3 at the chosen design
wm       N_STUDY panels of WM (the published PAPER2 preset) at the chosen design; not part of all (decided
         1 Nov 2026), kept for later. Its P-M fails by construction (worlds3d2.BY_CONSTRUCTION)
thin     WM thinning: THIN_SETTINGS (shares x rules x kappa) on paired market days, THIN_PANELS panels, P-M at
         the chosen design; not part of all, kept for later
analyse  rates before and after the one Holm family (primary3d.combine at the chosen design), FWER in W-null,
         ranges used, pairs that left, bin occupancy, and WM and thinning when present; writes report/*.md, .json
all      design, power, analyse in order: STUDY_WORLDS only (each step resumes where it stopped); --panels N
         caps every world's panels at N (a smoke run; PHASE3D2_OUT sets the output directory)
topup    more panels for every world of STUDY_WORLDS (--panels N), then analyse
check    integrity of every saved study chunk (loads, complete, one file per job, no panel twice); a bad or
         duplicated file, or a leftover .tmp, moves to quarantine/<time>/ (nothing is deleted)
trim     keeps only worlds3d2.PM_DESIGNS in chunks saved with the 18-design grid (decided 1 Nov 2026); each
         original moves unchanged to study_18designs/, nothing is recomputed

Every step takes an exclusive lock on <output>/phase3d2.lock (an flock, dropped by the OS when the process
ends however it ends). A second instance on the same output directory refuses to start (1 Nov 2026: two
instances of all once ran together on the same chunks and log).

Jobs are chunks of CHUNK_PANELS panels of one (world, tick bin), saved as JSON under
~/CODE/itch-data/phase3d2/{study,thin}/ and skipped when present (resumable). Seeds: panel i of world w and
tick bin tb uses default_rng([SEED, world index, tick bin index, i]).

  uv run --group analysis tools/est/phase3d2.py pilot|all|design|power|wm|thin|analyse|topup \
      [--panels N] [--workers 3]
"""
from __future__ import annotations

import argparse
import copy
import fcntl
import filecmp
import json
import os
import re
import shutil
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

import explore as E
import prereg3d as C
import primary3d as P
import worlds3d2 as W

OUT = Path(os.environ.get("PHASE3D2_OUT", E.data_dir() / "phase3d2"))  # override for smoke runs
SEED = 20261031
PILOT_PANELS = 10
N_STUDY = {"W-null": 200, "W0": 200, "W0-psi3": 100, "WB": 100, "WM": 100}
CHUNK_PANELS = 2
DEFAULT_PM = P.PMParams(3, 16, 128, 50)  # the pilot's design (test_3d1's PM)
PM_CHUNK = P.CHUNK  # run_primary's chunk (a smaller one saved no memory and cost 1.6x the time, measured)
THIN_PANELS = 50
THIN_T1_WORKERS = 2
THIN_SHARES = (0.46, 0.26, 0.17)
THIN_H = 0.17
THIN_KAPPAS = (1.0, 4.0, 16.0)  # CHOICE: thinning.py has only kappa = 4 (default); a grid around it
THIN_SETTINGS = (("full",),) + tuple(("random", s) for s in THIN_SHARES) + tuple(
    ("router", s, THIN_H, k) for s in THIN_SHARES for k in THIN_KAPPAS)
# Settings run in groups (one job each, memory); every group regenerates the SAME market days from the
# panel's market seed and thins with the setting's own seed, so all settings stay paired.
THIN_GROUPS = ((0, 1, 2, 3), (4, 5, 6), (7, 8, 9), (10, 11, 12))
DESIGN_WORLDS = ("W-null", "W0", "WB")  # stage A: every P-M design (the design rule's worlds)
POWER_WORLDS = ("W0-psi3",)  # stage B: the chosen design only
STUDY_WORLDS = DESIGN_WORLDS + POWER_WORLDS  # what all, topup and the pilot run (decided 1 Nov 2026)
LATER_WORLDS = ("WM",)  # the wm command, at the chosen design


def setting_name(s: tuple) -> str:
    if s[0] == "full":
        return "full 1.0"
    if s[0] == "random":
        return f"random {s[1]:g}"
    return f"router {s[1]:g} kappa {s[3]:g}"


# ---- one tick bin of one panel ----------------------------------------------------------------------------


def _f(x):
    x = np.asarray(x, dtype=float)
    if x.ndim == 0:
        return float(x) if np.isfinite(x) else None
    return [float(v) if np.isfinite(v) else None for v in x.ravel()]


def test_record(t: P.TestResult) -> dict:
    """The compact record of one TestResult: p, conditions, range, the statistics the report uses."""
    rec = {"p": float(t.p), "conditions": {k: bool(v) for k, v in t.conditions.items()},
           "range_used": list(t.range_used) if t.range_used else None, "range_rule": t.range_rule}
    st = t.statistic
    if t.name == "P-V":
        rec.update(p_S1=st.get("p_S1"), p_S2=st.get("p_S2"),
                   slope_S1=_f(st["slope_S1"]["estimate"]) if st else None,
                   slope_S2=_f(st["slope_S2"]["estimate"]) if st else None,
                   errors={k: _f(v) for k, v in st["loo_errors"].items() if k in ("M1", "M2", "M3")} if st else None,
                   by_range={r: v["p"] for r, v in t.by_range.items()})
    elif t.name == "P-R":
        rec.update(a_star=_f(st["a_star_T1_T2"]["estimate"]) if st else None,
                   a_lo=_f(st["a_star_T1_T2"]["lo"]) if st else None, a_hi=_f(st["a_star_T1_T2"]["hi"]) if st else None,
                   by_range={r: v["p"] for r, v in t.by_range.items()},
                   psi_hat=_f(t.secondary["psi_hat"]["estimate"]) if t.secondary else None)
    elif t.name == "P-F":
        rec.update(slope=_f(st["pf"]["estimate"]) if st else None, by_range={r: v["p"] for r, v in t.by_range.items()})
    elif t.name == "P-M":
        rec.update(p_slope=st.get("p_slope"), p_top=st.get("p_top"),
                   slope=_f(st["slope"]["estimate"]) if st else None, top=_f(st["top"]["estimate"]) if st else None,
                   beta=_f(st["beta_by_bin"]["estimate"]) if st else None,
                   n_stocks=len(t.stocks_used.get("stocks", [])), stock_days=t.stocks_used.get("stock_days", 0))
    return rec


def run_tick_bin(flows: dict, tb: str, rng: np.random.Generator, designs, n_boot: int = C.N_BOOT,
                 log_time: dict | None = None) -> dict:
    """As run_primary does for one tick bin (unit, evaluability by the reach rule alone, ONLY evaluable tests
    computed by run_unit), with P-M evaluated at every design on the same draws as run_unit's (the generator
    state is copied before run_unit, whose first draws are P.draws)."""
    t0 = time.time()
    u = P.tick_bin_unit(flows, tb)
    t1 = time.time()
    ev0 = P.evaluability(u.Ns, u.rho, designs[0])
    rng_dr = copy.deepcopy(rng)
    range_tests = tuple(t for t in ("P-V", "P-R", "P-F") if ev0[t]["evaluable"])
    res = P.run_unit(u, designs[0], rng, tests=range_tests, n_boot=n_boot)
    t2 = time.time()
    S_, D_ = u.flow.Nsd.shape
    dr = P.draws(S_, D_, rng_dr, n_boot)
    _, by = P.select_range(u.Ns, u.rho)
    pm, done = {}, {}
    for d in designs:
        ev = P.evaluability(u.Ns, u.rho, d)["P-M"]
        rec = {"evaluable": ev["evaluable"], "reason": ev["reason"]}
        if ev["evaluable"]:
            # designs that differ only in min runs per bin and keep the same stock-days have identical moments,
            # draws and so identical results: the first one's is reused (exact)
            same = (d.n_bins, d.l_max, d.L, P.pm_moments(u.flow, u.rho, d)["used"].tobytes())
            if same not in done:
                tests: dict = {}
                P._pm_test(u, d, dr, PM_CHUNK, by, tests)
                done[same] = test_record(tests["P-M"])
            rec.update(done[same])
        pm[W.design_key(d)] = rec
    t3 = time.time()
    if log_time is not None:
        log_time.update(unit=t1 - t0, range_tests=t2 - t1, pm=t3 - t2)
    return {"tick_bin": tb, "n_stocks": S_, "Ns": u.Ns.tolist(), "rho": u.rho.tolist(),
            "evaluability": {t: ev0[t] for t in ("P-V", "P-R", "P-F")},
            "tests": {k: test_record(v) for k, v in res["tests"].items()}, "pm": pm}


# ---- jobs -------------------------------------------------------------------------------------------------


def panel_rng(world: str, tb: str, i: int, salt: int = 0) -> np.random.Generator:
    return np.random.default_rng([SEED, salt, W.WORLDS.index(world), C.TICK_BINS.index(tb), i])


def chunk_path(kind: str, world: str, tb: str, start: int):
    return OUT / kind / f"{world}_{tb}_{start:03d}.json"


def study_job(args) -> str:
    world, tb, start, count, design_keys = args
    path = chunk_path("study", world, tb, start)
    if path.exists():
        return str(path)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    designs = [W.design_by_key(k) for k in design_keys]
    recs = []
    for i in range(start, start + count):
        rng = panel_rng(world, tb, i)
        tg = time.time()
        flows = W.tick_bin_flows(world, tb, rng)
        gen_s = time.time() - tg
        tm: dict = {}
        r = run_tick_bin(flows, tb, rng, designs, log_time=tm)
        r.update(panel=i, seconds={"generate": gen_s, **tm})
        recs.append(r)
        del flows
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(recs))
    tmp.rename(path)
    return str(path)


def thin_path(tb: str, i: int, g: int):
    return OUT / "thin" / f"WM_{tb}_{i:03d}_g{g}.json"


def thin_job(args) -> str:
    """One panel's tick bin under one group of thinning settings; the market days come from the panel's
    market seed (identical in every group), each setting's thinning from its own seed."""
    tb, i, g, design_key = args
    path = thin_path(tb, i, g)
    if path.exists():
        return str(path)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    ks = THIN_GROUPS[g]
    t0 = time.time()
    per_cell = {}
    for c in W.cells_of(tb):  # one market generator per cell, so a cell's markets never depend on the group
        rng = np.random.default_rng([SEED, 7, C.TICK_BINS.index(tb), i, sorted(C.CONFIRMATORY).index(c)])
        thin_rngs = [np.random.default_rng([SEED, 8, C.TICK_BINS.index(tb), i, k, sorted(C.CONFIRMATORY).index(c)])
                     for k in range(len(THIN_SETTINGS))]
        per_cell[c] = W.wm_cell_flows(c, rng, THIN_SETTINGS, thin_rngs=thin_rngs, keep=ks)
    gen_s = time.time() - t0
    recs = []
    for j, k in enumerate(ks):
        flows = {c: W.to_flow(per_cell[c][j]) for c in per_cell}
        for c in per_cell:
            per_cell[c][j] = None
        tm: dict = {}
        r = run_tick_bin(flows, tb, panel_rng("WM", tb, i, salt=2 + k), [W.design_by_key(design_key)], log_time=tm)
        r.update(setting=setting_name(THIN_SETTINGS[k]), panel=i, seconds=tm)
        recs.append(r)
        del flows
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps({"generate_s": gen_s, "total_s": time.time() - t0, "records": recs}))
    tmp.rename(path)
    return str(path)


def study_jobs(n_from: dict, n_to: dict, keys, worlds=W.WORLDS) -> list:
    js = [(w, tb, s, min(CHUNK_PANELS, n_to[w] - s), tuple(keys)) for w in worlds for tb in C.TICK_BINS
          for s in range(n_from.get(w, 0), n_to[w], CHUNK_PANELS)]
    return sorted(js, key=lambda j: (j[1] != "Q-T1", j[0]))  # the largest tick bin first


def run_jobs(fn, js, workers: int, done, log=print) -> None:
    todo = [j for j in js if not done(j)]
    log(f"{len(js)} jobs, {len(todo)} to run on {workers} workers", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for k, _ in enumerate(ex.map(fn, todo), 1):
            if k % 10 == 0 or k == len(todo):
                log(f"  {k}/{len(todo)}, {(time.time() - t0) / 60:.1f} min", flush=True)


# ---- one instance at a time, integrity of the saved chunks, trimming to the retained designs -------------

LOCK_NAME = "phase3d2.lock"


def acquire_lock(out: Path):
    """An exclusive, non-blocking flock on out/phase3d2.lock, held while the returned file stays open (keep a
    reference). The OS drops it when the process ends, so a crash or kill leaves no stale lock. A second
    instance on the same directory exits with the holder's pid and start time."""
    out.mkdir(parents=True, exist_ok=True)
    f = open(out / LOCK_NAME, "a+")
    try:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        f.seek(0)
        held = f.read().strip() or "pid unknown"
        f.close()
        raise SystemExit(f"another phase3d2 run holds {out / LOCK_NAME} ({held}); not starting") from None
    f.truncate(0)
    f.write(f"pid {os.getpid()}, started {time.strftime('%Y-%m-%d %H:%M:%S')}, argv {' '.join(sys.argv[1:])}\n")
    f.flush()
    return f


CHUNK_RE = re.compile(r"^(?P<world>.+)_(?P<tb>Q-T[123])_(?P<start>\d{3})\.json$")
REC_KEYS = ("tick_bin", "n_stocks", "Ns", "rho", "evaluability", "tests", "pm", "panel", "seconds")


def _p_ok(p) -> bool:
    return isinstance(p, (int, float)) and 0.0 <= p <= 1.0


def verify_chunk(path: Path) -> list[str]:
    """Problems with one saved study chunk (empty when it is complete): its name is a job of STUDY_WORLDS or
    LATER_WORLDS, it parses, it holds panels start, start + 1, ... of its tick bin (CHUNK_PANELS of them unless
    the world ends sooner), every record has every field, every evaluable range test a p in [0, 1], and P-M
    records for exactly the retained designs, all 18 of the old grid (design worlds, before trim) or one
    retained design (the other worlds)."""
    m = CHUNK_RE.match(path.name)
    if not m or m["world"] not in STUDY_WORLDS + LATER_WORLDS:
        return ["not a study job name"]
    world, tb, start = m["world"], m["tb"], int(m["start"])
    if start % CHUNK_PANELS:
        return [f"start {start} is not a chunk start"]
    try:
        recs = json.loads(path.read_text())
    except (ValueError, UnicodeDecodeError) as e:
        return [f"does not parse: {e}"]
    if not isinstance(recs, list) or not recs:
        return ["not a non-empty list of records"]
    probs = []
    # a chunk past N_STUDY (topup) may end at any panel count
    want = min(CHUNK_PANELS, N_STUDY[world] - start) if start < N_STUDY[world] else min(len(recs), CHUNK_PANELS)
    if len(recs) != want:
        probs.append(f"{len(recs)} records, expected {want}")
    retained = {W.design_key(d) for d in W.PM_DESIGNS}
    old = {W.design_key(d) for d in W.PM_DESIGNS_18}
    for j, r in enumerate(recs):
        if not isinstance(r, dict) or any(k not in r for k in REC_KEYS):
            probs.append(f"record {j}: missing fields")
            continue
        if r["panel"] != start + j or r["tick_bin"] != tb:
            probs.append(f"record {j}: panel {r['panel']} {r['tick_bin']}, expected {start + j} {tb}")
        if len(r["Ns"]) != r["n_stocks"] or len(r["rho"]) != r["n_stocks"]:
            probs.append(f"record {j}: Ns or rho not n_stocks long")
        for t in ("P-V", "P-R", "P-F"):
            if r["evaluability"].get(t, {}).get("evaluable") and not _p_ok(r["tests"].get(t, {}).get("p")):
                probs.append(f"record {j}: {t} evaluable without a p in [0, 1]")
        keys = set(r["pm"])
        ok = keys in (retained, old) if world in DESIGN_WORLDS else (len(keys) == 1 and keys <= retained)
        if not ok:
            probs.append(f"record {j}: P-M designs {sorted(keys)}")
        for k, v in r["pm"].items():
            if "evaluable" not in v or (v["evaluable"] and (not _p_ok(v.get("p")) or "conditions" not in v)):
                probs.append(f"record {j}: P-M {k} incomplete")
    return probs


def check_study(out: Path = OUT, log=print) -> dict:
    """Verifies every file in out/study and moves each bad one, each duplicate (a panel already held by an
    earlier file of the same world and tick bin) and each leftover temporary file to
    out/quarantine/<time>/study/, with the reason in that folder's reasons.json. Returns the counts."""
    d = out / "study"
    files = sorted(d.iterdir()) if d.exists() else []
    bad: dict[str, list[str]] = {}
    seen: dict[tuple, str] = {}
    n_ok = 0
    for f in files:
        if f.suffix != ".json":
            bad[f.name] = ["leftover temporary or foreign file"]
            continue
        probs = verify_chunk(f)
        if not probs:
            m = CHUNK_RE.match(f.name)
            for r in json.loads(f.read_text()):
                key = (m["world"], m["tb"], r["panel"])
                if key in seen:
                    probs.append(f"panel {r['panel']} already in {seen[key]} (duplicate)")
                seen.setdefault(key, f.name)
        if probs:
            bad[f.name] = probs
        else:
            n_ok += 1
    q = None
    if bad:
        q = out / "quarantine" / time.strftime("%Y%m%d-%H%M%S") / "study"
        q.mkdir(parents=True)
        for name in bad:
            shutil.move(str(d / name), str(q / name))
        (q / "reasons.json").write_text(json.dumps(bad, indent=1))
    res = {"files": len(files), "complete": n_ok, "quarantined": len(bad),
           "duplicates": sum(any("duplicate" in p for p in v) for v in bad.values()),
           "quarantine": str(q) if q else None, "reasons": bad}
    log(f"check: {res['files']} files, {res['complete']} complete, {res['quarantined']} quarantined "
        f"({res['duplicates']} duplicates){', in ' + res['quarantine'] if q else ''}", flush=True)
    return res


def trim_designs(out: Path = OUT, log=print) -> int:
    """Keeps only worlds3d2.PM_DESIGNS in the design worlds' chunks saved with the 18-design grid (decided
    1 Nov 2026). P-M at one design never depends on the other designs in the job (same draws; reuse only
    between designs with the same moments, which is exact), so the retained records are what a 6-design job
    gives. The original moves unchanged to out/study_18designs/ first; nothing is recomputed."""
    keep = [W.design_key(d) for d in W.PM_DESIGNS]
    arch = out / "study_18designs"
    n = 0
    for f in sorted((out / "study").glob("*.json")):
        m = CHUNK_RE.match(f.name)
        if not m or m["world"] not in DESIGN_WORLDS:
            continue
        recs = json.loads(f.read_text())
        if all(set(r["pm"]) == set(keep) for r in recs):
            continue
        for r in recs:
            r["pm"] = {k: r["pm"][k] for k in keep}
        arch.mkdir(parents=True, exist_ok=True)
        if not (arch / f.name).exists():
            shutil.copy2(f, arch / f.name)
        elif not filecmp.cmp(f, arch / f.name, shallow=False):  # a rerun after a stop between copy and rename
            raise SystemExit(f"{arch / f.name} differs from {f}; not overwriting an archived original")
        tmp = f.with_suffix(".tmp")
        tmp.write_text(json.dumps(recs))
        tmp.rename(f)
        n += 1
    log(f"trim: {n} chunks kept to {len(keep)} designs; originals in {arch}", flush=True)
    return n


# ---- pilot ------------------------------------------------------------------------------------------------


def pilot_job(args) -> dict:
    world, tb, i = args
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    rng = panel_rng(world, tb, 10_000 + i)
    t0 = time.time()
    flows = W.tick_bin_flows(world, tb, rng)
    gen_s = time.time() - t0
    tm: dict = {}
    r = run_tick_bin(flows, tb, rng, [DEFAULT_PM], log_time=tm)
    t1 = time.time()
    # one design's P-M alone, to project the cost of every design
    return {"world": world, "tick_bin": tb, "panel": i, "generate": gen_s, **tm, "total": t1 - t0,
            "pm_evaluable": r["pm"][W.design_key(DEFAULT_PM)]["evaluable"], "result": r}


def pilot(workers: int, n: int = PILOT_PANELS, log=print) -> None:
    js = [(w, tb, i) for w in STUDY_WORLDS for tb in C.TICK_BINS for i in range(n)]
    js.sort(key=lambda j: (j[1] != "Q-T1", j[0]))
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(pilot_job, js))
    wall = time.time() - t0
    (OUT / "report").mkdir(parents=True, exist_ok=True)
    (OUT / "report" / "pilot.json").write_text(json.dumps(res))
    log(f"pilot: {len(js)} panel-tick-bins in {wall / 60:.1f} min on {workers} workers", flush=True)
    proj = 0.0
    lines = ["| world | tick bin | generate s | unit s | P-V/P-R/P-F s | P-M (1 design) s | total s |",
             "|---|---|---|---|---|---|---|"]
    for w in STUDY_WORLDS:
        for tb in C.TICK_BINS:
            rs = [r for r in res if r["world"] == w and r["tick_bin"] == tb]
            m = {k: float(np.mean([r[k] for r in rs])) for k in ("generate", "unit", "range_tests", "pm", "total")}
            lines.append(f"| {w} | {tb} | {m['generate']:.1f} | {m['unit']:.1f} | {m['range_tests']:.1f} | "
                         f"{m['pm']:.1f} | {m['total']:.1f} |")
            extra = len(W.PM_DESIGNS) - 1 if w in DESIGN_WORLDS else 0
            proj += (m["total"] + extra * m["pm"]) * N_STUDY[w]
    lines.append(f"\nProjected study (all: {len(W.PM_DESIGNS)} designs in {', '.join(DESIGN_WORLDS)}, the chosen "
                 f"one in {', '.join(POWER_WORLDS)}, N_STUDY panels): {proj / 3600:.1f} worker-hours, "
                 f"{proj / 3600 / workers:.1f} h wall on {workers} workers (P-M cost per design taken as the "
                 f"default design's; WM and thinning not included).")
    (OUT / "report" / "pilot.md").write_text("\n".join(lines) + "\n")
    log("\n".join(lines), flush=True)


def stage_a(workers: int, n_to: dict) -> None:
    keys = [W.design_key(d) for d in W.PM_DESIGNS]
    run_jobs(study_job, study_jobs({}, n_to, keys, DESIGN_WORLDS), workers,
             lambda j: chunk_path("study", *j[:3]).exists())


def stage_b(workers: int, n_to: dict, worlds=POWER_WORLDS) -> None:
    import analyse3d2 as A

    key = A.chosen_design(write=True)
    print(f"P-M design: {key}", flush=True)
    run_jobs(study_job, study_jobs({}, n_to, [key], worlds), workers,
             lambda j: chunk_path("study", *j[:3]).exists())


def thin(workers: int, n: int) -> None:
    import analyse3d2 as A

    key = A.chosen_design()
    js = [(tb, i, g, key) for tb in C.TICK_BINS for i in range(n) for g in range(len(THIN_GROUPS))]
    done = lambda j: thin_path(*j[:3]).exists()  # noqa: E731
    # a Q-T1 thinning job peaks near 1.0 GB (measured), so at most 2 run at once (RAM < 3 GB)
    run_jobs(thin_job, [j for j in js if j[0] == "Q-T1"], min(workers, THIN_T1_WORKERS), done)
    run_jobs(thin_job, [j for j in js if j[0] != "Q-T1"], workers, done)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("step", choices=["pilot", "all", "design", "power", "wm", "thin", "topup", "analyse", "check",
                                     "trim"])
    ap.add_argument("--panels", type=int, default=None)
    ap.add_argument("--workers", type=int, default=min(3, os.cpu_count() or 1))
    a = ap.parse_args(argv)
    if a.workers > 3:
        raise SystemExit("at most 3 workers")
    lock = acquire_lock(OUT)  # noqa: F841 (held until the process ends)
    t0 = time.time()
    if a.step == "check":
        check_study(OUT)
    if a.step == "trim":
        trim_designs(OUT)
    n_to = dict(N_STUDY) if a.panels is None else {w: min(a.panels, N_STUDY[w]) for w in W.WORLDS}
    if a.step == "topup":
        if a.panels is None:
            raise SystemExit("topup needs --panels")
        n_to = {w: max(a.panels, N_STUDY[w]) for w in STUDY_WORLDS}
    if a.step == "pilot":
        pilot(a.workers, a.panels or PILOT_PANELS)
    elif a.step in ("all", "design", "topup"):
        stage_a(a.workers, n_to)
    if a.step in ("all", "power", "topup"):
        stage_b(a.workers, n_to)
    if a.step == "wm":
        stage_b(a.workers, n_to, LATER_WORLDS)
    if a.step == "thin":
        thin(a.workers, min(a.panels, THIN_PANELS) if a.panels else THIN_PANELS)
    if a.step in ("all", "analyse", "topup"):
        import analyse3d2 as A

        A.analyse()
    print(f"{a.step}: {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
