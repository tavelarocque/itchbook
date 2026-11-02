"""Phase 4 reproduction check: the confirmatory run recomputed on the confirmation data and compared bit for
bit with data/phase4/results.json. Not a second confirmatory analysis: it cannot change any reported result,
never opens data/phase4/ for writing, and writes only to data/phase4_repro/.

From phase4 it imports only the confirmation-only loader (build_flows, split 'confirm') and the frozen design
constants (SEED, FROZEN_PM, NASDAQ_CELLS, PLAN_CONFIRM_COUNTS). The ceilings are prereg3d's, applied by
primary3d.make_unit. It calls the frozen primary3d.run_primary once over the nine NASDAQ-listed cells with
FROZEN_PM, default_rng(SEED) and its defaults (B = 200, four surrogate draws), as phase4.run_family did, without
its observe() recorder (which calls only frozen helpers that take no draws). Nothing frozen is wrapped or
patched.

Comparison. Every leaf of run_primary's return (dataclasses as dicts, tuples and arrays as lists, tuple keys
written as phase4 wrote them) against results.json["result"]: the same keys and lengths; None, bools, ints and
strings equal; floats equal by float.hex (a NaN stored as null, an infinity as "inf" or "-inf"). Holm is
recomputed with the frozen primary3d.holm from the reproduced family p-values and compared with the stored and
the reproduced Holm records, and passed with rejected and every condition. The run's record (seed, P-M
design, B, surrogate draws), the runs files' SHA-256 and each cell's symbols and counts are checked too.
The 'methods' record of results.json came from observe() and is not a field of run_primary; it is not compared.

    uv run --group analysis python tools/est/phase4_repro.py    # exit 0 only if nothing differs
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import explore as E
import prereg3d as C
import primary3d as P
from phase4 import FROZEN_PM, NASDAQ_CELLS, PLAN_CONFIRM_COUNTS, SEED, build_flows

REPO = Path(__file__).resolve().parents[2]
STORED = REPO / "data" / "phase4" / "results.json"
OUT = REPO / "data" / "phase4_repro"


def key(k) -> str:
    """A dict key as results.json holds it."""
    if isinstance(k, tuple):
        if len(k) == 2 and all(isinstance(x, (int, np.integer)) for x in k):
            return f"[{k[0]},{k[1]}]"
        return "|".join(str(x) for x in k)
    return str(k)


def show(x) -> str:
    if isinstance(x, (float, np.floating)):
        return f"{float(x)!r} ({float(x).hex()})"
    return repr(x)


class Compare:
    """Leaf-by-leaf comparison of a reproduced object with its stored JSON; every leaf counted once."""

    def __init__(self):
        self.n, self.bad, self.by_section = 0, [], {}

    def leaf(self, path: str, ok: bool, stored, reproduced) -> None:
        self.n += 1
        sec = path.split(".")[0]
        self.by_section[sec] = self.by_section.get(sec, 0) + 1
        if not ok:
            self.bad.append({"path": path, "stored": stored, "reproduced": show(reproduced)})

    def __call__(self, a, b, path: str) -> None:
        if dataclasses.is_dataclass(a) and not isinstance(a, type):
            a = {f.name: getattr(a, f.name) for f in dataclasses.fields(a)}
        if isinstance(a, np.ndarray):
            a = a.tolist()
        if isinstance(a, dict):
            ka = {key(k): v for k, v in a.items()}
            if len(ka) != len(a):
                raise RuntimeError(f"{path}: two keys of the reproduced dict write the same JSON key")
            if not isinstance(b, dict):
                return self.leaf(path, False, b, f"dict of {len(ka)} keys")
            if not ka and not b:
                return self.leaf(path, True, b, a)
            for k in list(ka) + [k for k in b if k not in ka]:
                if k not in ka or k not in b:
                    self.leaf(f"{path}.{k}", False, b.get(k, "<missing>"), ka.get(k, "<missing>"))
                else:
                    self(ka[k], b[k], f"{path}.{k}")
            return
        if isinstance(a, (list, tuple)):
            if not isinstance(b, list) or len(a) != len(b):
                return self.leaf(path, False, b, a)
            if not a:
                return self.leaf(path, True, b, a)
            for i, (x, y) in enumerate(zip(a, b)):
                self(x, y, f"{path}[{i}]")
            return
        if a is None:
            ok = b is None
        elif isinstance(a, (bool, np.bool_)):
            ok = isinstance(b, bool) and bool(a) == b
        elif isinstance(a, (int, np.integer)):  # b is a JSON int, or a numpy int in the in-memory Holm check
            ok = type(b) is int or isinstance(b, np.integer)
            ok = ok and int(a) == int(b)
        elif isinstance(a, (float, np.floating)):  # likewise a JSON float or a numpy float
            f = float(a)
            if math.isnan(f):
                ok = b is None
            elif math.isinf(f):
                ok = b == ("inf" if f > 0 else "-inf")
            else:
                ok = (type(b) is float or isinstance(b, np.floating)) and f.hex() == float(b).hex()
        elif isinstance(a, str):
            ok = isinstance(b, str) and a == b
        else:
            raise TypeError(f"{path}: unexpected type {type(a)}")
        self.leaf(path, ok, b, a)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(REPO), *args], check=True, capture_output=True, text=True).stdout.strip()


def cell_counts(flows: dict, syms: dict) -> dict:
    return {c: {"stocks": f.n_stocks, "stock_days": len(f.lengths), "runs": int(f.lengths.sum()),
                "invalid_mids": int((~f.ok).sum()), "N_s_min": int(f.Nsd.min()), "symbols": syms[c]}
            for c, f in flows.items()}


def holm_check(res: dict, stored: dict, cmp: Compare) -> list[dict]:
    """Holm from the reproduced family p-values (frozen primary3d.holm), against the stored and the reproduced
    records; passed against rejected and every condition."""
    ps = {k: res["by_tick_bin"][k[1]]["tests"][k[0]].p for k in res["family"]}
    h = P.holm(ps)
    rows = []
    for k, rec in h.items():
        cmp(rec, stored["holm"][key(k)], f"holm_recomputed.{key(k)}")
        cmp(rec, res["holm"][k], f"holm_recomputed_vs_reproduced.{key(k)}")
        t = res["by_tick_bin"][k[1]]["tests"][k[0]]
        cmp(bool(rec["rejected"] and t.all_conditions),
            stored["by_tick_bin"][k[1]]["tests"][k[0]]["passed"], f"holm_recomputed.{key(k)}.passed")
        rows.append({"pair": key(k), "p": float(rec["p"]), "rank": int(rec["rank"]),
                     "threshold": float(rec["threshold"]), "p_adj": float(rec["p_adj"]),
                     "rejected": bool(rec["rejected"]), "passed": bool(rec["rejected"] and t.all_conditions)})
    return sorted(rows, key=lambda r: r["rank"])


def main() -> int:
    stored = json.loads(STORED.read_text())
    head, status = git("rev-parse", "HEAD"), git("status", "--porcelain")
    t0 = time.time()
    data = E.data_dir()
    days, _ = E.read_universe()
    rh = {d: sha256(data / "runs" / f"{d}.runs.npz") for d in days}
    flows, syms = build_flows("confirm", cells=NASDAQ_CELLS, data=data)
    t1 = time.time()
    res = P.run_primary(flows, FROZEN_PM, np.random.default_rng(SEED))
    t2 = time.time()

    cmp = Compare()
    cmp(res, stored["result"], "result")
    cmp({"seed": SEED, "pm": FROZEN_PM, "n_boot": C.N_BOOT, "n_sur": P.N_SUR},
        {k: stored[k] for k in ("seed", "pm", "n_boot", "n_sur")}, "run")
    cmp(rh, stored["runs_npz_sha256"], "inputs_sha256")
    counts = cell_counts(flows, syms)
    cmp(counts, {c: {k: s[k] for k in counts[c]} if c in counts else s for c, s in stored["stock_counts"].items()},
        "inputs_cells")
    cmp({c: s["stocks"] for c, s in counts.items()}, PLAN_CONFIRM_COUNTS, "inputs_plan_counts")
    holm_rows = holm_check(res, stored["result"], cmp)

    OUT.mkdir(parents=True, exist_ok=True)
    summary = {"identical": not cmp.bad, "fields_compared": cmp.n, "mismatches": len(cmp.bad),
               "fields_by_section": cmp.by_section, "git_head": head, "git_status": status,
               "phase4_repro_py_sha256": sha256(Path(__file__)), "seed": SEED, "pm": dataclasses.asdict(FROZEN_PM),
               "n_boot": C.N_BOOT, "n_sur": P.N_SUR, "python": sys.version.split()[0], "numpy": np.__version__,
               "load_s": t1 - t0, "run_primary_s": t2 - t1, "wall_s": time.time() - t0,
               "stored_wall_s": stored["wall_s"], "stored_git_head": stored["preflight"]["git_head"],
               "holm_recomputed": holm_rows, "mismatch_list": cmp.bad}
    (OUT / "compare.json").write_text(json.dumps(summary, indent=1))
    L = ["# Phase 4 reproduction check", "",
         f"Result: **{'bit-identical' if not cmp.bad else 'MISMATCH'}**. {cmp.n} fields compared, "
         f"{len(cmp.bad)} mismatches.", "",
         "A reproduction of the confirmatory run, not a second confirmatory analysis; it changes no reported "
         "result. tools/est/phase4_repro.py loaded the confirmation stocks with phase4.build_flows('confirm') and "
         "called the frozen primary3d.run_primary once over the nine NASDAQ-listed cells with the frozen P-M "
         "design and seed, with no wrapper, then compared every field it returns with data/phase4/results.json "
         "(floats by float.hex).", "",
         f"- git HEAD: {head} (stored run: {stored['preflight']['git_head']}); working tree: "
         f"{status or 'clean'}",
         f"- tools/est/phase4_repro.py SHA-256 {summary['phase4_repro_py_sha256']}",
         f"- seed {SEED}; P-M design {summary['pm']}; B = {C.N_BOOT}; surrogate draws {P.N_SUR}",
         f"- python {summary['python']}, numpy {summary['numpy']}",
         f"- runtime: {summary['wall_s']:.1f} s (loading and hashing {summary['load_s']:.1f} s, run_primary "
         f"{summary['run_primary_s']:.1f} s); the stored run took {stored['wall_s']:.1f} s", "",
         "## Fields compared", "", "| section | fields |", "|---|---|"]
    L += [f"| {k} | {v} |" for k, v in cmp.by_section.items()]
    L += ["", "result = every leaf of run_primary's return (evaluability, reach, p-values, Holm, sub-conditions, "
          "estimates, intervals, P-R secondaries, descriptive P-K); run = seed, P-M design, B, surrogate draws; "
          "inputs_* = runs files' SHA-256, per-cell symbols and counts, plan stock counts; holm_recomputed* = "
          "frozen primary3d.holm on the reproduced family p-values against the stored and reproduced records, "
          "and passed. The 'methods' record (from phase4's observe() recorder) is not a field of run_primary "
          "and is not compared.", "",
          "## Holm recomputed from the reproduced p-values", "",
          "| rank | pair | p | threshold | adjusted p | rejects | passed |", "|---|---|---|---|---|---|---|"]
    L += [f"| {r['rank']} | {r['pair']} | {r['p']!r} | {r['threshold']!r} | {r['p_adj']!r} | "
          f"{'yes' if r['rejected'] else 'no'} | {'yes' if r['passed'] else 'no'} |" for r in holm_rows]
    if cmp.bad:
        L += ["", "## Mismatches", "", "| path | stored | reproduced |", "|---|---|---|"]
        L += [f"| {m['path']} | {m['stored']!r} | {m['reproduced']} |" for m in cmp.bad]
    (OUT / "compare.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    return 0 if not cmp.bad else 1


if __name__ == "__main__":
    sys.exit(main())
