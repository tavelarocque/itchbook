"""Phase 3d-2 analysis of the saved panels (phase3d2 study and thin chunks): the P-M design rule, rejection
rates before and after the one Holm family, family-wise error in W-null, ranges used, pairs that left the
family, bin occupancy of the real exploration flow, and the thinning changes. Writes
~/CODE/itch-data/phase3d2/report/*.md and results.json.

Definitions (CHOICE, fixed with the rule before any study result):
- "rejected before Holm": the pair is evaluable and its p <= 0.05;
- "rejected after Holm": primary3d.combine's one family at the chosen P-M design rejects it. Since 2 Nov 2026
  the family is Q-T1 and Q-T3 only (up to 8 pairs); Q-T2 is descriptive (before-Holm rates only). Beside it,
  for comparison, the former family of all three tick bins (up to 12 pairs, 11 here as P-M never reaches in
  Q-T2), recomputed from the same saved p-values with primary3d.step_down;
- "passed": rejected after Holm AND every sub-condition (primary3d.set_decisions); for P-R this is the
  equivalence established;
- a pair that left the family counts as not rejected; rates are over all panels, MC SE sqrt(r (1 - r) / n);
- the P-M design rule's rate is the unadjusted pass: p <= 0.05 and both conditions (worlds3d2.choose_design).
"""
from __future__ import annotations

import json

import numpy as np

import phase3d2 as D
import prereg3d as C
import primary3d as P
import worlds3d2 as W

REP = D.OUT / "report"
ALPHA = 0.05


def load_study(world: str) -> dict:
    """{tick bin: [record per panel, ordered by panel]}."""
    out = {}
    for tb in C.TICK_BINS:
        recs = []
        for f in sorted((D.OUT / "study").glob(f"{world}_{tb}_*.json")):
            recs += json.loads(f.read_text())
        out[tb] = sorted(recs, key=lambda r: r["panel"])
    return out


def common_panels(st: dict) -> list[int]:
    """Panels present in every tick bin (a panel is one draw of all three tick bins)."""
    sets = [{r["panel"] for r in st[tb]} for tb in C.TICK_BINS]
    return sorted(set.intersection(*sets))


def pm_unadjusted(rec: dict, key: str) -> bool:
    r = rec["pm"].get(key)
    return bool(r and r["evaluable"] and r["p"] <= ALPHA and all(r["conditions"].values()))


def design_rates(worlds=(*W.SIZE_WORLDS, W.POWER_WORLD)) -> tuple[dict, dict]:
    rates, n = {}, {}
    for w in worlds:
        st = load_study(w)
        pan = common_panels(st)
        n[w] = len(pan)
        for tb in C.TICK_BINS:
            by = {r["panel"]: r for r in st[tb]}
            for d in W.PM_DESIGNS:
                k = W.design_key(d)
                rates[(w, k, tb)] = float(np.mean([pm_unadjusted(by[i], k) for i in pan])) if pan else np.nan
    return rates, n


def chosen_design(write: bool = False) -> str:
    """The rule's design (its fallback, flagged, when none is admissible); cached in report/design.json."""
    path = REP / "design.json"
    if path.exists() and not write:
        return json.loads(path.read_text())["used"]
    rates, n = design_rates()
    out = W.choose_design(rates, n)
    out["used"] = out["chosen"] or out["fallback"]
    out["n_panels"] = n
    out["rates"] = {f"{w}|{k}|{tb}": v for (w, k, tb), v in rates.items()}
    REP.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1))
    return out["used"]


def family(rec_by_tb: dict, key: str) -> dict:
    """primary3d.combine over the panel's tick bins at P-M design `key`: evaluability from the reach rule as
    saved (P-V, P-R, P-F) and from the design's P-M evaluability; p-values of the evaluable pairs of Q-T1 and
    Q-T3 only (rejected, passed; None in Q-T2, which is descriptive). rejected_11, passed_11: the former
    family of every tick bin given, from the same p-values (primary3d.step_down), for comparison only."""
    ev, pv, recs = {}, {}, {}
    for tb, r in rec_by_tb.items():
        pm = r["pm"][key]
        ev[tb] = {**{t: {"evaluable": r["evaluability"][t]["evaluable"], "reason": r["evaluability"][t]["reason"]}
                     for t in ("P-V", "P-R", "P-F")},
                  "P-M": {"evaluable": pm["evaluable"], "reason": pm["reason"]}}
        ev[tb] = {t: ev[tb][t] for t in C.PRIMARY_TESTS}
        recs[tb] = {**r["tests"], "P-M": pm}
        for t in C.PRIMARY_TESTS:
            if ev[tb][t]["evaluable"]:
                pv[(t, tb)] = recs[tb][t]["p"]
    fam = P.combine(ev, {k: v for k, v in pv.items() if k[1] in C.FAMILY_TICK_BINS})
    former = P.step_down(pv, ALPHA)
    dec = {}
    for t in C.PRIMARY_TESTS:
        for tb in rec_by_tb:
            h, h11 = fam["holm"].get((t, tb)), former.get((t, tb))
            rec = recs[tb].get(t, {})
            conds = rec.get("conditions", {})
            ok = bool(conds) and all(conds.values())
            in_fam = tb in C.FAMILY_TICK_BINS
            dec[(t, tb)] = {"evaluable": ev[tb][t]["evaluable"],
                            "raw": bool(h11 is not None and rec["p"] <= ALPHA),
                            "rejected": bool(h and h["rejected"]) if in_fam else None,
                            "passed": bool(h and h["rejected"] and ok) if in_fam else None,
                            "rejected_11": bool(h11 and h11["rejected"]),
                            "passed_11": bool(h11 and h11["rejected"] and ok),
                            "range": tuple(rec["range_used"]) if rec.get("range_used") else None}
    return {"decisions": dec, "family_size": len(fam["family"]), "family_size_11": len(former)}


def rate(x) -> tuple[float, float, int]:
    x = np.asarray(x, dtype=float)
    n = len(x)
    r = float(x.mean()) if n else np.nan
    return r, float(np.sqrt(r * (1 - r) / n)) if n else np.nan, n


def fmt_rate(r) -> str:
    return "-" if not np.isfinite(r[0]) else f"{r[0]:.3f} ({r[1]:.3f})"


def correct_rate(world: str, test: str, decs: list, family: str = "") -> tuple | None:
    """Share of panels with the world's correct outcome; family "" (Q-T1 and Q-T3) or "_11" (former). None
    when unlabelled or descriptive (Q-T2 in the family)."""
    lab = W.CORRECT[world][test]
    rej, pas = "rejected" + family, "passed" + family
    if lab is None or any(d[rej] is None for d in decs):
        return None
    if lab == "no reject":
        return rate([not d[rej] for d in decs])
    if lab == "reject":
        return rate([d[pas] for d in decs])
    if lab == "pass":
        return rate([d[pas] for d in decs])
    return rate([not d[pas] for d in decs])  # "fail"


def world_table(world: str, key: str) -> dict:
    st = load_study(world)
    pan = common_panels(st)
    by = {tb: {r["panel"]: r for r in st[tb]} for tb in C.TICK_BINS}
    fams = [family({tb: by[tb][i] for tb in C.TICK_BINS}, key) for i in pan]
    out = {"panels": len(pan), "rows": {}, "fwer": None}
    for t in C.PRIMARY_TESTS:
        for tb in C.TICK_BINS:
            decs = [f["decisions"][(t, tb)] for f in fams]
            rng = [d["range"] for d in decs if d["evaluable"]]
            in_fam = tb in C.FAMILY_TICK_BINS
            out["rows"][(t, tb)] = {
                "evaluable": rate([d["evaluable"] for d in decs]),
                "raw": rate([d["raw"] for d in decs]),
                "holm": rate([d["rejected"] for d in decs]) if in_fam else None,
                "passed": rate([d["passed"] for d in decs]) if in_fam else None,
                "correct": correct_rate(world, t, decs),
                "holm_11": rate([d["rejected_11"] for d in decs]),
                "passed_11": rate([d["passed_11"] for d in decs]),
                "correct_11": correct_rate(world, t, decs, "_11"),
                "range_8_64": sum(r == (8, 64) for r in rng), "range_4_32": sum(r == (4, 32) for r in rng),
                "left": sum(not d["evaluable"] for d in decs),
            }
    nulls = [t for t in C.PRIMARY_TESTS if W.CORRECT[world][t] == "no reject"]
    if nulls:
        out["fwer"] = rate([any(f["decisions"][(t, tb)]["rejected"] for t in nulls for tb in C.FAMILY_TICK_BINS)
                            for f in fams])
        out["fwer_11"] = rate([any(f["decisions"][(t, tb)]["rejected_11"] for t in nulls for tb in C.TICK_BINS)
                               for f in fams])
        out["fwer_tests"] = nulls
    out["any_rejection"] = rate([any(d["rejected"] for d in f["decisions"].values()) for f in fams])
    out["any_rejection_11"] = rate([any(d["rejected_11"] for d in f["decisions"].values()) for f in fams])
    out["family_size"] = float(np.mean([f["family_size"] for f in fams])) if fams else np.nan
    out["family_size_11"] = float(np.mean([f["family_size_11"] for f in fams])) if fams else np.nan
    out["stats"] = stat_means(st, pan, key)
    return out


def stat_means(st: dict, pan: list, key: str) -> dict:
    """Mean (and SD over panels) of each test's statistic per tick bin."""
    out = {}
    for tb in C.TICK_BINS:
        by = {r["panel"]: r for r in st[tb]}
        rows = [by[i] for i in pan]
        vals = {"P-V slope S1": [r["tests"].get("P-V", {}).get("slope_S1") for r in rows],
                "P-V slope S2": [r["tests"].get("P-V", {}).get("slope_S2") for r in rows],
                "P-R a*(T1)": [(r["tests"].get("P-R", {}).get("a_star") or [None, None])[0] for r in rows],
                "P-R a*(T2)": [(r["tests"].get("P-R", {}).get("a_star") or [None, None])[1] for r in rows],
                "psi-hat": [r["tests"].get("P-R", {}).get("psi_hat") for r in rows],
                "P-F slope": [r["tests"].get("P-F", {}).get("slope") for r in rows],
                "P-M slope": [r["pm"][key].get("slope") for r in rows],
                "P-M top beta": [r["pm"][key].get("top") for r in rows]}
        out[tb] = {}
        for k, v in vals.items():
            x = np.array([np.nan if a is None else a for a in v], dtype=float)
            x = x[np.isfinite(x)]
            out[tb][k] = (float(x.mean()), float(x.std(ddof=1)) if len(x) > 1 else np.nan, len(x)) if len(x) else None
    return out


# ---- bin occupancy (real exploration flow, R2-clipped sizes, ties to the lower bin as pm_moments) -------------


def occupancy() -> list[dict]:
    rows = []
    for cell in C.CONFIRMATORY:
        fc = W.flow_cache(cell)
        lengths = fc["Nsd"].ravel()
        q = P.r2_clip(fc["q"], lengths, W.N_DAYS)
        starts = np.concatenate([[0], np.cumsum(lengths)[:-1]])
        rl = float(np.mean(fc["q"] % 100 == 0))
        for nb in C.PM_BINS_ALLOWED:
            mins, low_share, empty = [], [], 0
            for a, n in zip(starts, lengths, strict=True):
                x = q[a : a + n]
                b = np.searchsorted(np.quantile(x, np.linspace(0, 1, nb + 1)[1:-1]), x, side="left")
                cnt = np.bincount(b, minlength=nb)
                mins.append(cnt.min())
                low_share.append(cnt[0] / n)
                empty += int(cnt.min() == 0)
            mins = np.array(mins)
            rows.append({"cell": cell, "n_bins": nb, "stock_days": len(lengths), "round_lot_share": rl,
                         "lowest_bin_share_median": float(np.median(low_share)), "expected": 1 / nb,
                         "empty_bin_share": empty / len(lengths),
                         **{f"min_bin_ge_{m}": float(np.mean(mins >= m)) for m in (20, 50, 100)}})
    return rows


# ---- thinning ---------------------------------------------------------------------------------------------


def load_thin() -> dict:
    """{(setting, tick bin): {panel: record}}."""
    out: dict = {}
    for f in sorted((D.OUT / "thin").glob("WM_*.json")):
        for r in json.loads(f.read_text())["records"]:
            out.setdefault((r["setting"], r["tick_bin"]), {})[r["panel"]] = r
    return out


THIN_STATS = (("P-V", "slope_S1"), ("P-V", "slope_S2"), ("P-R", "a_star"), ("P-R", "psi_hat"), ("P-F", "slope"),
              ("P-M", "slope"), ("P-M", "top"))


def _stat(rec: dict, key: str, test: str, name: str):
    t = rec["pm"][key] if test == "P-M" else rec["tests"].get(test, {})
    v = t.get(name)
    if name == "a_star":
        v = None if v is None else v[1]  # a*(T2)
    return np.nan if v is None else float(v)


def thinning_table(key: str) -> list[dict]:
    """Per setting and tick bin: rejection rates after the setting's own Holm family (its tick bins of one
    panel), their change against share 1.0 on the same panels, and the paired mean change of each statistic."""
    th = load_thin()
    if not th:
        return []
    names = [D.setting_name(s) for s in D.THIN_SETTINGS]

    def fams(s):
        pans = set.intersection(*[set(th.get((s, tb), {})) for tb in C.TICK_BINS])
        return {i: family({tb: th[(s, tb)][i] for tb in C.TICK_BINS}, key)["decisions"] for i in pans}

    base = fams(names[0])
    rows = []
    for s in names:
        cur = fams(s)
        pan = sorted(set(cur) & set(base))
        for tb in C.TICK_BINS:
            if not pan:
                continue
            row = {"setting": s, "tick_bin": tb, "panels": len(pan)}
            for t in C.PRIMARY_TESTS:
                if tb not in C.FAMILY_TICK_BINS:  # descriptive: no rejection
                    row[f"{t} rej"] = row[f"{t} rej change"] = None
                    continue
                row[f"{t} rej"] = rate([cur[i][(t, tb)]["rejected"] for i in pan])
                row[f"{t} rej change"] = row[f"{t} rej"][0] - rate([base[i][(t, tb)]["rejected"] for i in pan])[0]
            for test, nm in THIN_STATS:
                d = np.array([_stat(th[(s, tb)][i], key, test, nm) - _stat(th[(names[0], tb)][i], key, test, nm)
                              for i in pan])
                d = d[np.isfinite(d)]
                row[f"{test} {nm} change"] = (float(d.mean()), float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1
                                              else np.nan, len(d)) if len(d) else None
            rows.append(row)
    return rows


# ---- report -----------------------------------------------------------------------------------------------


def _cell(r, descriptive: bool = False) -> str:
    return "descriptive" if descriptive else (fmt_rate(r) if r else "-")


def analyse(log=print) -> dict:
    REP.mkdir(parents=True, exist_ok=True)
    key = chosen_design(write=True)
    dj = json.loads((REP / "design.json").read_text())
    res = {"design": dj, "worlds": {}}
    L = ["# Phase 3d-2 tables (generated by analyse3d2)", "",
         f"P-M design used: **{key}** ({'chosen by the rule' if dj['chosen'] else 'FALLBACK: no design admissible'})",
         "", "## P-M design rule", "",
         "| design | WB objective | worst null rate | admissible | " + " | ".join(
             f"{w} {tb}" for w in (*W.SIZE_WORLDS, W.POWER_WORLD) for tb in C.TICK_BINS) + " |",
         "|---|---|---|---|" + "---|" * 9]
    for d in W.PM_DESIGNS:
        k = W.design_key(d)
        cols = [f"{dj['rates'][f'{w}|{k}|{tb}']:.3f}" for w in (*W.SIZE_WORLDS, W.POWER_WORLD) for tb in C.TICK_BINS]
        L.append(f"| {k} | {dj['objective'][k]:.3f} | {dj['worst_null'][k]:.3f} | {dj['admissible'][k]} | "
                 + " | ".join(cols) + " |")
    for w in W.WORLDS:
        if not any((D.OUT / "study").glob(f"{w}_*.json")):
            continue
        t = world_table(w, key)
        res["worlds"][w] = t
        L += ["", f"## {w} ({t['panels']} panels, mean family size {t['family_size']:.2f}; former family "
              f"{t['family_size_11']:.2f})", "",
              "Each after-Holm cell: the family of Q-T1 and Q-T3 (2 Nov 2026) / the former family of all three tick "
              "bins, on the same saved p-values. Q-T2 is descriptive in the family.", "",
              "| test | tick bin | evaluable | before Holm | after Holm | passed | correct outcome | rate correct | "
              "[8,64] / [4,32] | left |", "|---|---|---|---|---|---|---|---|---|---|"]
        for (te, tb), r in t["rows"].items():
            lab = W.CORRECT[w][te] or W.BY_CONSTRUCTION.get((w, te))
            rng_txt = "lags" if te == "P-M" else f"{r['range_8_64']} / {r['range_4_32']}"
            desc = tb not in C.FAMILY_TICK_BINS
            L.append(f"| {te} | {tb} | {r['evaluable'][0]:.2f} | {fmt_rate(r['raw'])} | "
                     f"{_cell(r['holm'], desc)} / {_cell(r['holm_11'])} | "
                     f"{_cell(r['passed'], desc)} / {_cell(r['passed_11'])} | {lab or 'none'} | "
                     f"{_cell(r['correct'], desc)} / {_cell(r['correct_11'])} | {rng_txt} | {r['left']} |")
        if t["fwer"]:
            L.append(f"\nFamily-wise error ({', '.join(t['fwer_tests'])}, after Holm): Q-T1 and Q-T3 family "
                     f"{fmt_rate(t['fwer'])}; former family, any tick bin {fmt_rate(t['fwer_11'])}.")
        L.append(f"\nAny rejection in the family incl. P-R: {fmt_rate(t['any_rejection'])}; former family "
                 f"{fmt_rate(t['any_rejection_11'])}.")
        L += ["", "| tick bin | statistic | mean | SD over panels | panels |", "|---|---|---|---|---|"]
        for tb, sv in t["stats"].items():
            for k, v in sv.items():
                if v:
                    name = f"{k} ({P.PSI_NOTE})" if k == "psi-hat" else k
                    L.append(f"| {tb} | {name} | {v[0]:+.4f} | {v[1]:.4f} | {v[2]} |")
    occ = occupancy()
    res["occupancy"] = occ
    L += ["", "## P-M bin occupancy, real exploration flow (R2-clipped, ties to the lower bin)", "",
          "| cell | bins | stock-days | round-lot share | median lowest-bin share (1/bins) | "
          "stock-days with an empty bin | "
          "min bin >= 20 | >= 50 | >= 100 |", "|---|---|---|---|---|---|---|---|---|"]
    for r in occ:
        L.append(f"| {r['cell']} | {r['n_bins']} | {r['stock_days']} | {r['round_lot_share']:.2f} | "
                 f"{r['lowest_bin_share_median']:.2f} ({r['expected']:.2f}) | {r['empty_bin_share']:.2f} | "
                 f"{r['min_bin_ge_20']:.2f} | {r['min_bin_ge_50']:.2f} | {r['min_bin_ge_100']:.2f} |")
    th = thinning_table(key)
    res["thinning"] = th
    if th:
        L += ["", "## Thinning (WM, paired with share 1.0 on the same market days)", "",
              "| setting | tick bin | panels | " + " | ".join(f"{t} after Holm (change)" for t in C.PRIMARY_TESTS)
              + " | " + " | ".join(f"{a} {b}" for a, b in THIN_STATS) + " |",
              "|---|---|---|" + "---|" * (len(C.PRIMARY_TESTS) + len(THIN_STATS))]
        for r in th:
            L.append(f"| {r['setting']} | {r['tick_bin']} | {r['panels']} | " + " | ".join(
                "descriptive" if r[f"{t} rej"] is None else f"{r[f'{t} rej'][0]:.2f} ({r[f'{t} rej change']:+.2f})"
                for t in C.PRIMARY_TESTS) + " | " + " | ".join(
                "-" if not r[f"{a} {b} change"] else f"{r[f'{a} {b} change'][0]:+.3f} ({r[f'{a} {b} change'][1]:.3f})"
                for a, b in THIN_STATS) + " |")
        L.append(f"\nP-R psi_hat: {P.PSI_NOTE}.")
    (REP / "tables.md").write_text("\n".join(L) + "\n")
    for t in res["worlds"].values():
        t["rows"] = {f"{te}|{tb}": r for (te, tb), r in t["rows"].items()}
    (REP / "results.json").write_text(json.dumps(res, default=str, indent=1))
    log(f"wrote {REP / 'tables.md'}")
    return res


def pilot_tables(log=print) -> str:
    """Rates on the pilot panels (default design, one family per panel across tick bins), indicative only."""
    res = json.loads((REP / "pilot.json").read_text())
    key = W.design_key(D.DEFAULT_PM)
    L = [f"Pilot: {len({r['panel'] for r in res})} panels per world, P-M design {key}. Rates over panels; "
         "MC SE in brackets.", "",
         "| world | test | tick bin | evaluable | before Holm | after Holm | passed | correct | rate correct |",
         "|---|---|---|---|---|---|---|---|---|"]
    out = {}
    for w in [w for w in W.WORLDS if any(r["world"] == w for r in res)]:
        by = {}
        for r in res:
            if r["world"] == w:
                by.setdefault(r["panel"], {})[r["tick_bin"]] = r["result"]
        fams = [family(v, key) for v in by.values() if len(v) == len(C.TICK_BINS)]
        out[w] = fams
        for t in C.PRIMARY_TESTS:
            for tb in C.TICK_BINS:
                decs = [f["decisions"][(t, tb)] for f in fams]
                cr = correct_rate(w, t, decs)
                ev, raw = (rate([d[k] for d in decs]) for k in ("evaluable", "raw"))
                hol, pas = ((rate([d[k] for d in decs]) if tb in C.FAMILY_TICK_BINS else (np.nan, np.nan, 0))
                            for k in ("rejected", "passed"))
                L.append(f"| {w} | {t} | {tb} | {ev[0]:.1f} | {fmt_rate(raw)} | {fmt_rate(hol)} | {fmt_rate(pas)} | "
                         f"{W.CORRECT[w][t] or W.BY_CONSTRUCTION.get((w, t)) or 'none'} | "
                         f"{fmt_rate(cr) if cr else '-'} |")
    txt = "\n".join(L) + "\n"
    (REP / "pilot_rates.md").write_text(txt)
    log(txt)
    return txt
