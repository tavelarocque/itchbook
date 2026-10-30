"""Phase 3c tables: TSVs in data/phase3/p3c_*.tsv and markdown fragments in PHASE3C/report/ for
docs/phase3c_report.md (tools/est/phase3c.py analyse | tmin | semi)."""

from __future__ import annotations

import json

import numpy as np

import explore as E
import gen3 as G
import panel3c as PC
import phase3c as P3
import pooled as PO

REP = P3.PHASE3C / "report"
RANGES = ["8_64", "16_256", "64_1024"]
COV_MC = 0.81  # 0.85 minus 2 Monte Carlo SE of one coverage over 200 panels (sqrt(0.9 x 0.1 / 200) = 0.021)


def f(x, d=3) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return "-" if not np.isfinite(v) else f"{v:.{d}f}"


def md(name: str, lines: list[str]) -> None:
    REP.mkdir(parents=True, exist_ok=True)
    (REP / f"{name}.md").write_text("\n".join(lines) + "\n")


def edge_str(e) -> str:
    if e == np.inf:
        return "undef"
    return "-" if not np.isfinite(e) else f"1/{round(1 / e)}"


def write_all(res: dict, n_panels: int) -> None:
    P3.write_tsv(
        P3.OUT / "p3c_truth.tsv", res["truth"], "Phase 3c truth with the R2 clip as a fixed quantile (a <= 1.5)"
    )
    P3.write_tsv(
        P3.OUT / "p3c_bands.tsv", res["bands"], f"Phase 3c bias and 90 % coverage by T/N band, {n_panels} panels"
    )
    P3.write_tsv(P3.OUT / "p3c_rho.tsv", res["rho"], "Phase 3c rho* per cell, generator and estimator")
    # truth table
    L = [
        "| cell | gen | days | domain SE | zeta_K(a) [8,64], a = 0..1.5 | zeta_V(a) [8,64] | P-F slope [8,64] | "
        "slope zV / slope zK [8,64] (a <= 0.5) | ratio [8,64] | ratio [16,256] | ratio [64,1024] | "
        "slope zV - 2 slope zK [8,64] | R1/2/R0 [8,64] |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in res["truth"]:
        L.append(
            f"| {r['cell']} | {r['gen']} | {r['n_days']} | {f(r['se_max'], 4)} | {r['zK_8_64']} | {r['zV_8_64']} | "
            f"{f(r['pf_8_64'])} ({f(r['pf_8_64_se'])}) | {f(r['slopeV_8_64'])} / {f(r['slopeK_8_64'])} | "
            f"{f(r['ratio_8_64'])} ({f(r['ratio_8_64_se'])}) | {f(r['ratio_16_256'])} ({f(r['ratio_16_256_se'])}) | "
            f"{f(r['ratio_64_1024'])} ({f(r['ratio_64_1024_se'])}) | "
            f"{f(r['vk_diff_8_64'])} ({f(r['vk_diff_8_64_se'])}) | "
            f"{f(r['Rratio_8_64'])} |"
        )
    md("truth", L)
    # rho* table
    rho = {}
    for r in res["rho"]:
        rho.setdefault((r["tag"], r["cell"]), {})[(r["gen"], r["estimator"])] = r["rho_star"]
    L = [
        "| sign model | cell | " + " | ".join(f"{g} {e}" for g in P3.GENS for e in P3.ESTIMATORS) + " | rho* (min) |",
        "|---|---|" + "---|" * (len(P3.GENS) * len(P3.ESTIMATORS) + 1),
    ]
    cell_rho = {}
    for (tag, c), d in rho.items():
        vals = [d[(g, e)] for g in P3.GENS for e in P3.ESTIMATORS]
        fin_vals = [x for x in vals if x != np.inf]
        mn = np.nan if np.any(np.isnan(fin_vals)) else min(fin_vals)
        cell_rho[(tag, c)] = mn
        L.append(f"| {tag} | {c} | " + " | ".join(edge_str(v) for v in vals) + f" | {edge_str(mn)} |")
    md("rho", L)
    # sensitivity (Choice, not a decision): rho* under looser pass rules, from the same band rows
    crit = {
        "draft": lambda r: r["pass"],
        "bias only": lambda r: r["max_abs_bias"] < P3.BIAS_MAX,
        "bias + coverage >= 0.81": lambda r: r["max_abs_bias"] < P3.BIAS_MAX and r["min_cov"] >= COV_MC,
    }
    sens = {}
    L = [
        "| rule | sign model | cell | " + " | ".join(P3.GENS) + " | rho* (min over generators and estimators) |",
        "|---|---|---|---|---|---|---|",
    ]
    for name, rule in crit.items():
        for tag, c in rho:
            per_gen = []
            for g in P3.GENS:
                vals = []
                for e in P3.ESTIMATORS:
                    if (g, e) in P3.UNDEFINED:
                        continue
                    rows = [
                        dict(r, **{"pass": bool(r["evaluable"] and rule(r))})
                        for r in res["bands"]
                        if (r["tag"], r["cell"], r["gen"], r["estimator"]) == (tag, c, g, e)
                    ]
                    vals.append(P3.rho_star(rows))
                per_gen.append(np.nan if np.any(np.isnan(vals)) else min(vals))
            mn = np.nan if np.any(np.isnan(per_gen)) else min(per_gen)
            sens[(name, tag, c)] = mn
            L.append(f"| {name} | {tag} | {c} | " + " | ".join(edge_str(x) for x in per_gen) + f" | {edge_str(mn)} |")
    md("rho_sens", L)
    # Choice for steps 3 and 4: G0 and GT only (GF's K-based estimators carry the K_a(1) exclusion bias) and
    # every estimator except the ratio (a ratio of two noisy slopes, biased beyond 0.02 in every band)
    L = [
        "| cell | rho*, draft bounds, G0 + GT, no ratio: main | tail | the same with coverage >= 0.81: main | tail |",
        "|---|---|---|---|---|",
    ]
    for c in P3.q_cells():
        cells = []
        for rname in ("draft", "bias + coverage >= 0.81"):
            for tag in ("main", "tail"):
                if (tag, c) not in rho:
                    cells.append("")
                    continue
                vals = []
                for g in ("G0", "GT"):
                    for e in P3.ESTIMATORS:
                        if e == "ratio":
                            continue
                        rows = [
                            dict(r, **{"pass": bool(r["evaluable"] and crit[rname](r))})
                            for r in res["bands"]
                            if (r["tag"], r["cell"], r["gen"], r["estimator"]) == (tag, c, g, e)
                        ]
                        vals.append(P3.rho_star(rows))
                mn = np.nan if np.any(np.isnan(vals)) else min(vals)
                sens[(rname, "G0GT", tag, c)] = mn
                cells.append(edge_str(mn))
        L.append(f"| {c} | " + " | ".join(cells) + " |")
        sens[("G0GT", "main", c)] = sens[("draft", "G0GT", "main", c)]
    md("rho_g0gt", L)
    # band detail: worst |bias| and min coverage per band, pooled over estimators, per cell and generator
    L = [
        "| sign model | cell | gen | estimator | " + " | ".join(f"T/N <= {edge_str(e)}" for e in P3.EDGES[1:]) + " |",
        "|---|---|---|---|" + "---|" * 5,
    ]
    for tag, c, g, e in sorted({(r["tag"], r["cell"], r["gen"], r["estimator"]) for r in res["bands"]}):
        cells = []
        for edge in P3.EDGES[1:]:
            r = next(
                x
                for x in res["bands"]
                if (x["tag"], x["cell"], x["gen"], x["estimator"], x["edge"]) == (tag, c, g, e, edge)
            )
            cells.append(
                "n/e"
                if not r["evaluable"]
                else f"{f(r['max_abs_bias'])} / {f(r['min_cov'], 2)}{'' if r['pass'] else ' x'}"
            )
        L.append(f"| {tag} | {c} | {g} | {e} | " + " | ".join(cells) + " |")
    md("bands", L)
    # step 3: coverage vs ESS at each cell's rho* (the ceiling mask)
    ess_rows, pts = [], {"V": [], "K": []}
    for (tag, c, g), (st, tq) in res["study"].items():
        if tag != "main":
            continue
        rs = sens[("G0GT", tag, c)]
        rs_used = rs if np.isfinite(rs) else P3.EDGES[1]
        m = P3.CEIL[rs_used]
        for key, kind in (("zeta_V", "V"), ("zeta_K", "K")):
            _, cov, fin = P3.bias_coverage(st, tq, key, m)
            ess = np.nanmean(st[f"ess{kind}"], axis=0)
            for i, a in enumerate(PO.A_GRID):
                for r in range(3):
                    if fin[i, r] >= P3.MIN_EVAL:
                        pts[kind].append((float(ess[i]), float(cov[i, r])))
                        ess_rows.append(
                            {
                                "cell": c,
                                "gen": g,
                                "weights": kind,
                                "a": a,
                                "range": RANGES[r],
                                "rho_used": rs_used,
                                "median_ess": ess[i],
                                "coverage": cov[i, r],
                            }
                        )
    P3.write_tsv(
        P3.OUT / "p3c_ess_coverage.tsv",
        ess_rows,
        "Phase 3c coverage of zeta(a) at each cell's rho* vs the panels' median stock-day Kish ESS (R2 clip)",
    )
    thr = {k: P3.ess_threshold(v) for k, v in pts.items()}
    thr81 = {k: P3.ess_threshold(v, COV_MC) for k, v in pts.items()}
    L = [
        f"Coverage at each cell's rho* (draft bounds, G0 and GT, no ratio; 1/80 where none). ESS "
        f"threshold (smallest ESS with coverage >= {P3.COV_MIN} at and above it): V weights {f(thr['V'], 1)}, "
        f"K weights {f(thr['K'], 1)}; with coverage >= {COV_MC}: V {f(thr81['V'], 1)}, K {f(thr81['K'], 1)}. "
        f"Points (generator x a x range): V {len(pts['V'])}, K {len(pts['K'])}. a_max below uses the 0.81 "
        "thresholds where the 0.85 ones do not exist.",
        "",
        "| cell | rho* used | a_max (V weights, real median ESS) | a_max (K weights) | a_max (both) |",
        "|---|---|---|---|---|",
    ]
    amax_rows = []
    for c in P3.q_cells():
        re = P3.real_ess(c)
        am = {}
        for k in ("V", "K"):
            t = thr[k] if np.isfinite(thr[k]) else thr81[k]
            ok = [a for a, e in sorted(re[k].items()) if a <= 1.5 and np.isfinite(t) and e >= t]
            am[k] = max(ok) if ok else np.nan
        both = np.nanmin([am["V"], am["K"]]) if np.all(np.isfinite([am["V"], am["K"]])) else np.nan
        amax_rows.append({"cell": c, "a_max_V": am["V"], "a_max_K": am["K"], "a_max": both})
        L.append(f"| {c} | {edge_str(sens[('G0GT', 'main', c)])} | {f(am['V'], 2)} | {f(am['K'], 2)} | {f(both, 2)} |")
    # binned coverage vs ESS for the report
    L += [
        "",
        "| median ESS bin | V points | V min coverage | V mean coverage | K points | K min coverage | "
        "K mean coverage |",
        "|---|---|---|---|---|---|---|",
    ]
    bins = [0, 5, 10, 20, 30, 50, 100, 300, 1e9]
    for lo, hi in zip(bins[:-1], bins[1:], strict=True):
        cells = []
        for k in ("V", "K"):
            c_ = [cv for e, cv in pts[k] if lo <= e < hi]
            cells += [str(len(c_)), f(min(c_), 2) if c_ else "-", f(np.mean(c_), 2) if c_ else "-"]
        L.append(f"| {lo:g} to {hi:g} | " + " | ".join(cells) + " |")
    md("ess", L)
    P3.write_tsv(
        P3.OUT / "p3c_amax.tsv",
        amax_rows,
        f"Phase 3c a_max per cell from the ESS thresholds V {thr['V']}, K {thr['K']}",
    )
    # step 4: N_floor table
    days, u = E.read_universe()
    head = None
    nmin = {}
    for line in open(E.UNIVERSE):
        if line.startswith("#"):
            continue
        row = line.rstrip("\n").split("\t")
        if head is None:
            head = {k: i for i, k in enumerate(row)}
            continue
        if row[head["split_v2"]] == "confirm" and row[head["cell"]].startswith("Q-"):
            nmin.setdefault(row[head["cell"]], []).append(int(row[head["n_min"]]))
    nf_rows = []
    L = [
        "| cell | ceiling rho | confirm stocks | "
        + " | ".join(f"N_floor {n}: stocks / T_max / ranges reached" for n in (50, 100, 200, 400, 800))
        + " |",
        "|---|---|---|" + "---|" * 5,
    ]
    for c, rs_used, label in [
        (cc, r, lab)
        for cc in P3.q_cells()
        for r, lab in ((sens[("G0GT", "main", cc)], "rho*"), (1 / 20, "1/20 (reference)"), (1 / 10, "1/10 (reference)"))
    ]:
        if not np.isfinite(rs_used):
            L.append(f"| {c} | rho*: none | {len(nmin[c])} | - | - | - | - | - |")
            continue
        label = f"rho* = {edge_str(rs_used)}" if label == "rho*" else label
        nm = np.array(nmin[c])
        cells = []
        for nfl in (50, 100, 200, 400, 800):
            keep = nm >= nfl
            tmax = rs_used * nfl
            reach = []
            for _, hi in PO.RANGES:
                n_reach = int(np.sum(nm[keep] >= hi / rs_used))
                reach.append(f"{n_reach}")
            nf_rows.append(
                {
                    "cell": c,
                    "rho": label,
                    "rho_used": rs_used,
                    "N_floor": nfl,
                    "stocks": int(keep.sum()),
                    "T_max": tmax,
                    **{f"reach_{lo}_{hi}": int(np.sum(nm[keep] >= hi / rs_used)) for lo, hi in PO.RANGES},
                }
            )
            cells.append(f"{int(keep.sum())} / {tmax:.0f} / {'+'.join(reach)}")
        L.append(f"| {c} | {label} | {len(nm)} | " + " | ".join(cells) + " |")
    L += ["", "ranges reached = remaining stocks whose N_min >= T2 / rho* for [8,64] + [16,256] + [64,1024]."]
    md("nfloor", L)
    P3.write_tsv(
        P3.OUT / "p3c_nfloor.tsv", nf_rows, "Phase 3c N_floor trade-off (confirmation N_min from universe.tsv)"
    )
    # sigma_l bias
    sg_rows = []
    L = [
        "| cell | " + " | ".join(f"N = {n}: bias (sd, stocks)" for n in (50, 100, 200, 400, 800)) + " |",
        "|---|" + "---|" * 5,
    ]
    for c in P3.q_cells():
        sb = P3.sigma_bias(c)
        L.append(f"| {c} | " + " | ".join(f"{f(b)} ({f(s)}, {n})" for b, s, n in sb.values()) + " |")
        sg_rows += [{"cell": c, "N": N, "bias": b, "sd": s, "stocks": n} for N, (b, s, n) in sb.items()]
    md("sigma", L)
    P3.write_tsv(
        P3.OUT / "p3c_sigma_bias.tsv",
        sg_rows,
        "Phase 3c sigma_l-hat minus full-sample sigma_l, N contiguous runs per day",
    )
    # timing
    L = ["| sign model | cell | gen | panels | seconds (one worker) |", "|---|---|---|---|---|"]
    for t in res["time"]:
        L.append(f"| {t['tag']} | {t['cell']} | {t['gen']} | {t['panels']} | {t['seconds']:.0f} |")
    md("time", L)
    (REP / "rho.json").write_text(json.dumps({f"{k[0]} {k[1]}": v for k, v in cell_rho.items()}))


def write_tmin() -> None:
    rows = []
    L = [
        "| cell | P(Delta_T != 0) at T = 1, 2, 4, 8, 16, 32, 64 | T_min | l_max (N_floor 100) | beta | "
        "beta-hat bias (sd) |",
        "|---|---|---|---|---|---|",
    ]
    for c in P3.q_cells():
        p = P3.tmin_cell(c)
        T = 2 ** np.arange(len(p))
        tmin = int(T[np.argmax(p >= 0.5)]) if np.any(p >= 0.5) else np.nan
        lmax = P3.lmax_cell(c)
        bb = P3.beta_hat_bias(c, lmax)
        rows.append(
            {"cell": c, "T_min": tmin, "l_max": lmax, **{f"p_T{t}": v for t, v in zip(T, p, strict=True)}, **bb}
        )
        L.append(
            f"| {c} | {' '.join(f(x, 2) for x in p)} | {tmin} | {lmax} | {f(bb['beta'])} | "
            f"{f(bb['bias'])} ({f(bb['sd'])}) |"
        )
    md("tmin", L)
    P3.write_tsv(
        P3.OUT / "p3c_tmin_lmax.tsv", rows, "Phase 3c T_min from exploration mids alone; l_max and the G0 beta-hat bias"
    )


def load_cell_flow(cell: str):
    """Real exploration signs and sizes of a cell through the guard, stock-major, and day lengths (S, 8)."""
    days, _ = E.read_universe()
    per = {d: E.load_explore_runs(d, cell=cell) for d in days}
    syms = sorted(set.intersection(*(set(r) for r in per.values())))
    eps, q, Nsd = [], [], []
    for s in syms:
        row = []
        for d in days:
            r = per[d][s]
            eps.append(r["sign"].astype(float))
            q.append(r["q"].astype(float))
            row.append(len(r))
        Nsd.append(row)
    return np.concatenate(eps), np.concatenate(q), np.array(Nsd)


def semi_smoke(cell: str = "Q-T3S1") -> None:
    eps, q, Nsd = load_cell_flow(cell)
    cl = G.load_cell(cell)
    lines = []
    for gen in P3.GENS:
        p = PC.semi_synthetic(eps.copy(), q.copy(), Nsd, gen, cl, np.random.default_rng(1))
        same = bool(np.array_equal(p.eps, eps) and np.array_equal(p.q, q))
        lines.append(
            f"- {gen}: flow identical {same}; runs {len(eps)}; share of runs clipped {np.mean(p.q_mom < p.q):.4f}"
        )
    md(
        "semi",
        [f"Semi-synthetic smoke check on the real exploration flow of {cell} ({Nsd.shape[0]} stocks x 8 days):"]
        + lines,
    )
