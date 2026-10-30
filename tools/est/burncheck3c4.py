#!/usr/bin/env python3
"""Phase 3c-4 burn-in check: on the independent-days path, G0 and GT at Q-T1S3, zeta_K(0) and zeta_V(0) on
[4, 32] and [8, 64] at day lengths N = 300 and 3000, for burn-ins B in BURNS, against the 10^6-run truth
(phase3c3, B = 20000). Days are drawn in GROUPS independent groups (each group a "stock" of independent
stationary days); values from the pooled day means, SE from a delete-one-group jackknife.
Picks the smallest B whose values are all within 2 SE of the truth (SE of the value and of the truth
combined) and within TOL of the B = 20000 value."""

from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import panel3c as PC
import phase3c as P3
import phase3c3 as C3
import phase3c4 as C4
import pooled as PO

CELL = "Q-T1S3"
BURNS = (1000, 2000, 5000, 20000)
PLAN = {300: 8000, 3000: 1280}  # days per N (SE <= 0.01 for every value)
GROUPS = 16
T = 2 ** np.arange(7)
RANGES = ((4, 32), (8, 64))
TOL = 0.01


def group(args):
    gen, B, N, g = args
    cl, params, _, _ = P3.setup(CELL, gen, "main")
    D = PLAN[N] // GROUPS
    p = PC.make_panel(
        gen,
        cl,
        np.full((1, D), N),
        np.random.default_rng([C4.SEED, 7, B, N, g, len(gen)]),
        params=dict(params),
        burn=B,
        clip=None,
    )
    pm = PO.panel_moments(
        p.eps, p.q, p.lm, np.ones(len(p.eps), bool), p.starts, p.lengths, T_grid=T, a_grid=np.array([0.0])
    )
    return pm.V[:, 0].mean(0), pm.K[:, 0].mean(0)


def zetas(V, K):
    out = []
    for lo, hi in RANGES:
        sel = (T >= lo) & (T <= hi)
        x = np.log(T[sel]) - np.log(T[sel]).mean()
        out += [1 + np.log(V[sel] / (T[sel] * V[0])) @ x / (x @ x), 1 + np.log(K[sel] / T[sel]) @ x / (x @ x)]
    return np.array(out)  # zV[4,32], zK[4,32], zV[8,64], zK[8,64]


def main() -> int:
    names = ["zeta_V(0) [4,32]", "zeta_K(0) [4,32]", "zeta_V(0) [8,64]", "zeta_K(0) [8,64]"]
    res = {}
    t0 = time.time()
    with ProcessPoolExecutor(3) as ex:
        for gen in ("G0", "GT"):
            for B in BURNS:
                for N in PLAN:
                    out = list(ex.map(group, [(gen, B, N, g) for g in range(GROUPS)]))
                    V, K = np.array([a for a, _ in out]), np.array([b for _, b in out])
                    est = zetas(V.mean(0), K.mean(0))
                    loo = np.array(
                        [zetas(np.delete(V, i, 0).mean(0), np.delete(K, i, 0).mean(0)) for i in range(GROUPS)]
                    )
                    se = np.sqrt((GROUPS - 1) / GROUPS * np.sum((loo - loo.mean(0)) ** 2, axis=0))
                    res[(gen, B, N)] = (est, se)
                    print(
                        f"{gen} B={B:>5} N={N:>4}: "
                        + "  ".join(f"{e:.4f} ({s:.4f})" for e, s in zip(est, se, strict=True))
                        + f"   {time.time() - t0:.0f}s",
                        flush=True,
                    )
    truth = {}
    for gen in ("G0", "GT"):
        tv = C3.truth_values(C3.truth_path("main", CELL, gen))
        truth[gen] = (
            np.array([tv["zeta_V"][0][0, 0], tv["zeta_K"][0][0, 0], tv["zeta_V"][0][0, 1], tv["zeta_K"][0][0, 1]]),
            np.array([tv["zeta_V"][1][0, 0], tv["zeta_K"][1][0, 0], tv["zeta_V"][1][0, 1], tv["zeta_K"][1][0, 1]]),
        )
    rows, ok_by_B = [], {}
    for B in BURNS:
        ok = True
        for gen in ("G0", "GT"):
            for N in PLAN:
                est, se = res[(gen, B, N)]
                ref = res[(gen, 20000, N)][0]
                t, tse = truth[gen]
                z = np.abs(est - t) / np.sqrt(se**2 + tse**2)
                d20 = np.abs(est - ref)
                good = (z <= 2) & (d20 <= TOL)
                ok &= bool(good.all())
                for i, nm in enumerate(names):
                    rows.append(
                        {
                            "gen": gen,
                            "B": B,
                            "N": N,
                            "estimand": nm,
                            "value": est[i],
                            "se": se[i],
                            "truth": t[i],
                            "truth_se": tse[i],
                            "z_vs_truth": z[i],
                            "diff_vs_B20000": d20[i],
                            "ok": int(good[i]),
                        }
                    )
        ok_by_B[B] = ok
    chosen = min([B for B, ok in ok_by_B.items() if ok], default=None)
    P3.write_tsv(
        P3.OUT / "p3c4_burnin_check.tsv",
        rows,
        "Phase 3c-4 burn-in check, Q-T1S3, independent stationary days; "
        f"ok = within 2 SE of the 10^6-run truth and within {TOL} of the B = 20000 value; chosen B = {chosen}",
    )
    (C4.PHASE3C4 / "report").mkdir(parents=True, exist_ok=True)
    (C4.PHASE3C4 / "report" / "burnin.json").write_text(
        json.dumps({"chosen": chosen, "ok_by_B": {str(k): v for k, v in ok_by_B.items()}})
    )
    print("passes by B:", ok_by_B, " chosen B:", chosen, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
