"""Tests for Phase 3d-2 (worlds3d2, phase3d2, analyse3d2): world labels and specs, the P-M design rule on the
right worlds, thinning shares applied once and paired settings, the burn-in drawn from the stock's other days,
the WB kernel, the response against a direct sum, and the driver's tick-bin run against run_primary. Every
input is synthetic; no data file is read (the WM test replaces panel3c.day_lengths).

    uv run --group dev --group analysis pytest tools/est/test_3d2.py
"""
from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pytest

import analyse3d2 as A
import diag3d2b as DB
import gen3 as G
import mbmodel as M
import panel3c as PC
import phase3d2 as D
import prereg3d as C
import primary3d as P
import worlds3d2 as W


@pytest.fixture(autouse=True)
def restore_grid():
    """Every test leaves prereg3d's grid and primary3d's KW as it found them (the planted grid-patch bug
    leaves them patched, which would leak into later tests and other files)."""
    saved = (C.T_GRID, C.RANGES, C.PRIMARY_RANGE, C.FALLBACK_RANGE, P.KW, P.RANGE_NAMES)
    yield
    C.T_GRID, C.RANGES, C.PRIMARY_RANGE, C.FALLBACK_RANGE, P.KW, P.RANGE_NAMES = saved


CELL = G.Cell("TEST", np.arange(100, 200), 3.0, 1.0, (0.3, 0.7), {"w1": 0.3, "H": 0.8, "w2": 0.0, "phi": 0.5})


# ---- world labels and specs ---------------------------------------------------------------------------------


def check_world_labels():
    """The correct outcome per world and test, as instructed; P-R in W-null and everything in WM unlabelled."""
    assert W.CORRECT["W-null"] == {"P-V": "no reject", "P-R": None, "P-M": "no reject", "P-F": "no reject"}
    assert W.CORRECT["W0"]["P-M"] == "no reject" and W.CORRECT["W0"]["P-R"] == "pass"
    assert W.CORRECT["W0-psi3"]["P-R"] == "fail"
    assert W.CORRECT["WB"]["P-M"] == "reject"
    assert all(v is None for v in W.CORRECT["WM"].values())
    assert W.SIZE_WORLDS == ("W0", "W-null") and W.POWER_WORLD == "WB"
    assert (W.PV_SIZE_WORLD, W.PV_POWER_WORLD) == ("W-null", "WM")
    d = [{"rejected": True, "passed": True}, {"rejected": True, "passed": False}, {"rejected": False, "passed": False},
         {"rejected": False, "passed": False}]
    assert A.correct_rate("W-null", "P-F", d)[0] == 0.5  # no reject: 2 of 4
    assert A.correct_rate("WB", "P-M", d)[0] == 0.25  # reject needs both conditions: 1 of 4
    assert A.correct_rate("W0", "P-R", d)[0] == 0.25  # pass
    assert A.correct_rate("W0-psi3", "P-R", d)[0] == 0.75  # fail = not passed
    assert A.correct_rate("WM", "P-V", d) is None


def check_world_specs():
    assert W.world_spec("W-null") == {"flow": "real_s1", "signs": "coin_informed", "gen": "GF", "psi": 0.5}
    assert W.world_spec("W0") == {"flow": "real", "signs": "real", "gen": "G0", "psi": 0.5}
    assert W.world_spec("W0-psi3") == {"flow": "real", "signs": "real", "gen": "G0", "psi": 0.3}
    assert W.world_spec("W-null-real")["signs"] == "real" and W.world_spec("W0-shuf")["signs"] == "day_shuffled"
    assert not set(W.DIAG_WORLDS) & set(W.WORLDS)
    assert W.world_spec("WB")["gen"] == "WB" and W.world_spec("WM")["flow"] == "mb"
    assert W.WB_BETA == (0.6, -0.3) and W.BURN == 20_000
    with pytest.raises(ValueError):
        W.world_spec("W1")


# ---- the P-M design rule ------------------------------------------------------------------------------------


def toy_rates(power: dict, null: dict, default_null=0.03, default_power=0.1) -> dict:
    """power / null: design key -> rate in every tick bin (null applies to both size worlds unless keyed by
    (world, key))."""
    r = {}
    for d in W.PM_DESIGNS:
        k = W.design_key(d)
        for tb in C.TICK_BINS:
            r[("WB", k, tb)] = power.get(k, default_power)
            for w in W.SIZE_WORLDS:
                r[(w, k, tb)] = null.get((w, k), null.get(k, default_null))
    return r


N200 = {"W0": 200, "W-null": 200, "WB": 100}


def check_design_rule():
    """Highest WB rate among designs within 0.05 + 2 MC SE in BOTH size worlds, in every tick bin. The worlds
    disagree on purpose: W0's highest rate is another design, W-null alone excludes one design and W0 alone
    (in one tick bin) another."""
    k_best, k_w0, k_null, k_one = "b5_l32_L256_m50", "b3_l16_L128_m50", "b4_l16_L128_m50", "b4_l32_L256_m50"
    rates = toy_rates({k_best: 0.6, k_null: 0.9, k_one: 0.8}, {k_best: 0.0, k_w0: 0.045},
                      default_null=0.01)
    for tb in C.TICK_BINS:
        rates[("W-null", k_null, tb)] = 0.2
    rates[("W0", k_one, "Q-T2")] = 0.12
    rates[("W0", k_best, "Q-T3")] = 0.08  # inside the band: 0.05 + 2 sqrt(0.05 0.95 / 200) = 0.081
    out = W.choose_design(rates, N200)
    assert not out["admissible"][k_null] and not out["admissible"][k_one] and out["admissible"][k_best]
    assert out["chosen"] == k_best and out["fallback"] is None
    assert out["objective"][k_best] == pytest.approx(0.6)


def check_design_rule_ties_and_fallback():
    keys = [W.design_key(d) for d in W.PM_DESIGNS]
    out = W.choose_design(toy_rates({k: 0.5 for k in keys}, {}), N200)
    assert out["chosen"] == "b3_l16_L128_m50"  # every rate equal: fewer bins, then l_max
    tie = {"b5_l16_L128_m50": 0.7, "b4_l32_L256_m50": 0.7}
    assert W.choose_design(toy_rates(tie, {}), N200)["chosen"] == "b4_l32_L256_m50"
    worst = toy_rates({}, {k: 0.3 for k in keys})
    for tb in C.TICK_BINS:
        for w in W.SIZE_WORLDS:
            worst[(w, "b5_l32_L256_m50", tb)] = 0.12
    out = W.choose_design(worst, N200)
    assert out["chosen"] is None and out["fallback"] == "b5_l32_L256_m50"


def check_design_grid():
    """The 6 designs of 1 Nov 2026: bins {3, 4, 5} x l_max {16, 32}, L = 8 l_max, min runs per bin 50; the
    rule and stage A run over exactly these."""
    got = {(d.n_bins, d.l_max, d.L, d.min_runs_per_bin) for d in W.PM_DESIGNS}
    assert len(W.PM_DESIGNS) == 6 and got == {(nb, lm, 8 * lm, 50) for nb in (3, 4, 5) for lm in (16, 32)}
    assert set(W.PM_DESIGNS) < set(W.PM_DESIGNS_18) and len(W.PM_DESIGNS_18) == 18
    keys = [W.design_key(d) for d in W.PM_DESIGNS]
    assert set(W.choose_design(toy_rates({}, {}), N200)["admissible"]) == set(keys)
    js = D.study_jobs({}, {"W0": 4}, keys, ("W0",))
    assert all(j[4] == tuple(keys) for j in js) and len(js) == 6


# ---- thinning ----------------------------------------------------------------------------------------------


def check_router_share_once():
    """The observed share of the router rule is the target share: h is applied once (inside
    router_thinning), so p = share / (1 - h)."""
    assert W.router_p(0.46, 0.17) == pytest.approx(0.46 / 0.83)
    with pytest.raises(ValueError):
        W.router_p(1.0, 0.17)
    day = W.mb_market_day(40_000, np.random.default_rng(1), burn=1000)
    for share in (0.46, 0.17):
        seen = []
        for s in range(6):
            o = W.thin_day(day, ("router", share, 0.17, 16.0), np.random.default_rng(s))
            seen.append(len(o.t) / len(day.t))
        assert abs(np.mean(seen) - share) < 0.01, (share, np.mean(seen))  # seed SD 0.003 (measured)
        o = W.thin_day(day, ("random", share), np.random.default_rng(9))
        assert abs(len(o.t) / len(day.t) - share) < 0.02
    assert W.thin_day(day, ("full",), np.random.default_rng(0)) is day


def check_thin_settings_paired():
    """A job that keeps a subset of the settings sees the same markets: its flows equal the corresponding
    flows of the job that keeps all of them."""
    real = PC.day_lengths
    PC.day_lengths = lambda cell, n, n_days, rng: np.full((n, n_days), 150)
    try:
        settings = (("full",), ("random", 0.26), ("router", 0.17, 0.17, 1.0))

        def run(keep):
            return W.wm_cell_flows("Q-T3S1", np.random.default_rng(3), settings, n_stocks=1,
                                   thin_rngs=[np.random.default_rng([5, k]) for k in range(3)], keep=keep)

        every = run(None)
        for k in range(3):
            one = run([k])[0]
            for f in ("eps", "q", "lm", "ts"):
                assert np.array_equal(one[f], every[k][f]), (k, f)
        assert not np.array_equal(every[0]["lm"], every[1]["lm"])
        assert every[0]["lengths"].sum() == len(every[2]["eps"])
    finally:
        PC.day_lengths = real


def check_setting_grid():
    names = [D.setting_name(s) for s in D.THIN_SETTINGS]
    assert names[0] == "full 1.0" and len(names) == 1 + 3 + 9 and len(set(names)) == len(names)
    assert sorted(k for g in D.THIN_GROUPS for k in g) == list(range(len(names)))
    assert {s[1] for s in D.THIN_SETTINGS[1:]} == {0.46, 0.26, 0.17}
    assert all(s[2] == 0.17 for s in D.THIN_SETTINGS if s[0] == "router")


# ---- responses ---------------------------------------------------------------------------------------------


def check_burn_from_other_days():
    Nd = np.array([5, 7, 3, 9])
    rng = np.random.default_rng(0)
    for d in range(4):
        b = W.burn_index(100, Nd, d, 50, rng)
        s0 = 100 + np.concatenate([[0], np.cumsum(Nd)[:-1]])
        own = set(range(s0[d], s0[d] + Nd[d]))
        assert len(b) == 50 and not own & set(b.tolist())
        assert b.min() >= 100 and b.max() < 100 + Nd.sum()
        # contiguous cyclic stretch of the other days
        other = np.concatenate([np.arange(s0[j], s0[j] + Nd[j]) for j in range(4) if j != d])
        pos = np.searchsorted(other, b)
        assert np.all(np.diff(pos) % len(other) == 1)


def direct(x, beta):
    n = len(x)
    out = np.zeros(n)
    for t in range(1, n):
        k = np.arange(t)
        out[t] = np.sum(x[k] * (1.0 + t - k) ** -np.broadcast_to(beta, (n,))[k])
    return out


def check_wb_beta_and_impact():
    q = np.exp(np.random.default_rng(1).normal(5, 1, 4000))
    m1, m5 = W.wb_anchors(q)
    assert m1 < np.log(np.median(q)) < m5
    b = W.wb_beta(np.exp(np.array([m1 - 1, m1, (m1 + m5) / 2, m5, m5 + 1])), m1, m5)
    assert np.allclose(b, [0.6, 0.6, 0.15, -0.3, -0.3])
    lq = np.sort(np.log(q))
    bb = W.wb_beta(np.exp(lq), m1, m5)
    assert np.all(np.diff(bb) <= 1e-15)  # falls with size
    x = np.random.default_rng(2).standard_normal(300)
    beta = W.wb_beta(q[:300], m1, m5)
    got = W.wb_impact(x, beta, np.array([0]), np.array([300]))
    assert np.allclose(got, direct(x, beta), rtol=1e-7, atol=1e-8)


def check_semi_response_direct():
    """G0 and WB mids on stationary days equal a direct sum over [burn-in from the other days, the day],
    burn-in dropped, plus the random walk; the same generator draws give the same burn-ins."""
    S_, D_, N, burn = 2, 3, 40, 60
    rng = np.random.default_rng(4)
    Nsd = np.full((S_, D_), N)
    eps = np.where(rng.random(S_ * D_ * N) < 0.5, -1.0, 1.0)
    q = np.exp(rng.normal(3, 1, S_ * D_ * N))
    for world in ("W0", "WB"):
        r1 = np.random.default_rng(7)
        r2 = copy.deepcopy(r1)
        lm = W.semi_response(world, CELL, eps, q, Nsd, r1, burn=burn)
        theta = 1 / G.mean_q_psi(CELL, 0.5, "lognormal")
        exp, ser = [], []
        for i in range(S_ * D_):
            s, d = divmod(i, D_)
            b = W.burn_index(s * D_ * N, Nsd[s], d, burn, r2)
            ii = np.concatenate([b, np.arange(i * N, (i + 1) * N)])
            beta = (G.g0_beta(CELL) if world == "W0" else
                    W.wb_beta(q[ii], *W.wb_anchors(q[s * D_ * N : (s + 1) * D_ * N])))
            exp.append(direct(eps[ii] * q[ii] ** 0.5 * theta, beta))
            ser.append(len(ii))
        ser = np.array(ser)
        sst = np.concatenate([[0], np.cumsum(ser)[:-1]])
        rw = PC.random_walk(sst, ser, 1.0, r2)
        want = np.concatenate([e[burn:] + rw[a + burn : a + n] for e, a, n in zip(exp, sst, ser, strict=True)])
        assert np.allclose(lm, want, rtol=1e-7, atol=1e-7), world


def check_semi_flow_unchanged():
    """The real flow passes through: resampled stocks keep their days whole; W-null permutes sizes within
    each stock-day only."""
    rng = np.random.default_rng(0)
    Nsd = np.array([[3, 4], [2, 5], [6, 1]])
    n = Nsd.sum()
    fc = {"eps": np.where(rng.random(n) < 0.5, -1, 1).astype(np.int8), "q": np.arange(1.0, n + 1),
          "ts": np.arange(n, dtype=np.int64), "Nsd": Nsd}
    f = W.resample_stocks(fc, 5, np.random.default_rng(1))
    per = Nsd.sum(1)
    s0 = np.concatenate([[0], np.cumsum(per)[:-1]])
    a = 0
    for s in f["stocks"]:
        assert np.array_equal(f["q"][a : a + per[s]], fc["q"][s0[s] : s0[s] + per[s]])
        assert np.array_equal(f["eps"][a : a + per[s]], fc["eps"][s0[s] : s0[s] + per[s]])
        a += per[s]
    assert np.array_equal(f["Nsd"], Nsd[f["stocks"]])
    lengths = f["Nsd"].ravel()
    sh = W.shuffle_sizes(f["q"], lengths, np.random.default_rng(2))
    st = np.concatenate([[0], np.cumsum(lengths)[:-1]])
    for b, m in zip(st, lengths, strict=True):
        assert np.array_equal(np.sort(sh[b : b + m]), np.sort(f["q"][b : b + m]))
    assert not np.array_equal(sh, f["q"])


def check_wm_burn():
    """mb_market_day drops the first `burn` trades but keeps their impact in the mid."""
    p = M.PAPER2
    d = W.mb_market_day(3000, np.random.default_rng(5), burn=2000)
    n_meta = int(1.25 * 5000 / G._trades_per_meta(p)) + 50
    full = M.simulate_day(p, n_meta, np.random.default_rng(5))
    assert len(d.t) >= 3000 and np.array_equal(d.mid, full.mid[2000:]) and np.array_equal(d.q, full.q[2000:])
    ts = W.clock_ts(d.t)
    assert ts[0] == W.DAY_NS[0] and ts[-1] < W.DAY_NS[1] and np.all(np.diff(ts) >= 0)


# ---- the driver's tick-bin run against run_primary ----------------------------------------------------------


def g0_cell_flow(S, Dn, N, seed) -> P.Flow:
    p = PC.make_panel("G0", CELL, np.full((S, Dn), N), np.random.default_rng(seed), params={"beta": 0.3}, clip=None)
    ts = np.concatenate([34200 * 10**9 + np.arange(n) * (23400 * 10**9 // n) for n in p.lengths]).astype(np.int64)
    return P.Flow(p.eps, p.q, p.lm, np.ones(len(p.eps), bool), ts, p.starts, p.lengths, Dn)


def check_tick_bin_matches_run_primary():
    """phase3d2.run_tick_bin (P-M at several designs on run_unit's draws, a result reused across min-runs
    values that keep the same stock-days) + analyse3d2.family give exactly run_primary's p-values and Holm
    decisions at each design."""
    flows = {"Q-T3S1": g0_cell_flow(6, 4, 900, 1), "Q-T3S2": g0_cell_flow(5, 4, 900, 2)}
    designs = [P.PMParams(3, 16, 128, 20), P.PMParams(3, 16, 128, 50), P.PMParams(4, 16, 128, 50),
               P.PMParams(4, 16, 128, 400)]  # 400 drops stock-days: no reuse; 20 and 50 keep the same ones
    rng = np.random.default_rng(11)
    rec = D.run_tick_bin(flows, "Q-T3", rng, designs, n_boot=30)
    for d in designs:
        ref = P.run_primary(flows, d, np.random.default_rng(11), n_boot=30)
        fam = A.family({"Q-T3": rec}, W.design_key(d))
        for t in C.PRIMARY_TESTS:
            got = rec["pm"][W.design_key(d)] if t == "P-M" else rec["tests"][t]
            want = ref["by_tick_bin"]["Q-T3"]["tests"][t]
            assert got["p"] == pytest.approx(want.p, abs=0), (t, d)
            if t == "P-M" and d.min_runs_per_bin < 400:  # p_slope is intermediate under G0: it sees the draws
                assert got["p_slope"] == want.statistic["p_slope"] and 0.05 < got["p_slope"] < 0.95
            assert fam["decisions"][(t, "Q-T3")]["rejected"] == ref["holm"][(t, "Q-T3")]["rejected"]
        assert fam["family_size"] == len(ref["family"])
    assert rec["pm"]["b4_l16_L128_m400"]["stock_days"] < rec["pm"]["b4_l16_L128_m50"]["stock_days"]


def check_m3_short_days():
    """primary3d.m3_V on days shorter than the largest horizon (fix of 3d-2): finite where the day has a
    window, NaN beyond; a day of 5 runs no longer crashes."""
    lengths = np.array([5, 40, 200])
    rng = np.random.default_rng(0)
    n = lengths.sum()
    f = P.Flow(np.where(rng.random(n) < 0.5, -1.0, 1.0), np.exp(rng.normal(3, 1, n)), np.zeros(n), np.ones(n, bool),
               np.arange(n, dtype=np.int64), np.concatenate([[0], np.cumsum(lengths)[:-1]]), lengths, 3)
    V = P.m3_V(f)
    for i, N in enumerate(lengths):
        assert np.all(np.isfinite(V[i][:, C.T_GRID <= N])) and np.all(np.isnan(V[i][:, C.T_GRID > N]))



# ---- Phase 3d-2b: tie-aware bins, the clean P-F null, the long-day unit ---------------------------------------


def tied_sizes(seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.permutation(np.concatenate([np.full(400, 100.0), np.full(150, 200.0), rng.integers(1, 99, 300),
                                           rng.integers(201, 2000, 150)]).astype(float))


def check_tie_bins():
    """Every size in one bin, every bin at least min_runs and none empty, the known answer on a tied sample,
    one bin for one size, equal bins without ties; pm_moments uses only stock-days at the target count."""
    q = tied_sizes()
    for nb, mr in [(3, 20), (5, 50), (10, 100), (5, 400)]:
        b, upper = P.tie_bins(q, nb, mr)
        k = len(upper) + 1
        assert k <= nb and np.array_equal(np.unique(b), np.arange(k))
        cnt = np.bincount(b, minlength=k)
        assert cnt.min() >= min(mr, len(q)) and cnt.min() > 0
        for v in np.unique(q):
            assert len(np.unique(b[q == v])) == 1, (nb, mr, v)
    b, upper = P.tie_bins(q, 5, 50)
    assert np.bincount(b).tolist() == [199, 101, 400, 150, 150]
    assert len(P.tie_bins(np.full(30, 7.0), 3, 1)[1]) == 0
    assert np.bincount(P.tie_bins(np.arange(9.0), 3, 1)[0]).tolist() == [3, 3, 3]
    n = 3 * len(q)
    lm = np.cumsum(np.random.default_rng(1).normal(size=n))  # any mid with moves (rms r > 0)
    f = P.Flow(np.ones(n), np.concatenate([q, np.full(len(q), 100.0), q]), lm, np.ones(n, bool),
               np.arange(n, dtype=np.int64), np.arange(3) * len(q), np.full(3, len(q)), 3)
    used = P.pm_moments(f, np.array([1.0]), P.PMParams(3, 16, 128, 50))["used"]
    assert used.tolist() == [True, False, True]  # the all-100 day has one realised bin
    mo = P.psi_moments(f)
    assert mo["n_real"][0] >= 3 and np.sum(np.isfinite(mo["lmed"][0])) == mo["n_real"][0]


def check_coin_informed():
    """W-null's informed runs: a random 10 % with fair-coin signs; every other run keeps its real sign."""
    eps = np.ones(200_000)
    out, inf = W.coin_informed(eps, np.random.default_rng(0))
    assert abs(inf.mean() - W.GF_F_INF) < 0.005
    assert np.all(out[~inf] == 1.0) and abs(np.mean(out[inf] == -1.0) - 0.5) < 0.01


def stock_flow(f: P.Flow, s: int) -> P.Flow:
    D_ = f.n_days
    a, b = int(f.starts[s * D_]), int(f.starts[(s + 1) * D_ - 1] + f.lengths[(s + 1) * D_ - 1])
    lengths = f.lengths[s * D_ : (s + 1) * D_]
    return P.Flow(f.eps[a:b], f.q[a:b], f.lm[a:b], f.ok[a:b], f.ts[a:b], np.concatenate([[0], np.cumsum(lengths)[:-1]]),
                  lengths, D_)


def check_long_unit_matches_run_unit():
    """diag3d2b's long-day unit (inputs computed one stock at a time, run_unit fed through patched_inputs)
    gives exactly run_unit's p-values and a*; the grid patch is restored."""
    flows = {"Q-T3S1": g0_cell_flow(6, 4, 900, 1)}
    u = P.cell_unit(flows, "Q-T3S1")
    pm = P.PMParams(3, 16, 128, 50)
    want = C.PRIMARY_TESTS
    parts = [DB.stock_parts(stock_flow(u.flow, s), C.CONFIRMATORY["Q-T3S1"], np.random.default_rng(s), [pm], n_sur=2)
             for s in range(u.flow.n_stocks)]
    u2, inputs = DB.assemble(parts, u.flow.n_days, C.CONFIRMATORY["Q-T3S1"])
    sur, P.surrogate_raws = P.surrogate_raws, lambda u, kind, n, rng: inputs["sur"][kind][:n]  # the same draws
    try:
        ref = P.run_unit(u, pm, np.random.default_rng(3), n_boot=30, tests=want)
    finally:
        P.surrogate_raws = sur
    with DB.patched_inputs(inputs):
        got = P.run_unit(u2, pm, np.random.default_rng(3), n_boot=30, tests=want)
    assert P.m3_V is not inputs and callable(P.m3_V) and P.m3_V.__name__ == "m3_V"
    for t in want:
        assert got["tests"][t].p == ref["tests"][t].p, t
    for k in ("estimate", "lo", "hi"):  # equal to rounding (per-stock FFT moments)
        for st in (lambda r: r["tests"]["P-R"].statistic["a_star_T1_T2"],
                   lambda r: r["tests"]["P-R"].secondary["psi_hat"],
                   lambda r: r["tests"]["P-V"].statistic["slope_S2"], lambda r: r["tests"]["P-M"].statistic["slope"]):
            assert np.allclose(st(got)[k], st(ref)[k], rtol=1e-9, atol=1e-12)
    assert np.allclose(DB.a_star_all_T(u2), DB.a_star_all_T(u), rtol=1e-9)
    g0 = C.T_GRID
    with DB.grid(2 ** np.arange(11), ((64, 512), (128, 1024))):
        assert C.T_GRID[-1] == 1024 and P.KW["ranges"] == ((64, 512), (128, 1024))
    assert C.T_GRID is g0 and C.RANGES == ((4, 32), (8, 64)) and P.KW["T_grid"] is g0


# ---- 1 Nov 2026: a* = the local vertex at the grid maximum, censored at a grid end ----------------------------

R_SKEW = np.array([0.857, 0.972, 1.0, 0.983, 0.935])  # 3d-2b's Q-T1 sizes-only R(a) / R(0.5), a = 0 .. 1


def local3(R, k, h=0.25):
    """The parabola vertex through grid points k - 1, k, k + 1, written out."""
    return k * h + 0.5 * h * (R[k - 1] - R[k + 1]) / (R[k - 1] - 2 * R[k] + R[k + 1])


def check_astar_local():
    """The known answer on 3d-2b's skewed curve (0.531, where the global quadratic gives 0.58), a maximum at
    a = 0.25 (neighbours on both sides), within h / 2 of the grid maximum, symmetric curves exact, NaN kept."""
    assert P.vertex(R_SKEW) == pytest.approx(local3(R_SKEW, 2), abs=1e-12)
    assert abs(P.vertex(R_SKEW) - 0.531) < 0.001 and P.vertex_ls(R_SKEW) - P.vertex(R_SKEW) > 0.04
    R1 = np.array([0.95, 1.0, 0.97, 0.9, 0.8])
    assert P.vertex(R1) == pytest.approx(local3(R1, 1), abs=1e-12) and P.censored(R1) == 0
    rng = np.random.default_rng(0)
    for _ in range(200):
        c = rng.uniform(0.2, 0.8)
        R = np.exp(-rng.uniform(0.3, 3) * (C.A_GRID - c) ** 2) + 0.01 * rng.standard_normal(5)
        k = int(np.argmax(R))
        if 0 < k < 4:
            assert abs(P.vertex(R) - C.A_GRID[k]) <= 0.125 + 1e-12 and P.vertex(R) == pytest.approx(local3(R, k))
    for c in (0.3, 0.5, 0.6):  # a symmetric curve's vertex is its centre when the stencil is the parabola itself
        assert P.vertex(1 - (C.A_GRID - c) ** 2) == pytest.approx(c, abs=1e-12)
    assert np.isnan(P.vertex(np.array([0.9, np.nan, 1.0, 0.9, 0.8])))
    Rb = np.stack([R_SKEW, R1, R_SKEW[::-1]])  # batched along leading axes
    assert np.allclose(P.vertex(Rb[None]), [[local3(R_SKEW, 2), local3(R1, 1), 1 - local3(R_SKEW, 2)]])


def check_astar_censored():
    """A maximum at a grid end: a* is that end's a, censored -1 / +1, and P-R fails at that end (p = 1, the
    condition False) even when the replicates sit inside the margin."""
    dec = np.array([1.0, 0.95, 0.85, 0.7, 0.5])
    assert P.vertex(dec) == 0.0 and P.censored(dec) == -1
    assert P.vertex(dec[::-1]) == 1.0 and P.censored(dec[::-1]) == 1
    assert P.censored(R_SKEW) == 0 and P.censored(np.array([np.nan, 1, 0.9, 0.8, 0.7])) == 0
    r = np.random.default_rng(1)
    ends = P.Reps(np.array([0.5, 0.5]), 0.5 + 0.01 * r.standard_normal((200, 2)),
                  0.5 + 0.003 * r.standard_normal((8, 2)), 0.5 + 0.01 * r.standard_normal((200, 2)))
    ok = P.decide_pr(ends, np.array([0, 0]))
    assert ok["p"] < 0.01 and all(ok["conditions"].values()) and ok["censored"] == [None, None]
    for cens, bad in (([1, 0], "T1_within_margin"), ([0, -1], "T2_within_margin")):
        d = P.decide_pr(ends, np.array(cens))
        assert d["p"] == 1.0 and not d["conditions"][bad] and sum(d["conditions"].values()) == 1, (cens, d)


def check_pr_uses_local_vertex():
    """run_unit's P-R a*(T1), a*(T2) are primary3d.vertex of core2's pooled R at those horizons (not core2's
    global a_star), and the censoring reaches the decision."""
    flows = {"Q-T3S1": g0_cell_flow(6, 4, 900, 1)}
    u = P.cell_unit(flows, "Q-T3S1")
    res = P.run_unit(u, P.PMParams(3, 16, 128, 50), np.random.default_rng(3), n_boot=20, tests=("P-R",))
    t = res["tests"]["P-R"]
    S_, D_ = u.flow.Nsd.shape
    o = P.PO.core2(u.raw, np.ones((1, D_)), np.ones((1, S_)), u.specs, "tilde", **P.KW)
    Ra = np.moveaxis(o["R"][0, ..., 0], 0, -1)  # (nT, na)
    for r in C.RANGES:
        ix = [int(np.flatnonzero(C.T_GRID == r[0])[0]), int(np.flatnonzero(C.T_GRID == r[1])[0])]
        rn = P.RANGE_NAMES[C.RANGES.index(r)]
        got = t.by_range[rn]["a_star_T1_T2"]["estimate"]
        assert np.allclose(got, P.vertex(Ra[ix]), rtol=0, atol=1e-12), r
        assert np.min(np.abs(np.asarray(got) - o["a_star"][0, ix, 0])) > 1e-6, r  # the check can tell them apart
        assert t.by_range[rn]["censored"] == [None, None]


# ---- the family: Q-T1 and Q-T3 only (2 Nov 2026) --------------------------------------------------------


def toy_tick_bin(ps: dict, pm_ok: bool, key: str) -> dict:
    """A saved tick-bin record: P-V, P-R, P-F evaluable with the given p and every condition met; P-M at
    design key evaluable or not (p = 0.5)."""
    return {"evaluability": {t: {"evaluable": True, "reason": "reach"} for t in ("P-V", "P-R", "P-F")},
            "tests": {t: {"p": ps.get(t, 0.5), "conditions": {"c": True}, "range_used": [8, 64]}
                      for t in ("P-V", "P-R", "P-F")},
            "pm": {key: {"evaluable": pm_ok, "reason": "reach", **({"p": 0.5, "conditions": {"c": True}}
                                                                     if pm_ok else {})}}}


def check_analysis_family_without_qt2():
    """analyse3d2.family on saved records: Q-T2 is descriptive (rejected, passed None; its raw rate kept),
    the family has the 8 pairs of Q-T1 and Q-T3, and P-R in Q-T1 at p = 0.006 is rejected and passed
    (<= 0.05/8). The former family (_11) holds Q-T2's three evaluable pairs too, m = 11, so the same p is
    not rejected there (0.006 > 0.05/11). correct_rate is None for Q-T2 and scored for the former family."""
    key = W.design_key(D.DEFAULT_PM)
    recs = {"Q-T1": toy_tick_bin({"P-R": 0.006}, True, key), "Q-T2": toy_tick_bin({"P-F": 0.01}, False, key),
            "Q-T3": toy_tick_bin({}, True, key)}
    f = A.family(recs, key)
    assert f["family_size"] == 8 and f["family_size_11"] == 11
    d = f["decisions"]
    assert d[("P-R", "Q-T1")]["rejected"] and d[("P-R", "Q-T1")]["passed"]
    assert not d[("P-R", "Q-T1")]["rejected_11"] and not d[("P-R", "Q-T1")]["passed_11"]
    for t in C.PRIMARY_TESTS:
        assert d[(t, "Q-T2")]["rejected"] is None and d[(t, "Q-T2")]["passed"] is None
        assert d[(t, "Q-T2")]["rejected_11"] is False
    assert d[("P-F", "Q-T2")]["raw"] and not d[("P-M", "Q-T2")]["evaluable"]
    assert A.correct_rate("W0", "P-R", [d[("P-R", "Q-T2")]]) is None
    assert A.correct_rate("W0", "P-R", [d[("P-R", "Q-T2")]], "_11") == (0.0, 0.0, 1)
    assert A.correct_rate("W0", "P-R", [d[("P-R", "Q-T1")]]) == (1.0, 0.0, 1)


# ---- one instance at a time, chunk integrity, trimming ---------------------------------------------------


def _refused(out: Path) -> bool:
    try:
        D.acquire_lock(out).close()
    except SystemExit as e:
        return "another phase3d2 run holds" in str(e)
    return False


def check_lock_refuses_second():
    """While one holder keeps the lock, a second acquire on the same directory exits; once the holder closes,
    or its process dies without releasing, the lock is free again."""
    with tempfile.TemporaryDirectory() as td:
        out = Path(td)
        held = D.acquire_lock(out)
        assert _refused(out)
        assert f"pid {os.getpid()}" in (out / D.LOCK_NAME).read_text()
        held.close()
        assert not _refused(out)
        code = f"import os, phase3d2 as D; from pathlib import Path; D.acquire_lock(Path({td!r})); os._exit(0)"
        subprocess.run([sys.executable, "-c", code], cwd=Path(D.__file__).parent, check=True)
        assert not _refused(out)  # a dead holder leaves no stale lock


def check_main_takes_lock():
    """main refuses to run any step while another instance holds the output directory's lock: in this
    process (D.main) and from a second process (the command line)."""
    with tempfile.TemporaryDirectory() as td:
        out, saved = Path(td), D.OUT
        held = D.acquire_lock(out)
        try:
            D.OUT = out
            try:
                D.main(["check"])
                refused = False
            except SystemExit as e:
                refused = "another phase3d2 run holds" in str(e)
            assert refused and not (out / "study").exists()
        finally:
            D.OUT = saved
        env = {**os.environ, "PHASE3D2_OUT": td}
        r = subprocess.run([sys.executable, "phase3d2.py", "check"], cwd=Path(D.__file__).parent, env=env,
                           capture_output=True, text=True)
        assert r.returncode != 0 and "another phase3d2 run holds" in r.stderr, r.stderr
        held.close()
        r = subprocess.run([sys.executable, "phase3d2.py", "check"], cwd=Path(D.__file__).parent, env=env,
                           capture_output=True, text=True)
        assert r.returncode == 0 and "0 quarantined" in r.stdout, r.stderr


def toy_chunk(start: int, keys, tb: str = "Q-T1") -> list[dict]:
    ev = {"P-V": {"evaluable": True}, "P-R": {"evaluable": False}, "P-F": {"evaluable": True}}
    return [{"tick_bin": tb, "n_stocks": 2, "Ns": [10, 20], "rho": [0.1, 0.2], "evaluability": ev,
             "tests": {"P-V": {"p": 0.3}, "P-F": {"p": 1.0}},
             "pm": {k: {"evaluable": True, "p": 0.01 * (j + 1), "conditions": {"slope": True}}
                    for j, k in enumerate(keys)}, "panel": i, "seconds": {"generate": 1.0}}
            for i in range(start, start + D.CHUNK_PANELS)]


def check_verify_chunk():
    """A complete chunk passes (6 designs, or all 18 of the old grid); a truncated one, one with a wrong
    panel, a missing design, an evaluable test without p, a name that is no job and a leftover .tmp are all
    caught, and check_study moves each of them, and only them, to quarantine."""
    keep = [W.design_key(d) for d in W.PM_DESIGNS]
    old = [W.design_key(d) for d in W.PM_DESIGNS_18]
    with tempfile.TemporaryDirectory() as td:
        st = Path(td) / "study"
        st.mkdir()
        good = {"W0_Q-T1_004.json": toy_chunk(4, keep), "W-null_Q-T1_000.json": toy_chunk(0, old),
                "W0-psi3_Q-T3_002.json": toy_chunk(2, keep[:1], "Q-T3")}
        for name, recs in good.items():
            (st / name).write_text(json.dumps(recs))
            assert D.verify_chunk(st / name) == [], name
        bad = {}
        txt = json.dumps(toy_chunk(6, keep))
        bad["W0_Q-T1_006.json"] = txt[: len(txt) // 2]
        r = toy_chunk(8, keep)
        r[1]["panel"] = 8
        bad["W0_Q-T1_008.json"] = json.dumps(r)
        r = toy_chunk(10, keep)
        del r[0]["pm"][keep[3]]
        bad["W0_Q-T1_010.json"] = json.dumps(r)
        r = toy_chunk(12, keep)
        r[1]["tests"]["P-F"]["p"] = None
        bad["W0_Q-T1_012.json"] = json.dumps(r)
        bad["W0_Q-T1_004 2.json"] = json.dumps(toy_chunk(4, keep))
        bad["W0_Q-T1_014.json"] = json.dumps(toy_chunk(14, keep)[:1])
        for name, t in bad.items():
            (st / name).write_text(t)
            assert D.verify_chunk(st / name), name
        (st / "W0_Q-T1_016.tmp").write_text("[")
        res = D.check_study(Path(td), log=lambda *a, **k: None)
        assert res["complete"] == len(good) and res["quarantined"] == len(bad) + 1
        assert sorted(p.name for p in st.iterdir()) == sorted(good)
        q = Path(res["quarantine"])
        assert sorted(p.name for p in q.iterdir() if p.name != "reasons.json") == sorted([*bad, "W0_Q-T1_016.tmp"])


def check_trim_keeps_retained():
    """trim keeps exactly the 6 retained designs of an 18-design chunk, every value unchanged, archives the
    original byte for byte, leaves a W0-psi3 chunk and an already trimmed chunk alone, and is idempotent."""
    keep = [W.design_key(d) for d in W.PM_DESIGNS]
    old = [W.design_key(d) for d in W.PM_DESIGNS_18]
    with tempfile.TemporaryDirectory() as td:
        out = Path(td)
        st = out / "study"
        st.mkdir()
        orig = toy_chunk(0, old)
        (st / "W0_Q-T1_000.json").write_text(json.dumps(orig))
        (st / "WB_Q-T1_002.json").write_text(json.dumps(toy_chunk(2, keep)))
        (st / "W0-psi3_Q-T1_000.json").write_text(json.dumps(toy_chunk(0, keep[:1])))
        before = {p.name: p.read_bytes() for p in st.iterdir()}
        assert D.trim_designs(out, log=lambda *a, **k: None) == 1
        got = json.loads((st / "W0_Q-T1_000.json").read_text())
        for g, o in zip(got, orig, strict=True):
            assert list(g["pm"]) == keep and all(g["pm"][k] == o["pm"][k] for k in keep)
            assert {k: v for k, v in g.items() if k != "pm"} == {k: v for k, v in o.items() if k != "pm"}
        assert (out / "study_18designs" / "W0_Q-T1_000.json").read_bytes() == before["W0_Q-T1_000.json"]
        for name in ("WB_Q-T1_002.json", "W0-psi3_Q-T1_000.json"):
            assert (st / name).read_bytes() == before[name]
        assert D.trim_designs(out, log=lambda *a, **k: None) == 0
        assert D.verify_chunk(st / "W0_Q-T1_000.json") == []


CHECKS = {f.__name__: f for f in [
    check_world_labels, check_world_specs, check_design_rule, check_design_rule_ties_and_fallback,
    check_router_share_once, check_thin_settings_paired, check_setting_grid, check_burn_from_other_days,
    check_wb_beta_and_impact, check_semi_response_direct, check_semi_flow_unchanged, check_wm_burn,
    check_tick_bin_matches_run_primary, check_m3_short_days, check_tie_bins, check_coin_informed,
    check_long_unit_matches_run_unit, check_astar_local, check_astar_censored, check_pr_uses_local_vertex,
    check_design_grid, check_analysis_family_without_qt2, check_lock_refuses_second, check_main_takes_lock,
    check_verify_chunk, check_trim_keeps_retained,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_3d2(name):
    CHECKS[name]()


class Bug(NamedTuple):
    id: str
    module: object
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("design rule's objective read from W0", W, "choose_design", "rates[(POWER_WORLD, k, tb)]",
        'rates[("W0", k, tb)]', ("check_design_rule",)),
    Bug("design rule's size constraint read from WB", W, "choose_design", "for wd in SIZE_WORLDS:",
        "for wd in (POWER_WORLD,):", ("check_design_rule",)),
    Bug("design rule's size constraint on W-null only", W, "choose_design", "for wd in SIZE_WORLDS:",
        'for wd in ("W-null",):', ("check_design_rule",)),
    Bug("design ties to more bins", W, "choose_design", "key=lambda k: (-obj[k], order[k])",
        "key=lambda k: (-obj[k], tuple(-x for x in order[k]))", ("check_design_rule_ties_and_fallback",)),
    Bug("hidden share applied twice (router_p ignores h)", W, "router_p", "p = share / (1.0 - h)", "p = share",
        ("check_router_share_once",)),
    Bug("hidden share applied twice (thin_day)", W, "thin_day", "router_p(share, h), rng, h, kappa",
        "router_p(share, h) * (1 - h), rng, h, kappa", ("check_router_share_once",)),
    Bug("only the kept settings thinned, from the market generator", W, "wm_cell_flows",
        "obs = [thin_day(day, s, r) for s, r in zip(settings, thin_rngs, strict=True)]",
        "obs = [thin_day(day, settings[j], rng) if j in keep else day for j in range(len(settings))]",
        ("check_thin_settings_paired",)),
    Bug("pass and fail swapped in the P-R label", A, "correct_rate", "return rate([not d[pas] for d in decs])",
        "return rate([d[pas] for d in decs])", ("check_world_labels",)),
    Bug("reject counted without the conditions", A, "correct_rate",
        'if lab == "reject":\n        return rate([d[pas] for d in decs])',
        'if lab == "reject":\n        return rate([d[rej] for d in decs])', ("check_world_labels",)),
    Bug("W0-psi3 generated at psi = 0.5", W, "world_spec",
        '"W0-psi3": {"flow": "real", "signs": "real", "gen": "G0", "psi": 0.3}',
        '"W0-psi3": {"flow": "real", "signs": "real", "gen": "G0", "psi": PSI}', ("check_world_specs",)),
    Bug("burn-in includes the day itself", W, "burn_index", "if j != d", "if True", ("check_burn_from_other_days",)),
    Bug("WB anchors reversed", W, "wb_beta", "b0, b1 = betas", "b1, b0 = betas", ("check_wb_beta_and_impact",)),
    Bug("burn-in kept in the day's mid", W, "semi_response", "m[sst[k] + burn : sst[k] + ser[k]]",
        "m[sst[k] : sst[k] + ser[k] - burn]", ("check_semi_response_direct",)),
    Bug("P-M draws not run_unit's", D, "run_tick_bin", "dr = P.draws(S_, D_, rng_dr, n_boot)",
        "dr = P.draws(S_, D_, rng, n_boot)", ("check_tick_bin_matches_run_primary",)),
    Bug("P-M result reused across designs whatever stock-days they keep", D, "run_tick_bin",
        'P.pm_moments(u.flow, u.rho, d)["used"].tobytes())', "b'')", ("check_tick_bin_matches_run_primary",)),
    Bug("W-null sizes shuffled across days", W, "shuffle_sizes", "out[a : a + n] = SU.s1(q[a : a + n], rng)",
        "out[a : a + n] = SU.s1(q, rng)[:n]", ("check_semi_flow_unchanged",)),
    Bug("m3_V ACF to lag 63 whatever the day length", P, "m3_V", "lag = min(int(T[-1]) - 1, int(flow.lengths[i]) - 1)",
        "lag = int(T[-1]) - 1", ("check_m3_short_days",)),
    Bug("a size value split across two bins (rank cuts)", P, "tie_bins",
        'return np.searchsorted(upper, q, side="left"), upper',
        'return np.minimum(np.argsort(np.argsort(q, kind="stable")) * (len(upper) + 1) // len(q), len(upper)), upper',
        ("check_tie_bins",)),
    Bug("empty bin kept (coinciding cuts kept, no merge)", P, "tie_bins",
        "cuts = set()\n    for k in range(1, n_target):\n        if len(cum):\n            cuts.add(int(np.argmin("
        "np.abs(cum - k * len(q) / n_target))))\n    ends = sorted(cuts) + [len(v) - 1]  # last distinct-value index"
        " of each bin\n    sizes = np.diff(np.concatenate([[0], np.cumsum(cnt)[ends]])).tolist()\n"
        "    while len(sizes) > 1 and min(sizes) < min_runs:",
        "cuts = []\n    for k in range(1, n_target):\n        if len(cum):\n            cuts.append(int(np.argmin("
        "np.abs(cum - k * len(q) / n_target))))\n    ends = sorted(cuts) + [len(v) - 1]\n    sizes = []\n"
        "    while False:",
        ("check_tie_bins",)),
    Bug("P-M uses a stock-day below the target bin count", P, "pm_moments",
        "if len(upper) != nb - 1:", "if len(upper) > nb - 1:", ("check_tie_bins",)),
    Bug("W-null informed runs keep their real signs", W, "coin_informed", "np.where(informed, coin, eps)",
        "np.where(informed, eps, eps)", ("check_coin_informed",)),
    Bug("W-null spec without the coin-flip signs", W, "world_spec",
        '"W-null": {"flow": "real_s1", "signs": "coin_informed"', '"W-null": {"flow": "real_s1", "signs": "real"',
        ("check_world_specs",)),
    Bug("long unit: surrogate draws for P-V from the wrong kind", DB, "patched_inputs",
        'inputs["sur"][kind][:n]', 'inputs["sur"]["S1"][:n]', ("check_long_unit_matches_run_unit",)),
    Bug("long unit: the grid patch not restored", DB, "grid",
        "C.T_GRID, C.RANGES, C.PRIMARY_RANGE, C.FALLBACK_RANGE, P.KW, P.RANGE_NAMES = saved", "pass",
        ("check_long_unit_matches_run_unit",)),
    Bug("a*: the global quadratic used (vertex)", P, "vertex", "return np.where(fin, ast, np.nan)",
        "return vertex_ls(R, a_grid)", ("check_astar_local",)),
    Bug("a*: the global quadratic used (P-R reads core2's a_star)", P, "run_unit",
        "ast, cen = vertex(Ra), censored(Ra)", 'ast, cen = o["a_star"][..., 0], censored(Ra)',
        ("check_pr_uses_local_vertex",)),
    Bug("a*: neighbours taken from one side of the maximum", P, "vertex", "c = np.clip(k, 1, len(a) - 2)",
        "c = np.clip(k + 1, 1, len(a) - 2)", ("check_astar_local",)),
    Bug("a*: neighbours swapped (vertex on the wrong side)", P, "vertex", "0.5 * h * (Rm - Rp) / den",
        "0.5 * h * (Rp - Rm) / den", ("check_astar_local",)),
    Bug("a*: censored end treated as interior (vertex)", P, "vertex",
        "ast = np.where(k == 0, a[0], np.where(k == len(a) - 1, a[-1], inner))", "ast = inner",
        ("check_astar_censored",)),
    Bug("a*: censored end treated as interior (P-R decision)", P, "decide_pr",
        "cens = np.zeros(2) if cens is None else np.asarray(cens)", "cens = np.zeros(2)",
        ("check_astar_censored",)),
    Bug("analysis family over every tick bin (Q-T2 in Holm)", A, "family",
        "fam = P.combine(ev, {k: v for k, v in pv.items() if k[1] in C.FAMILY_TICK_BINS})",
        'fam = {"holm": P.step_down(pv, ALPHA), "family": list(pv)}', ("check_analysis_family_without_qt2",)),
    Bug("Q-T2 scored as part of the family", A, "family",
        '"rejected": bool(h and h["rejected"]) if in_fam else None,', '"rejected": bool(h and h["rejected"]),',
        ("check_analysis_family_without_qt2",)),
    Bug("lock never taken (no flock)", D, "acquire_lock", "fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)", "pass",
        ("check_lock_refuses_second", "check_main_takes_lock")),
    Bug("lock taken blocking-free but shared", D, "acquire_lock", "fcntl.LOCK_EX | fcntl.LOCK_NB",
        "fcntl.LOCK_SH | fcntl.LOCK_NB", ("check_lock_refuses_second", "check_main_takes_lock")),
    Bug("main does not take the lock", D, "main", "lock = acquire_lock(OUT)", "lock = None",
        ("check_main_takes_lock",)),
    Bug("chunk's panels not checked", D, "verify_chunk", 'if r["panel"] != start + j or r["tick_bin"] != tb:',
        "if False:", ("check_verify_chunk",)),
    Bug("a bad chunk reported but left in place", D, "check_study",
        "shutil.move(str(d / name), str(q / name))", "pass", ("check_verify_chunk",)),
    Bug("trim keeps the first 6 of the old grid", D, "trim_designs", "keep = [W.design_key(d) for d in W.PM_DESIGNS]",
        "keep = [W.design_key(d) for d in W.PM_DESIGNS_18[:6]]", ("check_trim_keeps_retained",)),
    Bug("trim does not archive the original", D, "trim_designs", "shutil.copy2(f, arch / f.name)", "pass",
        ("check_trim_keeps_retained",)),
]
DATA_BUGS = [
    ("WB's P-M labelled no reject", ("WB", "P-M", "no reject")),
    ("W-null's P-F unlabelled", ("W-null", "P-F", None)),
    ("W0's P-R labelled fail", ("W0", "P-R", "fail")),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(bug.module, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()


@pytest.mark.parametrize("bug", DATA_BUGS, ids=[b[0] for b in DATA_BUGS])
def test_planted_label_bug_is_caught(bug, monkeypatch):
    world, test, label = bug[1]
    monkeypatch.setitem(W.CORRECT[world], test, label)
    with pytest.raises(AssertionError):
        check_world_labels()
