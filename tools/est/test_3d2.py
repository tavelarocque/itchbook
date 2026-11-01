"""Tests for Phase 3d-2 (worlds3d2, phase3d2, analyse3d2): world labels and specs, the P-M design rule on the
right worlds, thinning shares applied once and paired settings, the burn-in drawn from the stock's other days,
the WB kernel, the response against a direct sum, and the driver's tick-bin run against run_primary. Every
input is synthetic; no data file is read (the WM test replaces panel3c.day_lengths).

    uv run --group dev --group analysis pytest tools/est/test_3d2.py
"""
from __future__ import annotations

import copy
from typing import NamedTuple

import numpy as np
import pytest

import analyse3d2 as A
import gen3 as G
import mbmodel as M
import panel3c as PC
import phase3d2 as D
import prereg3d as C
import primary3d as P
import worlds3d2 as W

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
    assert W.world_spec("W-null") == {"flow": "real_s1", "gen": "GF", "psi": 0.5}
    assert W.world_spec("W0") == {"flow": "real", "gen": "G0", "psi": 0.5}
    assert W.world_spec("W0-psi3") == {"flow": "real", "gen": "G0", "psi": 0.3}
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
    k_best, k_w0, k_null, k_one = "b5_l32_L256_m20", "b3_l16_L128_m50", "b4_l16_L128_m20", "b4_l32_L256_m50"
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
    assert out["chosen"] == "b3_l16_L128_m20"  # every rate equal: fewer bins, then l_max, then min runs
    tie = {"b5_l16_L128_m20": 0.7, "b4_l32_L256_m100": 0.7}
    assert W.choose_design(toy_rates(tie, {}), N200)["chosen"] == "b4_l32_L256_m100"
    worst = toy_rates({}, {k: 0.3 for k in keys})
    for tb in C.TICK_BINS:
        for w in W.SIZE_WORLDS:
            worst[(w, "b5_l32_L256_m100", tb)] = 0.12
    out = W.choose_design(worst, N200)
    assert out["chosen"] is None and out["fallback"] == "b5_l32_L256_m100"


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



CHECKS = {f.__name__: f for f in [
    check_world_labels, check_world_specs, check_design_rule, check_design_rule_ties_and_fallback,
    check_router_share_once, check_thin_settings_paired, check_setting_grid, check_burn_from_other_days,
    check_wb_beta_and_impact, check_semi_response_direct, check_semi_flow_unchanged, check_wm_burn,
    check_tick_bin_matches_run_primary, check_m3_short_days,
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
    Bug("pass and fail swapped in the P-R label", A, "correct_rate", 'return rate([not d["passed"] for d in decs])',
        'return rate([d["passed"] for d in decs])', ("check_world_labels",)),
    Bug("reject counted without the conditions", A, "correct_rate",
        'if lab == "reject":\n        return rate([d["passed"] for d in decs])',
        'if lab == "reject":\n        return rate([d["rejected"] for d in decs])', ("check_world_labels",)),
    Bug("W0-psi3 generated at psi = 0.5", W, "world_spec", '"W0-psi3": {"flow": "real", "gen": "G0", "psi": 0.3}',
        '"W0-psi3": {"flow": "real", "gen": "G0", "psi": PSI}', ("check_world_specs",)),
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
