"""Tests for Phase 3d-1 (primary3d, prereg3d): known-answer inputs for each primary test, Holm over exactly the
four p-values in the right order, the [8, 64] to [4, 32] fallback by the reach rule only, descriptive cells
kept out, both surrogates for P-V (ii), both conditions for P-M, and P-R as the equivalence test of
31 Oct 2026. Planted bugs through the mutate fixture. Every input is synthetic; no data file is read.

    uv run --group dev --group analysis pytest tools/est/test_3d1.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import gen3 as G
import moments as MO
import panel3c as PC
import pooled as PO
import prereg3d as C
import primary3d as P

CELL = G.Cell("TEST", np.arange(100, 200), 3.0, 1.0, (0.3, 0.7), {"w1": 0.3, "H": 0.8, "w2": 0.0, "phi": 0.5})
PM = P.PMParams(3, 16, 128, 50)
NA, NT = len(C.A_GRID), len(C.T_GRID)
T = C.T_GRID.astype(float)


def even_ts(lengths) -> np.ndarray:
    return np.concatenate([34200 * 10**9 + np.arange(n) * (23400 * 10**9 // n) for n in lengths]).astype(np.int64)


def g0_flow(S: int, D: int, N: int, seed: int) -> P.Flow:
    """G0 flow (psi = 0.5): log-normal sizes independent of the signs."""
    p = PC.make_panel("G0", CELL, np.full((S, D), N), np.random.default_rng(seed), params={"beta": 0.3}, clip=None)
    return P.Flow(p.eps, p.q, p.lm, np.ones(len(p.eps), bool), even_ts(p.lengths), p.starts, p.lengths, D)


def make_flow(eps, q, lm, Nsd) -> P.Flow:
    lengths = np.asarray(Nsd).ravel().astype(np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)[:-1]]).astype(np.int64)
    return P.Flow(np.asarray(eps, float), np.asarray(q, float), np.asarray(lm, float), np.ones(len(eps), bool),
                  even_ts(lengths), starts, lengths, np.asarray(Nsd).shape[1])


def reps(est, sd, rng, B=200, D=8) -> P.Reps:
    """Replicates scattered around est (array) with standard deviation sd."""
    est = np.asarray(est, dtype=float)
    return P.Reps(est, est + sd * rng.standard_normal((B, *est.shape)),
                  est + sd / np.sqrt(D) * rng.standard_normal((D, *est.shape)),
                  est + sd * rng.standard_normal((B, *est.shape)))


# ---- config, pooling, clip ----------------------------------------------------------------------------------


def check_config_table():
    """The pre-registered table, in one module: six confirmatory cells with their ceilings, three descriptive."""
    assert C.CONFIRMATORY == {"Q-T1S1": 1 / 5, "Q-T1S2": 1 / 20, "Q-T2S3": 1 / 5, "Q-T3S1": 1 / 5,
                              "Q-T3S2": 1 / 5, "Q-T3S3": 1 / 10}
    assert set(C.DESCRIPTIVE_CELLS) == {"Q-T1S3", "Q-T2S1", "Q-T2S2"}
    assert not set(C.DESCRIPTIVE_CELLS) & set(C.CONFIRMATORY)
    assert C.PRIMARY_TESTS == ("P-V", "P-R", "P-M", "P-F") and C.ALPHA_FW == 0.05
    assert np.allclose(C.A_GRID, [0, 0.25, 0.5, 0.75, 1]) and C.RANGES == ((4, 32), (8, 64))


def check_descriptive_excluded():
    """A descriptive cell in the input never enters the tick bin's pool, even poisoned, and cannot be pooled
    on its own."""
    good = {"Q-T1S1": g0_flow(3, 2, 400, 1), "Q-T1S2": g0_flow(2, 2, 400, 2)}
    bad = g0_flow(4, 2, 400, 3)
    bad.lm[:] = np.nan
    u = P.tick_bin_unit({**good, "Q-T1S3": bad}, "Q-T1")
    assert u.flow.n_stocks == 5 and set(u.cells) == {"Q-T1S1", "Q-T1S2"}
    assert u.excluded_cells == ("Q-T1S3",)
    assert np.all(np.isfinite(u.flow.lm))
    ref = P.tick_bin_unit(good, "Q-T1")
    assert np.array_equal(u.raw.V, ref.raw.V, equal_nan=True)
    for c in C.DESCRIPTIVE_CELLS:
        with pytest.raises(ValueError):
            P.cell_unit({c: bad}, c)


def check_cell_ceilings():
    """In a tick-bin pool each stock keeps its own cell's ceiling: with N = 1000, a Q-T1S1 stock (1/5) is
    eligible up to T = 128 > 64, a Q-T1S2 stock (1/20) only up to 50, so T = 64 is out for it and T = 32 in."""
    u = P.tick_bin_unit({"Q-T1S1": g0_flow(2, 2, 1000, 1), "Q-T1S2": g0_flow(3, 2, 1000, 2)}, "Q-T1")
    assert np.allclose(u.rho, [1 / 5] * 2 + [1 / 20] * 3)
    M = PO.masks(u.Ns, np.array([32, 64]), u.specs)[:, :, 0]
    assert M[:2].all() and M[2:, 0].all() and not M[2:, 1].any()
    _, by = P.select_range(u.Ns, u.rho)
    assert by == {(4, 32): 5, (8, 64): 2}


def check_r2_clip():
    """The unit's sizes are clipped at each stock's own 99.9th percentile over its days; the raw sizes stay."""
    f = g0_flow(3, 2, 2000, 4)
    u = P.cell_unit({"Q-T3S1": f}, "Q-T3S1")
    for s in range(3):
        sel = slice(int(f.starts[2 * s]), int(f.starts[2 * s + 1] + f.lengths[2 * s + 1]))
        cap = np.quantile(f.q[sel], 0.999)
        assert np.isclose(u.flow.q[sel].max(), cap) and np.allclose(u.flow.q[sel], np.minimum(f.q[sel], cap))
        assert f.q[sel].max() > cap


# ---- range rule ---------------------------------------------------------------------------------------------


def check_fallback_by_reach():
    """[8, 64] with 5 stocks at 64; [4, 32] with 4 at 64 and 5 at 32; none with 4 at 32."""
    rho = np.full(10, 1 / 5)
    Ns = np.array([320] * 5 + [100] * 5)  # 320 / 5 = 64
    assert P.select_range(Ns, rho)[0] == (8, 64)
    Ns = np.array([320] * 4 + [160] * 1 + [100] * 5)  # 160 / 5 = 32
    assert P.select_range(Ns, rho)[0] == (4, 32)
    Ns = np.array([320] * 4 + [100] * 6)
    r, by = P.select_range(Ns, rho)
    assert r is None and by == {(4, 32): 4, (8, 64): 4}


def check_fallback_never_by_result():
    """The test reports the range the reach rule picked even when the other range has the smaller p; both
    ranges are kept in by_range. The rule's inputs are day lengths and ceilings only."""
    import inspect

    assert list(inspect.signature(P.select_range).parameters) == ["Ns", "rho"]
    by_r = {(4, 32): {"p": 0.001, "alternative": "two-sided", "conditions": {"interval_excludes_0": True}},
            (8, 64): {"p": 0.4, "alternative": "two-sided", "conditions": {"interval_excludes_0": False}}}
    t = P._assemble("P-F", by_r, (8, 64), {(4, 32): 9, (8, 64): 6}, {(4, 32): [], (8, 64): []},
                    ("interval_excludes_0",))
    assert t.p == 0.4 and t.range_used == (8, 64) and t.range_rule.startswith("primary")
    assert set(t.by_range) == {"[4,32]", "[8,64]"}
    t = P._assemble("P-F", by_r, (4, 32), {(4, 32): 9, (8, 64): 3}, {(4, 32): [], (8, 64): []},
                    ("interval_excludes_0",))
    assert t.p == 0.001 and t.range_rule.startswith("fallback")


# ---- replicates, intervals, p-values ----------------------------------------------------------------------------


def check_interval_matches_intervals3():
    """The generic widest-of-three interval equals pooled.intervals3 on the same draws (same seed)."""
    u = P.cell_unit({"Q-T3S1": g0_flow(6, 4, 800, 5)}, "Q-T3S1")
    S_, D_ = u.flow.Nsd.shape
    dr = P.draws(S_, D_, np.random.default_rng(7), 40)
    r = P.run_stats(lambda Wd, Ws: {"pf": PO.core2(u.raw, Wd, Ws, u.specs, "tilde", **P.KW)["pf"]}, dr)
    lo, hi, m = P.interval(r["pf"])
    ref = PO.intervals3(u.raw, u.specs, np.random.default_rng(7), n_boot=40, keys=("pf",), **P.KW)["pf"]
    fin = np.isfinite(ref[1])
    assert fin.any()
    assert np.allclose(lo[fin], ref[1][fin]) and np.allclose(hi[fin], ref[2][fin]) and np.array_equal(m, ref[3])


def check_p_value_conservative():
    """The p-value is the largest of the three methods', and 1 for an undefined estimate."""
    rng = np.random.default_rng(1)
    r = reps(np.array([-0.1]), 0.05, rng)
    r = P.Reps(r.est, r.ds, r.jk, r.est + 0.2 * rng.standard_normal((200, 1)))  # whole-stock reps much wider
    ps = [P._boot_p(r.ds, 0.0, "less"), P._jk_p(r.est, r.jk, 0.0, "less"), P._boot_p(r.ws, 0.0, "less")]
    assert np.isclose(P.p_value(r, 0.0, "less")[0], max(float(x[0]) for x in ps))
    assert P.p_value(r, 0.0, "less")[0] > 0.2
    assert P.p_value(P.Reps(np.array([np.nan]), r.ds, r.jk, r.ws), 0.0, "less")[0] == 1.0


# ---- Holm ----------------------------------------------------------------------------------------------------


def fam(ps: dict, tb: str = "Q-T1") -> dict:
    return {(k, tb): v for k, v in ps.items()}


def all_ev(left=()) -> dict:
    """Evaluability of the 12 pairs of the three tick bins (Q-T2's are descriptive), every pair evaluable
    except those in left."""
    return {tb: {t: {"evaluable": (t, tb) not in left, "reason": "reach" if (t, tb) not in left else "below reach"}
                 for t in C.PRIMARY_TESTS} for tb in C.TICK_BINS}


def check_holm_order():
    """Ascending p against 0.05/4, 0.05/3, 0.05/2, 0.05, stopping at the first failure: P-R (0.001) and P-F
    (0.015 <= 0.0167, above the Bonferroni 0.0125) are rejected; P-M (0.03 > 0.025) stops, so P-V (0.04)
    is not rejected although 0.04 <= 0.05."""
    h = P.holm(fam({"P-V": 0.04, "P-R": 0.001, "P-M": 0.03, "P-F": 0.015}))
    assert [k for k in C.PRIMARY_TESTS if h[(k, "Q-T1")]["rejected"]] == ["P-R", "P-F"]
    assert [h[(k, "Q-T1")]["rank"] for k in ("P-R", "P-F", "P-M", "P-V")] == [1, 2, 3, 4]
    assert np.isclose(h[("P-F", "Q-T1")]["threshold"], 0.05 / 3) and np.isclose(h[("P-V", "Q-T1")]["threshold"], 0.05)
    assert np.allclose([h[(k, "Q-T1")]["p_adj"] for k in ("P-R", "P-F", "P-M", "P-V")], [0.004, 0.045, 0.06, 0.06])


def check_holm_family_keys():
    """Members must be (primary test, Q-T1 or Q-T3) pairs: P-K, P-C, a bare test name, an unknown tick bin or
    a Q-T2 pair (descriptive since 2 Nov 2026) raise; NaN is 1."""
    four = fam({"P-V": 0.01, "P-R": 0.01, "P-M": 0.01, "P-F": 0.01})
    for bad in ({**four, ("P-K", "Q-T1"): 0.001}, {**four, ("P-C", "Q-T3"): 0.001}, {**four, "P-V": 0.001},
                {**four, ("P-V", "Q-T4"): 0.001}, {**four, ("P-V", "Q-T2"): 0.001}):
        try:
            P.holm(bad)
        except ValueError:
            continue
        raise AssertionError(f"Holm accepted the family {list(bad)}")
    h = P.holm({**four, ("P-V", "Q-T1"): np.nan})
    assert h[("P-V", "Q-T1")]["p"] == 1.0 and not h[("P-V", "Q-T1")]["rejected"]


def check_one_family_across_bins():
    """One Holm family over the 8 pairs of Q-T1 and Q-T3: one p = 0.01 per tick bin, the rest 0.5, rejects
    nothing (0.01 > 0.05/8), where per-tick-bin Holm (0.01 <= 0.05/4) would reject both."""
    ps = {(t, tb): (0.01 if t == "P-F" else 0.5) for tb in C.FAMILY_TICK_BINS for t in C.PRIMARY_TESTS}
    out = P.combine(all_ev(), ps)
    assert len(out["family"]) == 8 and out["left"] == []
    assert not any(r["rejected"] for r in out["holm"].values())
    assert np.isclose(out["holm"][("P-F", "Q-T1")]["threshold"], 0.05 / 8)
    ps[("P-F", "Q-T3")] = 0.006  # <= 0.05/8: rejected, then 0.01 > 0.05/7 stops
    out = P.combine(all_ev(), ps)
    assert [k for k, r in out["holm"].items() if r["rejected"]] == [("P-F", "Q-T3")]


def check_qt2_never_in_holm():
    """Q-T2 is descriptive (2 Nov 2026): Holm refuses a Q-T2 pair; combine never makes one a member, lists
    each under 'descriptive' with its evaluability, and refuses a Q-T2 p-value. With Q-T2's pairs evaluable,
    the family is still the 8 pairs: p = 0.006 in Q-T1 is rejected (<= 0.05/8), which a family also holding
    Q-T2's three evaluable pairs at p = 0.5 (m = 11, 0.05/11 = 0.0045) would not do."""
    try:
        P.holm({("P-V", "Q-T2"): 0.001})
        raise AssertionError("Holm accepted a Q-T2 pair")
    except ValueError:
        pass
    ev = all_ev(left=(("P-M", "Q-T2"),))
    ps = {(t, tb): 0.5 for tb in C.FAMILY_TICK_BINS for t in C.PRIMARY_TESTS}
    ps[("P-R", "Q-T1")] = 0.006
    out = P.combine(ev, ps)
    assert len(out["family"]) == 8 and not any(k[1] == "Q-T2" for k in out["family"])
    assert not any(k[1] == "Q-T2" for k in out["holm"]) and out["left"] == []
    assert out["holm"][("P-R", "Q-T1")]["rejected"]
    assert {(d["test"], d["tick_bin"]): d["evaluable"] for d in out["descriptive"]} == {
        ("P-V", "Q-T2"): True, ("P-R", "Q-T2"): True, ("P-M", "Q-T2"): False, ("P-F", "Q-T2"): True}
    try:
        P.combine(ev, {**ps, ("P-V", "Q-T2"): 0.0})
        raise AssertionError("combine accepted a Q-T2 p-value")
    except ValueError:
        pass


def check_nonevaluable_leaves_family():
    """Two pairs below the reach leave the family before any p-value (m = 6, logged with their reason):
    p = 0.008 <= 0.05/6 is rejected; kept with p = 1 (m = 8, 0.05/8 = 0.00625) it would not be. A p-value
    supplied for a pair that left is refused."""
    left = (("P-M", "Q-T1"), ("P-V", "Q-T1"))
    ev = all_ev(left)
    ps = {(t, tb): 0.5 for tb in C.FAMILY_TICK_BINS for t in C.PRIMARY_TESTS if (t, tb) not in left}
    ps[("P-R", "Q-T3")] = 0.008
    out = P.combine(ev, ps)
    assert len(out["family"]) == 6 and out["holm"][("P-R", "Q-T3")]["rejected"]
    assert {(r["test"], r["tick_bin"]) for r in out["left"]} == set(left)
    assert all(r["reason"] == "below reach" for r in out["left"]) and not set(left) & set(out["holm"])
    try:
        P.combine(ev, {**ps, ("P-M", "Q-T1"): 1.0})
    except ValueError:
        return
    raise AssertionError("a p-value for a pair outside the family was accepted")


def check_evaluability_reach_only():
    """Evaluability comes from day lengths and ceilings alone: an evaluable pair whose statistic turned out
    undefined (p = 1) stays in the family, so m = 8 and 0.0065 > 0.05/8 is not rejected (with m = 7 it would
    be). The evaluability function sees only Ns, rho and the P-M parameters; P-M needs MIN_REACH stocks with
    L <= rho N_s."""
    import inspect

    ps = {(t, tb): 0.5 for tb in C.FAMILY_TICK_BINS for t in C.PRIMARY_TESTS}
    ps[("P-M", "Q-T1")], ps[("P-F", "Q-T3")] = 1.0, 0.0065
    out = P.combine(all_ev(), ps)
    assert len(out["family"]) == 8 and not out["holm"][("P-F", "Q-T3")]["rejected"]
    assert list(inspect.signature(P.evaluability).parameters) == ["Ns", "rho", "pm"]
    rho = np.full(6, 0.2)
    e = P.evaluability(np.array([640] * 5 + [100]), rho, PM)  # 640 / 5 = 128 = L
    assert e["P-M"]["evaluable"] and e["P-V"]["evaluable"] and e["P-F"]["evaluable"]
    e = P.evaluability(np.array([640] * 4 + [320, 100]), rho, PM)
    assert not e["P-M"]["evaluable"] and e["P-R"]["evaluable"] and "4 stocks reach L" in e["P-M"]["reason"]
    e = P.evaluability(np.array([100] * 6), rho, PM)
    assert not any(v["evaluable"] for v in e.values())


def check_passed_needs_conditions():
    """passed = Holm rejection AND every sub-condition."""
    mk = lambda n, p, ok: P.TestResult(n, p, "x", {"c1": True, "c2": ok}, {}, (8, 64), "", {}, {})  # noqa: E731
    res = {"P-V": mk("P-V", 0.001, False), "P-R": mk("P-R", 0.001, True), "P-M": mk("P-M", 0.9, True),
           "P-F": mk("P-F", 0.002, True)}
    P.set_decisions({"Q-T1": {"tests": res}}, P.holm(fam({k: v.p for k, v in res.items()})))
    assert res["P-V"].rejected and not res["P-V"].passed
    assert res["P-R"].passed and res["P-F"].passed and not res["P-M"].passed


# ---- P-F --------------------------------------------------------------------------------------------------------


def raw_from_K(slope_fn, S=8, D=8, seed=0) -> PO.Raw:
    """Raw moments with K_a(T) = T^(1 + slope_sd(a)) per stock-day, V = S = T (so v = 1), long days."""
    rng = np.random.default_rng(seed)
    base = 0.1 * rng.standard_normal((S, D))
    K = np.stack([T[None, None, :] ** (1 + base[:, :, None] + slope_fn(a)) for a in C.A_GRID], axis=2)
    V = np.broadcast_to(T, (S, D, NA, NT)).copy()
    return PO.Raw(V, K, np.broadcast_to(T, (S, D, NT)).copy(), np.full((S, D, NA), np.nan), np.full((S, D), 10**6))


def pf_reps(raw, rho=0.2, n_boot=60):
    S_, D_ = raw.Nsd.shape
    specs = [("cellceil", np.full(S_, rho))]
    dr = P.draws(S_, D_, np.random.default_rng(3), n_boot)
    return P.run_stats(lambda Wd, Ws: {"pf": P.snap(PO.core2(raw, Wd, Ws, specs, "tilde", **P.KW)["pf"][..., 0])},
                       dr)["pf"]


def check_pf_flat_not_rejected():
    """A k-tilde whose slope in a is exactly 0 (zeta_K the same at every a, stock-days differing) must not
    reject flatness: p = 1 on both ranges."""
    r = pf_reps(raw_from_K(lambda a: 0.0))
    for j in range(2):
        d = P.decide_pf(r[j])
        assert d["p"] == 1.0 and not d["conditions"]["interval_excludes_0"], d


def check_pf_two_sided():
    """zeta_K rising in a and falling in a are both rejected (two-sided)."""
    for sign in (+1, -1):
        r = pf_reps(raw_from_K(lambda a, s=sign: s * 0.05 * a + 0.01 * np.sin(7 * a)))
        d = P.decide_pf(r[1])
        assert d["p"] < 0.01 and d["conditions"]["interval_excludes_0"], (sign, d)


# ---- P-R -----------------------------------------------------------------------------------------------------------


def eq74(a, T_, s2=0.81, lam=0.1, A=1.0, B=0.0):
    """Maitrier and Bouchaud Eq. 74."""
    return np.exp(-s2 * a**2 / 2) * (A * np.exp(s2 * a / 2) + B * np.exp(lam * s2 * a * np.log(T_)))


def astar_reps(R_fn, sd=0.002, seed=0) -> P.Reps:
    """Reps of (a*(T1), a*(T2)) on [8, 64] from noisy copies of R_a(T) through primary3d.vertex."""
    rng = np.random.default_rng(seed)
    R0 = np.array([R_fn(C.A_GRID, t) for t in (8.0, 64.0)])  # (2, na)
    noisy = lambda n, s=sd: R0[None] + s * rng.standard_normal((n, *R0.shape))  # noqa: E731
    # leave-one-day-out values scatter sd / sqrt(D) around the estimate, as in reps()
    return P.Reps(P.vertex(R0), P.vertex(noisy(200)), P.vertex(noisy(8, sd / np.sqrt(8))), P.vertex(noisy(200)))


def check_vertex_matches_core2():
    """primary3d.vertex_ls (the a* before 1 Nov 2026) is pooled.core2's a* (same quadratic, same concavity rule)."""
    u = P.cell_unit({"Q-T3S1": g0_flow(6, 3, 800, 6)}, "Q-T3S1")
    o = PO.estimate2(u.raw, u.specs, **P.KW)
    assert np.allclose(P.vertex_ls(np.moveaxis(o["R"][..., 0], 0, -1)), o["a_star"][..., 0], equal_nan=True)


def check_pr_eq74_known_answer():
    """R_a(T) from Eq. 74 with B = 0 (a* = 1/2 at every T, the diagonal regime) passes the equivalence;
    with a* at 0.25 at both ends it fails."""
    d = P.decide_pr(astar_reps(lambda a, t: eq74(a, t, A=0.3 * t**-0.1)))
    assert np.allclose(P.vertex(np.array([eq74(C.A_GRID, 8.0)])), 0.5)
    assert d["p"] < 0.01 and all(d["conditions"].values()), d
    d = P.decide_pr(astar_reps(lambda a, t: 0.3 * np.exp(-0.81 * (a - 0.25) ** 2 / 2)))
    assert d["p"] > 0.5 and not any(d["conditions"].values()), d


def check_pr_uses_both_ends():
    """a* inside the margin at one end only fails, whichever end it is."""
    for good_T in (8.0, 64.0):
        R = lambda a, t, g=good_T: np.exp(-0.81 * (a - (0.5 if t == g else 0.2)) ** 2 / 2)  # noqa: E731
        d = P.decide_pr(astar_reps(R))
        assert d["p"] > 0.5 and sum(d["conditions"].values()) == 1, (good_T, d)


# ---- P-V -----------------------------------------------------------------------------------------------------------


def check_m3_is_shuffle_expectation():
    """M3 equals the mean of V_a(T) over within-day size shuffles (persistent signs, so the sign ACF term
    matters), within 3 % at every a and T."""
    f = g0_flow(1, 1, 3000, 8)
    m3 = P.m3_V(f)[0]
    rng = np.random.default_rng(0)
    acc = np.zeros_like(m3)
    n = 300
    for _ in range(n):
        acc += MO.day_moments(f.eps, rng.permutation(f.q), f.lm, f.ok, C.T_GRID, C.A_GRID).V
    assert np.max(np.abs(acc / n / m3 - 1)) < 0.03, np.max(np.abs(acc / n / m3 - 1))
    assert np.max(m3 / (T * m3[:, :1])) > 1.5  # the ACF term is far from negligible here


def check_m2_fit_recovers():
    """Per-a slopes generated exactly by M2 (mu_m 1.6, lambda sigma_l^2 0.1, a_c in (0.5, 0.75]) are recovered."""
    s = np.where(C.A_GRID < 0.6, (3 - 1.6 - 2 * C.A_GRID * 0.1) - 1, 0.0)
    m = P.fit_m2(s)
    assert m["k_below"] == 3 and np.isclose(m["mu_m"], 1.6) and np.isclose(m["lambda_sigma2"], 0.1)
    assert m["a_c"] == (0.5, 0.75) and m["rss"] < 1e-20


def check_loo_prefers_true_model():
    """Leave-one-day-out: with curves from M2 plus day noise, M2 beats M1 and a wrong M3; with curves from
    M3 itself, M3 beats M2."""
    r = (8, 64)
    x, _ = P.range_wx(r)
    zm1 = np.where(C.A_GRID < 0.6, 0.4 - 0.4 * C.A_GRID, 0.0)
    rng = np.random.default_rng(2)
    days = zm1[None, :, None] * x[None, None, :] + 0.01 * rng.standard_normal((8, NA, len(x)))
    days[:, :, 0] = 0.0
    Dtr = np.stack([days[np.arange(8) != d].mean(axis=0) for d in range(8)])
    e = P.model_errors(Dtr, days, np.zeros_like(days), r)
    assert e["M2"] < e["M1"] and e["M2"] < e["M3"] and e["days_used"] == 8, e
    e = P.model_errors(Dtr, days, days + 1e-4, r)
    assert e["M3"] < e["M2"], e


def check_pv_null_independent_sizes():
    """On flow whose sizes are independent of the signs, P-V (ii) does not reject: neither surrogate slope
    is significantly negative, so P-V fails."""
    f = {"Q-T3S1": g0_flow(10, 6, 1500, 11)}
    res = P.run_tick_bin(f, "Q-T3", PM, np.random.default_rng(1), n_boot=60, n_sur=2)
    t = res["tests"]["P-V"]
    st = t.statistic
    assert not (t.conditions["ii_S1_negative"] or t.conditions["ii_S2_negative"]), st
    assert min(st["p_S1"], st["p_S2"]) > 0.05, st


def dependent_flow(S, D, N, seed) -> P.Flow:
    """Half the runs carry a persistent sign and a small size, half a random sign and a large size: zeta_V(a)
    falls with a in the real flow, not in its surrogates."""
    rng = np.random.default_rng(seed)
    import synth as SY

    n = S * D * N
    pers = np.concatenate([SY.fgn_signs(N, 0.9, rng) for _ in range(S * D)])
    noise = rng.random(n) < 0.5
    eps = np.where(noise, np.where(rng.random(n) < 0.5, 1.0, -1.0), pers)
    q = np.exp(np.where(noise, 6.0, 3.0) + 0.3 * rng.standard_normal(n))
    lm = np.cumsum(0.01 * rng.standard_normal(n) + 0.002 * eps)
    return make_flow(eps, q, lm, np.full((S, D), N))


def check_pv_ii_detects_dependence():
    """Positive control: large runs with random signs make the real-minus-surrogate slope of zeta_V(a)
    negative for S1 and for S2."""
    res = P.run_tick_bin({"Q-T3S1": dependent_flow(10, 6, 1500, 12)}, "Q-T3", PM, np.random.default_rng(1),
                         n_boot=60, n_sur=2)
    t = res["tests"]["P-V"]
    assert t.conditions["ii_S1_negative"] and t.conditions["ii_S2_negative"], t.statistic
    assert max(t.statistic["p_S1"], t.statistic["p_S2"]) < 0.05


def check_pv_both_surrogates():
    """(ii) needs S1 AND S2: S1 clearly negative with S2 at 0 gives a large p and a failed S2 condition."""
    rng = np.random.default_rng(4)
    errs = {"M1": 2.0, "M2": 1.0, "M3": 3.0}
    d = P.decide_pv(errs, reps(np.array(-0.2), 0.02, rng), reps(np.array(0.0), 0.02, rng))
    assert d["p"] > 0.2 and d["conditions"]["ii_S1_negative"] and not d["conditions"]["ii_S2_negative"], d
    d = P.decide_pv(errs, reps(np.array(-0.2), 0.02, rng), reps(np.array(-0.2), 0.02, rng))
    assert d["p"] < 0.01 and all(d["conditions"].values()), d


def check_pv_gate():
    """(i) is required: M3 beating M2 (or M1 beating M2) gives p = 1 whatever (ii) says."""
    rng = np.random.default_rng(5)
    s = reps(np.array(-0.2), 0.02, rng)
    for errs in ({"M1": 2.0, "M2": 1.0, "M3": 0.5}, {"M1": 0.5, "M2": 1.0, "M3": 3.0}):
        d = P.decide_pv(errs, s, s)
        assert d["p"] == 1.0 and not (d["conditions"]["i_M2_beats_M1"] and d["conditions"]["i_M2_beats_M3"])


# ---- P-M -------------------------------------------------------------------------------------------------------------


def binned_flow(betas, S=6, D=4, N=4000, seed=0) -> P.Flow:
    """iid signs and log-normal sizes; run k moves the mid by G_b(l) = (1 + l)^-beta_b at lag l, b its size
    tercile within its stock-day (the estimator's bins), plus small noise."""
    rng = np.random.default_rng(seed)
    nb = len(betas)
    eps, q, lm = [], [], []
    lag = np.arange(N + 1, dtype=float)
    for _ in range(S * D):
        e = np.where(rng.random(N) < 0.5, 1.0, -1.0)
        qq = np.exp(3 + rng.standard_normal(N))
        b = np.searchsorted(np.quantile(qq, np.linspace(0, 1, nb + 1)[1:-1]), qq, side="left")
        m = np.zeros(N)
        for k, beta in enumerate(betas):
            g = np.concatenate([[0.0], (1 + lag[1:N]) ** -beta])
            m += np.convolve(np.where(b == k, e, 0.0), g)[:N]
        eps.append(e)
        q.append(qq)
        lm.append(m + 0.05 * rng.standard_normal(N).cumsum())
    return make_flow(np.concatenate(eps), np.concatenate(q), np.concatenate(lm), np.full((S, D), N))


def check_pm_recovers_beta():
    """Known kernels per size tercile (beta 0.6, 0.2, -0.3): beta-hat within 0.1 of each, the top bin's is the
    top tercile's, and P-M's two conditions hold."""
    truth = (0.6, 0.2, -0.3)
    res = P.run_tick_bin({"Q-T3S1": binned_flow(truth)}, "Q-T3", PM, np.random.default_rng(2), n_boot=40, n_sur=1)
    t = res["tests"]["P-M"]
    beta = np.array(t.statistic["beta_by_bin"]["estimate"])
    assert np.all(np.abs(beta - truth) < 0.1), beta
    assert t.statistic["top"]["estimate"] == beta[-1]
    assert all(t.conditions.values()) and t.p < 0.05, t.statistic


def check_pm_both_required():
    """Slope negative but top bin above 0 fails; top bin below 0 but slope positive fails; both pass."""
    rng = np.random.default_rng(6)
    neg, pos = reps(np.array(-0.3), 0.03, rng), reps(np.array(0.3), 0.03, rng)
    d = P.decide_pm(neg, pos)
    c = d["conditions"]
    assert d["p"] > 0.5 and c["slope_declines_with_log_q"] and not c["top_bin_upper_below_0"]
    d = P.decide_pm(pos, neg)
    assert d["p"] > 0.5 and not d["conditions"]["slope_declines_with_log_q"]
    d = P.decide_pm(neg, neg)
    assert d["p"] < 0.01 and all(d["conditions"].values())


def check_pm_params_and_min_runs():
    """P-M parameters must be among the allowed values with L >= 8 l_max; a stock-day with a bin below
    min_runs_per_bin, or a stock whose ceiling does not reach L, is not used."""
    for bad in ((6, 16, 128, 10), (3, 24, 192, 10), (3, 32, 128, 10), (3, 16, 128, 0)):
        with pytest.raises(ValueError):
            P.PMParams(*bad)
    f = binned_flow((0.3, 0.3, 0.3), S=2, D=2, N=1000)
    rho = np.array([0.2, 0.2])
    assert P.pm_moments(f, rho, P.PMParams(3, 16, 128, 300))["used"].all()
    assert not P.pm_moments(f, rho, P.PMParams(3, 16, 128, 400))["used"].any()  # 333 runs per tercile
    assert not P.pm_moments(f, np.array([0.1, 0.2]), P.PMParams(3, 16, 128, 300))["used"][:2].any()  # 128 > 100


# ---- the whole unit -------------------------------------------------------------------------------------------------


def check_run_primary():
    """End to end: Q-T3 (three cells, long days) has every test evaluable; Q-T2 (Q-T2S3 alone, 200-run days)
    reaches [4, 32] only and not L = 128, so P-M of Q-T2 is not computed; Q-T2 is descriptive (2 Nov 2026):
    its P-V, P-R and P-F are computed and never enter Holm (rejected and passed None). Q-T1 has no cell in the
    input and leaves whole. Holm sees exactly the 4 pairs of Q-T3; P-K is descriptive."""
    flows = {"Q-T3S1": g0_flow(6, 4, 1500, 21), "Q-T3S2": g0_flow(4, 4, 1500, 22), "Q-T3S3": g0_flow(3, 4, 1500, 23),
             "Q-T2S3": g0_flow(6, 4, 200, 24)}
    res = P.run_primary(flows, PM, np.random.default_rng(3), n_boot=30, n_sur=1)
    fam_ = {(t, "Q-T3") for t in C.PRIMARY_TESTS}
    assert set(res["family"]) == fam_ and set(res["holm"]) == fam_
    left = {(r["test"], r["tick_bin"]): r["reason"] for r in res["left"]}
    assert set(left) == {(t, "Q-T1") for t in C.PRIMARY_TESTS} and "no confirmatory cell" in left[("P-V", "Q-T1")]
    desc = {(d["test"], d["tick_bin"]): d for d in res["descriptive"]}
    assert set(desc) == {(t, "Q-T2") for t in C.PRIMARY_TESTS} and "stocks reach L" in desc[("P-M", "Q-T2")]["reason"]
    assert "P-M" not in res["by_tick_bin"]["Q-T2"]["tests"] and "Q-T1" not in res["by_tick_bin"]
    assert res["by_tick_bin"]["Q-T2"]["tests"]["P-F"].range_used == (4, 32)
    for tb, tests in (("Q-T3", C.PRIMARY_TESTS), ("Q-T2", ("P-V", "P-R", "P-F"))):
        assert "P-K" in res["by_tick_bin"][tb]["descriptive"]
        for k in tests:
            t = res["by_tick_bin"][tb]["tests"][k]
            assert 0.0 <= t.p <= 1.0
            if tb == "Q-T2":
                assert t.rejected is None and t.passed is None
            else:
                assert t.rejected is not None and t.passed == bool(t.rejected and t.all_conditions)
    t = res["by_tick_bin"]["Q-T3"]["tests"]["P-R"]
    assert set(t.by_range) == {"[4,32]", "[8,64]"} and t.range_used == (8, 64)
    assert {"Rratio_minus_exp_sig2_over_8", "dastar_dlogT", "sigma_l2", "psi_hat"} <= set(t.secondary)


def psi_flow(psi, S=8, D=4, N=4000, seed=0) -> P.Flow:
    """iid signs, permanent impact eps q^psi: E[Delta_1 eps | q] = q^psi exactly (plus zero-mean noise)."""
    rng = np.random.default_rng(seed)
    n = S * D * N
    eps = np.where(rng.random(n) < 0.5, 1.0, -1.0)
    q = np.exp(3 + rng.standard_normal(n))
    x = eps * q**psi + 0.5 * rng.standard_normal(n)
    lm = np.concatenate([[0.0], np.cumsum(x)[:-1]])  # the mid read just before run n
    return make_flow(eps, q, lm, np.full((S, D), N))


def check_psi_known_answer():
    """psi-hat (per stock OLS slope of log E[Delta_1 eps | q-bin] on log median q of the bin, mean over
    stocks) recovers psi = 0.3 and 0.6 within 0.03 when E[Delta_1 eps | q] = q^psi."""
    for psi in (0.3, 0.6):
        f = psi_flow(psi)
        mo = P.psi_moments(f)
        S_ = f.n_stocks
        est, n = P.psi_hat(mo, np.ones((1, f.n_days)), np.ones((1, S_)), f.n_days)
        assert abs(est[0] - psi) < 0.03 and n[0] == S_, (psi, est, n)


CHECKS = {f.__name__: f for f in [
    check_config_table, check_descriptive_excluded, check_cell_ceilings, check_r2_clip, check_fallback_by_reach,
    check_fallback_never_by_result, check_interval_matches_intervals3, check_p_value_conservative,
    check_holm_order, check_holm_family_keys, check_one_family_across_bins, check_qt2_never_in_holm,
    check_nonevaluable_leaves_family,
    check_evaluability_reach_only, check_passed_needs_conditions, check_pf_flat_not_rejected,
    check_pf_two_sided, check_vertex_matches_core2, check_pr_eq74_known_answer, check_pr_uses_both_ends,
    check_m3_is_shuffle_expectation, check_m2_fit_recovers, check_loo_prefers_true_model,
    check_pv_null_independent_sizes, check_pv_ii_detects_dependence, check_pv_both_surrogates, check_pv_gate,
    check_pm_recovers_beta, check_pm_both_required, check_pm_params_and_min_runs, check_run_primary,
    check_psi_known_answer,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_3d1(name):
    CHECKS[name]()


class Bug(NamedTuple):
    id: str
    module: object
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("Holm sorted descending", P, "step_down", 'order = np.argsort(p, kind="stable")',
        'order = np.argsort(-p, kind="stable")', ("check_holm_order",)),
    Bug("Holm as plain Bonferroni", P, "step_down", "thr = alpha / (m - i)", "thr = alpha / m", ("check_holm_order",)),
    Bug("Holm not stopping at the first failure", P, "step_down", "stop = True", "stop = False", ("check_holm_order",)),
    Bug("Holm accepting a member that is not a (primary test, tick bin) pair", P, "holm",
        "if not (isinstance(k, tuple) and len(k) == 2 and k[0] in C.PRIMARY_TESTS and k[1] in C.FAMILY_TICK_BINS):",
        "if False:", ("check_holm_family_keys", "check_qt2_never_in_holm")),
    Bug("Holm accepting a Q-T2 pair (family over every tick bin)", P, "holm", "k[1] in C.FAMILY_TICK_BINS):",
        "k[1] in C.TICK_BINS):", ("check_holm_family_keys", "check_qt2_never_in_holm")),
    Bug("combine making Q-T2 pairs members", P, "combine",
        "pairs = [(t, tb) for tb in ev if tb in C.FAMILY_TICK_BINS for t in C.PRIMARY_TESTS]",
        "pairs = [(t, tb) for tb in ev for t in C.PRIMARY_TESTS]", ("check_qt2_never_in_holm",)),
    Bug("run_primary passing Q-T2's p-values to the family", P, "run_primary",
        "for tb in out if tb in C.FAMILY_TICK_BINS for t in out[tb]", "for tb in out for t in out[tb]",
        ("check_run_primary",)),
    Bug("Holm per tick bin instead of one family", P, "combine", "h = holm({k: pvals[k] for k in members}, alpha)",
        "h = {k: v for b in {x[1] for x in members} for k, v in "
        "holm({k: pvals[k] for k in members if k[1] == b}, alpha).items()}",
        ("check_one_family_across_bins",)),
    Bug("non-evaluable pair kept in the family with p = 1", P, "combine",
        "h = holm({k: pvals[k] for k in members}, alpha)",
        "h = holm({**{k: 1.0 for k in pairs}, **{k: pvals[k] for k in members}}, alpha)",
        ("check_nonevaluable_leaves_family",)),
    Bug("evaluability decided from the computed p-value", P, "combine",
        'members = [k for k in pairs if ev[k[1]][k[0]]["evaluable"]]',
        "members = [k for k in pairs if k in pvals and pvals[k] < 1]",
        ("check_evaluability_reach_only",)),
    Bug("every test computed before evaluability is applied", P, "run_primary",
        'tests=tuple(t for t in C.PRIMARY_TESTS if ev[tb][t]["evaluable"])', "tests=C.PRIMARY_TESTS",
        ("check_run_primary",)),
    Bug("fallback range picked by the smaller p", P, "_assemble", "d = by_r[r_used]",
        'd = min(by_r.values(), key=lambda v: v["p"])', ("check_fallback_never_by_result",)),
    Bug("[8, 64] used whatever its reach", P, "select_range", "if by[C.PRIMARY_RANGE] >= PO.MIN_REACH:",
        "if True:", ("check_fallback_by_reach",)),
    Bug("descriptive cells pooled into the tick bin", P, "tick_bin_unit",
        "C.tick_bin(c) == tb and c in C.CONFIRMATORY", "C.tick_bin(c) == tb", ("check_descriptive_excluded",)),
    Bug("one ceiling for every stock of the tick bin", P, "make_unit",
        "np.full(flows[c].n_stocks, C.CONFIRMATORY[c])", "np.full(flows[c].n_stocks, max(C.CONFIRMATORY.values()))",
        ("check_cell_ceilings",)),
    Bug("R2 clip skipped", P, "make_unit", "r2_clip(fl.q, fl.lengths, fl.n_days)", "fl.q", ("check_r2_clip",)),
    Bug("p-value from the least conservative method", P, "p_value", "p = np.nanmax(ps, axis=0)",
        "p = np.nanmin(ps, axis=0)", ("check_p_value_conservative",)),
    Bug("P-F tested one-sided", P, "decide_pf", 'p_value(pf, 0.0, "two-sided")', 'p_value(pf, 0.0, "less")',
        ("check_pf_two_sided",)),
    Bug("P-F slope not snapped at numerical zero", P, "snap", "np.abs(x) < SNAP", "np.abs(x) < 0",
        ("check_pf_flat_not_rejected",)),
    Bug("P-R one-sided on drift instead of the equivalence", P, "decide_pr", '"p": float(np.max(p_end))',
        '"p": float(p_value(Reps(ends.est[1] - ends.est[0], ends.ds[:, 1] - ends.ds[:, 0], '
        'ends.jk[:, 1] - ends.jk[:, 0], ends.ws[:, 1] - ends.ws[:, 0]), 0.0, "greater"))',
        ("check_pr_eq74_known_answer",)),
    Bug("P-R at T2 only", P, "decide_pr", '"p": float(np.max(p_end))', '"p": float(p_end[1])',
        ("check_pr_uses_both_ends",)),
    Bug("M3 missing the factor 2 of the sign ACF term", P, "m3_V", "2 * np.mean(qa) ** 2 * tri",
        "np.mean(qa) ** 2 * tri", ("check_m3_is_shuffle_expectation",)),
    Bug("P-V (ii) on S1 only", P, "decide_pv", '"p": max(p1, p2) if gate else 1.0', '"p": p1 if gate else 1.0',
        ("check_pv_both_surrogates",)),
    Bug("P-V without the model-comparison gate", P, "decide_pv", '"p": max(p1, p2) if gate else 1.0',
        '"p": max(p1, p2)', ("check_pv_gate",)),
    Bug("P-M on the slope condition only", P, "decide_pm", '"p": max(p_s, p_t)', '"p": p_s',
        ("check_pm_both_required",)),
    Bug("P-M top bin taken as the smallest sizes", P, "pm_stats", '"top": beta[:, -1]', '"top": beta[:, 0]',
        ("check_pm_recovers_beta",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(bug.module, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()


def test_planted_qt2_in_holm_both_layers_is_caught(mutate):
    """Q-T2 entering Holm with nothing to stop it: holm's family check AND combine's tick-bin filter widened
    together, so no ValueError is raised and Q-T2's p-values reach the step-down."""
    mutate(P, "holm", "k[1] in C.FAMILY_TICK_BINS):", "k[1] in C.TICK_BINS):")
    mutate(P, "combine", "pairs = [(t, tb) for tb in ev if tb in C.FAMILY_TICK_BINS for t in C.PRIMARY_TESTS]",
           "pairs = [(t, tb) for tb in ev for t in C.PRIMARY_TESTS]")
    with pytest.raises(Exception):  # noqa: B017
        check_qt2_never_in_holm()
