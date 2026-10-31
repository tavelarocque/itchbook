"""Tests for Phase 3c-4: independent stationary days (panel3c), the whole-stock bootstrap and the widest of
three intervals (pooled), the pass rule with the Monte Carlo band and the all-ceilings-pass rho* rule
(phase3c4); and the rules decided 31 Oct 2026 from docs/phase3c4_diagnosis.md: range reach (pooled.range_rule),
the 90 % evaluability threshold, the P-F slope by coverage only, and the tick-bin rho*. No data files are read.

    uv run --group dev --group analysis pytest tools/est/test_3c3.py tools/est/test_3c4.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import gen3 as G
import panel3c as PC
import phase3c4 as C4
import pooled as PO

STRONG = {"w1": 0.6, "H": 0.99, "w2": 0.0, "phi": 0.5}
CELL = G.Cell("TEST", np.arange(100, 200), 3.0, 1.0, (0.3, 0.7), STRONG)
T = PO.T_GRID2.astype(float)
NA, NT = len(PO.A_GRID2), len(T)


def check_days_independent():
    """New path: a stock's days are separate series, so with a strong slowly varying sign component the mean
    signs of its day 1 and day 2 are uncorrelated across 300 stocks (|r| < 0.25); the continuous path
    correlates them (> 0.5)."""
    Nsd = np.full((300, 2), 400)
    for cont, lo, hi in ((False, -0.25, 0.25), (True, 0.5, 1.0)):
        p = PC.make_panel("G0", CELL, Nsd, np.random.default_rng(2), params={"beta": 0.4}, burn=2000, clip=None,
                          continuous=cont)
        assert len(p.eps) == Nsd.sum()
        m = p.eps.reshape(300, 2, 400).mean(axis=2)
        r = np.corrcoef(m[:, 0], m[:, 1])[0, 1]
        assert lo < r < hi, (cont, r)


def check_burn_per_day():
    """Every stock-day gets its own burn-in: the series lengths are burn + N for every day of every stock."""
    Nsd = np.array([[50, 70, 20], [30, 90, 10]])
    s = PC.independent_series(Nsd, 500)
    assert s.shape == (6, 1) and np.array_equal(s.ravel(), 500 + Nsd.ravel())


def raw_from(slopes):
    """Raw moments with K = T^(1 + slope_sd) per stock-day (slopes (S, D)), V = T, S = T."""
    slopes = np.asarray(slopes, float)
    S_, D_ = slopes.shape
    K = T[None, None, :] ** (1 + slopes[:, :, None])
    rep = lambda x: np.repeat(x[:, :, None, :], NA, axis=2)  # noqa: E731
    V = np.broadcast_to(T, (S_, D_, NT)).copy()
    return PO.Raw(rep(V), rep(K), V.copy(), np.full((S_, D_, NA), np.nan), np.full((S_, D_), 10**6))


def check_whole_stock_bootstrap():
    """Whole-stock bootstrap: stock multiplicities sum to S and every day keeps weight 1; its replicates vary
    with stock differences only, and are constant with day differences only."""
    Wd, Ws = PO.whole_stock_weights(6, 8, 50, np.random.default_rng(3))
    assert np.all(Wd == 1) and np.allclose(Ws.sum(axis=1), 6) and Ws.std() > 0
    rng = np.random.default_rng(4)
    stock_only = raw_from(np.repeat(0.2 + 0.1 * rng.standard_normal((6, 1)), 8, axis=1))
    day_only = raw_from(np.repeat(0.2 + 0.1 * rng.standard_normal((1, 8)), 6, axis=0))
    for raw, varies in ((stock_only, True), (day_only, False)):
        Wd, Ws = PO.whole_stock_weights(6, 8, 40, np.random.default_rng(5))
        z = PO.core2(raw, Wd, Ws, [("all",)])["zeta_K"][:, 0, 0, 0]
        assert (z.std() > 1e-3) == varies, (varies, z.std())


def check_widest_of_three():
    """intervals3 reports, per element, the widest of its three candidates."""
    rng = np.random.default_rng(6)
    raw = raw_from(0.2 + 0.05 * rng.standard_normal((6, 8)))
    iv = PO.intervals3(raw, [("all",)], np.random.default_rng(7), n_boot=40, keys=("zeta_K",))
    e, lo, hi, meth = iv["zeta_K"]
    assert np.all(np.isin(meth, (1, 2, 3))) and np.all(hi > lo)
    se, n = PO.jackknife2(raw, [("all",)], "day")
    jlo, jhi = PO.jk_interval(e, se["zeta_K"], n, PO.LEVEL)
    assert np.all(hi - lo >= jhi - jlo - 1e-12)


def check_pass_rule_mc_band():
    """Coverage fails only below 0.85 by more than 2 MC SE; a bin within 2 MC SE of 0.85 is in the band."""
    assert C4.passes(0.0, 0.83, 0.02) and C4.passes(0.0, 0.90, 0.02)
    assert not C4.passes(0.0, 0.80, 0.02) and not C4.passes(0.03, 0.95, 0.01)
    assert C4.in_mc_band(0.83, 0.02) and C4.in_mc_band(0.88, 0.02) and not C4.in_mc_band(0.95, 0.02)


def check_rho_star_all_pass():
    """rho* = the largest ceiling with every evaluable ceiling up to it passing; not-evaluable ones are skipped."""
    assert C4.rho_star({1 / 80: True, 1 / 40: False, 1 / 20: True, 1 / 10: True}) == 1 / 80
    assert C4.rho_star({1 / 80: None, 1 / 40: True, 1 / 20: True, 1 / 10: False}) == 1 / 20
    assert np.isnan(C4.rho_star({1 / 80: False, 1 / 40: True}))


# ---- the 31 Oct 2026 rules (docs/phase3c4_diagnosis.md) ----------------------------------------------

CEIL20 = [("ceil", 1 / 20)]


def raw_ns(Ns, c=None, D=2):
    """Raw moments with K = c_s T^1.25, V = S = T for stocks with fewest-day lengths Ns (eligible at T under
    the 1/20 ceiling iff T <= Ns / 20)."""
    Ns = np.asarray(Ns)
    c = np.ones(len(Ns)) if c is None else np.asarray(c, float)
    K = c[:, None, None] * np.broadcast_to(T**1.25, (len(Ns), D, NT))
    rep = lambda x: np.repeat(x[:, :, None, :], NA, axis=2)  # noqa: E731
    V = np.broadcast_to(T, (len(Ns), D, NT)).copy()
    return PO.Raw(rep(V), rep(K), V.copy(), np.full((len(Ns), D, NA), np.nan), np.repeat(Ns[:, None], D, axis=1))


def check_reach_at_top():
    """At least 5 stocks must reach T2 under the ceiling. 5 stocks reach T = 32 (N = 640) and 2 more only T = 4:
    [4, 32] is defined (zeta_K = 1.25, zeta_V = 1) and [8, 64] is not. 4 stocks reach 32 and 3 more reach 8:
    [4, 32] is undefined for every estimand, although 7 stocks reach T1 = 4 (and the old estimator, with every
    horizon populated, is finite there)."""
    r = raw_ns([640] * 5 + [100] * 2)
    assert np.array_equal(PO.range_reach(r.Nsd.min(axis=1), CEIL20), [[5], [0], [0]])
    e = PO.estimate2(r, CEIL20)
    assert np.allclose(e["zeta_K"][:, 0, 0], 1.25) and np.allclose(e["zeta_V"][:, 0, 0], 1.0)
    assert np.isfinite(e["pf"][0, 0]) and np.isfinite(e["Rratio"][0, 0])
    for key in PO.ESTIMANDS2:
        assert np.all(np.isnan(e[key][..., 1:, 0])), key
    r = raw_ns([640] * 4 + [200] * 3)
    e = PO.estimate2(r, CEIL20)
    for key in PO.ESTIMANDS2:
        assert np.all(np.isnan(e[key][..., 0, 0])), key
    assert np.all(np.isfinite(PO.estimate2(r, CEIL20, **PO.LEGACY_RULE)["zeta_K"][:, 0, 0]))


def drop_case():
    """6 stocks, every horizon reached; stock 0 has k~ = -100 at T = 8, so the pooled k~ there is negative."""
    r = raw_ns([10**6] * 6, c=[1, 3, 3, 3, 3, 3])
    r.K[0, :, :, 3] = -100 * 8.0
    return r


def check_full_range_only():
    """zeta_K is never a slope over a shortened range: with T = 8 dropped (non-positive pooled k~) zeta_K and the
    P-F slope are undefined on [4, 32] and [8, 64] (the old estimator joins the other horizons) and defined on
    [16, 256]; zeta_V is defined on every range. With no stock at T = 32 (6 stocks reaching T = 16), [4, 32] is
    undefined although the old estimator fits it from T = 4, 8, 16."""
    e = PO.estimate2(drop_case(), [("all",)])
    assert e["dropped_K"][0, 3, 0]
    assert np.all(np.isnan(e["zeta_K"][:, :2, 0])) and np.all(np.isnan(e["pf"][:2, 0]))
    assert np.allclose(e["zeta_K"][:, 2, 0], 1.25) and np.isfinite(e["pf"][2, 0])
    assert np.allclose(e["zeta_V"][:, :, 0], 1.0)
    assert np.all(np.isfinite(PO.estimate2(drop_case(), [("all",)], **PO.LEGACY_RULE)["zeta_K"][:, 0, 0]))
    r = raw_ns([400] * 6)
    assert np.all(np.isnan(PO.estimate2(r, CEIL20)["zeta_K"][:, 0, 0]))
    assert np.all(np.isfinite(PO.estimate2(r, CEIL20, **PO.LEGACY_RULE)["zeta_K"][:, 0, 0]))


def check_saved_panels_rule():
    """phase3c4 applies the rule to saved (old-estimator) panels from their reach and dropped_K: masking the old
    estimates reproduces the new estimator exactly."""
    cases = [(raw_ns([640] * 5 + [100] * 2), CEIL20), (raw_ns([640] * 4 + [200] * 3), CEIL20),
             (raw_ns([400] * 6), CEIL20), (drop_case(), [("all",)])]
    for r, sp in cases:
        old, new = PO.estimate2(r, sp, **PO.LEGACY_RULE), PO.estimate2(r, sp)
        kept = C4.kept_full_from_dropped(old["dropped_K"], T_grid=PO.T_GRID2, ranges=PO.RANGES2)
        rule = PO.range_rule(PO.range_reach(r.Nsd.min(axis=1), sp), kept)
        for key in PO.ESTIMANDS2:
            assert np.allclose(np.where(rule[key], old[key], np.nan), new[key], equal_nan=True), key


def saved(P, defined, est=0.05, half=1.0):
    """A saved-study dict for pf and Rratio (truth 0): estimate est, interval est -+ half, defined in the first
    `defined` of P panels."""
    st, dm = {}, {}
    for key in ("pf", "Rratio"):
        e = np.full((P, 2, 5), est, dtype=np.float32)
        st[f"{key}_est"], st[f"{key}_lo"], st[f"{key}_hi"] = e, e - half, e + half
        dm[key] = np.arange(P)[:, None, None] < defined
    return st, {"pf": (np.zeros(2), None), "Rratio": (np.zeros(2), None)}, dm


def check_pf_coverage_only():
    """Bias 0.05 with full coverage: the P-F slope passes (coverage only, its bias still reported), R_{1/2}/R_0
    fails on bias."""
    st, tr, dm = saved(40, 40)
    o = C4.evaluate(st, tr, 0, ("pf", "Rratio"), dm)
    assert o["pf"]["pass"] and np.isclose(o["pf"]["max_abs_bias"], 0.05)
    assert not o["Rratio"]["pass"]
    assert C4.passes(0.05, 0.9, 0.01, coverage_only=True) and not C4.passes(0.05, 0.9, 0.01)


def check_evaluable_90():
    """An element is evaluable only if defined in at least 90 % of panels (36 of 40, not 35); coverage is over
    the defined panels (all of them cover: 1.0, not 36 / 40)."""
    st, tr, dm = saved(40, 36, est=0.0)
    o = C4.evaluate(st, tr, 0, ("pf",), dm)["pf"]
    assert o["evaluable"] and o["min_cov"] == 1.0 and np.isclose(o["min_defined"], 0.9)
    st, tr, dm = saved(40, 35, est=0.0)
    assert not C4.evaluate(st, tr, 0, ("pf",), dm)["pf"]["evaluable"]


def check_tick_bin_rho():
    """A tick bin's rho* is the minimum of its cells' rho*, and none if any cell has none."""
    assert C4.tick_bin_rho({"a": 1 / 5, "b": 1 / 20, "c": 1 / 10}) == 1 / 20
    assert np.isnan(C4.tick_bin_rho({"a": 1 / 5, "b": np.nan, "c": 1 / 10}))


CHECKS = {f.__name__: f for f in [
    check_days_independent, check_burn_per_day, check_whole_stock_bootstrap, check_widest_of_three,
    check_pass_rule_mc_band, check_rho_star_all_pass, check_reach_at_top, check_full_range_only,
    check_saved_panels_rule, check_pf_coverage_only, check_evaluable_90, check_tick_bin_rho,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_3c4(name):
    CHECKS[name]()


class Bug(NamedTuple):
    id: str
    module: object
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("days of a stock drawn from one series on the new path", PC, "make_panel",
        "return _independent_panel(gen, cell, Nsd, rng, sign_model, sm, clip, params, int(burn))",
        "return _stationary_panel(gen, cell, Nsd, rng, sign_model, sm, clip, params, int(burn))",
        ("check_days_independent",)),
    Bug("burn-in shared across a stock's days", PC, "independent_series",
        "return (burn + Nsd.reshape(-1, 1)).astype(np.int64)",
        "return (Nsd.reshape(-1, 1) + burn * (np.arange(Nsd.size) % Nsd.shape[1] == 0)[:, None]).astype(np.int64)",
        ("check_burn_per_day",)),
    Bug("whole-stock bootstrap resampling days instead of stocks", PO, "whole_stock_weights",
        "return np.ones((B, n_d)), Ws",
        "return np.stack([np.bincount(rng.integers(0, n_d, n_d), minlength=n_d) for _ in range(B)]).astype(float), "
        "np.ones((B, n_s))",
        ("check_whole_stock_bootstrap",)),
    Bug("rho* taken as the largest passing ceiling", C4, "rho_star",
        "        if not status[e]:\n            break", "        if not status[e]:\n            continue",
        ("check_rho_star_all_pass",)),
    Bug("Monte Carlo band ignored", C4, "passes",
        "~(np.asarray(cov) + 2 * np.asarray(mc_se) < COV_MIN)", "~(np.asarray(cov) < COV_MIN)",
        ("check_pass_rule_mc_band",)),
    Bug("slope computed from fewer horizons than the range", PO, "range_rule",
        "zK = zV & np.asarray(kept_full) if full_range else zV", "zK = zV",
        ("check_full_range_only",)),
    Bug("reach threshold applied to the lowest T instead of T2", PO, "range_reach",
        "np.array([hi for _, hi in ranges])", "np.array([lo for lo, _ in ranges])",
        ("check_reach_at_top",)),
    Bug("P-F slope judged on bias too", C4, "evaluate",
        "coverage_only=key in COVERAGE_ONLY", "coverage_only=False",
        ("check_pf_coverage_only",)),
    Bug("evaluable at 50 % of panels", C4, "evaluate",
        "ev = dfrac >= MIN_EVAL", "ev = dfrac >= 0.5",
        ("check_evaluable_90",)),
    Bug("coverage over all panels instead of the defined ones", C4, "evaluate",
        "axis=0).ravel() / nd)", "axis=0).ravel() / len(e))",
        ("check_evaluable_90",)),
    Bug("tick-bin rho* ignoring a cell with none", C4, "tick_bin_rho",
        "return float(np.min(v)) if np.all(np.isfinite(v)) else np.nan", "return float(np.nanmin(v))",
        ("check_tick_bin_rho",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(bug.module, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
