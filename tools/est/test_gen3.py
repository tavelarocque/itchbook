"""Tests for the Phase 3a calibration (calib3.py), generators (gen3.py) and
ground-truth runner (truth3.py), on synthetic cells built in memory (no data
files are read, except the temporary calibration files written by one check).

1. Calibration pieces: sign ACF pooling never crosses a stock-day boundary,
   power-law fit, fGn-sign calibration reproduces a target C(l).
2. Generators: per-bin FFT convolution against the direct sum, GT with c = 0
   is G0 (identical arrays), G0's ground-truth zeta_K(a) is non-increasing,
   GF's is 1 for every a, day lengths come from the right cell.
3. Ground truth: the runner stops only when the SE target is met, and resumes.
4. Planted bugs (every one must be caught).

    uv run --group dev --group analysis pytest tools/est/test_gen3.py
"""
from __future__ import annotations

import json
from typing import NamedTuple

import numpy as np
import pytest

import calib3 as C3
import gen3 as G
import truth3 as TR

SM = {"w1": 0.1, "H": 0.85, "w2": 0.6, "phi": 0.5}


def cell(sigma: float = 1.2, Ns=None, name: str = "TEST") -> G.Cell:
    return G.Cell(name, np.arange(100, 200) if Ns is None else np.asarray(Ns), 3.0, sigma, (0.3, 0.7), SM)


def quiet(*_):
    pass


# ---- calibration ------------------------------------------------------------------


def check_sign_acf_pooling_within_day():
    """Two stock-days of constant opposite signs: every within-day pair agrees, so C(l) = 1 at every
    lag; pooling across the boundary would mix in disagreeing pairs. Then random series against the
    direct sum over pairs of each series."""
    lags = np.array([1, 5, 10, 20, 39])
    C, num, cnt = C3.pooled_acf([np.ones(40), -np.ones(40)], lags)
    assert np.allclose(C, 1.0) and np.allclose(cnt.sum(0), 2 * (40 - lags))
    rng = np.random.default_rng(3)
    xs = [rng.choice([-1.0, 1.0], n) for n in (7, 30, 55)]
    C, _, _ = C3.pooled_acf(xs, lags)
    want = [sum(float(x[:-l] @ x[l:]) for x in xs if len(x) > l) / sum(max(len(x) - l, 0) for x in xs) for l in lags]
    assert np.allclose(C, want)


def check_power_fit():
    l = C3.sign_lags(500)
    assert {10, 100} <= set(l.tolist()) and l.max() == 500
    c, g = C3.power_fit(l[l >= 10], 0.2 * l[l >= 10] ** -0.4)
    assert np.isclose(c, 0.2) and np.isclose(g, 0.4)


def check_sign_calibration():
    """fit_sign_model on an exact target C(l) (lags 1..500), then 8 x 2^16 signs from the fit: their
    pooled C(l) within 0.01 of the target at lags 1..100 (measured gap below 0.004 at lags <= 10,
    0.003 beyond)."""
    lags = C3.sign_lags(500)
    target = G.sign_model_C(lags, **SM)
    fit = G.fit_sign_model(lags, target)
    rng = np.random.default_rng(5)
    Cm, _, _ = C3.pooled_acf([G.signs(1 << 16, fit, rng).astype(float) for _ in range(8)], lags)
    pick = lags <= 100
    assert np.max(np.abs(Cm[pick] - target[pick])) < 0.01, np.round(Cm[pick] - target[pick], 4)


# ---- generators -------------------------------------------------------------------


def check_binned_fft_equals_direct():
    rng = np.random.default_rng(6)
    n = 300
    x = rng.standard_normal(n)
    b = rng.integers(0, G.N_BINS, n)
    for betas in (rng.uniform(0, 1, G.N_BINS), np.repeat([0.2, 0.6, 0.9], G.N_BINS // 3), np.full(G.N_BINS, 0.4)):
        assert np.allclose(G.binned_impact(x, b, betas), G.direct_binned_impact(x, b, betas), atol=1e-10)


def check_gt_c0_is_g0():
    """Same seed, c = 0: GT and G0 give identical arrays (bit for bit), for both size models."""
    c = cell()
    for copula in (False, True):
        a = G.simulate_g0(c, 5000, np.random.default_rng(7), beta=0.45, psi=0.5, copula=copula)
        b = G.simulate_gt(c, 5000, np.random.default_rng(7), beta_max=0.45, c=0.0, l_star=4.0, psi=0.5,
                          copula=copula)
        for f in ("eps", "q", "lm", "ok"):
            assert np.array_equal(getattr(a, f), getattr(b, f)), f


def check_g0_monotone():
    """G0's ground-truth zeta_K(a) on [8, 64] falls at every step of the grid up to a = 1.5 (measured
    steps -0.014 to -0.016; 12 days of 2 x 10^5 runs, sigma_l = 1.2). Larger a is left out: its SE
    reaches 0.1 at this size."""
    r = TR.run_truth("G0", cell(1.2), {"beta": 0.4}, N=200_000, min_days=12, max_days=12, chunk=12, workers=1,
                     out_dir=None, n_collapse=0, log=quiet)
    z = r.curves["zK_8_64"][TR.A_GRID <= 1.5]
    assert np.all(np.diff(z) < 0), np.round(np.diff(z), 4)


def check_gf_zeta_one():
    """GF: zeta_K(a) = 1 within 0.02 for every a in A_GRID on [8, 64] and [16, 256] (8 days of
    2 x 10^5 runs, sigma_l = 0.4, informed share 0.3, no noise; measured 0.987 to 1.006)."""
    r = TR.run_truth("GF", cell(0.4), {"sign_model": SM, "f_inf": 0.3, "noise": 0.0}, N=200_000, min_days=8,
                     max_days=8, chunk=8, workers=1, out_dir=None, n_collapse=0, log=quiet)
    for k in ("zK_8_64", "zK_16_256"):
        assert np.all(np.abs(r.curves[k] - 1) < 0.02), (k, np.round(r.curves[k], 3))


def _write_calibration(path, name, Ns):
    row = {"mu_logq": "3", "sigma_l": "1", "copula_w": "0.3", "copula_H": "0.7"} | {f"sgn_{k}": str(v)
                                                                                     for k, v in SM.items()}
    np.savez(path / f"calibration_{name}.npz", N=np.asarray(Ns), row=json.dumps(row))


def check_day_lengths_from_cell(tmp):
    """Two calibrated cells with disjoint N supports: every day of a cell's panel has a length from that
    cell's stock-day N values, and both halves of the support turn up."""
    _write_calibration(tmp, "A", np.arange(100, 200))
    _write_calibration(tmp, "B", np.arange(5000, 5100))
    saved = G.phase3_dir
    G.phase3_dir = lambda: tmp
    G.load_cell.cache_clear()
    try:
        for name, lo, hi in (("B", 5000, 5100), ("A", 100, 200)):
            p = G.panel("G0", name, 3, 6, seed=9)
            n = np.array([len(d.eps) for row in p for d in row])
            assert np.all((n >= lo) & (n < hi)), (name, n)
            assert n.min() < (lo + hi) / 2 <= n.max()
            assert all(len(d.q) == len(d.lm) == len(d.ok) == len(d.eps) for row in p for d in row)
    finally:
        G.phase3_dir = saved
        G.load_cell.cache_clear()


# ---- ground truth -----------------------------------------------------------------


def check_truth_se_target(tmp):
    """target 0.036 on the SE domain (a <= 1, T <= 128): 4 days give 0.058, 32 days 0.032 (G0, 2 x 10^4
    runs per day). The run must stop with the target met, and the saved per-day moments must agree."""
    target = 0.036
    r = TR.run_truth("G0", cell(1.0), {"beta": 0.4}, N=20_000, target_se=target, min_days=4, chunk=4, max_days=64,
                     workers=1, out_dir=tmp, n_collapse=0, log=quiet)
    assert r.se_met and r.se_max < target and r.n_days < 64
    assert TR.max_slope_se(TR.jackknife(*r.days)) < target
    z = np.load(tmp / "G0_TEST.npz")
    assert len(z["Vd"]) == r.n_days and json.loads(str(z["meta"]))["se_met"]


def check_truth_resumes(tmp):
    kw = dict(N=5000, min_days=4, max_days=4, chunk=4, workers=1, out_dir=tmp, n_collapse=0)
    a = TR.run_truth("G0", cell(1.0), {"beta": 0.4}, log=quiet, **kw)
    seen = []
    b = TR.run_truth("G0", cell(1.0), {"beta": 0.4}, log=seen.append, **kw)
    assert any("resuming from 4 days" in s for s in seen)
    assert np.array_equal(a.curves["k"], b.curves["k"], equal_nan=True)


CHECKS = {f.__name__: f for f in [
    check_sign_acf_pooling_within_day, check_power_fit, check_sign_calibration, check_binned_fft_equals_direct,
    check_gt_c0_is_g0, check_g0_monotone, check_gf_zeta_one, check_day_lengths_from_cell, check_truth_se_target,
    check_truth_resumes,
]}
NEEDS_TMP = {"check_day_lengths_from_cell", "check_truth_se_target", "check_truth_resumes"}


def run_check(name, tmp):
    return CHECKS[name](tmp) if name in NEEDS_TMP else CHECKS[name]()


@pytest.mark.parametrize("name", list(CHECKS))
def test_gen3(name, tmp_path):
    run_check(name, tmp_path)


# ---- planted bugs -----------------------------------------------------------------


class Bug(NamedTuple):
    id: str
    module: object
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("sign ACF pooled across day boundaries", C3, "pooled_acf",
        "for i, x in enumerate(series):",
        "for i, x in enumerate([np.concatenate(series)] + [np.zeros(0)] * (len(series) - 1)):",
        ("check_sign_acf_pooling_within_day",)),
    Bug("fGn calibration without the arcsine inversion", G, "fit_sign_model",
        "l, r = lags[use], np.sin(np.pi * C[use] / 2)", "l, r = lags[use], C[use]",
        ("check_sign_calibration",)),
    Bug("G0 kernel made size-dependent", G, "simulate_g0",
        "lambda x, q: causal_conv(x, M.trade_kernel(len(x), beta))",
        "lambda x, q: binned_impact(x, np.searchsorted(size_bins(cell)[0], np.log(q), side='right'), "
        "np.linspace(0.8, 0.0, N_BINS))",
        ("check_g0_monotone",)),
    Bug("per-bin kernel off by one lag", G, "binned_impact",
        "M.trade_kernel(n, beta))", "np.roll(M.trade_kernel(n, beta), 1))",
        ("check_binned_fft_equals_direct", "check_gt_c0_is_g0")),
    Bug("bins sharing a beta convolved once each", G, "binned_impact",
        "for beta in np.unique(betas):", "for beta in betas:",
        ("check_binned_fft_equals_direct", "check_gt_c0_is_g0")),
    Bug("GF informed runs carry the herding signs", G, "simulate_gf",
        "eps = np.where(informed, u, herd).astype(np.int8)", "eps = herd",
        ("check_gf_zeta_one",)),
    Bug("day lengths from the wrong cell", G, "load_cell",
        'z = np.load(phase3_dir() / f"calibration_{name}.npz")',
        'z = np.load(sorted(phase3_dir().glob("calibration_*.npz"))[0])',
        ("check_day_lengths_from_cell",)),
    Bug("ground truth stops after the first days", TR, "run_truth",
        "if se_max < target_se or n >= max_days:", "if True:",
        ("check_truth_se_target",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate, tmp_path):
    mutate(bug.module, bug.func, bug.old, bug.new)
    for i, name in enumerate(bug.caught_by):
        d = tmp_path / str(i)
        d.mkdir()
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            run_check(name, d)
