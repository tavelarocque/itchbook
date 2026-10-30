"""Tests for Phase 3c-3: stationary days (panel3c), v pooled as the mean of stock logs (pooled), the
rho* criterion (phase3c3), the semi-synthetic response across days, and the stable fGn autocorrelation.
No data files are read.

    uv run --group dev --group analysis pytest tools/est/test_3c3.py
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np
import pytest

import gen3 as G
import panel3c as PC
import phase3c3 as P3C3
import pooled as PO

STRONG = {"w1": 0.6, "H": 0.99, "w2": 0.0, "phi": 0.5}  # a strong, slowly varying sign component
CELL = G.Cell("TEST", np.arange(100, 200), 3.0, 1.0, (0.3, 0.7), STRONG)


def check_burnin_discarded():
    """Every stock loses exactly its burn-in: the panel holds sum(Nsd) runs, and cut_days drops the first
    `burn` runs of each stock's series."""
    x = np.arange(3 + 5 + 3 + 2)
    out = PC.cut_days(x, np.array([3 + 5, 3 + 2]), 3)
    assert out.tolist() == [3, 4, 5, 6, 7, 11, 12]
    Nsd = np.array([[50, 70], [30, 90]])
    p = PC.make_panel("G0", CELL, Nsd, np.random.default_rng(1), params={"beta": 0.4}, burn=500, clip=None)
    assert len(p.eps) == len(p.q) == len(p.lm) == Nsd.sum() and p.lengths.tolist() == [50, 70, 30, 90]


def check_days_from_one_series():
    """With a strong slowly varying sign component, a stock's two consecutive days share it: their mean signs
    correlate across 300 stocks (above 0.5); fresh-start days do not (below 0.25 in absolute value)."""
    Nsd = np.full((300, 2), 400)
    for burn, lo, hi in ((2000, 0.5, 1.0), (None, -0.25, 0.25)):
        p = PC.make_panel("G0", CELL, Nsd, np.random.default_rng(2), params={"beta": 0.4}, burn=burn, clip=None,
                          continuous=True)
        m = p.eps.reshape(300, 2, 400).mean(axis=2)
        r = np.corrcoef(m[:, 0], m[:, 1])[0, 1]
        assert lo < r < hi, (burn, r)


def check_v_mean_of_logs():
    """Two stocks with v = 2^(0.2 j) and 2^(0.6 j): zeta_V = 1 + the mean of the stocks' log slopes = 1.4
    exactly (the log of the cross-stock mean is not a power law)."""
    T = PO.T_GRID2.astype(float)
    nA, nT = len(PO.A_GRID2), len(T)
    j = np.arange(nT)
    V = np.stack([np.stack([T * 2.0 ** (s * j)] * 2) for s in (0.2, 0.6)])
    V = np.repeat(V[:, :, None, :], nA, axis=2)
    K = V.copy()
    Sm = np.broadcast_to(T, (2, 2, nT)).copy()
    raw = PO.Raw(V, K, Sm, np.full((2, 2, nA), np.nan), np.full((2, 2), 10**6))
    e = PO.estimate2(raw, [("all",)])
    assert np.allclose(e["zeta_V"][:, :, 0], 1.4), e["zeta_V"][:, :, 0]


def check_criterion():
    assert "Rratio" not in P3C3.criterion("GF")
    assert "Rratio" in P3C3.criterion("G0") and "Rratio" in P3C3.criterion("GT")
    assert {"zeta_V", "zeta_K", "pf"} <= set(P3C3.criterion("GF"))


def check_semi_carries_across_days():
    """The semi-synthetic mid on day 2 depends on day 1's runs (one convolution over the stock's days, with
    the prepended copy as burn-in); with fresh=True it does not."""
    rng = np.random.default_rng(3)
    Nsd = np.array([[200, 300]])
    eps, q = rng.choice([-1.0, 1.0], 500), np.exp(3 + rng.standard_normal(500))
    eps2 = eps.copy()
    eps2[:200] = -eps2[:200]  # change day 1 only
    kw = dict(params={"beta": 0.4, "noise": 0.0})
    d2 = lambda p: np.diff(p.lm[200:])  # noqa: E731  day 2's mid increments
    a = PC.semi_synthetic(eps, q, Nsd, "G0", CELL, np.random.default_rng(4), **kw)
    b = PC.semi_synthetic(eps2, q, Nsd, "G0", CELL, np.random.default_rng(4), **kw)
    assert not np.allclose(d2(a), d2(b))
    a = PC.semi_synthetic(eps, q, Nsd, "G0", CELL, np.random.default_rng(4), fresh=True, **kw)
    b = PC.semi_synthetic(eps2, q, Nsd, "G0", CELL, np.random.default_rng(4), fresh=True, **kw)
    assert np.allclose(d2(a), d2(b))


def check_fgn_rho_stable():
    """The stable form equals the direct formula at short lags and the asymptote H (2H - 1) l^(2H - 2) at
    l = 10^7, where the direct form has lost its digits."""
    l = np.array([0.0, 0.5, 1, 2, 7, 100, 1000])
    for H in (0.55, 0.75, 0.98):
        direct = 0.5 * ((l + 1) ** (2 * H) - 2 * l ** (2 * H) + np.abs(l - 1) ** (2 * H))
        assert np.allclose(G.fgn_rho(l, H), direct, rtol=1e-8, atol=1e-12)
        big = 1e7
        assert np.isclose(G.fgn_rho(np.array([big]), H)[0], H * (2 * H - 1) * big ** (2 * H - 2), rtol=1e-6)


CHECKS = {f.__name__: f for f in [
    check_burnin_discarded, check_days_from_one_series, check_v_mean_of_logs, check_criterion,
    check_semi_carries_across_days, check_fgn_rho_stable,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_3c3(name):
    CHECKS[name]()


class Bug(NamedTuple):
    id: str
    module: object
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("burn-in not discarded", PC, "cut_days",
        "keep[s0 : s0 + burn] = False", "keep[s0 : s0] = False", ("check_burnin_discarded",)),
    Bug("days cut from separate series", PC, "_stationary_group",
        "series = burn + Nsd.sum(axis=1).astype(np.int64)", "series = burn + Nsd.ravel().astype(np.int64)",
        ("check_days_from_one_series",)),
    Bug("v pooled as the log of the cross-stock mean", PO, "core2",
        'out["zeta_V"] = zeta_mean_of_logs(v, ws, Ns, specs, T_grid, ranges)',
        'out["zeta_V"], _, out["kept_V"] = chain_zeta(v, ws, Mt, T_grid, ranges)',
        ("check_v_mean_of_logs",)),
    Bug("R ratio counted under GF", P3C3, "criterion",
        '(() if gen == "GF" else ("Rratio",))', '("Rratio",)', ("check_criterion",)),
    Bug("semi-synthetic fresh start per day", PC, "semi_synthetic",
        "fresh: bool = False", "fresh: bool = True", ("check_semi_carries_across_days",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(bug.module, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
