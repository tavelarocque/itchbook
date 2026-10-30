"""Tests for describe.py on synthetic inputs only (no market data): the Hurwitz zeta,
the discrete power-law MLE, sampler, KS distance and CSN fit; eta on hand-made and
simulated mid paths; per-day run statistics; hidden-print counts; eligibility; strata,
merges and the split. Then planted bugs.

    uv run --group dev --group analysis pytest tools/est/test_describe.py
"""
from __future__ import annotations

import gzip
import math
import tempfile
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pytest

import describe as DS
import runs as R
import tape

LO, HI = R.SESSION_LO, R.SESSION_HI
MID = (LO + HI) // 2


def runs_array(rows):
    a = np.zeros(len(rows), dtype=R.RUN_DTYPE)
    for i, r in enumerate(rows):
        d = dict(locate=1, sign=1, mid_valid=1, n_fills=1, ts=MID + i, q=100, bid_px=99_900, ask_px=100_100,
                 mid2=200_000, first_px=100_000, last_px=100_000)
        d.update(r)
        for k, v in d.items():
            a[i][k] = v
    return a


def tape_array(rows):
    a = np.zeros(len(rows), dtype=tape.numpy_dtype())
    for i, r in enumerate(rows):
        d = dict(ts=MID, type=b"P", phase=tape.REGULAR, tstate=b"T", locate=1, price=100_000, bid_px=99_900,
                 ask_px=100_100, shares=100, side=b"B", sign=0, seq=i)
        d.update(r)
        for k, v in d.items():
            a[i][k] = v
    return a


# ---- power law ----------------------------------------------------------------------


def check_hurwitz():
    assert abs(DS.hurwitz(2.0, 1.0) - math.pi**2 / 6) < 1e-12
    assert abs(DS.hurwitz(1.5, 1.0) - 2.612375348685488) < 1e-12
    assert abs(DS.hurwitz(3.0, 0.5) - 7 * 1.2020569031595942) < 1e-12  # zeta(s, 1/2) = (2^s - 1) zeta(s)
    a = 37.0
    direct = np.sum((a + np.arange(2_000_000)) ** -2.5) + (a + 2e6 - 0.5) ** -1.5 / 1.5
    assert abs(DS.hurwitz(2.5, a) / direct - 1) < 1e-12
    assert np.allclose(DS.hurwitz(np.array([2.0, 3.0]), np.array([1.0, 4.0])),
                       [math.pi**2 / 6, 1.2020569031595942 - 1 - 1 / 8 - 1 / 27])


def check_alpha_mle_is_exact():
    x = DS.sample_powerlaw(5000, 2.3, 7, np.random.default_rng(1))
    al = DS.alpha_mle(np.array([7.0]), np.array([len(x)]), np.array([np.log(x).sum()]))[0]
    g = np.linspace(2.0, 2.6, 6001)
    L = -len(x) * np.log(DS.hurwitz(g, 7.0)) - g * np.log(x).sum()
    assert abs(al - g[np.argmax(L)]) < 2e-4


def check_sampler_is_exact():
    """pmf of 400000 draws within 4 binomial SE of x^-alpha / zeta(alpha, x_min) at small x,
    where the continuous approximation used as a start is worst."""
    rng = np.random.default_rng(2)
    for alpha, xmin in ((2.0, 1), (3.0, 5)):
        x = DS.sample_powerlaw(400_000, alpha, xmin, rng)
        z = DS.hurwitz(alpha, xmin)
        for k in (xmin, xmin + 1, xmin + 2, xmin + 9):
            p = k**-alpha / z
            assert abs(np.mean(x == k) - p) < 4 * math.sqrt(p * (1 - p) / len(x)), (alpha, xmin, k)


def check_ks_distance_brute_force():
    rng = np.random.default_rng(3)
    # a round-lot-like lattice (multiples of 10) in the tail: the model CDF rises between observed values
    q = np.concatenate([rng.integers(1, 6, 200), DS.sample_powerlaw(300, 2.4, 1, rng) * 10])
    u, c = np.unique(q.astype(float), return_counts=True)
    for k in (0, 2, 5):
        alpha = 2.1 + 0.1 * k
        xmin = u[k]
        tail = q[q >= xmin]
        xs = np.arange(int(xmin), int(u[-1]) + 1)
        pmf = xs ** -alpha / DS.hurwitz(alpha, xmin)
        F = np.cumsum(pmf)
        Femp = np.searchsorted(np.sort(tail), xs, side="right") / len(tail)
        want = np.max(np.abs(Femp - F))
        got = DS.ks_distance(u, c, np.array([k]), np.array([alpha]))[0]
        assert abs(got - want) < 1e-10, (k, got, want)


def check_csn_recovers_pareto():
    """Measured over 10 seeds: pure discrete Pareto alpha = 2.5, x_min = 1, n = 20000: alpha
    2.4947, sd 0.012, x_min = 1 always, Vuong +20.5 to +21.6. Uniform body below 100 plus a
    power-law tail alpha = 2.2 from 100 (6000 + 6000): alpha 2.195, sd 0.021, x_min 95 to 161
    over 11 seeds (the uniform body continues the tail's density near 100). Bootstrap
    p-values on true power laws were uniform (2 of 30 below 0.1). Tolerances:
    alpha 0.05 and 0.08 (4 sd), x_min 90 to 200, p > 0.01."""
    r = np.random.default_rng(4)
    f = DS.csn(DS.sample_powerlaw(20000, 2.5, 1, r), r, 50)
    assert abs(f["csn_alpha"] - 2.5) < 0.05 and f["csn_xmin"] == 1 and f["vuong"] > 2, f
    r = np.random.default_rng(5)
    f = DS.csn(np.concatenate([r.integers(1, 100, 6000), DS.sample_powerlaw(6000, 2.2, 100, r)]), r, 50)
    assert abs(f["csn_alpha"] - 2.2) < 0.08 and 90 <= f["csn_xmin"] <= 200 and f["csn_p"] > 0.01, f
    small = DS.csn(DS.sample_powerlaw(1500, 2.5, 1, r), r, 50)
    assert math.isnan(small["csn_p"]) and not math.isnan(small["csn_alpha"])  # p only from 2000 runs
    assert all(math.isnan(v) for v in DS.csn(np.arange(1, 40), r, 50).values())  # below 50 runs: blank


def check_csn_rejects_lognormal():
    """Log-normal q (mu = 4, sigma = 1, n = 50000). The KS-chosen x_min often lands far in the
    tail, where a log-normal is locally power-law-like, so the KS p-value alone rejects in only
    about half the seeds (measured 5 of 10 with p < 0.1), as CSN warn; the Vuong statistic
    favoured the log-normal in 20 of 20 (sigma 0.7 and 1, max -0.08). Required on 8 seeds:
    Vuong < 0 on every seed and at least 2 KS rejections (chance of fewer at rate 0.5: 3.5%)."""
    rej = 0
    for s in range(8):
        r = np.random.default_rng(600 + s)
        f = DS.csn(np.maximum(1, np.round(np.exp(4 + r.standard_normal(50000)))), r, 50)
        assert f["vuong"] < 0, (s, f)
        rej += f["csn_p"] < 0.1
    assert rej >= 2, rej


# ---- eta, run statistics, hidden prints ---------------------------------------------


def check_eta_hand_example():
    # locate 1 mid2: +1 tick, +1 (continues), half tick (skipped), -1 (reverses the last one-tick +),
    # 0, -1 (continues), +3 ticks (skipped), +1 (reverses); locate 2 starts fresh
    m2 = [1_020_000 + x for x in (0, 200, 400, 300, 100, 100, -100, 500, 700)] + [3_000_000, 3_000_200]
    loc = [1] * 9 + [2] * 2
    nc, na, n1, nm = DS.eta_counts(np.array(loc), np.array(m2), 3)
    assert (nc[1], na[1], n1[1], nm[1]) == (2, 2, 5, 7)
    assert (nc[2], na[2], n1[2], nm[2]) == (0, 0, 1, 1)
    assert DS.eta(nc[1], na[1]) == 0.5
    assert np.isnan(DS.eta(3, 0)) and np.isnan(DS.eta(0, 0))  # undefined without an alternation, never inf
    sub = DS.eta_counts(np.array([1, 1, 1]), np.array([4000, 4002, 4004]), 2)  # below $1 the tick is 1 unit
    assert (sub[0][1], sub[1][1]) == (1, 0)


def check_eta_on_simulated_paths():
    """One-tick moves continue the previous one with probability c, mixed with zero and half-tick
    moves; eta = c / (2 (1 - c)). Tolerance 4 delta-method SE, sqrt(c (1 - c) / N) / (2 (1 - c)^2)."""
    rng = np.random.default_rng(7)
    for c in (0.2, 0.5, 0.7):
        n = 20000
        step = np.empty(n)
        step[0] = 1
        for i in range(1, n):
            step[i] = step[i - 1] if rng.random() < c else -step[i - 1]
        # each one-tick move is followed by a zero or half-tick move, which eta must skip
        seq = np.empty(2 * n)
        seq[0::2] = step * 200
        seq[1::2] = np.where(rng.random(n) < 0.5, 0, rng.choice([-100, 100], n))
        m2 = 1_000_000 + np.cumsum(seq).astype(np.int64)
        nc, na, n1, _ = DS.eta_counts(np.ones(len(m2), dtype=int), m2, 2)
        want = c / (2 * (1 - c))
        se = math.sqrt(c * (1 - c) / n) / (2 * (1 - c) ** 2)
        assert abs(DS.eta(nc[1], na[1]) - want) < 4 * se, (c, DS.eta(nc[1], na[1]), want)


def check_day_stats():
    r = runs_array([
        {"q": 100, "first_px": 100_000, "mid2": 200_000},
        {"q": 300, "first_px": 100_200, "mid2": 200_200, "bid_px": 100_000, "ask_px": 100_200},
        {"q": 50, "first_px": 100_000, "mid_valid": 0, "mid2": 0, "bid_px": 100_000, "ask_px": 100_000},
        {"locate": 3, "q": 10, "first_px": 50_000, "mid2": 100_000, "bid_px": 49_900, "ask_px": 50_100},
    ])
    st = DS.day_stats(r)
    assert st["n"][1] == 3 and st["sum_q"][1] == 450 and st["n"][3] == 1
    assert np.isclose(st["sum_pq"][1] / st["sum_q"][1], (100 * 100_000 + 300 * 100_200 + 50 * 100_000) / 450)
    assert st["n_quote"][1] == 2 and st["n_spread1"][1] == 0 and np.isclose(st["sum_spread_ticks"][1], 4.0)
    assert st["nm"][1] == 1 and st["n1"][1] == 1  # mid 10.00 to 10.01; the invalid mid is skipped
    assert np.isclose(st["sum_spread_ticks"][3], 2.0)


def check_hidden_prints():
    t = tape_array([
        {"price": 100_000},  # at its own mid, no run there: unsigned
        {"price": 100_100},  # above mid: counted, signed
        {"price": 100_000, "ts": MID + 7},  # at mid, shares a run's (locate, ts): signed by the run
        {"ts": LO - 1},  # before 10:00
        {"ts": HI},  # 15:30 is outside
        {"phase": tape.PRE},
        {"tstate": b"H"},
        {"type": b"E", "sign": 1},
        {"price": 100_000, "bid_px": 100_000, "ask_px": 100_000},  # locked quote: not at a mid, unsigned
    ])
    runs = runs_array([{"ts": MID + 7, "sign": -1}])
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "x.tape.gz"
        with gzip.open(path, "wb") as f:
            f.write(tape.HEADER.pack(tape.MAGIC, tape.VERSION, tape.RECORD.size, 0x01020304, 0))
            f.write(t.tobytes())
        chunks = list(DS.iter_tape_chunks(path, records_per_chunk=4))
    assert [len(c) for c in chunks] == [4, 4, 1]
    p = DS.hidden_prints(chunks)
    assert len(p) == 4
    hs = DS.hidden_stats(p, runs, 2)
    assert (hs["hid_n"][1], hs["hid_mid"][1], hs["hid_uns"][1]) == (4, 2, 2)


# ---- eligibility, strata, split -----------------------------------------------------


def check_eligibility_filters():
    good = {"mc": "Q", "ic": "C", "etp": "N", "n": 500, "vwap": 200_000.0}
    days = [dict() for _ in range(8)]
    for i, d in enumerate(days):
        d["GOOD"] = dict(good)
        d["NYSE"] = dict(good, mc="N")
        d["MOVE"] = dict(good, mc="G" if i < 4 else "S")  # NASDAQ tier change is still NASDAQ-listed
        d["ARCA"] = dict(good, mc="P")
        d["PREF"] = dict(good, ic="P" if i == 3 else "C")
        d["ETF"] = dict(good, etp="Y")
        if i != 5:
            d["MISS"] = dict(good)  # absent from the symbol table one day
        d["IDLE"] = dict(good, n=0 if i == 2 else 500, vwap=0.0 if i == 2 else 200_000.0)
        d["PENNY"] = dict(good, vwap=49_999.0 if i == 7 else 200_000.0)
        d["EDGE"] = dict(good, vwap=50_000.0)  # exactly $5 passes
    f = DS.eligibility(days)
    assert f["GOOD"] == [] and f["NYSE"] == [] and f["MOVE"] == [] and f["EDGE"] == []
    assert f["ARCA"] == ["listing"] and f["PREF"] == ["issue"] and f["ETF"] == ["etp"]
    assert f["IDLE"] == ["presence", "price"] and f["PENNY"] == ["price"]
    assert set(f["MISS"]) == {"listing", "issue", "etp", "presence", "price"}


def check_strata_edges():
    assert DS.eta_bin(np.array([0.0, 0.3499, 0.35, 0.4999, 0.5, 2.0])).tolist() == [1, 1, 2, 2, 3, 3]
    t = DS.terciles(np.arange(9.0))
    assert t.tolist() == [1, 1, 1, 2, 2, 2, 3, 3, 3]
    assert [DS.n_explore(n) for n in (10, 12, 15, 17)] == [3, 4, 5, 5]


def check_merges():
    sizes = {(1, 1): 40, (1, 2): 12, (1, 3): 30, (2, 1): 5, (2, 2): 6, (2, 3): 50, (3, 1): 20, (3, 2): 20, (3, 3): 20}
    eb = np.concatenate([[e] * n for (e, t), n in sizes.items()])
    st = np.concatenate([[t] * n for (e, t), n in sizes.items()])
    cells, merges = DS.merge_cells(eb, st)
    lab = {(e, t): cells[(eb == e) & (st == t)][0] for e, t in sizes}
    assert lab[(1, 2)] == lab[(1, 3)] == "E1S2+3" and lab[(1, 1)] == "E1S1"  # 12 stocks: 8 to confirm
    assert lab[(2, 1)] == lab[(2, 2)] == lab[(2, 3)] == "E2S1+2+3"  # 5 into 6, then 11 into 50
    assert {lab[(3, t)] for t in (1, 2, 3)} == {"E3S1", "E3S2", "E3S3"}
    assert len(merges) == 3


def check_stratify_leaves_out_undefined_eta():
    rng = np.random.default_rng(8)
    rows = [{"symbol": f"N{i:03d}", "eligible": 1, "listing": "QGS"[i % 3], "eta": float(rng.uniform(0.2, 0.8)),
             "sigma_l": float(rng.uniform(1.2, 2.0))} for i in range(120)]
    rows += [{"symbol": "UNDEF", "eligible": 1, "listing": "Q", "eta": float("nan"), "sigma_l": 9.0},
             {"symbol": "NYSE1", "eligible": 1, "listing": "N", "eta": 0.4, "sigma_l": 1.5},
             {"symbol": "INELG", "eligible": 0, "listing": "Q", "eta": 0.4, "sigma_l": 1.5}]
    meta = DS.stratify(rows)
    by = {r["symbol"]: r for r in rows}
    assert meta["unstratified"] == ["UNDEF"]
    for s in ("UNDEF", "NYSE1", "INELG"):
        assert "cell" not in by[s] and "split" not in by[s]
    nas = [r for r in rows if r["listing"] != "N" and r["eligible"]]
    assert np.allclose(meta["sigma_terciles"], np.quantile([r["sigma_l"] for r in nas], [1 / 3, 2 / 3]))
    for r in rows[:120]:
        assert r["eta_bin"] == DS.eta_bin(np.array([r["eta"]]))[0] and r["split"] in ("explore", "confirm")


def check_split_is_deterministic_and_stratified():
    syms = [f"S{i:03d}" for i in range(150)]
    cells = np.array(["E1S1"] * 40 + ["E1S2"] * 61 + ["E2S1+2"] * 49)
    a = DS.split(syms, cells)
    assert a == DS.split(list(reversed(syms)), cells[::-1])  # input order does not matter
    assert a == DS.split(syms, cells, seed=DS.SPLIT_SEED)
    assert a != DS.split(syms, cells, seed=DS.SPLIT_SEED + 1)
    for cell, n in (("E1S1", 40), ("E1S2", 61), ("E2S1+2", 49)):
        members = [s for s, c in zip(syms, cells, strict=True) if c == cell]
        n_exp = sum(a[s] == "explore" for s in members)
        assert n_exp == DS.n_explore(n) and len(members) == n
        assert [a[s] for s in sorted(members)][: n_exp] != ["explore"] * n_exp  # not the first names


CHECKS = {f.__name__: f for f in [
    check_hurwitz, check_alpha_mle_is_exact, check_sampler_is_exact, check_ks_distance_brute_force,
    check_csn_recovers_pareto, check_csn_rejects_lognormal, check_eta_hand_example, check_eta_on_simulated_paths,
    check_day_stats, check_hidden_prints, check_eligibility_filters, check_strata_edges, check_merges,
    check_stratify_leaves_out_undefined_eta,
    check_split_is_deterministic_and_stratified,
]}


@pytest.mark.parametrize("name", list(CHECKS))
def test_describe(name):
    CHECKS[name]()


# ---- planted bugs -----------------------------------------------------------------


class Bug(NamedTuple):
    id: str
    func: str
    old: str
    new: str
    caught_by: tuple[str, ...]


BUGS = [
    Bug("eta without the factor 2", "eta",
        "nc / (2 * na), np.nan", "nc / na, np.nan",
        ("check_eta_hand_example", "check_eta_on_simulated_paths")),
    Bug("eta over all mid moves, not one-tick moves", "eta_counts",
        "p = np.flatnonzero(one)", "p = np.flatnonzero(moves)",
        ("check_eta_hand_example", "check_eta_on_simulated_paths")),
    Bug("one tick of the mid taken as one tick of mid2", "eta_counts",
        "np.abs(dm) == 2 * tick_of(m2[:-1])", "np.abs(dm) == tick_of(m2[:-1])",
        ("check_eta_hand_example", "check_eta_on_simulated_paths")),
    Bug("eta pairs across stocks", "eta_counts",
        "pair = lp[1:] == lp[:-1]", "pair = np.ones(len(lp) - 1, dtype=bool)",
        ("check_eta_hand_example",)),
    Bug("Euler-Maclaurin terms with the wrong sign", "hurwitz",
        "out += B / fact * rising * pw", "out -= B / fact * rising * pw",
        ("check_hurwitz",)),
    Bug("KS only at observed values", "ks_distance",
        "np.maximum(np.abs(F_emp - F_at), np.abs(F_emp_prev - F_below))", "np.abs(F_emp - F_at)",
        ("check_ks_distance_brute_force",)),
    Bug("continuous approximation instead of the exact discrete MLE", "fit_powerlaw",
        "al = alpha_mle(u[i0], n_t, s_l)", "al = 1 + n_t / (s_l - n_t * np.log(u[i0] - 0.5))",
        ("check_csn_recovers_pareto",)),
    Bug("sampler without the exact correction", "sample_powerlaw",
        "for _ in range(200):", "for _ in range(0):",
        ("check_sampler_is_exact",)),
    Bug("bootstrap counts fits at least as good, not as bad", "ks_pvalue",
        "worse += fit_powerlaw(us, cs)[3] >= ks", "worse += fit_powerlaw(us, cs)[3] <= ks",
        ("check_csn_rejects_lognormal",)),
    Bug("Vuong ratio with the sign flipped", "vuong",
        "d = l_pl - l_ln", "d = l_ln - l_pl",
        ("check_csn_recovers_pareto", "check_csn_rejects_lognormal")),
    Bug("presence on any day instead of every day", "eligibility",
        'if not all(r is not None and r.get("n", 0) > 0 for r in rows):',
        'if not any(r is not None and r.get("n", 0) > 0 for r in rows):',
        ("check_eligibility_filters",)),
    Bug("price filter strict at $5", "eligibility",
        'r["vwap"] >= PRICE_MIN', 'r["vwap"] > PRICE_MIN',
        ("check_eligibility_filters",)),
    Bug("NYSE-listed dropped", "eligibility",
        'r["mc"] in NASDAQ + NYSE', 'r["mc"] in NASDAQ',
        ("check_eligibility_filters",)),
    Bug("eta of a stock without alternations is infinite", "eta",
        "np.where(na > 0, nc / (2 * na), np.nan)", "nc / (2 * na)",
        ("check_eta_hand_example",)),
    Bug("stocks with undefined eta stratified (NaN falls in the first eta bin)", "stratify",
        'strat = [r for r in nas if np.isfinite(r["eta"])]', "strat = nas",
        ("check_stratify_leaves_out_undefined_eta",)),
    Bug("split ignores the seed", "split",
        "rng = np.random.default_rng(seed)", "rng = np.random.default_rng()",
        ("check_split_is_deterministic_and_stratified",)),
    Bug("split not randomized (alphabetical)", "split",
        "perm = rng.permutation(len(members))", "perm = np.arange(len(members))",
        ("check_split_is_deterministic_and_stratified",)),
    Bug("explore count truncated, not rounded", "n_explore",
        "int(math.floor(EXPLORE_FRAC * n + 0.5))", "int(EXPLORE_FRAC * n)",
        ("check_strata_edges",)),
    Bug("merge on cell size instead of confirmation count", "merge_cells",
        "size[g] - n_explore(size[g]) < MIN_CONFIRM", "size[g] < MIN_CONFIRM",
        ("check_merges",)),
    Bug("at-mid test without the factor 2", "hidden_stats",
        '2 * p["price"].astype(np.int64) == bid + ask', 'p["price"].astype(np.int64) == bid + ask',
        ("check_hidden_prints",)),
    Bug("hidden prints not filtered on trading state", "hidden_prints",
        ' & (t["tstate"] == b"T")', "",
        ("check_hidden_prints",)),
]


@pytest.mark.parametrize("bug", BUGS, ids=[b.id for b in BUGS])
def test_planted_bug_is_caught(bug, mutate):
    mutate(DS, bug.func, bug.old, bug.new)
    for name in bug.caught_by:
        with pytest.raises(Exception):  # noqa: B017 (an assertion or a crash both count as caught)
            CHECKS[name]()
