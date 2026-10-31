# Phase 3c-4 report: independent stationary days, the widest interval, the rules of 31 Oct 2026, and rho*

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken here. Real data were read
only through `explore.load_explore_runs`, and only signs and sizes (step 1). NASDAQ cells only.

## Headlines

1. **Cross-day dependence (Measurement).** A stock's daily mean sign is uncorrelated across its 8 study
   days: the correlation is about 0, with an interval containing 0, in 8 of 9 cells (Q-T2S3, with 7 stocks,
   gives -0.11). The centred sign ACF at lag 10 is mostly 0. The daily mean log size is strongly persistent
   in every cell (+0.12 to +0.74), at every calendar distance: that is a stable stock-level size effect,
   not a slowly varying flow component.
2. **Independent stationary days work (Measurement).** Each stock-day is its own series after its own
   discarded burn-in. G0's zeta_K(0) and zeta_V(0) on [8, 64] at Q-T1S3, with N = 300 and 3000, sit within
   1.3 SE of the 10^6-run truth.
3. **The burn-in check keeps B = 20,000 (Measurement and Choice).** B = 1000, 2000 and 5000 all fail on
   zeta_K(0), which sits 2.4 to 6.8 SE above the truth on [8, 64]. zeta_V passes at every B.
4. **The full study ran (Measurement).** 200 panels per cell and generator; top-up to 400 panels for the 18
   (cell, generator) pairs with an element inside the Monte Carlo band; 100 panels for the continuous-series
   stress cell. 201 min of wall time on 3 workers, 12.2 worker-hours. Under the pre-registered rules rho* was
   **none in every cell**.
5. **Diagnosis (Measurement).** `docs/phase3c4_diagnosis.md` (read-only, saved panels). There were three
   failure modes. zeta_K was "estimated" at small ceilings from the lower horizons of a range that no stock
   reached (83 failures, all with 0 stocks at T2). zeta_V had small negative bias and too-narrow coverage at
   large ceilings. The GF P-F slope failed on bias alone, with coverage 0.85 to 1.00.
6. **Three rule changes (Choice, 31 Oct 2026, from the diagnosis).** (1) Range reach: an estimate over
   [T1, T2] is defined for a panel only if at least 5 stocks reach T2 under the ceiling. An element is
   evaluable only if defined in at least 90 % of panels. (2) The P-F slope is judged by coverage only.
   (3) rho* is reported per cell and per tick bin. The bias bound (0.02), the coverage bound (0.85), the
   Monte Carlo band and the panels are unchanged, and no simulation was rerun.
7. **rho* under the new rules (Measurement).** Per cell: Q-T1S1 1/5, Q-T1S2 1/20, Q-T1S3 none, Q-T2S1 none,
   Q-T2S2 none, Q-T2S3 1/5, Q-T3S1 1/5, Q-T3S2 1/5, Q-T3S3 1/10. Per tick bin: **Q-T1 none, Q-T2 none,
   Q-T3 1/10**. Every remaining failure is zeta_V(a) coverage (161) or R_{1/2}/R_0 (14). zeta_K and the P-F
   slope have none. The continuous-series stress cell gets 1/40, against 1/5 for the same cell on
   independent days.
8. **Widest interval (Measurement).** Over defined elements, the day-and-stock bootstrap was widest in 98 to
   99 % of element-panels and the day jackknife-t in the rest; the whole-stock bootstrap essentially never.

## 1. Cross-day dependence (exploration stocks, signs and sizes only)

Per stock-day: (a) mean sign, (b) mean log run size, (c) sign ACF at lag 10 centred on the day's own mean
sign (stock-days with at least 50 runs).

**Choice: the dependence measure.** For each pair of study days, the correlation across the cell's stocks
of the two days' values, after removing each day's cross-sectional mean, averaged over the 28 pairs and
within calendar-distance groups. A persistent stock effect counts as cross-day dependence here: with 8
days per stock it cannot be told apart from slow time variation. Intervals come from a 95 % stock
bootstrap (B = 999).

| cell | statistic | stocks | correlation [95 % CI] | excludes 0 | <= 45 days apart | 46 to 200 | 201 to 400 |
|---|---|---|---|---|---|---|---|
| Q-T1S1 | mean sign | 69 | -0.035 [-0.076, +0.019] | no | +0.029 (3 pairs) | -0.014 (16) | -0.093 (9) |
| Q-T1S1 | mean log size | 69 | +0.207 [+0.091, +0.314] | yes | +0.305 (3 pairs) | +0.240 (16) | +0.117 (9) |
| Q-T1S1 | centred sign ACF(10) | 69 | +0.031 [-0.038, +0.127] | no | +0.052 (3 pairs) | -0.004 (16) | +0.086 (9) |
| Q-T1S2 | mean sign | 92 | +0.012 [-0.029, +0.067] | no | -0.136 (3 pairs) | +0.056 (16) | -0.015 (9) |
| Q-T1S2 | mean log size | 92 | +0.349 [+0.198, +0.479] | yes | +0.483 (3 pairs) | +0.322 (16) | +0.352 (9) |
| Q-T1S2 | centred sign ACF(10) | 92 | +0.035 [-0.011, +0.084] | no | +0.092 (3 pairs) | +0.027 (16) | +0.030 (9) |
| Q-T1S3 | mean sign | 108 | +0.008 [-0.044, +0.079] | no | -0.040 (3 pairs) | -0.003 (16) | +0.044 (9) |
| Q-T1S3 | mean log size | 108 | +0.121 [+0.048, +0.205] | yes | +0.157 (3 pairs) | +0.150 (16) | +0.056 (9) |
| Q-T1S3 | centred sign ACF(10) | 108 | +0.081 [+0.011, +0.151] | yes | +0.023 (3 pairs) | +0.125 (16) | +0.021 (9) |
| Q-T2S1 | mean sign | 35 | -0.003 [-0.076, +0.059] | no | +0.267 (3 pairs) | -0.014 (16) | -0.073 (9) |
| Q-T2S1 | mean log size | 35 | +0.547 [+0.222, +0.691] | yes | +0.563 (3 pairs) | +0.570 (16) | +0.502 (9) |
| Q-T2S1 | centred sign ACF(10) | 35 | +0.104 [-0.057, +0.301] | no | +0.223 (3 pairs) | +0.060 (16) | +0.145 (9) |
| Q-T2S2 | mean sign | 21 | +0.012 [-0.066, +0.064] | no | +0.182 (3 pairs) | +0.081 (16) | -0.168 (9) |
| Q-T2S2 | mean log size | 21 | +0.342 [+0.180, +0.445] | yes | +0.553 (3 pairs) | +0.395 (16) | +0.178 (9) |
| Q-T2S2 | centred sign ACF(10) | 21 | -0.058 [-0.100, -0.010] | yes | -0.015 (3 pairs) | -0.038 (16) | -0.107 (9) |
| Q-T2S3 | mean sign | 7 | -0.111 [-0.137, -0.044] | yes | -0.264 (3 pairs) | -0.168 (16) | +0.042 (9) |
| Q-T2S3 | mean log size | 7 | +0.258 [-0.028, +0.481] | no | +0.497 (3 pairs) | +0.213 (16) | +0.259 (9) |
| Q-T2S3 | centred sign ACF(10) | 7 | -0.084 [-0.145, +0.635] | no | +0.022 (3 pairs) | -0.145 (16) | -0.011 (9) |
| Q-T3S1 | mean sign | 22 | -0.047 [-0.092, +0.022] | no | +0.195 (3 pairs) | -0.073 (16) | -0.080 (9) |
| Q-T3S1 | mean log size | 22 | +0.531 [+0.314, +0.706] | yes | +0.723 (3 pairs) | +0.533 (16) | +0.464 (9) |
| Q-T3S1 | centred sign ACF(10) | 22 | -0.039 [-0.089, +0.020] | no | +0.174 (3 pairs) | -0.135 (16) | +0.059 (9) |
| Q-T3S2 | mean sign | 11 | +0.039 [-0.073, +0.195] | no | -0.157 (3 pairs) | +0.110 (16) | -0.022 (9) |
| Q-T3S2 | mean log size | 11 | +0.571 [+0.186, +0.787] | yes | +0.508 (3 pairs) | +0.529 (16) | +0.668 (9) |
| Q-T3S2 | centred sign ACF(10) | 11 | -0.065 [-0.106, +0.003] | no | -0.162 (3 pairs) | -0.060 (16) | -0.040 (9) |
| Q-T3S3 | mean sign | 10 | +0.103 [-0.055, +0.248] | no | -0.118 (3 pairs) | +0.106 (16) | +0.172 (9) |
| Q-T3S3 | mean log size | 10 | +0.737 [+0.563, +0.883] | yes | +0.734 (3 pairs) | +0.756 (16) | +0.702 (9) |
| Q-T3S3 | centred sign ACF(10) | 10 | +0.365 [-0.088, +0.584] | no | +0.337 (3 pairs) | +0.282 (16) | +0.522 (9) |

**Measurement.**
- **Mean sign:** no dependence. Every interval contains 0 except Q-T2S3, at -0.11 with 7 stocks.
- **Centred ACF(10):** two cells exclude 0 weakly (Q-T1S3 +0.08, Q-T2S2 -0.06).
- **Mean log size:** +0.12 to +0.74, about the same at every calendar distance. This is stock
  heterogeneity in size level, which the synthetic stocks (one size distribution per cell) do not carry.
  It is a Choice of the generators, now measured.

## 2. Independent stationary days

**Choice:** each stock-day is its own continuous series of B + N runs with the first B discarded,
independent of the stock's other days (`panel3c.make_panel(burn=B)`). Stock-level parameters are the
cell's on every day. The 3c-3 path (one series per stock) is kept as `continuous=True`, the stress test.

**Generation speed-up.** The FFT and circulant size is now the smallest even 5-smooth number a 2^k
(a in 1, 3, 5, 9, 15) at least 2N, instead of the next power of 2. Most stock-days use 49,152 points
instead of 65,536, about 25 % faster. The draws change, the distributions do not, and every 3c and 3c-3
test still passes.

**Stationarity check on the new path** (G0, Q-T1S3, [8, 64], a = 0; 8 independent groups; SE from a
delete-one-group jackknife):

| N | zeta_V(0) (SE) | zeta_K(0) (SE) |
|---|---|---|
| 300 | 1.364 (0.005) | 0.831 (0.019) |
| 3000 | 1.347 (0.011) | 0.814 (0.007) |
| 10^6 truth (3c-3) | 1.379 (0.022) | 0.811 (0.002) |

Every value is within 1.3 SE of the truth.

### Burn-in check (the follow-up instruction, step 1)

G0 and GT at Q-T1S3, both ranges, N = 300 (8000 days) and 3000 (1280 days), in 16 independent groups.
Every SE is at most 0.009. A value is "ok" when it is within 2 SE of the 10^6-run truth (both SEs
combined) and within 0.01 of the B = 20,000 value; "x" marks a failure.

| gen | N | estimand | B = 1000 | B = 2000 | B = 5000 | B = 20000 | truth (10^6-run days) |
|---|---|---|---|---|---|---|---|
| G0 | 300 | zeta_V(0) [4,32] | 1.361 (0.002) | 1.357 (0.003) | 1.358 (0.002) | 1.360 (0.003) | 1.372 (0.015) |
| G0 | 300 | zeta_K(0) [4,32] | 0.961 (0.003) x | 0.949 (0.004) | 0.953 (0.004) | 0.953 (0.004) | 0.946 (0.001) |
| G0 | 300 | zeta_V(0) [8,64] | 1.363 (0.003) | 1.357 (0.005) | 1.360 (0.003) | 1.357 (0.004) | 1.379 (0.022) |
| G0 | 300 | zeta_K(0) [8,64] | 0.857 (0.006) x | 0.824 (0.008) | 0.832 (0.008) x | 0.825 (0.007) | 0.811 (0.002) |
| G0 | 3000 | zeta_V(0) [4,32] | 1.364 (0.004) | 1.363 (0.003) | 1.362 (0.004) | 1.364 (0.004) | 1.372 (0.015) |
| G0 | 3000 | zeta_K(0) [4,32] | 0.957 (0.003) x | 0.958 (0.003) x | 0.950 (0.003) | 0.949 (0.003) | 0.946 (0.001) |
| G0 | 3000 | zeta_V(0) [8,64] | 1.368 (0.006) | 1.363 (0.005) | 1.364 (0.006) | 1.367 (0.006) | 1.379 (0.022) |
| G0 | 3000 | zeta_K(0) [8,64] | 0.842 (0.005) x | 0.837 (0.004) x | 0.826 (0.005) x | 0.813 (0.005) | 0.811 (0.002) |
| GT | 300 | zeta_V(0) [4,32] | 1.361 (0.002) | 1.357 (0.003) | 1.358 (0.002) | 1.360 (0.003) | 1.372 (0.015) |
| GT | 300 | zeta_K(0) [4,32] | 0.938 (0.003) x | 0.925 (0.004) | 0.931 (0.004) | 0.932 (0.005) | 0.924 (0.001) |
| GT | 300 | zeta_V(0) [8,64] | 1.363 (0.003) | 1.357 (0.005) | 1.360 (0.003) | 1.357 (0.004) | 1.379 (0.022) |
| GT | 300 | zeta_K(0) [8,64] | 0.838 (0.007) x | 0.802 (0.009) | 0.815 (0.009) x | 0.807 (0.008) | 0.793 (0.002) |
| GT | 3000 | zeta_V(0) [4,32] | 1.364 (0.004) | 1.363 (0.003) | 1.362 (0.004) | 1.364 (0.004) | 1.372 (0.015) |
| GT | 3000 | zeta_K(0) [4,32] | 0.934 (0.003) x | 0.936 (0.003) x | 0.929 (0.003) | 0.927 (0.003) | 0.924 (0.001) |
| GT | 3000 | zeta_V(0) [8,64] | 1.368 (0.006) | 1.363 (0.005) | 1.364 (0.006) | 1.367 (0.006) | 1.379 (0.022) |
| GT | 3000 | zeta_K(0) [8,64] | 0.823 (0.005) x | 0.818 (0.005) x | 0.808 (0.006) x | 0.793 (0.006) | 0.793 (0.002) |

**Measurement and Choice.**
- The smallest B with every value ok is **B = 20,000**, the value already in use.
- Every failure is zeta_K(0): the shorter burn-ins leave it 2.4 to 6.8 SE above the truth, mostly on
  [8, 64] and, at B <= 2000, also on [4, 32]. zeta_V(0) is ok at every B.
- The reduced-run panels used B = 20,000, so they were kept.

## 3. The interval

`pooled.intervals3` reports the widest of three 90 % intervals:
- (i) the day-and-stock bootstrap, using `inference.cluster_bootstrap`'s draws;
- (ii) the delete-one-day jackknife with Student t at G - 1 df;
- (iii) a whole-stock bootstrap that resamples stocks with all their days, B = 200.

A method undefined for an element is skipped for it. Which was widest, over the full study and the elements
defined under the range-reach rule (section 5):

| days | estimator | elements x panels | day_stock_bootstrap | day_jackknife_t | whole_stock_bootstrap |
|---|---|---|---|---|---|
| cont | Rratio | 1726 | 1.00 | 0.00 | 0.00 |
| cont | pf | 1726 | 0.99 | 0.01 | 0.00 |
| cont | zeta_K | 8630 | 0.99 | 0.01 | 0.00 |
| cont | zeta_V | 8630 | 0.98 | 0.00 | 0.02 |
| ind | Rratio | 43553 | 0.99 | 0.01 | 0.00 |
| ind | pf | 57065 | 0.98 | 0.02 | 0.00 |
| ind | zeta_K | 217765 | 0.99 | 0.01 | 0.00 |
| ind | zeta_V | 217765 | 0.98 | 0.02 | 0.00 |

## 4. The rho* study

**Design (as pre-registered).**
- Panels: 200 per NASDAQ cell for G0 and GT, and for GF, which is used for the P-F slope only.
  Confirmation-sized panels with exploration day lengths.
- Paths: independent stationary days with B = 20,000; R2 clip.
- Grid: ranges [4, 32] and [8, 64]; ceilings 1/80 to 1/5.
- Interval: the widest of the three in section 3.
- Pass rule: |bias| < 0.02, and coverage fails only when coverage + 2 MC SE < 0.85.
- Top-up: every (cell, generator) with an element inside the Monte Carlo band (coverage within 2 MC SE of
  0.85) at some ceiling was topped up to 400 panels. That is 18 pairs: G0 and GT of every cell except
  Q-T1S1, and GF of Q-T3S2 and Q-T3S3.
- Stress test: the continuous-series path at Q-T3S1, G0 and GT, 100 panels.

**Compute.** The study took 115.7 min and the top-up 85.0 min, on 3 workers. Worker seconds per
combination (summed over its chunks, including the reduced run's):

| days | cell | gen | panels | worker seconds |
|---|---|---|---|---|
| ind | Q-T1S1 | G0 | 200 | 1146 |
| ind | Q-T1S1 | GT | 200 | 2535 |
| ind | Q-T1S1 | GF | 200 | 1050 |
| ind | Q-T1S2 | G0 | 400 | 3131 |
| ind | Q-T1S2 | GT | 400 | 6677 |
| ind | Q-T1S2 | GF | 200 | 1950 |
| ind | Q-T1S3 | G0 | 400 | 2281 |
| ind | Q-T1S3 | GT | 400 | 9766 |
| ind | Q-T1S3 | GF | 200 | 1050 |
| ind | Q-T2S1 | G0 | 400 | 984 |
| ind | Q-T2S1 | GT | 400 | 2862 |
| ind | Q-T2S1 | GF | 200 | 524 |
| ind | Q-T2S2 | G0 | 400 | 572 |
| ind | Q-T2S2 | GT | 400 | 1416 |
| ind | Q-T2S2 | GF | 200 | 298 |
| ind | Q-T2S3 | G0 | 400 | 182 |
| ind | Q-T2S3 | GT | 400 | 591 |
| ind | Q-T2S3 | GF | 200 | 99 |
| ind | Q-T3S1 | G0 | 400 | 787 |
| ind | Q-T3S1 | GT | 400 | 2046 |
| ind | Q-T3S1 | GF | 200 | 380 |
| ind | Q-T3S2 | G0 | 400 | 314 |
| ind | Q-T3S2 | GT | 400 | 1024 |
| ind | Q-T3S2 | GF | 400 | 302 |
| ind | Q-T3S3 | G0 | 400 | 304 |
| ind | Q-T3S3 | GT | 400 | 890 |
| ind | Q-T3S3 | GF | 400 | 291 |
| cont | Q-T3S1 | G0 | 100 | 119 |
| cont | Q-T3S1 | GT | 100 | 252 |

**Under the pre-registered rules (Measurement).** rho* was none in every cell, including the stress cell.
Below is the status of each ceiling from `docs/phase3c4_diagnosis.md`, with failing elements by family:
V = zeta_V(a), K = zeta_K(a), PF = P-F slope, R = R_{1/2}/R_0.

| days | cell | 1/80 | 1/40 | 1/20 | 1/10 | 1/5 | ceilings that pass |
|---|---|---|---|---|---|---|---|
| ind | Q-T1S1 | fail: PF 1 | fail: PF 1 | fail: PF 1 | fail: PF 1 | fail: PF 1 | none |
| ind | Q-T1S2 | fail: PF 1 | fail: PF 1 | pass | fail: V 9 | fail: V 14 R 3 | 1/20 |
| ind | Q-T1S3 | n/e (K 10 PF 3) | fail: V 1 K 10 PF 3 R 3 | fail: V 15 PF 1 R 4 | fail: V 18 PF 1 R 4 | fail: V 20 PF 1 R 4 | none |
| ind | Q-T2S1 | fail: K 10 | fail: V 12 | fail: V 17 | fail: V 18 | fail: V 18 | none |
| ind | Q-T2S2 | n/e (K 10 PF 3) | fail: K 10 PF 1 | fail: V 11 PF 1 | fail: V 14 | fail: V 14 | none |
| ind | Q-T2S3 | n/e | n/e | n/e (K 10 PF 3) | fail: PF 1 | fail: PF 2 | none |
| ind | Q-T3S1 | fail: K 3 | pass | fail: PF 1 | pass | pass | 1/40, 1/10, 1/5 |
| ind | Q-T3S2 | n/e (K 10 PF 3) | fail: PF 2 | pass | fail: PF 1 | fail: PF 2 | 1/20 |
| ind | Q-T3S3 | n/e (K 10 PF 3) | fail: PF 2 | fail: PF 2 | fail: PF 2 R 1 | fail: PF 2 R 1 | none |
| cont | Q-T3S1 | fail: K 4 | fail: V 2 | fail: V 5 | fail: V 7 | fail: V 7 | none |

The diagnosis rebuilt each panel's stock counts from its seeds. Day lengths are drawn before any
simulation, and the rebuild matches the saved estimates' finiteness in every case. It found:
- **zeta_K.** All 83 failures had a median of 0 stocks at the range's top horizon. The estimate was a slope
  over the lower horizons only, while the truth is over the whole range.
- **zeta_V.** 181 failures, all with negative bias (at most 0.029 in size), concentrated at the larger
  ceilings. Coverage falls as the number of contributing stocks rises.
- **GF P-F slope.** 36 failures, all on bias alone, with coverage 0.85 to 1.00.

## 5. The rule changes of 31 Oct 2026

Decided on 31 Oct 2026 from `docs/phase3c4_diagnosis.md`. The same provenance is recorded in the code
comments (`pooled.MIN_REACH`, `phase3c4` analysis block, `phase3c4.passes`). Unchanged: the bias bound
0.02, the coverage bound 0.85, the Monte Carlo band rule, the all-ceilings-pass rule and the panels.

1. **Range reach** (`pooled.range_reach`, `pooled.range_rule`, applied in `pooled.core2`, so Phase 4 uses
   it by default):
   - An estimate over [T1, T2] is defined for a panel only if at least `MIN_REACH = 5` of the panel's stocks
     reach T2 under the mask. Otherwise it is NaN, never a slope over a shortened range.
   - **Choice:** zeta_K is also undefined when a horizon inside the range was dropped because its pooled k~
     was not positive. Such a slope is over a shortened range too. Among saved point estimates with the reach
     met, this happens in 145 of 294,580 (panel, a, range, ceiling) cases.
   - The P-F slope needs zeta_K defined at every a <= 1, and vk_diff (descriptive) at every a <= 0.5.
   - **Choice:** R_{1/2}/R_0 over [T1, T2] is measured at the range's middle horizon, but is held to the
     same reach at T2.
   - **Choice:** reach is the panel's, the same in every bootstrap and jackknife replicate.
   - An element is evaluable at a ceiling only if it is defined in at least 90 % of panels (`MIN_EVAL`, was
     50 %).
   - **Choice:** bias, coverage and MC SE are taken over the defined panels. A defined estimate with no
     interval still counts as a miss, as before.
2. **P-F criterion** (`phase3c4.passes(coverage_only=True)` for `pf`): the P-F slope is judged by coverage
   only, under G0, GT and GF. Its bias is still computed and reported (`max_abs_bias` in
   `p3c4_ceilings.tsv`).
3. **rho* per cell and per tick bin.** A tick bin's rho* is the minimum of its cells' rho*: none if any cell
   has none. At each failing ceiling the binding estimator is named: the generator and estimator of the worst
   failing element, scored as in `rho.md`, without the bias term for the P-F slope. A ceiling is skipped, as
   before, when some criterion estimator has no evaluable element there. No evaluable element fails at any
   skipped ceiling, so skipping only ceilings where nothing at all is evaluable gives the same rho*.

**Applying the rule to the saved panels.**
- The chunks were computed with the old estimator. `analyse` applies rule 1 to them through the same
  `pooled.range_rule`, using two inputs:
  - each panel's reach, from its day lengths re-drawn from its seed (`phase3c4.panel_ns`, verified in the
    diagnosis);
  - the saved point estimate's `dropped_K`.
- Masking the old point estimates this way reproduces the new estimator exactly
  (`test_3c4.check_saved_panels_rule`).
- **Not measured:** the bootstrap and jackknife replicates inside the saved intervals also come from the
  old estimator. In the new code, a replicate that drops a k~ horizon inside the range is discarded; in the
  saved ones it was kept as a shortened slope. The raw moments are not saved, so this cannot be redone
  without simulating again. Point estimates drop such a horizon in 0.05 % of cases.
- **Reproducibility:** `pooled.LEGACY_RULE` (`min_reach=0, full_range=False`) is pinned in the 3c-2 and
  3c-3 drivers and in the 3c-4 `study` step. Their committed outputs stay reproducible, and the 3c-4 chunks
  still resume exactly.

## 6. rho* under the new rules

**Status per ceiling** (same layout as the diagnosis). n/e = some criterion estimator has no evaluable
element.

| days | cell | 1/80 | 1/40 | 1/20 | 1/10 | 1/5 | ceilings that pass | rho* |
|---|---|---|---|---|---|---|---|---|
| ind | Q-T1S1 | n/e | pass | pass | pass | pass | 1/40, 1/20, 1/10, 1/5 | 1/5 |
| ind | Q-T1S2 | n/e | pass | pass | fail: V 9 | fail: V 14 R 3 | 1/40, 1/20 | 1/20 |
| ind | Q-T1S3 | n/e | n/e | fail: V 8 R 2 | fail: V 18 R 4 | fail: V 20 R 4 | none | none |
| ind | Q-T2S1 | n/e | fail: V 6 | fail: V 17 | fail: V 18 | fail: V 18 | none | none |
| ind | Q-T2S2 | n/e | n/e | fail: V 5 | fail: V 14 | fail: V 14 | none | none |
| ind | Q-T2S3 | n/e | n/e | n/e | pass | pass | 1/10, 1/5 | 1/5 |
| ind | Q-T3S1 | n/e | pass | pass | pass | pass | 1/40, 1/20, 1/10, 1/5 | 1/5 |
| ind | Q-T3S2 | n/e | n/e | n/e | pass | pass | 1/10, 1/5 | 1/5 |
| ind | Q-T3S3 | n/e | n/e | n/e | pass | fail: R 1 | 1/10 | 1/10 |
| cont | Q-T3S1 | n/e | pass | fail: V 5 | fail: V 7 | fail: V 7 | 1/40 | 1/40 |

**Per tick bin**, with the cells at the tick bin's rho* and what binds them:

| tick bin | cells' rho* | rho* | binding |
|---|---|---|---|
| Q-T1 | Q-T1S1 1/5, Q-T1S2 1/20, Q-T1S3 none | none | Q-T1S3 fails at 1/20: GT Rratio (Rratio[4,32]) |
| Q-T2 | Q-T2S1 none, Q-T2S2 none, Q-T2S3 1/5 | none | Q-T2S1 fails at 1/40: G0 zeta_V (zeta_V(0)[4,32]); Q-T2S2 fails at 1/20: G0 zeta_V (zeta_V(0)[4,32]) |
| Q-T3 | Q-T3S1 1/5, Q-T3S2 1/5, Q-T3S3 1/10 | 1/10 | Q-T3S3 fails at 1/5: GT Rratio (Rratio[8,64]) |

**Per cell**: binding estimator and worst element at each failing ceiling. b = bias, c = coverage +- MC SE,
def = fraction of panels in which it is defined.

| days | cell | panels (G0 / GT / GF) | 1/80 | 1/40 | 1/20 | 1/10 | 1/5 | rho* |
|---|---|---|---|---|---|---|---|---|
| ind | Q-T1S1 | 200 / 200 / 200 | n/e | pass | pass | pass | pass | 1/5 |
| ind | Q-T1S2 | 400 / 400 / 200 | n/e | pass | pass | fail 9: binding GT zeta_V, worst zeta_V(0)[8,64] b -0.007 c 0.74+-0.02 def 1.00 | fail 17: binding GT zeta_V, worst zeta_V(0)[8,64] b -0.008 c 0.58+-0.02 def 1.00 | 1/20 |
| ind | Q-T1S3 | 400 / 400 / 200 | n/e | n/e | fail 10: binding GT Rratio, worst Rratio[4,32] b -0.016 c 0.47+-0.02 def 1.00 | fail 22: binding G0 zeta_V, worst zeta_V(0)[4,32] b -0.016 c 0.19+-0.02 def 1.00 | fail 24: binding GT zeta_V, worst zeta_V(0)[8,64] b -0.029 c 0.03+-0.01 def 1.00 | none |
| ind | Q-T2S1 | 400 / 400 / 200 | n/e | fail 6: binding G0 zeta_V, worst zeta_V(0)[4,32] b -0.013 c 0.75+-0.02 def 1.00 | fail 17: binding G0 zeta_V, worst zeta_V(0.25)[8,64] b -0.024 c 0.61+-0.02 def 1.00 | fail 18: binding G0 zeta_V, worst zeta_V(0)[8,64] b -0.026 c 0.39+-0.02 def 1.00 | fail 18: binding G0 zeta_V, worst zeta_V(0)[8,64] b -0.028 c 0.25+-0.02 def 1.00 | none |
| ind | Q-T2S2 | 400 / 400 / 200 | n/e | n/e | fail 5: binding G0 zeta_V, worst zeta_V(0)[4,32] b -0.015 c 0.75+-0.02 def 1.00 | fail 14: binding G0 zeta_V, worst zeta_V(0)[8,64] b -0.027 c 0.59+-0.02 def 1.00 | fail 14: binding G0 zeta_V, worst zeta_V(0)[8,64] b -0.028 c 0.46+-0.02 def 1.00 | none |
| ind | Q-T2S3 | 400 / 400 / 200 | n/e | n/e | n/e | pass | pass | 1/5 |
| ind | Q-T3S1 | 400 / 400 / 200 | n/e | pass | pass | pass | pass | 1/5 |
| ind | Q-T3S2 | 400 / 400 / 400 | n/e | n/e | n/e | pass | pass | 1/5 |
| ind | Q-T3S3 | 400 / 400 / 400 | n/e | n/e | n/e | pass | fail 1: binding GT Rratio, worst Rratio[8,64] b -0.021 c 0.94+-0.01 def 1.00 | 1/10 |
| cont | Q-T3S1 | 100 / 100 | n/e | pass | fail 5: binding G0 zeta_V, worst zeta_V(0.25)[8,64] b -0.016 c 0.70+-0.05 def 1.00 | fail 7: binding G0 zeta_V, worst zeta_V(0)[8,64] b -0.018 c 0.64+-0.05 def 1.00 | fail 7: binding G0 zeta_V, worst zeta_V(0)[8,64] b -0.018 c 0.65+-0.05 def 1.00 | 1/40 |

**rho\* of each estimator alone** (generator and estimator):

| days | cell | G0 zeta_V | G0 zeta_K | G0 pf | G0 Rratio | GT zeta_V | GT zeta_K | GT pf | GT Rratio | GF pf |
|---|---|---|---|---|---|---|---|---|---|---|
| ind | Q-T1S1 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 |
| ind | Q-T1S2 | 1/20 | 1/5 | 1/5 | 1/10 | 1/20 | 1/5 | 1/5 | 1/10 | 1/5 |
| ind | Q-T1S3 | none | 1/5 | 1/5 | none | none | 1/5 | 1/5 | none | 1/5 |
| ind | Q-T2S1 | none | 1/5 | 1/5 | 1/5 | none | 1/5 | 1/5 | 1/5 | 1/5 |
| ind | Q-T2S2 | none | 1/5 | 1/5 | 1/5 | none | 1/5 | 1/5 | 1/5 | 1/5 |
| ind | Q-T2S3 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 |
| ind | Q-T3S1 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 |
| ind | Q-T3S2 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 |
| ind | Q-T3S3 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/5 | 1/10 | 1/5 |
| cont | Q-T3S1 | 1/40 | 1/5 | 1/5 | 1/5 | 1/20 | 1/5 | 1/5 | 1/5 |

**Fraction of panels with at least 5 stocks at T2** (G0; GT and GF differ from G0 by at most 0.09). This is the
reach part of "defined"; an element is evaluable at 0.90 or more.

| cell | range | 1/80 | 1/40 | 1/20 | 1/10 | 1/5 |
|---|---|---|---|---|---|---|
| Q-T1S1 | [4,32] | 0.47 | 0.96 | 0.98 | 1.00 | 1.00 |
| Q-T1S1 | [8,64] | 0.47 | 0.47 | 0.96 | 0.98 | 1.00 |
| Q-T1S2 | [4,32] | 0.53 | 1.00 | 1.00 | 1.00 | 1.00 |
| Q-T1S2 | [8,64] | 0.10 | 0.53 | 1.00 | 1.00 | 1.00 |
| Q-T1S3 | [4,32] | 0.00 | 0.80 | 1.00 | 1.00 | 1.00 |
| Q-T1S3 | [8,64] | 0.00 | 0.00 | 0.80 | 1.00 | 1.00 |
| Q-T2S1 | [4,32] | 0.83 | 1.00 | 1.00 | 1.00 | 1.00 |
| Q-T2S1 | [8,64] | 0.00 | 0.83 | 1.00 | 1.00 | 1.00 |
| Q-T2S2 | [4,32] | 0.00 | 0.49 | 1.00 | 1.00 | 1.00 |
| Q-T2S2 | [8,64] | 0.00 | 0.00 | 0.49 | 1.00 | 1.00 |
| Q-T2S3 | [4,32] | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 |
| Q-T2S3 | [8,64] | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 |
| Q-T3S1 | [4,32] | 0.81 | 1.00 | 1.00 | 1.00 | 1.00 |
| Q-T3S1 | [8,64] | 0.00 | 0.81 | 1.00 | 1.00 | 1.00 |
| Q-T3S2 | [4,32] | 0.00 | 0.10 | 0.55 | 1.00 | 1.00 |
| Q-T3S2 | [8,64] | 0.00 | 0.00 | 0.10 | 0.55 | 1.00 |
| Q-T3S3 | [4,32] | 0.00 | 0.05 | 0.88 | 1.00 | 1.00 |
| Q-T3S3 | [8,64] | 0.00 | 0.00 | 0.05 | 0.88 | 1.00 |

**Measurement.**
- Independent-days cells: 1,175 evaluable (element, ceiling) pairs, 175 failing: 161 zeta_V(a) and 14
  R_{1/2}/R_0. zeta_K(a) has 470 evaluable pairs and no failure. The P-F slope has 141 evaluable pairs
  and no coverage failure. Its largest |bias| is 0.069 (GF, Q-T3S3, 1/5, coverage 0.98), reported but no
  longer binding.
- Every evaluable element has a defined fraction of at least 0.905.
- 1/80 is not evaluable in any cell.
- **Q-T1 and Q-T2 have no rho\*:**
  - Q-T1S3's first evaluable ceiling is 1/20, where zeta_V and R_{1/2}/R_0 fail on coverage (worst: GT
    R_{1/2}/R_0 [4, 32], bias -0.016, coverage 0.47).
  - Q-T2S1 fails at 1/40 and Q-T2S2 at 1/20, both on zeta_V(a) coverage under G0.
- **Q-T3's rho\* is 1/10,** bound by Q-T3S3, where GT R_{1/2}/R_0 [8, 64] fails on bias (-0.021,
  coverage 0.94) at 1/5.
- **Stress cell:** the continuous-series path gives 1/40, against 1/5 on independent days. zeta_V(a) on
  [8, 64] fails from 1/20 on, with coverage 0.64 to 0.76.

## 7. Reach at each tick bin's rho*

Confirmation stocks (N_min from `data/universe.tsv`) with N_min >= T2 / rho, per cell and pooled per tick
bin, at the tick bin's rho*. Q-T1 and Q-T2 have no rho*. Their rows use the earlier 1/20 fallback, marked
"none found", and are not a reach at rho*.

| group | rho used | confirm stocks | [4,32] | [8,64] |
|---|---|---|---|---|
| Q-T1S1 | 1/20 (none found) | 161 | 34 | 22 |
| Q-T1S2 | 1/20 (none found) | 216 | 52 | 23 |
| Q-T1S3 | 1/20 (none found) | 252 | 26 | 5 |
| **Q-T1** (pooled) | 1/20 (none found) | 629 | 112 | 50 |
| Q-T2S1 | 1/20 (none found) | 80 | 32 | 17 |
| Q-T2S2 | 1/20 (none found) | 50 | 12 | 3 |
| Q-T2S3 | 1/20 (none found) | 17 | 1 | 0 |
| **Q-T2** (pooled) | 1/20 (none found) | 147 | 45 | 20 |
| Q-T3S1 | 1/10 | 50 | 39 | 31 |
| Q-T3S2 | 1/10 | 27 | 23 | 12 |
| Q-T3S3 | 1/10 | 23 | 13 | 5 |
| **Q-T3** (pooled) | 1/10 | 100 | 75 | 48 |

## 8. Tests and planted bugs

`tools/est/test_3c4.py` has 12 checks and 11 planted bugs. 6 checks and 6 bugs are new for the rule
changes. Every bug is caught by its intended assertion; each new mutant was also rerun outside pytest. Full
suite: 344 passed, 28 skipped.

| planted bug (new) | mutated function | caught by |
|---|---|---|
| slope computed from fewer horizons than the range | `pooled.range_rule` | check_full_range_only (zeta_K on [4, 32] finite with T = 8 dropped) |
| reach threshold applied to the lowest T instead of T2 | `pooled.range_reach` | check_reach_at_top (7 stocks counted at T = 4 where 5 reach T2 = 32) |
| P-F slope judged on bias too | `phase3c4.evaluate` | check_pf_coverage_only |
| evaluable at 50 % of panels | `phase3c4.evaluate` | check_evaluable_90 |
| coverage over all panels instead of the defined ones | `phase3c4.evaluate` | check_evaluable_90 |
| tick-bin rho* ignoring a cell with none | `phase3c4.tick_bin_rho` | check_tick_bin_rho |

New checks:
- `check_reach_at_top`: exactly 5 stocks at T2 is defined, 4 is not, whatever the count at T1.
- `check_full_range_only`: a dropped horizon, and a range whose top no stock reaches, are both undefined;
  the old estimator gave a finite shortened slope for both.
- `check_saved_panels_rule`: masking old estimates by reach and `dropped_K` equals the new estimator. This
  check caught a stacking-axis bug in `phase3c4.kept_full_from_dropped` during development. The `analyse`
  run before the fix had used the earlier inline code, and its outputs are byte-identical to the final run.
- `check_pf_coverage_only`, `check_evaluable_90`, `check_tick_bin_rho`.

Three older tests needed `**PO.LEGACY_RULE` because they build 2-stock panels, below the reach of 5:
- `test_3c3.check_v_mean_of_logs`
- `test_pooled2.check_log_slope_equals_K_over_T`
- `test_pooled2.check_nonpositive_pooled_drops_T`, which tests the old horizon-dropping mechanism; its
  docstring now says so.

## 9. Wall time

| job | wall | notes |
|---|---|---|
| cross-day dependence (step 1) | 63 s | guard reads of 8 days, signs and sizes only |
| stationarity check on the new path | 8 s | |
| burn-in check, two passes | 50 s + 53 s | 3 workers |
| reduced study (33 combinations x 10 panels) | 20.4 min | machine swapping |
| full study | 115.7 min | 3 workers |
| top-up to 400 panels (18 pairs) | 85.0 min | 3 workers |
| diagnosis (saved panels) | about 2 min | |
| `analyse` under the new rules | 124 s | includes re-drawing day lengths for every panel |
| test suite | 67 s | |

Peak memory was not recorded for the full study and top-up. The stopped run's peak was 0.91 GB over all
processes.

## 10. Open (decisions not taken here)

- Q-T1 and Q-T2 have no rho* under the new rules. What binds them is coverage of zeta_V(a), whose bias is
  small (-0.004 to -0.029 in the failing elements) while its intervals shrink as more stocks contribute, and of R_{1/2}/R_0 in
  Q-T1S3.
- The effect of the old-rule replicates inside the saved intervals (section 5) is not measured.
