# Phase 3c-4 report: independent stationary days, the widest interval, and why rho* was not run

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
   1.3 SE of the 10^6-run truth. zeta_V's SE falls back to 0.005 to 0.011, against 0.05 to 0.11 with
   3c-3's continuous days.
3. **The burn-in check keeps B = 20,000 (Measurement and Choice).** B = 1000, 2000 and 5000 all fail on
   zeta_K(0), which sits 2.4 to 6.8 SE above the truth on [8, 64]. zeta_V passes at every B.
4. **The full rho* study was stopped before its 10 % mark (Measurement).** After 28 min it projected about
   5.5 h on 3 workers, over the 3 h stop rule. I stopped at 9 new chunks, short of the 10 % mark (53
   chunks), because the projection does not depend on the first chunks alone: it uses the reduced run's
   measured cost for all 33 combinations (373 min), scaled by the observed speed ratio (0.90). GT accounts
   for 61 % of the cost. No rho* was determined. The stress test was not run, and the reach table is at the
   fallback 1/20.
5. **Which interval is widest (Measurement, partial).** The day-and-stock bootstrap was widest in 95 % to
   98 % of element-panels, the day jackknife in 2 % to 4 %, and the whole-stock bootstrap essentially
   never. This comes from the 39 completed chunks, mostly the largest cells.

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

A method undefined for an element is skipped for it.

Which was widest, over the 39 completed chunks (mostly Q-T1S3, Q-T1S2 and the reduced chunks of every
cell; partial):

| estimator | elements x panels with an interval | day_stock_bootstrap | day_jackknife_t | whole_stock_bootstrap |
|---|---|---|---|---|
| zeta_V | 10360 | 0.96 | 0.04 | 0.00 |
| zeta_K | 12560 | 0.96 | 0.04 | 0.00 |
| pf | 3321 | 0.95 | 0.04 | 0.00 |
| Rratio | 2736 | 0.98 | 0.02 | 0.00 |

## 4. The rho* study (stopped)

Planned as specified: 200 panels per NASDAQ cell for G0 and GT, and GF for the P-F slope only;
confirmation-sized panels with exploration day lengths; independent stationary days with B = 20,000; R2
clip; ranges [4, 32] and [8, 64]; ceilings 1/80 to 1/5; the widest-of-three interval.

The pass rule is |bias| < 0.02, and coverage fails only when coverage + 2 MC SE < 0.85. rho* is the
largest ceiling such that every evaluable ceiling up to it passes. Bins inside the MC band would have been
topped up to 400 panels.

The reduced run (10 panels for each of the 33 combinations; machine swapping at the time) measured the
cost per combination:

| combination | stock-days per panel | seconds per 10-panel chunk (reduced run) | projected minutes for 20 chunks on 1 worker |
|---|---|---|---|
| cont_Q-T3S1_G0 | 400 | 23 | 4 |
| cont_Q-T3S1_GT | 400 | 64 | 11 |
| ind_Q-T1S1_G0 | 1288 | 104 | 35 |
| ind_Q-T1S1_GF | 1288 | 79 | 26 |
| ind_Q-T1S1_GT | 1288 | 413 | 138 |
| ind_Q-T1S2_G0 | 1728 | 210 | 70 |
| ind_Q-T1S2_GF | 1728 | 222 | 74 |
| ind_Q-T1S2_GT | 1728 | 527 | 176 |
| ind_Q-T1S3_G0 | 2016 | 105 | 35 |
| ind_Q-T1S3_GF | 2016 | 108 | 36 |
| ind_Q-T1S3_GT | 2016 | 477 | 159 |
| ind_Q-T2S1_G0 | 640 | 74 | 25 |
| ind_Q-T2S1_GF | 640 | 77 | 26 |
| ind_Q-T2S1_GT | 640 | 163 | 54 |
| ind_Q-T2S2_G0 | 400 | 56 | 19 |
| ind_Q-T2S2_GF | 400 | 56 | 19 |
| ind_Q-T2S2_GT | 400 | 63 | 21 |
| ind_Q-T2S3_G0 | 136 | 9 | 3 |
| ind_Q-T2S3_GF | 136 | 8 | 3 |
| ind_Q-T2S3_GT | 136 | 57 | 19 |
| ind_Q-T3S1_G0 | 400 | 46 | 15 |
| ind_Q-T3S1_GF | 400 | 43 | 14 |
| ind_Q-T3S1_GT | 400 | 107 | 36 |
| ind_Q-T3S2_G0 | 216 | 21 | 7 |
| ind_Q-T3S2_GF | 216 | 30 | 10 |
| ind_Q-T3S2_GT | 216 | 84 | 28 |
| ind_Q-T3S3_G0 | 184 | 36 | 12 |
| ind_Q-T3S3_GF | 184 | 32 | 11 |
| ind_Q-T3S3_GT | 184 | 110 | 37 |
| total | | | 1120 worker-minutes, 373 min on 3 workers |

**Measurement.**
- The full run's first 28 min confirmed these costs: Q-T1S3 GT chunks took 315 to 509 s, against 477 s in
  the reduced run, a speed ratio of 0.90.
- The projection was **about 5.5 h** on 3 workers before top-ups, so the run was stopped under the 3 h
  rule.
- GT costs 3 times G0 (its 12 size bins need 12 forward FFTs per series). Every stock-day carries its own
  20,000-run burn-in, so a panel generates 20,000 x stock-days extra runs: about 40 million for Q-T1S3.
- In the pool a Q-T1S3 GT panel took about 48 s, against 15 s alone in one process: 3 parallel FFT-heavy
  workers on a fanless 8 GB machine that was swapping (11.7 of 12.3 GB swap used at the stop).

All 39 completed chunks are kept under `$ITCH_DATA/phase3c4/study/` and resume exactly
(`phase3c4.py study`). **The continuous-series stress test (Q-T3S1, G0 and GT, 100 panels) was not run.**
Its reduced chunks exist.

## 5. Reach (no rho* determined: fallback 1/20)

Confirmation stocks (N_min from `data/universe.tsv`) with N_min >= T2 / rho, per cell and pooled per tick
bin. In this run every cell uses 1/20, so the "each at its cell's rho" in the pooled rows means 1/20 here.

| group | rho used | confirm stocks | [4,32] | [8,64] |
|---|---|---|---|---|
| Q-T1S1 | 1/20 (none found) | 161 | 34 | 22 |
| Q-T1S2 | 1/20 (none found) | 216 | 52 | 23 |
| Q-T1S3 | 1/20 (none found) | 252 | 26 | 5 |
| Q-T2S1 | 1/20 (none found) | 80 | 32 | 17 |
| Q-T2S2 | 1/20 (none found) | 50 | 12 | 3 |
| Q-T2S3 | 1/20 (none found) | 17 | 1 | 0 |
| Q-T3S1 | 1/20 (none found) | 50 | 31 | 21 |
| Q-T3S2 | 1/20 (none found) | 27 | 12 | 6 |
| Q-T3S3 | 1/20 (none found) | 23 | 5 | 1 |
| Q-T1 (S1+S2+S3, each at its cell's rho) | per cell | 629 | 112 | 50 |
| Q-T2 (S1+S2+S3, each at its cell's rho) | per cell | 147 | 45 | 20 |
| Q-T3 (S1+S2+S3, each at its cell's rho) | per cell | 100 | 48 | 28 |

## 6. Tests and planted bugs

`tools/est/test_3c4.py` has 6 checks and 5 planted bugs, all caught, each by its intended assertion (each
mutant rerun outside pytest). Full suite: 332 passed, 28 skipped.

| planted bug | mutated function | caught by |
|---|---|---|
| days of a stock drawn from one series on the new path | `panel3c.make_panel` | check_days_independent (day-to-day correlation of mean signs 0.97 instead of about 0) |
| burn-in shared across a stock's days | `panel3c.independent_series` | check_burn_per_day |
| whole-stock bootstrap resampling days instead of stocks | `pooled.whole_stock_weights` | check_whole_stock_bootstrap |
| rho* taken as the largest passing ceiling | `phase3c4.rho_star` | check_rho_star_all_pass |
| Monte Carlo band ignored | `phase3c4.passes` | check_pass_rule_mc_band |

Two earlier tests needed changes:
- The 3c-3 test of continuous days now asks for `continuous=True` explicitly.
- The 3c planted bug "batched kernel off by one lag" now mutates `panel3c.kernel_spectrum`, where the
  kernel now lives.

## 7. Wall time and peak memory

| job | wall | peak RAM | notes |
|---|---|---|---|
| cross-day dependence (step 1) | 63 s | | guard reads of 8 days, signs and sizes only |
| stationarity check on the new path | 8 s | 0.22 GB | |
| burn-in check, two passes (the second with N = 3000 doubled for SE <= 0.01) | 50 s + 53 s | | 3 workers |
| reduced study (33 combinations x 10 panels) | 20.4 min | | machine swapping |
| full study, stopped | 28 min | 0.91 GB all processes | 3 workers, 10 new chunks |
| test suite | 55 s | | |

Total compute in this phase is about 55 min. The study stopped under the time rule; nothing ran over.

## 8. To finish the study (a decision, not taken here)

- Run the remaining chunks unattended on a quiet machine (about 5.5 h, or less without swapping). They
  resume from the checkpoints.
- Or cut GT's cost: fewer size bins, or a Chebyshev interpolation of its kernel in beta as in
  `mbmodel.impact_trade_clock` (approximate, about 5 nodes instead of 12 FFTs).
- Or lower the panel count for GT only.

A burn-in below 20,000 is ruled out by the check above.
