# Phase 3d-2b report: four diagnostics before the 3d-2 study

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken here. Real data were read
only through `explore.load_explore_runs` (exploration stocks of the six confirmatory cells, 8 study days, signs,
sizes and timestamps only, via `worlds3d2.build_flow_cache`). No replication file was opened. Every price
response is synthetic; nothing combines real signs or sizes with real price changes. numpy only, 3 workers.

Code: `tools/est/diag3d2b.py` (driver, resumable, one JSON per panel under
`~/CODE/itch-data/phase3d2/diag3d2b/`, tables in `tables.md` there), `primary3d.tie_bins` and its callers,
`worlds3d2` (clean W-null, W-null-real, W-null coin signs, W0-shuf). The full 3d-2 study was **not** started.

Decisions recorded from your message: the m3_V fix is approved and recorded in `docs/phase3d_tests.md` section
7.1 as a pre-freeze bug fix with its reason; the kappa grid {1, 4, 16} is approved.

**Corrections (1 Nov 2026, after this report).**
- P-R's a* was the least-squares quadratic over the primary a grid {0, 0.25, 0.5, 0.75, 1} (`C.A_GRID`,
  `a_max = 1`), not over [0, 1.5]. The relevant sizes-only vertex is therefore the "vertex on [0, 1]" column
  of section 4 (0.556 to 0.579). It matches the panels' a*(1) (0.57 to 0.58), so the gap that section 4 left
  open does not exist. Option (b) of section 5 was already the estimator in use.
- Nothing is frozen before the end of Phase 3, so the a* change of section 5 is a pre-freeze decision. It was
  taken on 1 Nov 2026 as option (a); see `docs/phase3d_tests.md` 7.3.

## Headlines

1. **Tie-aware bins make psi-hat and P-M usable on real sizes (Measurement; the rule is a Choice).** A psi-hat
   slope is now defined by the bins for 139 of 161 Q-T1 stocks, 6 of 7 Q-T2 and 43 of 43 Q-T3 (old deciles: 3,
   0, 0). P-M usable stock-days rise in every design, most at 5 bins (Q-T3 b5_m50: 0.68 against 0.32).
2. **The clean P-F null holds (Measurement).** With fair-coin signs on the informed runs, P-F rejects in 0 of 20
   panels in each tick bin, and zeta_K(a) is flat. With their real signs (the pilot's W-null, now
   "W-null-real") it rejects 19/20, 2/20, 19/20: the pilot's false rejections came from the real signs.
3. **Partly pre-asymptotic (Measurement).** In the Paper II model, P-V and P-F pass at long horizons (N =
   100000, ranges up to 1024: 30/30 each) but not at ours (P-V 0/30, P-F at most 11/30). P-R and P-M fail at
   both. Even at T = 1024 the measured zeta_V(a) falls with a at about 0.15 per unit a, against Eq. 30's 0.25,
   and a*(T) sits above Eq. 74's band at every T >= 2.
4. **The a* bias is not sign memory (Measurement).** With signs shuffled within the day, a*(1) is still 0.57 to
   0.58 (psi = 0.5). The source is the a* estimator itself: a* is the vertex of a quadratic fitted to R(a)
   over a in [0, 1] (corrected; first written as [0, 1.5]). On real within-day size distributions R(a) peaks
   at psi but falls more slowly to the right than to the left. A sizes-only calculation puts the argmax at
   0.51 and the quadratic vertex at 0.56 to 0.58.
5. **Runtime (Measurement).** The diagnostic run took 114.4 min of wall time on 3 workers (6864 s; peak RSS of
   the parent 0.84 GB, workers about 0.2 GB each), inside the 2 h target. One long panel averaged 456 s. The
   code, tests and the sizes-only step were on top of that.

## 1. Tie-aware size bins (Choice; occupancy is Measurement)

**Rule (Choice, fixed from the pilot's bin occupancy before any power result; `primary3d.tie_bins`).**

- Bins are built on a stock's sorted run sizes. Quantile cuts snap to a boundary between distinct size values,
  so all runs of one size fall in one bin (the nearer boundary; on equal distance, the lower one).
- Coinciding cuts are dropped.
- The smallest bin (the first on ties) is merged with its smaller neighbour (the lower on ties), until every bin
  holds the minimum number of runs or one bin is left. The realised count is recorded.
- psi-hat: 10 target bins, at least PSI_MIN_RUNS = 100 runs per bin, and a slope needs at least
  PSI_MIN_BINS = 3 realised bins (both Choice).
- P-M: a stock-day enters only if its realised count equals the design's target.

Tests: the known answer, plus two planted bugs (a size value split across two bins; an empty bin kept), both
caught (`test_3d2.py`).

**Occupancy on exploration sizes (sizes only, R2-clipped as `make_unit` does).**

The round-lot share of runs (a multiple of 100 shares) is 0.28 in Q-T1, 0.38 in Q-T2 and 0.49 in Q-T3.

psi-hat, per stock (all steps of its 8 days):

| tick bin | stocks | slope defined (>= 3 realised bins) | old rule (all 10 deciles non-empty) | realised bins (bins: stocks) |
|---|---|---|---|---|
| Q-T1 | 161 | 139 | 3 | 1: 14, 2: 8, 3: 4, 4: 6, 5: 6, 6: 3, 7: 8, 8: 35, 9: 49, 10: 28 |
| Q-T2 | 7 | 6 | 0 | 2: 1, 9: 6 |
| Q-T3 | 43 | 43 | 0 | 5: 1, 6: 1, 7: 3, 8: 8, 9: 23, 10: 7 |

Of the 22 Q-T1 stocks without a slope, 18 have fewer than 300 steps over 8 days (58 to 287). The other 4
(306 to 573) lose bins to ties.

In the panels, a stock also loses its slope when a bin's mean signed response is not positive. Stocks with a
slope, mean per panel:

| world | Q-T1 | Q-T2 | Q-T3 |
|---|---|---|---|
| W0 | 237 of 377 | 13 of 17 | 58 of 100 |
| W0-shuf | 304 of 377 | 14 of 17 | 91 of 100 |

P-M, share of stock-days usable:

| design | Q-T1 new (old) | Q-T2 new (old) | Q-T3 new (old) |
|---|---|---|---|
| b3_m20 | 0.79 (0.72) | 0.89 (0.88) | 0.99 (0.92) |
| b3_m50 | 0.64 (0.54) | 0.84 (0.73) | 0.93 (0.80) |
| b3_m100 | 0.46 (0.40) | 0.75 (0.52) | 0.81 (0.68) |
| b4_m50 | 0.52 (0.39) | 0.77 (0.57) | 0.83 (0.63) |
| b5_m20 | 0.61 (0.38) | 0.80 (0.43) | 0.86 (0.44) |
| b5_m50 | 0.40 (0.24) | 0.62 (0.34) | 0.68 (0.32) |
| b5_m100 | 0.26 (0.16) | 0.20 (0.12) | 0.48 (0.21) |

All 27 designs, with the realised-bin counts per stock-day, are in `tables.md`. Most of Q-T1's remaining
losses are short days: 215 of its 1288 stock-days have fewer than 60 runs. Among days with at least 60 runs,
only 55 realise fewer than 3 bins at b3_m20.

**psi-hat bias (Measurement).** psi-hat is 0.553 (MC SE 0.001), 0.532 (0.004) and 0.521 (0.002) in W0-shuf
(true psi 0.5). In W0 with real signs it is 0.666, 0.618 and 0.705. So sign memory adds about 0.1 to 0.2, and
a residual of +0.02 to +0.05 remains without it. Whether that residual is a Jensen effect of the wide merged
bins (log of a bin mean against the log of its median size) was not isolated (open).

## 2. The clean P-F null (Measurement)

In W-null, the informed runs (a random GF_F_INF share, as in GF of Phase 3a) get independent fair-coin signs.
The rest of the flow stays real. W-null-real is the pilot's W-null, with the informed runs keeping their real
signs. Both use 20 panels per tick bin and P-F alone. A rejection here means p <= 0.05 unadjusted with every
condition (the most lenient Holm step).

| world | tick bin | P-F rejects, k/n = rate (MC SE) | slope of zeta_K on a (MC SE) | zeta_K(a), a = 0, 0.25, 0.5, 0.75, 1 |
|---|---|---|---|---|
| W-null | Q-T1 | 0/20 = 0.00 (0.00) | 0.011 (0.013) | 0.981, 0.987, 0.991, 0.993, 0.993 |
| W-null | Q-T2 | 0/20 = 0.00 (0.00) | 0.083 (0.042) | 0.929, 0.913, 0.959, 0.984, 1.005 |
| W-null | Q-T3 | 0/20 = 0.00 (0.00) | -0.001 (0.009) | 0.977, 0.979, 0.979, 0.978, 0.977 |
| W-null-real | Q-T1 | 19/20 = 0.95 (0.05) | -0.041 (0.002) | 1.312, 1.298, 1.287, 1.278, 1.270 |
| W-null-real | Q-T2 | 2/20 = 0.10 (0.07) | -0.039 (0.013) | 1.246, 1.232, 1.223, 1.215, 1.205 |
| W-null-real | Q-T3 | 19/20 = 0.95 (0.05) | -0.083 (0.004) | 1.333, 1.308, 1.285, 1.267, 1.250 |

zeta_K is flat in W-null within 2 MC SE in each tick bin (Q-T2's 0.083 is 2.0 SE, on 17 stocks per panel). A
0 of 20 rate is consistent with a size of 0.05 (probability 0.36). Its MC SE shows 0.00 because k = 0; with 20
panels the rate is only bounded below about 0.14. With their real signs, the informed runs (permanent impact)
carry the real flow's sign memory, so zeta_K sits near 1.3 and falls with a. That is the pilot's 8 in 10.

## 3. Is the paper's world pre-asymptotic? (Measurement)

Setup: the Paper II model (`mbmodel.PAPER2`: mu_m = 1.5, lambda = 1/8, sigma_l = 1, beta_m = 0.25,
lambda' = 0.25, beta_q clipped to [0, 0.99]), with 30 panels in each setting.

- **(a) Ours.** Day lengths from the exploration stocks of each cell, as in the pilot. The ranges are [4,32]
  and [8,64].
- **(b) Long.** 50 stocks, 8 days, N = 100000 runs per day, rho = 0.2 (Choice: the ceiling, and any
  rho >= 0.0103 reaches 1024). The ranges are [64,512] and [128,1024], with T up to 1024.

Long panels use the same `primary3d.run_unit`, with only the T grid and the ranges patched. A test checks that
p-values are identical to the unpatched path and statistics agree to 1e-9. Intervals below are the panel means
of each panel's interval ends.

**Paper predictions (PAPER2).**

- Eq. 30: zeta_V(a) = 3 - mu_m - 2 a lambda sigma_l^2 = 1.50, 1.44, 1.38, 1.31, 1.25 at a = 0, 0.25, 0.5, 0.75,
  1 (a_c = 2, so all five points lie below it).
- Eq. 29, the finite-T form:
  - T = 64: 1.53, 1.47, 1.41, 1.35, 1.28.
  - T = 1024: 1.55, 1.49, 1.43, 1.37, 1.30.
- Eq. 74: a*(T) lies between 1/2 and lambda log T, so the upper edge is 0.52 at T = 64, 0.78 at T = 512 and
  0.87 at T = 1024. B/A grows like T^0.21.
- Generator: beta_q = 0.25 - 0.25 (log q - 3), clipped at 0.

### P-V

| setting | range | M1 err | M2 err | M3 err | M2 < M1 | M2 < M3 | real - S1 slope [lo, hi] | real - S2 slope [lo, hi] | P-V passes |
|---|---|---|---|---|---|---|---|---|---|
| ours Q-T1 | [4,32] | 0.0456 | 0.0122 | 0.0112 | 30/30 | 9/30 | +0.071 [0.052, 0.078] | +0.066 [0.052, 0.073] | 0/30 |
| ours Q-T1 | [8,64] | 0.0275 | 0.0081 | 0.0029 | 30/30 | 0/30 | +0.015 [-0.008, 0.026] | +0.031 [0.012, 0.040] | 0/30 |
| ours Q-T2 | [4,32] | 0.0633 | 0.0293 | 0.0177 | 30/30 | 3/30 | +0.073 [0.023, 0.113] | +0.069 [0.029, 0.103] | 0/30 |
| ours Q-T2 | [8,64] | 0.0495 | 0.0288 | 0.0085 | 30/30 | 0/30 | +0.015 [-0.045, 0.063] | +0.029 [-0.022, 0.069] | 0/30 |
| ours Q-T3 | [4,32] | 0.0428 | 0.0111 | 0.0144 | 30/30 | 28/30 | +0.074 [0.052, 0.086] | +0.067 [0.050, 0.078] | 0/30 |
| ours Q-T3 | [8,64] | 0.0294 | 0.0101 | 0.0034 | 30/30 | 0/30 | +0.017 [-0.011, 0.032] | +0.031 [0.009, 0.044] | 0/30 |
| long | [64,512] | 0.0027 | 0.0003 | 0.0035 | 30/30 | 30/30 | -0.098 [-0.117, -0.084] | -0.060 [-0.077, -0.048] | 30/30 |
| long | [128,1024] | 0.0016 | 0.0002 | 0.0026 | 30/30 | 30/30 | -0.119 [-0.143, -0.102] | -0.080 [-0.102, -0.065] | 30/30 |

Measured zeta_V(a) at a = 0, 0.25, 0.5, 0.75, 1:

| setting | measured zeta_V(a) | paper |
|---|---|---|
| ours, [4,32] | 1.62, 1.59, 1.56, 1.52, 1.49 | Eq. 30: 1.50 to 1.25 |
| ours, [8,64] | 1.63, 1.60, 1.57, 1.53, 1.48 | Eq. 30: 1.50 to 1.25 |
| long, [64,512] | 1.585, 1.549, 1.511, 1.473, 1.434 | Eq. 29 at T = 512: 1.55 to 1.30 |
| long, [128,1024] | 1.566, 1.528, 1.490, 1.450, 1.409 | Eq. 29 at T = 1024: 1.55 to 1.30 |

Our rows were the same in every tick bin to 0.01.

M2 fitted on all days gives:

- ours: mu_m 1.37 to 1.39 and lambda sigma^2 0.06 to 0.07.
- long: mu_m 1.41 to 1.43 and lambda sigma^2 0.074 to 0.078.

The preset's values are 1.5 and 0.125.

At long horizons the level at a = 0 is close to Eq. 29 (1.566 against 1.554 at T = 1024). The slope in a is
-0.15 to -0.16, against the paper's -0.25, and it has not moved toward -0.25 between [64,512] and [128,1024].

**P-V: passes at long horizons but not at ours.** At ours, M2 beats M1 always but loses to M3 at [8,64]. The
real-minus-surrogate slopes are positive rather than negative.

### P-R

a*(T), all stocks and days pooled:

| setting | T = 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 | 256 | 512 | 1024 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Eq. 74 upper edge | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.50 | 0.52 | 0.61 | 0.69 | 0.78 | 0.87 |
| ours Q-T1 | 0.481 | 0.541 | 0.607 | 0.673 | 0.735 | 0.796 | 0.862 | | | | |
| ours Q-T2 | 0.480 | 0.541 | 0.608 | 0.675 | 0.743 | 0.814 | 0.890 | | | | |
| ours Q-T3 | 0.475 | 0.537 | 0.604 | 0.670 | 0.734 | 0.800 | 0.870 | | | | |
| long | 0.451 | 0.516 | 0.583 | 0.643 | 0.697 | 0.747 | 0.800 | 0.858 | 0.928 | 1.014 | 1.122 |

MC SEs are at most 0.006. P-R needs a*(T1) and a*(T2) both within 0.5 +- 0.15.

| setting | range | a*(T1), a*(T2) | P-R passes |
|---|---|---|---|
| ours | [4,32] | 0.61, 0.80 | 0/30 |
| ours | [8,64] | 0.67, 0.86 to 0.89 | 0/30 |
| long | [64,512] | 0.80, 1.01 | 0/30 |
| long | [128,1024] | 0.86, 1.12 | 0/30 |

a* is above Eq. 74's band at every T >= 2 in both settings. On the paper's band, P-R would pass at our
horizons (upper edge 0.52 at T = 64).

**P-R: fails at both.** Caveat: section 4 shows the a* estimator reads high on real size distributions.
PAPER2 sizes are lognormal (log q ~ N(3, 1), redrawn below q = 1). For a lognormal, the curve behind a* is
symmetric and the quadratic vertex is nearly unbiased (0.500 to 0.507, section 4), so that bias should be
small here.
It was not computed on the model's own flow (open).

### P-F

P-F rejects flatness:

| setting | range | Q-T1 | Q-T2 | Q-T3 |
|---|---|---|---|---|
| ours | [4,32] | 1/30 | 0/30 | 0/30 |
| ours | [8,64] | 11/30 | 0/30 | 7/30 |
| long | [64,512] | 30/30 | | |
| long | [128,1024] | 30/30 | | |

**P-F: passes at long horizons but not at ours** (power at most 0.37 at ours).

### P-M

The tie-aware bins use 3 bins, at least 50 runs per bin, and lag windows [16,128] and [32,256].

| setting | design | mean log q by bin | beta-hat by bin | generator beta_q at those log q | slope [lo, hi] | slope upper < 0 | top-bin upper < 0 | P-M passes |
|---|---|---|---|---|---|---|---|---|
| ours Q-T1 | l16 L128 | 1.73, 2.80, 3.91 | 0.485, 0.349, 0.243 | 0.57, 0.30, 0.02 | -0.111 [-0.140, -0.080] | 30/30 | 0/30 | 0/30 |
| ours Q-T1 | l32 L256 | 1.73, 2.80, 3.92 | 0.468, 0.347, 0.217 | 0.57, 0.30, 0.02 | -0.115 [-0.138, -0.093] | 30/30 | 0/30 | 0/30 |
| ours Q-T3 | l16 L128 | 1.74, 2.80, 3.91 | 0.473, 0.350, 0.241 | 0.56, 0.30, 0.02 | -0.107 [-0.137, -0.078] | 30/30 | 0/30 | 0/30 |
| ours Q-T3 | l32 L256 | 1.73, 2.80, 3.91 | 0.461, 0.348, 0.217 | 0.57, 0.30, 0.02 | -0.112 [-0.139, -0.088] | 30/30 | 0/30 | 0/30 |
| long | l16 L128 | 1.70, 2.80, 3.92 | 0.495, 0.360, 0.250 | 0.57, 0.30, 0.02 | -0.110 [-0.119, -0.102] | 30/30 | 0/30 | 0/30 |
| long | l32 L256 | 1.70, 2.80, 3.92 | 0.477, 0.353, 0.223 | 0.57, 0.30, 0.02 | -0.115 [-0.123, -0.105] | 30/30 | 0/30 | 0/30 |

Q-T2 was not evaluable in any panel: 0 stocks reach L = 128 under their ceilings, and 5 are needed.

The estimated slope of beta on log q is -0.11, against the generator's -0.25. It is attenuated by more than
half, and N = 100000 does not reduce the attenuation. The top-bin condition (beta below 0) cannot hold in
PAPER2, because the generator clips beta_q at 0.

**P-M: fails at both.** The slope condition holds at both. The failure is the top-bin condition, which is a
property of the preset, not of the horizon.

### Summary of item 3

| test | ours | long | verdict |
|---|---|---|---|
| P-V | 0/30 | 30/30 | passes at long horizons but not at ours |
| P-R | 0/30 | 0/30 | at neither |
| P-F | 0 to 11/30 | 30/30 | passes at long horizons but not at ours |
| P-M | 0/30 | 0/30 | at neither (top-bin condition impossible under PAPER2) |

## 4. The a* bias source (Measurement)

W0 has psi = 0.5, with 30 panels per tick bin. a*(T) at every T, with MC SEs of at most 0.005 except Q-T2 at
T >= 32:

| world | tick bin | T = 1 | 2 | 4 | 8 | 16 | 32 | 64 | P-R passes |
|---|---|---|---|---|---|---|---|---|---|
| W0 (real signs) | Q-T1 | 0.631 | 0.622 | 0.613 | 0.600 | 0.595 | 0.592 | 0.585 | 14/30 = 0.47 (0.09) |
| W0 (real signs) | Q-T2 | 0.613 | 0.610 | 0.608 | 0.592 | 0.593 | 0.596 | 0.575 | 0/30 = 0.00 (0.00) |
| W0 (real signs) | Q-T3 | 0.614 | 0.610 | 0.605 | 0.602 | 0.596 | 0.582 | 0.553 | 24/30 = 0.80 (0.07) |
| W0-shuf (no sign memory) | Q-T1 | 0.583 | 0.588 | 0.596 | 0.604 | 0.619 | 0.642 | 0.660 | 0/30 = 0.00 (0.00) |
| W0-shuf (no sign memory) | Q-T2 | 0.581 | 0.589 | 0.598 | 0.597 | 0.616 | 0.635 | 0.644 | 0/30 = 0.00 (0.00) |
| W0-shuf (no sign memory) | Q-T3 | 0.570 | 0.576 | 0.581 | 0.594 | 0.611 | 0.628 | 0.634 | 0/30 = 0.00 (0.00) |

**Answer: removing sign memory does not bring a* to 0.5.** It lowers a*(1) by about 0.04 to 0.05. a*(T) then
rises with T, and P-R fails everywhere, because a*(T2) leaves the margin.

**Where the bias comes from (`diag3d2b.py vertex`, sizes only).** Without sign memory, the expected T = 1
correlation, given the sizes, is the following for each stock-day:

R(a) proportional to mean(qc^a q^psi) / sqrt(mean(qc^2a) mean(q^2psi))

Here q is the raw size the response uses and qc is the R2-clipped size the estimator uses. By Cauchy-Schwarz
this peaks at a = psi on every day. Pooled as `core2` pools it (Fisher z, stock means), on the exploration sizes:

| tick bin | R(a) / R(0.5) at a = 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5 | argmax (fine grid) | quadratic vertex on [0, 1.5] | vertex on [0, 1] (= a* until 1 Nov 2026) |
|---|---|---|---|---|
| Q-T1 | 0.857, 0.972, 1, 0.983, 0.935, 0.869, 0.797 | 0.51 | 0.635 | 0.579 |
| Q-T2 | 0.839, 0.969, 1, 0.979, 0.921, 0.842, 0.756 | 0.51 | 0.623 | 0.572 |
| Q-T3 | 0.839, 0.964, 1, 0.974, 0.905, 0.816, 0.725 | 0.51 | 0.597 | 0.556 |

The curve falls more slowly to the right of 0.5 than to the left. The quadratic fitted over the grid
therefore puts its vertex to the right of the true maximum: 0.05 to 0.07 on P-R's grid [0, 1] (the a* in use;
corrected, first written as [0, 1.5]), 0.09 to 0.13 on [0, 1.5].
Within-day size distributions drive this. On the same Q-T1 sizes without the clip, pooled over each stock's 8
days rather than day by day, the vertex is 0.53. For lognormal sizes the curve is symmetric: the analytic curve
for sigma of log q from 0.8 to 1.6 gives a vertex of 0.500 to 0.507.

Corrected: on P-R's grid the sizes-only vertex is 0.556 to 0.579, against the panels' a*(1) of 0.57 to 0.58.
As first written, this paragraph compared the [0, 1.5] vertex (0.60 to 0.64) and reported a gap of 0.03 to
0.05 that does not exist.

A scratch check (not committed) generated W0's response on clipped sizes too, with 5 panels each in Q-T1 and
Q-T3. a*(1) moved by 0.005, so the raw/clipped mismatch is not the source.

**Sign memory then shapes a*(T) (Measurement).** With real signs, a*(1) is 0.04 to 0.05 higher and a* falls
with T (to 0.55 to 0.59 at T = 64). Without sign memory, a* rises with T (to 0.63 to 0.66). W0's P-R passes
(0.47, 0.00, 0.80) therefore depend on two effects of opposite sign at T2: the vertex bias and the sign-memory
drift.

## 5. Open, for your decision (none acted on)

1. **The a* estimator (`pooled.core2`).** On real sizes, a* reads 0.06 to 0.13 above psi through the
   quadratic vertex alone (corrected: 0.05 to 0.07 on P-R's grid [0, 1]). Options:
   - (a) A local estimate: a quadratic through the grid maximum and its two neighbours, or a finer a grid near
     the maximum.
   - (b) Fit over [0, 1] instead of [0, 1.5]. This roughly halves the bias, to 0.56 to 0.58. (Corrected:
     this was already the estimator in use.)
   - (c) Keep a* and centre P-R's margin on the W0 calibration instead of on 1/2.
   
   (a) is closest to "the argmax of R(a)", which is what a* is meant to be. (Decided 1 Nov 2026: (a), as a
   pre-freeze decision, `docs/phase3d_tests.md` 7.3.)
2. **P-M under PAPER2.** The top-bin condition cannot hold, because beta_q is clipped at 0. Either WM is not a
   P-M power world, or PAPER2 needs a beta clip that allows negative beta_q. Separately, beta-hat's slope is
   attenuated to about 45 % of the generator's at both horizons.
3. **psi-hat's residual +0.02 to +0.05 without sign memory.** A Jensen effect of the merged bins is a
   candidate. It was not isolated.
4. **The 3d-2 pilot numbers predate the redefined W-null and the tie-aware bins.** `docs/phase3d2_report.md`
   now says so. The full study was not started.

## 6. Reproduce

```
cd ~/CODE/itchbook/tools/est
uv run --group analysis diag3d2b.py all --workers 3     # occupancy, vertex, then 420 panel jobs (resumable)
uv run --group analysis diag3d2b.py summarise           # tables.md
```

Tests: `test_3d2.py`, with checks for tie_bins, coin_informed and the long unit's equality with run_unit, and
planted bugs for:

- a size split across bins
- an empty bin kept
- P-M using a stock-day below the target count
- W-null informed runs keeping their real signs
- a W-null spec without coin signs
- long-unit surrogates of the wrong kind
- the grid patch not restored
