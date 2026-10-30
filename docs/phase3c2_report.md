# Phase 3c-2 report: the k-tilde normalisation, the t jackknife, and rho* again (NASDAQ cells)

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken here. No real data were
read in this step: every panel is synthetic, and the truth is the Phase 3c truth (3b for GT), recomputed on
the new ranges and a grid. The decisions fixed since 3c are implemented as given:
- P-K is descriptive, and the slope difference is computed but not in the criterion;
- the a grid is {0, 0.25, 0.5, 0.75, 1};
- [64, 1024] is dropped; the ranges are [4, 32], [8, 64], [16, 256].

New code: `pooled.py` (norm="tilde" path: `stock_curves_tilde`, `chain_zeta`, `core2`, `estimate2`,
`intervals2`, `jackknife2`, `jk_interval`; the 3c path kept as norm="old"), `phase3c2.py` (driver),
`test_pooled2.py`. Tables: `data/phase3/p3c2_bins.tsv`, `p3c2_rho.tsv`, `p3c2_reach.tsv`,
`p3c2_finite_truth.tsv`. Per-panel results are in `$ITCH_DATA/phase3c2/`.

## Headlines

1. **Still no rho* in any cell (Measurement).** No ceiling from T/N <= 1/160 to <= 1/10 gives |bias| < 0.02
   and coverage >= 0.85 for every criterion estimator under all three generators. Per estimator it is much
   better than in 3c:
   - zeta_V passes up to 1/10 in 4 of 9 cells for G0 and GT;
   - the P-F slope passes in 7 of 9 for G0 and GT;
   - the binding estimators are G0 and GT zeta_K and GF's R_{1/2}/R_0.
2. **G0's zeta_K target depends on day length (Measurement).** The 3c truth uses days of 10^6 runs, but
   K_a(T) depends on N itself: each day starts with no past impact, and with long sign memory and a slowly
   decaying kernel the reversal from the day's history keeps growing with N.
   - G0's population zeta_K(0) on [8, 64] in Q-T1S3 is 0.948 at N = 300, 0.889 at 1000, 0.855 at 3000,
     0.836 at 10^4 and 0.812 at 10^6. On [16, 256] it is 1.167 at N = 300 against 0.762.
   - This is a property of the target, not of the estimator, and no T/N ceiling removes it. It accounts
     for most of the zeta_K "bias" against the 10^6 truth.
3. **The new pooling (Choice) biases zeta_V where few stocks contribute (Measurement).** On the same panels
   G0's zeta_V at ceiling 1/40 in Q-T3S1 has bias 0.103 and coverage 0.02 with the new path, against 0.005
   and 0.94 with the 3c path.
   - At a = 0 the two normalisations coincide (V_0(1) = 1), so the cause is the pooling. To handle a
     non-positive pooled k~ as specified, I pooled both v and k~ as logs of matched means across stocks;
     3c took the mean of each stock's own log.
   - With only a few heavy-tailed stocks above the ceiling at the top of the range, the ratio of means is
     strongly biased (an independent check on 30 fresh G0 panels gives -0.107 on [16, 256] at ceiling 1/40).
   - The normalisation fix itself does what it should for zeta_K: same bias as 3c, with no stock-day
     excluded.
4. **The interval rule was already right (Measurement).** `inference.jackknife_days` and
   `pooled.intervals` both already used Student t quantiles with G - 1 df. Nothing was changed except
   factoring the quantile into `jk_interval`, which is now tested with a planted normal-quantile bug.
5. **Reach is the hard limit (Measurement).** At ceiling 1/40 (no rho* found) at most 23 confirmation
   stocks per cell reach [4, 32], at most 6 reach [8, 64], and at most 2 reach [16, 256]. Every one of
   them is far above N_floor = 800, so the floor does not bind; the ceiling does.

## 1. The normalisation fix

- k~_a(T) = Kbar_a(T) / (T sqrt(Sbar(1) Vbar_a(1))) and v_a(T) = Vbar_a(T) / (T Vbar_a(1)), with Xbar the
  stock's day average under the replicate's day weights.
- No stock-day is excluded for the sign of K_a(1). Days are averaged within a stock, stocks get equal
  weight, sets are matched, and logs are chained.
- **Choice:** both v and k~ are pooled as logs of matched-set means (logs of curves already normalised by
  each stock's own constants, never of raw averages). This is what makes "a pooled k~ that is not
  positive" well defined: such a horizon is dropped for that estimate and the chain joins the neighbouring
  kept horizons with a matched set. A range needs at least 3 kept horizons.
- The 3c path is kept as `norm="old"`.

**Measurement (dropped horizons):**
- G0: no panel drops a horizon.
- GT: 2 of 4400 panels.
- GF: 2 % to 26 % of panels drop at least one horizon at ceiling 1/10. GF's K is carried by its 10 %
  informed runs alone, so the pooled k~ at large T is often near 0.

| sign model | cell | gen | panels with a dropped T (ceiling 1/10) | panels with a dropped T at any ceiling |
|---|---|---|---|---|
| main | Q-T1S1 | G0 | 0 | 0 |
| main | Q-T1S1 | GT | 0 | 0 |
| main | Q-T1S1 | GF | 87 | 119 |
| main | Q-T1S2 | G0 | 0 | 0 |
| main | Q-T1S2 | GT | 0 | 0 |
| main | Q-T1S2 | GF | 86 | 129 |
| main | Q-T1S3 | G0 | 0 | 0 |
| main | Q-T1S3 | GT | 0 | 0 |
| main | Q-T1S3 | GF | 73 | 74 |
| main | Q-T2S1 | G0 | 0 | 0 |
| main | Q-T2S1 | GT | 0 | 0 |
| main | Q-T2S1 | GF | 55 | 55 |
| main | Q-T2S2 | G0 | 0 | 0 |
| main | Q-T2S2 | GT | 1 | 1 |
| main | Q-T2S2 | GF | 63 | 63 |
| main | Q-T2S3 | G0 | 0 | 0 |
| main | Q-T2S3 | GT | 0 | 0 |
| main | Q-T2S3 | GF | 9 | 9 |
| main | Q-T3S1 | G0 | 0 | 0 |
| main | Q-T3S1 | GT | 0 | 0 |
| main | Q-T3S1 | GF | 53 | 53 |
| main | Q-T3S2 | G0 | 0 | 0 |
| main | Q-T3S2 | GT | 0 | 0 |
| main | Q-T3S2 | GF | 95 | 100 |
| main | Q-T3S3 | G0 | 0 | 0 |
| main | Q-T3S3 | GT | 0 | 0 |
| main | Q-T3S3 | GF | 104 | 108 |
| tail | Q-T1S2 | G0 | 0 | 0 |
| tail | Q-T1S2 | GT | 1 | 1 |
| tail | Q-T1S2 | GF | 74 | 109 |
| tail | Q-T3S1 | G0 | 0 | 0 |
| tail | Q-T3S1 | GT | 0 | 0 |
| tail | Q-T3S1 | GF | 41 | 41 |

## 2. The interval

- The wider of the percentile bootstrap (B = 200, day-and-stock clusters) and the delete-one-day
  jackknife with Student t at G - 1 = 7 df.
- Diagnostic only: a delete-one-stock jackknife with t at S - 1 df, on the first 50 panels of each
  combination.

## 3. rho* (400 panels per cell and generator)

Panels are as in 3c: confirmation-sized, exploration day lengths, 3a sign model, log-normal sizes, R2 clip
on the moments' sizes only.

**Choices.**
- rho* is evaluated on the ceiling estimator (stock s at T iff T <= rho N_s), the rule Phase 4 applies,
  with a ceiling at each bin edge 1/160, 1/80, 1/40, 1/20, 1/10.
- The 3c-style bands are a diagnostic on the first 100 panels. In a band the matched set changes at every
  increment, so the small-sample bias of each log of a matched mean does not telescope. In the 10-panel
  pilot this showed up as a +0.05 bias in the narrowest band, absent from the 10^6 truth.
- An element is evaluable when at least 50 % of panels give an estimate; an undefined interval counts as
  not covering.
- An estimator passes a bin when it is evaluable there and every evaluable element passes.
- rho* = the largest bin where all four criterion estimators pass, for each generator. The cell's rho* =
  the largest bin where all 12 (generator, estimator) pairs pass.
- The Monte Carlo SE of every coverage is sqrt(c (1 - c) / 400), about 0.018 at 0.85; it is reported with
  the lowest coverage.

Per cell, generator and estimator (largest passing ceiling; "none" = no ceiling passes):

| sign model | cell | G0 zeta_V | G0 zeta_K | G0 pf | G0 Rratio | GT zeta_V | GT zeta_K | GT pf | GT Rratio | GF zeta_V | GF zeta_K | GF pf | GF Rratio | rho* (cell) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| main | Q-T1S1 | 1/10 | 1/40 | 1/10 | 1/40 | 1/10 | 1/20 | 1/10 | 1/20 | 1/10 | none | none | none | none |
| main | Q-T1S2 | 1/10 | none | 1/20 | 1/160 | 1/10 | 1/40 | 1/20 | 1/160 | 1/10 | none | none | 1/40 | none |
| main | Q-T1S3 | 1/80 | none | none | 1/160 | 1/80 | none | none | 1/160 | none | 1/80 | 1/80 | 1/40 | none |
| main | Q-T2S1 | 1/160 | none | 1/40 | 1/80 | 1/160 | none | 1/40 | 1/80 | 1/160 | 1/80 | 1/40 | 1/10 | none |
| main | Q-T2S2 | 1/80 | none | none | 1/40 | 1/80 | none | none | 1/40 | 1/80 | 1/80 | 1/80 | 1/40 | none |
| main | Q-T2S3 | 1/10 | none | 1/10 | 1/10 | 1/10 | none | 1/10 | 1/10 | 1/10 | none | none | 1/40 | none |
| main | Q-T3S1 | 1/10 | 1/40 | 1/10 | 1/40 | 1/10 | 1/20 | 1/10 | 1/10 | 1/10 | 1/80 | 1/80 | 1/80 | none |
| main | Q-T3S2 | none | none | 1/10 | 1/20 | none | none | 1/10 | 1/20 | 1/80 | 1/80 | 1/80 | 1/10 | none |
| main | Q-T3S3 | none | none | 1/10 | 1/40 | none | none | 1/10 | none | none | none | none | none | none |
| tail | Q-T1S2 | none | none | 1/160 | none | none | none | none | none | none | none | 1/10 | 1/160 | none |
| tail | Q-T3S1 | none | none | 1/80 | none | none | none | 1/80 | none | none | 1/20 | 1/10 | 1/20 | none |

Binding constraint at each ceiling (failing elements of the evaluable ones, the worst bias, and the lowest
coverage with its MC SE):

| sign model | cell | T/N bin | passes | not evaluable | failing elements (of evaluable) | worst bias | lowest coverage (MC SE) |
|---|---|---|---|---|---|---|---|
| main | Q-T1S1 | 1/160 | no | - | 32 of 108 | GF Rratio[16,256] +0.214 | GT zeta_V(1)[16,256] 0.058 (0.012) |
| main | Q-T1S1 | 1/80 | no | - | 25 of 108 | GF Rratio[16,256] -1.011 | GF zeta_V(1)[16,256] 0.427 (0.025) |
| main | Q-T1S1 | 1/40 | no | - | 10 of 108 | GF Rratio[16,256] -1.011 | GF zeta_K(0)[16,256] 0.927 (0.013) |
| main | Q-T1S1 | 1/20 | no | - | 13 of 108 | GF zeta_K(0)[16,256] -0.100 | G0 Rratio[4,32] 0.848 (0.018) |
| main | Q-T1S1 | 1/10 | no | - | 27 of 108 | GF zeta_K(0)[16,256] -0.103 | G0 Rratio[4,32] 0.618 (0.024) |
| main | Q-T1S2 | 1/160 | no | - | 37 of 72 | GF Rratio[8,64] +0.055 | G0 zeta_K(0)[8,64] 0.448 (0.025) |
| main | Q-T1S2 | 1/80 | no | - | 40 of 108 | GF Rratio[16,256] -0.351 | GT zeta_V(0.75)[16,256] 0.065 (0.012) |
| main | Q-T1S2 | 1/40 | no | - | 31 of 108 | GF zeta_K(0)[16,256] -0.098 | GT zeta_V(0.75)[16,256] 0.312 (0.023) |
| main | Q-T1S2 | 1/20 | no | - | 32 of 108 | GF Rratio[16,256] +0.096 | G0 zeta_K(0)[4,32] 0.318 (0.023) |
| main | Q-T1S2 | 1/10 | no | - | 45 of 108 | GF zeta_K(0)[16,256] -0.058 | G0 zeta_K(0)[4,32] 0.142 (0.017) |
| main | Q-T1S3 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 0 of 3 | - | - |
| main | Q-T1S3 | 1/80 | no | - | 18 of 39 | G0 zeta_K(0)[4,32] +0.099 | G0 zeta_K(0)[4,32] 0.000 (0.000) |
| main | Q-T1S3 | 1/40 | no | - | 43 of 72 | G0 zeta_K(0)[8,64] +0.079 | G0 Rratio[4,32] 0.110 (0.016) |
| main | Q-T1S3 | 1/20 | no | - | 81 of 108 | GF zeta_V(0.5)[16,256] -0.156 | GF zeta_V(0.25)[16,256] 0.005 (0.004) |
| main | Q-T1S3 | 1/10 | no | - | 89 of 108 | GF Rratio[16,256] +0.126 | G0 Rratio[4,32] 0.000 (0.000) |
| main | Q-T2S1 | 1/160 | no | - | 12 of 39 | GT zeta_K(0)[4,32] +0.076 | G0 zeta_K(0)[4,32] 0.000 (0.000) |
| main | Q-T2S1 | 1/80 | no | - | 25 of 72 | G0 zeta_K(0)[8,64] +0.055 | G0 zeta_K(0.25)[8,64] 0.237 (0.021) |
| main | Q-T2S1 | 1/40 | no | - | 45 of 108 | GF zeta_V(0)[16,256] -0.158 | GT zeta_V(1)[16,256] 0.003 (0.002) |
| main | Q-T2S1 | 1/20 | no | - | 62 of 108 | GF zeta_V(0)[16,256] -0.102 | GF zeta_V(0.25)[16,256] 0.060 (0.012) |
| main | Q-T2S1 | 1/10 | no | - | 68 of 108 | G0 zeta_K(0)[16,256] +0.100 | G0 zeta_K(0)[8,64] 0.720 (0.022) |
| main | Q-T2S2 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 0 of 3 | - | - |
| main | Q-T2S2 | 1/80 | no | - | 12 of 39 | G0 zeta_K(0)[4,32] +0.078 | G0 zeta_K(0)[4,32] 0.048 (0.011) |
| main | Q-T2S2 | 1/40 | no | - | 29 of 72 | G0 zeta_K(0)[8,64] +0.058 | GF zeta_V(0)[8,64] 0.482 (0.025) |
| main | Q-T2S2 | 1/20 | no | - | 56 of 108 | GF Rratio[16,256] -5.531 | GF zeta_V(0.5)[16,256] 0.020 (0.007) |
| main | Q-T2S2 | 1/10 | no | - | 58 of 108 | GF zeta_V(0)[16,256] -0.103 | GF zeta_V(0.25)[16,256] 0.122 (0.016) |
| main | Q-T2S3 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, G0 Rratio, GT zeta_V, GT zeta_K, GT pf, GT Rratio, GF zeta_V, GF zeta_K, GF pf, GF Rratio | 0 of 0 | - | - |
| main | Q-T2S3 | 1/80 | no | G0 zeta_V, G0 zeta_K, G0 pf, G0 Rratio, GT zeta_V, GT zeta_K, GT pf, GT Rratio, GF zeta_V, GF zeta_K, GF pf, GF Rratio | 0 of 0 | - | - |
| main | Q-T2S3 | 1/40 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 0 of 3 | - | - |
| main | Q-T2S3 | 1/20 | no | - | 18 of 39 | GF Rratio[8,64] +0.172 | G0 zeta_K(0)[4,32] 0.168 (0.019) |
| main | Q-T2S3 | 1/10 | no | - | 18 of 72 | GF zeta_K(0)[8,64] -0.064 | G0 zeta_K(0)[8,64] 0.887 (0.016) |
| main | Q-T3S1 | 1/160 | no | - | 12 of 39 | GT zeta_K(0)[4,32] +0.061 | G0 zeta_K(0)[4,32] 0.000 (0.000) |
| main | Q-T3S1 | 1/80 | no | - | 25 of 72 | G0 zeta_V(0.25)[8,64] -0.029 | G0 zeta_V(0.75)[8,64] 0.440 (0.025) |
| main | Q-T3S1 | 1/40 | no | - | 18 of 108 | G0 zeta_V(0)[16,256] -0.103 | GF zeta_V(0.5)[16,256] 0.010 (0.005) |
| main | Q-T3S1 | 1/20 | no | - | 25 of 108 | G0 zeta_V(0.25)[16,256] -0.057 | G0 zeta_V(0.75)[16,256] 0.223 (0.021) |
| main | Q-T3S1 | 1/10 | no | - | 27 of 108 | GF zeta_K(0)[16,256] -0.095 | G0 zeta_K(0)[8,64] 0.718 (0.023) |
| main | Q-T3S2 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 1 of 3 | GF Rratio[4,32] +0.039 | GF Rratio[4,32] 0.868 (0.017) |
| main | Q-T3S2 | 1/80 | no | - | 16 of 39 | GT zeta_K(0)[4,32] +0.058 | G0 zeta_K(0)[4,32] 0.245 (0.022) |
| main | Q-T3S2 | 1/40 | no | - | 33 of 72 | GF Rratio[8,64] +0.098 | GT zeta_V(0.5)[8,64] 0.672 (0.023) |
| main | Q-T3S2 | 1/20 | no | - | 36 of 108 | GF Rratio[16,256] -0.271 | GT zeta_V(0.75)[16,256] 0.258 (0.022) |
| main | Q-T3S2 | 1/10 | no | - | 40 of 108 | GF zeta_K(0)[16,256] -0.114 | GT zeta_V(0.5)[16,256] 0.588 (0.025) |
| main | Q-T3S3 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 3 of 3 | GF Rratio[4,32] +0.111 | GT Rratio[4,32] 0.820 (0.019) |
| main | Q-T3S3 | 1/80 | no | - | 29 of 39 | GF Rratio[8,64] +0.351 | G0 zeta_K(0)[4,32] 0.180 (0.019) |
| main | Q-T3S3 | 1/40 | no | - | 31 of 72 | GF Rratio[8,64] +0.075 | G0 zeta_V(0.75)[8,64] 0.682 (0.023) |
| main | Q-T3S3 | 1/20 | no | - | 41 of 108 | GF Rratio[16,256] -1.336 | GT zeta_V(0.75)[16,256] 0.253 (0.022) |
| main | Q-T3S3 | 1/10 | no | - | 38 of 108 | GF Rratio[16,256] +0.560 | GT zeta_V(0.75)[16,256] 0.657 (0.024) |
| tail | Q-T1S2 | 1/160 | no | - | 57 of 72 | G0 zeta_V(0.25)[8,64] -0.078 | G0 Rratio[4,32] 0.115 (0.016) |
| tail | Q-T1S2 | 1/80 | no | - | 88 of 108 | G0 zeta_V(0.5)[16,256] -0.153 | G0 Rratio[4,32] 0.003 (0.002) |
| tail | Q-T1S2 | 1/40 | no | - | 92 of 108 | G0 zeta_V(0.5)[16,256] -0.104 | G0 Rratio[4,32] 0.000 (0.000) |
| tail | Q-T1S2 | 1/20 | no | - | 94 of 108 | G0 zeta_K(0)[16,256] +0.129 | G0 zeta_K(0)[4,32] 0.000 (0.000) |
| tail | Q-T1S2 | 1/10 | no | - | 95 of 108 | G0 zeta_K(0)[16,256] +0.169 | G0 zeta_K(0)[4,32] 0.000 (0.000) |
| tail | Q-T3S1 | 1/160 | no | - | 29 of 39 | GT zeta_V(0)[4,32] -0.053 | GT zeta_V(0.75)[4,32] 0.085 (0.014) |
| tail | Q-T3S1 | 1/80 | no | - | 45 of 72 | GT zeta_V(0)[8,64] -0.073 | GT zeta_V(0.75)[8,64] 0.113 (0.016) |
| tail | Q-T3S1 | 1/40 | no | - | 74 of 108 | GT zeta_V(0.5)[16,256] -0.150 | G0 zeta_V(0.75)[16,256] 0.000 (0.000) |
| tail | Q-T3S1 | 1/20 | no | - | 79 of 108 | GT zeta_V(0.5)[16,256] -0.099 | GT zeta_V(0.75)[16,256] 0.050 (0.011) |
| tail | Q-T3S1 | 1/10 | no | - | 84 of 108 | G0 zeta_K(0)[16,256] +0.150 | G0 zeta_K(0)[8,64] 0.142 (0.017) |

Band diagnostic (3c-style, first 100 panels):

| sign model | cell | band <= 1/160: failing / evaluable, worst bias | band <= 1/80: failing / evaluable, worst bias | band <= 1/40: failing / evaluable, worst bias | band <= 1/20: failing / evaluable, worst bias | band <= 1/10: failing / evaluable, worst bias |
|---|---|---|---|---|---|---|
| main | Q-T1S1 | 28 / 108, GF Rratio[16,256] +0.350 | 98 / 108, GT zeta_K(0)[4,32] +0.086 | 41 / 105, GF zeta_K(0.5)[16,256] -0.095 | 30 / 108, GF Rratio[16,256] -0.210 | 49 / 108, GF Rratio[16,256] -2.803 |
| main | Q-T1S2 | 30 / 72, GF Rratio[8,64] +0.051 | 51 / 108, GF zeta_V(0.5)[16,256] -0.092 | 41 / 108, GF Rratio[16,256] +0.120 | 33 / 108, GF Rratio[16,256] +0.358 | 36 / 108, GF Rratio[16,256] +0.080 |
| main | Q-T1S3 | 0 / 3 | 21 / 39, G0 zeta_K(0)[4,32] +0.101 | 41 / 72, G0 zeta_K(0)[8,64] +0.082 | 82 / 108, GF Rratio[16,256] +0.545 | 85 / 108, GF Rratio[16,256] -8.829 |
| main | Q-T2S1 | 12 / 39, GT zeta_K(0)[4,32] +0.077 | 31 / 72, G0 zeta_K(0)[8,64] +0.054 | 51 / 108, GF zeta_V(0)[16,256] -0.162 | 61 / 108, GF Rratio[16,256] +3.224 | 86 / 108, GF Rratio[16,256] -0.599 |
| main | Q-T2S2 | 0 / 3 | 13 / 39, G0 zeta_K(0)[4,32] +0.078 | 33 / 72, GF zeta_K(0)[8,64] -0.068 | 59 / 108, GF Rratio[16,256] -0.445 | 67 / 108, GF Rratio[8,64] -0.638 |
| main | Q-T2S3 | 0 / 0 | 0 / 0 | 0 / 3 | 2 / 6, GF Rratio[8,64] +0.094 | 1 / 3, GF Rratio[8,64] +0.332 |
| main | Q-T3S1 | 13 / 39, GT zeta_K(0)[4,32] +0.061 | 25 / 72, GT zeta_K(0)[8,64] +0.032 | 24 / 108, G0 zeta_V(0)[16,256] -0.105 | 39 / 108, GF Rratio[16,256] +0.348 | 55 / 108, GF Rratio[8,64] -0.414 |
| main | Q-T3S2 | 2 / 3, GF Rratio[4,32] +0.028 | 36 / 39, GF Rratio[8,64] -0.064 | 38 / 72, GF Rratio[8,64] -0.456 | 73 / 108, GF Rratio[16,256] -0.270 | 42 / 72, GF Rratio[8,64] +0.503 |
| main | Q-T3S3 | 2 / 3, GF Rratio[4,32] +0.067 | 28 / 39, GF Rratio[8,64] +0.202 | 38 / 72, GF Rratio[8,64] +0.160 | 53 / 108, GF Rratio[16,256] -2.010 | 37 / 72, GF Rratio[16,256] +3.696 |
| tail | Q-T1S2 | 52 / 72, G0 zeta_V(0.25)[8,64] -0.079 | 93 / 108, G0 zeta_V(0.5)[16,256] -0.153 | 93 / 108, GF Rratio[16,256] +2.165 | 91 / 108, G0 zeta_K(0)[16,256] +0.135 | 90 / 108, G0 zeta_K(0)[16,256] +0.176 |
| tail | Q-T3S1 | 29 / 39, GT zeta_V(0)[4,32] -0.053 | 43 / 72, GT zeta_V(0)[8,64] -0.073 | 69 / 108, GT zeta_V(0.5)[16,256] -0.150 | 72 / 108, G0 zeta_K(0)[16,256] +0.105 | 93 / 108, GF Rratio[4,32] -10.397 |

**Measurement.**
- **G0 and GT.** The worst biases are zeta_K(0) on the short ranges, +0.03 to +0.10, the size of the
  finite-N shift of the target in the table below. Next come zeta_V at the smallest ceilings, -0.03 to
  -0.10, from the pooling (headline 3).
- **GF.** zeta_K's truth is 1 within 0.002; its estimator passes at the smaller ceilings in 5 of 9 cells.
  R_{1/2}/R_0 fails in most cells (biases up to 0.35 at the ceilings, several units in narrow bands): R_0
  is near 0 under GF, so the ratio is ill-conditioned.
- **Tail-lag sign model** (Q-T1S2, Q-T3S1). More estimators fail than under the main model. rho* does not
  improve.

### Old against new, same panels (Q-T1S2 and Q-T3S1)

For each ceiling, estimator and generator: max |bias| and min coverage, new path then old (3c) path, and
the stock-cluster jackknife's minimum coverage (new path, first 50 panels).

| cell | gen | T/N bin | estimator | new: max abs bias / min coverage | old: max abs bias / min coverage | stock jackknife (new, first 50 panels): min coverage |
|---|---|---|---|---|---|---|
| Q-T1S2 | G0 | 1/160 | zeta_V | 0.019 / 0.61 | 0.001 / 0.86 | 0.32 |
| Q-T1S2 | G0 | 1/160 | zeta_K | 0.032 / 0.45 | 0.006 / 0.86 | 0.12 |
| Q-T1S2 | G0 | 1/160 | pf | 0.012 / 0.79 | 0.002 / 0.88 | 0.36 |
| Q-T1S2 | G0 | 1/160 | Rratio | 0.005 / 0.89 | 0.005 / 0.91 | 0.68 |
| Q-T1S2 | G0 | 1/80 | zeta_V | 0.088 / 0.10 | 0.003 / 0.85 | 0.00 |
| Q-T1S2 | G0 | 1/80 | zeta_K | 0.019 / 0.78 | 0.011 / 0.84 | 0.46 |
| Q-T1S2 | G0 | 1/80 | pf | 0.005 / 0.88 | 0.005 / 0.86 | 0.56 |
| Q-T1S2 | G0 | 1/80 | Rratio | 0.010 / 0.76 | 0.010 / 0.73 | 0.36 |
| Q-T1S2 | G0 | 1/40 | zeta_V | 0.052 / 0.39 | 0.004 / 0.94 | 0.30 |
| Q-T1S2 | G0 | 1/40 | zeta_K | 0.017 / 0.66 | 0.016 / 0.79 | 0.34 |
| Q-T1S2 | G0 | 1/40 | pf | 0.006 / 0.93 | 0.006 / 0.95 | 0.70 |
| Q-T1S2 | G0 | 1/40 | Rratio | 0.013 / 0.55 | 0.013 / 0.55 | 0.16 |
| Q-T1S2 | G0 | 1/20 | zeta_V | 0.007 / 0.87 | 0.009 / 0.80 | 0.72 |
| Q-T1S2 | G0 | 1/20 | zeta_K | 0.034 / 0.32 | 0.030 / 0.65 | 0.08 |
| Q-T1S2 | G0 | 1/20 | pf | 0.011 / 0.86 | 0.011 / 0.87 | 0.48 |
| Q-T1S2 | G0 | 1/20 | Rratio | 0.016 / 0.35 | 0.016 / 0.33 | 0.04 |
| Q-T1S2 | G0 | 1/10 | zeta_V | 0.003 / 0.95 | 0.011 / 0.81 | 0.84 |
| Q-T1S2 | G0 | 1/10 | zeta_K | 0.054 / 0.14 | 0.041 / 0.61 | 0.00 |
| Q-T1S2 | G0 | 1/10 | pf | 0.016 / 0.72 | 0.014 / 0.89 | 0.32 |
| Q-T1S2 | G0 | 1/10 | Rratio | 0.020 / 0.22 | 0.020 / 0.27 | 0.04 |
| Q-T1S2 | GT | 1/160 | zeta_V | 0.019 / 0.56 | 0.001 / 0.82 | 0.30 |
| Q-T1S2 | GT | 1/160 | zeta_K | 0.033 / 0.55 | 0.004 / 0.83 | 0.24 |
| Q-T1S2 | GT | 1/160 | pf | 0.013 / 0.79 | 0.002 / 0.84 | 0.40 |
| Q-T1S2 | GT | 1/160 | Rratio | 0.004 / 0.96 | 0.004 / 0.96 | 0.78 |
| Q-T1S2 | GT | 1/80 | zeta_V | 0.087 / 0.07 | 0.002 / 0.81 | 0.02 |
| Q-T1S2 | GT | 1/80 | zeta_K | 0.021 / 0.80 | 0.006 / 0.84 | 0.46 |
| Q-T1S2 | GT | 1/80 | pf | 0.011 / 0.82 | 0.005 / 0.84 | 0.62 |
| Q-T1S2 | GT | 1/80 | Rratio | 0.008 / 0.82 | 0.008 / 0.83 | 0.54 |
| Q-T1S2 | GT | 1/40 | zeta_V | 0.053 / 0.31 | 0.003 / 0.93 | 0.34 |
| Q-T1S2 | GT | 1/40 | zeta_K | 0.013 / 0.86 | 0.008 / 0.96 | 0.56 |
| Q-T1S2 | GT | 1/40 | pf | 0.007 / 0.94 | 0.006 / 0.97 | 0.80 |
| Q-T1S2 | GT | 1/40 | Rratio | 0.012 / 0.70 | 0.012 / 0.68 | 0.34 |
| Q-T1S2 | GT | 1/20 | zeta_V | 0.009 / 0.85 | 0.009 / 0.79 | 0.76 |
| Q-T1S2 | GT | 1/20 | zeta_K | 0.020 / 0.73 | 0.011 / 0.85 | 0.38 |
| Q-T1S2 | GT | 1/20 | pf | 0.010 / 0.90 | 0.008 / 0.84 | 0.60 |
| Q-T1S2 | GT | 1/20 | Rratio | 0.017 / 0.51 | 0.017 / 0.50 | 0.22 |
| Q-T1S2 | GT | 1/10 | zeta_V | 0.003 / 0.95 | 0.012 / 0.78 | 0.82 |
| Q-T1S2 | GT | 1/10 | zeta_K | 0.041 / 0.56 | 0.011 / 0.96 | 0.26 |
| Q-T1S2 | GT | 1/10 | pf | 0.020 / 0.85 | 0.009 / 0.98 | 0.50 |
| Q-T1S2 | GT | 1/10 | Rratio | 0.022 / 0.36 | 0.022 / 0.36 | 0.12 |
| Q-T1S2 | GF | 1/160 | zeta_V | 0.019 / 0.61 | 0.001 / 0.85 | 0.34 |
| Q-T1S2 | GF | 1/160 | zeta_K | 0.038 / 0.88 | 0.112 / 0.84 | 0.48 |
| Q-T1S2 | GF | 1/160 | pf | 0.027 / 0.88 | 0.087 / 0.89 | 0.56 |
| Q-T1S2 | GF | 1/160 | Rratio | 0.055 / 0.98 | 0.055 / 0.98 | 0.90 |
| Q-T1S2 | GF | 1/80 | zeta_V | 0.092 / 0.08 | 0.003 / 0.83 | 0.04 |
| Q-T1S2 | GF | 1/80 | zeta_K | 0.070 / 0.84 | 0.152 / 0.86 | 0.50 |
| Q-T1S2 | GF | 1/80 | pf | 0.047 / 0.84 | 0.113 / 0.86 | 0.50 |
| Q-T1S2 | GF | 1/80 | Rratio | 0.351 / 0.89 | 0.351 / 0.88 | 0.52 |
| Q-T1S2 | GF | 1/40 | zeta_V | 0.054 / 0.38 | 0.004 / 0.93 | 0.22 |
| Q-T1S2 | GF | 1/40 | zeta_K | 0.098 / 0.94 | 0.153 / 0.80 | 0.72 |
| Q-T1S2 | GF | 1/40 | pf | 0.066 / 0.94 | 0.079 / 0.96 | 0.72 |
| Q-T1S2 | GF | 1/40 | Rratio | 0.018 / 0.98 | 0.018 / 0.97 | 0.82 |
| Q-T1S2 | GF | 1/20 | zeta_V | 0.008 / 0.88 | 0.010 / 0.83 | 0.76 |
| Q-T1S2 | GF | 1/20 | zeta_K | 0.088 / 0.97 | 0.126 / 0.63 | 0.78 |
| Q-T1S2 | GF | 1/20 | pf | 0.050 / 0.97 | 0.103 / 0.64 | 0.84 |
| Q-T1S2 | GF | 1/20 | Rratio | 0.096 / 0.98 | 0.096 / 0.99 | 0.82 |
| Q-T1S2 | GF | 1/10 | zeta_V | 0.004 / 0.95 | 0.013 / 0.77 | 0.84 |
| Q-T1S2 | GF | 1/10 | zeta_K | 0.058 / 0.98 | 0.137 / 0.77 | 0.84 |
| Q-T1S2 | GF | 1/10 | pf | 0.032 / 0.99 | 0.199 / 0.85 | 0.90 |
| Q-T1S2 | GF | 1/10 | Rratio | 0.021 / 0.98 | 0.021 / 0.99 | 0.82 |
| Q-T3S1 | G0 | 1/160 | zeta_V | 0.005 / 0.89 | - / - | 0.70 |
| Q-T3S1 | G0 | 1/160 | zeta_K | 0.056 / 0.00 | - / - | 0.00 |
| Q-T3S1 | G0 | 1/160 | pf | 0.020 / 0.28 | - / - | 0.04 |
| Q-T3S1 | G0 | 1/160 | Rratio | 0.005 / 0.95 | 0.005 / 0.94 | 0.76 |
| Q-T3S1 | G0 | 1/80 | zeta_V | 0.029 / 0.44 | 0.003 / 0.95 | 0.20 |
| Q-T3S1 | G0 | 1/80 | zeta_K | 0.026 / 0.62 | 0.008 / 0.93 | 0.30 |
| Q-T3S1 | G0 | 1/80 | pf | 0.008 / 0.95 | 0.004 / 0.97 | 0.78 |
| Q-T3S1 | G0 | 1/80 | Rratio | 0.006 / 0.95 | 0.006 / 0.95 | 0.70 |
| Q-T3S1 | G0 | 1/40 | zeta_V | 0.103 / 0.02 | 0.005 / 0.94 | 0.00 |
| Q-T3S1 | G0 | 1/40 | zeta_K | 0.017 / 0.88 | 0.015 / 0.93 | 0.66 |
| Q-T3S1 | G0 | 1/40 | pf | 0.008 / 0.96 | 0.007 / 0.97 | 0.82 |
| Q-T3S1 | G0 | 1/40 | Rratio | 0.014 / 0.88 | 0.014 / 0.90 | 0.62 |
| Q-T3S1 | G0 | 1/20 | zeta_V | 0.057 / 0.22 | 0.006 / 0.94 | 0.08 |
| Q-T3S1 | G0 | 1/20 | zeta_K | 0.022 / 0.82 | 0.018 / 0.91 | 0.52 |
| Q-T3S1 | G0 | 1/20 | pf | 0.010 / 0.94 | 0.009 / 0.94 | 0.64 |
| Q-T3S1 | G0 | 1/20 | Rratio | 0.017 / 0.85 | 0.017 / 0.85 | 0.48 |
| Q-T3S1 | G0 | 1/10 | zeta_V | 0.006 / 0.97 | 0.015 / 0.85 | 0.84 |
| Q-T3S1 | G0 | 1/10 | zeta_K | 0.054 / 0.72 | 0.034 / 0.93 | 0.32 |
| Q-T3S1 | G0 | 1/10 | pf | 0.018 / 0.90 | 0.013 / 0.94 | 0.62 |
| Q-T3S1 | G0 | 1/10 | Rratio | 0.020 / 0.83 | 0.020 / 0.84 | 0.48 |
| Q-T3S1 | GT | 1/160 | zeta_V | 0.007 / 0.94 | - / - | 0.72 |
| Q-T3S1 | GT | 1/160 | zeta_K | 0.061 / 0.01 | - / - | 0.00 |
| Q-T3S1 | GT | 1/160 | pf | 0.022 / 0.37 | - / - | 0.16 |
| Q-T3S1 | GT | 1/160 | Rratio | 0.004 / 0.96 | 0.004 / 0.96 | 0.76 |
| Q-T3S1 | GT | 1/80 | zeta_V | 0.026 / 0.52 | 0.000 / 0.97 | 0.20 |
| Q-T3S1 | GT | 1/80 | zeta_K | 0.029 / 0.78 | 0.005 / 0.98 | 0.48 |
| Q-T3S1 | GT | 1/80 | pf | 0.009 / 0.95 | 0.004 / 0.97 | 0.82 |
| Q-T3S1 | GT | 1/80 | Rratio | 0.005 / 0.97 | 0.005 / 0.98 | 0.80 |
| Q-T3S1 | GT | 1/40 | zeta_V | 0.099 / 0.01 | 0.002 / 0.96 | 0.00 |
| Q-T3S1 | GT | 1/40 | zeta_K | 0.016 / 0.96 | 0.007 / 0.98 | 0.78 |
| Q-T3S1 | GT | 1/40 | pf | 0.009 / 0.97 | 0.006 / 0.98 | 0.80 |
| Q-T3S1 | GT | 1/40 | Rratio | 0.015 / 0.94 | 0.015 / 0.94 | 0.80 |
| Q-T3S1 | GT | 1/20 | zeta_V | 0.053 / 0.28 | 0.004 / 0.96 | 0.12 |
| Q-T3S1 | GT | 1/20 | zeta_K | 0.018 / 0.94 | 0.007 / 0.98 | 0.78 |
| Q-T3S1 | GT | 1/20 | pf | 0.010 / 0.96 | 0.006 / 0.99 | 0.80 |
| Q-T3S1 | GT | 1/20 | Rratio | 0.016 / 0.94 | 0.016 / 0.94 | 0.74 |
| Q-T3S1 | GT | 1/10 | zeta_V | 0.002 / 0.97 | 0.011 / 0.90 | 0.84 |
| Q-T3S1 | GT | 1/10 | zeta_K | 0.039 / 0.92 | 0.026 / 0.99 | 0.68 |
| Q-T3S1 | GT | 1/10 | pf | 0.018 / 0.97 | 0.005 / 0.99 | 0.80 |
| Q-T3S1 | GT | 1/10 | Rratio | 0.018 / 0.93 | 0.018 / 0.93 | 0.68 |
| Q-T3S1 | GF | 1/160 | zeta_V | 0.006 / 0.94 | - / - | 0.74 |
| Q-T3S1 | GF | 1/160 | zeta_K | 0.007 / 0.97 | - / - | 0.90 |
| Q-T3S1 | GF | 1/160 | pf | 0.005 / 0.98 | - / - | 0.96 |
| Q-T3S1 | GF | 1/160 | Rratio | 0.015 / 0.98 | 0.015 / 0.98 | 0.82 |
| Q-T3S1 | GF | 1/80 | zeta_V | 0.027 / 0.49 | 0.002 / 0.96 | 0.22 |
| Q-T3S1 | GF | 1/80 | zeta_K | 0.012 / 0.98 | 0.099 / 0.94 | 0.90 |
| Q-T3S1 | GF | 1/80 | pf | 0.008 / 0.99 | 0.067 / 0.99 | 0.92 |
| Q-T3S1 | GF | 1/80 | Rratio | 0.012 / 0.98 | 0.012 / 0.98 | 0.86 |
| Q-T3S1 | GF | 1/40 | zeta_V | 0.103 / 0.01 | 0.004 / 0.96 | 0.00 |
| Q-T3S1 | GF | 1/40 | zeta_K | 0.030 / 0.99 | 0.140 / 0.94 | 0.88 |
| Q-T3S1 | GF | 1/40 | pf | 0.022 / 0.98 | 0.071 / 1.00 | 0.92 |
| Q-T3S1 | GF | 1/40 | Rratio | 0.065 / 0.98 | 0.065 / 0.99 | 0.84 |
| Q-T3S1 | GF | 1/20 | zeta_V | 0.056 / 0.26 | 0.006 / 0.95 | 0.12 |
| Q-T3S1 | GF | 1/20 | zeta_K | 0.054 / 0.98 | 0.119 / 0.93 | 0.86 |
| Q-T3S1 | GF | 1/20 | pf | 0.043 / 0.99 | 0.045 / 0.99 | 0.90 |
| Q-T3S1 | GF | 1/20 | Rratio | 0.054 / 0.99 | 0.054 / 0.98 | 0.84 |
| Q-T3S1 | GF | 1/10 | zeta_V | 0.004 / 0.96 | 0.013 / 0.89 | 0.82 |
| Q-T3S1 | GF | 1/10 | zeta_K | 0.095 / 0.98 | 0.106 / 0.92 | 0.84 |
| Q-T3S1 | GF | 1/10 | pf | 0.067 / 0.99 | 0.103 / 0.95 | 0.94 |
| Q-T3S1 | GF | 1/10 | Rratio | 0.028 / 0.99 | 0.028 / 0.98 | 0.84 |

**Measurement.**
- zeta_K: the new path matches the old on bias (Q-T1S2 G0 at ceiling 1/40: 0.017 against 0.016), with
  no exclusion. Its coverage is lower there (0.66 against 0.79).
- zeta_V: the new path is worse at small ceilings (headline 3).
- The stock jackknife's coverage is below the primary interval's in every row of the table, often below
  0.5. It should stay a diagnostic.

### Finite-N truth (diagnostic, G0 and GT)

zeta_K(0) of the population at day length N, against the 10^6-run truth: about 2 x 10^6 runs per entry,
cap fixed at the population 99.9 % quantile.

| cell | gen | N = 300: zK(0) [4,32] / [8,64] / [16,256] | N = 1000: zK(0) [4,32] / [8,64] / [16,256] | N = 3000: zK(0) [4,32] / [8,64] / [16,256] | N = 10000: zK(0) [4,32] / [8,64] / [16,256] | N = 1000000: zK(0) [4,32] / [8,64] / [16,256] |
|---|---|---|---|---|---|---|
| Q-T1S1 | G0 | 1.188 / 1.102 / 1.172 | 1.175 / 1.076 / 1.088 | 1.165 / 1.057 / 1.055 | 1.167 / 1.055 / 1.034 | 1.161 / 1.049 / 1.018 |
| Q-T1S1 | GT | 1.075 / 0.971 / 1.056 | 1.060 / 0.940 / 0.944 | 1.051 / 0.921 / 0.902 | 1.055 / 0.922 / 0.879 | 1.049 / 0.917 / 0.867 |
| Q-T1S2 | G0 | 1.066 / 1.002 / 1.134 | 1.047 / 0.968 / 1.034 | 1.037 / 0.946 / 0.989 | 1.037 / 0.941 / 0.958 | 1.030 / 0.932 / 0.931 |
| Q-T1S2 | GT | 0.984 / 0.908 / 1.054 | 0.963 / 0.870 / 0.934 | 0.956 / 0.851 / 0.888 | 0.957 / 0.849 / 0.850 | 0.950 / 0.840 / 0.831 |
| Q-T1S3 | G0 | 1.009 / 0.948 / 1.167 | 0.980 / 0.889 / 1.006 | 0.966 / 0.855 / 0.930 | 0.959 / 0.836 / 0.841 | 0.945 / 0.812 / 0.762 |
| Q-T1S3 | GT | 0.985 / 0.925 / 1.154 | 0.956 / 0.866 / 0.990 | 0.944 / 0.836 / 0.928 | 0.938 / 0.818 / 0.834 | 0.924 / 0.794 / 0.765 |
| Q-T2S1 | G0 | 0.872 / 0.801 / 0.962 | 0.853 / 0.765 / 0.846 | 0.848 / 0.748 / 0.784 | 0.846 / 0.736 / 0.701 | 0.836 / 0.722 / 0.653 |
| Q-T2S1 | GT | 0.835 / 0.760 / 0.925 | 0.817 / 0.724 / 0.806 | 0.814 / 0.710 / 0.749 | 0.811 / 0.700 / 0.657 | 0.802 / 0.686 / 0.621 |
| Q-T2S2 | G0 | 0.874 / 0.805 / 0.973 | 0.854 / 0.766 / 0.854 | 0.849 / 0.749 / 0.789 | 0.846 / 0.737 / 0.704 | 0.836 / 0.722 / 0.654 |
| Q-T2S2 | GT | 0.847 / 0.777 / 0.951 | 0.827 / 0.739 / 0.831 | 0.824 / 0.724 / 0.776 | 0.822 / 0.714 / 0.684 | 0.812 / 0.699 / 0.641 |
| Q-T2S3 | G0 | 1.023 / 0.972 / 1.074 | 1.010 / 0.948 / 1.003 | 1.004 / 0.932 / 0.969 | 1.006 / 0.935 / 0.955 | 1.001 / 0.929 / 0.938 |
| Q-T2S3 | GT | 0.922 / 0.856 / 0.969 | 0.907 / 0.829 / 0.882 | 0.903 / 0.816 / 0.847 | 0.906 / 0.822 / 0.828 | 0.900 / 0.815 / 0.816 |
| Q-T3S1 | G0 | 0.963 / 0.924 / 1.086 | 0.944 / 0.889 / 0.987 | 0.936 / 0.869 / 0.935 | 0.934 / 0.864 / 0.899 | 0.928 / 0.854 / 0.870 |
| Q-T3S1 | GT | 0.876 / 0.824 / 0.995 | 0.856 / 0.786 / 0.879 | 0.851 / 0.769 / 0.825 | 0.851 / 0.766 / 0.778 | 0.844 / 0.758 / 0.762 |
| Q-T3S2 | G0 | 1.019 / 0.991 / 1.149 | 0.999 / 0.955 / 1.053 | 0.991 / 0.935 / 1.008 | 0.989 / 0.930 / 0.979 | 0.983 / 0.922 / 0.953 |
| Q-T3S2 | GT | 0.925 / 0.883 / 1.057 | 0.904 / 0.844 / 0.940 | 0.898 / 0.826 / 0.894 | 0.898 / 0.825 / 0.858 | 0.891 / 0.817 / 0.840 |
| Q-T3S3 | G0 | 1.186 / 1.157 / 1.266 | 1.168 / 1.125 / 1.181 | 1.156 / 1.104 / 1.150 | 1.158 / 1.104 / 1.132 | 1.151 / 1.096 / 1.114 |
| Q-T3S3 | GT | 1.081 / 1.042 / 1.176 | 1.059 / 1.003 / 1.064 | 1.050 / 0.985 / 1.034 | 1.053 / 0.988 / 1.011 | 1.046 / 0.979 / 0.996 |

**Measurement.**
- zeta_K's target moves with N in every cell, by +0.02 to +0.06 on [4, 32] and +0.04 to +0.14 on
  [8, 64] between N = 300 and 10^6. On [16, 256] it moves by +0.14 to +0.40.
- The shift is largest where the sign model has a strong near-constant day component (Q-T1S3, Q-T2S1,
  Q-T2S2).
- V is stationary within a day and does not move.

## 4. Reach (ceiling 1/40, no rho* found)

Confirmation stocks (N_min from `data/universe.tsv`, a descriptive column) contributing to each range
(N_min >= T2 / rho), and how many of them are above each N_floor:

| cell | rho used | confirm stocks | [4,32]: contributing / N_floor 100, 200, 400, 800 | [8,64]: contributing / N_floor 100, 200, 400, 800 | [16,256]: contributing / N_floor 100, 200, 400, 800 |
|---|---|---|---|---|---|
| Q-T1S1 | 1/40 (none found) | 161 | 22 / 22, 22, 22, 22 | 4 / 4, 4, 4, 4 | 0 / 0, 0, 0, 0 |
| Q-T1S2 | 1/40 (none found) | 216 | 23 / 23, 23, 23, 23 | 5 / 5, 5, 5, 5 | 0 / 0, 0, 0, 0 |
| Q-T1S3 | 1/40 (none found) | 252 | 5 / 5, 5, 5, 5 | 0 / 0, 0, 0, 0 | 0 / 0, 0, 0, 0 |
| Q-T2S1 | 1/40 (none found) | 80 | 17 / 17, 17, 17, 17 | 6 / 6, 6, 6, 6 | 2 / 2, 2, 2, 2 |
| Q-T2S2 | 1/40 (none found) | 50 | 3 / 3, 3, 3, 3 | 0 / 0, 0, 0, 0 | 0 / 0, 0, 0, 0 |
| Q-T2S3 | 1/40 (none found) | 17 | 0 / 0, 0, 0, 0 | 0 / 0, 0, 0, 0 | 0 / 0, 0, 0, 0 |
| Q-T3S1 | 1/40 (none found) | 50 | 21 / 21, 21, 21, 21 | 5 / 5, 5, 5, 5 | 1 / 1, 1, 1, 1 |
| Q-T3S2 | 1/40 (none found) | 27 | 6 / 6, 6, 6, 6 | 3 / 3, 3, 3, 3 | 0 / 0, 0, 0, 0 |
| Q-T3S3 | 1/40 (none found) | 23 | 1 / 1, 1, 1, 1 | 1 / 1, 1, 1, 1 | 0 / 0, 0, 0, 0 |

## 5. Tests and planted bugs

`tools/est/test_pooled2.py` has 7 checks and 6 planted bugs, all caught, each by its intended assertion
(each mutant rerun outside pytest). The 3c tests (`test_pooled.py`) still pass on the old path.

| planted bug | mutated function | caught by |
|---|---|---|
| per-stock-day denominators instead of day-averaged | `pooled.stock_curves_tilde` | check_day_averaged_denominators |
| stock-days excluded on the sign of K_a(1) | `pooled.stock_curves_tilde` | check_no_k1_exclusion |
| denominator at T instead of T = 1 (log slope of k~ differs from that of K/T) | `pooled.stock_curves_tilde` | check_log_slope_equals_K_over_T |
| jackknife with the normal quantile | `pooled.jk_interval` | check_jackknife_t_quantile |
| stock jackknife deleting days | `pooled.jackknife2` | check_stock_jackknife |
| non-positive pooled k~ kept (logs of a negative mean) | `pooled.chain_zeta` | check_nonpositive_pooled_drops_T |

check_old_path_is_3c shows that `norm="old"` reproduces the 3c core on the same moments. Full suite:
310 passed, 28 skipped.

## 6. Wall time and peak memory

| job | wall | peak RAM | notes |
|---|---|---|---|
| reduced study, band masks (discarded after the band finding) | 1.8 min | 0.43 GB main | |
| reduced study, ceiling masks (kept: same seeds as the full run) | 2.2 min | 0.43 GB main | estimate: at most 85 min, realistic 60 |
| study, 33 combinations x 400 panels | 68.5 min | main 0.51 GB, 1.37 GB all processes | 4 workers |
| finite-N truth (diagnostic, 72 runs) | 45.3 min | main 0.16 GB (workers not sampled) | estimated 12 to 15 min, ran 3 times longer |
| analyse | 4 s | 0.06 GB | |
| test suite | about 45 s | | |

Total compute about 2 h 15 min, over the 2 h target: the finite-N diagnostic, added after the study
pointed to it, overran its estimate. Every job checkpoints under `$ITCH_DATA/phase3c2/` and resumes.

## 7. What this leaves open (for a decision, not taken here)

- **Pooling.** Pool v as the mean of stock logs (3c, unbiased here), and handle k~ separately (log of
  matched means only where the pooled k~ is positive, or a positive-part rule).
- **Target.** Define the estimand at the stock's own day length rather than at N = 10^6, or accept
  zeta_K's N dependence as part of the result.
- **GF's R_{1/2}/R_0** cannot meet the criterion while R_0 is near 0.
