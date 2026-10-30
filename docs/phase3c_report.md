# Phase 3c report: operating characteristics of the pooled estimator (NASDAQ cells)

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken here. Real data were read
only through `explore.load_explore_runs` (exploration stocks, 8 study days). On real data, only signs and
sizes (steps 4 and 6) or mids alone (step 5) were used; nothing combined them. NASDAQ-listed cells only.

New code:
- `tools/est/pooled.py`: the primary pooled estimator; Phase 4 uses it.
- `tools/est/panel3c.py`: batched panels and the semi-synthetic generator.
- `tools/est/phase3c.py`: the driver.
- `tools/est/phase3c_report.py`: the tables.
- Tests: `test_pooled.py`, `test_panel3c.py`.
- Two small additions elsewhere: `gen3.signs` accepts the 3b two-component model, and `pooled.core`
  reports `vk_diff`.

Tables: `data/phase3/p3c_truth.tsv`, `p3c_bands.tsv`, `p3c_rho.tsv`, `p3c_ess_coverage.tsv`,
`p3c_amax.tsv`, `p3c_nfloor.tsv`, `p3c_sigma_bias.tsv`, `p3c_tmin_lmax.tsv`. Checkpoints and per-panel
results are in `$ITCH_DATA/phase3c/`.

## Headlines

1. **The candidate primary cannot discriminate (Measurement).** Under G0, the monotone null, slope
   zeta_V(a) / slope zeta_K(a) over [0, 0.5] on [8, 64] is 2.14 (SE 0.15) in Q-T2S1 and 2.19 (0.15) in
   Q-T2S2: within 0.2 of 2. There the difference slope zeta_V - 2 slope zeta_K is -0.006 +- 0.006 and
   -0.010 +- 0.008, zero within error. GT differs little from G0 in every cell: on [8, 64] its ratio is
   0.02 to 0.66 away from G0's, and its difference 0.004 to 0.039 away. The equality holds for GT in Q-T1S3 (diff 0.008), Q-T2S1 (-0.010) and, with the tail-lag sign model, Q-T3S1 (0.009). Under GF
   the zeta_K slope is 0, so the ratio has no finite truth.
2. **No rho* exists under the draft bounds (Measurement).** No cell has a T/N band at which every
   estimator, for every generator, has |bias| < 0.02 and 90 % coverage >= 0.85. The reasons:
   - the ratio misses the bias bound in every band from T/N 1/20 up, and in 4 of its 5 evaluable
     elements below 1/20 (G0 and GT);
   - GF's K-based estimators are biased by the K_a(1) exclusion, which drops 35 % to 75 % of GF
     stock-days;
   - at small T/N, coverage of single elements falls below 0.85 in up to half the elements;
   - at T/N >= 1/20, bias grows.
   The most permissive variant tried (G0 and GT only, ratio left out, coverage >= 0.81) gives a rho*
   only in Q-T1S1 (1/20), Q-T1S2 (1/80) and Q-T3S1 (1/40).
3. **Coverage does not improve with ESS (Measurement).** Coverage at median stock-day Kish ESS above 300
   (mean 0.91 to 0.94, minimum 0.81 to 0.82) is no better than at ESS 5 to 20 (mean 0.92 to 0.97).
4. **The long range is unreachable (Measurement).** Even at a ceiling of 1/10, at most 2 confirmation
   stocks per NASDAQ cell reach [64, 1024], and at most 6 reach [16, 256].

## 0. The primary pooled estimator

No existing function implemented it: `moments.day_moments` is per day, and `inference.py` is generic. So
`tools/est/pooled.py` was written to the spec:
- moments per stock-day with stride max(1, floor(T / 4)), windows inside a day, uncentred (centred arm
  available);
- v and k normalised by T times their T = 1 value;
- stock-days with K_a(1) not positive at the one-sided 95 % level excluded from k (the count is
  returned);
- days averaged within a stock first, then equal weight per stock;
- matched increments of each stock's own log curve, eligibility decided at the upper horizon, chained;
- zeta = 1 + OLS slope of the chained curve;
- R pooled as the mean Fisher z;
- R_{1/2} / R_0, and a*(T) as the vertex of a quadratic on [0, min(1.5, a_max)];
- inference: day-and-stock cluster bootstrap (B = 200, the same draws as `inference.cluster_bootstrap`)
  and the delete-one-day jackknife, the wider interval reported.

One weighted core runs the point estimate, every bootstrap replicate and every jackknife replicate, so all
three use the same code.

**Choices.**
- A bootstrap replicate undefined for one element is dropped for that element only; `inference.py` drops
  the whole replicate.
- The K_a(1) exclusion applies to k only; v and z keep every stock-day.
- A stock whose day-averaged curve is not positive at a horizon drops out of that increment.
- `vk_diff` = slope zeta_V - 2 slope zeta_K is added next to the ratio. It is the candidate's equality,
  and it stays defined when the zeta_K slope is 0.

## 1. Truth with the R2 clip (fixed 99.9 % quantile of the size model)

Days of 10^6 runs, a in {0, 0.25, ..., 1.5}, T = 1..4096 (log2), SE target 0.005 for a <= 1 and T <= 128
(cap 64 days, 48 for GM).

**Choices.**
- GT's truth is the 3b run: same cap, seed and settings.
- GM's cap is the 99.9 % quantile of trade sizes in a pilot day.
- The beta-floor variant is GM_floor_G01, which differs from PAPER2 only in beta_clip.
- The sign-model sensitivity uses the 3b tail-lag fit (H, tau) at Q-T1S2 and Q-T3S1, with G0's beta and
  GT's tent kept at the main values.

**Measurement.** The SE target was met in 2 of 35 runs (Q-T3S3 G0 and GT); every SE is in the table.

| cell | gen | days | domain SE | zeta_K(a) [8,64], a = 0..1.5 | zeta_V(a) [8,64] | P-F slope [8,64] | slope zV / slope zK [8,64] (a <= 0.5) | ratio [8,64] | ratio [16,256] | ratio [64,1024] | slope zV - 2 slope zK [8,64] | R1/2/R0 [8,64] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Q-T1S1 | G0 | 64 | 0.0054 | 1.049 1.022 0.993 0.964 0.937 0.914 0.894 | 1.385 1.373 1.332 1.260 1.175 1.105 1.062 | -0.105 (0.001) | -0.105 / -0.113 | 0.929 (0.005) | 0.940 (0.009) | 0.816 (0.014) | 0.121 (0.001) | 1.128 |
| Q-T1S1 | GT | 64 | 0.0054 | 0.917 0.888 0.859 0.833 0.813 0.799 0.790 | 1.385 1.373 1.332 1.260 1.175 1.105 1.062 | -0.087 (0.001) | -0.105 / -0.116 | 0.907 (0.007) | 0.920 (0.015) | 0.768 (0.030) | 0.126 (0.001) | 1.152 |
| Q-T1S1 | GF | 64 | 0.0296 | 0.999 0.998 0.998 0.999 1.000 1.001 1.002 | 1.398 1.385 1.341 1.265 1.176 1.105 1.062 | 0.002 (0.007) | -0.114 / -0.002 | 51.008 (3114.674) | 32.103 (1915.946) | 1.464 (2.198) | -0.110 (0.020) | 1.629 |
| Q-T1S2 | G0 | 64 | 0.0102 | 0.932 0.904 0.874 0.847 0.823 0.805 0.791 | 1.336 1.321 1.271 1.191 1.113 1.062 1.035 | -0.096 (0.001) | -0.130 / -0.116 | 1.114 (0.012) | 1.091 (0.016) | 0.827 (0.018) | 0.103 (0.002) | 1.193 |
| Q-T1S2 | GT | 64 | 0.0102 | 0.840 0.811 0.786 0.769 0.760 0.757 0.759 | 1.336 1.321 1.271 1.191 1.113 1.062 1.035 | -0.054 (0.001) | -0.130 / -0.109 | 1.191 (0.016) | 1.138 (0.027) | 0.775 (0.030) | 0.088 (0.002) | 1.232 |
| Q-T1S2 | GF | 64 | 0.0271 | 0.998 0.998 0.998 0.999 1.000 1.001 1.002 | 1.345 1.329 1.275 1.192 1.112 1.061 1.035 | 0.003 (0.006) | -0.140 / -0.001 | 208.826 (5532.478) | 18.110 (112.089) | 1.614 (1.681) | -0.139 (0.021) | 1.738 |
| Q-T1S3 | G0 | 64 | 0.0327 | 0.812 0.786 0.760 0.737 0.720 0.707 0.698 | 1.380 1.360 1.296 1.197 1.109 1.058 1.033 | -0.077 (0.001) | -0.166 / -0.104 | 1.592 (0.050) | 1.734 (0.069) | 0.849 (0.079) | 0.043 (0.005) | 1.269 |
| Q-T1S3 | GT | 64 | 0.0327 | 0.794 0.768 0.750 0.744 0.747 0.754 0.762 | 1.380 1.360 1.296 1.197 1.109 1.058 1.033 | -0.018 (0.001) | -0.166 / -0.087 | 1.904 (0.056) | 1.934 (0.115) | 0.720 (0.073) | 0.008 (0.005) | 1.324 |
| Q-T1S3 | GF | 64 | 0.0377 | 1.000 0.999 0.999 0.999 1.000 1.001 1.001 | 1.395 1.373 1.305 1.201 1.110 1.058 1.033 | 0.001 (0.007) | -0.179 / -0.003 | 56.145 (973.613) | 8.815 (18.477) | 1.250 (1.396) | -0.173 (0.024) | 1.905 |
| Q-T2S1 | G0 | 64 | 0.0403 | 0.722 0.711 0.700 0.691 0.683 0.676 0.671 | 1.237 1.225 1.190 1.137 1.084 1.048 1.028 | -0.034 (0.001) | -0.094 / -0.044 | 2.141 (0.146) | 3.648 (0.482) | 2.342 (0.837) | -0.006 (0.006) | 1.184 |
| Q-T2S1 | GT | 64 | 0.0403 | 0.686 0.673 0.665 0.664 0.668 0.675 0.684 | 1.237 1.225 1.190 1.137 1.084 1.048 1.028 | 0.000 (0.001) | -0.094 / -0.042 | 2.244 (0.155) | 3.769 (0.648) | 2.235 (0.969) | -0.010 (0.006) | 1.193 |
| Q-T2S1 | GF | 64 | 0.0406 | 0.997 0.997 0.997 0.998 0.999 1.001 1.002 | 1.246 1.233 1.194 1.138 1.084 1.047 1.027 | 0.003 (0.005) | -0.103 / 0.001 | -171.076 (115858.578) | 30.396 (489.328) | 3.077 (5.011) | -0.104 (0.015) | 1.489 |
| Q-T2S2 | G0 | 64 | 0.0405 | 0.722 0.708 0.695 0.684 0.675 0.669 0.664 | 1.243 1.228 1.184 1.121 1.067 1.036 1.020 | -0.039 (0.001) | -0.119 / -0.054 | 2.188 (0.148) | 3.679 (0.458) | 2.247 (0.744) | -0.010 (0.008) | 1.241 |
| Q-T2S2 | GT | 64 | 0.0405 | 0.699 0.684 0.678 0.681 0.690 0.701 0.712 | 1.243 1.228 1.184 1.121 1.067 1.036 1.020 | 0.013 (0.001) | -0.119 / -0.042 | 2.849 (0.193) | 5.137 (1.091) | 1.993 (0.701) | -0.035 (0.008) | 1.272 |
| Q-T2S2 | GF | 64 | 0.0408 | 0.998 0.998 0.998 0.999 1.000 1.001 1.003 | 1.252 1.235 1.187 1.121 1.066 1.035 1.019 | 0.003 (0.005) | -0.130 / 0.000 | -360.011 (8871.480) | 25.267 (987.939) | 3.091 (4.641) | -0.130 (0.018) | 1.633 |
| Q-T2S3 | G0 | 64 | 0.0052 | 0.929 0.897 0.865 0.838 0.817 0.801 0.790 | 1.270 1.253 1.200 1.125 1.066 1.034 1.019 | -0.093 (0.001) | -0.141 / -0.126 | 1.120 (0.010) | 1.101 (0.015) | 0.926 (0.017) | 0.111 (0.002) | 1.247 |
| Q-T2S3 | GT | 64 | 0.0052 | 0.815 0.784 0.761 0.751 0.751 0.755 0.762 | 1.270 1.253 1.200 1.125 1.066 1.034 1.019 | -0.032 (0.001) | -0.141 / -0.107 | 1.327 (0.018) | 1.241 (0.033) | 0.897 (0.038) | 0.072 (0.003) | 1.314 |
| Q-T2S3 | GF | 64 | 0.0216 | 0.998 0.997 0.997 0.998 1.000 1.001 1.002 | 1.278 1.259 1.202 1.125 1.065 1.033 1.018 | 0.003 (0.006) | -0.152 / -0.001 | 126.096 (67168.481) | 27.176 (220.580) | 3.151 (4.221) | -0.149 (0.020) | 1.814 |
| Q-T3S1 | G0 | 64 | 0.0130 | 0.854 0.834 0.815 0.797 0.781 0.768 0.758 | 1.284 1.271 1.231 1.169 1.105 1.060 1.035 | -0.065 (0.001) | -0.106 / -0.079 | 1.349 (0.026) | 1.265 (0.029) | 0.891 (0.028) | 0.051 (0.002) | 1.177 |
| Q-T3S1 | GT | 64 | 0.0130 | 0.758 0.737 0.721 0.710 0.706 0.707 0.710 | 1.284 1.271 1.231 1.169 1.105 1.060 1.035 | -0.031 (0.001) | -0.106 / -0.074 | 1.427 (0.034) | 1.324 (0.052) | 0.902 (0.057) | 0.043 (0.003) | 1.193 |
| Q-T3S1 | GF | 64 | 0.0226 | 0.999 0.998 0.997 0.998 0.998 0.999 0.999 | 1.291 1.277 1.234 1.170 1.105 1.060 1.034 | 0.000 (0.006) | -0.114 / -0.003 | 33.403 (151.431) | 15.422 (40.324) | 2.134 (2.234) | -0.107 (0.015) | 1.511 |
| Q-T3S2 | G0 | 64 | 0.0089 | 0.922 0.892 0.862 0.836 0.814 0.797 0.784 | 1.327 1.310 1.256 1.175 1.101 1.054 1.031 | -0.093 (0.001) | -0.141 / -0.119 | 1.192 (0.013) | 1.081 (0.014) | 0.816 (0.016) | 0.096 (0.002) | 1.221 |
| Q-T3S2 | GT | 64 | 0.0089 | 0.817 0.787 0.763 0.749 0.743 0.743 0.746 | 1.327 1.310 1.256 1.175 1.101 1.054 1.031 | -0.046 (0.001) | -0.141 / -0.107 | 1.320 (0.020) | 1.136 (0.027) | 0.772 (0.029) | 0.073 (0.003) | 1.262 |
| Q-T3S2 | GF | 64 | 0.0245 | 0.999 0.997 0.997 0.998 0.998 0.999 0.999 | 1.336 1.318 1.261 1.176 1.100 1.054 1.030 | 0.001 (0.006) | -0.151 / -0.003 | 50.576 (374.481) | 22.848 (132.958) | 2.306 (2.962) | -0.145 (0.019) | 1.681 |
| Q-T3S3 | G0 | 64 | 0.0042 | 1.096 1.041 0.984 0.934 0.896 0.869 0.850 | 1.431 1.407 1.325 1.202 1.105 1.054 1.031 | -0.167 (0.001) | -0.213 / -0.223 | 0.955 (0.004) | 0.880 (0.005) | 0.748 (0.008) | 0.233 (0.001) | 1.284 |
| Q-T3S3 | GT | 64 | 0.0042 | 0.979 0.923 0.874 0.842 0.828 0.823 0.824 | 1.431 1.407 1.325 1.202 1.105 1.054 1.031 | -0.101 (0.001) | -0.213 / -0.209 | 1.019 (0.006) | 0.864 (0.009) | 0.619 (0.010) | 0.205 (0.002) | 1.398 |
| Q-T3S3 | GF | 64 | 0.0301 | 1.000 0.999 0.999 1.000 1.001 1.002 1.002 | 1.440 1.414 1.328 1.203 1.104 1.053 1.030 | 0.002 (0.007) | -0.224 / -0.002 | 99.192 (7177.044) | 79.845 (1425.949) | 1.944 (3.055) | -0.220 (0.028) | 2.130 |
| (all) | GM | 48 | 0.0067 | 1.390 1.391 1.390 1.388 1.384 1.379 1.373 | 1.635 1.606 1.574 1.538 1.502 1.470 1.443 | -0.012 (0.001) | -0.122 / 0.001 | -102.338 (482.471) | 7.054 (0.779) | 3.980 (0.334) | -0.124 (0.010) | 1.219 |
| (all) | GM_floor_G01 | 48 | 0.0114 | 1.461 1.467 1.475 1.486 1.497 1.510 1.522 | 1.635 1.606 1.574 1.538 1.502 1.470 1.443 | 0.041 (0.003) | -0.122 / 0.028 | -4.361 (0.767) | 22.395 (82.141) | 3.069 (1.782) | -0.178 (0.009) | 1.259 |
| Q-T1S2 (tail) | G0 | 64 | 0.0329 | 0.907 0.865 0.829 0.800 0.779 0.763 0.753 | 1.452 1.423 1.335 1.216 1.118 1.061 1.034 | -0.102 (0.001) | -0.234 / -0.157 | 1.495 (0.031) | 1.300 (0.038) | 0.819 (0.085) | 0.079 (0.005) | 1.368 |
| Q-T1S2 (tail) | GT | 64 | 0.0329 | 0.811 0.769 0.742 0.730 0.729 0.733 0.739 | 1.452 1.423 1.335 1.216 1.118 1.061 1.034 | -0.043 (0.002) | -0.234 / -0.136 | 1.717 (0.046) | 1.475 (0.076) | 0.980 (0.202) | 0.039 (0.007) | 1.430 |
| Q-T1S2 (tail) | GF | 64 | 0.0329 | 0.991 0.991 0.993 0.995 0.996 0.996 0.996 | 1.416 1.385 1.298 1.186 1.099 1.051 1.028 | 0.004 (0.005) | -0.235 / 0.004 | -57.355 (179.172) | -10.812 (14.275) | -1.494 (1.896) | -0.244 (0.022) | 1.619 |
| Q-T3S1 (tail) | G0 | 64 | 0.0340 | 0.845 0.817 0.793 0.773 0.758 0.746 0.737 | 1.384 1.361 1.294 1.201 1.118 1.065 1.036 | -0.072 (0.001) | -0.180 / -0.104 | 1.741 (0.058) | 1.613 (0.061) | 1.081 (0.149) | 0.027 (0.006) | 1.288 |
| Q-T3S1 (tail) | GT | 64 | 0.0340 | 0.746 0.718 0.699 0.690 0.689 0.693 0.698 | 1.384 1.361 1.294 1.201 1.118 1.065 1.036 | -0.029 (0.002) | -0.180 / -0.095 | 1.905 (0.080) | 1.826 (0.145) | 1.993 (1.307) | 0.009 (0.008) | 1.309 |
| Q-T3S1 (tail) | GF | 64 | 0.0340 | 0.992 0.992 0.993 0.994 0.996 0.997 0.997 | 1.348 1.324 1.259 1.173 1.099 1.054 1.030 | 0.004 (0.005) | -0.178 / 0.002 | -91.149 (2385.637) | -17.106 (36.110) | -2.215 (3.080) | -0.182 (0.018) | 1.437 |

**Measurement (truth).**
- P-F slope: G0 gives -0.034 to -0.167, GT -0.101 to +0.013, GF and GM PAPER2 about 0 (0.000 to
  0.004, and -0.012), and the beta-floor variant +0.041.
- GT is as flat as GF in Q-T2S1 (0.000) and Q-T2S2 (+0.013), so the flatness test cannot separate them
  there.
- On [64, 1024] every ratio is far from its [8, 64] value and every SE is large.

## 2. Panels and rho*

200 panels per NASDAQ cell and generator, each sized like the cell's confirmation set (17 to 252 stocks x
8 days). Each stock's 8 day lengths are copied from a random exploration stock of the cell. Signs follow
the 3a fGn + AR(1) model (GF: herding signs refitted, 10 % informed coins). Sizes are log-normal, and the
moments' sizes are clipped at each stock's own 99.9th percentile (R2); the impact uses raw sizes. The
estimator ran without a ceiling.

**Choices.**
- For a T/N bin b, every increment uses only the stocks whose T_j+1 / N_s falls in b, so an estimate
  belongs to one band.
- Bias = mean over panels minus truth. Coverage = share of 90 % intervals (the wider rule) containing the
  truth; an undefined interval counts as not covering.
- An element is evaluable when at least 50 % of panels give an estimate.
- rho* = the largest band edge whose band, and every evaluable band below it, passes for every evaluable
  element (all a <= 1.5 and every range).
- The ratio under GF is left out (no finite truth).

Per cell, generator and estimator ("-" = the first evaluable band already fails):

| sign model | cell | G0 zeta_V | G0 zeta_K | G0 pf | G0 ratio | G0 vk_diff | G0 Rratio | GT zeta_V | GT zeta_K | GT pf | GT ratio | GT vk_diff | GT Rratio | GF zeta_V | GF zeta_K | GF pf | GF ratio | GF vk_diff | GF Rratio | rho* (min) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| main | Q-T1S1 | 1/10 | 1/80 | 1/5 | 1/80 | 1/5 | 1/80 | 1/80 | 1/80 | 1/5 | 1/80 | 1/5 | 1/5 | 1/80 | - | - | undef | - | - | - |
| main | Q-T1S2 | - | - | 1/40 | - | - | - | 1/80 | 1/80 | 1/80 | 1/80 | 1/80 | 1/80 | - | - | - | undef | - | - | - |
| main | Q-T1S3 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | undef | - | - | - |
| main | Q-T2S1 | - | - | 1/20 | - | - | 1/80 | - | 1/20 | 1/20 | - | - | 1/80 | - | - | - | undef | - | 1/80 | - |
| main | Q-T2S2 | - | - | 1/20 | - | - | 1/40 | - | - | - | - | - | 1/80 | - | - | 1/20 | undef | - | - | - |
| main | Q-T2S3 | - | - | - | - | - | 1/10 | - | - | - | - | - | 1/10 | - | - | - | undef | - | 1/20 | - |
| main | Q-T3S1 | 1/20 | 1/40 | 1/20 | - | 1/20 | 1/20 | 1/10 | 1/20 | 1/10 | - | 1/10 | 1/20 | 1/20 | - | - | undef | - | 1/80 | - |
| main | Q-T3S2 | - | - | - | - | - | 1/20 | - | - | - | - | - | 1/80 | - | - | - | undef | - | - | - |
| main | Q-T3S3 | - | - | 1/5 | - | - | 1/40 | - | - | - | - | - | - | - | - | - | undef | - | - | - |
| tail | Q-T1S2 | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | undef | - | - | - |
| tail | Q-T3S1 | - | - | - | - | - | - | - | - | 1/20 | - | - | - | - | - | - | undef | - | 1/80 | - |

The same with looser pass rules (sensitivity, not decisions):

| rule | sign model | cell | G0 | GT | GF | rho* (min over generators and estimators) |
|---|---|---|---|---|---|---|
| draft | main | Q-T1S1 | 1/80 | 1/80 | - | - |
| draft | main | Q-T1S2 | - | 1/80 | - | - |
| draft | main | Q-T1S3 | - | - | - | - |
| draft | main | Q-T2S1 | - | - | - | - |
| draft | main | Q-T2S2 | - | - | - | - |
| draft | main | Q-T2S3 | - | - | - | - |
| draft | main | Q-T3S1 | - | - | - | - |
| draft | main | Q-T3S2 | - | - | - | - |
| draft | main | Q-T3S3 | - | - | - | - |
| draft | tail | Q-T1S2 | - | - | - | - |
| draft | tail | Q-T3S1 | - | - | - | - |
| bias only | main | Q-T1S1 | 1/80 | 1/80 | - | - |
| bias only | main | Q-T1S2 | - | 1/40 | - | - |
| bias only | main | Q-T1S3 | - | - | - | - |
| bias only | main | Q-T2S1 | - | - | - | - |
| bias only | main | Q-T2S2 | - | - | - | - |
| bias only | main | Q-T2S3 | - | - | - | - |
| bias only | main | Q-T3S1 | - | - | - | - |
| bias only | main | Q-T3S2 | - | - | - | - |
| bias only | main | Q-T3S3 | - | - | - | - |
| bias only | tail | Q-T1S2 | - | - | - | - |
| bias only | tail | Q-T3S1 | - | - | - | - |
| bias + coverage >= 0.81 | main | Q-T1S1 | 1/80 | 1/80 | - | - |
| bias + coverage >= 0.81 | main | Q-T1S2 | - | 1/80 | - | - |
| bias + coverage >= 0.81 | main | Q-T1S3 | - | - | - | - |
| bias + coverage >= 0.81 | main | Q-T2S1 | - | - | - | - |
| bias + coverage >= 0.81 | main | Q-T2S2 | - | - | - | - |
| bias + coverage >= 0.81 | main | Q-T2S3 | - | - | - | - |
| bias + coverage >= 0.81 | main | Q-T3S1 | - | - | - | - |
| bias + coverage >= 0.81 | main | Q-T3S2 | - | - | - | - |
| bias + coverage >= 0.81 | main | Q-T3S3 | - | - | - | - |
| bias + coverage >= 0.81 | tail | Q-T1S2 | - | - | - | - |
| bias + coverage >= 0.81 | tail | Q-T3S1 | - | - | - | - |

**Choice for steps 3 and 4:** G0 and GT only, ratio left out:

| cell | rho*, draft bounds, G0 + GT, no ratio: main | tail | the same with coverage >= 0.81: main | tail |
|---|---|---|---|---|
| Q-T1S1 | 1/80 |  | 1/20 |  |
| Q-T1S2 | - | - | 1/80 | - |
| Q-T1S3 | - |  | - |  |
| Q-T2S1 | - |  | - |  |
| Q-T2S2 | - |  | - |  |
| Q-T2S3 | - |  | - |  |
| Q-T3S1 | 1/40 | - | 1/40 | - |
| Q-T3S2 | - |  | - |  |
| Q-T3S3 | - |  | - |  |

Where it fails, by generator, estimator and band (pooled over cells):

| gen | estimator | T/N band upper edge | elements evaluated | share with abs bias >= 0.02 | share with coverage < 0.85 | median coverage | 10th pct coverage |
|---|---|---|---|---|---|---|---|
| G0 | zeta_V | 1/80 | 14 | 0.00 | 0.29 | 0.912 | 0.835 |
| G0 | zeta_V | 1/40 | 21 | 0.10 | 0.48 | 0.875 | 0.795 |
| G0 | zeta_V | 1/20 | 63 | 0.13 | 0.43 | 0.860 | 0.722 |
| G0 | zeta_V | 1/10 | 77 | 0.18 | 0.34 | 0.895 | 0.657 |
| G0 | zeta_V | 1/5 | 105 | 0.27 | 0.55 | 0.835 | 0.376 |
| G0 | zeta_K | 1/80 | 14 | 0.00 | 0.29 | 0.895 | 0.833 |
| G0 | zeta_K | 1/40 | 21 | 0.05 | 0.19 | 0.965 | 0.795 |
| G0 | zeta_K | 1/20 | 63 | 0.37 | 0.44 | 0.850 | 0.756 |
| G0 | zeta_K | 1/10 | 77 | 0.61 | 0.19 | 0.945 | 0.808 |
| G0 | zeta_K | 1/5 | 105 | 0.76 | 0.17 | 0.960 | 0.810 |
| G0 | pf | 1/80 | 2 | 0.00 | 0.00 | 0.905 | 0.881 |
| G0 | pf | 1/40 | 3 | 0.00 | 0.00 | 0.960 | 0.884 |
| G0 | pf | 1/20 | 9 | 0.11 | 0.22 | 0.870 | 0.789 |
| G0 | pf | 1/10 | 11 | 0.55 | 0.27 | 0.920 | 0.805 |
| G0 | pf | 1/5 | 15 | 0.73 | 0.33 | 0.875 | 0.784 |
| G0 | ratio | 1/80 | 2 | 0.50 | 0.50 | 0.880 | 0.824 |
| G0 | ratio | 1/40 | 3 | 1.00 | 0.33 | 0.960 | 0.860 |
| G0 | ratio | 1/20 | 9 | 1.00 | 0.22 | 0.885 | 0.795 |
| G0 | ratio | 1/10 | 11 | 1.00 | 0.00 | 0.930 | 0.880 |
| G0 | ratio | 1/5 | 15 | 1.00 | 0.27 | 0.945 | 0.807 |
| G0 | vk_diff | 1/80 | 2 | 0.00 | 0.50 | 0.883 | 0.845 |
| G0 | vk_diff | 1/40 | 3 | 0.33 | 0.33 | 0.960 | 0.856 |
| G0 | vk_diff | 1/20 | 9 | 0.78 | 0.22 | 0.885 | 0.775 |
| G0 | vk_diff | 1/10 | 11 | 0.73 | 0.00 | 0.950 | 0.855 |
| G0 | vk_diff | 1/5 | 15 | 0.67 | 0.13 | 0.975 | 0.825 |
| G0 | Rratio | 1/80 | 10 | 0.10 | 0.10 | 0.893 | 0.847 |
| G0 | Rratio | 1/40 | 12 | 0.17 | 0.33 | 0.897 | 0.779 |
| G0 | Rratio | 1/20 | 19 | 0.26 | 0.37 | 0.910 | 0.779 |
| G0 | Rratio | 1/10 | 20 | 0.60 | 0.20 | 0.895 | 0.732 |
| G0 | Rratio | 1/5 | 23 | 0.74 | 0.17 | 0.885 | 0.741 |
| GT | zeta_V | 1/80 | 14 | 0.00 | 0.00 | 0.922 | 0.877 |
| GT | zeta_V | 1/40 | 21 | 0.14 | 0.48 | 0.850 | 0.725 |
| GT | zeta_V | 1/20 | 63 | 0.13 | 0.51 | 0.845 | 0.725 |
| GT | zeta_V | 1/10 | 77 | 0.18 | 0.36 | 0.890 | 0.650 |
| GT | zeta_V | 1/5 | 105 | 0.26 | 0.61 | 0.830 | 0.404 |
| GT | zeta_K | 1/80 | 14 | 0.00 | 0.00 | 0.930 | 0.872 |
| GT | zeta_K | 1/40 | 21 | 0.00 | 0.29 | 0.970 | 0.820 |
| GT | zeta_K | 1/20 | 63 | 0.10 | 0.33 | 0.915 | 0.795 |
| GT | zeta_K | 1/10 | 77 | 0.27 | 0.06 | 0.960 | 0.855 |
| GT | zeta_K | 1/5 | 105 | 0.60 | 0.23 | 0.955 | 0.790 |
| GT | pf | 1/80 | 2 | 0.00 | 0.00 | 0.938 | 0.907 |
| GT | pf | 1/40 | 3 | 0.00 | 0.33 | 0.970 | 0.838 |
| GT | pf | 1/20 | 9 | 0.22 | 0.44 | 0.865 | 0.788 |
| GT | pf | 1/10 | 11 | 0.55 | 0.09 | 0.960 | 0.865 |
| GT | pf | 1/5 | 15 | 0.80 | 0.20 | 0.940 | 0.744 |
| GT | ratio | 1/80 | 2 | 0.00 | 0.00 | 0.943 | 0.908 |
| GT | ratio | 1/40 | 3 | 0.67 | 0.33 | 0.970 | 0.866 |
| GT | ratio | 1/20 | 9 | 1.00 | 0.22 | 0.880 | 0.818 |
| GT | ratio | 1/10 | 11 | 1.00 | 0.00 | 0.975 | 0.855 |
| GT | ratio | 1/5 | 15 | 1.00 | 0.33 | 0.885 | 0.759 |
| GT | vk_diff | 1/80 | 2 | 0.00 | 0.00 | 0.950 | 0.918 |
| GT | vk_diff | 1/40 | 3 | 0.33 | 0.33 | 0.990 | 0.870 |
| GT | vk_diff | 1/20 | 9 | 0.33 | 0.22 | 0.900 | 0.819 |
| GT | vk_diff | 1/10 | 11 | 0.36 | 0.00 | 0.990 | 0.885 |
| GT | vk_diff | 1/5 | 15 | 0.67 | 0.13 | 0.995 | 0.825 |
| GT | Rratio | 1/80 | 10 | 0.10 | 0.10 | 0.910 | 0.857 |
| GT | Rratio | 1/40 | 12 | 0.25 | 0.25 | 0.917 | 0.827 |
| GT | Rratio | 1/20 | 19 | 0.26 | 0.16 | 0.945 | 0.820 |
| GT | Rratio | 1/10 | 20 | 0.60 | 0.15 | 0.923 | 0.809 |
| GT | Rratio | 1/5 | 23 | 0.74 | 0.17 | 0.935 | 0.816 |
| GF | zeta_V | 1/80 | 14 | 0.00 | 0.43 | 0.912 | 0.820 |
| GF | zeta_V | 1/40 | 21 | 0.14 | 0.38 | 0.875 | 0.750 |
| GF | zeta_V | 1/20 | 63 | 0.14 | 0.57 | 0.840 | 0.697 |
| GF | zeta_V | 1/10 | 77 | 0.17 | 0.35 | 0.860 | 0.624 |
| GF | zeta_V | 1/5 | 105 | 0.30 | 0.63 | 0.820 | 0.344 |
| GF | zeta_K | 1/80 | 14 | 1.00 | 0.14 | 0.910 | 0.848 |
| GF | zeta_K | 1/40 | 21 | 1.00 | 0.14 | 0.985 | 0.840 |
| GF | zeta_K | 1/20 | 63 | 0.97 | 0.44 | 0.970 | 0.665 |
| GF | zeta_K | 1/10 | 77 | 0.82 | 0.31 | 0.965 | 0.740 |
| GF | zeta_K | 1/5 | 105 | 0.87 | 0.44 | 0.885 | 0.447 |
| GF | pf | 1/80 | 2 | 1.00 | 0.50 | 0.902 | 0.852 |
| GF | pf | 1/40 | 3 | 1.00 | 0.33 | 0.985 | 0.797 |
| GF | pf | 1/20 | 9 | 0.44 | 0.44 | 0.935 | 0.449 |
| GF | pf | 1/10 | 11 | 0.91 | 0.55 | 0.840 | 0.570 |
| GF | pf | 1/5 | 12 | 0.83 | 0.42 | 0.903 | 0.388 |
| GF | vk_diff | 1/80 | 2 | 1.00 | 0.00 | 0.920 | 0.868 |
| GF | vk_diff | 1/40 | 3 | 0.67 | 0.33 | 0.985 | 0.801 |
| GF | vk_diff | 1/20 | 9 | 1.00 | 0.44 | 0.950 | 0.558 |
| GF | vk_diff | 1/10 | 11 | 1.00 | 0.27 | 0.960 | 0.710 |
| GF | vk_diff | 1/5 | 15 | 1.00 | 0.47 | 0.910 | 0.332 |
| GF | Rratio | 1/80 | 10 | 0.70 | 0.00 | 0.975 | 0.860 |
| GF | Rratio | 1/40 | 12 | 0.83 | 0.08 | 0.978 | 0.861 |
| GF | Rratio | 1/20 | 19 | 0.68 | 0.16 | 0.985 | 0.841 |
| GF | Rratio | 1/10 | 20 | 0.90 | 0.15 | 0.975 | 0.841 |
| GF | Rratio | 1/5 | 23 | 1.00 | 0.13 | 0.975 | 0.844 |

Tail-lag sign model (Q-T1S2, Q-T3S1), same layout:

| gen | estimator | T/N band upper edge | elements evaluated | share with abs bias >= 0.02 | share with coverage < 0.85 | median coverage | 10th pct coverage |
|---|---|---|---|---|---|---|---|
| G0 | zeta_V | 1/80 | 7 | 0.57 | 1.00 | 0.665 | 0.623 |
| G0 | zeta_V | 1/40 | 14 | 0.57 | 0.79 | 0.667 | 0.396 |
| G0 | zeta_V | 1/20 | 21 | 0.62 | 0.76 | 0.510 | 0.020 |
| G0 | zeta_V | 1/10 | 28 | 0.71 | 0.79 | 0.435 | 0.039 |
| G0 | zeta_V | 1/5 | 35 | 0.74 | 0.89 | 0.345 | 0.005 |
| G0 | zeta_K | 1/80 | 7 | 0.71 | 1.00 | 0.635 | 0.553 |
| G0 | zeta_K | 1/40 | 14 | 0.93 | 0.79 | 0.738 | 0.404 |
| G0 | zeta_K | 1/20 | 21 | 0.95 | 0.67 | 0.685 | 0.115 |
| G0 | zeta_K | 1/10 | 28 | 0.89 | 0.46 | 0.855 | 0.304 |
| G0 | zeta_K | 1/5 | 35 | 0.91 | 0.54 | 0.830 | 0.447 |
| G0 | pf | 1/80 | 1 | 1.00 | 1.00 | 0.735 | 0.735 |
| G0 | pf | 1/40 | 2 | 1.00 | 0.50 | 0.778 | 0.675 |
| G0 | pf | 1/20 | 3 | 1.00 | 0.67 | 0.545 | 0.389 |
| G0 | pf | 1/10 | 4 | 1.00 | 0.50 | 0.748 | 0.355 |
| G0 | pf | 1/5 | 5 | 1.00 | 0.60 | 0.685 | 0.564 |
| G0 | ratio | 1/80 | 1 | 1.00 | 1.00 | 0.485 | 0.485 |
| G0 | ratio | 1/40 | 2 | 1.00 | 1.00 | 0.623 | 0.477 |
| G0 | ratio | 1/20 | 3 | 1.00 | 0.67 | 0.550 | 0.358 |
| G0 | ratio | 1/10 | 4 | 1.00 | 0.50 | 0.853 | 0.653 |
| G0 | ratio | 1/5 | 5 | 1.00 | 0.40 | 0.975 | 0.779 |
| G0 | vk_diff | 1/80 | 1 | 1.00 | 1.00 | 0.535 | 0.535 |
| G0 | vk_diff | 1/40 | 2 | 1.00 | 1.00 | 0.667 | 0.545 |
| G0 | vk_diff | 1/20 | 3 | 1.00 | 0.67 | 0.555 | 0.395 |
| G0 | vk_diff | 1/10 | 4 | 1.00 | 0.50 | 0.867 | 0.706 |
| G0 | vk_diff | 1/5 | 5 | 1.00 | 0.20 | 0.985 | 0.793 |
| G0 | Rratio | 1/80 | 3 | 1.00 | 1.00 | 0.430 | 0.110 |
| G0 | Rratio | 1/40 | 4 | 1.00 | 1.00 | 0.542 | 0.190 |
| G0 | Rratio | 1/20 | 5 | 1.00 | 1.00 | 0.570 | 0.074 |
| G0 | Rratio | 1/10 | 6 | 1.00 | 0.83 | 0.667 | 0.070 |
| G0 | Rratio | 1/5 | 6 | 1.00 | 1.00 | 0.470 | 0.212 |
| GT | zeta_V | 1/80 | 7 | 0.57 | 1.00 | 0.680 | 0.613 |
| GT | zeta_V | 1/40 | 14 | 0.57 | 0.79 | 0.660 | 0.448 |
| GT | zeta_V | 1/20 | 21 | 0.62 | 0.76 | 0.500 | 0.040 |
| GT | zeta_V | 1/10 | 28 | 0.68 | 0.82 | 0.475 | 0.059 |
| GT | zeta_V | 1/5 | 35 | 0.74 | 0.89 | 0.435 | 0.005 |
| GT | zeta_K | 1/80 | 7 | 0.43 | 0.57 | 0.835 | 0.782 |
| GT | zeta_K | 1/40 | 14 | 0.43 | 0.29 | 0.922 | 0.789 |
| GT | zeta_K | 1/20 | 21 | 0.57 | 0.33 | 0.960 | 0.755 |
| GT | zeta_K | 1/10 | 28 | 0.64 | 0.04 | 0.982 | 0.875 |
| GT | zeta_K | 1/5 | 35 | 0.66 | 0.43 | 0.950 | 0.669 |
| GT | pf | 1/80 | 1 | 1.00 | 1.00 | 0.785 | 0.785 |
| GT | pf | 1/40 | 2 | 0.50 | 0.50 | 0.883 | 0.841 |
| GT | pf | 1/20 | 3 | 0.67 | 0.33 | 0.870 | 0.782 |
| GT | pf | 1/10 | 4 | 1.00 | 0.00 | 0.925 | 0.883 |
| GT | pf | 1/5 | 5 | 1.00 | 0.40 | 0.910 | 0.669 |
| GT | ratio | 1/80 | 1 | 1.00 | 1.00 | 0.735 | 0.735 |
| GT | ratio | 1/40 | 2 | 1.00 | 0.50 | 0.907 | 0.857 |
| GT | ratio | 1/20 | 3 | 1.00 | 0.33 | 0.985 | 0.829 |
| GT | ratio | 1/10 | 4 | 1.00 | 0.00 | 0.970 | 0.924 |
| GT | ratio | 1/5 | 5 | 1.00 | 0.40 | 0.940 | 0.708 |
| GT | vk_diff | 1/80 | 1 | 1.00 | 1.00 | 0.775 | 0.775 |
| GT | vk_diff | 1/40 | 2 | 1.00 | 0.00 | 0.915 | 0.871 |
| GT | vk_diff | 1/20 | 3 | 0.67 | 0.33 | 0.990 | 0.842 |
| GT | vk_diff | 1/10 | 4 | 0.75 | 0.00 | 0.995 | 0.941 |
| GT | vk_diff | 1/5 | 5 | 0.60 | 0.20 | 0.990 | 0.724 |
| GT | Rratio | 1/80 | 3 | 1.00 | 1.00 | 0.650 | 0.374 |
| GT | Rratio | 1/40 | 4 | 1.00 | 0.50 | 0.810 | 0.391 |
| GT | Rratio | 1/20 | 5 | 1.00 | 0.60 | 0.740 | 0.368 |
| GT | Rratio | 1/10 | 6 | 1.00 | 0.50 | 0.853 | 0.435 |
| GT | Rratio | 1/5 | 6 | 1.00 | 0.67 | 0.770 | 0.547 |
| GF | zeta_V | 1/80 | 7 | 0.57 | 0.86 | 0.715 | 0.637 |
| GF | zeta_V | 1/40 | 14 | 0.50 | 0.79 | 0.660 | 0.409 |
| GF | zeta_V | 1/20 | 21 | 0.62 | 0.76 | 0.535 | 0.020 |
| GF | zeta_V | 1/10 | 28 | 0.68 | 0.75 | 0.488 | 0.030 |
| GF | zeta_V | 1/5 | 35 | 0.71 | 0.77 | 0.395 | 0.000 |
| GF | zeta_K | 1/80 | 7 | 1.00 | 0.00 | 0.890 | 0.874 |
| GF | zeta_K | 1/40 | 14 | 1.00 | 0.50 | 0.887 | 0.777 |
| GF | zeta_K | 1/20 | 21 | 1.00 | 0.33 | 0.970 | 0.610 |
| GF | zeta_K | 1/10 | 28 | 0.86 | 0.25 | 0.972 | 0.680 |
| GF | zeta_K | 1/5 | 33 | 0.82 | 0.48 | 0.905 | 0.388 |
| GF | pf | 1/80 | 1 | 1.00 | 0.00 | 0.850 | 0.850 |
| GF | pf | 1/40 | 2 | 1.00 | 0.50 | 0.847 | 0.749 |
| GF | pf | 1/20 | 3 | 0.67 | 0.33 | 0.980 | 0.520 |
| GF | pf | 1/10 | 4 | 1.00 | 0.25 | 0.883 | 0.640 |
| GF | pf | 1/5 | 4 | 0.75 | 0.50 | 0.823 | 0.433 |
| GF | vk_diff | 1/80 | 1 | 1.00 | 0.00 | 0.880 | 0.880 |
| GF | vk_diff | 1/40 | 2 | 1.00 | 0.50 | 0.870 | 0.774 |
| GF | vk_diff | 1/20 | 3 | 1.00 | 0.33 | 0.990 | 0.570 |
| GF | vk_diff | 1/10 | 4 | 1.00 | 0.25 | 0.948 | 0.694 |
| GF | vk_diff | 1/5 | 4 | 1.00 | 0.25 | 0.912 | 0.726 |
| GF | Rratio | 1/80 | 3 | 0.33 | 0.00 | 0.970 | 0.886 |
| GF | Rratio | 1/40 | 4 | 0.50 | 0.25 | 0.988 | 0.863 |
| GF | Rratio | 1/20 | 5 | 0.40 | 0.20 | 0.975 | 0.867 |
| GF | Rratio | 1/10 | 6 | 0.83 | 0.17 | 0.982 | 0.860 |
| GF | Rratio | 1/5 | 6 | 0.83 | 0.00 | 0.965 | 0.922 |

**Measurement.**
- **Bias.** For G0 and GT, the share of elements with |bias| >= 0.02 rises with T/N, for zeta_K from 0
  to 0.05 at T/N <= 1/40 up to 0.6 to 0.76 at 1/5. zeta_V's bias is smaller but stays at -0.02 to -0.03
  on [8, 64] in Q-T1S3 and Q-T2S1 even in small bands. Under GF, zeta_K's bias is >= 0.02 in nearly every
  element (the exclusion selection).
- **Coverage.** At the smallest bands, median coverage is 0.85 to 0.97, but the 10th percentile is 0.72
  to 0.91. One element's coverage over 200 panels has a Monte Carlo SE of 0.021, so part of the "every
  element >= 0.85" failure is noise. Below that, R_{1/2} / R_0 reaches coverage 0.24 at band 1/20 in
  Q-T1S3: its intervals are too narrow.
- **Stability.** rho* is not stable across cells: per estimator it runs from none to 1/5.
- **Tail sign model.** With the 3b tail-lag model, rho* moves down: none in Q-T1S2 and Q-T3S1 under
  every rule, where the main model gave 1/80 and 1/40 for G0 and GT.

## 3. a_max threshold (coverage against Kish ESS)

Coverage at each cell's rho* (draft bounds, G0 and GT, no ratio; 1/80 where none). ESS threshold (smallest ESS with coverage >= 0.85 at and above it): V weights 694.9, K weights 694.9; with coverage >= 0.81: V 2.7, K 10.4. Points (generator x a x range): V 63, K 63. a_max below uses the 0.81 thresholds where the 0.85 ones do not exist.

| cell | rho* used | a_max (V weights, real median ESS) | a_max (K weights) | a_max (both) |
|---|---|---|---|---|
| Q-T1S1 | 1/80 | - | - | - |
| Q-T1S2 | - | - | - | - |
| Q-T1S3 | - | - | - | - |
| Q-T2S1 | - | - | - | - |
| Q-T2S2 | - | 0.00 | 0.00 | 0.00 |
| Q-T2S3 | - | - | - | - |
| Q-T3S1 | 1/40 | 0.25 | 0.50 | 0.25 |
| Q-T3S2 | - | - | - | - |
| Q-T3S3 | - | - | - | - |

| median ESS bin | V points | V min coverage | V mean coverage | K points | K min coverage | K mean coverage |
|---|---|---|---|---|---|---|
| 0 to 5 | 9 | 0.86 | 0.93 | 0 | - | - |
| 5 to 10 | 12 | 0.83 | 0.92 | 0 | - | - |
| 10 to 20 | 6 | 0.94 | 0.97 | 6 | 0.96 | 0.97 |
| 20 to 30 | 9 | 0.82 | 0.92 | 6 | 0.84 | 0.91 |
| 30 to 50 | 0 | - | - | 3 | 0.85 | 0.86 |
| 50 to 100 | 9 | 0.83 | 0.93 | 12 | 0.84 | 0.94 |
| 100 to 300 | 8 | 0.83 | 0.93 | 17 | 0.83 | 0.93 |
| 300 to 1e+09 | 10 | 0.81 | 0.91 | 19 | 0.82 | 0.94 |

**Measurement.**
- Coverage barely depends on the median stock-day ESS, from 5 to above 300.
- The 0.85 rule has a threshold only because the highest-ESS points pass (695). The 0.81 rule gives 2.7
  (V weights) and 10.4 (K weights).
- So the coverage shortfall is not an ESS problem. It comes from the interval itself, with 8 days and
  17 to 252 stocks.
- The a_max column uses the 0.85 threshold of 695. Only Q-T3S1 (a_max 0.25 on V weights, 0.5 on K
  weights) and Q-T2S2 (0) have a real median ESS that high; elsewhere even a = 0 falls short (median N
  below 695).

## 4. N_floor

Confirmation stocks remaining (N_min from `data/universe.tsv`, a descriptive column; no confirmation run
was read), the horizon ceiling T_max = rho x N_floor, and the confirmation stocks whose N_min reaches
each fitting range's top (N_min >= T2 / rho). The table is shown at rho* (G0 + GT, no ratio) where one
exists, and at 1/20 and 1/10 for reference (display only).

| cell | ceiling rho | confirm stocks | N_floor 50: stocks / T_max / ranges reached | N_floor 100: stocks / T_max / ranges reached | N_floor 200: stocks / T_max / ranges reached | N_floor 400: stocks / T_max / ranges reached | N_floor 800: stocks / T_max / ranges reached |
|---|---|---|---|---|---|---|---|
| Q-T1S1 | rho* = 1/80 | 161 | 108 / 1 / 0+0+0 | 79 / 1 / 0+0+0 | 51 / 2 / 0+0+0 | 38 / 5 / 0+0+0 | 29 / 10 / 0+0+0 |
| Q-T1S1 | 1/20 (reference) | 161 | 108 / 2 / 22+0+0 | 79 / 5 / 22+0+0 | 51 / 10 / 22+0+0 | 38 / 20 / 22+0+0 | 29 / 40 / 22+0+0 |
| Q-T1S1 | 1/10 (reference) | 161 | 108 / 5 / 34+4+0 | 79 / 10 / 34+4+0 | 51 / 20 / 34+4+0 | 38 / 40 / 34+4+0 | 29 / 80 / 29+4+0 |
| Q-T1S2 | rho*: none | 216 | - | - | - | - | - |
| Q-T1S2 | 1/20 (reference) | 216 | 173 / 2 / 23+2+0 | 141 / 5 / 23+2+0 | 114 / 10 / 23+2+0 | 81 / 20 / 23+2+0 | 42 / 40 / 23+2+0 |
| Q-T1S2 | 1/10 (reference) | 216 | 173 / 5 / 52+5+0 | 141 / 10 / 52+5+0 | 114 / 20 / 52+5+0 | 81 / 40 / 52+5+0 | 42 / 80 / 42+5+0 |
| Q-T1S3 | rho*: none | 252 | - | - | - | - | - |
| Q-T1S3 | 1/20 (reference) | 252 | 163 / 2 / 5+0+0 | 142 / 5 / 5+0+0 | 98 / 10 / 5+0+0 | 55 / 20 / 5+0+0 | 11 / 40 / 5+0+0 |
| Q-T1S3 | 1/10 (reference) | 252 | 163 / 5 / 26+0+0 | 142 / 10 / 26+0+0 | 98 / 20 / 26+0+0 | 55 / 40 / 26+0+0 | 11 / 80 / 11+0+0 |
| Q-T2S1 | rho*: none | 80 | - | - | - | - | - |
| Q-T2S1 | 1/20 (reference) | 80 | 76 / 2 / 17+2+0 | 72 / 5 / 17+2+0 | 63 / 10 / 17+2+0 | 38 / 20 / 17+2+0 | 28 / 40 / 17+2+0 |
| Q-T2S1 | 1/10 (reference) | 80 | 76 / 5 / 32+6+2 | 72 / 10 / 32+6+2 | 63 / 20 / 32+6+2 | 38 / 40 / 32+6+2 | 28 / 80 / 28+6+2 |
| Q-T2S2 | rho*: none | 50 | - | - | - | - | - |
| Q-T2S2 | 1/20 (reference) | 50 | 46 / 2 / 3+0+0 | 42 / 5 / 3+0+0 | 37 / 10 / 3+0+0 | 24 / 20 / 3+0+0 | 9 / 40 / 3+0+0 |
| Q-T2S2 | 1/10 (reference) | 50 | 46 / 5 / 12+0+0 | 42 / 10 / 12+0+0 | 37 / 20 / 12+0+0 | 24 / 40 / 12+0+0 | 9 / 80 / 9+0+0 |
| Q-T2S3 | rho*: none | 17 | - | - | - | - | - |
| Q-T2S3 | 1/20 (reference) | 17 | 16 / 2 / 0+0+0 | 14 / 5 / 0+0+0 | 11 / 10 / 0+0+0 | 5 / 20 / 0+0+0 | 0 / 40 / 0+0+0 |
| Q-T2S3 | 1/10 (reference) | 17 | 16 / 5 / 1+0+0 | 14 / 10 / 1+0+0 | 11 / 20 / 1+0+0 | 5 / 40 / 1+0+0 | 0 / 80 / 0+0+0 |
| Q-T3S1 | rho* = 1/40 | 50 | 48 / 1 / 5+1+0 | 46 / 2 / 5+1+0 | 43 / 5 / 5+1+0 | 38 / 10 / 5+1+0 | 26 / 20 / 5+1+0 |
| Q-T3S1 | 1/20 (reference) | 50 | 48 / 2 / 21+2+0 | 46 / 5 / 21+2+0 | 43 / 10 / 21+2+0 | 38 / 20 / 21+2+0 | 26 / 40 / 21+2+0 |
| Q-T3S1 | 1/10 (reference) | 50 | 48 / 5 / 31+5+1 | 46 / 10 / 31+5+1 | 43 / 20 / 31+5+1 | 38 / 40 / 31+5+1 | 26 / 80 / 26+5+1 |
| Q-T3S2 | rho*: none | 27 | - | - | - | - | - |
| Q-T3S2 | 1/20 (reference) | 27 | 26 / 2 / 6+1+0 | 26 / 5 / 6+1+0 | 25 / 10 / 6+1+0 | 18 / 20 / 6+1+0 | 11 / 40 / 6+1+0 |
| Q-T3S2 | 1/10 (reference) | 27 | 26 / 5 / 12+3+0 | 26 / 10 / 12+3+0 | 25 / 20 / 12+3+0 | 18 / 40 / 12+3+0 | 11 / 80 / 11+3+0 |
| Q-T3S3 | rho*: none | 23 | - | - | - | - | - |
| Q-T3S3 | 1/20 (reference) | 23 | 22 / 2 / 1+1+0 | 22 / 5 / 1+1+0 | 19 / 10 / 1+1+0 | 10 / 20 / 1+1+0 | 1 / 40 / 1+1+0 |
| Q-T3S3 | 1/10 (reference) | 23 | 22 / 5 / 5+1+0 | 22 / 10 / 5+1+0 | 19 / 20 / 5+1+0 | 10 / 40 / 5+1+0 | 1 / 80 / 1+1+0 |

ranges reached = remaining stocks whose N_min >= T2 / rho* for [8,64] + [16,256] + [64,1024].

sigma_l-hat at small N: N contiguous runs from each exploration stock-day with at least N runs (random
offset), pooled over the stock's days, minus its full-sample sigma_l (sizes only):

| cell | N = 50: bias (sd, stocks) | N = 100: bias (sd, stocks) | N = 200: bias (sd, stocks) | N = 400: bias (sd, stocks) | N = 800: bias (sd, stocks) |
|---|---|---|---|---|---|
| Q-T1S1 | -0.034 (0.153, 62) | 0.017 (0.066, 54) | -0.003 (0.105, 45) | -0.013 (0.132, 31) | 0.009 (0.103, 15) |
| Q-T1S2 | -0.014 (0.105, 91) | -0.002 (0.065, 89) | 0.010 (0.072, 83) | 0.021 (0.080, 70) | 0.031 (0.065, 57) |
| Q-T1S3 | -0.049 (0.170, 100) | -0.015 (0.080, 91) | 0.001 (0.101, 82) | -0.011 (0.129, 67) | 0.029 (0.068, 41) |
| Q-T2S1 | 0.028 (0.090, 34) | 0.020 (0.073, 32) | 0.006 (0.050, 31) | 0.026 (0.076, 26) | -0.005 (0.076, 17) |
| Q-T2S2 | -0.006 (0.106, 20) | 0.026 (0.078, 19) | 0.034 (0.037, 18) | 0.043 (0.034, 16) | 0.023 (0.051, 13) |
| Q-T2S3 | 0.008 (0.098, 7) | -0.013 (0.057, 7) | 0.035 (0.052, 6) | 0.030 (0.018, 6) | -0.089 (0.129, 6) |
| Q-T3S1 | 0.045 (0.111, 22) | 0.024 (0.076, 22) | 0.052 (0.111, 22) | 0.045 (0.065, 21) | 0.019 (0.096, 18) |
| Q-T3S2 | 0.010 (0.188, 11) | 0.063 (0.057, 11) | 0.010 (0.057, 11) | -0.026 (0.124, 11) | -0.049 (0.186, 9) |
| Q-T3S3 | -0.060 (0.211, 10) | -0.098 (0.123, 10) | -0.054 (0.154, 9) | -0.008 (0.148, 9) | -0.111 (0.145, 7) |

**Measurement.**
- [64, 1024] is out of reach for NASDAQ confirmation stocks at any ceiling up to 1/10: 0 to 2 stocks per
  cell. [16, 256] is reached by at most 2 stocks per cell at 1/20 and at most 6 at 1/10; [8, 64] by up to
  23 at 1/20 and up to 52 at 1/10.
- sigma_l-hat at N = 50 to 800 is unbiased within about +-0.065 in most cells (larger, -0.05 to -0.11, in
  Q-T3S3, the widest sizes), with a stock-to-stock SD of 0.04 to 0.2.

## 5. T_min and l_max

T_min = smallest T with P(Delta_T != 0) >= 0.5, pooled over the windows of the exploration stock-days
(mids alone). l_max = the largest power of 2 with 8 l_max <= the 25th percentile of N_min among
exploration stocks with N_min > 100 (provisional N_floor). beta-hat comes from G0 panels (exploration
day lengths with N_min > 100, 20 panels per cell), with the pooled sign ACF and response,
`propagator.calibrate_G`, and `fit_power` on [1, l_max] with l0 fixed at its true value 1 (Choice).

| cell | P(Delta_T != 0) at T = 1, 2, 4, 8, 16, 32, 64 | T_min | l_max (N_floor 100) | beta | beta-hat bias (sd) |
|---|---|---|---|---|---|
| Q-T1S1 | 0.60 0.79 0.89 0.95 0.97 0.98 0.99 | 1 | 16 | 0.365 | 0.042 (0.003) |
| Q-T1S2 | 0.65 0.83 0.92 0.96 0.97 0.98 0.99 | 1 | 16 | 0.405 | 0.042 (0.003) |
| Q-T1S3 | 0.68 0.85 0.93 0.96 0.98 0.99 0.99 | 1 | 16 | 0.480 | 0.052 (0.003) |
| Q-T2S1 | 0.52 0.69 0.81 0.88 0.92 0.94 0.96 | 1 | 16 | 0.495 | 0.033 (0.005) |
| Q-T2S2 | 0.53 0.70 0.81 0.88 0.92 0.94 0.96 | 1 | 32 | 0.495 | 0.029 (0.009) |
| Q-T2S3 | 0.55 0.71 0.82 0.88 0.92 0.95 0.96 | 1 | 32 | 0.365 | 0.020 (0.013) |
| Q-T3S1 | 0.46 0.61 0.72 0.81 0.86 0.90 0.93 | 2 | 32 | 0.420 | 0.035 (0.004) |
| Q-T3S2 | 0.38 0.53 0.65 0.75 0.83 0.88 0.91 | 2 | 16 | 0.395 | 0.041 (0.009) |
| Q-T3S3 | 0.30 0.44 0.57 0.68 0.77 0.83 0.88 | 4 | 32 | 0.345 | 0.035 (0.008) |

**Measurement.**
- T_min is 1 in the T1 and T2 cells, 2 in Q-T3S1 and Q-T3S2, and 4 in Q-T3S3: large-tick mids often do
  not move over a run.
- beta-hat at l_max = 16 or 32 is biased up by 0.020 to 0.052, with SD 0.003 to 0.013, so the bias is
  systematic.

## 6. Semi-synthetic generator (for Phase 3d; not used for rho*)

`panel3c.semi_synthetic`: the real exploration signs and sizes per stock-day pass through unchanged, and
the mid is the G0, GT or GF response to the raw sizes; the moments see R2-clipped sizes.
`phase3c_report.load_cell_flow` reads a cell's flow through the guard. Tests:
- the flow passes through bit-identical;
- GT with c = 0 equals G0;
- the response uses raw sizes while the moments see clipped ones.

Semi-synthetic smoke check on the real exploration flow of Q-T3S1 (22 stocks x 8 days):
- G0: flow identical True; runs 368203; share of runs clipped 0.0010
- GT: flow identical True; runs 368203; share of runs clipped 0.0010
- GF: flow identical True; runs 368203; share of runs clipped 0.0010

**Choice / caveat.** On real flow, GF's informed runs are a random 10 % of real runs, whose signs are not
independent of the rest of the flow. GF's exact zeta_K = 1 therefore does not carry over to the
semi-synthetic panels; Phase 3d must take this into account.

## 7. Tests and planted bugs

`test_pooled.py` has 8 checks and 6 planted bugs; `test_panel3c.py` has 7 checks and 6 planted bugs. All
12 are caught, each by its intended assertion (each mutant rerun outside pytest). Full suite: 297 passed,
28 skipped.

| planted bug | mutated function | caught by |
|---|---|---|
| logs of raw cross-stock averages | `pooled.core` | check_mean_of_logs |
| days and stocks weighted by run count | `pooled.core` | check_equal_weights (1.666 for 1.733) |
| a stock used above its ceiling | `pooled.increment_masks` | check_ceiling_and_matching (1.400 for 1.460) |
| non-matched slopes (per-horizon means differenced) | `pooled._increments` | check_ceiling_and_matching (1.402 for 1.460) |
| jackknife deleting stocks instead of days | `pooled.jackknife_weights` | check_jackknife_deletes_days |
| K_a(1) exclusion ignored | `pooled.normalize` | check_k1_exclusion |
| batched kernel off by one lag | `panel3c.batched_causal_conv` | check_batched_conv_equals_gen3 |
| semi-synthetic returns the clipped sizes as the flow | `panel3c.semi_synthetic` | check_flow_unchanged |
| GT bins sharing a beta convolved once each | `panel3c.response` | check_gt_c0_is_g0 |
| response on the clipped sizes | `panel3c.semi_synthetic` | check_response_uses_raw_sizes |
| R2 cap over the whole panel | `panel3c.stock_caps` | check_r2_caps_per_stock |
| day lengths drawn day by day across stocks | `panel3c.day_lengths` | check_day_lengths_copied |

Other checks cover:
- panel moments equal `moments.day_moments`, and the centred arm equals direct centred means;
- the increment weights equal the OLS slope of the chained curve;
- the weighted bootstrap reproduces `inference.cluster_bootstrap`'s replicates;
- batched signs match the arcsine-law C(l) of both sign models.

## 8. Wall time and peak memory

| job | wall | peak RAM | notes |
|---|---|---|---|
| truth (26 runs; 9 GT reused from 3b) | 10.7 min | main 1.16 GB; workers not sampled, about 0.3 GB each in 3a | 4 workers |
| reduced study (10 panels per combination) | 2.1 min | | gave the 42-min estimate |
| study (33 combinations x 200 panels) | 43.8 min | main 0.49 GB, 1.34 GB all processes | 4 workers, 10,459 worker-seconds |
| analyse | 3 s | 0.12 GB | |
| T_min, l_max, beta-hat | 1.7 min | | guard reads mids only |
| semi-synthetic smoke check | 6 s | 0.46 GB | |
| test suite | 44 s | | 297 passed, 28 skipped |

Total compute about 62 min, inside the 2.5 h target. Every job checkpoints under `$ITCH_DATA/phase3c/` and
resumes.
