# Phase 3b report: sign model, clipping, usable a range

Labels: **Measurement** = what the data or a simulation showed; **Choice** = a decision taken here.
Real data are read only through `explore.load_explore_runs` (exploration stocks, 8 study days), once, into
`$ITCH_DATA/phase3/p3b/raw_<cell>.npz` (per stock-day signs and sizes). Nothing combines signs or sizes with
price changes. Code: `tools/est/signs3b.py`, `clip3b.py`, `phase3b.py` (driver, per-cell checkpoints under
`$ITCH_DATA/phase3/p3b/`), tests `tools/est/test_phase3b.py`. Tables: `data/phase3/p3b_signs.tsv`,
`p3b_gamma_estimators.tsv`, `p3b_clip.tsv`, `p3b_r1_slope.tsv`, `p3b_ess.tsv`, `p3b_gt_clip.tsv`.

## 1. Sign model

### 1a. Uncentred and centred C(l), and the daily mean sign

Cu is the uncentred pooled C(l) (as in 3a, the same lags and per-stock-day N / 10 rule; the code asserts it
equals the calibration's). Cc is the same with every stock-day demeaned by its own mean sign. var(m_d) is
the variance across stock-days of the daily mean sign. The slope is the OLS of log var(m_d) on log
(median N) over 8 quantile bins of log N (bins with at least 10 stock-days).

| cell | stock-days | Cu(1) / Cc(1) | Cu(10) / Cc(10) | Cu(100) / Cc(100) | var(m_d) | slope log var on log N | var bins |
|---|---|---|---|---|---|---|---|
| N-T1S1 | 352 | 0.483 / 0.459 | 0.022 / 0.001 | 0.025 / 0.007 | 0.290 | -0.54 | 8 |
| N-T1S2 | 632 | 0.485 / 0.463 | 0.018 / -0.003 | 0.022 / 0.005 | 0.050 | -0.56 | 8 |
| N-T1S3 | 952 | 0.482 / 0.456 | 0.019 / -0.006 | 0.028 / 0.006 | 0.065 | -0.53 | 8 |
| N-T2S1 | 424 | 0.425 / 0.411 | 0.012 / -0.002 | 0.016 / 0.005 | 0.096 | -0.60 | 8 |
| N-T2S2 | 320 | 0.405 / 0.390 | 0.009 / -0.006 | 0.024 / 0.011 | 0.051 | -0.60 | 8 |
| N-T2S3 | 112 | 0.402 / 0.369 | 0.023 / -0.008 |  /  | 0.073 | -0.63 | 8 |
| N-T3S1 | 384 | 0.358 / 0.347 | 0.014 / 0.003 | 0.015 / 0.004 | 0.019 | -0.50 | 8 |
| N-T3S2 | 208 | 0.383 / 0.368 | 0.026 / 0.011 | 0.012 / 0.002 | 0.024 | -0.54 | 8 |
| N-T3S3 | 104 | 0.460 / 0.442 | 0.049 / 0.031 | 0.024 / 0.011 | 0.023 | -0.50 | 8 |
| Q-T1S1 | 552 | 0.586 / 0.559 | 0.088 / 0.064 | 0.031 / 0.014 | 0.151 | -0.58 | 8 |
| Q-T1S2 | 736 | 0.552 / 0.524 | 0.046 / 0.019 | 0.034 / 0.011 | 0.059 | -0.50 | 8 |
| Q-T1S3 | 864 | 0.537 / 0.489 | 0.052 / 0.006 | 0.038 / 0.005 | 0.128 | -0.50 | 8 |
| Q-T2S1 | 280 | 0.463 / 0.444 | 0.013 / -0.006 | 0.018 / 0.003 | 0.046 | -0.51 | 8 |
| Q-T2S2 | 168 | 0.467 / 0.447 | 0.014 / -0.005 | 0.021 / 0.005 | 0.038 | -0.48 | 8 |
| Q-T2S3 | 56 | 0.460 / 0.439 | 0.033 / 0.014 |  /  | 0.048 |  | 0 |
| Q-T3S1 | 176 | 0.406 / 0.392 | 0.031 / 0.016 | 0.022 / 0.008 | 0.020 | -0.32 | 8 |
| Q-T3S2 | 88 | 0.443 / 0.422 | 0.038 / 0.017 | 0.026 / 0.010 | 0.023 | -0.24 | 8 |
| Q-T3S3 | 80 | 0.517 / 0.478 | 0.089 / 0.050 | 0.034 / 0.010 | 0.059 | -0.32 | 8 |

**Measurement.**
- Demeaning removes most of the tail. Cc(100) is 0.002 to 0.014, against Cu(100) of 0.012 to 0.038. At lag
  10, Cc is negative in 7 cells (-0.002 to -0.008); that is the usual downward bias of a demeaned ACF, of
  order 1/N.
- The fast part is untouched: Cc(1) is within 0.05 of Cu(1).
- var(m_d) scales as N^-0.24 to N^-0.63. That is between iid signs (slope -1) and a fixed per-day
  component (slope 0), and slower than iid in every cell: the daily mean sign carries more spread than its
  N would give. The NASDAQ large-tick cells (Q-T3*) are flattest (-0.24 to -0.32).
- In Q-T2S3 no log-N bin has 10 stock-days (56 stock-days in all), so no slope is reported.

### 1b. Two-component generator eps = sign(X + m_d), X = fGn(H), m_d ~ N(0, tau^2)

Its uncentred C(l) is exactly (2/pi) arcsin((rho_H(l) + tau^2) / (1 + tau^2)). The fit is least squares on
C over the log-spaced lags present (12 per decade; lags 1 to 10 thin out to 1, 2, 3, 4, 5, 6, 7, 9, 10).
The H grid is 0.505..0.995 (step 0.005), widened to 0.01..0.999 (step 0.001) when a fit sits on an edge;
the tau grid is 0..0.6 (step 0.0025). Fit quality is chi2/df = mean of ((C_model - C) / SE)^2 over the
same lags, SE from the calibration's stock-day bootstrap. **Choice:** a model "fails" in a cell when
chi2/df > 4 (a typical residual above 2 SE). The one-component columns are the 3a model (fGn + AR(1),
4 parameters) on the same lags and metric. The tail fit is a **sensitivity** (Choice): the same
two-component fit on lags >= 10 only.

| cell | two-comp H | tau | edge (first/final) | RMSE C | chi2/df | fails | worst lag (resid) | one-comp RMSE C | chi2/df | fails | tail fit H | tail tau | tail edge (first/final) | tail RMSE |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N-T1S1 | 0.755 | 0.0000 | 0/0 | 0.0554 | 28.0 | yes | 1 (-0.204) | 0.0225 | 1.9 | no | 0.320 | 0.2050 | 1/0 | 0.0041 |
| N-T1S2 | 0.750 | 0.0000 | 0/0 | 0.0552 | 203.2 | yes | 1 (-0.213) | 0.0298 | 43.0 | yes | 0.223 | 0.1975 | 1/0 | 0.0039 |
| N-T1S3 | 0.750 | 0.0000 | 0/0 | 0.0567 | 254.1 | yes | 1 (-0.210) | 0.0349 | 72.7 | yes | 0.346 | 0.2175 | 1/0 | 0.0037 |
| N-T2S1 | 0.730 | 0.0000 | 0/0 | 0.0466 | 124.1 | yes | 1 (-0.180) | 0.0277 | 21.0 | yes | 0.510 | 0.1650 | 0/0 | 0.0043 |
| N-T2S2 | 0.720 | 0.0000 | 0/0 | 0.0447 | 91.8 | yes | 1 (-0.173) | 0.0279 | 17.6 | yes | 0.307 | 0.1775 | 1/0 | 0.0036 |
| N-T2S3 | 0.735 | 0.0000 | 0/0 | 0.0463 | 16.2 | yes | 1 (-0.151) | 0.0204 | 3.0 | no | 0.530 | 0.2325 | 0/0 | 0.0081 |
| N-T3S1 | 0.710 | 0.0000 | 0/0 | 0.0338 | 110.4 | yes | 1 (-0.139) | 0.0091 | 4.1 | yes | 0.545 | 0.1525 | 0/0 | 0.0026 |
| N-T3S2 | 0.720 | 0.0000 | 0/0 | 0.0379 | 40.6 | yes | 1 (-0.151) | 0.0106 | 1.7 | no | 0.620 | 0.1375 | 0/0 | 0.0025 |
| N-T3S3 | 0.770 | 0.0000 | 0/0 | 0.0430 | 15.3 | yes | 1 (-0.160) | 0.0125 | 0.8 | no | 0.680 | 0.1750 | 0/0 | 0.0038 |
| Q-T1S1 | 0.810 | 0.0000 | 0/0 | 0.0607 | 32.9 | yes | 1 (-0.225) | 0.0173 | 1.9 | no | 0.740 | 0.1425 | 0/0 | 0.0035 |
| Q-T1S2 | 0.790 | 0.0000 | 0/0 | 0.0595 | 219.7 | yes | 1 (-0.223) | 0.0243 | 28.0 | yes | 0.645 | 0.2050 | 0/0 | 0.0047 |
| Q-T1S3 | 0.790 | 0.0000 | 0/0 | 0.0583 | 136.6 | yes | 1 (-0.207) | 0.0187 | 12.1 | yes | 0.605 | 0.2625 | 0/0 | 0.0032 |
| Q-T2S1 | 0.750 | 0.0000 | 0/0 | 0.0517 | 61.7 | yes | 1 (-0.192) | 0.0291 | 11.4 | yes | 0.228 | 0.1850 | 1/0 | 0.0045 |
| Q-T2S2 | 0.745 | 0.0000 | 0/0 | 0.0576 | 66.3 | yes | 1 (-0.202) | 0.0317 | 12.3 | yes | 0.343 | 0.1950 | 1/0 | 0.0047 |
| Q-T2S3 | 0.755 | 0.0000 | 0/0 | 0.0564 | 26.3 | yes | 1 (-0.181) | 0.0213 | 3.4 | no | 0.615 | 0.1850 | 0/0 | 0.0076 |
| Q-T3S1 | 0.740 | 0.0000 | 0/0 | 0.0367 | 24.7 | yes | 1 (-0.148) | 0.0087 | 0.9 | no | 0.615 | 0.1725 | 0/0 | 0.0026 |
| Q-T3S2 | 0.755 | 0.0000 | 0/0 | 0.0446 | 20.9 | yes | 1 (-0.164) | 0.0135 | 1.5 | no | 0.635 | 0.1975 | 0/0 | 0.0043 |
| Q-T3S3 | 0.800 | 0.0000 | 0/0 | 0.0469 | 9.8 | yes | 1 (-0.172) | 0.0174 | 1.2 | no | 0.745 | 0.1700 | 0/0 | 0.0036 |

**Measurement.**
- The two-component model fails in all 18 cells (chi2/df 9.8 to 254). Least squares puts tau = 0 in every
  cell and H = 0.71 to 0.81. The worst residual is always at lag 1 (-0.14 to -0.23).
- Two parameters cannot fit both the steep fall from lag 1 to lag 10 and the plateau. Lags 1 to 5 carry
  most of the absolute error, so the plateau is given up.
- The 3a one-component model fails in 9 cells (chi2/df 4.1 to 73) and passes in 9 (0.8 to 3.4).
- No all-lags fit sat on an H edge, so the widening was not needed there.
- The tail fit gets the plateau right: tau = 0.14 to 0.26 (tau^2 = 0.02 to 0.07, the size of var(m_d) at
  large N), RMSE 0.0025 to 0.008.
- In 6 cells (N-T1S1, N-T1S2, N-T1S3, N-T2S2, Q-T2S1, Q-T2S2) the tail fit first hit the H = 0.505 edge.
  The widened grid then gave H = 0.22 to 0.35, an anti-persistent within-day remainder; none of these
  final fits is on an edge.

### 1c. gamma estimators on the fitted generators

50 panels per cell, each with the cell's actual stock-day N values; true exponent gamma_true = 2 - 2H.
- **LW:** local Whittle on stock-days with N >= 500, pooled by inverse variance.
- **power:** the calibration's power-law fit of the uncentred C on [10, l_hi] (bootstrap B = 199 for the
  l_hi cut; Choice).
- **v0:** 2 - zeta_V(0), with zeta_V(0) from the OLS slope of log v_0(T) on T = 8..64, on stock-days with
  N >= 256 (a fixed set at every T; Choice).

This is characterization only; no estimator is picked.

Generator: all-lags fit (the spec's generator). gamma_true = 2 - 2H; bias = mean over 50 panels minus gamma_true; sd over panels.

| cell | gamma_true | LW bias (sd) | power bias (sd) | v0 bias (sd) |
|---|---|---|---|---|
| N-T1S1 | 0.49 | 0.076 (0.011) | 0.022 (0.063) | 0.058 (0.014) |
| N-T1S2 | 0.50 | 0.082 (0.005) | -0.002 (0.037) | 0.055 (0.006) |
| N-T1S3 | 0.50 | 0.086 (0.006) | 0.009 (0.048) | 0.057 (0.005) |
| N-T2S1 | 0.54 | 0.075 (0.006) | -0.006 (0.045) | 0.058 (0.006) |
| N-T2S2 | 0.56 | 0.076 (0.006) | -0.008 (0.056) | 0.058 (0.009) |
| N-T2S3 | 0.53 | 0.089 (0.019) | -0.018 (0.142) | 0.057 (0.025) |
| N-T3S1 | 0.58 | 0.073 (0.005) | -0.022 (0.036) | 0.058 (0.005) |
| N-T3S2 | 0.56 | 0.077 (0.010) | 0.001 (0.061) | 0.056 (0.011) |
| N-T3S3 | 0.46 | 0.087 (0.013) | 0.010 (0.087) | 0.060 (0.015) |
| Q-T1S1 | 0.38 | 0.079 (0.007) | 0.002 (0.054) | 0.052 (0.012) |
| Q-T1S2 | 0.42 | 0.087 (0.005) | 0.018 (0.056) | 0.056 (0.007) |
| Q-T1S3 | 0.42 | 0.094 (0.007) | 0.011 (0.058) | 0.056 (0.007) |
| Q-T2S1 | 0.50 | 0.078 (0.009) | 0.000 (0.050) | 0.058 (0.009) |
| Q-T2S2 | 0.51 | 0.086 (0.011) | 0.015 (0.080) | 0.058 (0.014) |
| Q-T2S3 | 0.49 | 0.084 (0.021) | 0.008 (0.093) | 0.051 (0.023) |
| Q-T3S1 | 0.52 | 0.072 (0.008) | -0.009 (0.053) | 0.056 (0.009) |
| Q-T3S2 | 0.49 | 0.088 (0.013) | 0.025 (0.105) | 0.058 (0.016) |
| Q-T3S3 | 0.40 | 0.088 (0.015) | 0.016 (0.069) | 0.057 (0.020) |

Generator: tail-lags fit (sensitivity). gamma_true = 2 - 2H; bias = mean over 50 panels minus gamma_true; sd over panels.

| cell | gamma_true | LW bias (sd) | power bias (sd) | v0 bias (sd) |
|---|---|---|---|---|
| N-T1S1 | 1.36 | -0.231 (0.011) | -1.393 (0.035) | -0.806 (0.049) |
| N-T1S2 | 1.55 | -0.413 (0.005) | -1.566 (0.050) | -1.017 (0.023) |
| N-T1S3 | 1.31 | -0.173 (0.005) | -1.334 (0.046) | -0.771 (0.016) |
| N-T2S1 | 0.98 | 0.013 (0.007) | -0.957 (0.041) | -0.274 (0.020) |
| N-T2S2 | 1.39 | -0.244 (0.007) | -1.429 (0.051) | -0.762 (0.034) |
| N-T2S3 | 0.94 | 0.027 (0.021) | -0.896 (0.194) | -0.367 (0.048) |
| N-T3S1 | 0.91 | 0.033 (0.005) | -0.862 (0.038) | -0.188 (0.016) |
| N-T3S2 | 0.76 | 0.060 (0.010) | -0.521 (0.079) | -0.070 (0.017) |
| N-T3S3 | 0.64 | 0.077 (0.013) | -0.375 (0.083) | -0.060 (0.026) |
| Q-T1S1 | 0.52 | 0.072 (0.008) | -0.168 (0.061) | 0.009 (0.015) |
| Q-T1S2 | 0.71 | 0.068 (0.005) | -0.542 (0.068) | -0.129 (0.016) |
| Q-T1S3 | 0.79 | 0.062 (0.006) | -0.728 (0.045) | -0.272 (0.014) |
| Q-T2S1 | 1.54 | -0.412 (0.007) | -1.558 (0.047) | -0.972 (0.037) |
| Q-T2S2 | 1.31 | -0.175 (0.012) | -1.359 (0.062) | -0.714 (0.044) |
| Q-T2S3 | 0.77 | 0.061 (0.023) | -0.557 (0.292) | -0.139 (0.041) |
| Q-T3S1 | 0.77 | 0.055 (0.008) | -0.623 (0.056) | -0.121 (0.023) |
| Q-T3S2 | 0.73 | 0.070 (0.012) | -0.536 (0.090) | -0.137 (0.033) |
| Q-T3S3 | 0.51 | 0.081 (0.014) | -0.150 (0.077) | -0.004 (0.024) |

The same estimators on the exploration signs (for reference; power uses B = 199 for its l_hi cut):

| cell | LW | power | v0 |
|---|---|---|---|
| N-T1S1 | 0.903 | 0.000 | 0.745 |
| N-T1S2 | 0.906 | 0.037 | 0.748 |
| N-T1S3 | 0.915 | -0.019 | 0.733 |
| N-T2S1 | 0.926 | 0.090 | 0.785 |
| N-T2S2 | 0.937 | -0.011 | 0.777 |
| N-T2S3 | 0.897 | 0.047 | 0.687 |
| N-T3S1 | 0.906 | 0.093 | 0.781 |
| N-T3S2 | 0.884 | 0.245 | 0.756 |
| N-T3S3 | 0.767 | 0.318 | 0.656 |
| Q-T1S1 | 0.773 | 0.308 | 0.620 |
| Q-T1S2 | 0.802 | 0.194 | 0.684 |
| Q-T1S3 | 0.805 | 0.077 | 0.639 |
| Q-T2S1 | 0.904 | 0.063 | 0.763 |
| Q-T2S2 | 0.892 | -0.080 | 0.763 |
| Q-T2S3 | 0.803 | -0.172 | 0.718 |
| Q-T3S1 | 0.852 | 0.172 | 0.720 |
| Q-T3S2 | 0.791 | 0.201 | 0.662 |
| Q-T3S3 | 0.698 | 0.351 | 0.576 |

**Measurement.**
- **All-lags generator** (pure fGn, tau = 0):
  - LW is biased up by +0.07 to +0.09 (sd 0.005 to 0.02).
  - v0 is biased up by +0.05 to +0.06 (sd 0.005 to 0.025).
  - The power fit is unbiased within its spread (-0.02 to +0.03), but that spread is several times larger
    (sd 0.04 to 0.14).
- **Tail generator** (a per-day component tau of 0.14 to 0.26):
  - The power fit is biased down by -0.15 to -1.57: the plateau flattens C(l) and pulls the fitted exponent
    toward 0, as on the real signs.
  - v0 is biased by +0.01 to -1.02.
  - LW is the least affected (-0.41 to +0.08), because it ignores frequency 0. Its bias is largest in the
    anti-persistent cells (gamma_true 1.3 to 1.55).
- On the real signs the three give 0.70 to 0.94 (LW), -0.17 to 0.35 (power) and 0.58 to 0.79 (v0). That
  ordering is the tail generator's ordering, not the all-lags generator's.

## 2. Clipping (sizes only)

R1 = 1 % of the stock-day's own total run volume (`moments.clip_volume`); R2 / R3 = the stock's 99.9th /
99th percentile of run size over its 8 days (numpy linear quantile; Choice). Shares are pooled over all
runs of the cell. At a = 0, q^0 = 1 for clipped and unclipped runs alike, so the removed share is 0.

| cell | rule | runs clipped | sum q removed | q^0.5 removed | q^1 removed | q^1.5 removed | q^2 removed | q^2.5 removed | q^3 removed | q^3.5 removed | q^4 removed |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N-T1S1 | R1 | 0.0235 | 0.0320 | 0.014 | 0.032 | 0.078 | 0.193 | 0.357 | 0.493 | 0.584 | 0.647 |
| N-T1S1 | R2 | 0.0011 | 0.0094 | 0.001 | 0.009 | 0.051 | 0.179 | 0.373 | 0.540 | 0.652 | 0.727 |
| N-T1S1 | R3 | 0.0086 | 0.0442 | 0.007 | 0.044 | 0.172 | 0.449 | 0.748 | 0.911 | 0.971 | 0.991 |
| N-T1S2 | R1 | 0.0115 | 0.0149 | 0.005 | 0.015 | 0.043 | 0.153 | 0.475 | 0.820 | 0.952 | 0.986 |
| N-T1S2 | R2 | 0.0010 | 0.0059 | 0.001 | 0.006 | 0.033 | 0.153 | 0.498 | 0.847 | 0.967 | 0.993 |
| N-T1S2 | R3 | 0.0091 | 0.0246 | 0.005 | 0.025 | 0.085 | 0.260 | 0.621 | 0.908 | 0.986 | 0.998 |
| N-T1S3 | R1 | 0.0176 | 0.0213 | 0.008 | 0.021 | 0.048 | 0.127 | 0.327 | 0.594 | 0.768 | 0.854 |
| N-T1S3 | R2 | 0.0010 | 0.0054 | 0.001 | 0.005 | 0.028 | 0.116 | 0.356 | 0.666 | 0.852 | 0.928 |
| N-T1S3 | R3 | 0.0091 | 0.0228 | 0.005 | 0.023 | 0.076 | 0.221 | 0.514 | 0.812 | 0.949 | 0.988 |
| N-T2S1 | R1 | 0.0066 | 0.0132 | 0.004 | 0.013 | 0.045 | 0.185 | 0.548 | 0.846 | 0.951 | 0.983 |
| N-T2S1 | R2 | 0.0010 | 0.0080 | 0.001 | 0.008 | 0.048 | 0.229 | 0.635 | 0.910 | 0.982 | 0.996 |
| N-T2S1 | R3 | 0.0090 | 0.0290 | 0.005 | 0.029 | 0.110 | 0.347 | 0.740 | 0.949 | 0.992 | 0.999 |
| N-T2S2 | R1 | 0.0065 | 0.0121 | 0.003 | 0.012 | 0.046 | 0.184 | 0.492 | 0.751 | 0.874 | 0.930 |
| N-T2S2 | R2 | 0.0010 | 0.0092 | 0.001 | 0.009 | 0.055 | 0.244 | 0.617 | 0.874 | 0.961 | 0.986 |
| N-T2S2 | R3 | 0.0087 | 0.0313 | 0.006 | 0.031 | 0.121 | 0.372 | 0.744 | 0.942 | 0.990 | 0.998 |
| N-T2S3 | R1 | 0.0384 | 0.0968 | 0.024 | 0.097 | 0.323 | 0.678 | 0.889 | 0.962 | 0.986 | 0.995 |
| N-T2S3 | R2 | 0.0012 | 0.0092 | 0.001 | 0.009 | 0.058 | 0.185 | 0.326 | 0.439 | 0.529 | 0.603 |
| N-T2S3 | R3 | 0.0096 | 0.0650 | 0.009 | 0.065 | 0.285 | 0.657 | 0.884 | 0.963 | 0.988 | 0.996 |
| N-T3S1 | R1 | 0.0026 | 0.0103 | 0.002 | 0.010 | 0.075 | 0.485 | 0.911 | 0.988 | 0.998 | 1.000 |
| N-T3S1 | R2 | 0.0010 | 0.0103 | 0.001 | 0.010 | 0.087 | 0.520 | 0.929 | 0.993 | 0.999 | 1.000 |
| N-T3S1 | R3 | 0.0096 | 0.0360 | 0.006 | 0.036 | 0.162 | 0.610 | 0.953 | 0.996 | 1.000 | 1.000 |
| N-T3S2 | R1 | 0.0094 | 0.0180 | 0.005 | 0.018 | 0.050 | 0.114 | 0.224 | 0.380 | 0.553 | 0.702 |
| N-T3S2 | R2 | 0.0011 | 0.0076 | 0.001 | 0.008 | 0.035 | 0.109 | 0.256 | 0.463 | 0.672 | 0.825 |
| N-T3S2 | R3 | 0.0097 | 0.0348 | 0.006 | 0.035 | 0.113 | 0.260 | 0.465 | 0.677 | 0.839 | 0.932 |
| N-T3S3 | R1 | 0.0086 | 0.0425 | 0.007 | 0.043 | 0.154 | 0.360 | 0.600 | 0.789 | 0.902 | 0.958 |
| N-T3S3 | R2 | 0.0010 | 0.0076 | 0.001 | 0.008 | 0.032 | 0.087 | 0.180 | 0.307 | 0.448 | 0.584 |
| N-T3S3 | R3 | 0.0099 | 0.0474 | 0.008 | 0.047 | 0.145 | 0.305 | 0.493 | 0.667 | 0.801 | 0.890 |
| Q-T1S1 | R1 | 0.0221 | 0.0283 | 0.011 | 0.028 | 0.065 | 0.135 | 0.214 | 0.278 | 0.323 | 0.352 |
| Q-T1S1 | R2 | 0.0011 | 0.0183 | 0.002 | 0.018 | 0.117 | 0.397 | 0.706 | 0.881 | 0.956 | 0.985 |
| Q-T1S1 | R3 | 0.0097 | 0.0685 | 0.010 | 0.068 | 0.275 | 0.646 | 0.899 | 0.979 | 0.996 | 0.999 |
| Q-T1S2 | R1 | 0.0111 | 0.0204 | 0.005 | 0.020 | 0.076 | 0.269 | 0.588 | 0.807 | 0.909 | 0.955 |
| Q-T1S2 | R2 | 0.0010 | 0.0125 | 0.001 | 0.012 | 0.079 | 0.329 | 0.712 | 0.916 | 0.978 | 0.994 |
| Q-T1S2 | R3 | 0.0095 | 0.0450 | 0.008 | 0.045 | 0.178 | 0.506 | 0.853 | 0.975 | 0.997 | 1.000 |
| Q-T1S3 | R1 | 0.0298 | 0.0523 | 0.017 | 0.052 | 0.151 | 0.409 | 0.739 | 0.910 | 0.968 | 0.987 |
| Q-T1S3 | R2 | 0.0011 | 0.0107 | 0.001 | 0.011 | 0.064 | 0.264 | 0.583 | 0.792 | 0.886 | 0.928 |
| Q-T1S3 | R3 | 0.0094 | 0.0410 | 0.007 | 0.041 | 0.161 | 0.463 | 0.805 | 0.953 | 0.990 | 0.997 |
| Q-T2S1 | R1 | 0.0086 | 0.0202 | 0.005 | 0.020 | 0.076 | 0.259 | 0.560 | 0.767 | 0.868 | 0.919 |
| Q-T2S1 | R2 | 0.0010 | 0.0121 | 0.001 | 0.012 | 0.076 | 0.316 | 0.690 | 0.899 | 0.966 | 0.987 |
| Q-T2S1 | R3 | 0.0095 | 0.0414 | 0.007 | 0.041 | 0.164 | 0.474 | 0.823 | 0.964 | 0.993 | 0.999 |
| Q-T2S2 | R1 | 0.0102 | 0.0276 | 0.007 | 0.028 | 0.090 | 0.245 | 0.484 | 0.690 | 0.811 | 0.876 |
| Q-T2S2 | R2 | 0.0011 | 0.0097 | 0.001 | 0.010 | 0.050 | 0.178 | 0.402 | 0.610 | 0.737 | 0.807 |
| Q-T2S2 | R3 | 0.0091 | 0.0363 | 0.006 | 0.036 | 0.132 | 0.347 | 0.624 | 0.814 | 0.898 | 0.932 |
| Q-T2S3 | R1 | 0.0127 | 0.0477 | 0.010 | 0.048 | 0.188 | 0.511 | 0.798 | 0.921 | 0.966 | 0.984 |
| Q-T2S3 | R2 | 0.0010 | 0.0209 | 0.002 | 0.021 | 0.137 | 0.475 | 0.804 | 0.940 | 0.982 | 0.994 |
| Q-T2S3 | R3 | 0.0093 | 0.0572 | 0.009 | 0.057 | 0.240 | 0.620 | 0.896 | 0.979 | 0.996 | 0.999 |
| Q-T3S1 | R1 | 0.0025 | 0.0113 | 0.002 | 0.011 | 0.064 | 0.255 | 0.544 | 0.736 | 0.832 | 0.884 |
| Q-T3S1 | R2 | 0.0010 | 0.0139 | 0.001 | 0.014 | 0.089 | 0.349 | 0.700 | 0.893 | 0.963 | 0.986 |
| Q-T3S1 | R3 | 0.0099 | 0.0471 | 0.007 | 0.047 | 0.192 | 0.523 | 0.843 | 0.965 | 0.993 | 0.999 |
| Q-T3S2 | R1 | 0.0094 | 0.0312 | 0.005 | 0.031 | 0.142 | 0.464 | 0.810 | 0.950 | 0.987 | 0.996 |
| Q-T3S2 | R2 | 0.0011 | 0.0178 | 0.002 | 0.018 | 0.119 | 0.447 | 0.808 | 0.952 | 0.989 | 0.997 |
| Q-T3S2 | R3 | 0.0098 | 0.0528 | 0.008 | 0.053 | 0.214 | 0.573 | 0.880 | 0.977 | 0.996 | 0.999 |
| Q-T3S3 | R1 | 0.0122 | 0.0943 | 0.016 | 0.094 | 0.248 | 0.414 | 0.559 | 0.677 | 0.768 | 0.835 |
| Q-T3S3 | R2 | 0.0010 | 0.0072 | 0.001 | 0.007 | 0.025 | 0.054 | 0.095 | 0.146 | 0.206 | 0.273 |
| Q-T3S3 | R3 | 0.0100 | 0.0433 | 0.007 | 0.043 | 0.127 | 0.241 | 0.369 | 0.495 | 0.610 | 0.709 |

R1: share of runs clipped per stock-day regressed on log N.

| cell | slope [95% stock-cluster CI] | mean share clipped, N <= 100 | N > 1000 | stock-days with N <= 100 having a clipped run |
|---|---|---|---|---|
| N-T1S1 | -0.162 [-0.176, -0.151] | 0.870 | 0.001 | 1.000 |
| N-T1S2 | -0.142 [-0.157, -0.117] | 0.610 | 0.000 | 1.000 |
| N-T1S3 | -0.154 [-0.163, -0.143] | 0.626 | 0.000 | 1.000 |
| N-T2S1 | -0.164 [-0.178, -0.150] | 0.821 | 0.000 | 1.000 |
| N-T2S2 | -0.137 [-0.155, -0.104] | 0.699 | 0.000 | 1.000 |
| N-T2S3 | -0.192 [-0.215, -0.153] | 0.634 | 0.001 | 1.000 |
| N-T3S1 | -0.106 [-0.141, -0.041] | 0.723 | 0.000 | 1.000 |
| N-T3S2 | -0.085 [-0.131, -0.028] | 0.554 | 0.001 | 1.000 |
| N-T3S3 | -0.058 [-0.085, -0.021] | 0.316 | 0.001 | 1.000 |
| Q-T1S1 | -0.171 [-0.202, -0.146] | 0.735 | 0.000 | 1.000 |
| Q-T1S2 | -0.132 [-0.150, -0.109] | 0.615 | 0.001 | 1.000 |
| Q-T1S3 | -0.167 [-0.176, -0.159] | 0.658 | 0.001 | 1.000 |
| Q-T2S1 | -0.123 [-0.157, -0.072] | 0.655 | 0.001 | 1.000 |
| Q-T2S2 | -0.153 [-0.175, -0.082] | 0.598 | 0.001 | 1.000 |
| Q-T2S3 | -0.164 [-0.192, -0.014] | 0.557 | 0.001 | 1.000 |
| Q-T3S1 | -0.038 [-0.069, -0.008] | 0.630 | 0.001 | 1.000 |
| Q-T3S2 | -0.034 [-0.056, -0.020] |  | 0.002 |  |
| Q-T3S3 | -0.071 [-0.102, -0.026] | 0.374 | 0.001 | 1.000 |

**Measurement.**
- R1 clips 0.25 % to 3.8 % of runs and 1.0 % to 9.7 % of sum q.
- R2 clips 0.10 % to 0.12 % of runs (the quantile rule's own 0.1 %, plus ties) and 0.5 % to 2.1 % of sum q.
- R3 clips 0.86 % to 1.0 % of runs and 2.3 % to 6.9 % of sum q.
- At a = 1 (q^2) the rules remove 5 % to 68 % of the moment, and at a = 2 (q^4) 27 % to 100 %. In most
  cells a quantile cap on the top 0.1 % of runs already removes most of q^4 (R2 above 80 % in 14 of 18).
- **R1 makes the clip depend on N.** The share of runs clipped falls with log N in all 18 cells, slope
  -0.03 to -0.19 per unit log N, and every 95 % stock-cluster interval excludes 0.
- Every stock-day with N <= 100 has at least one clipped run. This is forced: the largest run is at least
  V / N >= V / 100. Such stock-days lose 32 % to 87 % of their runs to the cap on average, against about
  0.1 % above N = 1000.

## 3. Usable a range (Kish ESS)

ESS = (sum w)^2 / sum w^2 per stock-day, with w = clipped q^(2a) for V moments and q^a for K moments. At
a = 0, ESS = N. a_max(threshold) is the largest grid a whose median stock-day ESS meets the threshold; "-"
means even a = 0 fails (median N below the threshold). The pooled a_max uses the ESS of each stock's 8
days together, then the median over stocks. No threshold is chosen.

Weights q^(2a) (V moments): median stock-day ESS at a = 0, 0.5, 1, 1.5, 2 (10th percentile at a = 1), a_max at 10 / 30 / 100, pooled-over-8-days a_max at 10 / 30 / 100.

| cell | rule | median ESS | p10 at a=1 | a_max 10/30/100 | pooled a_max 10/30/100 |
|---|---|---|---|---|---|
| N-T1S1 | none | 18 / 10 / 5 / 3 / 2 | 1.3 | 0.25/-/- | 1/0.75/0.25 |
| N-T1S1 | R1 | 18 / 18 / 17 / 13 / 10 | 3.0 | 1.75/-/- | 2/1.75/0.5 |
| N-T1S1 | R2 | 18 / 10 / 5 / 4 / 3 | 1.3 | 0.5/-/- | 1.25/0.75/0.25 |
| N-T1S1 | R3 | 18 / 11 / 6 / 4 / 3 | 1.4 | 0.5/-/- | 2/0.75/0.25 |
| N-T1S2 | none | 542 / 245 / 59 / 16 / 7 | 6.8 | 1.5/1/0.75 | 1.5/1.25/1 |
| N-T1S2 | R1 | 542 / 267 / 103 / 49 / 27 | 35.4 | 2/1.75/1 | 2/2/1.25 |
| N-T1S2 | R2 | 542 / 253 / 76 / 24 / 10 | 10.3 | 2/1.25/0.75 | 2/2/1.5 |
| N-T1S2 | R3 | 542 / 266 / 115 / 51 / 26 | 13.3 | 2/1.75/1 | 2/2/2 |
| N-T1S3 | none | 432 / 185 / 63 / 20 / 8 | 5.3 | 1.75/1.25/0.75 | 1.75/1.25/1 |
| N-T1S3 | R1 | 432 / 205 / 104 / 53 / 30 | 21.9 | 2/2/1 | 2/2/1.5 |
| N-T1S3 | R2 | 432 / 191 / 72 / 25 / 10 | 6.8 | 2/1.25/0.75 | 2/2/1.5 |
| N-T1S3 | R3 | 432 / 200 / 96 / 46 / 25 | 8.1 | 2/1.75/0.75 | 2/2/2 |
| N-T2S1 | none | 557 / 252 / 44 / 9 / 4 | 2.4 | 1.25/1/0.75 | 1.25/1/0.75 |
| N-T2S1 | R1 | 557 / 281 / 90 / 31 / 17 | 8.4 | 2/1.5/0.75 | 2/2/1.25 |
| N-T2S1 | R2 | 557 / 268 / 62 / 17 / 7 | 2.7 | 1.75/1/0.75 | 2/2/1.5 |
| N-T2S1 | R3 | 557 / 278 / 100 / 37 / 18 | 2.9 | 2/1.5/0.75 | 2/2/2 |
| N-T2S2 | none | 738 / 346 / 70 / 14 / 6 | 3.3 | 1.5/1/0.75 | 1.5/1.25/1 |
| N-T2S2 | R1 | 738 / 360 / 100 / 42 / 19 | 24.8 | 2/1.5/1 | 2/1.75/1.25 |
| N-T2S2 | R2 | 738 / 355 / 101 / 25 / 10 | 4.8 | 2/1.25/1 | 2/2/1.5 |
| N-T2S2 | R3 | 738 / 368 / 151 / 57 / 29 | 9.1 | 2/1.75/1 | 2/2/2 |
| N-T2S3 | none | 306 / 124 / 29 / 10 / 5 | 3.4 | 1.25/0.75/0.5 | 1.25/0.75/0.75 |
| N-T2S3 | R1 | 306 / 154 / 86 / 45 / 30 | 26.8 | 2/1.75/0.75 | 2/2/1.5 |
| N-T2S3 | R2 | 306 / 132 / 35 / 11 / 5 | 4.4 | 1.5/1/0.5 | 2/1.25/1 |
| N-T2S3 | R3 | 306 / 136 / 57 / 22 / 13 | 5.9 | 2/1.25/0.5 | 2/2/1.5 |
| N-T3S1 | none | 1557 / 617 / 83 / 16 / 6 | 8.2 | 1.5/1.25/0.75 | 1.75/1.25/1 |
| N-T3S1 | R1 | 1557 / 632 / 96 / 35 / 16 | 46.8 | 2/1.5/0.75 | 2/1.75/1.25 |
| N-T3S1 | R2 | 1557 / 658 / 129 / 34 / 15 | 12.2 | 2/1.5/1 | 2/2/2 |
| N-T3S1 | R3 | 1557 / 711 / 218 / 92 / 55 | 20.0 | 2/2/1.25 | 2/2/2 |
| N-T3S2 | none | 601 / 246 / 49 / 14 / 6 | 9.3 | 1.5/1/0.75 | 1.5/1.25/1 |
| N-T3S2 | R1 | 601 / 275 / 86 / 43 / 26 | 63.4 | 2/1.75/0.75 | 2/2/1.5 |
| N-T3S2 | R2 | 601 / 258 / 61 / 19 / 9 | 13.7 | 1.75/1.25/0.75 | 2/2/1.5 |
| N-T3S2 | R3 | 601 / 276 / 91 / 41 / 24 | 24.9 | 2/1.75/0.75 | 2/2/2 |
| N-T3S3 | none | 857 / 240 / 43 / 13 / 7 | 13.1 | 1.5/1/0.5 | 1.75/1.5/1 |
| N-T3S3 | R1 | 857 / 262 / 80 / 43 / 28 | 62.9 | 2/1.75/0.75 | 2/2/1.5 |
| N-T3S3 | R2 | 857 / 246 / 57 / 21 / 10 | 15.6 | 2/1.25/0.75 | 2/2/1.75 |
| N-T3S3 | R3 | 857 / 270 / 86 / 42 / 27 | 24.7 | 2/1.75/0.75 | 2/2/2 |
| Q-T1S1 | none | 163 / 71 / 18 / 7 / 4 | 2.5 | 1.25/0.75/0.25 | 1.25/1/0.75 |
| Q-T1S1 | R1 | 163 / 106 / 72 / 43 / 30 | 10.8 | 2/2/0.5 | 2/2/1.5 |
| Q-T1S1 | R2 | 163 / 74 / 21 / 9 / 5 | 3.0 | 1.25/0.75/0.25 | 2/1.25/1 |
| Q-T1S1 | R3 | 163 / 83 / 32 / 15 / 9 | 3.5 | 1.75/1/0.25 | 2/2/1.5 |
| Q-T1S2 | none | 580 / 222 / 39 / 9 / 4 | 4.8 | 1.25/1/0.5 | 1.25/1/0.75 |
| Q-T1S2 | R1 | 580 / 269 / 95 / 42 / 23 | 38.5 | 2/1.75/0.75 | 2/2/1.5 |
| Q-T1S2 | R2 | 580 / 244 / 53 / 15 / 7 | 8.1 | 1.75/1/0.75 | 2/1.75/1.25 |
| Q-T1S2 | R3 | 580 / 270 / 110 / 46 / 25 | 16.0 | 2/1.75/1 | 2/2/2 |
| Q-T1S3 | none | 297 / 103 / 24 / 8 / 4 | 2.5 | 1.25/0.75/0.5 | 1.25/1/0.75 |
| Q-T1S3 | R1 | 297 / 143 / 86 / 46 / 27 | 12.8 | 2/1.75/0.75 | 2/2/1.5 |
| Q-T1S3 | R2 | 297 / 116 / 32 / 11 / 5 | 3.0 | 1.5/1/0.5 | 2/1.5/1 |
| Q-T1S3 | R3 | 297 / 126 / 55 / 23 / 13 | 3.8 | 2/1.25/0.5 | 2/2/2 |
| Q-T2S1 | none | 472 / 196 / 36 / 8 / 4 | 3.5 | 1.25/1/0.5 | 1.25/1/0.75 |
| Q-T2S1 | R1 | 472 / 241 / 88 / 43 / 25 | 29.1 | 2/1.75/0.75 | 2/2/1.25 |
| Q-T2S1 | R2 | 472 / 215 / 48 / 15 / 7 | 5.1 | 1.5/1/0.75 | 2/2/1.25 |
| Q-T2S1 | R3 | 472 / 244 / 85 / 35 / 18 | 8.3 | 2/1.5/0.75 | 2/2/2 |
| Q-T2S2 | none | 714 / 270 / 36 / 8 / 4 | 3.2 | 1.25/1/0.5 | 1.25/1/0.75 |
| Q-T2S2 | R1 | 714 / 335 / 90 / 38 / 22 | 42.8 | 2/1.5/0.75 | 2/2/1.5 |
| Q-T2S2 | R2 | 714 / 295 / 52 / 15 / 6 | 3.9 | 1.5/1/0.75 | 2/2/1.25 |
| Q-T2S2 | R3 | 714 / 340 / 97 / 44 / 23 | 6.1 | 2/1.75/0.75 | 2/2/2 |
| Q-T2S3 | none | 559 / 222 / 44 / 8 / 4 | 4.3 | 1.25/1/0.5 | 1/0.75/0.75 |
| Q-T2S3 | R1 | 559 / 253 / 89 / 41 / 25 | 35.6 | 2/1.75/0.75 | 2/2/1.5 |
| Q-T2S3 | R2 | 559 / 229 / 46 / 13 / 5 | 7.0 | 1.5/1/0.75 | 2/1.5/1 |
| Q-T2S3 | R3 | 559 / 248 / 87 / 39 / 21 | 7.2 | 2/1.5/0.75 | 2/2/2 |
| Q-T3S1 | none | 1154 / 375 / 45 / 10 / 4 | 4.7 | 1.5/1/0.75 | 1.25/1/0.75 |
| Q-T3S1 | R1 | 1154 / 473 / 84 / 30 / 15 | 52.0 | 2/1.25/0.75 | 2/1.5/1.25 |
| Q-T3S1 | R2 | 1154 / 433 / 77 / 20 / 9 | 18.4 | 1.75/1.25/0.75 | 2/2/1.5 |
| Q-T3S1 | R3 | 1154 / 493 / 139 / 62 / 36 | 45.6 | 2/2/1 | 2/2/2 |
| Q-T3S2 | none | 642 / 204 / 31 / 9 / 4 | 4.8 | 1.25/1/0.5 | 1.5/1/0.75 |
| Q-T3S2 | R1 | 642 / 258 / 81 / 40 / 25 | 62.9 | 2/1.75/0.75 | 2/2/1.75 |
| Q-T3S2 | R2 | 642 / 218 / 47 / 13 / 7 | 14.8 | 1.5/1/0.5 | 2/2/1.25 |
| Q-T3S2 | R3 | 642 / 256 / 86 / 40 / 23 | 33.0 | 2/1.5/0.75 | 2/2/2 |
| Q-T3S3 | none | 614 / 145 / 24 / 8 / 5 | 7.0 | 1.25/0.75/0.5 | 1.75/1.25/0.75 |
| Q-T3S3 | R1 | 614 / 196 / 75 / 43 / 31 | 41.4 | 2/2/0.75 | 2/2/1.5 |
| Q-T3S3 | R2 | 614 / 148 / 31 / 12 / 6 | 7.6 | 1.5/1/0.5 | 2/2/1.25 |
| Q-T3S3 | R3 | 614 / 172 / 48 / 24 / 16 | 13.5 | 2/1.25/0.5 | 2/2/2 |

Weights q^a (K moments): median stock-day ESS at a = 0, 0.5, 1, 1.5, 2 (10th percentile at a = 1), a_max at 10 / 30 / 100, pooled-over-8-days a_max at 10 / 30 / 100.

| cell | rule | median ESS | p10 at a=1 | a_max 10/30/100 | pooled a_max 10/30/100 |
|---|---|---|---|---|---|
| N-T1S1 | none | 18 / 15 / 10 / 7 / 5 | 2.0 | 0.75/-/- | 2/1.5/0.5 |
| N-T1S1 | R1 | 18 / 18 / 18 / 17 / 17 | 3.0 | 2/-/- | 2/2/1 |
| N-T1S1 | R2 | 18 / 15 / 10 / 7 / 5 | 2.0 | 1/-/- | 2/1.5/0.5 |
| N-T1S1 | R3 | 18 / 15 / 11 / 8 / 6 | 2.0 | 1/-/- | 2/1.75/0.5 |
| N-T1S2 | none | 542 / 396 / 245 / 128 / 59 | 26.8 | 2/2/1.5 | 2/2/2 |
| N-T1S2 | R1 | 542 / 401 / 267 / 172 / 103 | 42.5 | 2/2/2 | 2/2/2 |
| N-T1S2 | R2 | 542 / 396 / 253 / 144 / 76 | 27.2 | 2/2/1.75 | 2/2/2 |
| N-T1S2 | R3 | 542 / 402 / 266 / 179 / 115 | 27.5 | 2/2/2 | 2/2/2 |
| N-T1S3 | none | 432 / 303 / 185 / 108 / 63 | 14.4 | 2/2/1.5 | 2/2/2 |
| N-T1S3 | R1 | 432 / 307 / 205 / 146 / 104 | 27.1 | 2/2/2 | 2/2/2 |
| N-T1S3 | R2 | 432 / 303 / 191 / 122 / 72 | 15.3 | 2/2/1.5 | 2/2/2 |
| N-T1S3 | R3 | 432 / 305 / 200 / 140 / 96 | 16.5 | 2/2/1.75 | 2/2/2 |
| N-T2S1 | none | 557 / 417 / 252 / 111 / 44 | 5.8 | 2/2/1.5 | 2/2/1.75 |
| N-T2S1 | R1 | 557 / 418 / 281 / 166 / 90 | 9.5 | 2/2/1.75 | 2/2/2 |
| N-T2S1 | R2 | 557 / 417 / 268 / 134 / 62 | 5.8 | 2/2/1.5 | 2/2/2 |
| N-T2S1 | R3 | 557 / 418 / 278 / 173 / 100 | 5.8 | 2/2/1.75 | 2/2/2 |
| N-T2S2 | none | 738 / 550 / 346 / 167 / 70 | 14.5 | 2/2/1.75 | 2/2/2 |
| N-T2S2 | R1 | 738 / 554 / 360 / 206 / 100 | 35.8 | 2/2/2 | 2/2/2 |
| N-T2S2 | R2 | 738 / 550 / 355 / 202 / 101 | 14.5 | 2/2/2 | 2/2/2 |
| N-T2S2 | R3 | 738 / 555 / 368 / 244 / 151 | 24.2 | 2/2/2 | 2/2/2 |
| N-T2S3 | none | 306 / 206 / 124 / 61 / 29 | 13.9 | 2/1.75/1 | 2/1.75/1.5 |
| N-T2S3 | R1 | 306 / 220 / 154 / 115 / 86 | 30.5 | 2/2/1.75 | 2/2/2 |
| N-T2S3 | R2 | 306 / 206 / 132 / 71 / 35 | 14.7 | 2/2/1 | 2/2/2 |
| N-T2S3 | R3 | 306 / 206 / 136 / 85 / 57 | 16.8 | 2/2/1.25 | 2/2/2 |
| N-T3S1 | none | 1557 / 1198 / 617 / 229 / 83 | 55.4 | 2/2/1.75 | 2/2/2 |
| N-T3S1 | R1 | 1557 / 1200 / 632 / 253 / 96 | 111.4 | 2/2/1.75 | 2/2/2 |
| N-T3S1 | R2 | 1557 / 1202 / 658 / 288 / 129 | 66.6 | 2/2/2 | 2/2/2 |
| N-T3S1 | R3 | 1557 / 1214 / 711 / 386 / 218 | 69.6 | 2/2/2 | 2/2/2 |
| N-T3S2 | none | 601 / 452 / 246 / 112 / 49 | 72.7 | 2/2/1.5 | 2/2/2 |
| N-T3S2 | R1 | 601 / 460 / 275 / 150 / 86 | 119.8 | 2/2/1.75 | 2/2/2 |
| N-T3S2 | R2 | 601 / 452 / 258 / 121 / 61 | 75.4 | 2/2/1.5 | 2/2/2 |
| N-T3S2 | R3 | 601 / 460 / 276 / 153 / 91 | 89.0 | 2/2/1.75 | 2/2/2 |
| N-T3S3 | none | 857 / 540 / 240 / 96 / 43 | 68.1 | 2/2/1.25 | 2/2/2 |
| N-T3S3 | R1 | 857 / 552 / 262 / 132 / 80 | 111.6 | 2/2/1.5 | 2/2/2 |
| N-T3S3 | R2 | 857 / 547 / 246 / 110 / 57 | 69.2 | 2/2/1.5 | 2/2/2 |
| N-T3S3 | R3 | 857 / 550 / 270 / 142 / 86 | 82.7 | 2/2/1.75 | 2/2/2 |
| Q-T1S1 | none | 163 / 119 / 71 / 34 / 18 | 6.7 | 2/1.5/0.5 | 2/2/1.5 |
| Q-T1S1 | R1 | 163 / 132 / 106 / 92 / 72 | 11.2 | 2/2/1 | 2/2/2 |
| Q-T1S1 | R2 | 163 / 122 / 74 / 39 / 21 | 6.7 | 2/1.5/0.5 | 2/2/2 |
| Q-T1S1 | R3 | 163 / 123 / 83 / 51 / 32 | 7.2 | 2/2/0.75 | 2/2/2 |
| Q-T1S2 | none | 580 / 405 / 222 / 90 / 39 | 32.8 | 2/2/1.25 | 2/2/1.75 |
| Q-T1S2 | R1 | 580 / 420 / 269 / 157 / 95 | 56.1 | 2/2/1.75 | 2/2/2 |
| Q-T1S2 | R2 | 580 / 411 / 244 / 119 / 53 | 33.6 | 2/2/1.5 | 2/2/2 |
| Q-T1S2 | R3 | 580 / 421 / 270 / 176 / 110 | 35.7 | 2/2/2 | 2/2/2 |
| Q-T1S3 | none | 297 / 193 / 103 / 52 / 24 | 6.3 | 2/1.75/1 | 2/2/1.75 |
| Q-T1S3 | R1 | 297 / 204 / 143 / 114 / 86 | 14.0 | 2/2/1.5 | 2/2/2 |
| Q-T1S3 | R2 | 297 / 193 / 116 / 61 / 32 | 6.3 | 2/2/1 | 2/2/2 |
| Q-T1S3 | R3 | 297 / 197 / 126 / 83 / 55 | 7.1 | 2/2/1.25 | 2/2/2 |
| Q-T2S1 | none | 472 / 359 / 196 / 90 / 36 | 24.6 | 2/2/1.25 | 2/2/1.75 |
| Q-T2S1 | R1 | 472 / 369 / 241 / 146 / 88 | 61.0 | 2/2/1.75 | 2/2/2 |
| Q-T2S1 | R2 | 472 / 364 / 215 / 103 / 48 | 25.2 | 2/2/1.5 | 2/2/2 |
| Q-T2S1 | R3 | 472 / 369 / 244 / 142 / 85 | 27.5 | 2/2/1.75 | 2/2/2 |
| Q-T2S2 | none | 714 / 525 / 270 / 96 / 36 | 16.9 | 2/2/1.25 | 2/2/1.75 |
| Q-T2S2 | R1 | 714 / 533 / 335 / 168 / 90 | 48.6 | 2/2/1.75 | 2/2/2 |
| Q-T2S2 | R2 | 714 / 526 / 295 / 125 / 52 | 16.9 | 2/2/1.5 | 2/2/2 |
| Q-T2S2 | R3 | 714 / 531 / 340 / 181 / 97 | 22.3 | 2/2/1.75 | 2/2/2 |
| Q-T2S3 | none | 559 / 394 / 222 / 89 / 44 | 18.6 | 2/2/1.25 | 2/1.75/1.5 |
| Q-T2S3 | R1 | 559 / 400 / 253 / 149 / 89 | 43.4 | 2/2/1.75 | 2/2/2 |
| Q-T2S3 | R2 | 559 / 396 / 229 / 110 / 46 | 18.6 | 2/2/1.5 | 2/2/2 |
| Q-T2S3 | R3 | 559 / 398 / 248 / 147 / 87 | 18.6 | 2/2/1.75 | 2/2/2 |
| Q-T3S1 | none | 1154 / 850 / 375 / 125 / 45 | 109.2 | 2/2/1.5 | 2/2/1.75 |
| Q-T3S1 | R1 | 1154 / 865 / 473 / 194 / 84 | 163.4 | 2/2/1.75 | 2/2/2 |
| Q-T3S1 | R2 | 1154 / 855 / 433 / 190 / 77 | 129.7 | 2/2/1.75 | 2/2/2 |
| Q-T3S1 | R3 | 1154 / 876 / 493 / 250 / 139 | 141.5 | 2/2/2 | 2/2/2 |
| Q-T3S2 | none | 642 / 445 / 204 / 76 / 31 | 88.7 | 2/2/1.25 | 2/2/1.75 |
| Q-T3S2 | R1 | 642 / 456 / 258 / 136 / 81 | 150.8 | 2/2/1.75 | 2/2/2 |
| Q-T3S2 | R2 | 642 / 445 / 218 / 93 / 47 | 106.6 | 2/2/1.25 | 2/2/2 |
| Q-T3S2 | R3 | 642 / 456 / 256 / 143 / 86 | 112.3 | 2/2/1.75 | 2/2/2 |
| Q-T3S3 | none | 614 / 337 / 145 / 53 / 24 | 34.7 | 2/1.75/1 | 2/2/1.75 |
| Q-T3S3 | R1 | 614 / 355 / 196 / 115 / 75 | 83.9 | 2/2/1.5 | 2/2/2 |
| Q-T3S3 | R2 | 614 / 338 / 148 / 61 / 31 | 34.7 | 2/2/1 | 2/2/2 |
| Q-T3S3 | R3 | 614 / 343 / 172 / 80 / 48 | 44.8 | 2/2/1.25 | 2/2/2 |

**Measurement.**
- V weights without clipping: the median stock-day ESS falls from N at a = 0 to 3 to 20 at a = 1.5
  (below 10 in 10 of 18 cells).
- a_max at threshold 30 is 0.75 to 1.25 (threshold 10: 0.25 to 1.75). At thresholds 30 and 100 it does
  not exist in N-T1S1, whose median N is 18.
- At a = 1.5, R1 gives the largest median ESS in 11 of 18 cells and R3 in the other 7; R3 is above R2 in
  every cell. R1's gain is partly the N dependence of section 2: the short stock-days are clipped hardest.
- K weights (q^a) reach a_max(30) = 1.5 to 2 (2 is the grid top) against 0.75 to 1.25 for V weights, as
  expected from q^a against q^(2a).
- Pooling a stock's 8 days lowers the unclipped a_max(10) in 2 cells (Q-T2S3, Q-T3S1: 1.5 per day, 1.25
  pooled), where the 8 days share one dominant run. With R2 or R3 the pooled a_max is at least the
  per-day one in all 18 cells.

## 4. GT under R2 clipping (NASDAQ cells)

**Choice.**
- In the simulation, R2 is the cap exp(mu + 3.0902 sigma_l), the 99.9 % quantile of the cell's log-normal
  size model. Sizes are i.i.d., so this is what a stock's 8-day 99.9th percentile estimates, and it does
  not depend on N.
- The cap applies only to the sizes the moments see. The generator's impact uses the raw q.
- The same seed as 3a is used, so the clipped and unclipped runs share their days exactly.
- a in {0, 0.25, ..., 1.5} is reported, and the depth is min over a <= 1.5 of zeta_K minus the smaller of
  zeta_K(0) and zeta_K(1.5). The 3a depth on the full grid (ends 0 and 3) is listed too.
- Same stopping rule as 3a: SE < 0.005 on local slopes with a <= 1 and T <= 128, cap 64 days.

| cell | cap (shares) | days | domain SE (met) | clipped zK(a) [8,64], a = 0..1.5 | depth (SE), argmin | unclipped zK(a) [8,64] | depth (SE), argmin | 3a full-grid depth |
|---|---|---|---|---|---|---|---|---|
| Q-T1S1 | 4378 | 64 | 0.0054 (no) | 0.917 0.888 0.859 0.833 0.813 0.799 0.790 | 0.000 (0.000), 1.50 | 0.917 0.888 0.859 0.833 0.814 0.807 0.809 | -0.003 (0.001), 1.25 | -0.070 |
| Q-T1S2 | 4851 | 64 | 0.0102 (no) | 0.840 0.811 0.786 0.769 0.760 0.757 0.759 | -0.001 (0.000), 1.25 | 0.840 0.811 0.787 0.774 0.776 0.793 0.817 | -0.043 (0.002), 0.75 | -0.066 |
| Q-T1S3 | 5003 | 64 | 0.0327 (no) | 0.794 0.768 0.750 0.744 0.747 0.754 0.762 | -0.018 (0.001), 0.75 | 0.794 0.769 0.755 0.760 0.785 0.821 0.859 | -0.039 (0.001), 0.50 | -0.039 |
| Q-T2S1 | 4945 | 64 | 0.0403 (no) | 0.686 0.673 0.665 0.664 0.668 0.675 0.684 | -0.020 (0.001), 0.75 | 0.686 0.673 0.667 0.669 0.682 0.702 0.727 | -0.019 (0.001), 0.50 | -0.019 |
| Q-T2S2 | 6803 | 64 | 0.0405 (no) | 0.699 0.684 0.678 0.681 0.690 0.701 0.712 | -0.021 (0.001), 0.50 | 0.699 0.685 0.682 0.693 0.717 0.751 0.787 | -0.017 (0.001), 0.50 | -0.017 |
| Q-T2S3 | 9129 | 64 | 0.0052 (no) | 0.815 0.784 0.761 0.751 0.751 0.755 0.762 | -0.011 (0.000), 1.00 | 0.815 0.784 0.765 0.764 0.783 0.815 0.851 | -0.050 (0.001), 0.75 | -0.050 |
| Q-T3S1 | 7953 | 64 | 0.0130 (no) | 0.758 0.737 0.721 0.710 0.706 0.707 0.710 | -0.003 (0.001), 1.00 | 0.758 0.737 0.722 0.714 0.716 0.727 0.746 | -0.032 (0.002), 0.75 | -0.044 |
| Q-T3S2 | 15957 | 64 | 0.0089 (no) | 0.817 0.787 0.763 0.749 0.743 0.743 0.746 | -0.003 (0.000), 1.00 | 0.817 0.787 0.765 0.755 0.761 0.780 0.806 | -0.051 (0.002), 0.75 | -0.061 |
| Q-T3S3 | 48543 | 64 | 0.0042 (yes) | 0.979 0.923 0.874 0.842 0.828 0.823 0.824 | -0.001 (0.000), 1.25 | 0.979 0.923 0.876 0.855 0.865 0.892 0.924 | -0.069 (0.003), 0.75 | -0.123 |

| cell | clipped zK(a) [16,256] | depth (SE), argmin | unclipped zK(a) [16,256] | depth (SE), argmin | 3a full-grid depth | max zK SE clipped / unclipped |
|---|---|---|---|---|---|---|
| Q-T1S1 | 0.867 0.847 0.825 0.804 0.786 0.773 0.764 | 0.000 (0.000), 1.50 | 0.867 0.847 0.825 0.804 0.788 0.782 0.786 | -0.004 (0.002), 1.25 | -0.058 | 0.0035 / 0.0037 |
| Q-T1S2 | 0.831 0.806 0.781 0.761 0.749 0.743 0.742 | 0.000 (0.000), 1.50 | 0.831 0.806 0.782 0.766 0.766 0.781 0.806 | -0.040 (0.003), 1.00 | -0.065 | 0.0043 / 0.0043 |
| Q-T1S3 | 0.765 0.743 0.726 0.720 0.723 0.730 0.738 | -0.018 (0.003), 0.75 | 0.765 0.743 0.731 0.740 0.769 0.810 0.851 | -0.034 (0.003), 0.50 | -0.034 | 0.0084 / 0.0084 |
| Q-T2S1 | 0.621 0.611 0.606 0.608 0.616 0.627 0.638 | -0.015 (0.003), 0.50 | 0.621 0.611 0.608 0.616 0.635 0.663 0.695 | -0.012 (0.003), 0.50 | -0.012 | 0.0088 / 0.0088 |
| Q-T2S2 | 0.641 0.630 0.627 0.634 0.648 0.662 0.676 | -0.014 (0.003), 0.50 | 0.641 0.631 0.633 0.651 0.685 0.726 0.768 | -0.010 (0.002), 0.25 | -0.010 | 0.0086 / 0.0086 |
| Q-T2S3 | 0.816 0.785 0.760 0.744 0.739 0.741 0.746 | -0.006 (0.001), 1.00 | 0.816 0.786 0.763 0.759 0.776 0.809 0.846 | -0.057 (0.003), 0.75 | -0.057 | 0.0039 / 0.0043 |
| Q-T3S1 | 0.762 0.740 0.721 0.706 0.697 0.694 0.695 | -0.000 (0.001), 1.25 | 0.762 0.740 0.721 0.709 0.707 0.715 0.732 | -0.025 (0.003), 1.00 | -0.055 | 0.0056 / 0.0056 |
| Q-T3S2 | 0.840 0.809 0.781 0.759 0.746 0.740 0.739 | 0.000 (0.000), 1.50 | 0.840 0.809 0.782 0.764 0.762 0.777 0.801 | -0.038 (0.003), 1.00 | -0.077 | 0.0045 / 0.0045 |
| Q-T3S3 | 0.996 0.950 0.901 0.862 0.838 0.828 0.824 | 0.000 (0.000), 1.50 | 0.996 0.950 0.901 0.871 0.872 0.896 0.925 | -0.054 (0.003), 0.75 | -0.117 | 0.0029 / 0.0035 |

**Measurement.**
- **Clipping at the 99.9 % quantile removes most of the V.** On [8, 64] the depth shrinks in 7 of 9
  cells, e.g. Q-T3S3 from -0.069 to -0.001, Q-T3S2 from -0.051 to -0.003, Q-T1S2 from -0.043 to -0.001.
- A V of at least 0.01 survives in 4 cells: Q-T1S3 -0.018, Q-T2S1 -0.020, Q-T2S2 -0.021, Q-T2S3 -0.011.
  In Q-T2S1 and Q-T2S2 it is slightly deeper than unclipped. Where a V remains, its minimum moves later by
  0 to 0.25.
- On [16, 256] the same 4 cells keep a V (-0.006 to -0.018). In the other 5 the clipped curve is flat or
  still falling at a = 1.5.
- The clipped and unclipped curves agree to 0.005 for a <= 0.5 and separate after. Clipping flattens the
  rising branch, which the oracle note attributes to the largest runs (`docs/phase3_oracle.md`, steps 3
  and 5). Capping them removes that branch's weight.
- The depth SEs are 0.000 to 0.003, so even depths of -0.001 to -0.003 are "significant" by the 2-SE rule
  while being negligible in size.
- The SE stopping rule was met in 1 of 9 cells (Q-T3S3). Elsewhere the largest domain slope SE at 64 days
  is 0.005 to 0.040. The zeta_K SE is at most 0.009 for a <= 1.5 on both ranges.
- On the a <= 1.5 grid, the unclipped Q-T1S1 V nearly disappears (depth -0.003 against -0.070 on the full
  grid): its rising branch lies mostly beyond a = 1.5.

## 5. Wall time and peak memory

| step | wall | peak RAM (main process) | notes |
|---|---|---|---|
| collect (guard, 8 days) | 4.8 s | 0.64 GB | one process |
| signs (18 cells, run twice) | 0.5 s | 0.09 GB | 4 workers |
| est (18 cells x 2 generators x 50 panels) | 84 s | 0.13 GB | 4 workers, each about the same size |
| clip and ESS | 3.7 s | 0.16 GB | 4 workers |
| GT under clipping (9 cells x 64 days of 10^6) | 6.6 min | 0.39 GB | 4 workers, about 0.3 GB each as measured in 3a |
| test suite | 44 s | | 270 passed, 28 skipped |

Total compute about 9 min; nothing over 10 min.

## 6. Tests and planted bugs

`tools/est/test_phase3b.py`: 6 checks and 7 planted bugs, all caught. Each mutant was rerun outside pytest
and fails on the intended assertion, not a crash. The 3a calibration refactor (its power-law rule factored
into `calib3.sign_acf_fit`) reproduces `calibration.tsv` byte for byte.

| planted bug | mutated function | caught by |
|---|---|---|
| Kish ESS without the square of the sum | `clip3b.kish` | check_kish_ess (equal weights no longer give n) |
| R1 on the whole sample's volume | `clip3b.clip_r1` | check_r1_own_volume |
| R2 quantile per cell, not per stock | `clip3b.stock_caps` | check_r2_per_stock |
| centred C(l) on the global mean | `signs3b.centred_acf` | check_centred_acf_within_day (constant days keep C = 1) |
| day mean drawn per run, not per stock-day | `signs3b.simulate_two_component` | check_two_component_recovery (tau 0.06 for 0.2) |
| fit model without the 1 + tau^2 normalisation | `signs3b._grid_fit` | check_two_component_recovery (H 0.67, tau 0.42 for 0.7, 0.45) |
| edge flag tests only the upper edge | `signs3b.fit_two_component` | check_edge_flag (H = 0.3 source fitted at 0.505, no flag) |
