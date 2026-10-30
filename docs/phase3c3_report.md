# Phase 3c-3 report: stationary days and the final rho* (NASDAQ cells)

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken here. Real data entered
only the semi-synthetic smoke check (exploration signs and sizes of Q-T3S1 through the guard; no prices).
The fixed decisions are implemented:
- stationary days;
- truth with burn-in;
- v pooled as the mean of stock logs and k~ as the log of the matched mean;
- R_{1/2}/R_0 counted under G0 and GT only.

New and changed code:
- `panel3c.py`: `make_panel(burn=B)`, `_stationary_panel`, `_stationary_group`, `cut_days`,
  `stationary_stock`, `batched_multi_conv`, and the semi-synthetic `fresh` flag.
- `pooled.py`: `zeta_mean_of_logs`, now the default for v.
- `gen3.py`: a stable `fgn_rho`, and `signs` accepts the two-component model.
- `truth3.py`: a `burn` parameter.
- `phase3c3.py` (driver) and `test_3c3.py`.

Tables: `data/phase3/p3c3_*.tsv`. Checkpoints are in `$ITCH_DATA/phase3c3/`.

## Headlines

1. **Stationary days fix the finite-day problem of zeta_K (Measurement).**
   - The step 1 check passes: G0's zeta_K(0) and zeta_V(0) agree across N = 300, 3000 and 10^6 within
     2 SE.
   - In the panels, G0 and GT zeta_K, the P-F slope and R_{1/2}/R_0 now pass |bias| < 0.02 and coverage
     >= 0.85 up to the largest ceiling (1/5) in most cells. Their remaining failures are at the smallest
     ceilings, where few stocks contribute.
   - Same Q-T1S2 panels, G0, ceiling 1/10: zeta_K has bias 0.004 and coverage 0.97 with stationary days,
     against 0.053 and 0.10 with fresh starts.
2. **Still no cell-level rho* (Measurement). The binding estimator is now zeta_V, through coverage.**
   - Across all cells and ceilings, the median of zeta_V's lowest coverage is 0.65 for every generator,
     and 88 % to 90 % of (cell, ceiling) points fall below 0.85. Its bias is borderline (median 0.02).
   - With continuous days, the slowly varying sign component carries across a stock's days, so the days
     are no longer independent clusters, and the day-and-stock intervals understate zeta_V's sampling
     variance. The step 1 check shows the same: zeta_V's SE grows about 20-fold once days are continuous.
   - Under GF, zeta_K is biased (median of the largest |bias| 0.066; the worst elements are negative, on
     [16, 256]), and the P-F slope too (0.039).
3. **The specified burn-in rule is not finite (Measurement), so B = 20,000 (Choice).**
   - The rule's B is infinite for every GT and GF case: their kernels have permanent parts (beta = 0
     bins; the informed runs).
   - For G0 it is infinite in 4 cells, because the calibrated sign ACF never falls below 1 % of C(1) (the
     per-day component).
   - Where the G0 rule is finite it asks 1.2 x 10^7 to 9.9 x 10^7 runs per synthetic stock, beyond this
     machine.
   - The floor B = 20,000 passes the stationarity check, so it is used everywhere.
4. **Run cost and an incident.**
   - You chose 200 panels per cell and generator, not 400, after the reduced run projected 3 h 15 min for
     400.
   - Midway the machine began swapping hard (17.7 of 18.4 GB swap; other applications held about 11 GB).
     I stopped the study at 441 of 780 chunks and, as you chose, resumed with 2 workers. Every chunk is
     checkpointed, so the panels are identical to an uninterrupted run.
   - Total compute for the phase was about 3 h 10 min, over the 2 h target.

## 1. Stationarity check (G0, Q-T1S3, zeta at [8, 64], a = 0)

Each entry comes from 8 independent replicates. A replicate is one series of B + D x N runs, the first B
discarded and the rest cut into D days of N runs (D = 1000, 200 and 4 for N = 300, 3000 and 10^6). Day
means are pooled within the replicate, and the SE comes from a delete-one-replicate jackknife. "Fresh" is
the 3c path: every day its own series, no history.

| B | N | zeta_K(0) (SE) | zeta_V(0) (SE) |
|---|---|---|---|
| fresh | 300 | 0.949 (0.005) | 1.359 (0.004) |
| fresh | 3000 | 0.868 (0.005) | 1.365 (0.005) |
| fresh | 10^6 (3c truth) | 0.812 | |
| 20,000 | 300 | 0.830 (0.015) | 1.401 (0.112) |
| 20,000 | 3000 | 0.815 (0.004) | 1.370 (0.046) |
| 20,000 | 10^6 | 0.807 (0.002) | 1.330 (0.073) |
| 200,000 | 300 | 0.814 (0.009) | 1.323 (0.049) |
| 200,000 | 3000 | 0.816 (0.003) | 1.404 (0.082) |
| 200,000 | 10^6 | 0.808 (0.001) | 1.457 (0.062) |

**Measurement.**
- With burn-in, zeta_K(0) agrees across N within 2 SE at both B. The largest gap, 0.023 between N = 300 and
  10^6 at B = 20,000, is 1.5 SE.
- zeta_V(0) agrees too, but its SE is 0.05 to 0.11 against 0.004 with fresh starts, because consecutive
  days share the slowly varying sign component.
- Fresh starts move zeta_K(0) by +0.14 between N = 10^6 and N = 300.

The first attempt at N = 10^6 failed: the circulant embedding of a series of 4 x 10^6 runs had 14,859
negative eigenvalues. The cause was numerical. `gen3.fgn_rho` computed 0.5 ((l+1)^2H - 2 l^2H + (l-1)^2H),
which cancels catastrophically at l = 10^6 with H near 1. It now uses the equivalent stable form
0.5 l^2H [expm1(2H log1p(1/l)) + expm1(2H log1p(-1/l))]. That form agrees with the old one to 6 x 10^-11
at short lags, and `data/phase3/calibration.tsv` is unchanged byte for byte.

## 2. The burn-in rule

B = max(20,000, 20 x the longest lag at which the kernel or the sign ACF is still above 1 % of its lag-1
value), with lags scanned on a log grid up to 10^12:

| cell | gen | sign ACF 1 % lag | kernel 1 % lag | B by the rule | B used |
|---|---|---|---|---|---|
| Q-T1S1 | G0 | 1.12e+05 | 5.96e+05 | 1.19e+07 | 20000 |
| Q-T1S1 | GT | 1.12e+05 | never (inf) | never (inf) | 20000 |
| Q-T1S1 | GF | 4.77e+05 | never (inf) | never (inf) | 20000 |
| Q-T1S2 | G0 | 1.7e+06 | 1.72e+05 | 3.41e+07 | 20000 |
| Q-T1S2 | GT | 1.7e+06 | never (inf) | never (inf) | 20000 |
| Q-T1S2 | GF | 1.18e+07 | never (inf) | never (inf) | 20000 |
| Q-T1S3 | G0 | never (inf) | 2.93e+04 | never (inf) | 20000 |
| Q-T1S3 | GT | never (inf) | never (inf) | never (inf) | 20000 |
| Q-T1S3 | GF | never (inf) | never (inf) | never (inf) | 20000 |
| Q-T2S1 | G0 | never (inf) | 2.19e+04 | never (inf) | 20000 |
| Q-T2S1 | GT | never (inf) | never (inf) | never (inf) | 20000 |
| Q-T2S1 | GF | never (inf) | never (inf) | never (inf) | 20000 |
| Q-T2S2 | G0 | never (inf) | 2.19e+04 | never (inf) | 20000 |
| Q-T2S2 | GT | never (inf) | never (inf) | never (inf) | 20000 |
| Q-T2S2 | GF | never (inf) | never (inf) | never (inf) | 20000 |
| Q-T2S3 | G0 | 2.07e+04 | 5.96e+05 | 1.19e+07 | 20000 |
| Q-T2S3 | GT | 2.07e+04 | never (inf) | never (inf) | 20000 |
| Q-T2S3 | GF | 7.81e+04 | never (inf) | never (inf) | 20000 |
| Q-T3S1 | G0 | 4.94e+06 | 1.15e+05 | 9.87e+07 | 20000 |
| Q-T3S1 | GT | 4.94e+06 | never (inf) | never (inf) | 20000 |
| Q-T3S1 | GF | 2.35e+07 | never (inf) | never (inf) | 20000 |
| Q-T3S2 | G0 | 7.74e+05 | 2.29e+05 | 1.55e+07 | 20000 |
| Q-T3S2 | GT | 7.74e+05 | never (inf) | never (inf) | 20000 |
| Q-T3S2 | GF | 2.37e+06 | never (inf) | never (inf) | 20000 |
| Q-T3S3 | G0 | 1.15e+05 | 1.24e+06 | 2.48e+07 | 20000 |
| Q-T3S3 | GT | 1.15e+05 | never (inf) | never (inf) | 20000 |
| Q-T3S3 | GF | 3.99e+05 | never (inf) | never (inf) | 20000 |

**Choice:** B = 20,000 for every generator and cell (the rule's floor), justified by the stationarity
check. The truth uses the same B. Each truth day is the last 10^6 runs of its own series of B + 10^6
runs. The days are independent, which gives the same expectations as one series cut into days, with
independent days for the SE.

**Measurement:** the stationary truth equals the 3c fresh truth within 0.002 on zeta_K [8, 64]; a 10^6-run
day already carries a long history.

| sign model | cell | gen | days | stationary zeta_K(a) [8,64], a = 0..1 | fresh (3c) zeta_K(a) [8,64] | stationary zeta_V(0) [8,64] (SE) | P-F slope [8,64] | R1/2/R0 [8,64] |
|---|---|---|---|---|---|---|---|---|
| main | Q-T1S1 | G0 | 64 | 1.049 1.021 0.992 0.963 0.936 | 1.049 1.022 0.993 0.964 0.937 | 1.387 (0.002) | -0.114 | 1.128 |
| main | Q-T1S1 | GT | 64 | 0.916 0.887 0.858 0.832 0.811 | 0.917 0.888 0.859 0.833 0.813 | 1.387 (0.002) | -0.106 | 1.152 |
| main | Q-T1S1 | GF | 64 | 1.000 0.999 0.998 0.997 0.996 | 0.999 0.998 0.998 0.999 1.000 | 1.400 (0.003) | -0.004 | 1.630 |
| main | Q-T1S2 | G0 | 64 | 0.931 0.902 0.873 0.845 0.822 | 0.932 0.904 0.874 0.847 0.823 | 1.338 (0.005) | -0.110 | 1.193 |
| main | Q-T1S2 | GT | 64 | 0.839 0.810 0.785 0.767 0.759 | 0.840 0.811 0.786 0.769 0.760 | 1.338 (0.005) | -0.082 | 1.231 |
| main | Q-T1S2 | GF | 64 | 1.001 1.000 0.999 0.997 0.996 | 0.998 0.998 0.998 0.999 1.000 | 1.347 (0.006) | -0.005 | 1.734 |
| main | Q-T1S3 | G0 | 64 | 0.811 0.784 0.758 0.736 0.718 | 0.812 0.786 0.760 0.737 0.720 | 1.379 (0.022) | -0.094 | 1.267 |
| main | Q-T1S3 | GT | 64 | 0.793 0.767 0.749 0.743 0.746 | 0.794 0.768 0.750 0.744 0.747 | 1.379 (0.022) | -0.047 | 1.322 |
| main | Q-T1S3 | GF | 64 | 1.009 1.006 1.004 1.001 0.999 | 1.000 0.999 0.999 0.999 1.000 | 1.394 (0.026) | -0.009 | 1.897 |
| main | Q-T2S1 | G0 | 64 | 0.722 0.711 0.700 0.690 0.682 | 0.722 0.711 0.700 0.691 0.683 | 1.236 (0.021) | -0.040 | 1.183 |
| main | Q-T2S1 | GT | 64 | 0.686 0.673 0.665 0.663 0.667 | 0.686 0.673 0.665 0.664 0.668 | 1.236 (0.021) | -0.019 | 1.192 |
| main | Q-T2S1 | GF | 64 | 1.007 1.006 1.004 1.003 1.002 | 0.997 0.997 0.997 0.998 0.999 | 1.245 (0.022) | -0.005 | 1.488 |
| main | Q-T2S2 | G0 | 64 | 0.722 0.708 0.695 0.684 0.675 | 0.722 0.708 0.695 0.684 0.675 | 1.243 (0.022) | -0.048 | 1.239 |
| main | Q-T2S2 | GT | 64 | 0.699 0.684 0.678 0.681 0.690 | 0.699 0.684 0.678 0.681 0.690 | 1.243 (0.022) | -0.009 | 1.270 |
| main | Q-T2S2 | GF | 64 | 1.007 1.005 1.004 1.002 1.001 | 0.998 0.998 0.998 0.999 1.000 | 1.251 (0.023) | -0.006 | 1.631 |
| main | Q-T2S3 | G0 | 64 | 0.928 0.896 0.865 0.837 0.816 | 0.929 0.897 0.865 0.838 0.817 | 1.272 (0.002) | -0.113 | 1.247 |
| main | Q-T2S3 | GT | 64 | 0.814 0.783 0.760 0.750 0.750 | 0.815 0.784 0.761 0.751 0.751 | 1.272 (0.002) | -0.064 | 1.313 |
| main | Q-T2S3 | GF | 64 | 1.001 1.000 0.999 0.998 0.997 | 0.998 0.997 0.997 0.998 1.000 | 1.280 (0.003) | -0.004 | 1.806 |
| main | Q-T3S1 | G0 | 64 | 0.854 0.833 0.814 0.796 0.780 | 0.854 0.834 0.815 0.797 0.781 | 1.286 (0.006) | -0.074 | 1.177 |
| main | Q-T3S1 | GT | 64 | 0.757 0.736 0.719 0.709 0.705 | 0.758 0.737 0.721 0.710 0.706 | 1.286 (0.006) | -0.053 | 1.192 |
| main | Q-T3S1 | GF | 64 | 1.002 1.001 1.001 1.000 0.999 | 0.999 0.998 0.997 0.998 0.998 | 1.293 (0.007) | -0.004 | 1.510 |
| main | Q-T3S2 | G0 | 64 | 0.921 0.891 0.861 0.835 0.813 | 0.922 0.892 0.862 0.836 0.814 | 1.330 (0.004) | -0.109 | 1.220 |
| main | Q-T3S2 | GT | 64 | 0.816 0.786 0.762 0.748 0.742 | 0.817 0.787 0.763 0.749 0.743 | 1.330 (0.004) | -0.075 | 1.261 |
| main | Q-T3S2 | GF | 64 | 1.001 1.000 1.000 0.999 0.998 | 0.999 0.997 0.997 0.998 0.998 | 1.338 (0.005) | -0.002 | 1.678 |
| main | Q-T3S3 | G0 | 48 | 1.094 1.039 0.982 0.933 0.894 | 1.096 1.041 0.984 0.934 0.896 | 1.434 (0.002) | -0.203 | 1.284 |
| main | Q-T3S3 | GT | 48 | 0.976 0.921 0.872 0.840 0.826 | 0.979 0.923 0.874 0.842 0.828 | 1.434 (0.002) | -0.152 | 1.398 |
| main | Q-T3S3 | GF | 64 | 1.000 0.999 0.998 0.997 0.997 | 1.000 0.999 0.999 1.000 1.001 | 1.443 (0.003) | -0.003 | 2.125 |
| tail | Q-T1S2 | G0 | 64 | 0.909 0.867 0.831 0.802 0.780 | 0.907 0.865 0.829 0.800 0.779 | 1.414 (0.024) | -0.129 | 1.352 |
| tail | Q-T1S2 | GT | 64 | 0.816 0.774 0.746 0.733 0.730 | 0.811 0.769 0.742 0.730 0.729 | 1.414 (0.024) | -0.085 | 1.411 |
| tail | Q-T1S2 | GF | 64 | 0.988 0.991 0.993 0.995 0.997 | 0.991 0.991 0.993 0.995 0.996 | 1.379 (0.024) | 0.008 | 1.608 |
| tail | Q-T3S1 | G0 | 64 | 0.848 0.820 0.795 0.775 0.759 | 0.845 0.817 0.793 0.773 0.758 | 1.348 (0.024) | -0.089 | 1.277 |
| tail | Q-T3S1 | GT | 64 | 0.752 0.723 0.703 0.693 0.691 | 0.746 0.718 0.699 0.690 0.689 | 1.348 (0.024) | -0.061 | 1.297 |
| tail | Q-T3S1 | GF | 64 | 0.991 0.993 0.995 0.997 0.998 | 0.992 0.992 0.993 0.994 0.996 | 1.313 (0.022) | 0.007 | 1.429 |

## 3. Semi-synthetic generator (for Phase 3d)

**Choice:** a stock's response is one convolution over its study days concatenated in calendar order, after
a burn-in made by prepending a copy of the same days, which is then discarded. The impact of a run
therefore reaches later days.
- `semi_synthetic(fresh=True)` keeps the 3c path.
- Tests: day 2's mid depends on day 1's runs (and not with `fresh=True`); the flow passes through
  unchanged; GT with c = 0 equals G0, both ways.
- **Measurement:** on the real Q-T3S1 exploration flow (22 stocks x 8 days, 368,203 runs, read through the
  guard), signs and sizes pass through bit-identical for G0, GT and GF.

## 4. rho* (200 stationary panels per cell and generator)

Panels are as in 3c-2: confirmation-sized, with exploration day lengths, the 3a sign model (the 3b tail-lag
model for the sensitivity cells), log-normal sizes, and the R2 clip on the moments' sizes only. Each stock
is one series with B = 20,000 burn-in, cut into its 8 days.

The ceiling estimator (T <= rho N_s) is evaluated at 1/160, 1/80, 1/40, 1/20, 1/10 and 1/5, on ranges
[4, 32], [8, 64] and [16, 256] with a = 0..1. The interval is the wider of the percentile bootstrap
(B = 200) and the t day jackknife, at 90 %.

**Choices.**
- An element is evaluable when at least 50 % of panels give an estimate.
- An undefined interval counts as not covering.
- rho* = the largest ceiling at which every evaluable element of every criterion estimator passes (per
  estimator, and per cell over all generators). With this definition a large ceiling can pass while a
  smaller one fails; the binding table shows every ceiling.
- The Monte Carlo SE of each coverage is sqrt(c (1 - c) / 200), about 0.025 at 0.85.

Largest passing ceiling per estimator ("none" = no ceiling passes):

| sign model | cell | G0 zeta_V | G0 zeta_K | G0 pf | G0 Rratio | GT zeta_V | GT zeta_K | GT pf | GT Rratio | GF zeta_V | GF zeta_K | GF pf | GF Rratio | rho* (cell) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| main | Q-T1S1 | 1/80 | 1/5 | 1/5 | 1/5 | 1/80 | 1/5 | 1/5 | 1/5 | none | none | none | none (not counted) | none |
| main | Q-T1S2 | none | 1/5 | 1/5 | 1/20 | none | 1/5 | 1/5 | 1/20 | none | none | none | 1/160 (not counted) | none |
| main | Q-T1S3 | none | 1/5 | 1/5 | none | none | 1/5 | 1/5 | 1/160 | none | 1/80 | 1/80 | 1/40 (not counted) | none |
| main | Q-T2S1 | none | 1/5 | 1/5 | 1/40 | none | 1/5 | 1/5 | 1/80 | none | 1/80 | 1/80 | 1/80 (not counted) | none |
| main | Q-T2S2 | none | 1/5 | 1/5 | 1/20 | none | 1/10 | 1/5 | 1/40 | none | 1/80 | 1/80 | 1/40 (not counted) | none |
| main | Q-T2S3 | 1/5 | 1/10 | 1/5 | 1/5 | 1/5 | none | 1/10 | 1/5 | 1/5 | 1/20 | 1/20 | 1/40 (not counted) | none |
| main | Q-T3S1 | 1/80 | 1/5 | 1/5 | 1/5 | none | 1/5 | 1/5 | 1/5 | 1/80 | 1/160 | 1/10 | 1/80 (not counted) | none |
| main | Q-T3S2 | none | 1/5 | 1/5 | 1/5 | none | 1/5 | 1/5 | 1/5 | none | none | none | none (not counted) | none |
| main | Q-T3S3 | 1/10 | 1/5 | 1/5 | 1/5 | 1/10 | 1/5 | 1/5 | 1/40 | 1/10 | none | none | none (not counted) | none |
| tail | Q-T1S2 | none | 1/20 | 1/5 | none | none | 1/5 | 1/5 | none | none | 1/160 | 1/5 | 1/20 (not counted) | none |
| tail | Q-T3S1 | none | 1/10 | 1/5 | 1/40 | none | 1/5 | 1/5 | 1/10 | none | 1/20 | 1/5 | 1/10 (not counted) | none |

Binding constraint at each ceiling (failing elements among the evaluable ones, the worst bias, and the
lowest coverage with its MC SE):

| sign model | cell | ceiling | passes | not evaluable | failing elements | worst bias | lowest coverage (MC SE) |
|---|---|---|---|---|---|---|---|
| main | Q-T1S1 | 1/160 | no | - | 29 of 90 | GF zeta_K(0)[8,64] -0.042 | GF zeta_V(0.5)[8,64] 0.695 (0.033) |
| main | Q-T1S1 | 1/80 | no | - | 9 of 90 | GF zeta_K(0)[8,64] -0.049 | GF zeta_V(0.5)[8,64] 0.835 (0.026) |
| main | Q-T1S1 | 1/40 | no | - | 22 of 105 | GF zeta_K(0)[16,256] -0.074 | G0 zeta_V(0)[16,256] 0.685 (0.033) |
| main | Q-T1S1 | 1/20 | no | - | 18 of 105 | GF pf[16,256] +0.044 | G0 zeta_V(0.25)[16,256] 0.795 (0.029) |
| main | Q-T1S1 | 1/10 | no | - | 18 of 105 | GF zeta_K(0)[16,256] -0.066 | G0 zeta_V(0)[16,256] 0.800 (0.028) |
| main | Q-T1S1 | 1/5 | no | - | 35 of 105 | GF zeta_K(0.25)[16,256] -0.100 | G0 zeta_V(0.25)[16,256] 0.735 (0.031) |
| main | Q-T1S2 | 1/160 | no | - | 30 of 55 | GF zeta_K(0)[8,64] -0.034 | G0 zeta_K(0)[8,64] 0.550 (0.035) |
| main | Q-T1S2 | 1/80 | no | - | 40 of 90 | GF zeta_K(0)[16,256] -0.042 | GT zeta_V(0.25)[8,64] 0.635 (0.034) |
| main | Q-T1S2 | 1/40 | no | - | 23 of 90 | GF zeta_K(0)[16,256] -0.041 | GT zeta_V(0.25)[8,64] 0.745 (0.031) |
| main | Q-T1S2 | 1/20 | no | - | 43 of 105 | GF zeta_K(0.25)[16,256] -0.057 | GT zeta_V(0)[16,256] 0.560 (0.035) |
| main | Q-T1S2 | 1/10 | no | - | 50 of 105 | GF zeta_K(0.25)[16,256] -0.051 | GT zeta_V(0.25)[8,64] 0.490 (0.035) |
| main | Q-T1S2 | 1/5 | no | - | 49 of 105 | GF zeta_K(0)[16,256] -0.041 | GF zeta_V(0)[16,256] 0.125 (0.023) |
| main | Q-T1S3 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 1 of 2 | G0 Rratio[4,32] -0.009 | G0 Rratio[4,32] 0.840 (0.026) |
| main | Q-T1S3 | 1/80 | no | G0 zeta_V, GT zeta_V, GF zeta_V | 14 of 22 | G0 zeta_K(0)[4,32] +0.086 | G0 zeta_K(0)[4,32] 0.000 (0.000) |
| main | Q-T1S3 | 1/40 | no | - | 33 of 55 | G0 zeta_K(0)[8,64] +0.049 | GF zeta_V(0)[4,32] 0.535 (0.035) |
| main | Q-T1S3 | 1/20 | no | - | 44 of 90 | GF zeta_K(0)[16,256] -0.095 | GF zeta_V(0)[4,32] 0.250 (0.031) |
| main | Q-T1S3 | 1/10 | no | - | 45 of 90 | GF zeta_K(0.25)[16,256] -0.080 | GF zeta_V(0)[8,64] 0.055 (0.016) |
| main | Q-T1S3 | 1/5 | no | - | 59 of 105 | GF zeta_V(0)[16,256] -0.143 | GF zeta_V(0)[8,64] 0.010 (0.007) |
| main | Q-T2S1 | 1/160 | no | G0 zeta_V, GT zeta_V, GF zeta_V | 12 of 22 | GT zeta_K(0)[4,32] +0.071 | G0 zeta_K(0)[4,32] 0.000 (0.000) |
| main | Q-T2S1 | 1/80 | no | - | 25 of 55 | GT zeta_K(0)[8,64] +0.042 | G0 zeta_K(0.25)[8,64] 0.505 (0.035) |
| main | Q-T2S1 | 1/40 | no | - | 48 of 90 | GF zeta_K(0)[16,256] -0.048 | G0 zeta_V(0)[8,64] 0.545 (0.035) |
| main | Q-T2S1 | 1/20 | no | - | 39 of 90 | GF zeta_K(0)[16,256] -0.063 | GT zeta_V(0)[8,64] 0.355 (0.034) |
| main | Q-T2S1 | 1/10 | no | - | 55 of 105 | GT zeta_V(0)[16,256] -0.114 | GT zeta_V(0)[16,256] 0.210 (0.029) |
| main | Q-T2S1 | 1/5 | no | - | 56 of 105 | GT zeta_V(0)[16,256] -0.119 | GT zeta_V(0)[16,256] 0.055 (0.016) |
| main | Q-T2S2 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 0 of 2 | - | - |
| main | Q-T2S2 | 1/80 | no | G0 zeta_V, GT zeta_V, GF zeta_V | 12 of 22 | GT zeta_K(0)[4,32] +0.069 | G0 zeta_K(0)[4,32] 0.080 (0.019) |
| main | Q-T2S2 | 1/40 | no | - | 33 of 55 | GF zeta_K(0)[8,64] -0.057 | GF zeta_V(0)[4,32] 0.605 (0.035) |
| main | Q-T2S2 | 1/20 | no | - | 49 of 90 | GF zeta_K(0)[16,256] -0.076 | GF zeta_V(0)[8,64] 0.540 (0.035) |
| main | Q-T2S2 | 1/10 | no | - | 41 of 90 | GF zeta_K(0.25)[16,256] -0.110 | G0 zeta_V(0)[8,64] 0.415 (0.035) |
| main | Q-T2S2 | 1/5 | no | - | 56 of 105 | GF zeta_K(0)[16,256] -0.131 | GF zeta_V(0)[16,256] 0.185 (0.027) |
| main | Q-T2S3 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, G0 Rratio, GT zeta_V, GT zeta_K, GT pf, GT Rratio, GF zeta_V, GF zeta_K, GF pf | 0 of 0 | - | - |
| main | Q-T2S3 | 1/80 | no | G0 zeta_V, G0 zeta_K, G0 pf, G0 Rratio, GT zeta_V, GT zeta_K, GT pf, GT Rratio, GF zeta_V, GF zeta_K, GF pf | 0 of 0 | - | - |
| main | Q-T2S3 | 1/40 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 0 of 2 | - | - |
| main | Q-T2S3 | 1/20 | no | G0 zeta_V, GT zeta_V, GF zeta_V | 12 of 22 | GT zeta_K(0)[4,32] +0.057 | G0 zeta_K(0)[4,32] 0.310 (0.033) |
| main | Q-T2S3 | 1/10 | no | - | 5 of 55 | GF zeta_K(0)[8,64] -0.046 | GT zeta_K(0)[8,64] 0.950 (0.015) |
| main | Q-T2S3 | 1/5 | no | - | 19 of 90 | GF zeta_K(0)[8,64] -0.077 | GF pf[16,256] 0.810 (0.028) |
| main | Q-T3S1 | 1/160 | no | G0 zeta_V, GT zeta_V, GF zeta_V | 12 of 22 | GT zeta_K(0)[4,32] +0.059 | G0 zeta_K(0)[4,32] 0.000 (0.000) |
| main | Q-T3S1 | 1/80 | no | - | 10 of 55 | GT zeta_K(0)[8,64] +0.023 | G0 zeta_K(0.5)[8,64] 0.835 (0.026) |
| main | Q-T3S1 | 1/40 | no | - | 21 of 90 | GF zeta_K(0)[16,256] -0.038 | GF zeta_V(0)[8,64] 0.790 (0.029) |
| main | Q-T3S1 | 1/20 | no | - | 21 of 90 | GF zeta_K(0)[16,256] -0.090 | GF zeta_V(0.25)[8,64] 0.740 (0.031) |
| main | Q-T3S1 | 1/10 | no | - | 39 of 105 | GF zeta_K(0.25)[16,256] -0.076 | GF zeta_V(0)[16,256] 0.610 (0.034) |
| main | Q-T3S1 | 1/5 | no | - | 41 of 105 | GF zeta_K(0)[16,256] -0.090 | GT zeta_V(0)[16,256] 0.430 (0.035) |
| main | Q-T3S2 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 0 of 2 | - | - |
| main | Q-T3S2 | 1/80 | no | G0 zeta_V, GT zeta_V, GF zeta_V | 18 of 22 | GT zeta_K(0)[4,32] +0.052 | G0 zeta_K(0)[4,32] 0.420 (0.035) |
| main | Q-T3S2 | 1/40 | no | - | 24 of 55 | GF zeta_K(0.5)[8,64] -0.074 | GF zeta_V(0.25)[4,32] 0.700 (0.032) |
| main | Q-T3S2 | 1/20 | no | - | 39 of 90 | GF zeta_K(0)[8,64] -0.075 | GF pf[16,256] 0.635 (0.034) |
| main | Q-T3S2 | 1/10 | no | - | 30 of 90 | GF zeta_K(0.5)[16,256] -0.085 | GF zeta_V(0.25)[8,64] 0.795 (0.029) |
| main | Q-T3S2 | 1/5 | no | - | 35 of 105 | GF zeta_K(0)[8,64] -0.090 | GF zeta_V(0)[16,256] 0.590 (0.035) |
| main | Q-T3S3 | 1/160 | no | G0 zeta_V, G0 zeta_K, G0 pf, GT zeta_V, GT zeta_K, GT pf, GF zeta_V, GF zeta_K, GF pf | 1 of 2 | G0 Rratio[4,32] -0.005 | G0 Rratio[4,32] 0.785 (0.029) |
| main | Q-T3S3 | 1/80 | no | G0 zeta_V, GT zeta_V, GF zeta_V | 16 of 22 | GF zeta_K(0)[4,32] -0.074 | G0 zeta_K(0)[4,32] 0.230 (0.030) |
| main | Q-T3S3 | 1/40 | no | - | 24 of 55 | GF zeta_K(0.25)[8,64] -0.080 | GF zeta_K(0)[8,64] 0.715 (0.032) |
| main | Q-T3S3 | 1/20 | no | - | 40 of 90 | GF zeta_K(0)[8,64] -0.094 | GF zeta_K(0)[16,256] 0.535 (0.035) |
| main | Q-T3S3 | 1/10 | no | - | 22 of 90 | GF zeta_K(0)[8,64] -0.100 | GF zeta_K(0)[16,256] 0.780 (0.029) |
| main | Q-T3S3 | 1/5 | no | - | 31 of 105 | GF zeta_K(0)[8,64] -0.166 | G0 zeta_V(0)[16,256] 0.655 (0.034) |
| tail | Q-T1S2 | 1/160 | no | - | 26 of 55 | GT zeta_V(0)[4,32] -0.025 | GT zeta_V(0)[4,32] 0.440 (0.035) |
| tail | Q-T1S2 | 1/80 | no | - | 33 of 90 | GT zeta_V(0)[8,64] -0.047 | GT zeta_V(0)[8,64] 0.410 (0.035) |
| tail | Q-T1S2 | 1/40 | no | - | 36 of 90 | GF zeta_K(0)[16,256] -0.053 | G0 zeta_V(0)[8,64] 0.505 (0.035) |
| tail | Q-T1S2 | 1/20 | no | - | 52 of 105 | GT zeta_V(0)[16,256] -0.093 | G0 zeta_V(0)[8,64] 0.315 (0.033) |
| tail | Q-T1S2 | 1/10 | no | - | 52 of 105 | G0 zeta_V(0)[16,256] -0.089 | G0 zeta_V(0)[8,64] 0.105 (0.022) |
| tail | Q-T1S2 | 1/5 | no | - | 54 of 105 | GT zeta_V(0)[16,256] -0.089 | G0 zeta_V(0)[16,256] 0.000 (0.000) |
| tail | Q-T3S1 | 1/160 | no | G0 zeta_V, GT zeta_V, GF zeta_V | 10 of 22 | GT zeta_K(0)[4,32] +0.024 | GT zeta_K(0.75)[4,32] 0.500 (0.035) |
| tail | Q-T3S1 | 1/80 | no | - | 10 of 55 | GT zeta_V(0)[4,32] -0.013 | GT zeta_V(0)[4,32] 0.725 (0.032) |
| tail | Q-T3S1 | 1/40 | no | - | 22 of 90 | GT zeta_V(0)[8,64] -0.031 | GT zeta_V(0)[8,64] 0.665 (0.033) |
| tail | Q-T3S1 | 1/20 | no | - | 23 of 90 | GT zeta_V(0)[8,64] -0.033 | GT zeta_V(0)[8,64] 0.635 (0.034) |
| tail | Q-T3S1 | 1/10 | no | - | 41 of 105 | GT zeta_V(0)[16,256] -0.082 | GT zeta_V(0)[16,256] 0.410 (0.035) |
| tail | Q-T3S1 | 1/5 | no | - | 43 of 105 | GT zeta_V(0)[16,256] -0.085 | GT zeta_V(0)[16,256] 0.210 (0.029) |

Horizons dropped because the pooled k~ was not positive (panels with at least one, of 200):

| sign model | cell | gen | panels with a dropped T (any ceiling) | panels |
|---|---|---|---|---|
| main | Q-T1S1 | G0 | 0 | 200 |
| main | Q-T1S1 | GT | 0 | 200 |
| main | Q-T1S1 | GF | 66 | 200 |
| main | Q-T1S2 | G0 | 0 | 200 |
| main | Q-T1S2 | GT | 0 | 200 |
| main | Q-T1S2 | GF | 74 | 200 |
| main | Q-T1S3 | G0 | 1 | 200 |
| main | Q-T1S3 | GT | 0 | 200 |
| main | Q-T1S3 | GF | 54 | 200 |
| main | Q-T2S1 | G0 | 0 | 200 |
| main | Q-T2S1 | GT | 0 | 200 |
| main | Q-T2S1 | GF | 47 | 200 |
| main | Q-T2S2 | G0 | 5 | 200 |
| main | Q-T2S2 | GT | 11 | 200 |
| main | Q-T2S2 | GF | 71 | 200 |
| main | Q-T2S3 | G0 | 0 | 200 |
| main | Q-T2S3 | GT | 0 | 200 |
| main | Q-T2S3 | GF | 14 | 200 |
| main | Q-T3S1 | G0 | 0 | 200 |
| main | Q-T3S1 | GT | 0 | 200 |
| main | Q-T3S1 | GF | 49 | 200 |
| main | Q-T3S2 | G0 | 1 | 200 |
| main | Q-T3S2 | GT | 4 | 200 |
| main | Q-T3S2 | GF | 85 | 200 |
| main | Q-T3S3 | G0 | 0 | 200 |
| main | Q-T3S3 | GT | 1 | 200 |
| main | Q-T3S3 | GF | 105 | 200 |
| tail | Q-T1S2 | G0 | 0 | 200 |
| tail | Q-T1S2 | GT | 4 | 200 |
| tail | Q-T1S2 | GF | 60 | 200 |
| tail | Q-T3S1 | G0 | 0 | 200 |
| tail | Q-T3S1 | GT | 1 | 200 |
| tail | Q-T3S1 | GF | 31 | 200 |

### Fresh-start against stationary (Q-T1S2 and Q-T3S1, 200 panels each)

For each ceiling: max |bias| / min coverage of stationary panels against the stationary truth, of
fresh-start panels against the fresh (3c) truth, and of fresh-start panels against the stationary truth.

| cell | gen | ceiling | estimator | stationary panels vs stationary truth: max abs bias / min cov | fresh panels vs fresh truth | fresh panels vs stationary truth |
|---|---|---|---|---|---|---|
| Q-T1S2 | G0 | 1/160 | zeta_V | 0.002 / 0.64 | 0.002 / 0.85 | 0.003 / 0.84 |
| Q-T1S2 | G0 | 1/160 | zeta_K | 0.028 / 0.55 | 0.031 / 0.45 | 0.032 / 0.42 |
| Q-T1S2 | G0 | 1/160 | pf | 0.010 / 0.85 | 0.013 / 0.79 | 0.012 / 0.79 |
| Q-T1S2 | G0 | 1/160 | Rratio | 0.002 / 0.94 | 0.005 / 0.90 | 0.005 / 0.94 |
| Q-T1S2 | G0 | 1/80 | zeta_V | 0.006 / 0.65 | 0.004 / 0.86 | 0.006 / 0.84 |
| Q-T1S2 | G0 | 1/80 | zeta_K | 0.026 / 0.73 | 0.021 / 0.79 | 0.017 / 0.82 |
| Q-T1S2 | G0 | 1/80 | pf | 0.007 / 0.89 | 0.006 / 0.88 | 0.005 / 0.88 |
| Q-T1S2 | G0 | 1/80 | Rratio | 0.006 / 0.85 | 0.011 / 0.78 | 0.010 / 0.84 |
| Q-T1S2 | G0 | 1/40 | zeta_V | 0.008 / 0.76 | 0.004 / 0.91 | 0.005 / 0.91 |
| Q-T1S2 | G0 | 1/40 | zeta_K | 0.020 / 0.90 | 0.016 / 0.62 | 0.017 / 0.63 |
| Q-T1S2 | G0 | 1/40 | pf | 0.004 / 0.96 | 0.008 / 0.94 | 0.007 / 0.94 |
| Q-T1S2 | G0 | 1/40 | Rratio | 0.005 / 0.89 | 0.013 / 0.54 | 0.012 / 0.58 |
| Q-T1S2 | G0 | 1/20 | zeta_V | 0.021 / 0.57 | 0.011 / 0.82 | 0.014 / 0.81 |
| Q-T1S2 | G0 | 1/20 | zeta_K | 0.003 / 0.97 | 0.031 / 0.28 | 0.035 / 0.29 |
| Q-T1S2 | G0 | 1/20 | pf | 0.001 / 0.96 | 0.013 / 0.84 | 0.013 / 0.87 |
| Q-T1S2 | G0 | 1/20 | Rratio | 0.008 / 0.89 | 0.017 / 0.31 | 0.017 / 0.41 |
| Q-T1S2 | G0 | 1/10 | zeta_V | 0.026 / 0.50 | 0.011 / 0.81 | 0.014 / 0.71 |
| Q-T1S2 | G0 | 1/10 | zeta_K | 0.004 / 0.97 | 0.053 / 0.10 | 0.058 / 0.09 |
| Q-T1S2 | G0 | 1/10 | pf | 0.003 / 0.97 | 0.018 / 0.74 | 0.018 / 0.78 |
| Q-T1S2 | G0 | 1/10 | Rratio | 0.009 / 0.80 | 0.021 / 0.22 | 0.021 / 0.26 |
| Q-T1S2 | G0 | 1/5 | zeta_V | 0.033 / 0.21 | 0.015 / 0.51 | 0.018 / 0.37 |
| Q-T1S2 | G0 | 1/5 | zeta_K | 0.008 / 0.96 | 0.074 / 0.01 | 0.078 / 0.01 |
| Q-T1S2 | G0 | 1/5 | pf | 0.003 / 0.98 | 0.022 / 0.69 | 0.022 / 0.70 |
| Q-T1S2 | G0 | 1/5 | Rratio | 0.011 / 0.77 | 0.025 / 0.20 | 0.025 / 0.24 |
| Q-T1S2 | GT | 1/160 | zeta_V | 0.004 / 0.66 | 0.001 / 0.85 | 0.002 / 0.84 |
| Q-T1S2 | GT | 1/160 | zeta_K | 0.030 / 0.63 | 0.034 / 0.55 | 0.035 / 0.53 |
| Q-T1S2 | GT | 1/160 | pf | 0.010 / 0.83 | 0.013 / 0.82 | 0.013 / 0.83 |
| Q-T1S2 | GT | 1/160 | Rratio | 0.002 / 0.96 | 0.004 / 0.95 | 0.003 / 0.96 |
| Q-T1S2 | GT | 1/80 | zeta_V | 0.009 / 0.64 | 0.003 / 0.84 | 0.005 / 0.84 |
| Q-T1S2 | GT | 1/80 | zeta_K | 0.025 / 0.80 | 0.020 / 0.81 | 0.013 / 0.83 |
| Q-T1S2 | GT | 1/80 | pf | 0.019 / 0.83 | 0.010 / 0.82 | 0.009 / 0.83 |
| Q-T1S2 | GT | 1/80 | Rratio | 0.003 / 0.87 | 0.011 / 0.84 | 0.010 / 0.86 |
| Q-T1S2 | GT | 1/40 | zeta_V | 0.009 / 0.74 | 0.003 / 0.94 | 0.005 / 0.91 |
| Q-T1S2 | GT | 1/40 | zeta_K | 0.020 / 0.97 | 0.012 / 0.91 | 0.013 / 0.92 |
| Q-T1S2 | GT | 1/40 | pf | 0.013 / 0.94 | 0.007 / 0.94 | 0.007 / 0.96 |
| Q-T1S2 | GT | 1/40 | Rratio | 0.005 / 0.92 | 0.010 / 0.68 | 0.009 / 0.76 |
| Q-T1S2 | GT | 1/20 | zeta_V | 0.023 / 0.56 | 0.008 / 0.82 | 0.011 / 0.81 |
| Q-T1S2 | GT | 1/20 | zeta_K | 0.013 / 0.98 | 0.021 / 0.73 | 0.025 / 0.76 |
| Q-T1S2 | GT | 1/20 | pf | 0.008 / 0.99 | 0.012 / 0.90 | 0.012 / 0.92 |
| Q-T1S2 | GT | 1/20 | Rratio | 0.006 / 0.89 | 0.018 / 0.53 | 0.018 / 0.62 |
| Q-T1S2 | GT | 1/10 | zeta_V | 0.028 / 0.49 | 0.010 / 0.83 | 0.014 / 0.75 |
| Q-T1S2 | GT | 1/10 | zeta_K | 0.002 / 0.97 | 0.045 / 0.60 | 0.052 / 0.56 |
| Q-T1S2 | GT | 1/10 | pf | 0.001 / 0.98 | 0.025 / 0.80 | 0.026 / 0.82 |
| Q-T1S2 | GT | 1/10 | Rratio | 0.010 / 0.80 | 0.024 / 0.41 | 0.024 / 0.49 |
| Q-T1S2 | GT | 1/5 | zeta_V | 0.035 / 0.17 | 0.014 / 0.56 | 0.017 / 0.41 |
| Q-T1S2 | GT | 1/5 | zeta_K | 0.006 / 0.99 | 0.069 / 0.47 | 0.077 / 0.44 |
| Q-T1S2 | GT | 1/5 | pf | 0.002 / 0.97 | 0.030 / 0.74 | 0.031 / 0.74 |
| Q-T1S2 | GT | 1/5 | Rratio | 0.013 / 0.77 | 0.027 / 0.30 | 0.027 / 0.38 |
| Q-T1S2 | GF | 1/160 | zeta_V | 0.006 / 0.69 | 0.000 / 0.83 | 0.001 / 0.81 |
| Q-T1S2 | GF | 1/160 | zeta_K | 0.034 / 0.86 | 0.029 / 0.85 | 0.032 / 0.86 |
| Q-T1S2 | GF | 1/160 | pf | 0.026 / 0.86 | 0.019 / 0.88 | 0.026 / 0.88 |
| Q-T1S2 | GF | 1/160 | Rratio | 0.009 / 0.96 | 0.028 / 0.96 | 0.032 / 0.96 |
| Q-T1S2 | GF | 1/80 | zeta_V | 0.012 / 0.68 | 0.002 / 0.83 | 0.003 / 0.83 |
| Q-T1S2 | GF | 1/80 | zeta_K | 0.042 / 0.82 | 0.050 / 0.83 | 0.053 / 0.83 |
| Q-T1S2 | GF | 1/80 | pf | 0.027 / 0.81 | 0.035 / 0.83 | 0.042 / 0.83 |
| Q-T1S2 | GF | 1/80 | Rratio | 0.044 / 0.90 | 0.158 / 0.87 | 0.169 / 0.88 |
| Q-T1S2 | GF | 1/40 | zeta_V | 0.011 / 0.78 | 0.004 / 0.90 | 0.005 / 0.88 |
| Q-T1S2 | GF | 1/40 | zeta_K | 0.041 / 0.95 | 0.089 / 0.93 | 0.073 / 0.93 |
| Q-T1S2 | GF | 1/40 | pf | 0.023 / 0.94 | 0.049 / 0.94 | 0.038 / 0.94 |
| Q-T1S2 | GF | 1/40 | Rratio | 0.336 / 0.96 | 0.055 / 0.96 | 0.043 / 0.97 |
| Q-T1S2 | GF | 1/20 | zeta_V | 0.028 / 0.60 | 0.009 / 0.81 | 0.011 / 0.80 |
| Q-T1S2 | GF | 1/20 | zeta_K | 0.057 / 0.97 | 0.069 / 0.96 | 0.047 / 0.96 |
| Q-T1S2 | GF | 1/20 | pf | 0.033 / 0.97 | 0.027 / 0.96 | 0.016 / 0.96 |
| Q-T1S2 | GF | 1/20 | Rratio | 0.038 / 0.98 | 0.042 / 0.97 | 0.053 / 0.97 |
| Q-T1S2 | GF | 1/10 | zeta_V | 0.032 / 0.49 | 0.012 / 0.74 | 0.015 / 0.65 |
| Q-T1S2 | GF | 1/10 | zeta_K | 0.051 / 0.98 | 0.047 / 0.97 | 0.034 / 0.97 |
| Q-T1S2 | GF | 1/10 | pf | 0.026 / 0.98 | 0.011 / 0.98 | 0.022 / 0.98 |
| Q-T1S2 | GF | 1/10 | Rratio | 0.032 / 0.98 | 0.023 / 0.99 | 0.020 / 0.99 |
| Q-T1S2 | GF | 1/5 | zeta_V | 0.038 / 0.12 | 0.016 / 0.48 | 0.019 / 0.36 |
| Q-T1S2 | GF | 1/5 | zeta_K | 0.041 / 0.98 | 0.053 / 0.98 | 0.032 / 0.98 |
| Q-T1S2 | GF | 1/5 | pf | 0.036 / 0.98 | 0.031 / 0.99 | 0.020 / 0.99 |
| Q-T1S2 | GF | 1/5 | Rratio | 0.024 / 0.98 | 0.024 / 0.97 | 0.022 / 0.98 |
| Q-T3S1 | G0 | 1/160 | zeta_V | - / - | - / - | - / - |
| Q-T3S1 | G0 | 1/160 | zeta_K | 0.051 / 0.00 | 0.056 / 0.00 | 0.056 / 0.00 |
| Q-T3S1 | G0 | 1/160 | pf | 0.017 / 0.45 | 0.020 / 0.23 | 0.019 / 0.25 |
| Q-T3S1 | G0 | 1/160 | Rratio | 0.000 / 0.94 | 0.005 / 0.95 | 0.004 / 0.96 |
| Q-T3S1 | G0 | 1/80 | zeta_V | 0.004 / 0.86 | 0.001 / 0.96 | 0.003 / 0.95 |
| Q-T3S1 | G0 | 1/80 | zeta_K | 0.017 / 0.83 | 0.026 / 0.62 | 0.027 / 0.61 |
| Q-T3S1 | G0 | 1/80 | pf | 0.003 / 0.97 | 0.006 / 0.97 | 0.006 / 0.98 |
| Q-T3S1 | G0 | 1/80 | Rratio | 0.001 / 0.98 | 0.006 / 0.93 | 0.006 / 0.93 |
| Q-T3S1 | G0 | 1/40 | zeta_V | 0.011 / 0.83 | 0.004 / 0.95 | 0.006 / 0.95 |
| Q-T3S1 | G0 | 1/40 | zeta_K | 0.031 / 0.84 | 0.017 / 0.86 | 0.018 / 0.87 |
| Q-T3S1 | G0 | 1/40 | pf | 0.014 / 0.94 | 0.007 / 0.97 | 0.008 / 0.97 |
| Q-T3S1 | G0 | 1/40 | Rratio | 0.005 / 0.96 | 0.013 / 0.89 | 0.012 / 0.89 |
| Q-T3S1 | G0 | 1/20 | zeta_V | 0.012 / 0.81 | 0.004 / 0.94 | 0.006 / 0.93 |
| Q-T3S1 | G0 | 1/20 | zeta_K | 0.017 / 0.94 | 0.023 / 0.77 | 0.024 / 0.78 |
| Q-T3S1 | G0 | 1/20 | pf | 0.006 / 0.96 | 0.009 / 0.94 | 0.008 / 0.95 |
| Q-T3S1 | G0 | 1/20 | Rratio | 0.007 / 0.96 | 0.017 / 0.85 | 0.016 / 0.88 |
| Q-T3S1 | G0 | 1/10 | zeta_V | 0.031 / 0.66 | 0.011 / 0.89 | 0.014 / 0.86 |
| Q-T3S1 | G0 | 1/10 | zeta_K | 0.007 / 0.97 | 0.055 / 0.69 | 0.060 / 0.68 |
| Q-T3S1 | G0 | 1/10 | pf | 0.002 / 0.97 | 0.017 / 0.91 | 0.016 / 0.93 |
| Q-T3S1 | G0 | 1/10 | Rratio | 0.008 / 0.96 | 0.021 / 0.81 | 0.020 / 0.83 |
| Q-T3S1 | G0 | 1/5 | zeta_V | 0.037 / 0.45 | 0.013 / 0.80 | 0.016 / 0.74 |
| Q-T3S1 | G0 | 1/5 | zeta_K | 0.010 / 0.98 | 0.072 / 0.60 | 0.077 / 0.57 |
| Q-T3S1 | G0 | 1/5 | pf | 0.002 / 0.97 | 0.022 / 0.89 | 0.020 / 0.90 |
| Q-T3S1 | G0 | 1/5 | Rratio | 0.008 / 0.95 | 0.024 / 0.82 | 0.023 / 0.85 |
| Q-T3S1 | GT | 1/160 | zeta_V | - / - | - / - | - / - |
| Q-T3S1 | GT | 1/160 | zeta_K | 0.059 / 0.01 | 0.063 / 0.01 | 0.062 / 0.01 |
| Q-T3S1 | GT | 1/160 | pf | 0.019 / 0.49 | 0.023 / 0.34 | 0.022 / 0.38 |
| Q-T3S1 | GT | 1/160 | Rratio | 0.002 / 0.97 | 0.006 / 0.96 | 0.005 / 0.96 |
| Q-T3S1 | GT | 1/80 | zeta_V | 0.005 / 0.84 | 0.002 / 0.96 | 0.004 / 0.96 |
| Q-T3S1 | GT | 1/80 | zeta_K | 0.023 / 0.83 | 0.029 / 0.73 | 0.029 / 0.72 |
| Q-T3S1 | GT | 1/80 | pf | 0.003 / 0.98 | 0.009 / 0.97 | 0.008 / 0.97 |
| Q-T3S1 | GT | 1/80 | Rratio | 0.003 / 0.96 | 0.006 / 0.94 | 0.006 / 0.95 |
| Q-T3S1 | GT | 1/40 | zeta_V | 0.013 / 0.81 | 0.005 / 0.95 | 0.008 / 0.94 |
| Q-T3S1 | GT | 1/40 | zeta_K | 0.022 / 0.94 | 0.016 / 0.96 | 0.013 / 0.96 |
| Q-T3S1 | GT | 1/40 | pf | 0.017 / 0.94 | 0.012 / 0.97 | 0.012 / 0.96 |
| Q-T3S1 | GT | 1/40 | Rratio | 0.005 / 0.96 | 0.014 / 0.93 | 0.013 / 0.94 |
| Q-T3S1 | GT | 1/20 | zeta_V | 0.015 / 0.78 | 0.006 / 0.94 | 0.008 / 0.92 |
| Q-T3S1 | GT | 1/20 | zeta_K | 0.018 / 0.97 | 0.017 / 0.97 | 0.018 / 0.96 |
| Q-T3S1 | GT | 1/20 | pf | 0.010 / 0.98 | 0.010 / 0.96 | 0.009 / 0.97 |
| Q-T3S1 | GT | 1/20 | Rratio | 0.008 / 0.98 | 0.017 / 0.92 | 0.016 / 0.93 |
| Q-T3S1 | GT | 1/10 | zeta_V | 0.034 / 0.64 | 0.014 / 0.88 | 0.017 / 0.83 |
| Q-T3S1 | GT | 1/10 | zeta_K | 0.003 / 0.97 | 0.039 / 0.93 | 0.048 / 0.93 |
| Q-T3S1 | GT | 1/10 | pf | 0.002 / 0.99 | 0.020 / 0.96 | 0.020 / 0.97 |
| Q-T3S1 | GT | 1/10 | Rratio | 0.009 / 0.97 | 0.021 / 0.91 | 0.019 / 0.92 |
| Q-T3S1 | GT | 1/5 | zeta_V | 0.038 / 0.43 | 0.017 / 0.77 | 0.020 / 0.69 |
| Q-T3S1 | GT | 1/5 | zeta_K | 0.003 / 0.96 | 0.057 / 0.93 | 0.066 / 0.92 |
| Q-T3S1 | GT | 1/5 | pf | 0.002 / 0.99 | 0.027 / 0.95 | 0.027 / 0.95 |
| Q-T3S1 | GT | 1/5 | Rratio | 0.010 / 0.96 | 0.024 / 0.90 | 0.023 / 0.91 |
| Q-T3S1 | GF | 1/160 | zeta_V | - / - | - / - | - / - |
| Q-T3S1 | GF | 1/160 | zeta_K | 0.011 / 0.98 | 0.007 / 0.97 | 0.003 / 0.97 |
| Q-T3S1 | GF | 1/160 | pf | 0.009 / 0.99 | 0.002 / 0.99 | 0.005 / 0.99 |
| Q-T3S1 | GF | 1/160 | Rratio | 0.012 / 0.95 | 0.009 / 0.97 | 0.010 / 0.97 |
| Q-T3S1 | GF | 1/80 | zeta_V | 0.006 / 0.85 | 0.002 / 0.98 | 0.003 / 0.97 |
| Q-T3S1 | GF | 1/80 | zeta_K | 0.021 / 0.97 | 0.007 / 0.98 | 0.003 / 0.98 |
| Q-T3S1 | GF | 1/80 | pf | 0.014 / 0.99 | 0.001 / 0.98 | 0.002 / 0.98 |
| Q-T3S1 | GF | 1/80 | Rratio | 0.008 / 0.98 | 0.011 / 0.98 | 0.012 / 0.98 |
| Q-T3S1 | GF | 1/40 | zeta_V | 0.014 / 0.79 | 0.005 / 0.94 | 0.007 / 0.94 |
| Q-T3S1 | GF | 1/40 | zeta_K | 0.038 / 0.98 | 0.019 / 0.98 | 0.009 / 0.98 |
| Q-T3S1 | GF | 1/40 | pf | 0.027 / 0.98 | 0.010 / 0.99 | 0.006 / 0.99 |
| Q-T3S1 | GF | 1/40 | Rratio | 0.142 / 0.99 | 0.041 / 0.98 | 0.043 / 0.97 |
| Q-T3S1 | GF | 1/20 | zeta_V | 0.016 / 0.74 | 0.007 / 0.92 | 0.009 / 0.90 |
| Q-T3S1 | GF | 1/20 | zeta_K | 0.090 / 0.98 | 0.044 / 0.98 | 0.030 / 0.98 |
| Q-T3S1 | GF | 1/20 | pf | 0.064 / 0.99 | 0.034 / 0.99 | 0.024 / 0.99 |
| Q-T3S1 | GF | 1/20 | Rratio | 0.039 / 0.99 | 0.043 / 0.98 | 0.045 / 0.98 |
| Q-T3S1 | GF | 1/10 | zeta_V | 0.037 / 0.61 | 0.016 / 0.84 | 0.018 / 0.82 |
| Q-T3S1 | GF | 1/10 | zeta_K | 0.076 / 0.97 | 0.062 / 0.98 | 0.048 / 0.97 |
| Q-T3S1 | GF | 1/10 | pf | 0.013 / 0.99 | 0.039 / 0.98 | 0.029 / 0.98 |
| Q-T3S1 | GF | 1/10 | Rratio | 0.032 / 0.99 | 0.029 / 0.99 | 0.031 / 0.99 |
| Q-T3S1 | GF | 1/5 | zeta_V | 0.041 / 0.45 | 0.018 / 0.76 | 0.020 / 0.73 |
| Q-T3S1 | GF | 1/5 | zeta_K | 0.090 / 0.98 | 0.093 / 0.99 | 0.080 / 0.98 |
| Q-T3S1 | GF | 1/5 | pf | 0.055 / 0.99 | 0.072 / 0.98 | 0.062 / 0.98 |
| Q-T3S1 | GF | 1/5 | Rratio | 0.026 / 0.99 | 0.055 / 0.99 | 0.057 / 0.99 |

**Measurement.**
- **zeta_K.** Fresh panels are badly off at the larger ceilings: Q-T1S2 G0 at 1/10 has bias 0.053 and
  coverage 0.10. Stationary panels are fine there (0.004, 0.97).
- **zeta_V.** Coverage is lower with stationary days (Q-T1S2 G0 at 1/10: 0.50 against 0.81). Bias is
  similar in size. The cross-day dependence makes zeta_V's true sampling variance larger than the day
  clusters can see.
- **Tail-lag sign model** (sensitivity, Q-T1S2 and Q-T3S1): no rho* either. Its zeta_K passes to 1/20 and
  1/10 for G0 at the two cells.

## 5. Reach (no rho* found: shown at 1/40 and 1/20)

Confirmation stocks (N_min from `data/universe.tsv`) contributing to each range (N_min >= T2 / rho), per
cell and pooled across the sigma_l cells of each tick bin:

| cell | rho | confirm stocks | [4,32] contributing | [8,64] contributing | [16,256] contributing |
|---|---|---|---|---|---|
| Q-T1S1 | 1/40 | 161 | 22 | 4 | 0 |
| Q-T1S2 | 1/40 | 216 | 23 | 5 | 0 |
| Q-T1S3 | 1/40 | 252 | 5 | 0 | 0 |
| Q-T2S1 | 1/40 | 80 | 17 | 6 | 2 |
| Q-T2S2 | 1/40 | 50 | 3 | 0 | 0 |
| Q-T2S3 | 1/40 | 17 | 0 | 0 | 0 |
| Q-T3S1 | 1/40 | 50 | 21 | 5 | 1 |
| Q-T3S2 | 1/40 | 27 | 6 | 3 | 0 |
| Q-T3S3 | 1/40 | 23 | 1 | 1 | 0 |
| Q-T1 (S1+S2+S3) | 1/40 | 629 | 50 | 9 | 0 |
| Q-T2 (S1+S2+S3) | 1/40 | 147 | 20 | 6 | 2 |
| Q-T3 (S1+S2+S3) | 1/40 | 100 | 28 | 9 | 1 |
| Q-T1S1 | 1/20 | 161 | 34 | 22 | 0 |
| Q-T1S2 | 1/20 | 216 | 52 | 23 | 2 |
| Q-T1S3 | 1/20 | 252 | 26 | 5 | 0 |
| Q-T2S1 | 1/20 | 80 | 32 | 17 | 2 |
| Q-T2S2 | 1/20 | 50 | 12 | 3 | 0 |
| Q-T2S3 | 1/20 | 17 | 1 | 0 | 0 |
| Q-T3S1 | 1/20 | 50 | 31 | 21 | 2 |
| Q-T3S2 | 1/20 | 27 | 12 | 6 | 1 |
| Q-T3S3 | 1/20 | 23 | 5 | 1 | 1 |
| Q-T1 (S1+S2+S3) | 1/20 | 629 | 112 | 50 | 2 |
| Q-T2 (S1+S2+S3) | 1/20 | 147 | 45 | 20 | 2 |
| Q-T3 (S1+S2+S3) | 1/20 | 100 | 48 | 28 | 4 |

**Measurement.**
- At 1/20, pooled by tick bin, [4, 32] has 112, 45 and 48 contributing stocks (T1, T2, T3), and [8, 64]
  has 50, 20 and 28.
- [16, 256] has at most 4 stocks in any tick bin.

## 6. Tests and planted bugs

`tools/est/test_3c3.py` has 6 checks and 5 planted bugs, all caught, each by its intended assertion (each
mutant rerun outside pytest). Full suite: 321 passed, 28 skipped.

| planted bug | mutated function | caught by |
|---|---|---|
| burn-in not discarded | `panel3c.cut_days` | check_burnin_discarded |
| days cut from separate series | `panel3c._stationary_group` | check_days_from_one_series (day-to-day correlation of mean signs -0.02 instead of > 0.5) |
| v pooled as the log of the cross-stock mean | `pooled.core2` | check_v_mean_of_logs (1.49 to 1.53 instead of 1.4) |
| R ratio counted under GF | `phase3c3.criterion` | check_criterion |
| semi-synthetic fresh start per day | `panel3c.semi_synthetic` | check_semi_carries_across_days |

Other checks: the stable `fgn_rho` equals the direct formula at short lags and the asymptote at 10^7.
`check_gt_c0_is_g0` and `check_response_uses_raw_sizes` in `test_panel3c.py` now run on both
semi-synthetic paths (the raw-size check on the fresh path).

## 7. Wall time and peak memory

| job | wall | peak RAM | notes |
|---|---|---|---|
| stationarity check (two passes; the first stopped at N = 10^6 by the embedding failure) | 3.2 min | 0.73 GB | 4 workers |
| truth, 33 runs | 24.7 min | 0.49 GB main | estimated 14 min |
| reduced study (10 panels per combination) | 5.3 min | 1.18 GB largest process | before stock grouping; chunks discarded |
| study, first part (441 of 780 chunks) | about 2 h | 1.24 GB all processes | stopped while the machine swapped |
| study, resumed (339 chunks, 2 workers) | 26.4 min | 0.91 GB all processes | |
| analyse | 9 s | 0.04 GB | |
| test suite | 48 s | | |

Grouping a panel's stocks in groups of 32 cut the largest worker from 1.18 GB to about 0.41 GB.

## 8. Open for a decision (not taken here)

- **zeta_V's interval.** The day clusters are not independent under stationary days. The candidates are a
  stock-level (or stock x block) bootstrap, or widening by the across-day dependence. The 3c-2
  stock-cluster jackknife covered worse with fresh days; it is untested here.
- **GF's zeta_K and P-F slope.** They are biased (about -0.04 to -0.09 on [16, 256]): the log of matched
  means of a k~ carried by 10 % of runs is noisy.
- **The "largest passing ceiling" definition** lets a ceiling pass above failing smaller ones. A
  contiguous rule (every ceiling up to rho passes) would give smaller values.
