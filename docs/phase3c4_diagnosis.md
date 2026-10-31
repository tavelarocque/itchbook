# Phase 3c-4 diagnosis: why rho* is none in every cell

Under the pre-registered rules. This diagnosis led to the rule changes of 31 Oct 2026 (range reach, P-F slope by coverage only, tick-bin rho*); see `docs/phase3c4_report.md` section 5.

Read-only diagnosis of the saved Phase 3c-4 study (`~/CODE/itch-data/phase3c4/study/*.npz`, the 3c-3 truth files that `analyse` uses). No simulation was run, and no rule, threshold or estimator was changed. Bias, coverage, MC SE, evaluability and the pass rule are recomputed exactly as `phase3c4.evaluate` computes them: |bias| < 0.02; coverage fails only when coverage + 2 MC SE < 0.85; an element counts only when at least half its panels give a finite estimate; criterion estimators are zeta_V(a), zeta_K(a), P-F slope and R_{1/2}/R_0 under G0 and GT, and the P-F slope alone under GF.

**Stock counts.** The chunks do not store how many stocks are eligible. Each panel's day lengths are therefore re-drawn from that panel's seeded generator (`PC.day_lengths` is the first draw in `study_job`, before any order flow is simulated), and eligibility is `pooled.masks`: stock s enters horizon T under ceiling c when T / min_d N_sd <= c. The reconstruction was checked against the saved estimates: zeta_V on a range is finite if and only if at least one stock is eligible at the range's top horizon, and over all 263,105 finite and 86,895 empty (panel, a, range, ceiling) cases the reconstructed counts agree with the saved finiteness without a single exception; no saved zeta_K is finite with fewer than 3 eligible horizons. "Stocks contributing to a range" below is the number eligible at the range's smallest horizon (every stock that enters any horizon of the range; for zeta_V a stock enters only from the second horizon, since an increment is decided at its upper horizon); "at the top" is the number eligible at T = 32 or 64. For R_{1/2}/R_0 both counts are taken at its single horizon (T = 8 on [4, 32], T = 16 on [8, 64]).

Panels: 400 for every (cell, generator) that was topped up, 200 otherwise, 100 in the continuous-series stress cell. Confirmation stocks per cell: Q-T1S1 161, Q-T1S2 216, Q-T1S3 252, Q-T2S1 80, Q-T2S2 50, Q-T2S3 17, Q-T3S1 50, Q-T3S2 27, Q-T3S3 23.

## 0. Headline

- Independent-days cells: 1917 evaluable (element, ceiling) pairs, 331 fail (138 on coverage only, 54 on bias only, 139 on both).

Status of each ceiling for the cell (all criterion estimators over G0, GT, GF), and the failing elements by family (V = zeta_V(a), K = zeta_K(a), PF = P-F slope, R = R ratio). `n/e` = some criterion estimator has no evaluable element, so `rho_star` skips the ceiling. `rho_star` stops at the first evaluable ceiling that fails, so a cell gets rho* = none whenever its smallest evaluable ceiling fails, whatever happens above it.

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

What the tables below show (independent-days cells):

- zeta_K(a): all 83 failures are elements whose median number of stocks at the range's top horizon is 0. Such an element is still evaluable (3 kept horizons below the top suffice), so the estimate is a slope over the lower part of the range while the truth is the slope over the whole range. All 83 fail at small ceilings only (1/80 or 1/40) and pass at every larger one; 37 other top-empty elements pass. They are the only failures at the smallest evaluable ceiling in Q-T2S1 and Q-T3S1, and most of them in Q-T1S3 and Q-T2S2; they also fail at the n/e ceilings of Q-T1S3, Q-T2S2, Q-T2S3, Q-T3S2 and Q-T3S3, which do not count toward rho*.
- zeta_V(a): 181 failures, every one with negative bias, largest |bias| 0.029; 121 fail on coverage alone. They begin at 1/40, grow with the ceiling (13, 43, 59, 66 failures at 1/40 ... 1/5), and fall off with a (49 at a = 0, 11 at a = 1). Patterns are large ceilings only (46 elements) or everywhere (20). Coverage falls as the number of stocks rises (section 4); that axis is confounded with the ceiling, since more stocks means a larger T / N.
- P-F slope under GF: 36 failures, all on bias alone (positive, up to +0.097) with coverage 0.85 to 1.00; 9 are top-empty truncations, the other 27 persist with 2 to 72 stocks at the top. GF alone fails Q-T1S1 at every ceiling and Q-T3S3 at every evaluable ceiling, and gives the first failure in Q-T1S2, Q-T2S3 and Q-T3S2.
- R ratio: 20 failures, all negative bias: 15 in Q-T1S3 from 1/40 up (coverage 0.27 to 0.81), 3 in Q-T1S2 at 1/5, 2 in Q-T3S3 (GT, [8, 64]) at 1/10 and 1/5.
- Widest interval: the day-and-stock bootstrap in 64 to 99 % of pairs (lowest when few stocks contribute), the day jackknife-t in the rest; the whole-stock bootstrap is widest in under 1 %.
- Stress cell: ind Q-T3S1 fails only at 1/80 (3 zeta_K elements with no stock at T = 64) and passes 1/40, 1/10, 1/5; cont Q-T3S1 fails at every ceiling. Its zeta_V bias is about twice the ind bias (about -0.017 against -0.008 on [8, 64]), and its coverage is 0.64 to 0.76.

## 1. Failing criterion elements

Per cell and ceiling, every evaluable criterion element that fails. `part` is which half of the rule failed. `pattern` is the element's status at 1/80, 1/40, 1/20, 1/10, 1/5 (F fail, P pass, - not evaluable). `larger ceilings` says whether it also fails at every larger evaluable ceiling.

### ind Q-T1S1

**ceiling 1/80**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0282 | 0.955 | 0.015 | 200 | bias | FFFFF | fails at all 4 larger |

**ceiling 1/40**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0293 | 0.990 | 0.007 | 200 | bias | FFFFF | fails at all 3 larger |

**ceiling 1/20**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0342 | 0.990 | 0.007 | 200 | bias | FFFFF | fails at all 2 larger |

**ceiling 1/10**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0373 | 0.995 | 0.005 | 200 | bias | FFFFF | fails at all 1 larger |

**ceiling 1/5**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0329 | 1.000 | 0.000 | 200 | bias | FFFFF | no larger evaluable |

### ind Q-T1S2

**ceiling 1/80**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0400 | 0.970 | 0.012 | 200 | bias | FFPPP | passes at 1/20, 1/10, 1/5 |

**ceiling 1/40**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0291 | 0.995 | 0.005 | 200 | bias | FFPPP | passes at 1/20, 1/10, 1/5 |

**ceiling 1/20**: 0 failing

**ceiling 1/10**: 9 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0039 | 0.805 | 0.020 | 400 | coverage | PPPFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0039 | 0.802 | 0.020 | 400 | coverage | PPPFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0038 | 0.802 | 0.020 | 400 | coverage | PPPFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0073 | 0.748 | 0.022 | 400 | coverage | PPPFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0070 | 0.755 | 0.022 | 400 | coverage | PPPFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0063 | 0.797 | 0.020 | 400 | coverage | PPPFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0071 | 0.738 | 0.022 | 400 | coverage | PPPFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0070 | 0.743 | 0.022 | 400 | coverage | PPPFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0064 | 0.782 | 0.021 | 400 | coverage | PPPFF | fails at all 1 larger |

**ceiling 1/5**: 17 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0044 | 0.770 | 0.021 | 400 | coverage | PPPFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0044 | 0.772 | 0.021 | 400 | coverage | PPPPF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0041 | 0.787 | 0.020 | 400 | coverage | PPPPF | no larger evaluable |
| [4,32] | G0 | Rratio[4,32] | -0.0058 | 0.802 | 0.020 | 400 | coverage | PPPPF | no larger evaluable |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0045 | 0.725 | 0.022 | 400 | coverage | PPPFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0044 | 0.752 | 0.022 | 400 | coverage | PPPFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0041 | 0.807 | 0.020 | 400 | coverage | PPPPF | no larger evaluable |
| [4,32] | GT | Rratio[4,32] | -0.0078 | 0.790 | 0.020 | 400 | coverage | PPPPF | no larger evaluable |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0084 | 0.603 | 0.024 | 400 | coverage | PPPFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0083 | 0.598 | 0.025 | 400 | coverage | PPPFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0077 | 0.645 | 0.024 | 400 | coverage | PPPFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0066 | 0.723 | 0.022 | 400 | coverage | PPPPF | no larger evaluable |
| [8,64] | G0 | Rratio[8,64] | -0.0068 | 0.805 | 0.020 | 400 | coverage | PPPPF | no larger evaluable |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0083 | 0.585 | 0.025 | 400 | coverage | PPPFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0082 | 0.600 | 0.024 | 400 | coverage | PPPFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0076 | 0.627 | 0.024 | 400 | coverage | PPPFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0066 | 0.713 | 0.023 | 400 | coverage | PPPPF | no larger evaluable |

### ind Q-T1S3

**ceiling 1/80**: 13 failing (33 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_K(0)[4,32] | +0.0856 | 0.000 | 0.000 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.25)[4,32] | +0.0725 | 0.003 | 0.002 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.5)[4,32] | +0.0613 | 0.010 | 0.005 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.75)[4,32] | +0.0527 | 0.018 | 0.007 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(1)[4,32] | +0.0466 | 0.083 | 0.014 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | pf[4,32] | -0.0391 | 0.182 | 0.019 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GF | pf[4,32] | +0.0245 | 0.995 | 0.005 | 200 | bias | FFPPP | passes at 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0)[4,32] | +0.0851 | 0.003 | 0.002 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.25)[4,32] | +0.0713 | 0.007 | 0.004 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.5)[4,32] | +0.0586 | 0.010 | 0.005 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.75)[4,32] | +0.0480 | 0.040 | 0.010 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(1)[4,32] | +0.0400 | 0.142 | 0.017 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | pf[4,32] | -0.0454 | 0.110 | 0.016 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |

**ceiling 1/40**: 17 failing (10 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | Rratio[4,32] | -0.0108 | 0.685 | 0.023 | 400 | coverage | PFFFF | fails at all 3 larger |
| [4,32] | GF | pf[4,32] | +0.0404 | 0.990 | 0.007 | 200 | bias | FFPPP | passes at 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0146 | 0.800 | 0.020 | 400 | coverage | -FFFF | fails at all 3 larger |
| [4,32] | GT | Rratio[4,32] | -0.0143 | 0.655 | 0.024 | 400 | coverage | PFFFF | fails at all 3 larger |
| [8,64] | G0 | zeta_K(0)[8,64] | +0.0479 | 0.565 | 0.025 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.25)[8,64] | +0.0427 | 0.562 | 0.025 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.5)[8,64] | +0.0379 | 0.590 | 0.025 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.75)[8,64] | +0.0339 | 0.682 | 0.023 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(1)[8,64] | +0.0308 | 0.780 | 0.021 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GF | pf[8,64] | +0.0605 | 0.975 | 0.011 | 200 | bias | -FFFF | fails at all 3 larger |
| [8,64] | GT | zeta_K(0)[8,64] | +0.0436 | 0.738 | 0.022 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.25)[8,64] | +0.0382 | 0.745 | 0.022 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.5)[8,64] | +0.0327 | 0.770 | 0.021 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.75)[8,64] | +0.0275 | 0.805 | 0.020 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(1)[8,64] | +0.0231 | 0.870 | 0.017 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | pf[8,64] | -0.0206 | 0.915 | 0.014 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | Rratio[8,64] | -0.0181 | 0.805 | 0.020 | 400 | coverage | PFFFF | fails at all 3 larger |

**ceiling 1/20**: 20 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0155 | 0.482 | 0.025 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0148 | 0.495 | 0.025 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0125 | 0.603 | 0.024 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | G0 | zeta_V(0.75)[4,32] | -0.0085 | 0.755 | 0.022 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | G0 | Rratio[4,32] | -0.0119 | 0.595 | 0.025 | 400 | coverage | PFFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0159 | 0.475 | 0.025 | 400 | coverage | -FFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0151 | 0.505 | 0.025 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0127 | 0.595 | 0.025 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0.75)[4,32] | -0.0086 | 0.760 | 0.021 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | GT | Rratio[4,32] | -0.0161 | 0.472 | 0.025 | 400 | coverage | PFFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0254 | 0.703 | 0.023 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0247 | 0.715 | 0.023 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0215 | 0.738 | 0.022 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | G0 | Rratio[8,64] | -0.0164 | 0.620 | 0.024 | 400 | coverage | PPFFF | fails at all 2 larger |
| [8,64] | GF | pf[8,64] | +0.0575 | 1.000 | 0.000 | 200 | bias | -FFFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0261 | 0.685 | 0.023 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0252 | 0.695 | 0.023 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0218 | 0.728 | 0.022 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0158 | 0.795 | 0.020 | 400 | coverage | --FFF | fails at all 2 larger |
| [8,64] | GT | Rratio[8,64] | -0.0206 | 0.652 | 0.024 | 400 | both | PFFFF | fails at all 2 larger |

**ceiling 1/10**: 23 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0162 | 0.193 | 0.020 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0155 | 0.200 | 0.020 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0132 | 0.310 | 0.023 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | G0 | zeta_V(0.75)[4,32] | -0.0091 | 0.590 | 0.025 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | G0 | Rratio[4,32] | -0.0135 | 0.492 | 0.025 | 400 | coverage | PFFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0164 | 0.200 | 0.020 | 400 | coverage | -FFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0156 | 0.235 | 0.021 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0132 | 0.320 | 0.023 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.75)[4,32] | -0.0091 | 0.585 | 0.025 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | GT | Rratio[4,32] | -0.0188 | 0.340 | 0.024 | 400 | coverage | PFFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0274 | 0.247 | 0.022 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0266 | 0.258 | 0.022 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0232 | 0.328 | 0.023 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0171 | 0.492 | 0.025 | 400 | coverage | --PFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(1)[8,64] | -0.0105 | 0.743 | 0.022 | 400 | coverage | --PFF | fails at all 1 larger |
| [8,64] | G0 | Rratio[8,64] | -0.0175 | 0.517 | 0.025 | 400 | coverage | PPFFF | fails at all 1 larger |
| [8,64] | GF | pf[8,64] | +0.0354 | 0.995 | 0.005 | 200 | bias | -FFFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0278 | 0.225 | 0.021 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0267 | 0.253 | 0.022 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0233 | 0.347 | 0.024 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0170 | 0.540 | 0.025 | 400 | coverage | --FFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(1)[8,64] | -0.0104 | 0.800 | 0.020 | 400 | coverage | --PFF | fails at all 1 larger |
| [8,64] | GT | Rratio[8,64] | -0.0227 | 0.453 | 0.025 | 400 | both | PFFFF | fails at all 1 larger |

**ceiling 1/5**: 25 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0167 | 0.113 | 0.016 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0161 | 0.128 | 0.017 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0139 | 0.198 | 0.020 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.75)[4,32] | -0.0097 | 0.475 | 0.025 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(1)[4,32] | -0.0057 | 0.810 | 0.020 | 400 | coverage | -PPPF | no larger evaluable |
| [4,32] | G0 | Rratio[4,32] | -0.0158 | 0.422 | 0.025 | 400 | coverage | PFFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0170 | 0.105 | 0.015 | 400 | coverage | -FFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0163 | 0.125 | 0.017 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0140 | 0.220 | 0.021 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.75)[4,32] | -0.0099 | 0.465 | 0.025 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | GT | zeta_V(1)[4,32] | -0.0057 | 0.810 | 0.020 | 400 | coverage | -PPPF | no larger evaluable |
| [4,32] | GT | Rratio[4,32] | -0.0221 | 0.265 | 0.022 | 400 | both | PFFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0289 | 0.037 | 0.009 | 400 | both | --FFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0280 | 0.052 | 0.011 | 400 | both | --FFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0247 | 0.092 | 0.014 | 400 | both | --FFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0185 | 0.193 | 0.020 | 400 | coverage | --PFF | no larger evaluable |
| [8,64] | G0 | zeta_V(1)[8,64] | -0.0117 | 0.568 | 0.025 | 400 | coverage | --PFF | no larger evaluable |
| [8,64] | G0 | Rratio[8,64] | -0.0192 | 0.495 | 0.025 | 400 | coverage | PPFFF | no larger evaluable |
| [8,64] | GF | pf[8,64] | +0.0282 | 0.990 | 0.007 | 200 | bias | -FFFF | no larger evaluable |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0291 | 0.033 | 0.009 | 400 | both | --FFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0281 | 0.052 | 0.011 | 400 | both | --FFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0247 | 0.105 | 0.015 | 400 | both | --FFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0183 | 0.253 | 0.022 | 400 | coverage | --FFF | no larger evaluable |
| [8,64] | GT | zeta_V(1)[8,64] | -0.0114 | 0.603 | 0.024 | 400 | coverage | --PFF | no larger evaluable |
| [8,64] | GT | Rratio[8,64] | -0.0256 | 0.372 | 0.024 | 400 | both | PFFFF | no larger evaluable |

### ind Q-T2S1

**ceiling 1/80**: 10 failing (10 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | G0 | zeta_K(0)[8,64] | +0.0415 | 0.522 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.25)[8,64] | +0.0383 | 0.505 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.5)[8,64] | +0.0354 | 0.540 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.75)[8,64] | +0.0328 | 0.603 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(1)[8,64] | +0.0305 | 0.720 | 0.022 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0)[8,64] | +0.0403 | 0.627 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.25)[8,64] | +0.0375 | 0.618 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.5)[8,64] | +0.0345 | 0.640 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.75)[8,64] | +0.0312 | 0.685 | 0.023 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(1)[8,64] | +0.0280 | 0.767 | 0.021 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |

**ceiling 1/40**: 12 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0132 | 0.750 | 0.022 | 400 | coverage | PFFFF | fails at all 3 larger |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0126 | 0.760 | 0.021 | 400 | coverage | PFFFF | fails at all 3 larger |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0108 | 0.775 | 0.021 | 400 | coverage | PFFFF | fails at all 3 larger |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0134 | 0.757 | 0.021 | 400 | coverage | PFFFF | fails at all 3 larger |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0130 | 0.752 | 0.022 | 400 | coverage | PFFFF | fails at all 3 larger |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0112 | 0.797 | 0.020 | 400 | coverage | PFFFF | fails at all 3 larger |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0241 | 0.750 | 0.022 | 400 | both | -FFFF | fails at all 3 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0232 | 0.765 | 0.021 | 400 | both | -FFFF | fails at all 3 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0200 | 0.782 | 0.021 | 400 | both | -FFFF | fails at all 3 larger |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0243 | 0.777 | 0.021 | 400 | both | -FFFF | fails at all 3 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0235 | 0.785 | 0.021 | 400 | both | -FFFF | fails at all 3 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0206 | 0.797 | 0.020 | 400 | both | -FFFF | fails at all 3 larger |

**ceiling 1/20**: 17 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0140 | 0.623 | 0.024 | 400 | coverage | PFFFF | fails at all 2 larger |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0134 | 0.632 | 0.024 | 400 | coverage | PFFFF | fails at all 2 larger |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0115 | 0.682 | 0.023 | 400 | coverage | PFFFF | fails at all 2 larger |
| [4,32] | G0 | zeta_V(0.75)[4,32] | -0.0086 | 0.787 | 0.020 | 400 | coverage | PPFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0139 | 0.623 | 0.024 | 400 | coverage | PFFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0134 | 0.632 | 0.024 | 400 | coverage | PFFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0117 | 0.703 | 0.023 | 400 | coverage | PFFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0.75)[4,32] | -0.0089 | 0.790 | 0.020 | 400 | coverage | PPFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0251 | 0.623 | 0.024 | 400 | both | -FFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0242 | 0.608 | 0.024 | 400 | both | -FFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0210 | 0.645 | 0.024 | 400 | both | -FFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0160 | 0.718 | 0.023 | 400 | coverage | -PFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(1)[8,64] | -0.0107 | 0.797 | 0.020 | 400 | coverage | -PFFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0245 | 0.627 | 0.024 | 400 | both | -FFFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0237 | 0.642 | 0.024 | 400 | both | -FFFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0207 | 0.690 | 0.023 | 400 | both | -FFFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0159 | 0.752 | 0.022 | 400 | coverage | -PFFF | fails at all 2 larger |

**ceiling 1/10**: 18 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0147 | 0.505 | 0.025 | 400 | coverage | PFFFF | fails at all 1 larger |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0141 | 0.517 | 0.025 | 400 | coverage | PFFFF | fails at all 1 larger |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0123 | 0.610 | 0.024 | 400 | coverage | PFFFF | fails at all 1 larger |
| [4,32] | G0 | zeta_V(0.75)[4,32] | -0.0094 | 0.728 | 0.022 | 400 | coverage | PPFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0146 | 0.492 | 0.025 | 400 | coverage | PFFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0141 | 0.510 | 0.025 | 400 | coverage | PFFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0124 | 0.578 | 0.025 | 400 | coverage | PFFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.75)[4,32] | -0.0096 | 0.700 | 0.023 | 400 | coverage | PPFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0261 | 0.390 | 0.024 | 400 | both | -FFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0251 | 0.407 | 0.025 | 400 | both | -FFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0220 | 0.505 | 0.025 | 400 | both | -FFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0170 | 0.637 | 0.024 | 400 | coverage | -PFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(1)[8,64] | -0.0116 | 0.760 | 0.021 | 400 | coverage | -PFFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0254 | 0.443 | 0.025 | 400 | both | -FFFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0246 | 0.470 | 0.025 | 400 | both | -FFFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0217 | 0.535 | 0.025 | 400 | both | -FFFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0168 | 0.620 | 0.024 | 400 | coverage | -PFFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(1)[8,64] | -0.0115 | 0.775 | 0.021 | 400 | coverage | -PPFF | fails at all 1 larger |

**ceiling 1/5**: 18 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0153 | 0.440 | 0.025 | 400 | coverage | PFFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0147 | 0.463 | 0.025 | 400 | coverage | PFFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0129 | 0.557 | 0.025 | 400 | coverage | PFFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.75)[4,32] | -0.0100 | 0.708 | 0.023 | 400 | coverage | PPFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0153 | 0.455 | 0.025 | 400 | coverage | PFFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0147 | 0.468 | 0.025 | 400 | coverage | PFFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0129 | 0.552 | 0.025 | 400 | coverage | PFFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.75)[4,32] | -0.0102 | 0.710 | 0.023 | 400 | coverage | PPFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0275 | 0.253 | 0.022 | 400 | both | -FFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0266 | 0.278 | 0.022 | 400 | both | -FFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0234 | 0.343 | 0.024 | 400 | both | -FFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0183 | 0.482 | 0.025 | 400 | coverage | -PFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(1)[8,64] | -0.0128 | 0.685 | 0.023 | 400 | coverage | -PFFF | no larger evaluable |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0268 | 0.275 | 0.022 | 400 | both | -FFFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0260 | 0.292 | 0.023 | 400 | both | -FFFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0231 | 0.380 | 0.024 | 400 | both | -FFFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0182 | 0.507 | 0.025 | 400 | coverage | -PFFF | no larger evaluable |
| [8,64] | GT | zeta_V(1)[8,64] | -0.0128 | 0.713 | 0.023 | 400 | coverage | -PPFF | no larger evaluable |

### ind Q-T2S2

**ceiling 1/80**: 13 failing (33 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_K(0)[4,32] | +0.0699 | 0.095 | 0.015 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.25)[4,32] | +0.0613 | 0.092 | 0.014 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.5)[4,32] | +0.0539 | 0.138 | 0.017 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.75)[4,32] | +0.0479 | 0.220 | 0.021 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(1)[4,32] | +0.0435 | 0.393 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | pf[4,32] | -0.0265 | 0.675 | 0.023 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GF | pf[4,32] | +0.0343 | 0.975 | 0.011 | 200 | bias | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0)[4,32] | +0.0694 | 0.142 | 0.017 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.25)[4,32] | +0.0609 | 0.155 | 0.018 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.5)[4,32] | +0.0527 | 0.203 | 0.020 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.75)[4,32] | +0.0452 | 0.297 | 0.023 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(1)[4,32] | +0.0391 | 0.463 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | pf[4,32] | -0.0306 | 0.682 | 0.023 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |

**ceiling 1/40**: 11 failing (10 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | G0 | zeta_K(0)[8,64] | +0.0428 | 0.790 | 0.020 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.25)[8,64] | +0.0380 | 0.792 | 0.020 | 400 | both | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.5)[8,64] | +0.0336 | 0.818 | 0.019 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(0.75)[8,64] | +0.0298 | 0.853 | 0.018 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | G0 | zeta_K(1)[8,64] | +0.0267 | 0.897 | 0.015 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GF | pf[8,64] | +0.0202 | 0.950 | 0.015 | 200 | bias | -FFPP | passes at 1/10, 1/5 |
| [8,64] | GT | zeta_K(0)[8,64] | +0.0363 | 0.887 | 0.016 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.25)[8,64] | +0.0329 | 0.897 | 0.015 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.5)[8,64] | +0.0292 | 0.902 | 0.015 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.75)[8,64] | +0.0255 | 0.915 | 0.014 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(1)[8,64] | +0.0221 | 0.930 | 0.013 | 400 | bias | -FPPP | passes at 1/20, 1/10, 1/5 |

**ceiling 1/20**: 12 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0148 | 0.752 | 0.022 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0143 | 0.767 | 0.021 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0122 | 0.800 | 0.020 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0136 | 0.772 | 0.021 | 400 | coverage | -PFFF | fails at all 2 larger |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0130 | 0.807 | 0.020 | 400 | coverage | -PFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0262 | 0.770 | 0.021 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0254 | 0.770 | 0.021 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0219 | 0.802 | 0.020 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | GF | pf[8,64] | +0.0389 | 0.995 | 0.005 | 200 | bias | -FFPP | passes at 1/10, 1/5 |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0243 | 0.775 | 0.021 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0234 | 0.770 | 0.021 | 400 | both | --FFF | fails at all 2 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0202 | 0.828 | 0.019 | 400 | bias | --FFF | fails at all 2 larger |

**ceiling 1/10**: 14 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0151 | 0.688 | 0.023 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0146 | 0.672 | 0.023 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0125 | 0.740 | 0.022 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0140 | 0.705 | 0.023 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0133 | 0.730 | 0.022 | 400 | coverage | -PFFF | fails at all 1 larger |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0112 | 0.795 | 0.020 | 400 | coverage | -PPFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0272 | 0.593 | 0.025 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0263 | 0.615 | 0.024 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0227 | 0.670 | 0.024 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0169 | 0.755 | 0.022 | 400 | coverage | --PFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0251 | 0.660 | 0.024 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0243 | 0.685 | 0.023 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0209 | 0.718 | 0.023 | 400 | both | --FFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0156 | 0.800 | 0.020 | 400 | coverage | --PFF | fails at all 1 larger |

**ceiling 1/5**: 14 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0156 | 0.627 | 0.024 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.25)[4,32] | -0.0151 | 0.642 | 0.024 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | G0 | zeta_V(0.5)[4,32] | -0.0131 | 0.728 | 0.022 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0)[4,32] | -0.0148 | 0.693 | 0.023 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.25)[4,32] | -0.0140 | 0.695 | 0.023 | 400 | coverage | -PFFF | no larger evaluable |
| [4,32] | GT | zeta_V(0.5)[4,32] | -0.0117 | 0.770 | 0.021 | 400 | coverage | -PPFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0280 | 0.455 | 0.025 | 400 | both | --FFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0270 | 0.485 | 0.025 | 400 | both | --FFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0235 | 0.547 | 0.025 | 400 | both | --FFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0177 | 0.680 | 0.023 | 400 | coverage | --PFF | no larger evaluable |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0264 | 0.517 | 0.025 | 400 | both | --FFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0254 | 0.540 | 0.025 | 400 | both | --FFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0220 | 0.632 | 0.024 | 400 | both | --FFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.75)[8,64] | -0.0165 | 0.745 | 0.022 | 400 | coverage | --PFF | no larger evaluable |

### ind Q-T2S3

**ceiling 1/80**: 0 failing (50 criterion elements not evaluable)

**ceiling 1/40**: 0 failing (48 criterion elements not evaluable)

**ceiling 1/20**: 13 failing (33 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_K(0)[4,32] | +0.0498 | 0.355 | 0.024 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.25)[4,32] | +0.0416 | 0.410 | 0.025 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.5)[4,32] | +0.0351 | 0.500 | 0.025 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.75)[4,32] | +0.0304 | 0.645 | 0.024 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | G0 | zeta_K(1)[4,32] | +0.0272 | 0.772 | 0.021 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | G0 | pf[4,32] | -0.0226 | 0.865 | 0.017 | 400 | bias | --FPP | passes at 1/10, 1/5 |
| [4,32] | GF | pf[4,32] | +0.0224 | 0.985 | 0.009 | 200 | bias | --FFF | fails at all 2 larger |
| [4,32] | GT | zeta_K(0)[4,32] | +0.0571 | 0.502 | 0.025 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.25)[4,32] | +0.0477 | 0.545 | 0.025 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.5)[4,32] | +0.0397 | 0.605 | 0.024 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.75)[4,32] | +0.0332 | 0.693 | 0.023 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | GT | zeta_K(1)[4,32] | +0.0281 | 0.810 | 0.020 | 400 | both | --FPP | passes at 1/10, 1/5 |
| [4,32] | GT | pf[4,32] | -0.0290 | 0.855 | 0.018 | 400 | bias | --FPP | passes at 1/10, 1/5 |

**ceiling 1/10**: 1 failing (10 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | GF | pf[4,32] | +0.0200 | 0.995 | 0.005 | 200 | bias | --FFF | fails at all 1 larger |

**ceiling 1/5**: 2 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | GF | pf[4,32] | +0.0207 | 0.995 | 0.005 | 200 | bias | --FFF | no larger evaluable |
| [8,64] | GF | pf[8,64] | +0.0633 | 0.990 | 0.007 | 200 | bias | ---PF | no larger evaluable |

### ind Q-T3S1

**ceiling 1/80**: 3 failing (10 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GT | zeta_K(0)[8,64] | +0.0225 | 0.865 | 0.017 | 400 | bias | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.25)[8,64] | +0.0215 | 0.845 | 0.018 | 400 | bias | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.5)[8,64] | +0.0206 | 0.833 | 0.019 | 400 | bias | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |

**ceiling 1/40**: 0 failing

**ceiling 1/20**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0237 | 0.995 | 0.005 | 200 | bias | PPFPP | passes at 1/10, 1/5 |

**ceiling 1/10**: 0 failing

**ceiling 1/5**: 0 failing

### ind Q-T3S2

**ceiling 1/80**: 13 failing (33 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_K(0)[4,32] | +0.0470 | 0.378 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.25)[4,32] | +0.0407 | 0.415 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.5)[4,32] | +0.0357 | 0.480 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.75)[4,32] | +0.0319 | 0.570 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(1)[4,32] | +0.0293 | 0.665 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | pf[4,32] | -0.0177 | 0.797 | 0.020 | 400 | coverage | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GF | pf[4,32] | +0.0450 | 0.900 | 0.015 | 400 | bias | FFPPF | passes at 1/20, 1/10 |
| [4,32] | GT | zeta_K(0)[4,32] | +0.0535 | 0.448 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.25)[4,32] | +0.0459 | 0.475 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.5)[4,32] | +0.0395 | 0.522 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.75)[4,32] | +0.0343 | 0.590 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(1)[4,32] | +0.0304 | 0.665 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | pf[4,32] | -0.0231 | 0.770 | 0.021 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |

**ceiling 1/40**: 2 failing (10 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | GF | pf[4,32] | +0.0523 | 0.970 | 0.009 | 400 | bias | FFPPF | passes at 1/20, 1/10 |
| [8,64] | GF | pf[8,64] | +0.0715 | 0.853 | 0.018 | 400 | bias | -FPFF | passes at 1/20 |

**ceiling 1/20**: 0 failing

**ceiling 1/10**: 1 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GF | pf[8,64] | +0.0472 | 0.988 | 0.006 | 400 | bias | -FPFF | fails at all 1 larger |

**ceiling 1/5**: 2 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | GF | pf[4,32] | +0.0237 | 0.983 | 0.007 | 400 | bias | FFPPF | no larger evaluable |
| [8,64] | GF | pf[8,64] | +0.0508 | 0.983 | 0.007 | 400 | bias | -FPFF | no larger evaluable |

### ind Q-T3S3

**ceiling 1/80**: 13 failing (33 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_K(0)[4,32] | +0.0426 | 0.250 | 0.022 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.25)[4,32] | +0.0300 | 0.417 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.5)[4,32] | +0.0213 | 0.600 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(0.75)[4,32] | +0.0164 | 0.713 | 0.023 | 400 | coverage | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | zeta_K(1)[4,32] | +0.0141 | 0.775 | 0.021 | 400 | coverage | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | G0 | pf[4,32] | -0.0282 | 0.537 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GF | pf[4,32] | +0.0356 | 0.880 | 0.016 | 400 | bias | FFFFF | fails at all 4 larger |
| [4,32] | GT | zeta_K(0)[4,32] | +0.0506 | 0.302 | 0.023 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.25)[4,32] | +0.0378 | 0.422 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.5)[4,32] | +0.0291 | 0.530 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(0.75)[4,32] | +0.0237 | 0.630 | 0.024 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | zeta_K(1)[4,32] | +0.0202 | 0.700 | 0.023 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [4,32] | GT | pf[4,32] | -0.0299 | 0.568 | 0.025 | 400 | both | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |

**ceiling 1/40**: 2 failing (10 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | GF | pf[4,32] | +0.0397 | 0.990 | 0.005 | 400 | bias | FFFFF | fails at all 3 larger |
| [8,64] | GF | pf[8,64] | +0.0529 | 0.855 | 0.018 | 400 | bias | -FFFF | fails at all 3 larger |

**ceiling 1/20**: 2 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | GF | pf[4,32] | +0.0389 | 0.990 | 0.005 | 400 | bias | FFFFF | fails at all 2 larger |
| [8,64] | GF | pf[8,64] | +0.0909 | 0.973 | 0.008 | 400 | bias | -FFFF | fails at all 2 larger |

**ceiling 1/10**: 3 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | GF | pf[4,32] | +0.0450 | 0.985 | 0.006 | 400 | bias | FFFFF | fails at all 1 larger |
| [8,64] | GF | pf[8,64] | +0.0841 | 0.985 | 0.006 | 400 | bias | -FFFF | fails at all 1 larger |
| [8,64] | GT | Rratio[8,64] | -0.0214 | 0.943 | 0.012 | 400 | bias | PPPFF | fails at all 1 larger |

**ceiling 1/5**: 3 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | GF | pf[4,32] | +0.0491 | 0.980 | 0.007 | 400 | bias | FFFFF | no larger evaluable |
| [8,64] | GF | pf[8,64] | +0.0966 | 0.978 | 0.007 | 400 | bias | -FFFF | no larger evaluable |
| [8,64] | GT | Rratio[8,64] | -0.0214 | 0.943 | 0.012 | 400 | bias | PPPFF | no larger evaluable |

### cont Q-T3S1

**ceiling 1/80**: 4 failing (10 criterion elements not evaluable)

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GT | zeta_K(0)[8,64] | +0.0245 | 0.870 | 0.034 | 100 | bias | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.25)[8,64] | +0.0232 | 0.850 | 0.036 | 100 | bias | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.5)[8,64] | +0.0220 | 0.840 | 0.037 | 100 | bias | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |
| [8,64] | GT | zeta_K(0.75)[8,64] | +0.0210 | 0.820 | 0.038 | 100 | bias | FPPPP | passes at 1/40, 1/20, 1/10, 1/5 |

**ceiling 1/40**: 2 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | GT | zeta_V(0)[8,64] | -0.0138 | 0.760 | 0.043 | 100 | coverage | -FPFF | passes at 1/20 |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0132 | 0.760 | 0.043 | 100 | coverage | -FPFF | passes at 1/20 |

**ceiling 1/20**: 5 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [4,32] | G0 | zeta_V(0)[4,32] | -0.0084 | 0.760 | 0.043 | 100 | coverage | PPFPP | passes at 1/10, 1/5 |
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0167 | 0.710 | 0.045 | 100 | coverage | -PFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0162 | 0.700 | 0.046 | 100 | coverage | -PFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0144 | 0.730 | 0.044 | 100 | coverage | -PFFF | fails at all 2 larger |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0114 | 0.760 | 0.043 | 100 | coverage | -PFFF | fails at all 2 larger |

**ceiling 1/10**: 7 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0176 | 0.640 | 0.048 | 100 | coverage | -PFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0171 | 0.670 | 0.047 | 100 | coverage | -PFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0151 | 0.650 | 0.048 | 100 | coverage | -PFFF | fails at all 1 larger |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0118 | 0.710 | 0.045 | 100 | coverage | -PFFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0159 | 0.700 | 0.046 | 100 | coverage | -FPFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0154 | 0.730 | 0.044 | 100 | coverage | -FPFF | fails at all 1 larger |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0134 | 0.750 | 0.043 | 100 | coverage | -PPFF | fails at all 1 larger |

**ceiling 1/5**: 7 failing

| range | gen | element | bias | coverage | MC SE | panels | part | pattern | larger ceilings |
|---|---|---|---|---|---|---|---|---|---|
| [8,64] | G0 | zeta_V(0)[8,64] | -0.0185 | 0.650 | 0.048 | 100 | coverage | -PFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.25)[8,64] | -0.0178 | 0.660 | 0.047 | 100 | coverage | -PFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.5)[8,64] | -0.0156 | 0.680 | 0.047 | 100 | coverage | -PFFF | no larger evaluable |
| [8,64] | G0 | zeta_V(0.75)[8,64] | -0.0123 | 0.710 | 0.045 | 100 | coverage | -PFFF | no larger evaluable |
| [8,64] | GT | zeta_V(0)[8,64] | -0.0163 | 0.690 | 0.046 | 100 | coverage | -FPFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.25)[8,64] | -0.0157 | 0.740 | 0.044 | 100 | coverage | -FPFF | no larger evaluable |
| [8,64] | GT | zeta_V(0.5)[8,64] | -0.0137 | 0.760 | 0.043 | 100 | coverage | -PPFF | no larger evaluable |

### 1b. Summary: failing elements and the worst one, per cell and ceiling

Worst = largest |bias| / 0.02 + max(0, 0.85 - coverage - 2 MC SE) / 0.05 (the score `analyse` uses for rho.md). n/e marks a ceiling at which some criterion estimator has no evaluable element.

| days | cell | 1/80 | 1/40 | 1/20 | 1/10 | 1/5 |
|---|---|---|---|---|---|---|
| ind | Q-T1S1 | 1: GF pf[8,64] b +0.028 c 0.95+-0.01 (bias) | 1: GF pf[8,64] b +0.029 c 0.99+-0.01 (bias) | 1: GF pf[8,64] b +0.034 c 0.99+-0.01 (bias) | 1: GF pf[8,64] b +0.037 c 0.99+-0.00 (bias) | 1: GF pf[8,64] b +0.033 c 1.00+-0.00 (bias) |
| ind | Q-T1S2 | 1: GF pf[8,64] b +0.040 c 0.97+-0.01 (bias) | 1: GF pf[8,64] b +0.029 c 0.99+-0.00 (bias) | 0 | 9: GT zeta_V(0)[8,64] b -0.007 c 0.74+-0.02 (coverage) | 17: GT zeta_V(0)[8,64] b -0.008 c 0.58+-0.02 (coverage) |
| ind | Q-T1S3 | n/e, 13: G0 zeta_K(0)[4,32] b +0.086 c 0.00+-0.00 (both) | 17: G0 zeta_K(0)[8,64] b +0.048 c 0.56+-0.02 (both) | 20: GT Rratio[4,32] b -0.016 c 0.47+-0.02 (coverage) | 23: G0 zeta_V(0)[4,32] b -0.016 c 0.19+-0.02 (coverage) | 25: GT zeta_V(0)[8,64] b -0.029 c 0.03+-0.01 (both) |
| ind | Q-T2S1 | 10: G0 zeta_K(0.25)[8,64] b +0.038 c 0.51+-0.02 (both) | 12: G0 zeta_V(0)[8,64] b -0.024 c 0.75+-0.02 (both) | 17: G0 zeta_V(0.25)[8,64] b -0.024 c 0.61+-0.02 (both) | 18: G0 zeta_V(0)[8,64] b -0.026 c 0.39+-0.02 (both) | 18: G0 zeta_V(0)[8,64] b -0.028 c 0.25+-0.02 (both) |
| ind | Q-T2S2 | n/e, 13: G0 zeta_K(0)[4,32] b +0.070 c 0.10+-0.01 (both) | 11: G0 zeta_K(0)[8,64] b +0.043 c 0.79+-0.02 (both) | 12: G0 zeta_V(0)[8,64] b -0.026 c 0.77+-0.02 (both) | 14: G0 zeta_V(0)[8,64] b -0.027 c 0.59+-0.02 (both) | 14: G0 zeta_V(0)[8,64] b -0.028 c 0.46+-0.02 (both) |
| ind | Q-T2S3 | n/e, 0 | n/e, 0 | n/e, 13: G0 zeta_K(0)[4,32] b +0.050 c 0.35+-0.02 (both) | 1: GF pf[4,32] b +0.020 c 0.99+-0.00 (bias) | 2: GF pf[8,64] b +0.063 c 0.99+-0.01 (bias) |
| ind | Q-T3S1 | 3: GT zeta_K(0)[8,64] b +0.022 c 0.86+-0.02 (bias) | 0 | 1: GF pf[8,64] b +0.024 c 0.99+-0.00 (bias) | 0 | 0 |
| ind | Q-T3S2 | n/e, 13: G0 zeta_K(0)[4,32] b +0.047 c 0.38+-0.02 (both) | 2: GF pf[8,64] b +0.072 c 0.85+-0.02 (bias) | 0 | 1: GF pf[8,64] b +0.047 c 0.99+-0.01 (bias) | 2: GF pf[8,64] b +0.051 c 0.98+-0.01 (bias) |
| ind | Q-T3S3 | n/e, 13: G0 zeta_K(0)[4,32] b +0.043 c 0.25+-0.02 (both) | 2: GF pf[8,64] b +0.053 c 0.85+-0.02 (bias) | 2: GF pf[8,64] b +0.091 c 0.97+-0.01 (bias) | 3: GF pf[8,64] b +0.084 c 0.98+-0.01 (bias) | 3: GF pf[8,64] b +0.097 c 0.98+-0.01 (bias) |
| cont | Q-T3S1 | 4: GT zeta_K(0)[8,64] b +0.024 c 0.87+-0.03 (bias) | 2: GT zeta_V(0)[8,64] b -0.014 c 0.76+-0.04 (coverage) | 5: G0 zeta_V(0.25)[8,64] b -0.016 c 0.70+-0.05 (coverage) | 7: G0 zeta_V(0)[8,64] b -0.018 c 0.64+-0.05 (coverage) | 7: G0 zeta_V(0)[8,64] b -0.018 c 0.65+-0.05 (coverage) |

## 2. Stocks contributing, across panels

Median (10th percentile) across panels of the number of stocks eligible. Range = eligible at the range's smallest horizon (T = 4 or 8); top = eligible at its largest horizon (T = 32 or 64). Day lengths are drawn per generator, so G0, GT and GF differ slightly; the continuous stress cell uses the same seeds, hence the same day lengths, as panels 0 to 99 of ind Q-T3S1.

### 2a. Stocks eligible at the range (smallest T)

| days | cell | gen | range | 1/80 | 1/40 | 1/20 | 1/10 | 1/5 |
|---|---|---|---|---|---|---|---|---|
| ind | Q-T1S1 | G0 | [4,32] | 24 (18) | 45 (37) | 75 (68) | 96 (88) | 110 (103) |
| ind | Q-T1S1 | G0 | [8,64] | 12 (8) | 24 (18) | 45 (37) | 75 (68) | 96 (88) |
| ind | Q-T1S1 | GT | [4,32] | 23 (18) | 44 (37) | 75 (67) | 96 (87) | 110 (102) |
| ind | Q-T1S1 | GT | [8,64] | 11 (7) | 23 (18) | 44 (37) | 75 (67) | 96 (87) |
| ind | Q-T1S1 | GF | [4,32] | 23 (17) | 44 (37) | 74 (65.9) | 95.5 (88) | 108 (101) |
| ind | Q-T1S1 | GF | [8,64] | 12 (7) | 23 (17) | 44 (37) | 74 (65.9) | 95.5 (88) |
| ind | Q-T1S2 | G0 | [4,32] | 103 (94) | 141 (132) | 171.5 (164) | 188 (182) | 200 (195) |
| ind | Q-T1S2 | G0 | [8,64] | 61 (53) | 103 (94) | 141 (132) | 171.5 (164) | 188 (182) |
| ind | Q-T1S2 | GT | [4,32] | 103.5 (95) | 140 (131) | 171 (164) | 189 (181) | 200 (194.9) |
| ind | Q-T1S2 | GT | [8,64] | 61 (52) | 103.5 (95) | 140 (131) | 171 (164) | 189 (181) |
| ind | Q-T1S2 | GF | [4,32] | 104 (94) | 142 (132) | 172 (165) | 188 (182) | 199 (195) |
| ind | Q-T1S2 | GF | [8,64] | 61 (53) | 104 (94) | 142 (132) | 172 (165) | 188 (182) |
| ind | Q-T1S3 | G0 | [4,32] | 73 (63) | 114 (104) | 146 (137) | 174 (166) | 195 (187) |
| ind | Q-T1S3 | G0 | [8,64] | 26 (20) | 73 (63) | 114 (104) | 146 (137) | 174 (166) |
| ind | Q-T1S3 | GT | [4,32] | 72 (63) | 115 (104) | 147 (137) | 175 (166) | 196 (188) |
| ind | Q-T1S3 | GT | [8,64] | 25 (19) | 72 (63) | 115 (104) | 147 (137) | 175 (166) |
| ind | Q-T1S3 | GF | [4,32] | 71.5 (63) | 114 (104) | 146 (138) | 175 (166.9) | 196 (187) |
| ind | Q-T1S3 | GF | [8,64] | 25 (20) | 71.5 (63) | 114 (104) | 146 (138) | 175 (166.9) |
| ind | Q-T2S1 | G0 | [4,32] | 39 (33) | 52 (47) | 71 (67) | 71 (67) | 71 (67) |
| ind | Q-T2S1 | G0 | [8,64] | 23 (18) | 39 (33) | 52 (47) | 71 (67) | 71 (67) |
| ind | Q-T2S1 | GT | [4,32] | 39 (33) | 52 (47) | 71 (67) | 71 (67) | 71 (67) |
| ind | Q-T2S1 | GT | [8,64] | 23 (18) | 39 (33) | 52 (47) | 71 (67) | 71 (67) |
| ind | Q-T2S1 | GF | [4,32] | 38 (33) | 53 (47) | 71 (67.9) | 71 (67.9) | 71 (67.9) |
| ind | Q-T2S1 | GF | [8,64] | 23 (18) | 38 (33) | 53 (47) | 71 (67.9) | 71 (67.9) |
| ind | Q-T2S2 | G0 | [4,32] | 28 (24) | 33 (29) | 41 (37) | 43 (40) | 48 (46) |
| ind | Q-T2S2 | G0 | [8,64] | 14 (10) | 28 (24) | 33 (29) | 41 (37) | 43 (40) |
| ind | Q-T2S2 | GT | [4,32] | 28 (24) | 33 (29) | 40 (37) | 43 (40) | 48 (45) |
| ind | Q-T2S2 | GT | [8,64] | 14 (10) | 28 (24) | 33 (29) | 40 (37) | 43 (40) |
| ind | Q-T2S2 | GF | [4,32] | 28 (23) | 33 (29) | 40 (36) | 42 (39) | 48 (45) |
| ind | Q-T2S2 | GF | [8,64] | 14 (10) | 28 (23) | 33 (29) | 40 (36) | 42 (39) |
| ind | Q-T2S3 | G0 | [4,32] | 12 (10) | 15 (13) | 15 (13) | 15 (13) | 17 (17) |
| ind | Q-T2S3 | G0 | [8,64] | 0 (0) | 12 (10) | 15 (13) | 15 (13) | 15 (13) |
| ind | Q-T2S3 | GT | [4,32] | 12 (10) | 15 (13) | 15 (13) | 15 (13) | 17 (17) |
| ind | Q-T2S3 | GT | [8,64] | 0 (0) | 12 (10) | 15 (13) | 15 (13) | 15 (13) |
| ind | Q-T2S3 | GF | [4,32] | 12 (10) | 15 (13) | 15 (13) | 15 (13) | 17 (17) |
| ind | Q-T2S3 | GF | [8,64] | 0 (0) | 12 (10) | 15 (13) | 15 (13) | 15 (13) |
| ind | Q-T3S1 | G0 | [4,32] | 41 (37) | 45 (42) | 48 (45) | 48 (45) | 50 (50) |
| ind | Q-T3S1 | G0 | [8,64] | 23 (19) | 41 (37) | 45 (42) | 48 (45) | 48 (45) |
| ind | Q-T3S1 | GT | [4,32] | 41 (37.9) | 46 (43) | 48 (46) | 48 (46) | 50 (50) |
| ind | Q-T3S1 | GT | [8,64] | 23 (18.9) | 41 (37.9) | 46 (43) | 48 (46) | 48 (46) |
| ind | Q-T3S1 | GF | [4,32] | 41 (37.9) | 45.5 (42) | 48 (46) | 48 (46) | 50 (50) |
| ind | Q-T3S1 | GF | [8,64] | 23 (18) | 41 (37.9) | 45.5 (42) | 48 (46) | 48 (46) |
| ind | Q-T3S2 | G0 | [4,32] | 17 (14) | 27 (27) | 27 (27) | 27 (27) | 27 (27) |
| ind | Q-T3S2 | G0 | [8,64] | 5 (2) | 17 (14) | 27 (27) | 27 (27) | 27 (27) |
| ind | Q-T3S2 | GT | [4,32] | 17 (14) | 27 (27) | 27 (27) | 27 (27) | 27 (27) |
| ind | Q-T3S2 | GT | [8,64] | 5 (3) | 17 (14) | 27 (27) | 27 (27) | 27 (27) |
| ind | Q-T3S2 | GF | [4,32] | 17 (14) | 27 (27) | 27 (27) | 27 (27) | 27 (27) |
| ind | Q-T3S2 | GF | [8,64] | 5 (2) | 17 (14) | 27 (27) | 27 (27) | 27 (27) |
| ind | Q-T3S3 | G0 | [4,32] | 14 (11) | 21 (19) | 21 (19) | 21 (19) | 23 (23) |
| ind | Q-T3S3 | G0 | [8,64] | 7 (4) | 14 (11) | 21 (19) | 21 (19) | 21 (19) |
| ind | Q-T3S3 | GT | [4,32] | 14 (11) | 21 (19) | 21 (19) | 21 (19) | 23 (23) |
| ind | Q-T3S3 | GT | [8,64] | 7 (4) | 14 (11) | 21 (19) | 21 (19) | 21 (19) |
| ind | Q-T3S3 | GF | [4,32] | 14 (11) | 21 (19) | 21 (19) | 21 (19) | 23 (23) |
| ind | Q-T3S3 | GF | [8,64] | 7 (4) | 14 (11) | 21 (19) | 21 (19) | 21 (19) |
| cont | Q-T3S1 | G0 | [4,32] | 41 (38) | 45 (42.9) | 48 (46) | 48 (46) | 50 (50) |
| cont | Q-T3S1 | G0 | [8,64] | 23 (19) | 41 (38) | 45 (42.9) | 48 (46) | 48 (46) |
| cont | Q-T3S1 | GT | [4,32] | 41 (38) | 46 (43) | 48 (46) | 48 (46) | 50 (50) |
| cont | Q-T3S1 | GT | [8,64] | 23 (19.9) | 41 (38) | 46 (43) | 48 (46) | 48 (46) |

### 2b. Stocks eligible at the top (largest T)

| days | cell | gen | range | 1/80 | 1/40 | 1/20 | 1/10 | 1/5 |
|---|---|---|---|---|---|---|---|---|
| ind | Q-T1S1 | G0 | [4,32] | 4 (2) | 9 (5) | 12 (8) | 24 (18) | 45 (37) |
| ind | Q-T1S1 | G0 | [8,64] | 4 (2) | 4 (2) | 9 (5) | 12 (8) | 24 (18) |
| ind | Q-T1S1 | GT | [4,32] | 5 (2) | 8.5 (5) | 11 (7) | 23 (18) | 44 (37) |
| ind | Q-T1S1 | GT | [8,64] | 5 (2) | 5 (2) | 8.5 (5) | 11 (7) | 23 (18) |
| ind | Q-T1S1 | GF | [4,32] | 5 (2) | 9 (6) | 12 (7) | 23 (17) | 44 (37) |
| ind | Q-T1S1 | GF | [8,64] | 5 (2) | 5 (2) | 9 (6) | 12 (7) | 23 (17) |
| ind | Q-T1S2 | G0 | [4,32] | 5 (2) | 23 (18) | 61 (53) | 103 (94) | 141 (132) |
| ind | Q-T1S2 | G0 | [8,64] | 2 (1) | 5 (2) | 23 (18) | 61 (53) | 103 (94) |
| ind | Q-T1S2 | GT | [4,32] | 5 (2) | 23 (18) | 61 (52) | 103.5 (95) | 140 (131) |
| ind | Q-T1S2 | GT | [8,64] | 2 (1) | 5 (2) | 23 (18) | 61 (52) | 103.5 (95) |
| ind | Q-T1S2 | GF | [4,32] | 4.5 (2) | 23 (18.9) | 61 (53) | 104 (94) | 142 (132) |
| ind | Q-T1S2 | GF | [8,64] | 2 (1) | 4.5 (2) | 23 (18.9) | 61 (53) | 104 (94) |
| ind | Q-T1S3 | G0 | [4,32] | 0 (0) | 7 (4) | 26 (20) | 73 (63) | 114 (104) |
| ind | Q-T1S3 | G0 | [8,64] | 0 (0) | 0 (0) | 7 (4) | 26 (20) | 73 (63) |
| ind | Q-T1S3 | GT | [4,32] | 0 (0) | 7 (4) | 25 (19) | 72 (63) | 115 (104) |
| ind | Q-T1S3 | GT | [8,64] | 0 (0) | 0 (0) | 7 (4) | 25 (19) | 72 (63) |
| ind | Q-T1S3 | GF | [4,32] | 0 (0) | 6 (4) | 25 (20) | 71.5 (63) | 114 (104) |
| ind | Q-T1S3 | GF | [8,64] | 0 (0) | 0 (0) | 6 (4) | 25 (20) | 71.5 (63) |
| ind | Q-T2S1 | G0 | [4,32] | 7 (4) | 16 (12) | 23 (18) | 39 (33) | 52 (47) |
| ind | Q-T2S1 | G0 | [8,64] | 0 (0) | 7 (4) | 16 (12) | 23 (18) | 39 (33) |
| ind | Q-T2S1 | GT | [4,32] | 7 (4) | 16 (12) | 23 (18) | 39 (33) | 52 (47) |
| ind | Q-T2S1 | GT | [8,64] | 0 (0) | 7 (4) | 16 (12) | 23 (18) | 39 (33) |
| ind | Q-T2S1 | GF | [4,32] | 7 (3.9) | 16 (12) | 23 (18) | 38 (33) | 53 (47) |
| ind | Q-T2S1 | GF | [8,64] | 0 (0) | 7 (3.9) | 16 (12) | 23 (18) | 38 (33) |
| ind | Q-T2S2 | G0 | [4,32] | 0 (0) | 4 (2) | 14 (10) | 28 (24) | 33 (29) |
| ind | Q-T2S2 | G0 | [8,64] | 0 (0) | 0 (0) | 4 (2) | 14 (10) | 28 (24) |
| ind | Q-T2S2 | GT | [4,32] | 0 (0) | 5 (2) | 14 (10) | 28 (24) | 33 (29) |
| ind | Q-T2S2 | GT | [8,64] | 0 (0) | 0 (0) | 5 (2) | 14 (10) | 28 (24) |
| ind | Q-T2S2 | GF | [4,32] | 0 (0) | 4.5 (2) | 14 (10) | 28 (23) | 33 (29) |
| ind | Q-T2S2 | GF | [8,64] | 0 (0) | 0 (0) | 4.5 (2) | 14 (10) | 28 (23) |
| ind | Q-T2S3 | G0 | [4,32] | 0 (0) | 0 (0) | 0 (0) | 12 (10) | 15 (13) |
| ind | Q-T2S3 | G0 | [8,64] | 0 (0) | 0 (0) | 0 (0) | 0 (0) | 12 (10) |
| ind | Q-T2S3 | GT | [4,32] | 0 (0) | 0 (0) | 0 (0) | 12 (10) | 15 (13) |
| ind | Q-T2S3 | GT | [8,64] | 0 (0) | 0 (0) | 0 (0) | 0 (0) | 12 (10) |
| ind | Q-T2S3 | GF | [4,32] | 0 (0) | 0 (0) | 0 (0) | 12 (10) | 15 (13) |
| ind | Q-T2S3 | GF | [8,64] | 0 (0) | 0 (0) | 0 (0) | 0 (0) | 12 (10) |
| ind | Q-T3S1 | G0 | [4,32] | 7 (4) | 13 (10) | 23 (19) | 41 (37) | 45 (42) |
| ind | Q-T3S1 | G0 | [8,64] | 0 (0) | 7 (4) | 13 (10) | 23 (19) | 41 (37) |
| ind | Q-T3S1 | GT | [4,32] | 7 (3.9) | 14 (10) | 23 (18.9) | 41 (37.9) | 46 (43) |
| ind | Q-T3S1 | GT | [8,64] | 0 (0) | 7 (3.9) | 14 (10) | 23 (18.9) | 41 (37.9) |
| ind | Q-T3S1 | GF | [4,32] | 7 (4) | 14 (9) | 23 (18) | 41 (37.9) | 45.5 (42) |
| ind | Q-T3S1 | GF | [8,64] | 0 (0) | 7 (4) | 14 (9) | 23 (18) | 41 (37.9) |
| ind | Q-T3S2 | G0 | [4,32] | 0 (0) | 2 (1) | 5 (2) | 17 (14) | 27 (27) |
| ind | Q-T3S2 | G0 | [8,64] | 0 (0) | 0 (0) | 2 (1) | 5 (2) | 17 (14) |
| ind | Q-T3S2 | GT | [4,32] | 0 (0) | 2 (1) | 5 (3) | 17 (14) | 27 (27) |
| ind | Q-T3S2 | GT | [8,64] | 0 (0) | 0 (0) | 2 (1) | 5 (3) | 17 (14) |
| ind | Q-T3S2 | GF | [4,32] | 0 (0) | 2 (1) | 5 (2) | 17 (14) | 27 (27) |
| ind | Q-T3S2 | GF | [8,64] | 0 (0) | 0 (0) | 2 (1) | 5 (2) | 17 (14) |
| ind | Q-T3S3 | G0 | [4,32] | 0 (0) | 2 (1) | 7 (4) | 14 (11) | 21 (19) |
| ind | Q-T3S3 | G0 | [8,64] | 0 (0) | 0 (0) | 2 (1) | 7 (4) | 14 (11) |
| ind | Q-T3S3 | GT | [4,32] | 0 (0) | 2 (0) | 7 (4) | 14 (11) | 21 (19) |
| ind | Q-T3S3 | GT | [8,64] | 0 (0) | 0 (0) | 2 (0) | 7 (4) | 14 (11) |
| ind | Q-T3S3 | GF | [4,32] | 0 (0) | 2 (1) | 7 (4) | 14 (11) | 21 (19) |
| ind | Q-T3S3 | GF | [8,64] | 0 (0) | 0 (0) | 2 (1) | 7 (4) | 14 (11) |
| cont | Q-T3S1 | G0 | [4,32] | 7 (4) | 13.5 (9.9) | 23 (19) | 41 (38) | 45 (42.9) |
| cont | Q-T3S1 | G0 | [8,64] | 0 (0) | 7 (4) | 13.5 (9.9) | 23 (19) | 41 (38) |
| cont | Q-T3S1 | GT | [4,32] | 7 (4) | 14 (10.9) | 23 (19.9) | 41 (38) | 46 (43) |
| cont | Q-T3S1 | GT | [8,64] | 0 (0) | 7 (4) | 14 (10.9) | 23 (19.9) | 41 (38) |

## 3. Do failures persist at larger ceilings?

Per element that fails at one or more evaluable ceilings, its pattern over the evaluable ceilings in order: fails everywhere (F only), fails at small ceilings only (F...FP...P), fails at large ceilings only (P...PF...F), mixed (anything else, e.g. FPF or PFP). Counted per element (each element once) and per failure (each failing (element, ceiling) pair, by its element's class). Independent-days cells only; the stress cell is in section 6.

| family | elements failing | small only | everywhere | large only | mixed | of which 'everywhere' with one evaluable ceiling |
|---|---|---|---|---|---|---|
| zeta_V(a) | 66 | 0 | 20 | 46 | 0 | 0 |
| zeta_K(a) | 83 | 83 | 0 | 0 | 0 | 0 |
| P-F slope | 24 | 15 | 5 | 1 | 3 | 0 |
| R ratio | 8 | 0 | 0 | 8 | 0 | 0 |
| all | 181 | 98 | 25 | 55 | 3 | 0 |

Per failing (element, ceiling) pair:

| family | failures | small only | everywhere | large only | mixed |
|---|---|---|---|---|---|
| zeta_V(a) | 181 | 0 | 67 | 114 | 0 |
| zeta_K(a) | 83 | 83 | 0 | 0 | 0 |
| P-F slope | 47 | 18 | 21 | 1 | 7 |
| R ratio | 20 | 0 | 0 | 20 | 0 |

Elements by generator and class:

| gen | family | small only | everywhere | large only | mixed |
|---|---|---|---|---|---|
| G0 | zeta_V(a) | 0 | 9 | 24 | 0 |
| G0 | zeta_K(a) | 40 | 0 | 0 | 0 |
| G0 | P-F slope | 5 | 0 | 0 | 0 |
| G0 | R ratio | 0 | 0 | 4 | 0 |
| GT | zeta_V(a) | 0 | 11 | 22 | 0 |
| GT | zeta_K(a) | 43 | 0 | 0 | 0 |
| GT | P-F slope | 6 | 0 | 0 | 0 |
| GT | R ratio | 0 | 0 | 4 | 0 |
| GF | P-F slope | 4 | 5 | 1 | 3 |

## 4. Coverage by median stocks contributing to the range

Evaluable criterion elements of the independent-days cells, pooled over cells, generators, ranges, a and ceilings, binned by the element's median (across panels) stock count. Coverage is the share of (element, panel) pairs whose interval holds the truth; MC SE is clustered on panels (elements from the same panel are not independent). `elem` = evaluable (element, ceiling) pairs in the bin; `fail` = how many of them fail the rule; `|bias|` = their mean absolute bias. Medians below 1 (fewer than half the panels have any eligible stock) are in bin 0.

| family | 0 | 1-2 | 3-5 | 6-10 | 11-20 | 21-50 | 51+ |
|---|---|---|---|---|---|---|---|
| zeta_V(a) | - | - | - | - | 0.967 +- 0.004 (elem 40, fail 0, abs bias 0.003) | 0.899 +- 0.003 (elem 360, fail 45, abs bias 0.006) | 0.724 +- 0.005 (elem 290, fail 136, abs bias 0.009) |
| zeta_K(a) | - | - | - | - | 0.816 +- 0.004 (elem 100, fail 30, abs bias 0.014) | 0.945 +- 0.001 (elem 400, fail 33, abs bias 0.005) | 0.941 +- 0.002 (elem 310, fail 20, abs bias 0.005) |
| P-F slope | - | - | - | - | 0.890 +- 0.004 (elem 30, fail 15, abs bias 0.019) | 0.977 +- 0.001 (elem 120, fail 19, abs bias 0.011) | 0.961 +- 0.001 (elem 93, fail 13, abs bias 0.009) |
| R ratio | - | 0.862 +- 0.009 (elem 4, fail 0, abs bias 0.006) | 0.958 +- 0.005 (elem 6, fail 0, abs bias 0.005) | 0.963 +- 0.004 (elem 8, fail 0, abs bias 0.007) | 0.972 +- 0.002 (elem 34, fail 0, abs bias 0.006) | 0.950 +- 0.002 (elem 74, fail 3, abs bias 0.007) | 0.783 +- 0.006 (elem 48, fail 17, abs bias 0.009) |

### 4b. The same, binned by median stocks at the range's top

| family | 0 | 1-2 | 3-5 | 6-10 | 11-20 | 21-50 | 51+ |
|---|---|---|---|---|---|---|---|
| zeta_V(a) | - | 0.858 +- 0.007 (elem 50, fail 0, abs bias 0.003) | 0.927 +- 0.004 (elem 90, fail 6, abs bias 0.005) | 0.906 +- 0.004 (elem 100, fail 14, abs bias 0.007) | 0.905 +- 0.003 (elem 150, fail 28, abs bias 0.007) | 0.790 +- 0.004 (elem 210, fail 74, abs bias 0.009) | 0.628 +- 0.008 (elem 90, fail 59, abs bias 0.009) |
| zeta_K(a) | 0.635 +- 0.004 (elem 120, fail 83, abs bias 0.031) | 0.957 +- 0.003 (elem 50, fail 0, abs bias 0.003) | 0.979 +- 0.002 (elem 90, fail 0, abs bias 0.001) | 0.985 +- 0.001 (elem 100, fail 0, abs bias 0.002) | 0.984 +- 0.001 (elem 150, fail 0, abs bias 0.001) | 0.985 +- 0.001 (elem 210, fail 0, abs bias 0.002) | 0.985 +- 0.002 (elem 90, fail 0, abs bias 0.002) |
| P-F slope | 0.827 +- 0.004 (elem 36, fail 20, abs bias 0.023) | 0.969 +- 0.002 (elem 15, fail 4, abs bias 0.018) | 0.984 +- 0.002 (elem 27, fail 5, abs bias 0.009) | 0.987 +- 0.001 (elem 30, fail 5, abs bias 0.011) | 0.988 +- 0.001 (elem 45, fail 8, abs bias 0.010) | 0.988 +- 0.001 (elem 63, fail 4, abs bias 0.006) | 0.985 +- 0.002 (elem 27, fail 1, abs bias 0.004) |
| R ratio | - | 0.862 +- 0.009 (elem 4, fail 0, abs bias 0.006) | 0.958 +- 0.005 (elem 6, fail 0, abs bias 0.005) | 0.963 +- 0.004 (elem 8, fail 0, abs bias 0.007) | 0.972 +- 0.002 (elem 34, fail 0, abs bias 0.006) | 0.950 +- 0.002 (elem 74, fail 3, abs bias 0.007) | 0.783 +- 0.006 (elem 48, fail 17, abs bias 0.009) |

## 5. Which interval was widest, by stock-count bin

Share of (element, panel) pairs with a defined interval in which each of the three intervals was the widest (1 = day-and-stock bootstrap, 2 = day jackknife with t, 3 = whole-stock bootstrap), same bins (median stocks contributing to the range) and same pool as section 4. n = pairs with an interval.

| family | bin | n | day-stock bootstrap | day jackknife t | whole-stock bootstrap |
|---|---|---|---|---|---|
| zeta_V(a) | 11-20 | 13985 | 0.962 | 0.038 | 0.000 |
| zeta_V(a) | 21-50 | 134550 | 0.944 | 0.055 | 0.000 |
| zeta_V(a) | 51+ | 105570 | 0.970 | 0.030 | 0.000 |
| zeta_K(a) | 11-20 | 36675 | 0.864 | 0.136 | 0.001 |
| zeta_K(a) | 21-50 | 151890 | 0.962 | 0.036 | 0.001 |
| zeta_K(a) | 51+ | 113965 | 0.981 | 0.018 | 0.001 |
| P-F slope | 11-20 | 9975 | 0.849 | 0.144 | 0.008 |
| P-F slope | 21-50 | 41141 | 0.939 | 0.057 | 0.004 |
| P-F slope | 51+ | 28990 | 0.973 | 0.024 | 0.002 |
| R ratio | 1-2 | 1469 | 0.640 | 0.356 | 0.004 |
| R ratio | 3-5 | 2384 | 0.914 | 0.083 | 0.003 |
| R ratio | 6-10 | 2800 | 0.969 | 0.029 | 0.003 |
| R ratio | 11-20 | 12800 | 0.986 | 0.014 | 0.000 |
| R ratio | 21-50 | 28000 | 0.989 | 0.011 | 0.000 |
| R ratio | 51+ | 18000 | 0.993 | 0.007 | 0.000 |

## 6. Stress cell: cont Q-T3S1 next to ind Q-T3S1

G0 and GT only (the stress test does not run GF). ind has 400 panels per generator, cont 100. Panels 0 to 99 of ind use the same seeds, hence the same day lengths, as the 100 cont panels (the simulated paths differ), so `ind cov 0-99` is the coverage on those matched panels.

### 6a. Failing elements per ceiling

| ceiling | ind fails | ind worst | cont fails | cont worst | fail in both | ind only | cont only |
|---|---|---|---|---|---|---|---|
| 1/80 | 3 | GT zeta_K(0)[8,64] b +0.022 c 0.86+-0.02 | 4 | GT zeta_K(0)[8,64] b +0.024 c 0.87+-0.03 | 3 | 0 | 1 |
| 1/40 | 0 | - | 2 | GT zeta_V(0)[8,64] b -0.014 c 0.76+-0.04 | 0 | 0 | 2 |
| 1/20 | 0 | - | 5 | G0 zeta_V(0.25)[8,64] b -0.016 c 0.70+-0.05 | 0 | 0 | 5 |
| 1/10 | 0 | - | 7 | G0 zeta_V(0)[8,64] b -0.018 c 0.64+-0.05 | 0 | 0 | 7 |
| 1/5 | 0 | - | 7 | G0 zeta_V(0)[8,64] b -0.018 c 0.65+-0.05 | 0 | 0 | 7 |

### 6b. Every criterion element, side by side

Bias and coverage (MC SE) per element and ceiling. Status: pass, or fail with the part that failed; - = not evaluable.

| gen | element | ceiling | ind bias | ind cov | ind status | ind cov 0-99 | cont bias | cont cov | cont status | stocks range / top (med) |
|---|---|---|---|---|---|---|---|---|---|---|
| G0 | zeta_V(0)[4,32] | 1/80 | -0.0029 | 0.950 (0.011) | pass | 0.950 | -0.0066 | 0.850 (0.036) | pass | 41 / 7 |
| G0 | zeta_V(0)[4,32] | 1/40 | -0.0033 | 0.953 (0.011) | pass | 0.950 | -0.0081 | 0.790 (0.041) | pass | 45 / 13.5 |
| G0 | zeta_V(0)[4,32] | 1/20 | -0.0039 | 0.938 (0.012) | pass | 0.960 | -0.0084 | 0.760 (0.043) | fail (coverage) | 48 / 23 |
| G0 | zeta_V(0)[4,32] | 1/10 | -0.0042 | 0.932 (0.013) | pass | 0.940 | -0.0090 | 0.800 (0.040) | pass | 48 / 41 |
| G0 | zeta_V(0)[4,32] | 1/5 | -0.0045 | 0.917 (0.014) | pass | 0.930 | -0.0090 | 0.800 (0.040) | pass | 50 / 45 |
| G0 | zeta_V(0)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| G0 | zeta_V(0)[8,64] | 1/40 | -0.0060 | 0.950 (0.011) | pass | 0.960 | -0.0138 | 0.840 (0.037) | pass | 41 / 7 |
| G0 | zeta_V(0)[8,64] | 1/20 | -0.0071 | 0.920 (0.014) | pass | 0.940 | -0.0167 | 0.710 (0.045) | fail (coverage) | 45 / 13.5 |
| G0 | zeta_V(0)[8,64] | 1/10 | -0.0081 | 0.890 (0.016) | pass | 0.890 | -0.0176 | 0.640 (0.048) | fail (coverage) | 48 / 23 |
| G0 | zeta_V(0)[8,64] | 1/5 | -0.0089 | 0.840 (0.018) | pass | 0.820 | -0.0185 | 0.650 (0.048) | fail (coverage) | 48 / 41 |
| G0 | zeta_V(0.25)[4,32] | 1/80 | -0.0028 | 0.953 (0.011) | pass | 0.950 | -0.0065 | 0.880 (0.032) | pass | 41 / 7 |
| G0 | zeta_V(0.25)[4,32] | 1/40 | -0.0032 | 0.955 (0.010) | pass | 0.940 | -0.0080 | 0.770 (0.042) | pass | 45 / 13.5 |
| G0 | zeta_V(0.25)[4,32] | 1/20 | -0.0037 | 0.945 (0.011) | pass | 0.950 | -0.0082 | 0.770 (0.042) | pass | 48 / 23 |
| G0 | zeta_V(0.25)[4,32] | 1/10 | -0.0040 | 0.940 (0.012) | pass | 0.970 | -0.0086 | 0.770 (0.042) | pass | 48 / 41 |
| G0 | zeta_V(0.25)[4,32] | 1/5 | -0.0044 | 0.917 (0.014) | pass | 0.930 | -0.0086 | 0.780 (0.041) | pass | 50 / 45 |
| G0 | zeta_V(0.25)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| G0 | zeta_V(0.25)[8,64] | 1/40 | -0.0058 | 0.938 (0.012) | pass | 0.950 | -0.0134 | 0.820 (0.038) | pass | 41 / 7 |
| G0 | zeta_V(0.25)[8,64] | 1/20 | -0.0068 | 0.922 (0.013) | pass | 0.920 | -0.0162 | 0.700 (0.046) | fail (coverage) | 45 / 13.5 |
| G0 | zeta_V(0.25)[8,64] | 1/10 | -0.0079 | 0.880 (0.016) | pass | 0.890 | -0.0171 | 0.670 (0.047) | fail (coverage) | 48 / 23 |
| G0 | zeta_V(0.25)[8,64] | 1/5 | -0.0087 | 0.848 (0.018) | pass | 0.840 | -0.0178 | 0.660 (0.047) | fail (coverage) | 48 / 41 |
| G0 | zeta_V(0.5)[4,32] | 1/80 | -0.0025 | 0.960 (0.010) | pass | 0.950 | -0.0056 | 0.880 (0.032) | pass | 41 / 7 |
| G0 | zeta_V(0.5)[4,32] | 1/40 | -0.0028 | 0.963 (0.009) | pass | 0.950 | -0.0071 | 0.790 (0.041) | pass | 45 / 13.5 |
| G0 | zeta_V(0.5)[4,32] | 1/20 | -0.0033 | 0.945 (0.011) | pass | 0.930 | -0.0072 | 0.810 (0.039) | pass | 48 / 23 |
| G0 | zeta_V(0.5)[4,32] | 1/10 | -0.0037 | 0.953 (0.011) | pass | 1.000 | -0.0075 | 0.810 (0.039) | pass | 48 / 41 |
| G0 | zeta_V(0.5)[4,32] | 1/5 | -0.0040 | 0.935 (0.012) | pass | 0.970 | -0.0074 | 0.810 (0.039) | pass | 50 / 45 |
| G0 | zeta_V(0.5)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| G0 | zeta_V(0.5)[8,64] | 1/40 | -0.0051 | 0.948 (0.011) | pass | 0.950 | -0.0118 | 0.850 (0.036) | pass | 41 / 7 |
| G0 | zeta_V(0.5)[8,64] | 1/20 | -0.0061 | 0.930 (0.013) | pass | 0.930 | -0.0144 | 0.730 (0.044) | fail (coverage) | 45 / 13.5 |
| G0 | zeta_V(0.5)[8,64] | 1/10 | -0.0072 | 0.900 (0.015) | pass | 0.920 | -0.0151 | 0.650 (0.048) | fail (coverage) | 48 / 23 |
| G0 | zeta_V(0.5)[8,64] | 1/5 | -0.0080 | 0.882 (0.016) | pass | 0.900 | -0.0156 | 0.680 (0.047) | fail (coverage) | 48 / 41 |
| G0 | zeta_V(0.75)[4,32] | 1/80 | -0.0020 | 0.968 (0.009) | pass | 0.980 | -0.0042 | 0.910 (0.029) | pass | 41 / 7 |
| G0 | zeta_V(0.75)[4,32] | 1/40 | -0.0023 | 0.970 (0.009) | pass | 0.960 | -0.0056 | 0.880 (0.032) | pass | 45 / 13.5 |
| G0 | zeta_V(0.75)[4,32] | 1/20 | -0.0028 | 0.958 (0.010) | pass | 0.960 | -0.0056 | 0.840 (0.037) | pass | 48 / 23 |
| G0 | zeta_V(0.75)[4,32] | 1/10 | -0.0032 | 0.943 (0.012) | pass | 0.990 | -0.0057 | 0.860 (0.035) | pass | 48 / 41 |
| G0 | zeta_V(0.75)[4,32] | 1/5 | -0.0035 | 0.940 (0.012) | pass | 0.970 | -0.0056 | 0.860 (0.035) | pass | 50 / 45 |
| G0 | zeta_V(0.75)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| G0 | zeta_V(0.75)[8,64] | 1/40 | -0.0041 | 0.960 (0.010) | pass | 0.960 | -0.0091 | 0.870 (0.034) | pass | 41 / 7 |
| G0 | zeta_V(0.75)[8,64] | 1/20 | -0.0050 | 0.948 (0.011) | pass | 0.930 | -0.0114 | 0.760 (0.043) | fail (coverage) | 45 / 13.5 |
| G0 | zeta_V(0.75)[8,64] | 1/10 | -0.0062 | 0.925 (0.013) | pass | 0.950 | -0.0118 | 0.710 (0.045) | fail (coverage) | 48 / 23 |
| G0 | zeta_V(0.75)[8,64] | 1/5 | -0.0070 | 0.907 (0.014) | pass | 0.940 | -0.0123 | 0.710 (0.045) | fail (coverage) | 48 / 41 |
| G0 | zeta_V(1)[4,32] | 1/80 | -0.0016 | 0.980 (0.007) | pass | 0.980 | -0.0030 | 0.960 (0.020) | pass | 41 / 7 |
| G0 | zeta_V(1)[4,32] | 1/40 | -0.0019 | 0.980 (0.007) | pass | 0.990 | -0.0042 | 0.950 (0.022) | pass | 45 / 13.5 |
| G0 | zeta_V(1)[4,32] | 1/20 | -0.0024 | 0.970 (0.009) | pass | 0.970 | -0.0041 | 0.920 (0.027) | pass | 48 / 23 |
| G0 | zeta_V(1)[4,32] | 1/10 | -0.0027 | 0.950 (0.011) | pass | 0.980 | -0.0041 | 0.910 (0.029) | pass | 48 / 41 |
| G0 | zeta_V(1)[4,32] | 1/5 | -0.0029 | 0.958 (0.010) | pass | 0.970 | -0.0039 | 0.920 (0.027) | pass | 50 / 45 |
| G0 | zeta_V(1)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| G0 | zeta_V(1)[8,64] | 1/40 | -0.0031 | 0.968 (0.009) | pass | 0.970 | -0.0063 | 0.880 (0.032) | pass | 41 / 7 |
| G0 | zeta_V(1)[8,64] | 1/20 | -0.0039 | 0.950 (0.011) | pass | 0.940 | -0.0082 | 0.820 (0.038) | pass | 45 / 13.5 |
| G0 | zeta_V(1)[8,64] | 1/10 | -0.0050 | 0.940 (0.012) | pass | 0.960 | -0.0085 | 0.840 (0.037) | pass | 48 / 23 |
| G0 | zeta_V(1)[8,64] | 1/5 | -0.0057 | 0.927 (0.013) | pass | 0.950 | -0.0088 | 0.780 (0.041) | pass | 48 / 41 |
| G0 | zeta_K(0)[4,32] | 1/80 | +0.0012 | 0.988 (0.006) | pass | 0.960 | +0.0001 | 0.970 (0.017) | pass | 41 / 7 |
| G0 | zeta_K(0)[4,32] | 1/40 | +0.0019 | 0.988 (0.006) | pass | 0.980 | -0.0001 | 0.990 (0.010) | pass | 45 / 13.5 |
| G0 | zeta_K(0)[4,32] | 1/20 | +0.0016 | 0.983 (0.007) | pass | 0.970 | -0.0002 | 0.980 (0.014) | pass | 48 / 23 |
| G0 | zeta_K(0)[4,32] | 1/10 | +0.0015 | 0.990 (0.005) | pass | 0.980 | -0.0008 | 0.960 (0.020) | pass | 48 / 41 |
| G0 | zeta_K(0)[4,32] | 1/5 | +0.0013 | 0.988 (0.006) | pass | 0.970 | -0.0006 | 0.970 (0.017) | pass | 50 / 45 |
| G0 | zeta_K(0)[8,64] | 1/80 | +0.0175 | 0.838 (0.018) | pass | 0.800 | +0.0159 | 0.880 (0.032) | pass | 23 / 0 |
| G0 | zeta_K(0)[8,64] | 1/40 | +0.0038 | 0.978 (0.007) | pass | 0.980 | +0.0002 | 0.960 (0.020) | pass | 41 / 7 |
| G0 | zeta_K(0)[8,64] | 1/20 | +0.0040 | 0.983 (0.007) | pass | 0.990 | +0.0000 | 0.990 (0.010) | pass | 45 / 13.5 |
| G0 | zeta_K(0)[8,64] | 1/10 | +0.0031 | 0.970 (0.009) | pass | 0.950 | -0.0003 | 0.990 (0.010) | pass | 48 / 23 |
| G0 | zeta_K(0)[8,64] | 1/5 | +0.0026 | 0.983 (0.007) | pass | 0.970 | -0.0001 | 0.980 (0.014) | pass | 48 / 41 |
| G0 | zeta_K(0.25)[4,32] | 1/80 | +0.0013 | 0.985 (0.006) | pass | 0.960 | -0.0000 | 0.970 (0.017) | pass | 41 / 7 |
| G0 | zeta_K(0.25)[4,32] | 1/40 | +0.0019 | 0.993 (0.004) | pass | 0.990 | -0.0001 | 0.990 (0.010) | pass | 45 / 13.5 |
| G0 | zeta_K(0.25)[4,32] | 1/20 | +0.0016 | 0.978 (0.007) | pass | 0.980 | -0.0002 | 0.980 (0.014) | pass | 48 / 23 |
| G0 | zeta_K(0.25)[4,32] | 1/10 | +0.0015 | 0.990 (0.005) | pass | 0.970 | -0.0006 | 0.970 (0.017) | pass | 48 / 41 |
| G0 | zeta_K(0.25)[4,32] | 1/5 | +0.0014 | 0.985 (0.006) | pass | 0.970 | -0.0005 | 0.970 (0.017) | pass | 50 / 45 |
| G0 | zeta_K(0.25)[8,64] | 1/80 | +0.0166 | 0.835 (0.019) | pass | 0.790 | +0.0148 | 0.870 (0.034) | pass | 23 / 0 |
| G0 | zeta_K(0.25)[8,64] | 1/40 | +0.0036 | 0.965 (0.009) | pass | 0.950 | +0.0000 | 0.960 (0.020) | pass | 41 / 7 |
| G0 | zeta_K(0.25)[8,64] | 1/20 | +0.0038 | 0.980 (0.007) | pass | 0.980 | -0.0001 | 0.980 (0.014) | pass | 45 / 13.5 |
| G0 | zeta_K(0.25)[8,64] | 1/10 | +0.0029 | 0.975 (0.008) | pass | 0.960 | -0.0004 | 0.990 (0.010) | pass | 48 / 23 |
| G0 | zeta_K(0.25)[8,64] | 1/5 | +0.0026 | 0.983 (0.007) | pass | 0.960 | -0.0001 | 0.980 (0.014) | pass | 48 / 41 |
| G0 | zeta_K(0.5)[4,32] | 1/80 | +0.0013 | 0.983 (0.007) | pass | 0.960 | -0.0001 | 0.970 (0.017) | pass | 41 / 7 |
| G0 | zeta_K(0.5)[4,32] | 1/40 | +0.0018 | 0.993 (0.004) | pass | 0.990 | -0.0001 | 0.970 (0.017) | pass | 45 / 13.5 |
| G0 | zeta_K(0.5)[4,32] | 1/20 | +0.0016 | 0.975 (0.008) | pass | 0.960 | -0.0001 | 0.980 (0.014) | pass | 48 / 23 |
| G0 | zeta_K(0.5)[4,32] | 1/10 | +0.0014 | 0.988 (0.006) | pass | 0.980 | -0.0004 | 0.980 (0.014) | pass | 48 / 41 |
| G0 | zeta_K(0.5)[4,32] | 1/5 | +0.0013 | 0.988 (0.006) | pass | 0.980 | -0.0003 | 0.980 (0.014) | pass | 50 / 45 |
| G0 | zeta_K(0.5)[8,64] | 1/80 | +0.0159 | 0.830 (0.019) | pass | 0.780 | +0.0140 | 0.890 (0.031) | pass | 23 / 0 |
| G0 | zeta_K(0.5)[8,64] | 1/40 | +0.0035 | 0.973 (0.008) | pass | 0.960 | +0.0000 | 0.970 (0.017) | pass | 41 / 7 |
| G0 | zeta_K(0.5)[8,64] | 1/20 | +0.0036 | 0.980 (0.007) | pass | 0.990 | -0.0001 | 0.980 (0.014) | pass | 45 / 13.5 |
| G0 | zeta_K(0.5)[8,64] | 1/10 | +0.0027 | 0.973 (0.008) | pass | 0.970 | -0.0004 | 0.990 (0.010) | pass | 48 / 23 |
| G0 | zeta_K(0.5)[8,64] | 1/5 | +0.0024 | 0.993 (0.004) | pass | 0.980 | +0.0002 | 0.990 (0.010) | pass | 48 / 41 |
| G0 | zeta_K(0.75)[4,32] | 1/80 | +0.0013 | 0.985 (0.006) | pass | 0.980 | -0.0003 | 0.970 (0.017) | pass | 41 / 7 |
| G0 | zeta_K(0.75)[4,32] | 1/40 | +0.0017 | 0.985 (0.006) | pass | 0.970 | -0.0002 | 0.990 (0.010) | pass | 45 / 13.5 |
| G0 | zeta_K(0.75)[4,32] | 1/20 | +0.0015 | 0.978 (0.007) | pass | 0.960 | +0.0000 | 1.000 (0.000) | pass | 48 / 23 |
| G0 | zeta_K(0.75)[4,32] | 1/10 | +0.0012 | 0.988 (0.006) | pass | 0.980 | -0.0001 | 0.980 (0.014) | pass | 48 / 41 |
| G0 | zeta_K(0.75)[4,32] | 1/5 | +0.0012 | 0.985 (0.006) | pass | 0.980 | -0.0000 | 0.980 (0.014) | pass | 50 / 45 |
| G0 | zeta_K(0.75)[8,64] | 1/80 | +0.0156 | 0.845 (0.018) | pass | 0.820 | +0.0134 | 0.870 (0.034) | pass | 23 / 0 |
| G0 | zeta_K(0.75)[8,64] | 1/40 | +0.0033 | 0.978 (0.007) | pass | 0.980 | -0.0000 | 0.980 (0.014) | pass | 41 / 7 |
| G0 | zeta_K(0.75)[8,64] | 1/20 | +0.0034 | 0.980 (0.007) | pass | 0.980 | -0.0001 | 1.000 (0.000) | pass | 45 / 13.5 |
| G0 | zeta_K(0.75)[8,64] | 1/10 | +0.0024 | 0.970 (0.009) | pass | 0.960 | -0.0002 | 0.990 (0.010) | pass | 48 / 23 |
| G0 | zeta_K(0.75)[8,64] | 1/5 | +0.0021 | 0.988 (0.006) | pass | 0.970 | +0.0005 | 0.990 (0.010) | pass | 48 / 41 |
| G0 | zeta_K(1)[4,32] | 1/80 | +0.0012 | 0.985 (0.006) | pass | 0.990 | -0.0005 | 0.990 (0.010) | pass | 41 / 7 |
| G0 | zeta_K(1)[4,32] | 1/40 | +0.0015 | 0.988 (0.006) | pass | 0.990 | -0.0004 | 0.990 (0.010) | pass | 45 / 13.5 |
| G0 | zeta_K(1)[4,32] | 1/20 | +0.0013 | 0.985 (0.006) | pass | 0.970 | +0.0001 | 0.990 (0.010) | pass | 48 / 23 |
| G0 | zeta_K(1)[4,32] | 1/10 | +0.0010 | 0.988 (0.006) | pass | 0.980 | +0.0003 | 0.980 (0.014) | pass | 48 / 41 |
| G0 | zeta_K(1)[4,32] | 1/5 | +0.0010 | 0.983 (0.007) | pass | 0.980 | +0.0003 | 0.990 (0.010) | pass | 50 / 45 |
| G0 | zeta_K(1)[8,64] | 1/80 | +0.0153 | 0.865 (0.017) | pass | 0.800 | +0.0129 | 0.860 (0.035) | pass | 23 / 0 |
| G0 | zeta_K(1)[8,64] | 1/40 | +0.0030 | 0.980 (0.007) | pass | 0.980 | -0.0003 | 0.970 (0.017) | pass | 41 / 7 |
| G0 | zeta_K(1)[8,64] | 1/20 | +0.0031 | 0.980 (0.007) | pass | 0.960 | -0.0002 | 1.000 (0.000) | pass | 45 / 13.5 |
| G0 | zeta_K(1)[8,64] | 1/10 | +0.0019 | 0.970 (0.009) | pass | 0.970 | -0.0000 | 0.990 (0.010) | pass | 48 / 23 |
| G0 | zeta_K(1)[8,64] | 1/5 | +0.0017 | 0.983 (0.007) | pass | 0.970 | +0.0008 | 1.000 (0.000) | pass | 48 / 41 |
| G0 | pf[4,32] | 1/80 | +0.0000 | 0.983 (0.007) | pass | 1.000 | -0.0006 | 0.970 (0.017) | pass | 41 / 7 |
| G0 | pf[4,32] | 1/40 | -0.0004 | 0.990 (0.005) | pass | 1.000 | -0.0003 | 0.950 (0.022) | pass | 45 / 13.5 |
| G0 | pf[4,32] | 1/20 | -0.0003 | 0.990 (0.005) | pass | 0.990 | +0.0003 | 0.980 (0.014) | pass | 48 / 23 |
| G0 | pf[4,32] | 1/10 | -0.0005 | 0.990 (0.005) | pass | 1.000 | +0.0011 | 0.980 (0.014) | pass | 48 / 41 |
| G0 | pf[4,32] | 1/5 | -0.0003 | 0.988 (0.006) | pass | 1.000 | +0.0009 | 0.990 (0.010) | pass | 50 / 45 |
| G0 | pf[8,64] | 1/80 | -0.0021 | 0.990 (0.005) | pass | 0.980 | -0.0030 | 0.970 (0.017) | pass | 23 / 0 |
| G0 | pf[8,64] | 1/40 | -0.0007 | 0.983 (0.007) | pass | 0.980 | -0.0004 | 0.980 (0.014) | pass | 41 / 7 |
| G0 | pf[8,64] | 1/20 | -0.0009 | 0.990 (0.005) | pass | 0.990 | -0.0002 | 0.980 (0.014) | pass | 45 / 13.5 |
| G0 | pf[8,64] | 1/10 | -0.0012 | 0.980 (0.007) | pass | 0.970 | +0.0003 | 1.000 (0.000) | pass | 48 / 23 |
| G0 | pf[8,64] | 1/5 | -0.0009 | 0.988 (0.006) | pass | 0.980 | +0.0010 | 0.990 (0.010) | pass | 48 / 41 |
| G0 | Rratio[4,32] | 1/80 | -0.0014 | 0.983 (0.007) | pass | 0.970 | -0.0018 | 0.970 (0.017) | pass | 23 / 23 |
| G0 | Rratio[4,32] | 1/40 | -0.0017 | 0.973 (0.008) | pass | 0.980 | -0.0023 | 1.000 (0.000) | pass | 41 / 41 |
| G0 | Rratio[4,32] | 1/20 | -0.0020 | 0.975 (0.008) | pass | 0.960 | -0.0024 | 0.990 (0.010) | pass | 45 / 45 |
| G0 | Rratio[4,32] | 1/10 | -0.0023 | 0.978 (0.007) | pass | 0.970 | -0.0026 | 1.000 (0.000) | pass | 48 / 48 |
| G0 | Rratio[4,32] | 1/5 | -0.0023 | 0.978 (0.007) | pass | 0.970 | -0.0026 | 1.000 (0.000) | pass | 48 / 48 |
| G0 | Rratio[8,64] | 1/80 | -0.0024 | 0.980 (0.007) | pass | 0.990 | -0.0027 | 0.970 (0.017) | pass | 13.5 / 13.5 |
| G0 | Rratio[8,64] | 1/40 | -0.0028 | 0.985 (0.006) | pass | 0.970 | -0.0028 | 0.960 (0.020) | pass | 23 / 23 |
| G0 | Rratio[8,64] | 1/20 | -0.0033 | 0.973 (0.008) | pass | 0.980 | -0.0034 | 0.980 (0.014) | pass | 41 / 41 |
| G0 | Rratio[8,64] | 1/10 | -0.0036 | 0.963 (0.009) | pass | 0.920 | -0.0035 | 0.980 (0.014) | pass | 45 / 45 |
| G0 | Rratio[8,64] | 1/5 | -0.0040 | 0.968 (0.009) | pass | 0.960 | -0.0040 | 0.970 (0.017) | pass | 48 / 48 |
| GT | zeta_V(0)[4,32] | 1/80 | -0.0032 | 0.955 (0.010) | pass | 0.940 | -0.0063 | 0.840 (0.037) | pass | 41 / 7 |
| GT | zeta_V(0)[4,32] | 1/40 | -0.0034 | 0.965 (0.009) | pass | 0.970 | -0.0072 | 0.870 (0.034) | pass | 46 / 14 |
| GT | zeta_V(0)[4,32] | 1/20 | -0.0039 | 0.943 (0.012) | pass | 0.960 | -0.0072 | 0.850 (0.036) | pass | 48 / 23 |
| GT | zeta_V(0)[4,32] | 1/10 | -0.0046 | 0.925 (0.013) | pass | 0.920 | -0.0070 | 0.860 (0.035) | pass | 48 / 41 |
| GT | zeta_V(0)[4,32] | 1/5 | -0.0049 | 0.940 (0.012) | pass | 0.940 | -0.0075 | 0.850 (0.036) | pass | 50 / 46 |
| GT | zeta_V(0)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| GT | zeta_V(0)[8,64] | 1/40 | -0.0063 | 0.927 (0.013) | pass | 0.930 | -0.0138 | 0.760 (0.043) | fail (coverage) | 41 / 7 |
| GT | zeta_V(0)[8,64] | 1/20 | -0.0068 | 0.940 (0.012) | pass | 0.940 | -0.0151 | 0.770 (0.042) | pass | 46 / 14 |
| GT | zeta_V(0)[8,64] | 1/10 | -0.0083 | 0.902 (0.015) | pass | 0.860 | -0.0159 | 0.700 (0.046) | fail (coverage) | 48 / 23 |
| GT | zeta_V(0)[8,64] | 1/5 | -0.0095 | 0.848 (0.018) | pass | 0.860 | -0.0163 | 0.690 (0.046) | fail (coverage) | 48 / 41 |
| GT | zeta_V(0.25)[4,32] | 1/80 | -0.0030 | 0.960 (0.010) | pass | 0.930 | -0.0060 | 0.880 (0.032) | pass | 41 / 7 |
| GT | zeta_V(0.25)[4,32] | 1/40 | -0.0033 | 0.965 (0.009) | pass | 0.980 | -0.0068 | 0.870 (0.034) | pass | 46 / 14 |
| GT | zeta_V(0.25)[4,32] | 1/20 | -0.0039 | 0.953 (0.011) | pass | 0.960 | -0.0068 | 0.810 (0.039) | pass | 48 / 23 |
| GT | zeta_V(0.25)[4,32] | 1/10 | -0.0045 | 0.930 (0.013) | pass | 0.920 | -0.0067 | 0.860 (0.035) | pass | 48 / 41 |
| GT | zeta_V(0.25)[4,32] | 1/5 | -0.0048 | 0.932 (0.013) | pass | 0.940 | -0.0071 | 0.870 (0.034) | pass | 50 / 46 |
| GT | zeta_V(0.25)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| GT | zeta_V(0.25)[8,64] | 1/40 | -0.0062 | 0.945 (0.011) | pass | 0.940 | -0.0132 | 0.760 (0.043) | fail (coverage) | 41 / 7 |
| GT | zeta_V(0.25)[8,64] | 1/20 | -0.0067 | 0.945 (0.011) | pass | 0.950 | -0.0145 | 0.770 (0.042) | pass | 46 / 14 |
| GT | zeta_V(0.25)[8,64] | 1/10 | -0.0082 | 0.892 (0.015) | pass | 0.840 | -0.0154 | 0.730 (0.044) | fail (coverage) | 48 / 23 |
| GT | zeta_V(0.25)[8,64] | 1/5 | -0.0094 | 0.858 (0.017) | pass | 0.830 | -0.0157 | 0.740 (0.044) | fail (coverage) | 48 / 41 |
| GT | zeta_V(0.5)[4,32] | 1/80 | -0.0026 | 0.965 (0.009) | pass | 0.930 | -0.0048 | 0.890 (0.031) | pass | 41 / 7 |
| GT | zeta_V(0.5)[4,32] | 1/40 | -0.0029 | 0.970 (0.009) | pass | 0.980 | -0.0055 | 0.880 (0.032) | pass | 46 / 14 |
| GT | zeta_V(0.5)[4,32] | 1/20 | -0.0035 | 0.965 (0.009) | pass | 0.970 | -0.0056 | 0.850 (0.036) | pass | 48 / 23 |
| GT | zeta_V(0.5)[4,32] | 1/10 | -0.0040 | 0.943 (0.012) | pass | 0.930 | -0.0056 | 0.910 (0.029) | pass | 48 / 41 |
| GT | zeta_V(0.5)[4,32] | 1/5 | -0.0043 | 0.925 (0.013) | pass | 0.920 | -0.0060 | 0.900 (0.030) | pass | 50 / 46 |
| GT | zeta_V(0.5)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| GT | zeta_V(0.5)[8,64] | 1/40 | -0.0057 | 0.945 (0.011) | pass | 0.930 | -0.0112 | 0.800 (0.040) | pass | 41 / 7 |
| GT | zeta_V(0.5)[8,64] | 1/20 | -0.0061 | 0.943 (0.012) | pass | 0.940 | -0.0125 | 0.780 (0.041) | pass | 46 / 14 |
| GT | zeta_V(0.5)[8,64] | 1/10 | -0.0074 | 0.905 (0.015) | pass | 0.880 | -0.0134 | 0.750 (0.043) | fail (coverage) | 48 / 23 |
| GT | zeta_V(0.5)[8,64] | 1/5 | -0.0086 | 0.858 (0.017) | pass | 0.850 | -0.0137 | 0.760 (0.043) | fail (coverage) | 48 / 41 |
| GT | zeta_V(0.75)[4,32] | 1/80 | -0.0022 | 0.978 (0.007) | pass | 0.950 | -0.0033 | 0.930 (0.026) | pass | 41 / 7 |
| GT | zeta_V(0.75)[4,32] | 1/40 | -0.0024 | 0.973 (0.008) | pass | 0.980 | -0.0038 | 0.960 (0.020) | pass | 46 / 14 |
| GT | zeta_V(0.75)[4,32] | 1/20 | -0.0029 | 0.968 (0.009) | pass | 0.970 | -0.0041 | 0.910 (0.029) | pass | 48 / 23 |
| GT | zeta_V(0.75)[4,32] | 1/10 | -0.0034 | 0.945 (0.011) | pass | 0.940 | -0.0043 | 0.910 (0.029) | pass | 48 / 41 |
| GT | zeta_V(0.75)[4,32] | 1/5 | -0.0036 | 0.940 (0.012) | pass | 0.940 | -0.0047 | 0.920 (0.027) | pass | 50 / 46 |
| GT | zeta_V(0.75)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| GT | zeta_V(0.75)[8,64] | 1/40 | -0.0047 | 0.958 (0.010) | pass | 0.950 | -0.0080 | 0.890 (0.031) | pass | 41 / 7 |
| GT | zeta_V(0.75)[8,64] | 1/20 | -0.0051 | 0.945 (0.011) | pass | 0.970 | -0.0092 | 0.830 (0.038) | pass | 46 / 14 |
| GT | zeta_V(0.75)[8,64] | 1/10 | -0.0062 | 0.922 (0.013) | pass | 0.900 | -0.0102 | 0.800 (0.040) | pass | 48 / 23 |
| GT | zeta_V(0.75)[8,64] | 1/5 | -0.0073 | 0.877 (0.016) | pass | 0.840 | -0.0106 | 0.860 (0.035) | pass | 48 / 41 |
| GT | zeta_V(1)[4,32] | 1/80 | -0.0018 | 0.983 (0.007) | pass | 0.990 | -0.0022 | 0.990 (0.010) | pass | 41 / 7 |
| GT | zeta_V(1)[4,32] | 1/40 | -0.0019 | 0.975 (0.008) | pass | 0.970 | -0.0025 | 0.960 (0.020) | pass | 46 / 14 |
| GT | zeta_V(1)[4,32] | 1/20 | -0.0024 | 0.970 (0.009) | pass | 0.980 | -0.0029 | 0.950 (0.022) | pass | 48 / 23 |
| GT | zeta_V(1)[4,32] | 1/10 | -0.0028 | 0.953 (0.011) | pass | 0.950 | -0.0032 | 0.940 (0.024) | pass | 48 / 41 |
| GT | zeta_V(1)[4,32] | 1/5 | -0.0029 | 0.953 (0.011) | pass | 0.940 | -0.0037 | 0.940 (0.024) | pass | 50 / 46 |
| GT | zeta_V(1)[8,64] | 1/80 | +nan | 0.000 (0.000) | - | 0.000 | +nan | 0.000 (0.000) | - | 23 / 0 |
| GT | zeta_V(1)[8,64] | 1/40 | -0.0036 | 0.968 (0.009) | pass | 0.970 | -0.0050 | 0.920 (0.027) | pass | 41 / 7 |
| GT | zeta_V(1)[8,64] | 1/20 | -0.0039 | 0.945 (0.011) | pass | 0.960 | -0.0061 | 0.880 (0.032) | pass | 46 / 14 |
| GT | zeta_V(1)[8,64] | 1/10 | -0.0048 | 0.945 (0.011) | pass | 0.940 | -0.0071 | 0.850 (0.036) | pass | 48 / 23 |
| GT | zeta_V(1)[8,64] | 1/5 | -0.0059 | 0.925 (0.013) | pass | 0.920 | -0.0077 | 0.880 (0.032) | pass | 48 / 41 |
| GT | zeta_K(0)[4,32] | 1/80 | +0.0009 | 0.975 (0.008) | pass | 0.980 | +0.0024 | 0.970 (0.017) | pass | 41 / 7 |
| GT | zeta_K(0)[4,32] | 1/40 | +0.0007 | 0.983 (0.007) | pass | 0.990 | +0.0030 | 0.970 (0.017) | pass | 46 / 14 |
| GT | zeta_K(0)[4,32] | 1/20 | +0.0005 | 0.985 (0.006) | pass | 1.000 | +0.0033 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(0)[4,32] | 1/10 | +0.0005 | 0.985 (0.006) | pass | 1.000 | +0.0031 | 0.990 (0.010) | pass | 48 / 41 |
| GT | zeta_K(0)[4,32] | 1/5 | +0.0003 | 0.988 (0.006) | pass | 1.000 | +0.0021 | 0.990 (0.010) | pass | 50 / 46 |
| GT | zeta_K(0)[8,64] | 1/80 | +0.0225 | 0.865 (0.017) | fail (bias) | 0.840 | +0.0245 | 0.870 (0.034) | fail (bias) | 23 / 0 |
| GT | zeta_K(0)[8,64] | 1/40 | +0.0006 | 0.980 (0.007) | pass | 0.980 | +0.0035 | 0.970 (0.017) | pass | 41 / 7 |
| GT | zeta_K(0)[8,64] | 1/20 | +0.0004 | 0.985 (0.006) | pass | 0.990 | +0.0059 | 1.000 (0.000) | pass | 46 / 14 |
| GT | zeta_K(0)[8,64] | 1/10 | +0.0004 | 0.985 (0.006) | pass | 1.000 | +0.0046 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(0)[8,64] | 1/5 | -0.0001 | 0.988 (0.006) | pass | 0.980 | +0.0040 | 0.990 (0.010) | pass | 48 / 41 |
| GT | zeta_K(0.25)[4,32] | 1/80 | +0.0011 | 0.975 (0.008) | pass | 0.980 | +0.0024 | 0.980 (0.014) | pass | 41 / 7 |
| GT | zeta_K(0.25)[4,32] | 1/40 | +0.0009 | 0.983 (0.007) | pass | 1.000 | +0.0033 | 0.970 (0.017) | pass | 46 / 14 |
| GT | zeta_K(0.25)[4,32] | 1/20 | +0.0008 | 0.983 (0.007) | pass | 1.000 | +0.0035 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(0.25)[4,32] | 1/10 | +0.0008 | 0.980 (0.007) | pass | 1.000 | +0.0032 | 0.990 (0.010) | pass | 48 / 41 |
| GT | zeta_K(0.25)[4,32] | 1/5 | +0.0006 | 0.983 (0.007) | pass | 1.000 | +0.0022 | 0.980 (0.014) | pass | 50 / 46 |
| GT | zeta_K(0.25)[8,64] | 1/80 | +0.0215 | 0.845 (0.018) | fail (bias) | 0.800 | +0.0232 | 0.850 (0.036) | fail (bias) | 23 / 0 |
| GT | zeta_K(0.25)[8,64] | 1/40 | +0.0010 | 0.985 (0.006) | pass | 0.990 | +0.0035 | 0.970 (0.017) | pass | 41 / 7 |
| GT | zeta_K(0.25)[8,64] | 1/20 | +0.0007 | 0.985 (0.006) | pass | 0.990 | +0.0060 | 0.990 (0.010) | pass | 46 / 14 |
| GT | zeta_K(0.25)[8,64] | 1/10 | +0.0007 | 0.983 (0.007) | pass | 1.000 | +0.0047 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(0.25)[8,64] | 1/5 | +0.0002 | 0.993 (0.004) | pass | 1.000 | +0.0038 | 0.980 (0.014) | pass | 48 / 41 |
| GT | zeta_K(0.5)[4,32] | 1/80 | +0.0012 | 0.980 (0.007) | pass | 0.990 | +0.0023 | 0.990 (0.010) | pass | 41 / 7 |
| GT | zeta_K(0.5)[4,32] | 1/40 | +0.0010 | 0.980 (0.007) | pass | 1.000 | +0.0034 | 0.980 (0.014) | pass | 46 / 14 |
| GT | zeta_K(0.5)[4,32] | 1/20 | +0.0009 | 0.988 (0.006) | pass | 1.000 | +0.0036 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(0.5)[4,32] | 1/10 | +0.0009 | 0.983 (0.007) | pass | 1.000 | +0.0034 | 0.990 (0.010) | pass | 48 / 41 |
| GT | zeta_K(0.5)[4,32] | 1/5 | +0.0007 | 0.980 (0.007) | pass | 1.000 | +0.0024 | 0.980 (0.014) | pass | 50 / 46 |
| GT | zeta_K(0.5)[8,64] | 1/80 | +0.0206 | 0.833 (0.019) | fail (bias) | 0.780 | +0.0220 | 0.840 (0.037) | fail (bias) | 23 / 0 |
| GT | zeta_K(0.5)[8,64] | 1/40 | +0.0011 | 0.978 (0.007) | pass | 0.980 | +0.0036 | 0.960 (0.020) | pass | 41 / 7 |
| GT | zeta_K(0.5)[8,64] | 1/20 | +0.0009 | 0.990 (0.005) | pass | 0.990 | +0.0060 | 0.990 (0.010) | pass | 46 / 14 |
| GT | zeta_K(0.5)[8,64] | 1/10 | +0.0008 | 0.983 (0.007) | pass | 1.000 | +0.0048 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(0.5)[8,64] | 1/5 | +0.0003 | 0.990 (0.005) | pass | 1.000 | +0.0037 | 0.980 (0.014) | pass | 48 / 41 |
| GT | zeta_K(0.75)[4,32] | 1/80 | +0.0010 | 0.980 (0.007) | pass | 0.990 | +0.0020 | 0.960 (0.020) | pass | 41 / 7 |
| GT | zeta_K(0.75)[4,32] | 1/40 | +0.0008 | 0.980 (0.007) | pass | 1.000 | +0.0033 | 0.950 (0.022) | pass | 46 / 14 |
| GT | zeta_K(0.75)[4,32] | 1/20 | +0.0008 | 0.985 (0.006) | pass | 1.000 | +0.0034 | 0.980 (0.014) | pass | 48 / 23 |
| GT | zeta_K(0.75)[4,32] | 1/10 | +0.0008 | 0.985 (0.006) | pass | 1.000 | +0.0034 | 0.990 (0.010) | pass | 48 / 41 |
| GT | zeta_K(0.75)[4,32] | 1/5 | +0.0005 | 0.985 (0.006) | pass | 1.000 | +0.0026 | 0.980 (0.014) | pass | 50 / 46 |
| GT | zeta_K(0.75)[8,64] | 1/80 | +0.0197 | 0.843 (0.018) | pass | 0.800 | +0.0210 | 0.820 (0.038) | fail (bias) | 23 / 0 |
| GT | zeta_K(0.75)[8,64] | 1/40 | +0.0009 | 0.980 (0.007) | pass | 0.980 | +0.0034 | 0.970 (0.017) | pass | 41 / 7 |
| GT | zeta_K(0.75)[8,64] | 1/20 | +0.0007 | 0.990 (0.005) | pass | 0.990 | +0.0057 | 0.990 (0.010) | pass | 46 / 14 |
| GT | zeta_K(0.75)[8,64] | 1/10 | +0.0007 | 0.988 (0.006) | pass | 1.000 | +0.0047 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(0.75)[8,64] | 1/5 | +0.0001 | 0.983 (0.007) | pass | 0.990 | +0.0036 | 0.980 (0.014) | pass | 48 / 41 |
| GT | zeta_K(1)[4,32] | 1/80 | +0.0008 | 0.985 (0.006) | pass | 1.000 | +0.0017 | 0.970 (0.017) | pass | 41 / 7 |
| GT | zeta_K(1)[4,32] | 1/40 | +0.0005 | 0.980 (0.007) | pass | 1.000 | +0.0029 | 0.950 (0.022) | pass | 46 / 14 |
| GT | zeta_K(1)[4,32] | 1/20 | +0.0006 | 0.980 (0.007) | pass | 1.000 | +0.0031 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(1)[4,32] | 1/10 | +0.0005 | 0.990 (0.005) | pass | 1.000 | +0.0033 | 0.990 (0.010) | pass | 48 / 41 |
| GT | zeta_K(1)[4,32] | 1/5 | +0.0002 | 0.980 (0.007) | pass | 1.000 | +0.0026 | 1.000 (0.000) | pass | 50 / 46 |
| GT | zeta_K(1)[8,64] | 1/80 | +0.0188 | 0.868 (0.017) | pass | 0.840 | +0.0200 | 0.830 (0.038) | pass | 23 / 0 |
| GT | zeta_K(1)[8,64] | 1/40 | +0.0004 | 0.973 (0.008) | pass | 0.990 | +0.0029 | 0.980 (0.014) | pass | 41 / 7 |
| GT | zeta_K(1)[8,64] | 1/20 | +0.0002 | 0.988 (0.006) | pass | 1.000 | +0.0050 | 0.980 (0.014) | pass | 46 / 14 |
| GT | zeta_K(1)[8,64] | 1/10 | +0.0003 | 0.988 (0.006) | pass | 1.000 | +0.0044 | 0.990 (0.010) | pass | 48 / 23 |
| GT | zeta_K(1)[8,64] | 1/5 | -0.0005 | 0.983 (0.007) | pass | 1.000 | +0.0032 | 0.980 (0.014) | pass | 48 / 41 |
| GT | pf[4,32] | 1/80 | -0.0002 | 0.985 (0.006) | pass | 0.990 | -0.0007 | 0.980 (0.014) | pass | 41 / 7 |
| GT | pf[4,32] | 1/40 | -0.0002 | 0.983 (0.007) | pass | 1.000 | -0.0000 | 0.990 (0.010) | pass | 46 / 14 |
| GT | pf[4,32] | 1/20 | +0.0001 | 0.985 (0.006) | pass | 1.000 | -0.0003 | 1.000 (0.000) | pass | 48 / 23 |
| GT | pf[4,32] | 1/10 | +0.0000 | 0.980 (0.007) | pass | 0.980 | +0.0003 | 1.000 (0.000) | pass | 48 / 41 |
| GT | pf[4,32] | 1/5 | -0.0001 | 0.988 (0.006) | pass | 1.000 | +0.0005 | 1.000 (0.000) | pass | 50 / 46 |
| GT | pf[8,64] | 1/80 | -0.0037 | 0.980 (0.007) | pass | 0.970 | -0.0045 | 0.980 (0.014) | pass | 23 / 0 |
| GT | pf[8,64] | 1/40 | -0.0002 | 0.978 (0.007) | pass | 0.990 | -0.0005 | 0.980 (0.014) | pass | 41 / 7 |
| GT | pf[8,64] | 1/20 | -0.0002 | 0.985 (0.006) | pass | 1.000 | -0.0008 | 0.990 (0.010) | pass | 46 / 14 |
| GT | pf[8,64] | 1/10 | -0.0001 | 0.985 (0.006) | pass | 1.000 | -0.0001 | 1.000 (0.000) | pass | 48 / 23 |
| GT | pf[8,64] | 1/5 | -0.0003 | 0.985 (0.006) | pass | 0.980 | -0.0008 | 0.990 (0.010) | pass | 48 / 41 |
| GT | Rratio[4,32] | 1/80 | -0.0017 | 0.978 (0.007) | pass | 0.970 | -0.0022 | 0.970 (0.017) | pass | 23 / 23 |
| GT | Rratio[4,32] | 1/40 | -0.0022 | 0.968 (0.009) | pass | 0.980 | -0.0024 | 0.980 (0.014) | pass | 41 / 41 |
| GT | Rratio[4,32] | 1/20 | -0.0024 | 0.960 (0.010) | pass | 1.000 | -0.0024 | 0.990 (0.010) | pass | 46 / 46 |
| GT | Rratio[4,32] | 1/10 | -0.0028 | 0.968 (0.009) | pass | 1.000 | -0.0027 | 0.990 (0.010) | pass | 48 / 48 |
| GT | Rratio[4,32] | 1/5 | -0.0028 | 0.968 (0.009) | pass | 1.000 | -0.0027 | 0.990 (0.010) | pass | 48 / 48 |
| GT | Rratio[8,64] | 1/80 | -0.0027 | 0.970 (0.009) | pass | 0.980 | -0.0027 | 0.980 (0.014) | pass | 14 / 14 |
| GT | Rratio[8,64] | 1/40 | -0.0029 | 0.970 (0.009) | pass | 0.970 | -0.0037 | 0.960 (0.020) | pass | 23 / 23 |
| GT | Rratio[8,64] | 1/20 | -0.0036 | 0.975 (0.008) | pass | 0.980 | -0.0039 | 0.990 (0.010) | pass | 41 / 41 |
| GT | Rratio[8,64] | 1/10 | -0.0041 | 0.975 (0.008) | pass | 1.000 | -0.0035 | 1.000 (0.000) | pass | 46 / 46 |
| GT | Rratio[8,64] | 1/5 | -0.0045 | 0.975 (0.008) | pass | 1.000 | -0.0039 | 1.000 (0.000) | pass | 48 / 48 |

