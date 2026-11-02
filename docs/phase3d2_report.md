# Phase 3d-2 report (final): operating characteristics of the four primary tests

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken in this phase. Real data
were read only through `worlds3d2`'s flow cache: the exploration stocks of the six confirmatory cells, 8
study days, and only their signs, sizes and timestamps. No replication or confirmation file was opened.
Every price response is synthetic.

Earlier write-ups:
- `docs/phase3d2a_report.md`: design, pilot and driver, before the study ran;
- `docs/phase3d2b_report.md`: diagnostics;
- `docs/phase3d2c_report.md`: the decisions of 1 Nov 2026 and the re-pilot.

The decisions are recorded in `docs/phase3d_tests.md` section 7. Tables: `~/CODE/itch-data/phase3d2/report/tables.md`
(8-test family, with the former 11-test family beside it) and `report/family11/` (the tables as first
written, 11-test family).

## Headlines

1. **Size holds, with room to spare (Measurement).**
   - In W-null (200 panels), no test rejects in any tick bin, before or after Holm. The family-wise error
     is 0 of 200 (95 % upper bound 0.015).
   - In W0, P-M rejects in 0 of 200 panels.
   - The tests are therefore conservative on these panels. Their size is well below 0.05.
2. **Power is low after Holm (Measurement, 8-test family).**
   - W0's P-R establishes the equivalence in 0.185 of panels in Q-T1 and 0.235 in Q-T3. Before Holm the
     rates are 0.770 and 0.750.
   - WB's P-M passes in 0.020 (Q-T1) and 0.220 (Q-T3). Before Holm: 0.250 and 0.900.
   - W0-psi3's P-R fails in every panel, which is its correct outcome.
   - P-V has no power world in the study.
3. **Q-T2 is descriptive (decision of 2 Nov 2026, before any real-data test).**
   - The Holm family is now the four tests in Q-T1 and Q-T3, up to 8 p-values.
   - Recomputed from the saved p-values, every labelled rate after Holm is the same as under the former
     11-test family or higher (up to +0.05), and the family-wise error stays 0. Section 5.
4. **The P-M design is b3_l16_L128_m50** (3 bins, l_max 16, L 128, 50 runs per bin), chosen by the
   pre-stated rule. Every one of the 6 designs was admissible. Section 2.
5. **psi-hat is far too high on real flow (Measurement).**
   - At psi = 0.5 (W0) it reads 0.63 to 0.70. At psi = 0.3 (W0-psi3) it reads 0.35 to 0.41.
   - Without sign memory (3d-2b's W0-shuf) the bias was only +0.02 to +0.05.
   - `primary3d.PSI_NOTE` states only that smaller figure. Section 6.
6. **Exposure: P-V uses signs and sizes only.** So its statistics in the semi-synthetic worlds were
   computed on resampled real exploration flow. They are, in effect, the exploration sample's own P-V
   statistics, and they have now been seen. Section 7.
7. **P-F rejects often in worlds with no P-F label (Measurement).** In W0 Q-T3 it rejects after Holm in 0.875
   of panels, and in W0-psi3 Q-T3 in 0.980. These are not scored. zeta_K falls with a in these worlds (slope
   -0.10 in W0 Q-T3), which have real signs and a size-independent transient kernel. Section 4.

## 1. What was run

**Worlds** (`worlds3d2`; `docs/phase3d2a_report.md` section 2):

| world | response | panels | correct outcome |
|---|---|---|---|
| W-null | real flow, sizes permuted within each stock-day; GF with fair-coin informed signs (3d-2b), psi = 0.5 | 200 | P-V, P-M, P-F no reject |
| W0 | real flow; G0 (one size-independent kernel), psi = 0.5 | 200 | P-R pass, P-M no reject |
| W0-psi3 | as W0 with psi = 0.3 | 100 | P-R fail |
| WB | real flow; beta(q) falls in log q from 0.6 to -0.3, psi = 0.5 | 100 | P-M reject |

- **Panels.** Each panel resamples each confirmatory cell's exploration stocks, with replacement, up to the
  cell's confirmation count. That gives Q-T1 377 stocks, Q-T2 17 and Q-T3 100.
- **What was not run.** WM (the published PAPER2 preset; its P-M fails by construction) and the thinning
  grid are separate commands and were not run (decision 1 Nov 2026, `docs/phase3d_tests.md` 7.5).
- **Run.**
  - Stage A: W-null, W0 and WB at the 6 P-M designs (the grid trimmed from 18 on 1 Nov 2026 for compute
    cost only, 7.6).
  - Then the design rule.
  - Stage B: W0-psi3 at the chosen design.
  - 900 chunks of 2 panels, all complete (`phase3d2.py check`, 2 Nov 2026).
  - The run resumed after the double-run incident of 1 Nov 2026 (7.6): 130 Q-T1 chunks were reused,
    trimmed to the 6 designs. The remaining 620 + 150 jobs and the analysis took 5.7 h on 3 workers.

**Rates (Choice, fixed with the rule before any result).**
- **Before Holm:** the pair is evaluable and p <= 0.05.
- **After Holm:** the panel's one family rejects it.
- **Passed:** rejected after Holm and every sub-condition holds. For P-R, that is the equivalence
  established.
- A pair that is not evaluable counts as not rejected.
- MC SE is sqrt(r (1 - r) / n). At r = 0 it prints 0.000; the 95 % upper bound is then 0.015 at 200 panels
  and 0.030 at 100.

## 2. The P-M design choice (Measurement, rule fixed beforehand)

The rule (`worlds3d2.choose_design`):
- A design is admissible when, in every tick bin, its unadjusted pass rate is at most 0.05 + 2 MC SE in W0
  and in W-null.
- Among admissible designs, it picks the highest WB pass rate averaged over the three tick bins.
- Ties go to fewer bins, then smaller l_max.

| design | WB Q-T1 | WB Q-T2 | WB Q-T3 | WB objective | worst rate in W0 and W-null |
|---|---|---|---|---|---|
| **b3_l16_L128_m50** | 0.250 | 0.000 | 0.900 | **0.383** | 0.000 |
| b4_l16_L128_m50 | 0.270 | 0.000 | 0.680 | 0.317 | 0.000 |
| b5_l16_L128_m50 | 0.040 | 0.000 | 0.430 | 0.157 | 0.000 |
| b4_l32_L256_m50 | 0.030 | 0.000 | 0.390 | 0.140 | 0.000 |
| b5_l32_L256_m50 | 0.000 | 0.000 | 0.420 | 0.140 | 0.000 |
| b3_l32_L256_m50 | 0.320 | 0.000 | 0.090 | 0.137 | 0.000 |

- All 6 designs are admissible, and none has a single unadjusted rejection in W0 or W-null.
- P-M is never evaluable in Q-T2 (no 5 stocks reach L under the ceilings), so Q-T2 adds 0 to every
  design's objective. The Q-T2 decision of 2 Nov 2026 therefore leaves the ranking unchanged, and
  `report/design.json` is byte-identical before and after.

## 3. Size (Measurement)

**W-null, 200 panels.**
- Every test, every tick bin: 0.000 before Holm and 0.000 after, in both families.
- Family-wise error (P-V, P-M, P-F, after Holm): 0.000 in the 8-test family and 0.000 in the former
  11-test family.
- Any rejection including P-R: 0.000.
- P-R has no label in W-null. Its a* sits near the grid end (0.85 to 0.99), and it never passes.

**W0, P-M (no reject):** 0.000 in every tick bin, before and after Holm.

Every null rejection count is 0. So the measured size is 0, with a 95 % upper bound of 0.015 per tick bin
and for the family. The study cannot say how far below 0.05 the true size lies, only that it is not above
about 0.015.

## 4. Power and correct outcomes, before and after Holm (Measurement)

Rates over panels, MC SE in brackets. "8" is the family of Q-T1 and Q-T3 (2 Nov 2026), "11" the former
family of all three tick bins. Both come from the same saved p-values.

| world | test | correct outcome | tick bin | before Holm | after Holm, 8 | after Holm, 11 | rate correct, 8 | rate correct, 11 |
|---|---|---|---|---|---|---|---|---|
| W0 | P-R | pass | Q-T1 | 0.770 (0.030) | 0.185 (0.027) | 0.145 (0.025) | 0.185 (0.027) | 0.145 (0.025) |
| W0 | P-R | pass | Q-T2 | 0.000 (0.000) | descriptive | 0.000 (0.000) | descriptive | 0.000 (0.000) |
| W0 | P-R | pass | Q-T3 | 0.750 (0.031) | 0.235 (0.030) | 0.220 (0.029) | 0.235 (0.030) | 0.220 (0.029) |
| W0-psi3 | P-R | fail | Q-T1 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) |
| W0-psi3 | P-R | fail | Q-T2 | 0.000 (0.000) | descriptive | 0.000 (0.000) | descriptive | 1.000 (0.000) |
| W0-psi3 | P-R | fail | Q-T3 | 0.030 (0.017) | 0.000 (0.000) | 0.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) |
| WB | P-M | reject | Q-T1 | 0.250 (0.043) | 0.020 (0.014) | 0.020 (0.014) | 0.020 (0.014) | 0.020 (0.014) |
| WB | P-M | reject | Q-T2 | not evaluable | descriptive | 0.000 (0.000) | descriptive | 0.000 (0.000) |
| WB | P-M | reject | Q-T3 | 0.900 (0.030) | 0.220 (0.041) | 0.170 (0.038) | 0.220 (0.041) | 0.170 (0.038) |

- For P-R and P-M, "after Holm" and "passed" agree in every row. Every rejected pair also met its
  sub-conditions.
- **P-V** has no power world in the study. WM was P-V's power world and is not part of `all` (7.5). P-V
  never rejected after Holm in any world.
- **W0-psi3's correct "fail" rests on width, not location.** Its a* is 0.38 / 0.36 / 0.41 at T1 and
  0.35 / 0.32 / 0.37 at T2 (Q-T1 / Q-T2 / Q-T3). That is at or just inside the band's lower edge, 0.35. P-R
  fails because the interval reaches outside the band.
- **Holm costs most of the power.** Before Holm, P-R passes in three quarters of W0's panels and P-M in 0.90
  of WB's Q-T3 panels. After Holm, about one panel in five passes.

**Rejections in pairs without a label (reported, not scored).**

| world | test | tick bin | before Holm | after Holm, 8 | after Holm, 11 | P-F slope of zeta_K on a, mean (SD) |
|---|---|---|---|---|---|---|
| W0 | P-F | Q-T1 | 0.330 (0.033) | 0.025 (0.011) | 0.020 (0.010) | -0.015 (0.006) |
| W0 | P-F | Q-T3 | 1.000 (0.000) | 0.875 (0.023) | 0.855 (0.025) | -0.100 (0.015) |
| W0-psi3 | P-F | Q-T1 | 0.370 (0.048) | 0.100 (0.030) | 0.080 (0.027) | -0.015 (0.006) |
| W0-psi3 | P-F | Q-T3 | 1.000 (0.000) | 0.980 (0.014) | 0.970 (0.017) | -0.119 (0.015) |
| WB | P-F | Q-T2 | 0.580 (0.049) | descriptive | 0.080 (0.027) | +0.084 (0.019) |
| WB | P-F | Q-T3 | 0.590 (0.049) | 0.010 (0.010) | 0.010 (0.010) | +0.063 (0.009) |
| W0 | P-V | Q-T3 | 0.010 (0.007) | 0.000 (0.000) | 0.000 (0.000) | |
| W0-psi3 | P-V | Q-T3 | 0.040 (0.020) | 0.000 (0.000) | 0.000 (0.000) | |

- W0 and W0-psi3 have real signs (with their memory) and a transient kernel, and zeta_K falls with a there.
  Under W-null's clean null, zeta_K is flat and P-F rejects nothing.
- So a P-F rejection on real data is consistent with transient impact on correlated real signs, without any
  size dependence of the kernel. It is not specific to the alternative P-F is meant to detect.
- Any rejection in the family, P-R included: W0 0.905 (8) against 0.875 (11); W0-psi3 0.980 against 0.970;
  WB 0.230 against 0.240.

## 5. The Q-T2 decision (2 Nov 2026; `docs/phase3d_tests.md` 7.7)

- **Decision.** Q-T2 leaves the confirmatory family and becomes descriptive. The Holm family is the four
  primary tests in Q-T1 and Q-T3 (`prereg3d.FAMILY_TICK_BINS`), up to 8 p-values.
- **When.** It was taken after this study and before any real-data test.
- **Reason.** Q-T2's only confirmatory cell (Q-T2S3) has 17 stocks, and this study measured near-zero power
  for every test there:
  - W0's P-R passes in 0 of 200 panels, even before Holm. a*(T2)'s SD over panels is 0.096, against 0.03
    in Q-T1 and Q-T3.
  - P-M is never evaluable at the chosen design.
  - P-V has no power world.
- **What still runs.** Q-T2's unit is still computed and reported, with before-Holm rates only.
- **Effect on the study, recomputed from the saved per-panel p-values (no new simulation):**

| rate after Holm | 11-test family | 8-test family |
|---|---|---|
| mean family size | 11.00 | 8.00 |
| W-null family-wise error | 0.000 | 0.000 |
| W0 P-M no reject (every tick bin) | 1.000 | 1.000 |
| W0 P-R pass, Q-T1 | 0.145 (0.025) | 0.185 (0.027) |
| W0 P-R pass, Q-T3 | 0.220 (0.029) | 0.235 (0.030) |
| W0-psi3 P-R fail, Q-T1 and Q-T3 | 1.000 | 1.000 |
| WB P-M pass, Q-T1 | 0.020 (0.014) | 0.020 (0.014) |
| WB P-M pass, Q-T3 | 0.170 (0.038) | 0.220 (0.041) |

- **The 11-test rebuild matches the archive.** `primary3d.step_down` rebuilt the 11-test column from the
  same p-values, and it reproduces the archived tables (`report/family11/results.json`) in every row of
  every world: before-Holm, after-Holm, passed and correct rates, the family-wise error and the family
  sizes.
- The gains come only from the smaller family. The same p-values meet looser Holm thresholds: 0.05/8 at
  the first step, against 0.05/11.

## 6. psi-hat on real flow (Measurement)

psi-hat is P-R's pre-registered secondary: per stock, the slope of log E[Delta_1 eps | q-bin] on log
median q, averaged over stocks. Mean over panels (SD over panels):

| world | true psi | Q-T1 | Q-T2 | Q-T3 |
|---|---|---|---|---|
| W0 (real signs, G0 kernel) | 0.5 | 0.665 (0.011) | 0.626 (0.032) | 0.702 (0.024) |
| W0-psi3 (as W0) | 0.3 | 0.382 (0.005) | 0.347 (0.008) | 0.413 (0.010) |
| WB (real signs, size-dependent decay) | 0.5 | 0.531 (0.010) | 0.637 (0.042) | 0.547 (0.016) |
| W0-shuf (3d-2b, signs shuffled within day) | 0.5 | 0.553 | 0.532 | 0.521 |

- **On real flow at psi = 0.5, psi-hat reads 0.63 to 0.70**, that is +0.13 to +0.20. At psi = 0.3 it reads
  0.35 to 0.41 (+0.05 to +0.11).
- Shuffling the signs within the day removes most of the bias (W0-shuf), so it comes from the real signs'
  memory.
- With sizes independent of signs, past runs' decaying impact adds only a q-independent term to
  E[Delta_1 eps | q] (`docs/phase3d_tests.md` section 6: 0.537 on synthetic G0 flow). A bias this large
  needs that term to depend on q in the real flow, that is, sizes related to sign memory. That mechanism was
  not measured here.
- psi-hat is reported, never tested (decision 7.4: no correction). `primary3d.PSI_NOTE`, printed beside
  every psi-hat, gives only the shuffled-sign figure (+0.02 to +0.05). On real flow it understates the bias
  (open item 1).
- For comparison, a*(T1) in W0 is 0.58 / 0.57 / 0.60 (Q-T1 / Q-T2 / Q-T3), +0.07 to +0.10. So P-R's a* is
  much less affected than psi-hat.

## 7. Exposure note: P-V on resampled exploration flow

- **P-V depends only on the flow.** Its statistic is built from v_a(T), the variance of the summed
  eps q^a over T runs. That uses signs and sizes only; no price enters. Its surrogates S1 and S2 permute
  sizes, so they are flow-only as well.
- **So the study computed it on real data.** In W0, W0-psi3 and WB, P-V's gate (i) and its surrogate slopes
  (ii) were computed on resampled real exploration flow. The synthetic response never enters them.
- **What has now been seen** is, in effect, the exploration sample's own P-V statistics. The slopes of
  zeta_V real minus surrogate, mean over W0's panels (S1 / S2):
  - Q-T1: +0.043 / +0.047;
  - Q-T2: +0.065 / +0.071;
  - Q-T3: -0.032 / -0.015.

  W0-psi3 and WB give the same to within 0.005. Their panels are other resamples of the same stocks.
- **P-V's outcome on that flow.**
  - In Q-T1 and Q-T2 the slopes are positive, where P-V predicts negative.
  - In Q-T3 they are negative. Yet P-V rejects before Holm in only 0.010 (W0) and 0.040 (W0-psi3) of
    panels, and after Holm in none.
- **What this means.**
  - The exploration flow, resampled to confirmation sizes, does not support P-V.
  - This is a property of the exploration sample, which was already exploratory.
  - The confirmation stocks' flow has not been read.
  - In W-null, sizes are permuted within each stock-day, so its P-V slopes are 0 to within 0.003 by
    construction.
- No other primary test is affected: P-R, P-M and P-F all use the synthetic price.

## 8. Tests (Measurement)

Full suite: 471 passed, 28 skipped (2 Nov 2026).

New for the Q-T2 decision:
- **Checks.**
  - `test_3d1.check_qt2_never_in_holm`: holm refuses a Q-T2 pair; combine lists Q-T2 as descriptive and
    refuses a Q-T2 p-value; the family of 8 rejects p = 0.006, where 11 members would not.
  - `check_run_primary`, updated: Q-T2 is computed and descriptive, with rejected and passed None.
  - `test_3d2.check_analysis_family_without_qt2`, on saved-record toys.
- **Planted bugs, all caught:**
  - holm's family check widened to every tick bin;
  - combine making Q-T2 pairs members;
  - run_primary passing Q-T2's p-values;
  - the analysis family over every tick bin;
  - Q-T2 scored as part of the family;
  - holm and combine widened together, so no guard raises.
- **Updated tests.** The Holm known answers in `test_3d1` (one family, pairs leaving it, a p = 1 pair) are
  restated for 8 members. Three Holm planted bugs now mutate `step_down`, where the step-down moved.

## 9. Open (none acted on)

1. **`PSI_NOTE` understates psi-hat's bias on real flow** (+0.13 to +0.20 at psi = 0.5, section 6). The
   note could cite both figures. psi-hat stays uncorrected either way.
2. **P-F's rejections without a label** (section 4). On real signs with a transient kernel, P-F rejects in
   Q-T3 in most panels. How to read a real-data P-F rejection needs a statement before confirmation.
3. **P-V has no power world**, and its exploration-flow outcome has now been seen (section 7).
4. **Power after Holm is about 0.2** for the labelled alternatives (section 4), against size well under
   0.05 (section 3).
5. **WM and thinning** were not run (`phase3d2.py wm`, `phase3d2.py thin`).

## 10. Reproduce

```
cd ~/CODE/itchbook/tools/est
uv run --group analysis phase3d2.py check     # integrity of the 900 saved chunks
uv run --group analysis phase3d2.py analyse   # tables from the saved p-values (seconds; no simulation)
uv run --group dev --group analysis pytest ~/CODE/itchbook/tools/est -q
```
