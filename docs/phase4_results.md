# Phase 4: the confirmatory run

Run once on 2026-11-02 10:50:01 -0700 with `tools/est/phase4.py confirm`. The plan is frozen at tag prereg-freeze-v1 (commit d126088fba041824950d7ffae4aa75a7a98b2703) and registered on OSF (osf.io/z4t7x; registered PDF SHA-256 79c8984ab832c0a9a0d858ac46959277399ce6394106fb526fe9bd48779ccafc).

- git HEAD at the run: d126088fba041824950d7ffae4aa75a7a98b2703; files added after the freeze: tools/est/phase4.py, tools/est/test_phase4.py
- data/manifest.tsv 9780239827c00763fd542ad43cd9a474fbfe815a0418847a95a6328804bc5a9d; data/universe.tsv 5fb7e1b7061ea8433d1916a3ca54c5f97a290e29b6c367f70881b248a3b42d25 (both as registered)
- bootstrap seed 20261102 (numpy default_rng, one generator through run_primary); B = 200; surrogate draws 4
- P-M design: {'n_bins': 3, 'l_max': 16, 'L': 128, 'min_runs_per_bin': 50}
- wall time 40 s (loading 6 s)
- outputs: data/phase4/results.json (every output), data/phase4/tables.md

Labels: **Confirmatory** = a member of the one Holm family; **Descriptive** = computed and reported, never in Holm. Secondary and robustness analyses are not part of this step.

## The family

One Holm-Bonferroni family at family-wise 0.05 over the evaluable (test, tick bin) pairs of Q-T1 and Q-T3: m = 8.


| rank | test | tick bin | range used | p | threshold | adjusted p | Holm rejects | passed |
|---|---|---|---|---|---|---|---|---|
| 1 | P-F | Q-T1 | [8,64] | 1.43e-04 | 0.0063 | 0.0011 | yes | yes |
| 2 | P-F | Q-T3 | [8,64] | 0.0048 | 0.0071 | 0.0335 | yes | yes |
| 3 | P-R | Q-T3 | [8,64] | 0.4450 | 0.0083 | 1.0000 | no | no |
| 4 | P-V | Q-T3 | [8,64] | 0.9700 | 0.0100 | 1.0000 | no | no |
| 5 | P-V | Q-T1 | [8,64] | 1.0000 | 0.0125 | 1.0000 | no | no |
| 6 | P-R | Q-T1 | [8,64] | 1.0000 | 0.0167 | 1.0000 | no | no |
| 7 | P-M | Q-T1 | lags [1,16], L = 128 | 1.0000 | 0.0250 | 1.0000 | no | no |
| 8 | P-M | Q-T3 | lags [1,16], L = 128 | 1.0000 | 0.0500 | 1.0000 | no | no |

## Each test (confirmatory)

### P-V, Q-T1

- range: primary [8,64]; stocks at T2 of each range: {'[4,32]': 107, '[8,64]': 64}
- p = 1.0000 (less (S1 and S2), gated by (i)); Holm: rank 5, threshold 0.0125, not rejected
- sub-condition i_M2_beats_M1: holds
- sub-condition i_M2_beats_M3: fails
- sub-condition ii_S1_negative: fails
- sub-condition ii_S2_negative: fails
- (i) leave-one-day-out weighted error: M1 0.0171, M2 0.0159, M3 0.0013 (days used 8); M2 fit on all days: {'mu_m': 1.753281322393341, 'lambda_sigma2': 0.016303476162487007, 'a_c': [1.0, 'inf'], 'k_below': 5}
- (ii) slope of zeta_V real minus S1: 0.0453 [0.0287, 0.0559] (day_stock_bootstrap); p_S1 1.0000
- (ii) slope of zeta_V real minus S2: 0.0480 [0.0327, 0.0569] (day_stock_bootstrap); p_S2 1.0000
- passed: no

Pre-registered interpretation: "No non-rejection is reported as evidence against the paper. The measured power at reachable horizons does not support that reading."

### P-R, Q-T1

- range: primary [8,64]; stocks at T2 of each range: {'[4,32]': 107, '[8,64]': 64}
- p = 1.0000 (equivalence (TOST, both ends)); Holm: rank 6, threshold 0.0167, not rejected
- sub-condition T1_within_margin: fails
- sub-condition T2_within_margin: fails
- a*(T1): 0.1271 [-0.1922, 0.4465] (day_jackknife_t); censored: None
- a*(T2): 0.0000 [0.0000, 0.0000] (day_stock_bootstrap); censored: low
- TOST p at T1, T2: 1.0000, 1.0000 (margin 0.35 to 0.65)
- passed: no

Pre-registered interpretation: "No non-rejection is reported as evidence against the paper. The measured power at reachable horizons does not support that reading."

### P-M, Q-T1

- range: lags [1,16], L = 128; stocks at T2 of each range: {'[4,32]': 107, '[8,64]': 64}
- p = 1.0000 (less (both)); Holm: rank 7, threshold 0.0250, not rejected
- sub-condition slope_declines_with_log_q: fails
- sub-condition top_bin_upper_below_0: fails
- beta-hat by size bin (small to large): 0.0300 [0.0100, 0.0500] (day_stock_bootstrap); -0.0300 [-0.0705, -0.0100] (day_stock_bootstrap); 0.0600 [0.0295, 1.0300] (day_stock_bootstrap)
- slope of beta-hat on log q: 0.0031 [-0.0115, 0.2562] (day_stock_bootstrap); p_slope 0.6867
- top bin beta-hat: 0.0600 [0.0295, 1.0300] (day_stock_bootstrap); p_top 1.0000
- stocks used 39, stock-days used 312
- passed: no

Pre-registered interpretation: "No non-rejection is reported as evidence against the paper. The measured power at reachable horizons does not support that reading."

### P-F, Q-T1

- range: primary [8,64]; stocks at T2 of each range: {'[4,32]': 107, '[8,64]': 64}
- p = 1.43e-04 (two-sided); Holm: rank 1, threshold 0.0063, rejected
- sub-condition interval_excludes_0: holds
- slope of zeta_K(a) on a over [0, 1]: -0.0398 [-0.0558, -0.0271] (day_stock_bootstrap)
- zeta_K(a) levels, a = 0, 0.25, 0.5, 0.75, 1 (reported, not tested): 1.0130 [0.9870, 1.0389] (day_jackknife_t); 0.9940 [0.9662, 1.0205] (day_stock_bootstrap); 0.9821 [0.9524, 1.0091] (day_stock_bootstrap); 0.9756 [0.9446, 1.0027] (day_stock_bootstrap); 0.9724 [0.9398, 1.0001] (day_stock_bootstrap)
- passed: yes

Pre-registered interpretation: "A P-F rejection is evidence against the informed-flow regime and nothing more. It does not indicate size-dependent impact, because a plain linear propagator also produces it."

### P-V, Q-T3

- range: primary [8,64]; stocks at T2 of each range: {'[4,32]': 82, '[8,64]': 67}
- p = 0.9700 (less (S1 and S2), gated by (i)); Holm: rank 4, threshold 0.0100, not rejected
- sub-condition i_M2_beats_M1: holds
- sub-condition i_M2_beats_M3: holds
- sub-condition ii_S1_negative: fails
- sub-condition ii_S2_negative: fails
- (i) leave-one-day-out weighted error: M1 0.0267, M2 0.0138, M3 0.0213 (days used 8); M2 fit on all days: {'mu_m': 1.7520478879559644, 'lambda_sigma2': 0.061183879445720464, 'a_c': [1.0, 'inf'], 'k_below': 5}
- (ii) slope of zeta_V real minus S1: 0.0040 [-0.0258, 0.0337] (day_jackknife_t); p_S1 0.6300
- (ii) slope of zeta_V real minus S2: 0.0187 [-0.0094, 0.0467] (day_jackknife_t); p_S2 0.9700
- passed: no

Pre-registered interpretation: "No non-rejection is reported as evidence against the paper. The measured power at reachable horizons does not support that reading."

### P-R, Q-T3

- range: primary [8,64]; stocks at T2 of each range: {'[4,32]': 82, '[8,64]': 67}
- p = 0.4450 (equivalence (TOST, both ends)); Holm: rank 3, threshold 0.0083, not rejected
- sub-condition T1_within_margin: holds
- sub-condition T2_within_margin: fails
- a*(T1): 0.4087 [0.3752, 0.4422] (day_jackknife_t); censored: None
- a*(T2): 0.3577 [0.2994, 0.4096] (day_stock_bootstrap); censored: None
- TOST p at T1, T2: 0.0064, 0.4450 (margin 0.35 to 0.65)
- passed: no

Pre-registered interpretation: "No non-rejection is reported as evidence against the paper. The measured power at reachable horizons does not support that reading."

### P-M, Q-T3

- range: lags [1,16], L = 128; stocks at T2 of each range: {'[4,32]': 82, '[8,64]': 67}
- p = 1.0000 (less (both)); Holm: rank 8, threshold 0.0500, not rejected
- sub-condition slope_declines_with_log_q: fails
- sub-condition top_bin_upper_below_0: fails
- beta-hat by size bin (small to large): 0.0200 [0.0100, 0.6700] (day_stock_bootstrap); 1.0600 [0.2500, 1.4905] (day_stock_bootstrap); -0.0100 [-0.6698, 0.6498] (day_jackknife_t)
- slope of beta-hat on log q: 0.0484 [-0.1478, 0.2446] (day_jackknife_t); p_slope 1.0000
- top bin beta-hat: -0.0100 [-0.6698, 0.6498] (day_jackknife_t); p_top 0.4889
- stocks used 44, stock-days used 352
- passed: no

Pre-registered interpretation: "No non-rejection is reported as evidence against the paper. The measured power at reachable horizons does not support that reading."

### P-F, Q-T3

- range: primary [8,64]; stocks at T2 of each range: {'[4,32]': 82, '[8,64]': 67}
- p = 0.0048 (two-sided); Holm: rank 2, threshold 0.0071, rejected
- sub-condition interval_excludes_0: holds
- slope of zeta_K(a) on a over [0, 1]: -0.0616 [-0.0988, -0.0381] (day_stock_bootstrap)
- zeta_K(a) levels, a = 0, 0.25, 0.5, 0.75, 1 (reported, not tested): 1.0403 [1.0127, 1.0773] (day_stock_bootstrap); 1.0103 [0.9834, 1.0381] (day_stock_bootstrap); 0.9923 [0.9641, 1.0192] (day_stock_bootstrap); 0.9816 [0.9513, 1.0083] (day_stock_bootstrap); 0.9776 [0.9411, 1.0039] (day_stock_bootstrap)
- passed: yes

Pre-registered interpretation: "A P-F rejection is evidence against the informed-flow regime and nothing more. It does not indicate size-dependent impact, because a plain linear propagator also produces it."

## Descriptive (never in Holm)

Q-T2 (its one cell, Q-T2S3, at its ceiling of 1/5) is computed by run_primary and is descriptive (decision of 2 Nov 2026). p-values are before Holm.

### P-V, Q-T2 (descriptive)

- range: primary [8,64]
- p before Holm = 1.0000 (less (S1 and S2), gated by (i))
- sub-condition i_M2_beats_M1: holds
- sub-condition i_M2_beats_M3: fails
- sub-condition ii_S1_negative: fails
- sub-condition ii_S2_negative: fails
- (i) leave-one-day-out weighted error: M1 0.0500, M2 0.0419, M3 0.0230 (days used 8); M2 fit on all days: {'mu_m': 1.6998362995495002, 'lambda_sigma2': 0.04352412220968709, 'a_c': [1.0, 'inf'], 'k_below': 5}
- (ii) slope of zeta_V real minus S1: 0.0354 [-0.0517, 0.1046] (day_stock_bootstrap); p_S1 0.9070
- (ii) slope of zeta_V real minus S2: 0.0588 [-0.0245, 0.1181] (day_stock_bootstrap); p_S2 0.9750

### P-R, Q-T2 (descriptive)

- range: primary [8,64]
- p before Holm = 0.8750 (equivalence (TOST, both ends))
- sub-condition T1_within_margin: fails
- sub-condition T2_within_margin: fails
- a*(T1): 0.3122 [0.2192, 0.3880] (day_stock_bootstrap); censored: None
- a*(T2): 0.2392 [0.0000, 1.0000] (day_stock_bootstrap); censored: None
- TOST p at T1, T2: 0.8750, 0.7450 (margin 0.35 to 0.65)

### P-M, Q-T2 (descriptive)

- not computed: 1 stocks reach L = 128 under their ceilings (need 5)

### P-F, Q-T2 (descriptive)

- range: primary [8,64]
- p before Holm = 0.2600 (two-sided)
- sub-condition interval_excludes_0: fails
- slope of zeta_K(a) on a over [0, 1]: -0.0554 [-0.1409, 0.0282] (day_stock_bootstrap)
- zeta_K(a) levels, a = 0, 0.25, 0.5, 0.75, 1 (reported, not tested): 0.9595 [0.8568, 1.0412] (day_stock_bootstrap); 0.9310 [0.8316, 1.0113] (day_stock_bootstrap); 0.9144 [0.8175, 0.9963] (day_stock_bootstrap); 0.9060 [0.8051, 0.9881] (day_stock_bootstrap); 0.9027 [0.7951, 0.9856] (day_stock_bootstrap)

### P-K (descriptive): vk_diff = slope zeta_V - 2 slope zeta_K, by tick bin and range

- Q-T1 [4,32]: 0.0646 [0.0395, 0.0914] (day_stock_bootstrap)
- Q-T1 [8,64]: 0.0869 [0.0488, 0.1215] (day_stock_bootstrap)
- Q-T2 [4,32]: 0.0114 [-0.1165, 0.1316] (day_stock_bootstrap)
- Q-T2 [8,64]: 0.0475 [-0.2094, 0.2603] (day_stock_bootstrap)
- Q-T3 [4,32]: -0.0129 [-0.0538, 0.0346] (day_stock_bootstrap)
- Q-T3 [8,64]: 0.0457 [-0.0093, 0.1171] (day_stock_bootstrap)

The descriptive cells Q-T1S3, Q-T2S1 and Q-T2S2 were loaded and dropped from the tick-bin units (logged by tick_bin_unit). Their per-cell results are a secondary analysis in the frozen plan (page 7) and are not computed here.

## Stock counts

| cell | role | stocks | stock-days | runs | invalid mids |
|---|---|---|---|---|---|
| Q-T1S1 | confirmatory | 161 | 1288 | 1076988 | 0 |
| Q-T1S2 | confirmatory | 216 | 1728 | 1645242 | 0 |
| Q-T1S3 | descriptive | 252 | 2016 | 965974 | 0 |
| Q-T2S1 | descriptive | 80 | 640 | 1299311 | 0 |
| Q-T2S2 | descriptive | 50 | 400 | 350111 | 0 |
| Q-T2S3 | descriptive | 17 | 136 | 76386 | 0 |
| Q-T3S1 | confirmatory | 50 | 400 | 980014 | 0 |
| Q-T3S2 | confirmatory | 27 | 216 | 368946 | 0 |
| Q-T3S3 | confirmatory | 23 | 184 | 230223 | 0 |

Full tables, each method's p-value and interval, and the evaluability log: data/phase4/tables.md.
