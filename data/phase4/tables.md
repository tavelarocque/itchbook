# Phase 4 confirmatory run: tables

Seed 20261102; P-M design {'n_bins': 3, 'l_max': 16, 'L': 128, 'min_runs_per_bin': 50}; B = 200, surrogate draws 4; git HEAD d126088fba041824950d7ffae4aa75a7a98b2703; wall time 40 s.

## Holm family (Q-T1 and Q-T3, family-wise 0.05)

Family size m = 8.

| rank | test | tick bin | range used | p | Holm threshold | adjusted p | Holm rejects | sub-conditions | passed |
|---|---|---|---|---|---|---|---|---|---|
| 1 | P-F | Q-T1 | [8,64] | 1.43e-04 | 0.0063 | 0.0011 | yes | interval_excludes_0 yes | yes |
| 2 | P-F | Q-T3 | [8,64] | 0.0048 | 0.0071 | 0.0335 | yes | interval_excludes_0 yes | yes |
| 3 | P-R | Q-T3 | [8,64] | 0.4450 | 0.0083 | 1.0000 | no | T1_within_margin yes, T2_within_margin no | no |
| 4 | P-V | Q-T3 | [8,64] | 0.9700 | 0.0100 | 1.0000 | no | i_M2_beats_M1 yes, i_M2_beats_M3 yes, ii_S1_negative no, ii_S2_negative no | no |
| 5 | P-V | Q-T1 | [8,64] | 1.0000 | 0.0125 | 1.0000 | no | i_M2_beats_M1 yes, i_M2_beats_M3 no, ii_S1_negative no, ii_S2_negative no | no |
| 6 | P-R | Q-T1 | [8,64] | 1.0000 | 0.0167 | 1.0000 | no | T1_within_margin no, T2_within_margin no | no |
| 7 | P-M | Q-T1 | lags [1,16], L = 128 | 1.0000 | 0.0250 | 1.0000 | no | slope_declines_with_log_q no, top_bin_upper_below_0 no | no |
| 8 | P-M | Q-T3 | lags [1,16], L = 128 | 1.0000 | 0.0500 | 1.0000 | no | slope_declines_with_log_q no, top_bin_upper_below_0 no | no |

## Pairs that left the family (reach rule, before any p-value)

- none

## Each method's p-value and 0.90 interval (range used)

Columns: p / [lo, hi] per method; the frozen p is the largest of the three. For P-R the frozen TOST p
at an end is the larger of its two rows, and 1 at a censored end (censoring is not applied in these rows).

| tick bin | test | quantity | estimate | null | day-and-stock bootstrap | day jackknife t | whole-stock bootstrap | frozen p |
|---|---|---|---|---|---|---|---|---|
| Q-T1 | P-F | slope | -0.0398 | two-sided 0 | 0.0000 / [-0.0558, -0.0271] | 1.43e-04 / [-0.0499, -0.0297] | 0.0000 / [-0.0502, -0.0300] | 1.43e-04 |
| Q-T1 | P-R | a_star_vs_0.35 (T1) | 0.1271 | greater 0.35 | 1.0000 / [0.0000, 0.1598] | 0.8862 / [-0.1922, 0.4465] | 1.0000 / [0.0000, 0.1444] | 1.0000 |
| Q-T1 | P-R | a_star_vs_0.35 (T2) | 0.0000 | greater 0.35 | 1.0000 / [0.0000, 0.0000] | 1.0000 / [0.0000, 0.0000] | 1.0000 / [0.0000, 0.0000] | 1.0000 |
| Q-T1 | P-R | a_star_vs_0.65 (T1) | 0.1271 | less 0.65 | 0.0000 / [0.0000, 0.1598] | 0.0086 / [-0.1922, 0.4465] | 0.0000 / [0.0000, 0.1444] | 0.0086 |
| Q-T1 | P-R | a_star_vs_0.65 (T2) | 0.0000 | less 0.65 | 0.0000 / [0.0000, 0.0000] | 0.0000 / [0.0000, 0.0000] | 0.0000 / [0.0000, 0.0000] | 0.0000 |
| Q-T1 | P-V | slope_S1 | 0.0453 | less 0 | 1.0000 / [0.0287, 0.0559] | 1.0000 / [0.0356, 0.0550] | 1.0000 / [0.0370, 0.0532] | 1.0000 |
| Q-T1 | P-V | slope_S2 | 0.0480 | less 0 | 1.0000 / [0.0327, 0.0569] | 1.0000 / [0.0386, 0.0573] | 1.0000 / [0.0397, 0.0551] | 1.0000 |
| Q-T1 | P-M | slope | 0.0031 | less 0 | 0.6200 / [-0.0115, 0.2562] | 0.6867 / [-0.0085, 0.0147] | 0.5850 / [-0.0066, 0.0439] | 0.6867 |
| Q-T1 | P-M | top | 0.0600 | less 0 | 1.0000 / [0.0295, 1.0300] | 0.9925 / [0.0246, 0.0954] | 1.0000 / [0.0300, 0.2010] | 1.0000 |
| Q-T2 | P-F | slope | -0.0554 | two-sided 0 | 0.2600 / [-0.1409, 0.0282] | 0.1503 / [-0.1204, 0.0096] | 0.1600 / [-0.1020, 0.0056] | 0.2600 |
| Q-T2 | P-R | a_star_vs_0.35 (T1) | 0.3122 | greater 0.35 | 0.8000 / [0.2192, 0.3880] | 0.8398 / [0.2452, 0.3792] | 0.8750 / [0.2659, 0.3652] | 0.8750 |
| Q-T2 | P-R | a_star_vs_0.35 (T2) | 0.2392 | greater 0.35 | 0.6850 / [0.0000, 1.0000] | 0.7424 / [-0.0670, 0.5455] | 0.7450 / [0.0000, 0.5726] | 0.7450 |
| Q-T2 | P-R | a_star_vs_0.65 (T1) | 0.3122 | less 0.65 | 0.0000 / [0.2192, 0.3880] | 1.44e-05 / [0.2452, 0.3792] | 0.0000 / [0.2659, 0.3652] | 1.44e-05 |
| Q-T2 | P-R | a_star_vs_0.65 (T2) | 0.2392 | less 0.65 | 0.1150 / [0.0000, 1.0000] | 0.0193 / [-0.0670, 0.5455] | 0.0300 / [0.0000, 0.5726] | 0.1150 |
| Q-T2 | P-V | slope_S1 | 0.0354 | less 0 | 0.6850 / [-0.0517, 0.1046] | 0.9070 / [-0.0103, 0.0812] | 0.7800 / [-0.0268, 0.0978] | 0.9070 |
| Q-T2 | P-V | slope_S2 | 0.0588 | less 0 | 0.8750 / [-0.0245, 0.1181] | 0.9748 / [0.0116, 0.1061] | 0.9750 / [0.0092, 0.1120] | 0.9750 |
| Q-T3 | P-F | slope | -0.0616 | two-sided 0 | 0.0000 / [-0.0988, -0.0381] | 0.0048 / [-0.0903, -0.0329] | 0.0000 / [-0.0768, -0.0496] | 0.0048 |
| Q-T3 | P-R | a_star_vs_0.35 (T1) | 0.4087 | greater 0.35 | 0.0000 / [0.3817, 0.4368] | 0.0064 / [0.3752, 0.4422] | 0.0000 / [0.3979, 0.4199] | 0.0064 |
| Q-T3 | P-R | a_star_vs_0.35 (T2) | 0.3577 | greater 0.35 | 0.4450 / [0.2994, 0.4096] | 0.3347 / [0.3250, 0.3904] | 0.3650 / [0.3156, 0.3929] | 0.4450 |
| Q-T3 | P-R | a_star_vs_0.65 (T1) | 0.4087 | less 0.65 | 0.0000 / [0.3817, 0.4368] | 1.35e-06 / [0.3752, 0.4422] | 0.0000 / [0.3979, 0.4199] | 1.35e-06 |
| Q-T3 | P-R | a_star_vs_0.65 (T2) | 0.3577 | less 0.65 | 0.0000 / [0.2994, 0.4096] | 3.07e-07 / [0.3250, 0.3904] | 0.0000 / [0.3156, 0.3929] | 3.07e-07 |
| Q-T3 | P-V | slope_S1 | 0.0040 | less 0 | 0.4850 / [-0.0308, 0.0248] | 0.5961 / [-0.0258, 0.0337] | 0.6300 / [-0.0150, 0.0172] | 0.6300 |
| Q-T3 | P-V | slope_S2 | 0.0187 | less 0 | 0.8300 / [-0.0151, 0.0385] | 0.8762 / [-0.0094, 0.0467] | 0.9700 / [0.0033, 0.0319] | 0.9700 |
| Q-T3 | P-M | slope | 0.0484 | less 0 | 0.9050 / [-0.0594, 0.2508] | 0.6727 / [-0.1478, 0.2446] | 1.0000 / [0.0192, 0.1901] | 1.0000 |
| Q-T3 | P-M | top | -0.0100 | less 0 | 0.4050 / [-0.0200, 0.7205] | 0.4889 / [-0.6698, 0.6498] | 0.3950 / [-0.0100, 0.4600] | 0.4889 |

## Evaluability log

- Q-T1 P-V: evaluable; primary [8,64] ([4,32]: 107 stocks at T2, [8,64]: 64 stocks at T2)
- Q-T1 P-R: evaluable; primary [8,64] ([4,32]: 107 stocks at T2, [8,64]: 64 stocks at T2)
- Q-T1 P-M: evaluable; 39 stocks reach L = 128 under their ceilings (need 5)
- Q-T1 P-F: evaluable; primary [8,64] ([4,32]: 107 stocks at T2, [8,64]: 64 stocks at T2)
- Q-T2 P-V: evaluable; primary [8,64] ([4,32]: 12 stocks at T2, [8,64]: 8 stocks at T2)
- Q-T2 P-R: evaluable; primary [8,64] ([4,32]: 12 stocks at T2, [8,64]: 8 stocks at T2)
- Q-T2 P-M: not evaluable; 1 stocks reach L = 128 under their ceilings (need 5)
- Q-T2 P-F: evaluable; primary [8,64] ([4,32]: 12 stocks at T2, [8,64]: 8 stocks at T2)
- Q-T3 P-V: evaluable; primary [8,64] ([4,32]: 82 stocks at T2, [8,64]: 67 stocks at T2)
- Q-T3 P-R: evaluable; primary [8,64] ([4,32]: 82 stocks at T2, [8,64]: 67 stocks at T2)
- Q-T3 P-M: evaluable; 44 stocks reach L = 128 under their ceilings (need 5)
- Q-T3 P-F: evaluable; primary [8,64] ([4,32]: 82 stocks at T2, [8,64]: 67 stocks at T2)

## Stock counts

| cell | role | stocks | stock-days | runs | invalid mids | min N_s |
|---|---|---|---|---|---|---|
| Q-T1S1 | confirmatory | 161 | 1288 | 1076988 | 0 | 1 |
| Q-T1S2 | confirmatory | 216 | 1728 | 1645242 | 0 | 1 |
| Q-T1S3 | descriptive | 252 | 2016 | 965974 | 0 | 1 |
| Q-T2S1 | descriptive | 80 | 640 | 1299311 | 0 | 1 |
| Q-T2S2 | descriptive | 50 | 400 | 350111 | 0 | 2 |
| Q-T2S3 | descriptive | 17 | 136 | 76386 | 0 | 46 |
| Q-T3S1 | confirmatory | 50 | 400 | 980014 | 0 | 13 |
| Q-T3S2 | confirmatory | 27 | 216 | 368946 | 0 | 37 |
| Q-T3S3 | confirmatory | 23 | 184 | 230223 | 0 | 14 |
