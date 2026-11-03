# Phase 5: secondary and robustness analyses (descriptive)

Secondary and descriptive. Nothing here is confirmatory, nothing here changes any Phase 4 result, and nothing here is a test of the paper. No multiplicity control is applied and no decision is made.

- **Measurement.** Run once on 2026-11-02 17:06:44 -0700 with `tools/est/phase5.py run`, at git HEAD c8187e2a3f66c9b6e4b6c2addca3b2ec5dad0799 (files under tools/ added after the freeze: tools/est/phase4.py, tools/est/phase4_repro.py, tools/est/phase5.py, tools/est/test_phase4.py, tools/est/test_phase5.py; tools/est/phase5.py SHA-256 e7c20eda05f5a16f36d71efeb01b72fa999ef01030c7110bf1b877935ef1f4b2), seed 20261102, B = 200, surrogate draws 4, P-M design {'n_bins': 3, 'l_max': 16, 'L': 128, 'min_runs_per_bin': 50}.
- **Measurement.** Registered plan read from /Users/tavelarocque/Downloads/Frozen_Preregistration_itchbook_2026-11-02.pdf, SHA-256 79c8984ab832c0a9a0d858ac46959277399ce6394106fb526fe9bd48779ccafc (as registered on OSF, osf.io/z4t7x). data/manifest.tsv and data/universe.tsv as registered.
- **Measurement.** Runs files: the SHA-256 of all 8 study days' runs files matched data/phase4/results.json (`runs_npz_sha256`) before anything was computed: 01302019 1207762557d59c7d; 03272019 9f92b33448294436; 07302019 a18e4ad88ca6ad1a; 08302019 c5656b8eed48dfda; 10182019 c1b4bc84bc2f95bb; 10302019 94d67cb87c620983; 12302019 7613e59dc9e808c8; 01302020 a3051e32469eea43 (full hashes in data/phase5/results.json).
- **Measurement.** Confirmation stocks only, loaded with `phase4.build_flows('confirm')`; every symbol checked against the universe's confirmation split; counts per cell as in the plan and as in Phase 4: Q-T1S1 161, Q-T1S2 216, Q-T1S3 252, Q-T2S1 80, Q-T2S2 50, Q-T2S3 17, Q-T3S1 50, Q-T3S2 27, Q-T3S3 23.
- **Measurement.** data/phase4/ and data/phase4_repro/ (5 files) and the tags were hashed before and after the run: unchanged.
- **Measurement.** Wall time 203 s (loading 6 s, anchor 34 s, arms 131 s, cells 32 s).

## Anchor

- **Measurement.** phase5's code path with every setting at the primary (clip R2 through clip3b.apply_rule, centered=False through pooled.panel_moments, seed 20261102) reproduces data/phase4/results.json["result"] bit for bit: 1790 fields compared (floats by float.hex), 0 mismatches.

## Choices made in Phase 5

- **Choice: Cells without a registered ceiling.** Q-T1S3, Q-T2S1 and Q-T2S2 have rho* 'none' in the plan's table (page 3), and the plan does not say how to estimate them; they are reported as "not estimable: no registered ceiling" and nothing is computed for them. No ceiling is invented.
- **Choice: Per-cell generator.** Each per-cell unit runs from its own numpy default_rng(20261102), so a cell's numbers do not depend on which other cells ran or in which order.
- **Choice: Per-cell tests computed.** Only the tests the frozen reach rule (primary3d.evaluability, day lengths and ceilings alone) finds evaluable in the cell are computed, as run_primary does for a tick bin.
- **Choice: Per-cell P-K.** Vk_diff is reported at the range the reach rule selects for the cell (both ranges are in data/phase5/results.json), evaluable when the reach rule gives a range, and with no p-value: the plan fixes no reference value or threshold for P-K.
- **Choice: Per-cell p-values.** Per statistic, the frozen p-value (the largest of the three methods') against the reference value the primary test uses: P-V slopes one-sided below 0, P-R the TOST p at that end against 0.35 and 0.65 (1 at a censored end), P-M slope and top bin one-sided below 0, P-F two-sided against 0; and per test, the frozen decision function's p-value. No multiplicity control is applied.
- **Choice: Q-T2S3 per cell.** Q-T2S3 is Q-T2's only cell, so its per-cell unit holds the same stocks as Phase 4's Q-T2 unit; its numbers here differ from Phase 4's descriptive Q-T2 numbers only through the draws (a fresh generator here, the generator after Q-T1 in Phase 4).
- **Choice: How an arm changes a setting.** The frozen make_unit hardcodes R2 and the frozen moments_raw does not expose centring, so two small wrappers in phase5.py are used for the duration of an arm: primary3d.r2_clip is replaced by clip_sizes, which calls the frozen clip3b.apply_rule(rule) on the same stock-days, and pooled.panel_moments is called with its frozen centered flag. The primary arm runs through the same wrappers (R2 by clip3b.apply_rule) and the anchor shows it bit-identical to Phase 4.
- **Choice: Clip arm scope.** The arm's rule replaces R2 wherever the primary applies R2 (every use of the unit's run sizes: window moments, surrogates, M3, sigma_l^2, the size bins of P-M and of the impact slope), because the plan caps the sizes entering any moment; nothing else changes.
- **Choice: R1 definition.** R1 is clip3b's frozen R1: each run capped at 1 % of the stock-day's summed run volume in this sample (NASDAQ displayed runs), the only R1 in the frozen code; consolidated volume is not in the data.
- **Choice: Centring scope.** Centred moments subtract each stock-day's window means in V, K and S (pooled.panel_moments centered=True) for the real flow and for the surrogates alike; the M3 prediction, sigma_l^2, the impact slope's inputs and the P-M propagator moments are unchanged, since the frozen code exposes centring only for the window moments. P-M therefore differs from the primary only if the draws differ.
- **Choice: Generator in the arms.** Each arm reruns the frozen run_primary over all nine cells from a fresh default_rng(20261102), in the frozen order Q-T1, Q-T2, Q-T3. Q-T1's bootstrap draws come before any surrogate draw and are identical in every arm (checked). A clip rule can change the S2 block lengths (surrogates.choose_L reads the sizes) and so the draws of later units; whether each unit's draws equal the primary's is measured and reported. No unit is re-seeded, which would change the primary path.
- **Choice: Q-T2 in the arms.** Q-T2 is computed inside each arm's run_primary (to keep the frozen generator order) and not reported; the arms' scope is Q-T1 and Q-T3.
- **Choice: Excludes the reference value.** An interval excludes v when v < lo or v > hi; for P-R both 0.35 and 0.65 are shown; P-M's top-bin statistic uses 0 like its slope. The column is descriptive; nothing is decided from it.
- **Choice: Primary column of section C.** The primary values are those of the anchor run, bit-identical to data/phase4/results.json.
- **Choice: Section B ranges.** Both ranges are shown for the range-dependent secondaries, the range used marked; values are read from data/phase4/results.json and shown to 4 decimals.

## A. Per-cell results, all nine NASDAQ-listed cells

**Measurement** (computed here; data/phase5/results.json, key `cells`). Each cell with a registered ceiling is one pooling unit; the reach rule decides which tests are computed. Intervals are 0.90.

### A.1 Cells and tests

| cell | rho* | stocks | test | evaluable (frozen reach rule) | reach rule's reason | range used | stocks at T2: [4,32] | stocks at T2: [8,64] | p (descriptive, no multiplicity control) | sub-conditions |
|---|---|---|---|---|---|---|---|---|---|---|
| Q-T1S1 | 1/5 | 161 | P-V | yes | primary [8,64] ([4,32]: 55 stocks at T2, [8,64]: 41 stocks at T2) | [8,64] | 55 | 41 | 1.0000 | i_M2_beats_M1 holds, i_M2_beats_M3 fails, ii_S1_negative fails, ii_S2_negative fails |
| Q-T1S1 | 1/5 | 161 | P-R | yes | primary [8,64] ([4,32]: 55 stocks at T2, [8,64]: 41 stocks at T2) | [8,64] | 55 | 41 | 1.0000 | T1_within_margin fails, T2_within_margin fails |
| Q-T1S1 | 1/5 | 161 | P-M | yes | 34 stocks reach L = 128 under their ceilings (need 5) | lags [1,16], L = 128 | 55 | 41 | 1.0000 | slope_declines_with_log_q fails, top_bin_upper_below_0 fails |
| Q-T1S1 | 1/5 | 161 | P-F | yes | primary [8,64] ([4,32]: 55 stocks at T2, [8,64]: 41 stocks at T2) | [8,64] | 55 | 41 | 0.0300 | interval_excludes_0 holds |
| Q-T1S1 | 1/5 | 161 | P-K | yes | primary [8,64] ([4,32]: 55 stocks at T2, [8,64]: 41 stocks at T2) | [8,64] | 55 | 41 | n/a (no reference value) | n/a |
| Q-T1S2 | 1/20 | 216 | P-V | yes | primary [8,64] ([4,32]: 52 stocks at T2, [8,64]: 23 stocks at T2) | [8,64] | 52 | 23 | 1.0000 | i_M2_beats_M1 holds, i_M2_beats_M3 fails, ii_S1_negative fails, ii_S2_negative fails |
| Q-T1S2 | 1/20 | 216 | P-R | yes | primary [8,64] ([4,32]: 52 stocks at T2, [8,64]: 23 stocks at T2) | [8,64] | 52 | 23 | 1.0000 | T1_within_margin fails, T2_within_margin fails |
| Q-T1S2 | 1/20 | 216 | P-M | yes | 5 stocks reach L = 128 under their ceilings (need 5) | lags [1,16], L = 128 | 52 | 23 | 0.9349 | slope_declines_with_log_q holds, top_bin_upper_below_0 fails |
| Q-T1S2 | 1/20 | 216 | P-F | yes | primary [8,64] ([4,32]: 52 stocks at T2, [8,64]: 23 stocks at T2) | [8,64] | 52 | 23 | 6.10e-04 | interval_excludes_0 holds |
| Q-T1S2 | 1/20 | 216 | P-K | yes | primary [8,64] ([4,32]: 52 stocks at T2, [8,64]: 23 stocks at T2) | [8,64] | 52 | 23 | n/a (no reference value) | n/a |
| Q-T1S3 | none | 252 | all | n/a | not estimable: no registered ceiling | n/a | n/a | n/a | n/a | n/a |
| Q-T2S1 | none | 80 | all | n/a | not estimable: no registered ceiling | n/a | n/a | n/a | n/a | n/a |
| Q-T2S2 | none | 50 | all | n/a | not estimable: no registered ceiling | n/a | n/a | n/a | n/a | n/a |
| Q-T2S3 | 1/5 | 17 | P-V | yes | primary [8,64] ([4,32]: 12 stocks at T2, [8,64]: 8 stocks at T2) | [8,64] | 12 | 8 | 1.0000 | i_M2_beats_M1 holds, i_M2_beats_M3 fails, ii_S1_negative fails, ii_S2_negative fails |
| Q-T2S3 | 1/5 | 17 | P-R | yes | primary [8,64] ([4,32]: 12 stocks at T2, [8,64]: 8 stocks at T2) | [8,64] | 12 | 8 | 0.9200 | T1_within_margin fails, T2_within_margin fails |
| Q-T2S3 | 1/5 | 17 | P-M | no | 1 stocks reach L = 128 under their ceilings (need 5) | n/a | 12 | 8 | n/a (not computed) | n/a |
| Q-T2S3 | 1/5 | 17 | P-F | yes | primary [8,64] ([4,32]: 12 stocks at T2, [8,64]: 8 stocks at T2) | [8,64] | 12 | 8 | 0.2300 | interval_excludes_0 fails |
| Q-T2S3 | 1/5 | 17 | P-K | yes | primary [8,64] ([4,32]: 12 stocks at T2, [8,64]: 8 stocks at T2) | [8,64] | 12 | 8 | n/a (no reference value) | n/a |
| Q-T3S1 | 1/5 | 50 | P-V | yes | primary [8,64] ([4,32]: 44 stocks at T2, [8,64]: 39 stocks at T2) | [8,64] | 44 | 39 | 1.0000 | i_M2_beats_M1 holds, i_M2_beats_M3 fails, ii_S1_negative fails, ii_S2_negative fails |
| Q-T3S1 | 1/5 | 50 | P-R | yes | primary [8,64] ([4,32]: 44 stocks at T2, [8,64]: 39 stocks at T2) | [8,64] | 44 | 39 | 0.9350 | T1_within_margin holds, T2_within_margin fails |
| Q-T3S1 | 1/5 | 50 | P-M | yes | 31 stocks reach L = 128 under their ceilings (need 5) | lags [1,16], L = 128 | 44 | 39 | 0.9750 | slope_declines_with_log_q fails, top_bin_upper_below_0 fails |
| Q-T3S1 | 1/5 | 50 | P-F | yes | primary [8,64] ([4,32]: 44 stocks at T2, [8,64]: 39 stocks at T2) | [8,64] | 44 | 39 | 0.0279 | interval_excludes_0 holds |
| Q-T3S1 | 1/5 | 50 | P-K | yes | primary [8,64] ([4,32]: 44 stocks at T2, [8,64]: 39 stocks at T2) | [8,64] | 44 | 39 | n/a (no reference value) | n/a |
| Q-T3S2 | 1/5 | 27 | P-V | yes | primary [8,64] ([4,32]: 25 stocks at T2, [8,64]: 23 stocks at T2) | [8,64] | 25 | 23 | 0.9000 | i_M2_beats_M1 holds, i_M2_beats_M3 holds, ii_S1_negative fails, ii_S2_negative fails |
| Q-T3S2 | 1/5 | 27 | P-R | yes | primary [8,64] ([4,32]: 25 stocks at T2, [8,64]: 23 stocks at T2) | [8,64] | 25 | 23 | 0.4000 | T1_within_margin holds, T2_within_margin fails |
| Q-T3S2 | 1/5 | 27 | P-M | yes | 12 stocks reach L = 128 under their ceilings (need 5) | lags [1,16], L = 128 | 25 | 23 | 0.8313 | slope_declines_with_log_q fails, top_bin_upper_below_0 fails |
| Q-T3S2 | 1/5 | 27 | P-F | yes | primary [8,64] ([4,32]: 25 stocks at T2, [8,64]: 23 stocks at T2) | [8,64] | 25 | 23 | 0.0452 | interval_excludes_0 holds |
| Q-T3S2 | 1/5 | 27 | P-K | yes | primary [8,64] ([4,32]: 25 stocks at T2, [8,64]: 23 stocks at T2) | [8,64] | 25 | 23 | n/a (no reference value) | n/a |
| Q-T3S3 | 1/10 | 23 | P-V | yes | primary [8,64] ([4,32]: 13 stocks at T2, [8,64]: 5 stocks at T2) | [8,64] | 13 | 5 | 0.2100 | i_M2_beats_M1 holds, i_M2_beats_M3 holds, ii_S1_negative fails, ii_S2_negative fails |
| Q-T3S3 | 1/10 | 23 | P-R | yes | primary [8,64] ([4,32]: 13 stocks at T2, [8,64]: 5 stocks at T2) | [8,64] | 13 | 5 | 0.1307 | T1_within_margin holds, T2_within_margin fails |
| Q-T3S3 | 1/10 | 23 | P-M | no | 1 stocks reach L = 128 under their ceilings (need 5) | n/a | 13 | 5 | n/a (not computed) | n/a |
| Q-T3S3 | 1/10 | 23 | P-F | yes | primary [8,64] ([4,32]: 13 stocks at T2, [8,64]: 5 stocks at T2) | [8,64] | 13 | 5 | 0.0603 | interval_excludes_0 holds |
| Q-T3S3 | 1/10 | 23 | P-K | yes | primary [8,64] ([4,32]: 13 stocks at T2, [8,64]: 5 stocks at T2) | [8,64] | 13 | 5 | n/a (no reference value) | n/a |

### A.2 Statistics, each method's interval

**Measurement** (computed here). Widest-of-three 0.90 interval and the method that gave it; then each method's 0.90 interval. P-R rows at a censored end show the end's a.

| cell | test | statistic | range | estimate | widest 0.90 interval (method) | day-and-stock bootstrap | day jackknife t | whole-stock bootstrap | p (descriptive, no multiplicity control) |
|---|---|---|---|---|---|---|---|---|---|
| Q-T1S1 | P-V | slope in a of zeta_V real minus S1 | [8,64] | 0.0567 | [0.0329, 0.0782] (day_stock_bootstrap) | [0.0329, 0.0782] | [0.0389, 0.0746] | [0.0444, 0.0695] | 1.0000 |
| Q-T1S1 | P-V | slope in a of zeta_V real minus S2 | [8,64] | 0.0602 | [0.0388, 0.0764] (day_stock_bootstrap) | [0.0388, 0.0764] | [0.0459, 0.0744] | [0.0492, 0.0734] | 1.0000 |
| Q-T1S1 | P-R | a*(T1) | [8,64] | 0.1624 | [0.0000, 0.2179] (day_stock_bootstrap) | [0.0000, 0.2179] | [0.1080, 0.2169] | [0.1321, 0.1938] | 1.0000 |
| Q-T1S1 | P-R | a*(T2), censored low | [8,64] | 0.0000 | [0.0000, 0.0000] (day_stock_bootstrap) | [0.0000, 0.0000] | [0.0000, 0.0000] | [0.0000, 0.0000] | 1.0000 |
| Q-T1S1 | P-M | slope of beta-hat on log q | lags [1,16], L = 128 | 0.0132 | [-0.0059, 0.3864] (day_stock_bootstrap) | [-0.0059, 0.3864] | [-0.0212, 0.0476] | [9.08e-04, 0.3710] | 0.9600 |
| Q-T1S1 | P-M | top size bin beta-hat | lags [1,16], L = 128 | 0.0900 | [0.0300, 1.4800] (day_stock_bootstrap) | [0.0300, 1.4800] | [-0.0249, 0.2049] | [0.0500, 1.4500] | 1.0000 |
| Q-T1S1 | P-F | slope of zeta_K(a) on a over [0, 1] | [8,64] | -0.0249 | [-0.0442, -0.0033] (day_stock_bootstrap) | [-0.0442, -0.0033] | [-0.0317, -0.0181] | [-0.0385, -0.0096] | 0.0300 |
| Q-T1S1 | P-K | vk_diff = slope zeta_V - 2 slope zeta_K | [8,64] | 0.0761 | [0.0319, 0.1199] (day_stock_bootstrap) | [0.0319, 0.1199] | [0.0438, 0.1084] | [0.0476, 0.1060] | n/a (no reference value) |
| Q-T1S2 | P-V | slope in a of zeta_V real minus S1 | [8,64] | 0.0305 | [0.0070, 0.0460] (day_stock_bootstrap) | [0.0070, 0.0460] | [0.0152, 0.0458] | [0.0192, 0.0421] | 1.0000 |
| Q-T1S2 | P-V | slope in a of zeta_V real minus S2 | [8,64] | 0.0360 | [0.0170, 0.0500] (day_stock_bootstrap) | [0.0170, 0.0500] | [0.0241, 0.0478] | [0.0266, 0.0467] | 1.0000 |
| Q-T1S2 | P-R | a*(T1), censored low | [8,64] | 0.0000 | [0.0000, 0.1436] (day_stock_bootstrap) | [0.0000, 0.1436] | [0.0000, 0.0000] | [0.0000, 0.0000] | 1.0000 |
| Q-T1S2 | P-R | a*(T2), censored low | [8,64] | 0.0000 | [0.0000, 0.0000] (day_stock_bootstrap) | [0.0000, 0.0000] | [0.0000, 0.0000] | [0.0000, 0.0000] | 1.0000 |
| Q-T1S2 | P-M | slope of beta-hat on log q | lags [1,16], L = 128 | -0.0225 | [-0.2852, -0.0062] (day_stock_bootstrap) | [-0.2852, -0.0062] | [-0.0330, -0.0119] | [-0.2840, -0.0083] | 0.0254 |
| Q-T1S2 | P-M | top size bin beta-hat | lags [1,16], L = 128 | 0.0300 | [-1.0000, 0.1000] (day_stock_bootstrap) | [-1.0000, 0.1000] | [-0.0032, 0.0632] | [-0.9700, 0.1100] | 0.9349 |
| Q-T1S2 | P-F | slope of zeta_K(a) on a over [0, 1] | [8,64] | -0.0557 | [-0.0813, -0.0379] (day_stock_bootstrap) | [-0.0813, -0.0379] | [-0.0737, -0.0378] | [-0.0711, -0.0432] | 6.10e-04 |
| Q-T1S2 | P-K | vk_diff = slope zeta_V - 2 slope zeta_K | [8,64] | 0.0943 | [0.0447, 0.1554] (day_stock_bootstrap) | [0.0447, 0.1554] | [0.0397, 0.1489] | [0.0682, 0.1272] | n/a (no reference value) |
| Q-T1S3 | all | not estimable: no registered ceiling | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| Q-T2S1 | all | not estimable: no registered ceiling | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| Q-T2S2 | all | not estimable: no registered ceiling | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| Q-T2S3 | P-V | slope in a of zeta_V real minus S1 | [8,64] | 0.0390 | [-0.0499, 0.1085] (day_stock_bootstrap) | [-0.0499, 0.1085] | [-0.0102, 0.0881] | [-0.0247, 0.1057] | 0.9117 |
| Q-T2S3 | P-V | slope in a of zeta_V real minus S2 | [8,64] | 0.0467 | [-0.0257, 0.1123] (day_stock_bootstrap) | [-0.0257, 0.1123] | [0.0014, 0.0921] | [-0.0037, 0.0941] | 0.9540 |
| Q-T2S3 | P-R | a*(T1) | [8,64] | 0.3122 | [0.2216, 0.3845] (day_stock_bootstrap) | [0.2216, 0.3845] | [0.2452, 0.3792] | [0.2692, 0.3581] | 0.9200 |
| Q-T2S3 | P-R | a*(T2) | [8,64] | 0.2392 | [0.0000, 0.7991] (day_stock_bootstrap) | [0.0000, 0.7991] | [-0.0670, 0.5455] | [0.0000, 0.5417] | 0.7500 |
| Q-T2S3 | P-M | not computed: 1 stocks reach L = 128 under their ceilings (need 5) | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| Q-T2S3 | P-F | slope of zeta_K(a) on a over [0, 1] | [8,64] | -0.0554 | [-0.1347, 0.0156] (day_stock_bootstrap) | [-0.1347, 0.0156] | [-0.1204, 0.0096] | [-0.1007, 0.0109] | 0.2300 |
| Q-T2S3 | P-K | vk_diff = slope zeta_V - 2 slope zeta_K | [8,64] | 0.0475 | [-0.1930, 0.2557] (day_stock_bootstrap) | [-0.1930, 0.2557] | [-0.1270, 0.2219] | [-0.0715, 0.1863] | n/a (no reference value) |
| Q-T3S1 | P-V | slope in a of zeta_V real minus S1 | [8,64] | 0.0097 | [-0.0211, 0.0364] (day_stock_bootstrap) | [-0.0211, 0.0364] | [-0.0183, 0.0377] | [-0.0087, 0.0240] | 0.8400 |
| Q-T3S1 | P-V | slope in a of zeta_V real minus S2 | [8,64] | 0.0206 | [-0.0080, 0.0417] (day_stock_bootstrap) | [-0.0080, 0.0417] | [-0.0018, 0.0429] | [0.0053, 0.0353] | 0.9800 |
| Q-T3S1 | P-R | a*(T1) | [8,64] | 0.3845 | [0.3518, 0.4136] (day_stock_bootstrap) | [0.3518, 0.4136] | [0.3547, 0.4144] | [0.3685, 0.3996] | 0.0350 |
| Q-T3S1 | P-R | a*(T2) | [8,64] | 0.3164 | [0.2282, 0.3741] (day_stock_bootstrap) | [0.2282, 0.3741] | [0.2789, 0.3540] | [0.2636, 0.3532] | 0.9350 |
| Q-T3S1 | P-M | slope of beta-hat on log q | lags [1,16], L = 128 | 0.1071 | [-0.1034, 0.3177] (day_jackknife_t) | [-0.0153, 0.3566] | [-0.1034, 0.3177] | [0.0056, 0.2259] | 0.9750 |
| Q-T3S1 | P-M | top size bin beta-hat | lags [1,16], L = 128 | 0.2300 | [-0.5082, 0.9682] (day_jackknife_t) | [-0.0200, 0.9225] | [-0.5082, 0.9682] | [-0.0100, 0.5810] | 0.7132 |
| Q-T3S1 | P-F | slope of zeta_K(a) on a over [0, 1] | [8,64] | -0.0540 | [-0.0959, -0.0206] (day_stock_bootstrap) | [-0.0959, -0.0206] | [-0.0911, -0.0170] | [-0.0727, -0.0369] | 0.0279 |
| Q-T3S1 | P-K | vk_diff = slope zeta_V - 2 slope zeta_K | [8,64] | 0.0707 | [0.0017, 0.1597] (day_stock_bootstrap) | [0.0017, 0.1597] | [-0.0042, 0.1457] | [0.0394, 0.1061] | n/a (no reference value) |
| Q-T3S2 | P-V | slope in a of zeta_V real minus S1 | [8,64] | 0.0114 | [-0.0618, 0.0582] (day_stock_bootstrap) | [-0.0618, 0.0582] | [-0.0408, 0.0636] | [-0.0176, 0.0427] | 0.7200 |
| Q-T3S2 | P-V | slope in a of zeta_V real minus S2 | [8,64] | 0.0228 | [-0.0367, 0.0659] (day_stock_bootstrap) | [-0.0367, 0.0659] | [-0.0225, 0.0682] | [-0.0084, 0.0526] | 0.9000 |
| Q-T3S2 | P-R | a*(T1) | [8,64] | 0.4004 | [0.3583, 0.4426] (day_jackknife_t) | [0.3540, 0.4358] | [0.3583, 0.4426] | [0.3727, 0.4231] | 0.0289 |
| Q-T3S2 | P-R | a*(T2) | [8,64] | 0.3667 | [0.2403, 0.4364] (day_stock_bootstrap) | [0.2403, 0.4364] | [0.3036, 0.4297] | [0.2964, 0.4144] | 0.4000 |
| Q-T3S2 | P-M | slope of beta-hat on log q | lags [1,16], L = 128 | 0.0322 | [-0.4493, 0.2904] (day_stock_bootstrap) | [-0.4493, 0.2904] | [-0.0271, 0.0916] | [-0.4061, 0.0596] | 0.8313 |
| Q-T3S2 | P-M | top size bin beta-hat | lags [1,16], L = 128 | -0.0200 | [-0.0400, 0.8010] (day_stock_bootstrap) | [-0.0400, 0.8010] | [-0.0500, 0.0100] | [-0.0300, -0.0100] | 0.1900 |
| Q-T3S2 | P-F | slope of zeta_K(a) on a over [0, 1] | [8,64] | -0.0577 | [-0.1171, -0.0173] (day_stock_bootstrap) | [-0.1171, -0.0173] | [-0.1026, -0.0128] | [-0.0837, -0.0342] | 0.0452 |
| Q-T3S2 | P-K | vk_diff = slope zeta_V - 2 slope zeta_K | [8,64] | 0.0476 | [-0.0201, 0.1430] (day_stock_bootstrap) | [-0.0201, 0.1430] | [-0.0078, 0.1030] | [-0.0085, 0.1105] | n/a (no reference value) |
| Q-T3S3 | P-V | slope in a of zeta_V real minus S1 | [8,64] | -0.0525 | [-0.1319, 0.0292] (day_stock_bootstrap) | [-0.1319, 0.0292] | [-0.0800, -0.0249] | [-0.1317, 0.0175] | 0.1407 |
| Q-T3S3 | P-V | slope in a of zeta_V real minus S2 | [8,64] | -0.0311 | [-0.1058, 0.0322] (day_stock_bootstrap) | [-0.1058, 0.0322] | [-0.0493, -0.0130] | [-0.1084, 0.0274] | 0.2100 |
| Q-T3S3 | P-R | a*(T1) | [8,64] | 0.4460 | [0.4077, 0.4848] (day_stock_bootstrap) | [0.4077, 0.4848] | [0.4084, 0.4836] | [0.4242, 0.4721] | 9.38e-04 |
| Q-T3S3 | P-R | a*(T2) | [8,64] | 0.4801 | [0.2729, 0.6145] (day_stock_bootstrap) | [0.2729, 0.6145] | [0.4209, 0.5392] | [0.3683, 0.5729] | 0.1307 |
| Q-T3S3 | P-M | not computed: 1 stocks reach L = 128 under their ceilings (need 5) | n/a | n/a | n/a | n/a | n/a | n/a | n/a |
| Q-T3S3 | P-F | slope of zeta_K(a) on a over [0, 1] | [8,64] | -0.0801 | [-0.1633, -0.0106] (day_stock_bootstrap) | [-0.1633, -0.0106] | [-0.1106, -0.0495] | [-0.1428, -0.0372] | 0.0603 |
| Q-T3S3 | P-K | vk_diff = slope zeta_V - 2 slope zeta_K | [8,64] | -0.0625 | [-0.2139, 0.0898] (day_stock_bootstrap) | [-0.2139, 0.0898] | [-0.1477, 0.0227] | [-0.1779, 0.0279] | n/a (no reference value) |

## B. P-R secondaries per tick bin (from Phase 4)

**Measurement** (computed in Phase 4 by primary3d.run_unit and stored in data/phase4/results.json; read here, not recomputed; each value's JSON path is given). Intervals are 0.90, widest of three.

- **Measurement** (Q-T1, `result.by_tick_bin.Q-T1.tests.P-R.secondary.psi_hat`): psi-hat 0.0535 [0.0464, 0.0694] (day_stock_bootstrap) over 330 stocks with a slope. Note: psi-hat reads 0.63 to 0.70 at psi = 0.5 and 0.35 to 0.41 at psi = 0.3 on real flow (Phase 3d-2); not corrected
- **Measurement** (Q-T2, `result.by_tick_bin.Q-T2.tests.P-R.secondary.psi_hat`): psi-hat 0.1305 [0.0786, 0.1760] (day_stock_bootstrap) over 17 stocks with a slope. Note: psi-hat reads 0.63 to 0.70 at psi = 0.5 and 0.35 to 0.41 at psi = 0.3 on real flow (Phase 3d-2); not corrected
- **Measurement** (Q-T3, `result.by_tick_bin.Q-T3.tests.P-R.secondary.psi_hat`): psi-hat 0.1794 [0.1566, 0.2022] (day_jackknife_t) over 99 stocks with a slope. Note: psi-hat reads 0.63 to 0.70 at psi = 0.5 and 0.35 to 0.41 at psi = 0.3 on real flow (Phase 3d-2); not corrected

| tick bin | quantity | range | estimate | widest 0.90 interval (method) | p (descriptive, no multiplicity control) | JSON path |
|---|---|---|---|---|---|---|
| Q-T1 | R_{1/2}/R_0 at T = 8 | [4,32] | 0.9691 | [0.9597, 0.9784] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T1.tests.P-R.secondary.Rratio.[4,32]` |
| Q-T1 | R_{1/2}/R_0 minus e^{sigma_l^2/8} at T = 8 | [4,32] | -0.3665 | [-0.3928, -0.3403] (day_jackknife_t) | 2.84e-08 (two-sided) | `result.by_tick_bin.Q-T1.tests.P-R.secondary.Rratio_minus_exp_sig2_over_8.[4,32]` |
| Q-T1 | da*/dlog T | [4,32] | -0.0829 | [-0.1413, -0.0245] (day_jackknife_t) | 0.2100 (two-sided) | `result.by_tick_bin.Q-T1.tests.P-R.secondary.dastar_dlogT.[4,32]` |
| Q-T1 | R_{1/2}/R_0 at T = 16 | [8,64] (range used) | 0.9533 | [0.9404, 0.9662] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T1.tests.P-R.secondary.Rratio.[8,64]` |
| Q-T1 | R_{1/2}/R_0 minus e^{sigma_l^2/8} at T = 16 | [8,64] (range used) | -0.3823 | [-0.4138, -0.3508] (day_jackknife_t) | 7.43e-08 (two-sided) | `result.by_tick_bin.Q-T1.tests.P-R.secondary.Rratio_minus_exp_sig2_over_8.[8,64]` |
| Q-T1 | da*/dlog T | [8,64] (range used) | -0.0550 | [-0.1932, 0.0832] (day_jackknife_t) | 0.9000 (two-sided) | `result.by_tick_bin.Q-T1.tests.P-R.secondary.dastar_dlogT.[8,64]` |
| Q-T1 | sigma_l^2 (pooled variance of log run size) | n/a | 2.3151 | [2.1986, 2.4316] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T1.tests.P-R.secondary.sigma_l2` |
| Q-T1 | a*(T1) | [8,64] (range used) | 0.1271 | [-0.1922, 0.4465] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T1.tests.P-R.statistic.a_star_T1_T2` |
| Q-T1 | a*(T2), censored low | [8,64] (range used) | 0.0000 | [0.0000, 0.0000] (day_stock_bootstrap) | n/a | `result.by_tick_bin.Q-T1.tests.P-R.statistic.a_star_T1_T2` |
| Q-T2 | R_{1/2}/R_0 at T = 8 | [4,32] | 1.0302 | [0.9885, 1.0727] (day_stock_bootstrap) | n/a | `result.by_tick_bin.Q-T2.tests.P-R.secondary.Rratio.[4,32]` |
| Q-T2 | R_{1/2}/R_0 minus e^{sigma_l^2/8} at T = 8 | [4,32] | -0.4018 | [-0.4710, -0.3326] (day_jackknife_t) | 1.13e-05 (two-sided) | `result.by_tick_bin.Q-T2.tests.P-R.secondary.Rratio_minus_exp_sig2_over_8.[4,32]` |
| Q-T2 | da*/dlog T | [4,32] | -0.0356 | [-0.1291, 0.0482] (day_stock_bootstrap) | 0.4900 (two-sided) | `result.by_tick_bin.Q-T2.tests.P-R.secondary.dastar_dlogT.[4,32]` |
| Q-T2 | R_{1/2}/R_0 at T = 16 | [8,64] (range used) | 1.0315 | [0.9832, 1.0885] (day_stock_bootstrap) | n/a | `result.by_tick_bin.Q-T2.tests.P-R.secondary.Rratio.[8,64]` |
| Q-T2 | R_{1/2}/R_0 minus e^{sigma_l^2/8} at T = 16 | [8,64] (range used) | -0.4006 | [-0.4778, -0.3233] (day_jackknife_t) | 2.41e-05 (two-sided) | `result.by_tick_bin.Q-T2.tests.P-R.secondary.Rratio_minus_exp_sig2_over_8.[8,64]` |
| Q-T2 | da*/dlog T | [8,64] (range used) | -0.0474 | [-0.1753, 0.2667] (day_stock_bootstrap) | 0.6500 (two-sided) | `result.by_tick_bin.Q-T2.tests.P-R.secondary.dastar_dlogT.[8,64]` |
| Q-T2 | sigma_l^2 (pooled variance of log run size) | n/a | 2.8727 | [2.4765, 3.2689] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T2.tests.P-R.secondary.sigma_l2` |
| Q-T2 | a*(T1) | [8,64] (range used) | 0.3122 | [0.2192, 0.3880] (day_stock_bootstrap) | n/a | `result.by_tick_bin.Q-T2.tests.P-R.statistic.a_star_T1_T2` |
| Q-T2 | a*(T2) | [8,64] (range used) | 0.2392 | [0.0000, 1.0000] (day_stock_bootstrap) | n/a | `result.by_tick_bin.Q-T2.tests.P-R.statistic.a_star_T1_T2` |
| Q-T3 | R_{1/2}/R_0 at T = 8 | [4,32] | 1.0835 | [1.0600, 1.1070] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T3.tests.P-R.secondary.Rratio.[4,32]` |
| Q-T3 | R_{1/2}/R_0 minus e^{sigma_l^2/8} at T = 8 | [4,32] | -0.2930 | [-0.3235, -0.2624] (day_jackknife_t) | 3.79e-07 (two-sided) | `result.by_tick_bin.Q-T3.tests.P-R.secondary.Rratio_minus_exp_sig2_over_8.[4,32]` |
| Q-T3 | da*/dlog T | [4,32] | 0.0111 | [-0.0062, 0.0250] (day_stock_bootstrap) | 0.3400 (two-sided) | `result.by_tick_bin.Q-T3.tests.P-R.secondary.dastar_dlogT.[4,32]` |
| Q-T3 | R_{1/2}/R_0 at T = 16 | [8,64] (range used) | 1.0871 | [1.0667, 1.1074] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T3.tests.P-R.secondary.Rratio.[8,64]` |
| Q-T3 | R_{1/2}/R_0 minus e^{sigma_l^2/8} at T = 16 | [8,64] (range used) | -0.2894 | [-0.3249, -0.2540] (day_jackknife_t) | 1.14e-06 (two-sided) | `result.by_tick_bin.Q-T3.tests.P-R.secondary.Rratio_minus_exp_sig2_over_8.[8,64]` |
| Q-T3 | da*/dlog T | [8,64] (range used) | -0.0246 | [-0.0533, -0.0022] (day_stock_bootstrap) | 0.0700 (two-sided) | `result.by_tick_bin.Q-T3.tests.P-R.secondary.dastar_dlogT.[8,64]` |
| Q-T3 | sigma_l^2 (pooled variance of log run size) | n/a | 2.5562 | [2.2491, 2.8634] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T3.tests.P-R.secondary.sigma_l2` |
| Q-T3 | a*(T1) | [8,64] (range used) | 0.4087 | [0.3752, 0.4422] (day_jackknife_t) | n/a | `result.by_tick_bin.Q-T3.tests.P-R.statistic.a_star_T1_T2` |
| Q-T3 | a*(T2) | [8,64] (range used) | 0.3577 | [0.2994, 0.4096] (day_stock_bootstrap) | n/a | `result.by_tick_bin.Q-T3.tests.P-R.statistic.a_star_T1_T2` |

### a* in Q-T1, with the secondaries beside it

**Measurement** (Phase 4 values, read from data/phase4/results.json). In Q-T1 over the range used [8,64], a*(T = 8) is 0.1271 [-0.1922, 0.4465], and a*(T = 64) is 0.0000 [0.0000, 0.0000], censored low (the grid maximum of pooled R_a(64) is at the lowest a of the grid). Beside it: R_{1/2}/R_0 is 0.9533 [0.9404, 0.9662] at T = 16 (below 1, its interval excludes 1) and 0.9691 [0.9597, 0.9784] at T = 8 (below 1, its interval excludes 1); da*/dlog T over [8,64] is -0.0550 [-0.1932, 0.0832], its interval includes 0; psi-hat is 0.0535 [0.0464, 0.0694]. Note: psi-hat reads 0.63 to 0.70 at psi = 0.5 and 0.35 to 0.41 at psi = 0.3 on real flow (Phase 3d-2); not corrected. These values are placed side by side as stored; no causal reading is drawn from them.

## C. Robustness arms in Q-T1 and Q-T3

**Measurement** (computed here; data/phase5/results.json, key `arms`). Each arm changes one setting from the primary (clip R2, uncentred moments), everything else at the primary setting, seed 20261102. Estimates at the range each test used (the reach rule does not depend on the arm). The excludes columns are descriptive; nothing is decided from them.

| tick bin | test | statistic | reference | primary: estimate [0.90 interval] (method) | primary interval excludes reference | arm | arm: estimate [0.90 interval] (method) | arm interval excludes reference |
|---|---|---|---|---|---|---|---|---|
| Q-T1 | P-V | slope in a of zeta_V real minus S1 | 0 | 0.0453 [0.0287, 0.0559] (day_stock_bootstrap) | yes | clip R1 | 0.0483 [0.0332, 0.0581] (day_stock_bootstrap) | yes |
| Q-T1 | P-V | slope in a of zeta_V real minus S1 | 0 | 0.0453 [0.0287, 0.0559] (day_stock_bootstrap) | yes | clip R3 | 0.0406 [0.0263, 0.0505] (day_stock_bootstrap) | yes |
| Q-T1 | P-V | slope in a of zeta_V real minus S1 | 0 | 0.0453 [0.0287, 0.0559] (day_stock_bootstrap) | yes | no clip | 0.0532 [0.0357, 0.0629] (day_stock_bootstrap) | yes |
| Q-T1 | P-V | slope in a of zeta_V real minus S1 | 0 | 0.0453 [0.0287, 0.0559] (day_stock_bootstrap) | yes | centred moments | 0.0177 [0.0012, 0.0284] (day_stock_bootstrap) | yes |
| Q-T1 | P-V | slope in a of zeta_V real minus S2 | 0 | 0.0480 [0.0327, 0.0569] (day_stock_bootstrap) | yes | clip R1 | 0.0516 [0.0374, 0.0608] (day_stock_bootstrap) | yes |
| Q-T1 | P-V | slope in a of zeta_V real minus S2 | 0 | 0.0480 [0.0327, 0.0569] (day_stock_bootstrap) | yes | clip R3 | 0.0451 [0.0313, 0.0539] (day_stock_bootstrap) | yes |
| Q-T1 | P-V | slope in a of zeta_V real minus S2 | 0 | 0.0480 [0.0327, 0.0569] (day_stock_bootstrap) | yes | no clip | 0.0532 [0.0384, 0.0615] (day_stock_bootstrap) | yes |
| Q-T1 | P-V | slope in a of zeta_V real minus S2 | 0 | 0.0480 [0.0327, 0.0569] (day_stock_bootstrap) | yes | centred moments | 0.0176 [0.0023, 0.0292] (day_stock_bootstrap) | yes |
| Q-T1 | P-R | a*(T1) | 0.35 and 0.65 | 0.1271 [-0.1922, 0.4465] (day_jackknife_t) | 0.35: no; 0.65: yes | clip R1 | 0.0000 [-0.3234, 0.3234] (day_jackknife_t), censored low | 0.35: yes; 0.65: yes |
| Q-T1 | P-R | a*(T1) | 0.35 and 0.65 | 0.1271 [-0.1922, 0.4465] (day_jackknife_t) | 0.35: no; 0.65: yes | clip R3 | 0.1263 [-0.2056, 0.4581] (day_jackknife_t) | 0.35: no; 0.65: yes |
| Q-T1 | P-R | a*(T1) | 0.35 and 0.65 | 0.1271 [-0.1922, 0.4465] (day_jackknife_t) | 0.35: no; 0.65: yes | no clip | 0.1271 [-0.1919, 0.4461] (day_jackknife_t) | 0.35: no; 0.65: yes |
| Q-T1 | P-R | a*(T1) | 0.35 and 0.65 | 0.1271 [-0.1922, 0.4465] (day_jackknife_t) | 0.35: no; 0.65: yes | centred moments | 0.1347 [-0.0942, 0.3637] (day_jackknife_t) | 0.35: no; 0.65: yes |
| Q-T1 | P-R | a*(T2) | 0.35 and 0.65 | 0.0000 [0.0000, 0.0000] (day_stock_bootstrap), censored low | 0.35: yes; 0.65: yes | clip R1 | 0.0000 [0.0000, 0.0000] (day_stock_bootstrap), censored low | 0.35: yes; 0.65: yes |
| Q-T1 | P-R | a*(T2) | 0.35 and 0.65 | 0.0000 [0.0000, 0.0000] (day_stock_bootstrap), censored low | 0.35: yes; 0.65: yes | clip R3 | 0.0000 [0.0000, 0.0000] (day_stock_bootstrap), censored low | 0.35: yes; 0.65: yes |
| Q-T1 | P-R | a*(T2) | 0.35 and 0.65 | 0.0000 [0.0000, 0.0000] (day_stock_bootstrap), censored low | 0.35: yes; 0.65: yes | no clip | 0.0000 [0.0000, 0.0000] (day_stock_bootstrap), censored low | 0.35: yes; 0.65: yes |
| Q-T1 | P-R | a*(T2) | 0.35 and 0.65 | 0.0000 [0.0000, 0.0000] (day_stock_bootstrap), censored low | 0.35: yes; 0.65: yes | centred moments | 0.0000 [0.0000, 0.0000] (day_stock_bootstrap), censored low | 0.35: yes; 0.65: yes |
| Q-T1 | P-M | slope of beta-hat on log q | 0 | 0.0031 [-0.0115, 0.2562] (day_stock_bootstrap) | no | clip R1 | 0.0031 [-0.0115, 0.2562] (day_stock_bootstrap) | no |
| Q-T1 | P-M | slope of beta-hat on log q | 0 | 0.0031 [-0.0115, 0.2562] (day_stock_bootstrap) | no | clip R3 | 0.0030 [-0.0116, 0.2562] (day_stock_bootstrap) | no |
| Q-T1 | P-M | slope of beta-hat on log q | 0 | 0.0031 [-0.0115, 0.2562] (day_stock_bootstrap) | no | no clip | 0.0031 [-0.0115, 0.2562] (day_stock_bootstrap) | no |
| Q-T1 | P-M | slope of beta-hat on log q | 0 | 0.0031 [-0.0115, 0.2562] (day_stock_bootstrap) | no | centred moments | 0.0031 [-0.0115, 0.2562] (day_stock_bootstrap) | no |
| Q-T1 | P-M | top size bin beta-hat | 0 | 0.0600 [0.0295, 1.0300] (day_stock_bootstrap) | yes | clip R1 | 0.0600 [0.0295, 1.0300] (day_stock_bootstrap) | yes |
| Q-T1 | P-M | top size bin beta-hat | 0 | 0.0600 [0.0295, 1.0300] (day_stock_bootstrap) | yes | clip R3 | 0.0600 [0.0295, 1.0300] (day_stock_bootstrap) | yes |
| Q-T1 | P-M | top size bin beta-hat | 0 | 0.0600 [0.0295, 1.0300] (day_stock_bootstrap) | yes | no clip | 0.0600 [0.0295, 1.0300] (day_stock_bootstrap) | yes |
| Q-T1 | P-M | top size bin beta-hat | 0 | 0.0600 [0.0295, 1.0300] (day_stock_bootstrap) | yes | centred moments | 0.0600 [0.0295, 1.0300] (day_stock_bootstrap) | yes |
| Q-T1 | P-F | slope of zeta_K(a) on a over [0, 1] | 0 | -0.0398 [-0.0558, -0.0271] (day_stock_bootstrap) | yes | clip R1 | -0.0411 [-0.0574, -0.0296] (day_stock_bootstrap) | yes |
| Q-T1 | P-F | slope of zeta_K(a) on a over [0, 1] | 0 | -0.0398 [-0.0558, -0.0271] (day_stock_bootstrap) | yes | clip R3 | -0.0453 [-0.0595, -0.0331] (day_stock_bootstrap) | yes |
| Q-T1 | P-F | slope of zeta_K(a) on a over [0, 1] | 0 | -0.0398 [-0.0558, -0.0271] (day_stock_bootstrap) | yes | no clip | -0.0400 [-0.0552, -0.0265] (day_stock_bootstrap) | yes |
| Q-T1 | P-F | slope of zeta_K(a) on a over [0, 1] | 0 | -0.0398 [-0.0558, -0.0271] (day_stock_bootstrap) | yes | centred moments | -0.0411 [-0.0570, -0.0294] (day_stock_bootstrap) | yes |
| Q-T3 | P-V | slope in a of zeta_V real minus S1 | 0 | 0.0040 [-0.0258, 0.0337] (day_jackknife_t) | no | clip R1 | 1.55e-04 [-0.0360, 0.0210] (day_stock_bootstrap) | no |
| Q-T3 | P-V | slope in a of zeta_V real minus S1 | 0 | 0.0040 [-0.0258, 0.0337] (day_jackknife_t) | no | clip R3 | -0.0108 [-0.0410, 0.0094] (day_stock_bootstrap) | no |
| Q-T3 | P-V | slope in a of zeta_V real minus S1 | 0 | 0.0040 [-0.0258, 0.0337] (day_jackknife_t) | no | no clip | 0.0061 [-0.0279, 0.0287] (day_stock_bootstrap) | no |
| Q-T3 | P-V | slope in a of zeta_V real minus S1 | 0 | 0.0040 [-0.0258, 0.0337] (day_jackknife_t) | no | centred moments | -0.0891 [-0.1245, -0.0595] (day_stock_bootstrap) | yes |
| Q-T3 | P-V | slope in a of zeta_V real minus S2 | 0 | 0.0187 [-0.0094, 0.0467] (day_jackknife_t) | no | clip R1 | 0.0120 [-0.0187, 0.0316] (day_stock_bootstrap) | no |
| Q-T3 | P-V | slope in a of zeta_V real minus S2 | 0 | 0.0187 [-0.0094, 0.0467] (day_jackknife_t) | no | clip R3 | 0.0017 [-0.0297, 0.0222] (day_stock_bootstrap) | no |
| Q-T3 | P-V | slope in a of zeta_V real minus S2 | 0 | 0.0187 [-0.0094, 0.0467] (day_jackknife_t) | no | no clip | 0.0176 [-0.0131, 0.0419] (day_stock_bootstrap) | no |
| Q-T3 | P-V | slope in a of zeta_V real minus S2 | 0 | 0.0187 [-0.0094, 0.0467] (day_jackknife_t) | no | centred moments | -0.0688 [-0.0962, -0.0413] (day_jackknife_t) | yes |
| Q-T3 | P-R | a*(T1) | 0.35 and 0.65 | 0.4087 [0.3752, 0.4422] (day_jackknife_t) | 0.35: yes; 0.65: yes | clip R1 | 0.4223 [0.3827, 0.4618] (day_jackknife_t) | 0.35: yes; 0.65: yes |
| Q-T3 | P-R | a*(T1) | 0.35 and 0.65 | 0.4087 [0.3752, 0.4422] (day_jackknife_t) | 0.35: yes; 0.65: yes | clip R3 | 0.4332 [0.3964, 0.4701] (day_jackknife_t) | 0.35: yes; 0.65: yes |
| Q-T3 | P-R | a*(T1) | 0.35 and 0.65 | 0.4087 [0.3752, 0.4422] (day_jackknife_t) | 0.35: yes; 0.65: yes | no clip | 0.3982 [0.3672, 0.4291] (day_jackknife_t) | 0.35: yes; 0.65: yes |
| Q-T3 | P-R | a*(T1) | 0.35 and 0.65 | 0.4087 [0.3752, 0.4422] (day_jackknife_t) | 0.35: yes; 0.65: yes | centred moments | 0.4065 [0.3725, 0.4405] (day_jackknife_t) | 0.35: yes; 0.65: yes |
| Q-T3 | P-R | a*(T2) | 0.35 and 0.65 | 0.3577 [0.2994, 0.4096] (day_stock_bootstrap) | 0.35: no; 0.65: yes | clip R1 | 0.3552 [0.3022, 0.4008] (day_stock_bootstrap) | 0.35: no; 0.65: yes |
| Q-T3 | P-R | a*(T2) | 0.35 and 0.65 | 0.3577 [0.2994, 0.4096] (day_stock_bootstrap) | 0.35: no; 0.65: yes | clip R3 | 0.3594 [0.3051, 0.4171] (day_stock_bootstrap) | 0.35: no; 0.65: yes |
| Q-T3 | P-R | a*(T2) | 0.35 and 0.65 | 0.3577 [0.2994, 0.4096] (day_stock_bootstrap) | 0.35: no; 0.65: yes | no clip | 0.3573 [0.3094, 0.4054] (day_stock_bootstrap) | 0.35: no; 0.65: yes |
| Q-T3 | P-R | a*(T2) | 0.35 and 0.65 | 0.3577 [0.2994, 0.4096] (day_stock_bootstrap) | 0.35: no; 0.65: yes | centred moments | 0.3895 [0.3359, 0.4356] (day_stock_bootstrap) | 0.35: no; 0.65: yes |
| Q-T3 | P-M | slope of beta-hat on log q | 0 | 0.0484 [-0.1478, 0.2446] (day_jackknife_t) | no | clip R1 | 0.0484 [-0.1478, 0.2446] (day_jackknife_t) | no |
| Q-T3 | P-M | slope of beta-hat on log q | 0 | 0.0484 [-0.1478, 0.2446] (day_jackknife_t) | no | clip R3 | 0.0494 [-0.1473, 0.2461] (day_jackknife_t) | no |
| Q-T3 | P-M | slope of beta-hat on log q | 0 | 0.0484 [-0.1478, 0.2446] (day_jackknife_t) | no | no clip | 0.0482 [-0.1479, 0.2444] (day_jackknife_t) | no |
| Q-T3 | P-M | slope of beta-hat on log q | 0 | 0.0484 [-0.1478, 0.2446] (day_jackknife_t) | no | centred moments | 0.0484 [-0.1478, 0.2446] (day_jackknife_t) | no |
| Q-T3 | P-M | top size bin beta-hat | 0 | -0.0100 [-0.6698, 0.6498] (day_jackknife_t) | no | clip R1 | -0.0100 [-0.6698, 0.6498] (day_jackknife_t) | no |
| Q-T3 | P-M | top size bin beta-hat | 0 | -0.0100 [-0.6698, 0.6498] (day_jackknife_t) | no | clip R3 | -0.0100 [-0.6698, 0.6498] (day_jackknife_t) | no |
| Q-T3 | P-M | top size bin beta-hat | 0 | -0.0100 [-0.6698, 0.6498] (day_jackknife_t) | no | no clip | -0.0100 [-0.6698, 0.6498] (day_jackknife_t) | no |
| Q-T3 | P-M | top size bin beta-hat | 0 | -0.0100 [-0.6698, 0.6498] (day_jackknife_t) | no | centred moments | -0.0100 [-0.6698, 0.6498] (day_jackknife_t) | no |
| Q-T3 | P-F | slope of zeta_K(a) on a over [0, 1] | 0 | -0.0616 [-0.0988, -0.0381] (day_stock_bootstrap) | yes | clip R1 | -0.0642 [-0.1004, -0.0420] (day_stock_bootstrap) | yes |
| Q-T3 | P-F | slope of zeta_K(a) on a over [0, 1] | 0 | -0.0616 [-0.0988, -0.0381] (day_stock_bootstrap) | yes | clip R3 | -0.0733 [-0.1033, -0.0440] (day_stock_bootstrap) | yes |
| Q-T3 | P-F | slope of zeta_K(a) on a over [0, 1] | 0 | -0.0616 [-0.0988, -0.0381] (day_stock_bootstrap) | yes | no clip | -0.0534 [-0.0862, -0.0277] (day_stock_bootstrap) | yes |
| Q-T3 | P-F | slope of zeta_K(a) on a over [0, 1] | 0 | -0.0616 [-0.0988, -0.0381] (day_stock_bootstrap) | yes | centred moments | -0.0600 [-0.0952, -0.0373] (day_stock_bootstrap) | yes |

### C.1 P-V (i): leave-one-day-out weighted error of log v_a(T), by arm

**Measurement** (computed here). The P-V statistics above are its part (ii); (i) is a gate without an interval.

| tick bin | arm | M1 | M2 | M3 | days used |
|---|---|---|---|---|---|
| Q-T1 | primary (R2, uncentred) | 0.0171 | 0.0159 | 0.0013 | 8 |
| Q-T1 | clip R1 | 0.0174 | 0.0164 | 0.0014 | 8 |
| Q-T1 | clip R3 | 0.0171 | 0.0162 | 0.0012 | 8 |
| Q-T1 | no clip | 0.0167 | 0.0150 | 0.0014 | 8 |
| Q-T1 | centred moments | 0.0167 | 0.0151 | 0.1014 | 8 |
| Q-T3 | primary (R2, uncentred) | 0.0267 | 0.0138 | 0.0213 | 8 |
| Q-T3 | clip R1 | 0.0253 | 0.0134 | 0.0208 | 8 |
| Q-T3 | clip R3 | 0.0258 | 0.0136 | 0.0230 | 8 |
| Q-T3 | no clip | 0.0272 | 0.0129 | 0.0207 | 8 |
| Q-T3 | centred moments | 0.0250 | 0.0084 | 0.1478 | 8 |

### C.2 Settings applied and draws

**Measurement** (computed here). The clip rule and centring flag the wrappers actually applied in each arm, and whether each unit's bootstrap draws (SHA-256 of the weights) equal the primary's.

| arm | clip rule applied | centred flag applied | seed | draws equal to the primary's: Q-T1 | Q-T2 | Q-T3 |
|---|---|---|---|---|---|---|
| primary (R2, uncentred) | R2 | False | 20261102 | yes | yes | yes |
| clip R1 | R1 | False | 20261102 | yes | no | no |
| clip R3 | R3 | False | 20261102 | yes | no | no |
| no clip | none | False | 20261102 | yes | no | no |
| centred moments | R2 | True | 20261102 | yes | yes | yes |

## Planned, not run in Phase 5

Every other secondary, descriptive or robustness item of the registered plan (pages 7 and 8):

- Planned, not run in Phase 5: P-K: the zeta_K(a) curve, unclipped and under each clip rule, with real minus S1 and S2 (plan page 7). Section A gives only the per-cell vk_diff summary.
- Planned, not run in Phase 5: P-C: the collapse exponent chi-hat beside 1/(1 + gamma-hat) and 1 - gamma-hat/2, gamma-hat by local Whittle (bandwidth n^0.65), a power-law fit to the sign autocorrelation and the slope of v_0(T), each with its measured bias (plan page 7).
- Planned, not run in Phase 5: Model comparison M0 (linear in T), M1, M2 and M3 for K_a(T) (plan page 7).
- Planned, not run in Phase 5: The P-R secondaries (R_{1/2}/R_0 against e^{sigma_l^2/8}, the impact slope beside a*(T), da*/dlog T) per cell; section B gives them per tick bin, from Phase 4.
- Planned, not run in Phase 5: The clip and centring arms for Q-T2, for the per-cell units and for P-K and P-C; section C covers the four primary statistics in Q-T1 and Q-T3.
- Planned, not run in Phase 5: Hidden-order prints signed by the recorded rule (side of the mid; same-timestamp prints at the mid take the run's sign; others unsigned and dropped, their share reported) (plan page 7).
- Planned, not run in Phase 5: Fill-level analysis, with sigma_l re-measured on fills (plan page 7).
- Planned, not run in Phase 5: Hierarchical pooling (DerSimonian-Laird with Hartung-Knapp-Sidik-Jonkman variance) (plan page 7).
- Planned, not run in Phase 5: The rescaled axis u = a sigma-hat_l^2 for P-V and P-K (plan page 7).
- Planned, not run in Phase 5: Deseasonalised mid returns (plan page 7).
- Planned, not run in Phase 5: Sub-windows 10:00 to 12:00, 12:00 to 14:00 and 14:00 to 15:30 (plan page 7).
- Planned, not run in Phase 5: The NYSE-listed secondary universe, with the same frozen code (plan pages 3 and 7).
- Planned, not run in Phase 5: The venue-thinning experiment: observed shares {1.0, 0.46, 0.26, 0.17}, hidden share 0.17, random and router-like thinning, kappa in {1, 4, 16} (plan page 8).
- Planned, not run in Phase 5: The NASDAQ-share dose-response check, after the consolidated-volume vendor is reconciled on about ten names (plan page 8).
- Planned, not run in Phase 5: The nine replication days, with the same frozen code (plan pages 2 and 8).
