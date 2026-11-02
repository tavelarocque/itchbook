# Phase 3d-1: the four primary tests as fixed functions

Nothing is frozen yet: the freeze comes after Phase 3. Every change recorded in section 7 is an ordinary
pre-freeze decision.

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken here. No real data were read:
no test ran on exploration or confirmation data, and there is no simulation study in this step.

Code: `tools/est/prereg3d.py` (every pre-registered constant), `tools/est/primary3d.py` (the tests, Holm, the
pooling unit, the range rule), `tools/est/test_3d1.py` (known answers and planted bugs). One addition to
`pooled.py`: the mask spec `("cellceil", rho_per_stock)`.

## 1. P-R: the decision of 31 Oct 2026 and why

### 1.1 The decision

The pre-registered P-R was a one-sided test that a*(T) drifts with T in the direction implied by Maitrier and
Bouchaud Eq. 74. The derivation in section 1.3 shows Eq. 74 does not fix that direction, so the step stopped. On
31 Oct 2026, after the derivation and before any data, P-R was replaced:

- **P-R primary (equivalence).** a*(T) lies within 1/2 +- 0.15 at BOTH T1 and T2 of the fitting range used. TOST
  at each end, intersection-union, so the p-value is the larger of the two ends' TOST p-values.
- **Secondary, not in Holm:** R_{1/2}/R_0 against e^{sigma_l^2/8}; psi-hat reported beside a*(T) (Proposition
  1: the diagonal peak sits at psi); a two-sided estimate of da*/dlog T with its interval, reported only.

This is the paper's stock regime, in which the diagonal term dominates. Source: Maitrier and Bouchaud, arXiv:2506.07711v6
[q-fin.TR], 4 Mar 2026 (PDF SHA-256 ad046872a8a20cf212487ca2b4eba97d24a9c7d80aa05d5e5a5518f421256831),
section 7.3, page 24, quoted exactly:

> "Finally, for point 3., Fig. 9 shows that the correlation between price returns and generalized order
> imbalance is non-monotonic as a function of a, reaching a maximum for a⋆ ≈ 0.5 for LLOY, TSCO and SPMINI, and
> a⋆ ≈ 1 for EUROSTOXX."

> "Such a non-monotonic behaviour is predicted by our theoretical analysis. Interestingly, when the diagonal
> contribution to Ra dominates, we expect that the maximum correlation is reached precisely for a⋆ = 1/2, with a
> peak amplitude that decreases with T, see Eq. (71). The data for the two stocks is therefore compatible with
> the fact that in the small a regime, Ra (T)d. is a decreasing function of T."

(Page 24 of the PDF. The quote was checked against the rendered page. SPMINI is a future; the paper groups it
with the two stocks at a⋆ about 0.5.)

### 1.2 Eq. 74

Section 6.3, page 21, for a < a_c:

    R_a(T) = exp(-s a^2 / 2) * ( A(T) exp(s a / 2) + B(T) exp(lambda s a tau) ),

with s = sigma_l^2, tau = log T, and A and B functions of T that the paper does not specify. The equation was
checked against the rendered page.

### 1.3 Why the one-sided drift test was dropped (the derivation, line by line)

**1. Completing the square.**

    R_a(T) = A(T) exp(s/8) exp(-s (a - 1/2)^2 / 2) + B(T) exp(s (lambda tau)^2 / 2) exp(-s (a - lambda tau)^2 / 2).

Each term is a Gaussian bump in a with width 1/sqrt(s). The A bump is centred at 1/2 for every T. The B bump is
centred at lambda log T.

**2. The log derivative.**

    d log R_a / da = s ( -a + (1 - w)/2 + w lambda tau ),
    w = B exp(lambda s a tau) / ( A exp(s a / 2) + B exp(lambda s a tau) ).

The log derivative of a sum is the share-weighted average of each term's log derivative. w is the B term's
share.

**3. The peak.**

    a*(T) = (1 - w*) / 2 + w* lambda tau.

The peak is a weighted average of the two centres.

**4. The pure regimes.**
- w* = 0, the diagonal term dominant: a* = 1/2 at every T, so there is no drift. This is the regime of section
  7.3, and the regime the new P-R tests.
- w* = 1: a* = lambda log T, which increases with T.

**5. The mixed regime.** With G = -a + 1/2 + w (lambda tau - 1/2), w = logistic(u),
u = log(B/A) + s a (lambda tau - 1/2), and g = d log(B/A) / d tau:

    da*/dtau = [ w lambda + (lambda tau - 1/2) w (1 - w) (g + s lambda a*) ] / [ 1 - s w (1 - w) (lambda tau - 1/2)^2 ].

This is implicit differentiation of G = 0. The denominator is positive at a maximum, because
d^2 log R / da^2 = s dG/da < 0 there.

**6. The sign is open.** The first term is positive. The second has the sign of (lambda tau - 1/2)(g + s lambda a*):
- g is not given by Eq. 74;
- with the paper's Fig. 9 lambda = 0.1, lambda log T is 0.14 to 0.42 on T = 4 to 64, below 1/2.

So a B/A ratio that grows with T pulls a* down while lambda log T pushes it up.

**7. Sufficient conditions for a sign, none implied by Eq. 74.** Any one of these gives a positive drift:
- g = 0 and s < 4;
- lambda log T >= 1/2 and g >= 0;
- w* = 1.

With w* = 0 the drift is zero.

**8. Numerical check (Measurement).** s = 0.81, lambda = 0.1, A = 1. Each cell is the exact argmax of Eq. 74 /
the quadratic vertex on a in {0, 0.25, 0.5, 0.75, 1}.

| B/A | T = 4 | 8 | 16 | 32 | 64 | direction |
|---|---|---|---|---|---|---|
| 0 (A only) | 0.500 / 0.500 | 0.500 / 0.500 | 0.500 / 0.500 | 0.500 / 0.500 | 0.500 / 0.500 | none |
| infinite (B only) | 0.139 / 0.109 | 0.208 / 0.196 | 0.277 / 0.275 | 0.347 / 0.348 | 0.416 / 0.418 | up |
| 1 | 0.328 / 0.329 | 0.360 / 0.362 | 0.393 / 0.395 | 0.425 / 0.427 | 0.459 / 0.460 | up |
| 0.3 T^0.26 | 0.400 / 0.402 | 0.407 / 0.409 | 0.419 / 0.421 | 0.437 / 0.438 | 0.461 / 0.462 | up |
| 0.05 T | 0.446 / 0.447 | 0.422 / 0.424 | 0.405 / 0.407 | 0.407 / 0.410 | 0.436 / 0.438 | down, then up |
| 20 / T | 0.202 / 0.189 | 0.296 / 0.295 | 0.380 / 0.382 | 0.443 / 0.445 | 0.481 / 0.481 | up |

The exponent 0.26 comes from the paper's defaults (mu_m = 3/2, lambda s = 1/8, gamma_x + beta_hat(0) = 0.8),
through the T powers of Eqs. 71 and 72.

## 2. API (`primary3d`)

**Input.**
- `Flow(eps, q, lm, ok, ts, starts, lengths, n_days)`: stocks' stock-days, stock-major, with RAW sizes q. The
  R2 clip is applied when a unit is built.
- `ts` is in ns since midnight ET, for S2's 30-minute bins.

**Units.**
- `tick_bin_unit(flows, "Q-T1")` (primary): pools the tick bin's confirmatory cells. Descriptive cells in
  `flows` are dropped and logged in `Unit.excluded_cells`.
- `cell_unit(flows, cell)` (secondary).
- `make_unit` raises for any non-confirmatory cell.
- Each stock carries its cell's ceiling from `prereg3d.CONFIRMATORY` (mask `("cellceil", rho)`).

**Running.**
- `run_primary(flows, pm, rng, n_boot=200, n_sur=4)` is the primary analysis. For each tick bin:
  - build the unit;
  - decide `evaluability(Ns, rho, pm)` by the reach rule alone;
  - run `run_unit` on the evaluable tests only.

  Then `combine` runs one Holm family over the evaluable (test, tick bin) pairs, and `set_decisions` sets
  `rejected` and `passed`. It returns `{"by_tick_bin", "evaluability", "family", "left", "holm"}`, where `left`
  logs every pair that left the family, with its reason.
- `run_tick_bin(flows, tb, ...)` gives one tick bin's results before the family's Holm. `run_cell(...)` is
  secondary and never in the family. Both return `{"tests": {name: TestResult}, "descriptive": {"P-K": ...},
  "unit": {...}}`.

**`TestResult`.**
- `p`, `alternative`, and `conditions` (every sub-condition as a bool).
- `statistic`: the range used, with estimates and widest intervals.
- `range_used`, `range_rule` (primary, fallback or not evaluable), `reach` (stocks at T2 of each range) and
  `stocks_used` (indices at T2 of the range used).
- `by_range`: both ranges, reported. `secondary`: never in Holm.
- `rejected` (the family's Holm) and `passed` = rejected AND every condition.

**Pure decision layers.** These take replicates (`Reps(est, ds, jk, ws)`): `decide_pv(errors, s1, s2)`,
`decide_pr(ends)`, `decide_pm(slope, top)`, `decide_pf(pf)`.

**Other functions.**
- `holm(pvals)`: keys must be (primary test, tick bin) pairs, so at most 12. P-K, P-C, a bare test name or an
  unknown tick bin raise.
- `combine(ev, pvals)`: refuses a p-value for a pair that left the family, and a family member without a
  p-value.
- `select_range(Ns, rho)`: the reach rule. Its only inputs are day lengths and ceilings.
- `PMParams(n_bins, l_max, L, min_runs_per_bin)`: validated against prereg3d (bins 3, 4 or 5; l_max 16 or 32;
  L >= 8 l_max). There is no default; Phase 3d-2 chooses.
- Building blocks: `m3_V`, `logv_rel`, `fit_m1`, `fit_m2`, `model_errors`, `pm_moments`, `pm_solve`,
  `fit_beta`, `psi_moments`, `psi_hat`, `vertex`, `interval`, `p_value`.

## 3. The tests as implemented

All four use the T grid 1..64, the a grid {0, 0.25, 0.5, 0.75, 1}, R2-clipped sizes, the tilde path of
`pooled.core2` with the reach rule, and the per-stock ceilings.

**Inference.** One set of draws per unit, shared by every test, in the order of `pooled.intervals3`:
- the day-and-stock bootstrap;
- the whole-stock bootstrap;
- the delete-one-day jackknife with Student t at D - 1 df.

The interval is the widest of the three at 0.90; `check_interval_matches_intervals3` checks equality with
`pooled.intervals3`. **Choice:** the p-value is the largest of the three methods' p-values, the counterpart of
the widest interval. A percentile bootstrap p-value is the share of replicates on the null side, and twice the
smaller share when two-sided.

**Range.** The range used is chosen by `select_range`: [8, 64] if at least 5 stocks reach 64, else [4, 32] if
at least 5 reach 32, else not evaluable (p = 1, every condition false).
- **Choice:** a statistic undefined on the chosen range (for example zeta_K dropping a horizon) gives p = 1. It
  never triggers a switch to the other range.

**P-V.**
- **(i)** Leave-one-day-out over the D days. For held-out day d:
  - M1 and M2 are fitted to the pooled curve of the other days;
  - M3 is built from day d's own stock-days;
  - each model's prediction of day d's pooled curve is scored.
  - **Choice:** the curve is log v_a(T) - log v_a(T1) over the range's horizons, so there are no free
    intercepts. The models predict (zeta_V(a) - 1)(log T - log T1).
  - **Choice:** the weights are w_T = 1/T, about the inverse variance of log v.
  - M1: zeta_V constant in a.
  - M2: zeta_V(a) = 3 - mu_m - 2 a lambda sigma_l^2 below a_c, 1 above. **Choice:** lambda sigma_l^2 >= 0; a_c
    leaves 2 to 5 grid points below it; ties go to the smaller a_c.
  - M3: V_a(T) = T E[q^2a] + 2 E[q^a]^2 sum_{l<T} (T - l) C(l), from the stock-day's uncentred sign ACF and
    size moments, pooled exactly like the real V. `check_m3_is_shuffle_expectation` checks it against the mean
    over 300 shuffles: within 3 %.
  - **Choice:** a day with any undefined value is skipped for every model.
  - (i) holds when M2's total error is below both M1's and M3's.
- **(ii)** The slope in a, over the grid points below M2's a_c, of zeta_V(a) real minus zeta_V(a) surrogate.
  - **Choice:** a_c is from M2 fitted on all days and held fixed across replicates.
  - **Choice:** the surrogate zeta_V is the mean over N_SUR = 4 draws. S2's block length comes from
    `surrogates.choose_L` per stock-day.
  - The test is one-sided (less), for S1 and S2 separately; each needs its upper bound below 0.
- **p** = max(p_S1, p_S2) when (i) holds, else 1. **Choice:** (i) is a gate without a p-value.

**P-R.**
- a*(T) (decided 1 Nov 2026, section 7.3): the vertex of the parabola through the grid maximum of the
  pooled Fisher-z R_a(T) (`pooled.core2`'s R) and its two neighbours (`primary3d.vertex`). A maximum at a grid
  end is censored there, and P-R fails at that T. Until 1 Nov 2026 a* was the vertex of the least-squares
  quadratic over the whole grid, NaN unless concave (`pooled.core2`'s a_star, now `primary3d.vertex_ls`).
  **Choice:** it is NaN below the range's reach.
- TOST at T1 and at T2: H0 a* <= 0.35 and H0 a* >= 0.65. p_end is the larger of the two, and
  p = max(p_end(T1), p_end(T2)).
- Conditions: the 0.90 interval (the TOST interval) lies inside (0.35, 0.65) at T1 and at T2.
- Secondary:
  - R_{1/2}/R_0 at the range's middle horizon, against e^{sigma_l^2/8}, as a difference with interval and
    two-sided p. **Choice:** sigma_l^2 is each stock-day's variance of log clipped size, equal weight per stock.
  - psi-hat, the pre-registered definition (decision of 31 Oct 2026; section 3.1).
  - da*/dlog T: the OLS slope over the range's horizons, two-sided p.

**P-M.**
- **Choice:** size bins are each stock-day's own tie-aware bins (`tie_bins`, section 7.2; until Phase 3d-2b,
  quantiles with ties to the lower bin).
- x_b(n) = eps_n 1[bin b]. r_n = lm_{n+1} - lm_n over valid steps, divided by the stock-day's rms (**Choice**).
- Moments: Cp[b, b', k] = mean x_b(n) x_b'(n + k) and Sr[b, l] = mean x_b(n) r_{n+l}, for k, l < L.
- They are pooled equal-weight per stock (days first), then solved for the size-binned propagator
  S_b(l) = sum_b' sum_{j<L} K_b'(j) C_bb'(l - j) (block Toeplitz). G_b(l) = sum_{j<l} K_b(j).
- beta_b comes from `propagator.fit_power`'s model G0 l0^beta / (l0 + l)^beta on lags [1, l_max].
  **Choice:** the beta grid is -1 to 1.5, because the top-bin condition needs beta < 0 to be reachable.
- The slope is OLS of beta_b on x_b, the bin's mean log size relative to the stock-day's.
- **Choice:** a stock-day is used if:
  - its realised bin count equals n_bins, so every bin holds >= min_runs_per_bin runs (section 7.2);
  - N > L + 1;
  - its stock reaches L under its cell's ceiling (L <= rho N_s).

  P-M needs at least 5 stocks.
- p = max(p_slope, p_top), both one-sided (less). Conditions: the slope's upper bound is below 0, and the top
  bin's upper bound is below 0.

**P-F.**
- The OLS slope of zeta_K(a) on a over [0, 1] (`pooled.core2` "pf"). Two-sided against 0; the condition is that
  the interval excludes 0.
- zeta_K(a) levels are reported with intervals, never tested.

**Numerical zero (Choice).** Every slope statistic within 1e-12 of 0 is set to 0. An exactly flat zeta_K gave
+-3e-17 from the OLS weights (Measurement), and a consistent sign across replicates would then reject flatness on
rounding alone.

**Holm (decision of 31 Oct 2026).** ONE family across all tick bins: the 4 tests x 3 tick bins give up to 12
p-values, at family-wise 0.05.
- Members are (test, tick bin) pairs. P-K and P-C cannot enter, and NaN counts as 1.
- Ascending p; the i-th smallest is compared with 0.05/(m - i), with m the family size; stop at the first
  failure. Adjusted p-values are reported.
- A pair that is not evaluable by the reach rule alone (`evaluability`: day lengths and ceilings, never an
  estimate) leaves the family before any p-value is computed. Its test is not run at all, and the pair is logged
  in `left` with its reason:
  - P-V, P-R, P-F need a range from `select_range`;
  - P-M needs at least 5 stocks with L <= rho N_s.
- A pair that is evaluable but whose statistic turns out undefined stays in the family with p = 1. That covers,
  for example, fewer than 5 P-M stocks after the min-runs rule, or zeta_K dropping a horizon.
- This replaces the earlier per-tick-bin Choice. Only the combining step changed; the four tests are as before.

### 3.1 psi-hat (decision of 31 Oct 2026)

psi-hat is exactly the pre-registered definition (`psi_moments`, `psi_hat`):
- per stock, take the quantile bins of run size over all its study days (R2-clipped sizes), made tie-aware in
  Phase 3d-2b (section 7.2);
- in each bin, E[Delta_1 eps | q-bin] is the mean of eps_n (lm_{n+1} - lm_n) over steps with ok at n and
  n + 1 inside a day;
- the stock's psi is the OLS slope of log E[Delta_1 eps | q-bin] on log(median q of the bin), over the bins;
- the tick bin's psi-hat pools the stocks: the equal-weight mean of the stocks' slopes. Replicates carry their
  own day and stock weights.

**How it differs from the earlier version.** The earlier psi-hat (Choice, now removed) had four differences:
- it used deciles within each stock-day, not per stock;
- it regressed the binned response on c0 + c1 mean (q/gm)^psi, which allowed an additive constant;
- it did not take logs;
- it fitted one psi common to all stocks, by minimising the summed 1 - R^2 over a grid.

**Choices left inside the definition:**
- 10 bins (the definition fixes quantile bins but not how many);
- ties go to the lower bin (replaced in Phase 3d-2b by tie-aware bins of at least 100 steps, and a slope needs
  at least 3 realised bins; section 7.2);
- the bin edges and medians are fixed per stock from all its days and held fixed across replicates, so only E
  is re-weighted;
- a stock with a bin whose E is not positive has no slope (log undefined) and is left out. The number of stocks
  with a slope is reported as `psi_n_stocks`.

**Descriptive.** P-K is reported as the slope difference vk_diff with its interval. P-C has no definition in
the repo; see section 6.

## 4. Known answers (Measurement, `test_3d1.py`)

| input | result |
|---|---|
| Eq. 74, B = 0, A(T) = 0.3 T^-0.1, replicate noise 0.002 (jackknife noise 0.002 / sqrt(8), from 1 Nov 2026) | a* = 0.500 at both ends; P-R p = 0.0002 (T1), 0.0004 (T2) with the local vertex; passes |
| R peaked at a* = 0.25 at both ends | P-R p = 1.0; both conditions fail |
| a* = 0.5 at one end, 0.2 at the other (either end) | p > 0.5; one condition fails |
| G0 flow, sizes independent of signs (10 stocks x 6 days x 1500 runs) | (ii) p_S1 = 1.0, p_S2 = 0.98; slopes +0.022, +0.016; M3 (the true model) beats M2; P-V fails |
| big sizes on random-sign runs, small on persistent ones | (ii) slopes -0.65 (S1), -0.64 (S2), p 4e-6; M2 beats M1 and M3 |
| k-tilde with slope in a exactly 0 | P-F p = 1 on both ranges, interval contains 0 |
| zeta_K rising or falling in a by 0.05 | P-F p < 0.01 both ways |
| size-binned kernels beta = (0.6, 0.2, -0.3) | beta-hat (0.62, 0.20, -0.30); p_slope 9e-5, p_top 2e-5; both conditions hold |
| iid signs, permanent impact eps q^psi (E[Delta_1 eps \| q] = q^psi), 8 stocks x 4 days x 4000 runs | psi-hat 0.314 for psi = 0.3, 0.628 for psi = 0.6 (within 0.03, the check's tolerance); the small upward bias comes from evaluating log E[q^psi] at the bin median |
| G0 flow (psi = 0.5, sign memory, transient kernel), 10 x 6 x 1500 | psi-hat 0.537 (reported, not checked: the past runs' decaying impact adds a q-independent term) |
| one Holm family, one p = 0.01 per tick bin, rest 0.5 | nothing rejected (0.01 > 0.05/12); per-bin Holm would reject three |
| two pairs below the reach, p = 0.0048 elsewhere | family of 10, 0.0048 <= 0.005 rejected; kept with p = 1 (m = 12) it would not be |
| an evaluable pair with p = 1 | stays: family of 12, so 0.0045 > 0.05/12 is not rejected |
| run_primary on Q-T3 (3 cells) + Q-T2S3 (200-run days), no Q-T1 cell | family of 7: Q-T3 x 4, Q-T2 P-V/P-R/P-F on [4,32]; P-M of Q-T2 and all of Q-T1 left, with reasons; Q-T2's P-M never computed |

## 5. Tests and planted bugs

`tools/est/test_3d1.py` has 31 checks and 23 planted bugs. Each mutant was also run outside pytest:
- 20 fail their check on an assertion;
- three are caught one layer earlier, by a guard raised inside the check:
  - "descriptive cells pooled into the tick bin": `make_unit`'s guard (ValueError: "Q-T1S3 is not a
    confirmatory cell");
  - the two "evaluability after the p-value" mutants: `combine`'s refusal of a p-value for a pair outside the
    family (ValueError: "p-values given for pairs outside the family"). That refusal is how "leaves before any
    p-value is computed" is enforced.

| planted bug | mutated function | caught by |
|---|---|---|
| Holm sorted descending | `holm` | check_holm_order |
| Holm as plain Bonferroni | `holm` | check_holm_order (P-F at 0.015: Holm 0.0167 rejects, Bonferroni 0.0125 does not) |
| Holm not stopping at the first failure | `holm` | check_holm_order (P-V at 0.04 must stay unrejected) |
| Holm accepting a member that is not a (primary test, tick bin) pair | `holm` | check_holm_family_keys (P-K, P-C, bare name, unknown tick bin) |
| Holm per tick bin instead of one family | `combine` | check_one_family_across_bins |
| non-evaluable pair kept in the family with p = 1 | `combine` | check_nonevaluable_leaves_family |
| evaluability decided from the computed p-value | `combine` | check_evaluability_reach_only (combine's guard) |
| every test computed before evaluability is applied | `run_primary` | check_run_primary (combine's guard) |
| fallback range picked by the smaller p | `_assemble` | check_fallback_never_by_result |
| [8, 64] used whatever its reach | `select_range` | check_fallback_by_reach |
| descriptive cells pooled into the tick bin | `tick_bin_unit` | check_descriptive_excluded (make_unit's guard) |
| one ceiling for every stock of the tick bin | `make_unit` | check_cell_ceilings |
| R2 clip skipped | `make_unit` | check_r2_clip |
| p-value from the least conservative method | `p_value` | check_p_value_conservative |
| P-F tested one-sided | `decide_pf` | check_pf_two_sided |
| P-F slope not snapped at numerical zero | `snap` | check_pf_flat_not_rejected |
| P-R one-sided on drift instead of the equivalence | `decide_pr` | check_pr_eq74_known_answer |
| P-R at T2 only | `decide_pr` | check_pr_uses_both_ends |
| M3 missing the factor 2 of the sign ACF term | `m3_V` | check_m3_is_shuffle_expectation |
| P-V (ii) on S1 only | `decide_pv` | check_pv_both_surrogates |
| P-V without the model-comparison gate | `decide_pv` | check_pv_gate |
| P-M on the slope condition only | `decide_pm` | check_pm_both_required |
| P-M top bin taken as the smallest sizes | `pm_stats` | check_pm_recovers_beta |

Other checks: config table, the cell ceilings in a pooled tick bin, the M2 fit recovering its parameters,
leave-one-day-out choosing the true model, vertex_ls = `pooled.core2`'s a*, P-M parameter validation and the
min-runs and ceiling rules, the whole primary run (family, log of pairs that left), and psi-hat's known answer.

Full suite (`uv run --group dev --group analysis pytest`): 398 passed, 28 skipped (was 344 passed before
Phase 3d-1). The new file runs in about 4 s.

## 6. Open (not decided here)

- **P-C** has no definition in the repo. Holm structurally refuses it, but no P-C output is produced here.
- **P-M's top-bin condition, read literally.** It is beta_top's upper bound below 0. With G(l) ~ l^-beta, that
  means impact grows with lag in the top size bin. It is implemented as written; if the intended reference
  value is not 0, it is one constant in `decide_pm`.
- **P-M bins on real sizes.** With a mass of runs at one size (round lots), quantile edges coincide; ties go to
  the lower bin, and min_runs_per_bin then drops the stock-day. Resolved in Phase 3d-2b by tie-aware bins
  (section 7.2).
- **Proposition 1** is not in the repo or in either Maitrier and Bouchaud paper. psi-hat follows the
  pre-registered definition (section 3.1). Under sign memory with a decaying kernel (G0), past runs add a
  q-independent term to E[Delta_1 eps | q], and psi-hat reads 0.537 against a true 0.5.
- **Compute cost** at confirmation tick-bin sizes, with N_SUR = 4 and B = 200, is not measured.

## 7. Pre-freeze decisions

Nothing is frozen until the end of Phase 3. Every change below was made before any confirmation data were read
and before any power result of Phase 3d-2.

### 7.1 `m3_V`: pre-freeze bug fix (Phase 3d-2; approved 1 Nov 2026)

- **The bug.** `m3_V` computed the sign ACF to lag 63 (the largest horizon minus 1) on every stock-day,
  whatever its length. `propagator.sign_acf` uses an FFT of size 2N, so:
  - on a day of fewer than 32 runs it crashed ("ValueError: operands could not be broadcast together with
    shapes (8,) (64,)");
  - on a day of 32 to 63 runs it returned lags that wrapped around the day, so M3 was wrong at every horizon.
- **Why it mattered.** Real Q-T1 data have many short days: 164 of the 552 exploration stock-days of Q-T1S1
  have fewer than 64 runs, the shortest has 1. `run_primary` would have crashed on the real Q-T1 tick bin.
- **The fix.** The ACF is computed only to lag N - 1, and M3 is NaN at a horizon T > N, where the day has no
  window and the real V is also undefined. A day with any undefined value was already skipped for every
  model (section 3, P-V (i)), so nothing else changed.
- **Check.** `test_3d2.check_m3_short_days` (days of 5, 40 and 200 runs: finite where the day has a window, NaN
  beyond), with the old line planted back as a bug. The 3d-1 tests are unchanged and pass.

### 7.2 Tie-aware size bins for P-M and psi-hat (Phase 3d-2b, Choice)

Made from the pilot's bin occupancy (docs/phase3d2_report.md), before any power result.
- **Rule** (`tie_bins(q, n_target, min_runs)`): on the sorted sizes, each of the n_target - 1 quantile cuts
  moves to the nearest boundary between distinct sizes (equal distance: the lower boundary). So every run of
  one size falls in one bin. Cuts that coincide drop out. Then the smallest bin (the first on ties) is merged
  with its smaller neighbour (the lower on ties) until every bin holds at least min_runs runs, or one bin is
  left. The realised bin count is recorded.
- **P-M**: bins per stock-day, target n_bins, min_runs = min_runs_per_bin. A stock-day is used only when its
  realised count equals n_bins, because the pooled moments need the same bins on every stock-day.
- **psi-hat**: bins per stock over all its study days' steps, target 10. **Choice:** at least 100 steps per
  bin (PSI_MIN_RUNS), and a stock's slope needs at least 3 realised bins (PSI_MIN_BINS). Bins beyond a
  stock's realised count are NaN and carry no weight.
- **Why.** Round lots are 28 to 49 % of runs in the confirmatory tick bins. Under the old rule, coinciding
  quantile edges left empty bins, so psi-hat was undefined for every Q-T3 stock and for all but a handful of
  Q-T1 stocks, and P-M lost up to half its stock-days at 5 bins.
- **Checks** (`test_3d2.check_tie_bins`): every size in one bin, no empty bin, every bin at least min_runs, a
  known answer on a tied sample, P-M using only stock-days at the target count, psi's realised count. Planted
  bugs: a size value split across two bins (rank cuts); an empty bin kept (coinciding cuts kept, no merge); a
  stock-day below the target count used by P-M.

### 7.3 a*: the local vertex at the grid maximum (decided 1 Nov 2026)

- **Rule** (`primary3d.vertex`, `primary3d.censored`). On the pooled Fisher-z R_a(T) at each T, take the grid
  maximum a_k (the first on ties). Then a* = a_k + (h / 2) (R_k-1 - R_k+1) / (R_k-1 - 2 R_k + R_k+1), with
  h = 0.25. So a* lies within h / 2 of a_k, and a* = a_k when the three values are equal.
- **Grid ends.** If the maximum is at a = 0 or a = 1, a* is censored at that end: its value is that end's a,
  and it is flagged low or high. A censored point estimate at T1 or T2 fails P-R at that end (p = 1 there,
  its condition False). Censored bootstrap and jackknife replicates keep the end value, which lies outside
  the margin, so they count on the null side. da*/dlog T uses the same values.
- **Why.** In 3d-2b, on real within-day size distributions without sign memory, the expected R(a) peaks at
  psi = 0.5 but falls more slowly to the right than to the left. The least-squares quadratic over the whole
  grid then puts its vertex to the right of the peak.
  - On the sizes-only curve, the old a* on {0, ..., 1} is 0.556 to 0.579, against an argmax of 0.51.
  - The panels in W0-shuf gave a*(1) of 0.57 to 0.58.
  - The local vertex on the same curves is 0.520 to 0.529.
  
  The local vertex is the estimator closest to the argmax of R(a), which is what a* stands for.
- **Check (Measurement, `diag3d2b.py abias` with `PHASE3D2B_OUT=.../diag3d2b_astar`).** These are the same 30
  panels per tick bin as 3d-2b (same seeds); the old a* reproduces 3d-2b exactly. Censored panels: none,
  except 2 high at T = 64 in W0-shuf Q-T2. a*(1), local (old):

  | world | Q-T1 | Q-T2 | Q-T3 |
  |---|---|---|---|
  | W0-shuf (no sign memory) | 0.540 (0.583) | 0.537 (0.581) | 0.534 (0.570) |
  | W0 (real signs) | 0.597 (0.631) | 0.577 (0.613) | 0.582 (0.614) |

  Without sign memory, a*(1) moves to within 0.04 of 0.5. Of the remaining +0.03 to +0.04, the sizes-only
  curve accounts for +0.02 to +0.03: the three-point stencil on a skewed curve at h = 0.25.
- **Checks** (`test_3d2.check_astar_local`, `check_astar_censored`, `check_pr_uses_local_vertex`). A known
  answer on 3d-2b's skewed curve (0.531; the global quadratic gives 0.58), a maximum at a = 0.25, within h / 2
  of the grid maximum, symmetric parabolas exact, NaN kept. Grid ends censored, and P-R failing there with the
  replicates inside the margin. run_unit's a* equal to vertex of core2's R and different from core2's a_star.
  Planted bugs:
  - the global quadratic used (in vertex; and P-R reading core2's a_star);
  - neighbours taken from one side of the maximum;
  - neighbours swapped;
  - a censored end treated as interior (in vertex; and in the P-R decision).
- **Unchanged.** `pooled.core2` still returns the old a_star, so Phases 3c to 3d-2b reproduce. Phase 3d reads
  a* only through `primary3d.vertex`.

### 7.4 psi-hat: no correction (decided 1 Nov 2026)

psi-hat keeps its definition. In 3d-2b's W0-shuf (psi = 0.5, no sign memory) it read 0.553, 0.532 and 0.521
in Q-T1, Q-T2 and Q-T3, so +0.02 to +0.05 high. Every report of psi-hat carries this (`primary3d.PSI_NOTE`:
P-R's secondary, the analyse3d2 tables, the diag3d2b tables).

### 7.5 WM and the scope of the 3d-2 study (decided 1 Nov 2026)

- WM stays the published PAPER2 preset, with beta_q clipped at 0. P-M's top-bin condition (beta_top below 0)
  cannot hold there. WM's P-M is reported as "fails by construction" (`worlds3d2.BY_CONSTRUCTION`), and it is
  not scored. WB remains the P-M power world.
- `phase3d2.py all` runs W-null (the clean coin-flip version), W0, W0-psi3 and WB, with the tie-aware bins:
  - design: W-null, W0, WB at every P-M design (18 when decided; 6 since the trim of 7.6);
  - power: W0-psi3 at the chosen design;
  - analyse.
- WM (`phase3d2.py wm`) and thinning (`phase3d2.py thin`) are separate commands, kept for later.
- With WM out of `all`, P-V has no power world in the study. WM was P-V's power world
  (`worlds3d2.PV_POWER_WORLD`).

### 7.6 The P-M design grid trimmed to 6 designs (decided 1 Nov 2026), and the double run

- **Decision.** The P-M design grid goes from 18 designs to 6: bins {3, 4, 5} x l_max {16, 32}, with
  L = 8 l_max and min runs per bin fixed at 50 (`worlds3d2.PM_DESIGNS`; the old grid is kept as
  `PM_DESIGNS_18`). It was taken before any study result was read.
- **Reason.** Compute cost only. No rejection rate from the finished jobs was read or summarised, before or
  after the decision.
- **Unchanged.** Worlds, panel counts, seeds, the four tests, the one Holm family, and the design rule
  (`worlds3d2.choose_design`: admissible in W-null and W0 in every tick bin, highest mean WB rate, ties to
  fewer bins then smaller l_max), which now runs over the 6. The min-runs tie-break no longer has anything to
  break.
- **The finished jobs are reused, not recomputed.** P-M at one design does not depend on which other
  designs share the job: every design uses the same draws, and a result is reused only between designs with
  the same moments, which is exact. So the 6 retained designs' records in an 18-design chunk are what a
  6-design job gives. `phase3d2.py trim` kept only the 6 in the 130 finished chunks. Each original moved
  unchanged to `~/CODE/itch-data/phase3d2/study_18designs/`; every other field and every retained record is
  identical. W0 chunk 058, recomputed alone at the 6 designs, matched the saved chunk in every field but the
  timings.
- **Incident: two instances of `phase3d2.py all` ran at once.** They ran from the start of the run
  (both logged "750 jobs, 750 to run") against the same chunk folder and the same `run.log`, until the run
  was stopped at 130 of 750 jobs (all Q-T1: W-null chunks 0 to 198, W0 chunks 0 to 58).
  - Both instances computed the same jobs from the same seeds and wrote them under the same names. So no
    job has two files, but a job's file may have been written twice, and both writers used the same
    temporary name.
  - **Integrity check (`phase3d2.py check`, Measurement).** All 130 chunks loaded and were complete: the
    right panels in order, every field, a p in [0, 1] for every evaluable test, all 18 designs. No panel was
    held twice and no temporary file was left. All 735 JSON files under `~/CODE/itch-data/phase3d2/` load.
    Quarantined: 0 (the folder `quarantine/<time>/` is created only when something moves).
  - Effect: the two instances only shared the machine, so the jobs ran slower. The results themselves are
    deterministic in the seeds, and the recomputation below confirms one chunk.
  - **Guard.** Every `phase3d2.py` step now takes an exclusive, non-blocking flock on
    `<output>/phase3d2.lock` (`phase3d2.acquire_lock`), held until the process ends; the OS drops it however
    the process ends, so a crash leaves no stale lock. A second instance on the same output directory exits
    with the holder's pid and start time.
- **Timing after the trim (Measurement, 1 Nov 2026).** Two jobs (2 panels each), one after the other, with
  no other study process running:
  - W0 Q-T1 chunk 058: 118 s (per panel: generate 3.7 / 2.7 s, range tests 21.6 / 17.3 s, P-M 36.8 / 32.8 s);
  - W0 Q-T3 chunk 000: 66 s (range tests 6.9 / 6.1 s, P-M 26.1 / 23.5 s).

  Remaining: stage A has 620 jobs (W0 Q-T1 70, WB Q-T1 50, every Q-T2 and Q-T3 job), stage B 150. WB, Q-T2
  and stage B are scaled from these two jobs with the pilot's and smoke run's ratios; they were not
  measured here. Total: about 11 worker-hours. That is 3.6 h wall on 3 workers if they scale as in the
  pilot, and about 7 h if they slow down as much as in the smoke run. Expect 4 to 7 h.
- **Checks** (`test_3d2`): `check_design_grid` (the 6 designs; the rule and stage A over exactly them),
  `check_lock_refuses_second` (second acquire refused; free again after release and after the holder's
  process dies), `check_main_takes_lock` (main refuses, in-process and from a second process),
  `check_verify_chunk` (truncated chunk, wrong panel, missing design, evaluable test without p, a name that
  is no job, short chunk and leftover .tmp all caught and quarantined, and only they), and
  `check_trim_keeps_retained`. Planted bugs: no flock; a shared lock; main not taking the lock; a chunk's
  panels not checked; a bad chunk left in place; trim keeping the wrong 6; trim not archiving the original.
  `check_design_rule` and `check_design_rule_ties_and_fallback` now use designs of the 6.
