# Phase 3d-2 report: operating characteristics of the four primary tests (design, pilot, driver)

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken here. Real data were read
only through `explore.load_explore_runs`, exploration stocks of the six confirmatory cells, 8 study days, and
only their signs, sizes and timestamps (`worlds3d2.build_flow_cache` keeps exactly those three fields). No
replication file was opened. Every price response is synthetic; nothing combines real signs or sizes with real
price changes.

**Superseded in part by Phase 3d-2b (`docs/phase3d2b_report.md`).** Since this report, W-null has been
redefined: its informed runs now get fair-coin signs. The old W-null is kept as W-null-real. psi-hat and P-M
now use tie-aware size bins (`primary3d.tie_bins`). The pilot's numbers below predate both changes. In
particular, the pilot's W-null P-F rejections (8 of 10) came from the real informed signs and do not occur
under the clean null (0 of 20 per tick bin).

**And by Phase 3d-2c (`docs/phase3d2c_report.md`, 1 Nov 2026).** a* is now the local vertex at the grid
maximum. `all` runs W-null, W0, W0-psi3 and WB only. This pilot's files moved to
`~/CODE/itch-data/phase3d2/report/pilot_3d2/`; `report/pilot.*` is the 3d-2c re-pilot.

## Headlines

1. **The study does not fit this session (Measurement).** The pilot (150 panel-tick-bins, 30.7 min on 2 workers)
   and single-job measurements put the full chain (stage A at 18 P-M designs, stage B, thinning) at about 45
   worker-hours, about 17 h of wall time within the 3-worker, 3 GB limits (section 8). Step 2 (and steps 3 to 5,
   which
   need its panels) was therefore not run here. The driver is resumable; the one-line command is in section 8.
2. **A defect in a frozen function (Measurement), fixed (Choice, needs your sign-off).**
   `primary3d.m3_V` took each stock-day's sign ACF to lag 63 regardless of its length. On a day of fewer than
   32 runs it crashed (the FFT buffer was shorter than 64 lags); on a day of 32 to 63 runs it divided
   wrapped-around products by zero or negative counts. Exploration days that short are common: Q-T1S1 has 164
   of 552 stock-days under 64 runs (minimum 1), Q-T1S2 62 of 736. So `run_primary` would have crashed on the
   real Q-T1 tick bin. The fix makes M3 undefined (NaN) at horizons above the day's length, where the real V
   is already undefined (no window of T runs exists). `_stock_mean` skips NaN days, so a stock's M3 at T comes
   from its days long enough for T, like its V. Nothing else in `primary3d` changed. Checked by
   `test_3d2.check_m3_short_days` with the old line as a planted bug; all 54 `test_3d1` tests still pass.
3. **P-M bins on real sizes (Measurement).** Round lots are 27 % to 51 % of runs. Quantile edges then
   coincide, and with ties to the lower bin a bin is empty on 31 % to 51 % of stock-days at 5 bins (3 % to
   11 % at 3 bins). Section 5.
4. **W-null is not a clean P-F null on real flow (Measurement, pilot).** In the pilot P-F rejects after Holm in
   8 of 10 panels in Q-T1 and in Q-T3 (pooled P-F slope -0.040 and -0.083, SD over panels 0.009 and 0.018). On
   real flow the informed runs keep their real signs, which carry the flow's sign memory, and with sizes
   independent of signs the cross terms weigh E[q^a] E[q^psi] against the diagonal's E[q^(a+psi)], so zeta_K(a)
   falls with a. The W-null family-wise error is therefore driven by P-F, not by the tests' calibration (section
   6).
5. **psi-hat is undefined on real sizes (Measurement).** Round-lot ties leave empty deciles: in a W0 panel every
   Q-T3 stock and 370 of 377 Q-T1 stocks had one, so psi-hat (secondary) was undefined in Q-T3 and rested on 6
   stocks in Q-T1 (section 6).
6. **Pilot (Measurement, 10 panels per world, default P-M design 3 bins, l_max 16, L 128, 50 runs; indicative
   only).**
   - P-V never rejects in any world, WM included: in WM M2 beats M1 in every panel but never M3, so the gate
     fails, and the (ii) slopes are positive (+0.02 S1, +0.03 S2), against the negative slope P-V predicts. In
     the semi-synthetic worlds P-V depends only on the real flow (V uses signs and sizes alone).
   - a* on real flow sits about 0.1 above psi: 0.58 to 0.60 in W0 (psi = 0.5) and 0.37 to 0.43 in W0-psi3 (psi =
     0.3). W0's P-R passes after Holm in 1 and 2 of 10 panels (Q-T1, Q-T3) and in none in Q-T2. W0-psi3 never
     passes, which is its correct outcome, but its a* lies inside the equivalence band, so that outcome rests on
     interval width, not on location.
   - WB's P-M passes after Holm in 4 of 10 panels in Q-T1 and 8 of 10 in Q-T3 at the default design; P-M is not
     evaluable in Q-T2 at L = 128 (no 5 stocks reach L under the ceilings).
   - In WM, a* is 0.67 (T1) and 0.86 to 0.89 (T2), outside the band, and P-M's top-bin beta is +0.22 to +0.24
     (Paper II clips beta_q at 0), so P-R and P-M cannot pass there.

## 1. Panels

**Composition.** Per confirmatory cell, the cell's exploration stocks are resampled with replacement up to its
confirmation count (Q-T1S1 161 from 69, Q-T1S2 216 from 92, Q-T2S3 17 from 7, Q-T3S1 50 from 22, Q-T3S2 27
from 11, Q-T3S3 23 from 10). Each drawn stock keeps its real 8 days of signs, raw sizes and timestamps. Tick
bins pool their cells through `primary3d.tick_bin_unit`, each stock at its own cell's ceiling: Q-T1 has 377
stocks, Q-T2 17, Q-T3 100.

**Synthetic response on independent stationary days** (3c-4, burn-in B = 20,000). Each stock-day is the end
of its own series of B + N runs, with the first B discarded.
- **Choice:** on real flow the B runs come from the same stock's other seven days, concatenated in calendar
  order, as one cyclic stretch with a random start. The history then has the stock's own size and sign
  statistics, and is independent of the day (3c-4 measured the cross-day mean-sign correlation at about 0). A
  copy of the same day (the 3c-3 semi-synthetic path) would correlate the history's impact with the day's
  own signs.
- Under GF (permanent impact) the burn-in only adds a constant to the day's mid, which no statistic sees, so
  GF skips it (exact).
- A stock whose other days are short repeats them cyclically to fill B runs.

**Seeds.** Panel i of world w and tick bin tb uses `default_rng([20261031, 0, w, tb, i])`. Tick bins of one
panel are separate draws of stocks; the Holm family of panel i joins the three tick bins with the same i.

## 2. Worlds

| world | flow | response | correct outcome |
|---|---|---|---|
| W-null | real, sizes permuted within each stock-day (`surrogates.s1`) | GF: a random 10 % of runs informed, permanent impact eps q^0.5 / E[q^0.5], plus a unit random walk | P-V no reject, P-F no reject, P-M no reject; P-R reported as observed |
| W0 | real | G0: one kernel (1 + l)^-beta, beta = `gen3.g0_beta(cell)`, impact eps q^0.5 | P-M no reject; P-R pass (Proposition 1: a* = psi = 0.5) |
| W0-psi3 | real | as W0 with psi = 0.3 | P-R fail (a* near 0.3) |
| WB | real | beta falls linearly in log q from 0.6 to -0.3 | P-M reject (both conditions) |
| WM | `mbmodel.PAPER2`, its own flow; day lengths of exploration stocks of each cell | the model's own mid | none imposed |

- "reject" means Holm rejects AND every condition holds (`passed`); "pass" for P-R is the same (the
  equivalence is established); "fail" is its negation; "no reject" means Holm does not reject.
- P-V power is read in WM and P-V size in W-null.
- **WB (Choice).** Per stock, over its 8 days, m1 and m5 are the median log size of its bottom and top
  quintile (raw sizes). beta(q) = 0.6 at log q = m1, -0.3 at m5, linear in between, clipped outside. The
  impact sum over per-run betas uses `mbmodel`'s Chebyshev interpolation in beta (20 nodes), one batched
  convolution per node; it matches a direct sum to 1e-7 (`check_wb_beta_and_impact`).
- **WM (Choice).** Each stock-day is the first N trades after B = 20,000 trades of an mbmodel day; the
  burn-in trades' impact stays in the mid. Timestamps spread the day's model time linearly over 9:30 to
  16:00 ET (for S2's 30-minute bins). Sizes are stored as float32.
- The pre-registered GF caveat (panel3c): on real flow, informed runs keep their real signs, which are
  correlated with the rest of the flow, so GF's exact zeta_K(a) = 1 does not carry over. See section 6.

## 3. How the tests are run

- `phase3d2.run_tick_bin` does what `run_primary` does for one tick bin: the unit, evaluability by the reach
  rule alone, then `run_unit` on the evaluable range tests only. It also runs P-M at every requested design
  through `primary3d._pm_test`, on the same draws as `run_unit` (the generator state is copied before
  `run_unit`, whose first draws are `P.draws`).
- `analyse3d2.family` then runs `primary3d.combine` (one Holm family of up to 12 pairs) at a chosen design.
- `check_tick_bin_matches_run_primary` shows that this gives exactly `run_primary`'s p-values (including
  P-M's p_slope) and Holm decisions at four designs, two of which differ only in min runs and share one reused
  result (designs that keep exactly the same stock-days have the same moments and draws, so the same result).
- Every `primary3d` function is used unchanged, except the `m3_V` fix (headline 2).
- **Rates (Choice).** "Before Holm" = evaluable and p <= 0.05. "After Holm" = rejected in the panel's family.
  A pair that left the family counts as not rejected. Rates are over all panels, with MC SE
  sqrt(r (1 - r) / n). The family-wise error in W-null is the share of panels in which any of P-V, P-F, P-M is
  rejected after Holm in any tick bin; the share including P-R is reported beside it.

## 4. The P-M design rule (stated before any study result)

**Trimmed on 1 Nov 2026 to 6 designs (min runs per bin fixed at 50), for compute cost only; see
`docs/phase3d_tests.md` 7.6.** The rule below is unchanged and now runs over those 6.

Grid: bins {3, 4, 5} x l_max {16, 32} x min runs per bin {20, 50, 100}, with **L = 8 l_max** (the smallest L
the pre-registration allows; Choice). The rule is coded in `worlds3d2.choose_design` and was committed to
before the study:
- **rate** of a design in a world and tick bin = the share of panels in which P-M is evaluable, has
  p <= 0.05 and both conditions hold (unadjusted; Choice: the same definition for size and power, as the
  instruction uses one word for both);
- **admissible**: in every tick bin, in W0 and in W-null, rate <= 0.05 + 2 sqrt(0.05 x 0.95 / panels)
  (0.081 at 200 panels; Choice: the MC SE at the nominal level, so the bound does not loosen with the
  observed rate);
- **objective**: the WB rate averaged over the three tick bins (a tick bin where P-M is not evaluable
  contributes 0, so reach counts);
- **ties**: fewer bins, then smaller l_max, then smaller min runs;
- if no design is admissible, none is chosen; the design with the smallest worst null rate is used for W0-psi3,
  WM and thinning, and flagged as a fallback.

Stage A runs W-null, W0 and WB at all 18 designs; the rule then fixes the design; stage B runs W0-psi3 and WM
at that design only.

## 5. P-M bin occupancy on real sizes (Measurement)

Exploration stock-days, R2-clipped sizes, each stock-day's own quantile edges, ties to the lower bin (exactly
`primary3d.pm_moments`). "Empty" = some bin holds no run. The last three columns are the share of stock-days
whose smallest bin holds at least 20 / 50 / 100 runs (that is, usable under each min-runs value, before
the ceiling rule).

| cell | bins | stock-days | round-lot share | median lowest-bin share (nominal) | empty bin | min bin >= 20 | >= 50 | >= 100 |
|---|---|---|---|---|---|---|---|---|
| Q-T1S1 | 3 | 552 | 0.30 | 0.34 (0.33) | 0.11 | 0.57 | 0.34 | 0.20 |
| Q-T1S1 | 4 | 552 | 0.30 | 0.26 (0.25) | 0.29 | 0.39 | 0.19 | 0.11 |
| Q-T1S1 | 5 | 552 | 0.30 | 0.21 (0.20) | 0.51 | 0.21 | 0.10 | 0.05 |
| Q-T1S2 | 3 | 736 | 0.27 | 0.34 (0.33) | 0.03 | 0.84 | 0.70 | 0.55 |
| Q-T1S2 | 4 | 736 | 0.27 | 0.26 (0.25) | 0.15 | 0.67 | 0.53 | 0.38 |
| Q-T1S2 | 5 | 736 | 0.27 | 0.21 (0.20) | 0.32 | 0.50 | 0.35 | 0.24 |
| Q-T2S3 | 3 | 56 | 0.38 | 0.34 (0.33) | 0.02 | 0.88 | 0.73 | 0.52 |
| Q-T2S3 | 4 | 56 | 0.38 | 0.26 (0.25) | 0.14 | 0.71 | 0.57 | 0.29 |
| Q-T2S3 | 5 | 56 | 0.38 | 0.21 (0.20) | 0.41 | 0.43 | 0.34 | 0.12 |
| Q-T3S1 | 3 | 176 | 0.49 | 0.35 (0.33) | 0.05 | 0.92 | 0.84 | 0.73 |
| Q-T3S1 | 4 | 176 | 0.49 | 0.25 (0.25) | 0.22 | 0.72 | 0.61 | 0.48 |
| Q-T3S1 | 5 | 176 | 0.49 | 0.20 (0.20) | 0.48 | 0.42 | 0.32 | 0.22 |
| Q-T3S2 | 3 | 88 | 0.51 | 0.41 (0.33) | 0.06 | 0.91 | 0.77 | 0.64 |
| Q-T3S2 | 4 | 88 | 0.51 | 0.25 (0.25) | 0.17 | 0.74 | 0.61 | 0.40 |
| Q-T3S2 | 5 | 88 | 0.51 | 0.20 (0.20) | 0.44 | 0.41 | 0.25 | 0.15 |
| Q-T3S3 | 3 | 80 | 0.49 | 0.41 (0.33) | 0.03 | 0.91 | 0.75 | 0.64 |
| Q-T3S3 | 4 | 80 | 0.49 | 0.26 (0.25) | 0.12 | 0.78 | 0.68 | 0.47 |
| Q-T3S3 | 5 | 80 | 0.49 | 0.20 (0.20) | 0.31 | 0.54 | 0.39 | 0.28 |

- Round lots (multiples of 100 shares) are 27 % to 51 % of runs.
- The lowest bin holds about its nominal share in the median (ties to the lower bin push it up to 0.41 at 3
  bins in Q-T3S2 and Q-T3S3), but an empty bin is frequent at 4 and 5 bins, so the min-runs rule removes those
  stock-days at any threshold.
- Q-T1S1 is the binding cell: at 5 bins only 10 % of its stock-days have every bin with 50 runs or more.

## 6. Pilot (step 1)

**Timing per tick bin** (mean over 10 panels, seconds on one worker, one P-M design; two single-job measurements
ran alongside part of the pilot, so these are slight overestimates):

| world | tick bin | generate s | unit s | P-V/P-R/P-F s | P-M (1 design) s | total s |
|---|---|---|---|---|---|---|
| W-null | Q-T1 | 0.2 | 1.8 | 20.1 | 1.8 | 23.9 |
| W-null | Q-T2 | 0.0 | 0.0 | 0.7 | 0.0 | 0.7 |
| W-null | Q-T3 | 0.1 | 0.7 | 7.4 | 1.1 | 9.3 |
| W0 | Q-T1 | 3.0 | 1.9 | 20.6 | 1.9 | 27.3 |
| W0 | Q-T2 | 0.1 | 0.0 | 0.6 | 0.0 | 0.8 |
| W0 | Q-T3 | 0.8 | 0.7 | 7.2 | 1.1 | 9.8 |
| W0-psi3 | Q-T1 | 3.5 | 1.8 | 20.2 | 1.9 | 27.3 |
| W0-psi3 | Q-T2 | 0.2 | 0.0 | 0.7 | 0.0 | 0.9 |
| W0-psi3 | Q-T3 | 1.0 | 0.7 | 7.4 | 1.2 | 10.3 |
| WB | Q-T1 | 42.0 | 2.1 | 22.6 | 2.1 | 68.8 |
| WB | Q-T2 | 1.5 | 0.0 | 0.7 | 0.0 | 2.2 |
| WB | Q-T3 | 9.1 | 0.7 | 7.4 | 1.2 | 18.4 |
| WM | Q-T1 | 101.8 | 1.8 | 19.2 | 1.8 | 124.6 |
| WM | Q-T2 | 4.6 | 0.0 | 0.6 | 0.0 | 5.3 |
| WM | Q-T3 | 29.2 | 0.7 | 6.9 | 1.1 | 37.9 |

The pilot's own projection (19.3 worker-hours for stages A and B) prices each of the 18 P-M designs at the
default design's cost, which is the cheapest. Measured on single jobs, the 18 designs take 50 s (Q-T1) and 66 s
(Q-T3), not 18 x 1.9 s and 18 x 1.1 s (section 8).

**Rates** (before Holm = evaluable and p <= 0.05; after Holm = the panel's one family; MC SE in brackets; 10
panels, so an MC SE is up to 0.16):

| world | test | tick bin | evaluable | before Holm | after Holm | passed | correct | rate correct |
|---|---|---|---|---|---|---|---|---|
| W-null | P-V | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W-null | P-V | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W-null | P-V | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W-null | P-R | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W-null | P-R | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W-null | P-R | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W-null | P-M | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W-null | P-M | Q-T2 | 0.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W-null | P-M | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W-null | P-F | Q-T1 | 1.0 | 0.900 (0.095) | 0.800 (0.126) | 0.800 (0.126) | no reject | 0.200 (0.126) |
| W-null | P-F | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W-null | P-F | Q-T3 | 1.0 | 1.000 (0.000) | 0.800 (0.126) | 0.800 (0.126) | no reject | 0.200 (0.126) |
| W0 | P-V | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0 | P-V | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0 | P-V | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0 | P-R | Q-T1 | 1.0 | 0.600 (0.155) | 0.100 (0.095) | 0.100 (0.095) | pass | 0.100 (0.095) |
| W0 | P-R | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | pass | 0.000 (0.000) |
| W0 | P-R | Q-T3 | 1.0 | 0.700 (0.145) | 0.200 (0.126) | 0.200 (0.126) | pass | 0.200 (0.126) |
| W0 | P-M | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W0 | P-M | Q-T2 | 0.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W0 | P-M | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | no reject | 1.000 (0.000) |
| W0 | P-F | Q-T1 | 1.0 | 0.300 (0.145) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0 | P-F | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0 | P-F | Q-T3 | 1.0 | 1.000 (0.000) | 0.800 (0.126) | 0.800 (0.126) | none | - |
| W0-psi3 | P-V | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0-psi3 | P-V | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0-psi3 | P-V | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0-psi3 | P-R | Q-T1 | 1.0 | 0.100 (0.095) | 0.000 (0.000) | 0.000 (0.000) | fail | 1.000 (0.000) |
| W0-psi3 | P-R | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | fail | 1.000 (0.000) |
| W0-psi3 | P-R | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | fail | 1.000 (0.000) |
| W0-psi3 | P-M | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0-psi3 | P-M | Q-T2 | 0.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0-psi3 | P-M | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0-psi3 | P-F | Q-T1 | 1.0 | 0.100 (0.095) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0-psi3 | P-F | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| W0-psi3 | P-F | Q-T3 | 1.0 | 1.000 (0.000) | 1.000 (0.000) | 1.000 (0.000) | none | - |
| WB | P-V | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WB | P-V | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WB | P-V | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WB | P-R | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WB | P-R | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WB | P-R | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WB | P-M | Q-T1 | 1.0 | 0.600 (0.155) | 0.400 (0.155) | 0.400 (0.155) | reject | 0.400 (0.155) |
| WB | P-M | Q-T2 | 0.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | reject | 0.000 (0.000) |
| WB | P-M | Q-T3 | 1.0 | 1.000 (0.000) | 0.800 (0.126) | 0.800 (0.126) | reject | 0.800 (0.126) |
| WB | P-F | Q-T1 | 1.0 | 0.200 (0.126) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WB | P-F | Q-T2 | 1.0 | 0.400 (0.155) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WB | P-F | Q-T3 | 1.0 | 0.700 (0.145) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-V | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-V | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-V | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-R | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-R | Q-T2 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-R | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-M | Q-T1 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-M | Q-T2 | 0.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-M | Q-T3 | 1.0 | 0.000 (0.000) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-F | Q-T1 | 1.0 | 0.300 (0.145) | 0.100 (0.095) | 0.100 (0.095) | none | - |
| WM | P-F | Q-T2 | 1.0 | 0.100 (0.095) | 0.000 (0.000) | 0.000 (0.000) | none | - |
| WM | P-F | Q-T3 | 1.0 | 0.300 (0.145) | 0.100 (0.095) | 0.100 (0.095) | none | - |

**Statistics** (mean over the 10 panels; SD over panels in brackets):

| world | tick bin | P-V slope S1 | P-V slope S2 | a*(T1) | a*(T2) | psi-hat | P-F slope | P-M slope | P-M top beta |
|---|---|---|---|---|---|---|---|---|---|
| W-null | Q-T1 | +0.001 (0.001) | +0.001 (0.002) | 0.523 (0.024) | 0.506 (0.064) | 0.484 (0.070) | -0.040 (0.009) | -0.08 (0.21) | -0.36 (0.42) |
| W-null | Q-T3 | +0.000 (0.002) | -0.000 (0.004) | 0.506 (0.032) | 0.473 (0.076) | - | -0.083 (0.018) | -0.08 (0.40) | -0.11 (0.56) |
| W0 | Q-T1 | +0.045 (0.004) | +0.049 (0.005) | 0.597 (0.007) | 0.582 (0.023) | 0.655 (0.119) | -0.015 (0.005) | -0.026 (0.023) | +0.366 (0.017) |
| W0 | Q-T2 | +0.062 (0.026) | +0.067 (0.019) | 0.602 (0.012) | 0.616 (0.076) | - | +0.014 (0.030) | - | - |
| W0 | Q-T3 | -0.035 (0.005) | -0.017 (0.004) | 0.604 (0.007) | 0.557 (0.021) | - | -0.097 (0.013) | -0.003 (0.017) | +0.419 (0.022) |
| W0-psi3 | Q-T1 | +0.046 (0.005) | +0.048 (0.004) | 0.430 (0.006) | 0.401 (0.033) | 0.367 (0.026) | -0.011 (0.005) | -0.008 (0.009) | +0.382 (0.023) |
| W0-psi3 | Q-T3 | -0.034 (0.010) | -0.017 (0.009) | 0.422 (0.009) | 0.374 (0.015) | - | -0.119 (0.015) | -0.004 (0.021) | +0.415 (0.024) |
| WB | Q-T1 | +0.043 (0.004) | +0.047 (0.003) | 0.641 (0.006) | 0.735 (0.012) | 0.446 (0.046) | +0.015 (0.005) | -0.464 (0.118) | -0.326 (0.008) |
| WB | Q-T3 | -0.036 (0.006) | -0.019 (0.005) | 0.735 (0.007) | 0.987 (0.027) | - | +0.062 (0.006) | -0.255 (0.028) | -0.288 (0.009) |
| WM | Q-T1 | +0.021 (0.006) | +0.036 (0.005) | 0.672 (0.002) | 0.858 (0.010) | 0.345 (0.005) | -0.009 (0.005) | -0.113 (0.009) | +0.237 (0.014) |
| WM | Q-T3 | +0.013 (0.005) | +0.027 (0.005) | 0.672 (0.003) | 0.874 (0.007) | 0.347 (0.009) | -0.012 (0.004) | -0.119 (0.012) | +0.223 (0.021) |

psi-hat is '-' where it is undefined: no stock has a slope. Q-T2 rows are omitted. **Measurement:** psi-hat's
deciles of run size are per stock over its 8 days with ties to the lower bin; with round lots, decile edges
coincide and some deciles are empty, and a stock with an empty (or non-positive) decile has no slope (the 3d-1
Choice). In one W0 panel, 100 of 100 Q-T3 stocks and 370 of 377 Q-T1 stocks had an empty decile: psi-hat was
undefined in Q-T3 and rested on 6 stocks in Q-T1. It is secondary (not in Holm), but as defined it does not
measure psi on real sizes.

**Measurement (pilot).**
- W-null: P-V, P-M never reject; P-F rejects in 8 of 10 panels in Q-T1 and Q-T3 (headline 4). The pilot's
  family-wise error in W-null is therefore at least 0.8.
- P-V's gate (i) holds only in Q-T3 of the semi-synthetic worlds (M2 beats M1 everywhere, M3 only there); (ii)
  there has S1 slope -0.035 and S2 -0.017, and the S2 upper bound is never below 0, so P-V never rejects. In
  Q-T1 and Q-T2 the slopes are positive. These are properties of the resampled exploration flow (signs and sizes
  only), the same in W0, W0-psi3 and WB up to the draw of stocks.
- WM (the paper's world): P-V's gate fails and its slopes are positive; a* is 0.67 to 0.89; P-M's top-bin beta
  is positive. Under the pre-registered tests the paper's own model would pass none of P-V, P-R, P-M.
- a* bias on real flow under G0: about +0.1 at both psi (W0: 0.58 to 0.60; W0-psi3: 0.37 to 0.43). psi-hat is
  close to psi in W0-psi3 (0.367) and above it in W0 (0.655, SD 0.12).

## 7. Thinning (step 5): design, not run

- WM flow at observed shares 1.0, 0.46, 0.26, 0.17; random rule, and router rule with hidden share h = 0.17.
- **Kappa (Choice, needs confirmation).** `thinning.py` has no kappa grid, only the default concentration
  kappa = 4. The driver uses kappa in {1, 4, 16}: 13 settings.
- **Share applied once.** The router's observed share is p (1 - h), so p = share / (1 - h) (`router_p`); h is
  applied only inside `thinning.router_thinning`. Share 1.0 is the unthinned market, not a router setting
  (p = 1 / 0.83 > 1 is not reachable). Measured observed shares over 40 seeds at kappa 16: 0.4598 for 0.46
  and 0.1694 for 0.17, seed SD 0.003.
- **Paired (Choice).** The observed day keeps the cell's day length N, so the market day behind it holds about
  N / share trades. Every setting thins the same market days: each setting has its own thinning seed, and every
  job thins all 13 settings to make the same market-redraw decisions, keeping only its own group (4 groups,
  memory). Changes are paired differences against share 1.0 on the same panels.
- 50 panels per setting; P-M at the chosen design; each setting is its own Holm family across tick bins.

## 8. Running the study

**Measured cost of the real jobs** (one worker, after the memory fixes; peak RSS of the job's process):

| job | seconds | peak RSS |
|---|---|---|
| stage A, W0, Q-T1, 18 designs | 75 (generate 3, range tests 21, P-M 50) | 0.67 GB |
| stage A, W0, Q-T3, 18 designs | 73 (P-M 66) | 0.48 GB |
| stage A, WB, Q-T1, 18 designs | 90 (generate 23, P-M 50) | 0.75 GB |
| stage B, WM, Q-T1, 1 design | 111 (generate 93) | 0.64 GB |
| thinning, Q-T1, one group of settings | 230 to 250 (generate 160) | 1.03 GB |
| thinning, Q-T3, one group | 83 (generate 51) | |

Totals: stage A about 21 worker-hours (W-null and W0 200 panels at about 145 s, WB 100 at about 170 s), stage B
about 6 (W0-psi3 100 at 39 s, WM 100 at 168 s), thinning about 18 (50 panels x 4 groups; each group regenerates
the same market days). About 45 worker-hours: about 9 h for stages A and B on 3 workers, then about 8 h of
thinning (its Q-T1 jobs at most 2 at a time, `THIN_T1_3`, since each peaks at 1.0 GB). Memory: at most 3 x 0.75
GB in stages A and B, 2 x 1.03 GB in Q-T1 thinning.

**Memory and speed changes found while measuring (no effect on results):**
- WB generates 8 stock-days at a time instead of 64 (its Chebyshev weights held 20 copies of a group's series:
  0.99 GB, now 0.58 GB).
- Designs that differ only in min runs per bin and keep exactly the same stock-days share one P-M result (same
  moments, same draws; exact). A smaller P-M solve chunk was tried and dropped: no memory saving and 1.6 times
  slower.

One line from Terminal (stage A, the design rule, stage B, thinning, analysis; every step resumes where it
stopped):

    cd ~/CODE/itchbook && caffeinate -i uv run --group analysis tools/est/phase3d2.py all --workers 3 2>&1 | tee -a ~/CODE/itch-data/phase3d2/run.log

Outputs land in `~/CODE/itch-data/phase3d2/report/`: `design.json` (the rule's inputs and choice),
`tables.md` (every (test, tick bin, world) rate before and after Holm with MC SE, FWER in W-null, ranges used,
pairs that left, statistics, occupancy, thinning) and `results.json`. `phase3d2.py analyse` rebuilds them from
the saved chunks; `phase3d2.py topup --panels N` adds panels to every world. `all --panels N` caps every
world's panels and the thinning panels at N, and `PHASE3D2_OUT` moves the output directory.

**Smoke run of the whole chain (Measurement).** `PHASE3D2_OUT=<scratch> phase3d2.py all --panels 2 --workers 3`
ran stage A (8.4 min), the design rule, stage B (4.1 min), thinning (Q-T1 16.6 min on 2 workers, the rest
5.2 min) and `analyse` in 34 min, peak 1.04 GB in one process. Its tables were checked by eye (design rule,
family, FWER, ranges, pairs that left, thinning changes). Two panels say nothing about rates; the smoke run's
outputs are not in the repo or in `~/CODE/itch-data/phase3d2`. The per-job times agree with section 8's
projection (a Q-T1 thinning group about 250 s).

## 9. Tests and planted bugs

`tools/est/test_3d2.py`: 14 checks, 17 planted function bugs (`mutate` fixture) and 3 planted label bugs
(`monkeypatch.setitem` on `worlds3d2.CORRECT`): 34 tests. Each function mutant was also run outside pytest:
16 fail their check on an assertion, and the m3_V mutant fails on the original crash.

| planted bug | mutated | caught by |
|---|---|---|
| design rule's objective read from W0 | `choose_design` | check_design_rule |
| design rule's size constraint read from WB | `choose_design` | check_design_rule |
| design rule's size constraint on W-null only | `choose_design` | check_design_rule |
| design ties to more bins | `choose_design` | check_design_rule_ties_and_fallback |
| hidden share applied twice (router_p ignores h) | `router_p` | check_router_share_once |
| hidden share applied twice (thin_day) | `thin_day` | check_router_share_once |
| only the kept settings thinned, from the market generator | `wm_cell_flows` | check_thin_settings_paired |
| pass and fail swapped in the P-R label | `correct_rate` | check_world_labels |
| reject counted without the conditions | `correct_rate` | check_world_labels |
| W0-psi3 generated at psi = 0.5 | `world_spec` | check_world_specs |
| burn-in includes the day itself | `burn_index` | check_burn_from_other_days |
| WB anchors reversed | `wb_beta` | check_wb_beta_and_impact |
| burn-in kept in the day's mid | `semi_response` | check_semi_response_direct |
| P-M draws not run_unit's | `run_tick_bin` | check_tick_bin_matches_run_primary |
| P-M result reused across designs whatever stock-days they keep | `run_tick_bin` | check_tick_bin_matches_run_primary |
| W-null sizes shuffled across days | `shuffle_sizes` | check_semi_flow_unchanged |
| m3_V ACF to lag 63 whatever the day length | `primary3d.m3_V` | check_m3_short_days |
| WB's P-M labelled no reject | `CORRECT` | check_world_labels |
| W-null's P-F unlabelled | `CORRECT` | check_world_labels |
| W0's P-R labelled fail | `CORRECT` | check_world_labels |

Other checks: the world specs; the router and random observed shares; the setting grid; the G0 and WB
responses against a direct sum over [burn-in, day] with the same draws; the WM burn-in equal to the full
simulation's later trades; real flow passed through unchanged.

Full suite (`uv run --group dev --group analysis pytest`): 432 passed, 28 skipped (was 398 passed before Phase
3d-2). The new file runs in about 9 s.

## 10. Open (not decided here)

- **The m3_V fix** changes a frozen function; it needs your sign-off (headline 2).
- **Kappa grid** for thinning (section 7).
- **W-null for P-F** (section 6): on real signs GF's informed runs carry the flow's sign memory, so zeta_K(a)
  is not flat. A clean P-F null on real flow would need the informed runs' signs replaced by fair coins (as in
  synthetic GF), which changes the flow; not done here.
- **psi-hat's bins on real sizes** (section 6): empty deciles from round-lot ties make it undefined for most
  stocks; fewer bins or ties split across bins would be a change to its definition.
- **P-V, P-R and P-M in WM** (pilot): the paper's own model passes none of them under the pre-registered
  tests; the full study will measure it with 100 panels.
- **P-M at 4 and 5 bins on real sizes** loses many stock-days to empty bins (section 5); the design rule
  weighs this only through evaluability and power.
