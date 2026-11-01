# Phase 3d-2c: the decisions of 1 Nov 2026, the a* rerun and the re-pilot

Labels: **Measurement** = what a computation showed; **Choice** = a decision taken here. Real data were read only
through `worlds3d2`'s flow cache (exploration stocks, signs, sizes and timestamps). No replication or
confirmation file was opened. Every price response is synthetic. The full 3d-2 study was **not** started.

Nothing is frozen until the end of Phase 3. The changes below are ordinary pre-freeze decisions, recorded in
`docs/phase3d_tests.md` section 7 (7.3 to 7.5).

## Headlines

1. **a* is now the local vertex at the grid maximum (decision 1).** Without sign memory (W0-shuf, psi = 0.5)
   a*(1) is 0.534 to 0.540, against 0.570 to 0.583 with the old global quadratic. That is closer to 0.5 but
   not at it: +0.03 to +0.04 remains. The sizes-only curve predicts +0.02 to +0.03 of that from the
   three-point stencil on a skewed curve.
2. **Correction to 3d-2b.** The old a* fitted over P-R's grid [0, 1], not [0, 1.5], as 3d-2b said. The
   sizes-only vertex on [0, 1] (0.556 to 0.579) matches the panels, so 3d-2b's "unexplained gap" of 0.03 to
   0.05 does not exist. `docs/phase3d2b_report.md` now carries the correction.
3. **WM stays PAPER2, and its P-M fails by construction (decision 2).** `all` runs W-null (clean), W0,
   W0-psi3 and WB only (decision 3). One consequence to note: P-V has no power world in `all`, because WM
   was it.
4. **Re-pilot and smoke run (Measurement).** The `all` study should take about 11 to 13 h wall on 3 workers.
   The pilot's own projection (18.5 worker-hours, 6.2 h) is too low: it prices each of the 18 P-M designs at
   the default design's cost. A 1-panel smoke run of `all` measured P-M at all 18 designs at 88 to 111 s per
   Q-T1 or Q-T3 job. The pilot projection assumed 27 to 49 s.
5. **psi-hat (decision 4).** No correction. `primary3d.PSI_NOTE` ("+0.02 to +0.05 high without sign memory")
   is printed next to every psi-hat the code reports.

## 1. a* (decision 1; Choice, docs/phase3d_tests.md 7.3)

**Rule.** a* = a_k + (h / 2) (R_k-1 - R_k+1) / (R_k-1 - 2 R_k + R_k+1). Here a_k is the grid maximum of the
pooled Fisher-z R_a(T) and h = 0.25.
- If the maximum is at a = 0 or a = 1, a* is censored at that end and P-R fails at that T (p = 1, the
  condition False).
- Censored replicates keep the end value, which lies outside the margin.
- `primary3d.vertex` and `primary3d.censored`; P-R reads a* only through them.
- `pooled.core2` still returns the old a_star (now also `primary3d.vertex_ls`), so Phases 3c to 3d-2b
  reproduce.

**Sizes-only curve (Measurement, `diag3d2b.vertex`, 3d-2b's record; the local vertex is computed from its
stored shape).**

| tick bin | argmax (fine grid) | global quadratic on [0, 1] (old a*) | local vertex (new a*) |
|---|---|---|---|
| Q-T1 | 0.51 | 0.579 | 0.529 |
| Q-T2 | 0.51 | 0.572 | 0.525 |
| Q-T3 | 0.51 | 0.556 | 0.520 |

**The rerun (Measurement).** W0 against W0-shuf, 30 panels per tick bin. These are the same seeds, hence the
same panels, as 3d-2b. The old a* recomputed on them reproduces 3d-2b exactly (max difference 0 over 180
panels).

- Output: `~/CODE/itch-data/phase3d2/diag3d2b_astar/` (`tables.md`); 3.4 min on 3 workers.
- Command: `PHASE3D2B_OUT=~/CODE/itch-data/phase3d2/diag3d2b_astar uv run --group analysis diag3d2b.py abias`.

a*(T), new (local vertex), mean over panels; MC SEs at most 0.007 except Q-T2 at T >= 32:

| world | tick bin | T = 1 | 2 | 4 | 8 | 16 | 32 | 64 | P-R passes (unadjusted) |
|---|---|---|---|---|---|---|---|---|---|
| W0 | Q-T1 | 0.597 | 0.596 | 0.592 | 0.580 | 0.568 | 0.557 | 0.548 | 19/30 = 0.63 (0.09) |
| W0 | Q-T2 | 0.577 | 0.580 | 0.583 | 0.568 | 0.565 | 0.566 | 0.521 | 0/30 = 0.00 (0.00) |
| W0 | Q-T3 | 0.582 | 0.588 | 0.593 | 0.596 | 0.595 | 0.585 | 0.548 | 25/30 = 0.83 (0.07) |
| W0-shuf | Q-T1 | 0.540 | 0.548 | 0.560 | 0.572 | 0.590 | 0.620 | 0.652 | 0/30 = 0.00 (0.00) |
| W0-shuf | Q-T2 | 0.537 | 0.548 | 0.561 | 0.562 | 0.588 | 0.615 | 0.636 | 0/30 = 0.00 (0.00) |
| W0-shuf | Q-T3 | 0.534 | 0.541 | 0.551 | 0.567 | 0.587 | 0.611 | 0.624 | 1/30 = 0.03 (0.03) |

Old a* on the same panels:
- W0 a*(1): 0.631, 0.613, 0.614.
- W0-shuf a*(1): 0.583, 0.581, 0.570.
- W0's P-R pass rates in 3d-2b: 0.47, 0.00, 0.80.

Censored panels: none, except 2 of 30 high at T = 64 in W0-shuf Q-T2.

**Answer: a*(1) moves toward 0.5, but does not reach it.** Without sign memory it is now 0.53 to 0.54, about
0.04 closer than before. With real signs it is 0.58 to 0.60.

Not changed by the new estimator: without sign memory a* still rises with T, to 0.62 to 0.65 at T = 64, so
P-R fails in W0-shuf at T2. With real signs it falls with T. Neither drift was chased here.

## 2. WM (decision 2; Choice, 7.5)

WM is the published PAPER2 preset, with beta_q clipped at 0. In every WM table, WM's P-M shows
"fails by construction (PAPER2 clips beta_q at 0: the top bin cannot be < 0)" (`worlds3d2.BY_CONSTRUCTION`)
instead of a scored label. `W.CORRECT["WM"]` stays unlabelled, and WB stays the P-M power world.

## 3. The trimmed study and the re-pilot (decision 3)

`phase3d2.py`:
- `all` runs design, power, analyse:
  - design: W-null (clean coin-flip version), W0 and WB at all 18 P-M designs, tie-aware bins;
  - power: W0-psi3 at the chosen design;
  - analyse.
- WM is the new `wm` command and thinning is `thin`. Neither is in `all`.
- `topup` and the pilot cover the four study worlds.
- Seeds are unchanged (`W.WORLDS` keeps its order).

The 3d-2 pilot's outputs were moved, not deleted, to `~/CODE/itch-data/phase3d2/report/pilot_3d2/`.

**Re-pilot timing (Measurement; 10 panels per world, the default design only).** Mean seconds per panel and
tick bin:

| world | Q-T1 | Q-T2 | Q-T3 |
|---|---|---|---|
| W-null | 29.7 | 0.8 | 10.5 |
| W0 | 31.5 | 1.0 | 11.2 |
| WB | 60.9 | 2.4 | 21.9 |
| W0-psi3 | 30.5 | 0.9 | 11.8 |

The pilot itself took 12.0 min. Its printed projection, 18.5 worker-hours (6.2 h wall on 3 workers), assumes
the 17 extra P-M designs cost what the default design does, and that is not so.

**Smoke run of `all` (Measurement).** `PHASE3D2_OUT=<scratch> phase3d2.py all --panels 1`: design (9 jobs),
then power (3 jobs, W0-psi3 only), then analyse, in 3.3 min. No WM. The psi-hat note appears in every
psi-hat row. The design chosen from 1 panel means nothing.

Seconds per job, P-M at all 18 designs:

| world | Q-T1 P-M | Q-T3 P-M | Q-T1 total | Q-T3 total |
|---|---|---|---|---|
| W-null | 90 | 96 | 145 | 104 |
| W0 | 88 | 102 | 145 | 112 |
| WB | 57 | 111 | 167 | 136 |

**Estimates of the full study.**
- Smoke P-M costs with the pilot's other costs: about 33 worker-hours, 11 h wall on 3 workers.
- The smoke run's own job times: 37.7 worker-hours, 12.6 h. These are 1 panel each, and the range tests ran
  slower than in the pilot with three Q-T1 jobs at once.

**Expect 11 to 13 h.** The run is resumable by chunk (2 panels).

**Indicative pilot rates (Measurement, 10 panels per world, default P-M design b3_l16_L128_m50, after the one
Holm family; `report/pilot_rates.md`).** These are not study results.
- **W-null (clean).** Nothing is rejected, in any test or tick bin.
- **W0, P-R pass.** Q-T1: 0.8 before Holm, 0.1 after. Q-T2: 0.0 before, 0.0 after. Q-T3: 0.7 before, 0.2
  after. Holm, not a*, is what holds P-R back here.
- **W0, P-F.** Rejects in Q-T3 in 0.8 of panels after Holm. There is no label (real signs; cf. W-null-real in
  3d-2b).
- **W0-psi3, P-R.** Fails everywhere (correct).
- **WB, P-M pass.** Q-T1: 0.3 before Holm, 0.0 after. Q-T3: 0.8 before, 0.1 after. Q-T2 is not evaluable.

## 4. psi-hat (decision 4)

No change to the estimator. `primary3d.PSI_NOTE` is attached to P-R's secondary psi-hat, to every psi-hat row
of the analyse3d2 tables (and below the thinning table) and to the diag3d2b a* tables.

## 5. Tests (Measurement)

- New checks: `check_astar_local`, `check_astar_censored`, `check_pr_uses_local_vertex`.
- Six new planted bugs, all caught:
  - the global quadratic, twice (in vertex; P-R reading core2's a_star);
  - neighbours from one side;
  - neighbours swapped;
  - a censored end treated as interior, twice (in vertex; in the P-R decision).
- `test_3d1`:
  - `check_vertex_matches_core2` now checks `vertex_ls`.
  - `astar_reps` gives its 8 jackknife replicates noise sd / sqrt(8), as `reps` does, instead of the full sd.
    With the full sd, the local vertex (3 of 5 points, noisier) gave P-R p = 0.040 on the Eq. 74 known
    answer, against the check's 0.01. With the fix it gives 0.0002 and 0.0004.
- `test_3d2` now restores prereg3d's grid after every test. The planted "grid patch not restored" bug used to
  leak a patched C.T_GRID into later tests: 13 test_3d1 failures when test_3d2 ran first. This was there
  before this step.

Full suite: 451 passed, 28 skipped.

## 6. Open (none acted on)

1. **P-V has no power world in `all`.** WM was it (`worlds3d2.PV_POWER_WORLD`).
2. **a*'s drift with T is unexplained.** Without sign memory it rises (to about 0.63 at T = 64); with real
   signs it falls.
3. **W0's P-R rate after Holm is low in the pilot** (0.1, 0.0, 0.2). The 200-panel study will measure it.
4. **The pilot's projection formula** prices every P-M design as the default one. It was left as it is; the
   smoke-run figures above replace it.
