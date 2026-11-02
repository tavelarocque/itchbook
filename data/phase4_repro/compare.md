# Phase 4 reproduction check

Result: **bit-identical**. 2823 fields compared, 0 mismatches.

A reproduction of the confirmatory run, not a second confirmatory analysis; it changes no reported result. tools/est/phase4_repro.py loaded the confirmation stocks with phase4.build_flows('confirm') and called the frozen primary3d.run_primary once over the nine NASDAQ-listed cells with the frozen P-M design and seed, with no wrapper, then compared every field it returns with data/phase4/results.json (floats by float.hex).

- git HEAD: 915cbd073386701eb383e96128a958c383f80ea7 (stored run: d126088fba041824950d7ffae4aa75a7a98b2703); working tree: ?? tools/est/phase4_repro.py
- tools/est/phase4_repro.py SHA-256 2ee4f35cc2ea0eb4c078a8d3e0385083255010ca5a6ac1109e645ca02975bf28
- seed 20261102; P-M design {'n_bins': 3, 'l_max': 16, 'L': 128, 'min_runs_per_bin': 50}; B = 200; surrogate draws 4
- python 3.13.16, numpy 2.5.4
- runtime: 55.2 s (loading and hashing 6.9 s, run_primary 48.3 s); the stored run took 39.8 s

## Fields compared

| section | fields |
|---|---|
| result | 1790 |
| run | 7 |
| inputs_sha256 | 8 |
| inputs_cells | 921 |
| inputs_plan_counts | 9 |
| holm_recomputed | 48 |
| holm_recomputed_vs_reproduced | 40 |

result = every leaf of run_primary's return (evaluability, reach, p-values, Holm, sub-conditions, estimates, intervals, P-R secondaries, descriptive P-K); run = seed, P-M design, B, surrogate draws; inputs_* = runs files' SHA-256, per-cell symbols and counts, plan stock counts; holm_recomputed* = frozen primary3d.holm on the reproduced family p-values against the stored and reproduced records, and passed. The 'methods' record (from phase4's observe() recorder) is not a field of run_primary and is not compared.

## Holm recomputed from the reproduced p-values

| rank | pair | p | threshold | adjusted p | rejects | passed |
|---|---|---|---|---|---|---|
| 1 | P-F|Q-T1 | 0.0001428966406721203 | 0.00625 | 0.0011431731253769624 | yes | yes |
| 2 | P-F|Q-T3 | 0.004780873334061808 | 0.0071428571428571435 | 0.033466113338432656 | yes | yes |
| 3 | P-R|Q-T3 | 0.445 | 0.008333333333333333 | 1.0 | no | no |
| 4 | P-V|Q-T3 | 0.97 | 0.01 | 1.0 | no | no |
| 5 | P-V|Q-T1 | 1.0 | 0.0125 | 1.0 | no | no |
| 6 | P-R|Q-T1 | 1.0 | 0.016666666666666666 | 1.0 | no | no |
| 7 | P-M|Q-T1 | 1.0 | 0.025 | 1.0 | no | no |
| 8 | P-M|Q-T3 | 1.0 | 0.05 | 1.0 | no | no |
