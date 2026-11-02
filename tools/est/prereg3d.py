"""Phase 3d pre-registered constants: one place for every value the primary tests take as fixed.

Fixed by the pre-registration (implemented, not revisited here):
- confirmatory cells and their ceilings rho (the Phase 3c-4 rho* per cell); Q-T1S3, Q-T2S1, Q-T2S2 are
  descriptive only and never enter a primary test;
- the pooling unit of the primary tests is the tick bin; per-cell results are secondary;
- a grid {0, 0.25, 0.5, 0.75, 1}; fitting ranges [4, 32] and [8, 64], [8, 64] primary where evaluable under
  the reach rule (pooled.MIN_REACH stocks at T2), else [4, 32];
- R2 clip: each stock's sizes clipped at its own 99.9th percentile over its study days;
- Holm-Bonferroni at family-wise 0.05 over exactly P-V, P-R, P-M, P-F; P-K and P-C descriptive.
Decided 31 Oct 2026 (docs/phase3d_tests.md): P-R is the equivalence test a*(T) in 1/2 +- 0.15 at both ends
of the range used.
Decided 2 Nov 2026, before any real-data test (docs/phase3d_tests.md 7.7): Q-T2 leaves the confirmatory
family and is descriptive. The Holm family is the four tests in Q-T1 and Q-T3 only (FAMILY_TICK_BINS, at
most 8 p-values). Q-T2's unit (its one cell Q-T2S3, which keeps its ceiling) is still computed and reported,
and never enters Holm.
"""
from __future__ import annotations

import numpy as np

CONFIRMATORY = {
    "Q-T1S1": 1 / 5,
    "Q-T1S2": 1 / 20,
    "Q-T2S3": 1 / 5,
    "Q-T3S1": 1 / 5,
    "Q-T3S2": 1 / 5,
    "Q-T3S3": 1 / 10,
}
DESCRIPTIVE_CELLS = ("Q-T1S3", "Q-T2S1", "Q-T2S2")
TICK_BINS = ("Q-T1", "Q-T2", "Q-T3")  # every tick bin computed (simulation, seeds, reports)
FAMILY_TICK_BINS = ("Q-T1", "Q-T3")  # the tick bins of the one Holm family (2 Nov 2026)
DESCRIPTIVE_TICK_BINS = tuple(tb for tb in TICK_BINS if tb not in FAMILY_TICK_BINS)  # Q-T2


def tick_bin(cell: str) -> str:
    """'Q-T2S3' -> 'Q-T2'."""
    return cell[:4]


A_GRID = np.round(np.arange(0.0, 1.001, 0.25), 2)
T_GRID = 2 ** np.arange(7)  # 1 .. 64
RANGES = ((4, 32), (8, 64))
PRIMARY_RANGE, FALLBACK_RANGE = (8, 64), (4, 32)
A_MAX = 1.0
CLIP_QUANTILE = 0.999  # R2

PRIMARY_TESTS = ("P-V", "P-R", "P-M", "P-F")
DESCRIPTIVE_TESTS = ("P-K", "P-C")
ALPHA_FW = 0.05

LEVEL = 0.90  # two-sided interval level; its one-sided bounds are 5 % tests, and it is the TOST interval
N_BOOT = 200

PR_CENTRE, PR_MARGIN = 0.5, 0.15  # P-R: a*(T1) and a*(T2) within 1/2 +- 0.15

# P-M parameters chosen in Phase 3d-2 (allowed values fixed here; no default)
PM_BINS_ALLOWED = (3, 4, 5)
PM_LMAX_ALLOWED = (16, 32)
PM_L_FACTOR = 8  # L >= 8 l_max
