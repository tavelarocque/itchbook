# itchbook

[![ci](https://github.com/tavelarocque/itchbook/actions/workflows/ci.yml/badge.svg)](https://github.com/tavelarocque/itchbook/actions/workflows/ci.yml)

A C++20 engine that rebuilds NASDAQ's full order-by-order (L3) book from real TotalView-ITCH 5.0 trading days.
It is then used as ground truth for two microstructure studies:

1. How wrong are standard L2 queue-position models, measured against the true L3 queue?
2. Do deep order-book models beat a simple order-flow-imbalance baseline after costs?

**At a glance**

- **Engine.** It replays a full NASDAQ day of 423 million messages at 92 ns per message on one core. That is 3.3× faster than CppTrader's equivalent order-by-order book, measured on the same file and machine. Its output is byte-identical to an independent Python implementation, checked at 35 times of day.
- **Queue study.** Across 7.5 million tracked orders, 85% of cancelled volume came from *behind* the order, not ahead of it. The pro-rata assumption behind common L2 backtest models is far too optimistic. Behind 10–24 orders, one in six of its simulated fills never happened.
- **Forecasting study.** Seven trading days, 20 symbols, strictly out-of-time test days. DeepLOB beats a ridge regression on multi-level order-flow imbalance (OFI) on macro-F1, 0.47–0.53 against 0.43–0.45. But both models' edge is 0.1–0.4 ticks per trade. Even on the 1-tick-spread stocks that make up most of the universe, a round trip costs 1.6 ticks: one tick of spread plus 0.6 ticks of taker fees. So neither model finds a single profitable trade after costs.

Each milestone has a walkthrough that explains the method and every design choice:

| Walkthrough | Covers |
|---|---|
| [M0 Setup](docs/m0_setup.pdf) | build, CI, sanitizers, test and fuzz tooling, data integrity |
| [M1 Parser](docs/m1_parser.pdf) | zero-copy ITCH 5.0 views, streaming gzip, verification against an independent parser |
| [M2 Book](docs/m2_book.pdf) | the L3 reference book and its whole-market cross-check |
| [M3 Performance](docs/m3_performance.pdf) | `FastBook`, each optimisation step, the CppTrader benchmark |
| [M4 Queue study](docs/m4_queue.pdf) | L2 queue-position models against the true L3 queue |
| [M5 Forecasting](docs/m5_forecasting.pdf) | DeepLOB and an order-flow model against an OFI baseline, after costs |

## Status

| Milestone | What | State |
|---|---|---|
| M0 | Build, CI, sanitizers, test/bench/fuzz tooling, data scripts | done |
| M1 | ITCH 5.0 parser: zero-copy views of all 23 message types, bounded-memory gzip streaming | done: 423,285,709 messages on 01302020, identical to the independent reference and the published counts |
| M2 | L3 book (`ReferenceBook`, `itch-book`) | done: 0 orphans on 01302020, empty book at close, whole-market checksums identical to the independent Python book at 35 times of day |
| M3 | Optimised engine (`FastBook`), identical output, benchmarked against CppTrader on the same file and machine | done: byte-identical to the M2 book on 01302020; 92 ns/message vs about 315 for CppTrader's L3 book on one pinned Xeon core |
| M4 | Research 1: how wrong are L2 queue-position models against the L3 truth? | done: cancellations come overwhelmingly from behind (15% of cancelled volume ahead; 4–5% at mid-queue vs 50% under the uniform model); same-price replaces move 6.6% of the queue ahead invisibly; behind 10–24 orders the uniform model's fills are 83% real, the pessimistic model's 99% real but it misses 56% of real fills |
| M5 | Research 2: deep LOB models vs an OFI baseline, after costs | done: DeepLOB beats ridge-on-OFI on macro-F1 (0.47–0.53 vs 0.43–0.45) mainly by predicting "no move"; ridge calls direction as well or better; a deep model on OFI does no better than DeepLOB; none of the three has a profitable trade after the spread and the $0.003 taker fee, against up to 7.6 ticks/trade available to perfect foresight. Both deep models were trained until validation performance stopped improving |

## Performance

Full day 01302020 (423,285,709 messages, 12.95 GB uncompressed). The time is total processing time divided by total messages, with input read on the measured thread. The AWS m7i.2xlarge runs (Xeon Platinum 8488C) were pinned to one core with interleaved rounds. Raw output is in `results/`.

| Engine | Book | Xeon 8488C, ns/msg | MacBook Air M2, ns/msg |
|---|---|---|---|
| CppTrader `market_manager` | L3, time priority | 309–321 | 224–254 |
| CppTrader `market_manager_optimized` | aggregated levels only; 16 GB order array | 75–82 | can't run (8 GB RAM) |
| itchbook `ReferenceBook` (M2 oracle) | L3, time priority | 336–356 (two runs) | |
| **itchbook `FastBook`** | L3, time priority | **91–97** | **89–97** |

Per-update latency on the Xeon (rdtsc, with the 12.5 ns timer cost included): p50 54 ns, p99 517 ns, p99.9 743 ns.
The order index peaks at 1.96 GB.
`docs/m3_performance.pdf` explains each design step, including why software prefetching ended up off by default.

```bash
scripts/bench_cpptrader.sh 01302020 3     # builds CppTrader at a pinned commit and runs both, interleaved
scripts/bench_prefetch_ab.sh 01302020 3   # what each prefetch stage is worth
```

## Research 1: L2 queue models against the L3 truth

`itch-queue` follows 1 in 16 orders that join the best bid or ask behind other orders during regular hours (4.6 million on 01302020).
At every trade print at their price, it compares each L2 model's estimate of the shares ahead with the true count from the L3 book.
The truth was checked against the book's own queue walk at all 1.4 million evaluation points over the two days, with no failures.
The models include hftbacktest's (`RiskAdverseQueueModel` and `ProbQueueModel` with its power, log and `PowerProbQueueFunc3` functions), plus one calibrated on the other trading day.

| Model (01302020) | Mean error / queue at entry | Fill precision, queue 10–24 | Fill recall, queue 10–24 |
|---|---|---|---|
| pessimistic | 0.076 | 98.9% | 43.7% |
| uniform (pro-rata) | 0.118 | 83.4% | 85.7% |
| hftbacktest power2 | 0.105 | 91.9% | 80.8% |
| hftbacktest front3 (r³) | 0.081 | 96.8% | 63.1% |
| calibrated on 12302019 | 0.075 | 97.5% | 57.1% |

```bash
scripts/queue_study.sh 16                                          # both days, both passes (~6 min on an M2 Air)
uv run --group analysis tools/queue/analyze.py --data ~/CODE/itch-data/m4 --days 12302019 01302020
```

`docs/m4_queue.pdf` has the method, the figures and the limits of the study.

## Research 2: deep order-book models against an order-flow baseline

`itch-lob` records a 250 ms grid of book features for 20 symbols over seven NASDAQ days (2019–20): the top 10 levels per side, multi-level OFI accumulated over every book event, and aggressor-signed volume.
The universe was fixed on the first day, before any other day was looked at.
Models are trained on four days, every choice (penalty, thresholds, early stopping) is made on a fifth, and they are scored once on the last two.

| Test days (12302019, 01302020) | Ridge on multi-level OFI | DeepLOB | Deep model on OFI |
|---|---|---|---|
| Macro-F1, 1 s / 5 s / 30 s | 0.43 / 0.45 / 0.43 | 0.47 / 0.53 / 0.51 | 0.43 / 0.53 / 0.49 |
| Direction right when it calls a move, 1 s | 62.7% | 55.7% | 53.6% |
| Edge per trade before costs (ticks), 1 s / 5 s / 30 s | 0.14 / 0.25 / 0.23 | 0.11 / 0.21 / 0.43 | 0.14 / 0.33 / 0.39 |
| Profitable taker trades after spread and fee | none | none | none |

Perfect foresight earns 1.1 / 2.6 / 7.6 ticks per trade after the same costs.

`docs/m5_forecasting.pdf` has the method, the results and the limits.

```bash
scripts/m5_prepare.sh          # rank day 1, fix the universe, extract features, fit the ridge baseline
scripts/m5_train.sh 150        # train DeepLOB and the order-flow model, then score everything
```

## Build

You need a C++20 compiler (GCC 13+, Clang 17+ or Apple clang 15+), CMake 3.24+, Ninja and zlib.
The test and benchmark dependencies (GoogleTest, RapidCheck and Google Benchmark) are fetched at pinned revisions.

```bash
cmake --preset release && cmake --build --preset release && ctest --preset release
cmake --preset asan    && cmake --build --preset asan    && ctest --preset asan    # ASan + UBSan
cmake --preset native  && cmake --build --preset native                          # -march/-mcpu=native for benchmarks
```

The fuzz preset needs a clang that ships libFuzzer. Apple clang does not, so use Linux or Homebrew LLVM:
`cmake --preset fuzz && cmake --build --preset fuzz && build/fuzz/fuzz/fuzz_stream -max_total_time=60`.

## Data

The sample days come from NASDAQ at `https://emi.nasdaq.com/ITCH/Nasdaq%20ITCH/`. They are 3.5–18 GB gzipped, and they never go in git.

```bash
scripts/fetch_data.sh 12302019 01302020     # into ~/CODE/itch-data (override with ITCH_DATA)
```

`data/manifest.tsv` lists each file's size from the server listing. The server's `.md5sum` files return 404, so integrity is checked two ways:

- `gzip -t`, which verifies every gzip member's CRC32;
- a SHA-256 recorded the first time a download verifies.

## Verifying the parser

The parser is verified against things it can't fake:

- **Unit tests.** Every field of all 23 message types is encoded by a sequential test writer and decoded by the offset-based views.
- **Tiling test.** Every byte of every message belongs to exactly one field, so the per-field fingerprint below misses nothing.
- **Property tests (RapidCheck).** Any message sequence, cut into chunks of any size, gzipped or raw, decodes to exactly the encoded messages. Any byte string is fully accounted for.
- **libFuzzer.** Arbitrary input runs through the splitter and every accessor under ASan and UBSan.
- **Independent Python reference** (`tools/itchref`). It is a separate implementation that uses `struct` format strings. It prints the same canonical counts and per-field sums as `itch-stats --fingerprint`, and `crosscheck.py` diffs the two.
- **Mutation tests of the cross-check.** Swapping any two adjacent fields in the Python layout must break the match. Without this, a broken check would pass silently.
- **Published counts.** On 01302020, per-type counts are compared with a published third-party table.

```bash
build/release/itch-book ~/CODE/itch-data/01302020.NASDAQ_ITCH50.gz      # L3 book: counters, checksums, top of book
build/release/itch-book --engine reference ~/CODE/itch-data/01302020.NASDAQ_ITCH50.gz   # same output from the M2 oracle
python3 tools/itchref/crosscheck_book.py --itch-book build/release/itch-book --file ~/CODE/itch-data/01302020.NASDAQ_ITCH50.gz
build/release/itch-stats --fingerprint ~/CODE/itch-data/01302020.NASDAQ_ITCH50.gz
python3 tools/itchref/crosscheck.py --itch-stats build/release/itch-stats \
    --file ~/CODE/itch-data/01302020.NASDAQ_ITCH50.gz --fingerprint --published 01302020
```

## Layout

```
include/itchbook/itch/   endian loads, message views, dispatch, field reflection, frame splitter
include/itchbook/io/     byte sources: file/stdin, streaming gzip, memory
src/                     source implementations
include/itchbook/book/   ReferenceBook (oracle), FastBook, order indexes, BookBuilder
include/itchbook/bench/  cycle-counter clock, latency histogram
include/itchbook/queue/  L2 queue models and the QueueStudy engine
include/itchbook/lob/    250 ms book and order-flow recorder for M5
tools/itch_stats.cpp     counts, fingerprint and throughput for a file
tools/itch_book.cpp      L3 book replay: counters, checksums, snapshots
tools/itch_bench.cpp     book throughput and per-update latency
tools/itch_queue.cpp     M4 queue-position study
tools/itch_lob.cpp       M5 feature extraction (book levels, multi-level OFI, signed volume)
tools/queue/             M4 analysis (pandas, matplotlib)
tools/m5/                M5 ridge baseline, DeepLOB and order-flow models (PyTorch), cost-aware evaluation
tools/itchref/           independent Python reference parser and cross-checks
tests/  bench/  fuzz/    GoogleTest + RapidCheck, Google Benchmark, libFuzzer
docs/                    walkthroughs (one PDF per milestone)
scripts/                 data download, benchmarks, public-repo sync, walkthrough renderer
results/                 raw benchmark and study results
```
