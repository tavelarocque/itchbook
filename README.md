# itchbook

A C++20 engine that rebuilds NASDAQ's full order-by-order (L3) book from real TotalView-ITCH 5.0 trading days.
It is then used as ground truth for two microstructure studies:

1. How wrong are standard L2 queue-position models, measured against the true L3 queue?
2. Do deep order-book models beat a simple order-flow-imbalance baseline after costs?

## Status

| Milestone | What | State |
|---|---|---|
| M0 | Build, CI, sanitizers, test/bench/fuzz tooling, data scripts | done |
| M1 | ITCH 5.0 parser: zero-copy views of all 23 message types, bounded-memory gzip streaming | done: 423,285,709 messages on 01302020, identical to the independent reference and the published counts |
| M2 | L3 book (`ReferenceBook`, `itch-book`) | done: 0 orphans on 01302020, empty book at close, whole-market checksums identical to the independent Python book at 35 times of day |
| M3 | Optimised engine (`FastBook`), identical output, benchmarked against CppTrader on the same file and machine | done: byte-identical to the M2 book on 01302020; 92 ns/message vs about 315 for CppTrader's L3 book on one pinned Xeon core |
| M4 | Research 1: how wrong are L2 queue-position models against the L3 truth? | done: cancellations come overwhelmingly from behind (15% of cancelled volume ahead; 4–5% at mid-queue vs 50% under the uniform model); same-price replaces move 6.6% of the queue ahead invisibly; behind 10–24 orders the uniform model's fills are 83% real, the pessimistic model's 99% real but it misses 56% of real fills |
| M5 | Research 2: forecasting audit | next |

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
tools/itch_stats.cpp     counts, fingerprint and throughput for a file
tools/itch_book.cpp      L3 book replay: counters, checksums, snapshots
tools/itch_bench.cpp     book throughput and per-update latency
tools/itch_queue.cpp     M4 queue-position study
tools/queue/             M4 analysis (pandas, matplotlib)
tools/itchref/           independent Python reference parser and cross-checks
tests/  bench/  fuzz/    GoogleTest + RapidCheck, Google Benchmark, libFuzzer
docs/                    walkthroughs (one PDF per milestone)
scripts/                 data download, benchmarks, public-repo sync, walkthrough renderer
results/                 raw benchmark output
```
