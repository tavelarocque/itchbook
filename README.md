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
| M2 | L3 book (`ReferenceBook`, `itch-book`) with zero-orphan, invariant and cross-implementation gates | in progress |
| M3 | Performance against CppTrader on the same file and machine | |
| M4 | Research 1: queue-position ground truth | |
| M5 | Research 2: forecasting audit | |

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
build/release/itch-stats --fingerprint ~/CODE/itch-data/01302020.NASDAQ_ITCH50.gz
python3 tools/itchref/crosscheck.py --itch-stats build/release/itch-stats \
    --file ~/CODE/itch-data/01302020.NASDAQ_ITCH50.gz --fingerprint --published 01302020
```

## Layout

```
include/itchbook/itch/   endian loads, message views, dispatch, field reflection, frame splitter
include/itchbook/io/     byte sources: file/stdin, streaming gzip, memory
src/                     source implementations
tools/itch_stats.cpp     counts, fingerprint and throughput for a file
tools/itchref/           independent Python reference parser and cross-checks
tests/  bench/  fuzz/    GoogleTest + RapidCheck, Google Benchmark, libFuzzer
docs/                    walkthroughs (one PDF per milestone)
scripts/                 data download, public-repo sync, walkthrough renderer
```
