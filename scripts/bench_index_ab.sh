#!/usr/bin/env bash
# A/B the FastBook order index (paged direct table vs open-addressing hash)
# against CppTrader's optimised book, interleaved, on one pinned core.
# Assumes scripts/bench_cpptrader.sh has already built CppTrader and
# decompressed the day.
#
#   scripts/bench_index_ab.sh [DAY] [ROUNDS]
set -euo pipefail
day="${1:-01302020}"
rounds="${2:-3}"
repo="$(cd "$(dirname "$0")/.." && pwd)"
data="${ITCH_DATA:-$HOME/CODE/itch-data}"
ct="${CPPTRADER_DIR:-$HOME/CODE/CppTrader}"
raw="$data/$day.NASDAQ_ITCH50"
pin=()
command -v taskset > /dev/null && pin=(taskset -c "${CPU:-2}")

cd "$repo"
cmake --preset release -DITCHBOOK_BUILD_TESTS=OFF -DITCHBOOK_BUILD_BENCH=OFF > /dev/null
cmake --build --preset release --target itch-bench > /dev/null
cmake -S . -B build/hash -G Ninja -DCMAKE_BUILD_TYPE=Release -DITCHBOOK_HASH_INDEX=ON \
  -DITCHBOOK_BUILD_TESTS=OFF -DITCHBOOK_BUILD_BENCH=OFF > /dev/null
cmake --build build/hash --target itch-bench > /dev/null

out="$repo/results/m3_index_ab_$(hostname -s)_$(date +%Y%m%d_%H%M).txt"
mkdir -p "$repo/results"
echo "# order index A/B, $(date '+%Y-%m-%dT%H:%M:%S%z'), itchbook $(git rev-parse --short HEAD), ${pin[*]:-no pinning}" | tee "$out"
run() {
  local label="$1"; shift
  echo "## $label" | tee -a "$out"
  ${pin[@]+"${pin[@]}"} "$@" 2>&1 | tee -a "$out" | grep -E "Processing time|message latency|ns_per_message|order_index" || true
}
for i in $(seq 1 "$rounds"); do
  echo "### round $i" | tee -a "$out"
  if [ -x "$ct/temp/cpptrader-performance-market_manager_optimized" ]; then
    run "cpptrader market_manager_optimized" "$ct/temp/cpptrader-performance-market_manager_optimized" -i "$raw"
  fi
  run "itchbook fast, paged index" build/release/itch-bench --inline-input "$raw"
  run "itchbook fast, hash index" build/hash/itch-bench --inline-input "$raw"
done
run "itchbook fast, paged index, per-update latency" build/release/itch-bench --inline-input --latency "$raw"
echo "results: $out"
