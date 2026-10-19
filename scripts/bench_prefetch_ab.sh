#!/usr/bin/env bash
# What is software prefetching worth? Interleaved runs of FastBook with both
# prefetch stages, each stage alone, none (lookahead cursors still running),
# and no lookahead at all, on one pinned core with input read inline.
#
#   scripts/bench_prefetch_ab.sh [DAY] [ROUNDS]
set -euo pipefail
day="${1:-01302020}"
rounds="${2:-3}"
repo="$(cd "$(dirname "$0")/.." && pwd)"
data="${ITCH_DATA:-$HOME/CODE/itch-data}"
raw="$data/$day.NASDAQ_ITCH50"
pin=()
command -v taskset > /dev/null && pin=(taskset -c "${CPU:-2}")

cd "$repo"
cmake --preset release -DITCHBOOK_BUILD_TESTS=OFF -DITCHBOOK_BUILD_BENCH=OFF > /dev/null
cmake --build --preset release --target itch-bench > /dev/null

out="$repo/results/m3_prefetch_ab_$(hostname -s)_$(date +%Y%m%d_%H%M).txt"
mkdir -p "$repo/results"
echo "# prefetch A/B, $(date '+%Y-%m-%dT%H:%M:%S%z'), itchbook $(git rev-parse --short HEAD), ${pin[*]:-no pinning}" | tee "$out"
run() {
  local label="$1"; shift
  echo "## $label" | tee -a "$out"
  ${pin[@]+"${pin[@]}"} build/release/itch-bench --inline-input "$@" "$raw" 2>&1 | tee -a "$out" | grep -E "ns_per_message" || true
}
for i in $(seq 1 "$rounds"); do
  echo "### round $i" | tee -a "$out"
  run "both stages" --prefetch-stages 3
  run "stage 0 only (index entry, book header)" --prefetch-stages 1
  run "stage 1 only (slot or best level)" --prefetch-stages 2
  run "lookahead, no prefetch" --prefetch-stages 0
  run "no lookahead (default)"
done
echo "results: $out"
