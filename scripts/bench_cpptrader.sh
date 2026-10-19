#!/usr/bin/env bash
# Head-to-head benchmark: itchbook vs CppTrader on the same NASDAQ day, same
# machine, same compiler flags, single-threaded, uncompressed input read on the
# measured thread (CppTrader's own setup). Linux or WSL2 with >= 24 GB RAM:
# CppTrader's optimised book pre-allocates a flat array indexed by order
# reference number.
#
#   scripts/bench_cpptrader.sh [DAY] [REPEATS]       e.g. scripts/bench_cpptrader.sh 01302020 3
#
# Environment: ITCH_DATA (default ~/itch-data), CPPTRADER_DIR (default
# ~/CppTrader), CPU (core to pin to, default 2).
set -euo pipefail

day="${1:-01302020}"
repeats="${2:-3}"
repo="$(cd "$(dirname "$0")/.." && pwd)"
data="${ITCH_DATA:-$HOME/itch-data}"
ct="${CPPTRADER_DIR:-$HOME/CppTrader}"
cpu="${CPU:-2}"
ct_rev=39421f50ae673be34235e24eac95827220763652  # CppTrader commit the published numbers were checked against
out="$repo/results/m3_$(hostname)_$(date +%Y%m%d_%H%M).txt"
mkdir -p "$repo/results"

gz="$data/$day.NASDAQ_ITCH50.gz"
raw="$data/$day.NASDAQ_ITCH50"
[ -f "$gz" ] || { echo "missing $gz (scripts/fetch_data.sh $day, with ITCH_DATA=$data)" >&2; exit 1; }
if [ ! -f "$raw" ]; then
  echo "decompressing $gz once (about 13 GB for 01302020)..."
  gzip -dc "$gz" > "$raw.tmp" && mv "$raw.tmp" "$raw"
fi

# Largest order reference in the file: CppTrader's optimised book indexes an
# array by it without bounds checks, sized for 300,000,000. 01302020 reaches
# 494,466,152, so the array is enlarged (one constant, nothing else changed).
cmake --preset native -S "$repo" > /dev/null
cmake --build --preset native --target itch-stats itch-bench > /dev/null
max_ref="$("$repo/build/native/itch-stats" "$raw" 2>&1 >/dev/null | sed -n 's/.*max order reference \([0-9]*\).*/\1/p')"
echo "max order reference: $max_ref"

if [ ! -d "$ct/.git" ]; then
  git clone https://github.com/chronoxor/CppTrader.git "$ct"
fi
git -C "$ct" fetch -q origin && git -C "$ct" checkout -q "$ct_rev"
export PATH="$HOME/.local/bin:$PATH"
# Ubuntu 24.04 marks the system Python as externally managed (PEP 668); gil is a
# small pure-Python tool, so a user-level install is safe.
command -v gil >/dev/null || python3 -m pip install --user --break-system-packages gil
(cd "$ct" && gil update > /dev/null)
size=$(( (max_ref / 100000000 + 1) * 100000000 ))
sed -i "s/_orders.resize([0-9]*);/_orders.resize($size);/" "$ct/performance/market_manager_optimized.cpp"
grep -n "_orders.resize" "$ct/performance/market_manager_optimized.cpp"
cmake -S "$ct" -B "$ct/temp" -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="-march=native" > /dev/null
ninja -C "$ct/temp" cpptrader-performance-itch_handler cpptrader-performance-market_manager \
  cpptrader-performance-market_manager_optimized > /dev/null

{
  echo "# itchbook vs CppTrader, $(date -Iseconds)"
  echo "# host: $(hostname); $(lscpu | sed -n 's/^Model name: *//p'); $(nproc) logical CPUs; $(free -g | awk '/Mem:/{print $2}') GB RAM"
  echo "# kernel: $(uname -r); compiler: $(c++ --version | head -1)"
  echo "# itchbook: $(git -C "$repo" rev-parse --short HEAD); CppTrader: $ct_rev (orders array resized to $size)"
  for m in "$ct/modules/"*/; do echo "#   module $(basename "$m"): $(git -C "$m" rev-parse --short HEAD 2>/dev/null || echo '?')"; done
  echo "# file: $raw ($(stat -c %s "$raw") bytes), pinned to CPU $cpu, $repeats repeats each"
} | tee "$out"

cat "$raw" > /dev/null  # warm the page cache so no run pays for disk reads
run() {
  local label="$1"; shift
  for i in $(seq 1 "$repeats"); do
    echo "## $label run $i" | tee -a "$out"
    taskset -c "$cpu" "$@" 2>&1 | tee -a "$out" | grep -Ei "processing time|latency|throughput|ns_per_message|messages_per_s|update_p|busy_s|orphans" || true
  done
}
run "cpptrader itch_handler (parse only)" "$ct/temp/cpptrader-performance-itch_handler" -i "$raw"
run "cpptrader market_manager (L3, time priority)" "$ct/temp/cpptrader-performance-market_manager" -i "$raw"
run "cpptrader market_manager_optimized (aggregated levels)" "$ct/temp/cpptrader-performance-market_manager_optimized" -i "$raw"
run "itchbook fast (L3, time priority)" "$repo/build/native/itch-bench" --inline-input --engine fast "$raw"
run "itchbook fast, no prefetch" "$repo/build/native/itch-bench" --inline-input --no-prefetch --engine fast "$raw"
run "itchbook fast, per-update latency" "$repo/build/native/itch-bench" --inline-input --latency --engine fast "$raw"
run "itchbook reference" "$repo/build/native/itch-bench" --inline-input --engine reference "$raw"
echo "results: $out"
