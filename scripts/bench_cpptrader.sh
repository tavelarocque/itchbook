#!/usr/bin/env bash
# Head-to-head benchmark: itchbook vs CppTrader on the same NASDAQ day, same
# machine, single-threaded, uncompressed input read on the measured thread
# (CppTrader's own setup). Runs on Linux (pinned to one core with taskset) and
# macOS (no pinning available).
#
#   scripts/bench_cpptrader.sh [DAY] [REPEATS]       e.g. scripts/bench_cpptrader.sh 01302020 3
#
# Environment:
#   ITCH_DATA       data directory (default ~/CODE/itch-data)
#   CPPTRADER_DIR   CppTrader checkout (default ~/CODE/CppTrader)
#   CPU             core to pin to on Linux (default 2)
#   OPTIMIZED       1 to include CppTrader's market_manager_optimized. It
#                   pre-allocates a flat array indexed by order reference
#                   number, about 16 GB for 01302020, so it needs >= 24 GB RAM.
#                   Default: 1 if the machine has >= 24 GB, else 0.
#
# Rounds are interleaved (every engine once per round) so thermal throttling on
# a fanless laptop hits all engines alike.
set -euo pipefail

day="${1:-01302020}"
repeats="${2:-3}"
repo="$(cd "$(dirname "$0")/.." && pwd)"
data="${ITCH_DATA:-$HOME/CODE/itch-data}"
ct="${CPPTRADER_DIR:-$HOME/CODE/CppTrader}"
cpu="${CPU:-2}"
ct_rev=39421f50ae673be34235e24eac95827220763652  # CppTrader commit whose README holds the published numbers

os="$(uname -s)"
if [ "$os" = Darwin ]; then
  ram_gb=$(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))
  cpu_name="$(sysctl -n machdep.cpu.brand_string)"
  ncpu="$(sysctl -n hw.ncpu)"
  fsize() { stat -f%z "$1"; }
  pin=()
else
  ram_gb=$(( $(awk '/MemTotal/{print $2}' /proc/meminfo) / 1024 / 1024 ))
  cpu_name="$(lscpu | sed -n 's/^Model name: *//p')"
  ncpu="$(nproc)"
  fsize() { stat -c%s "$1"; }
  pin=(taskset -c "$cpu")
fi
optimized="${OPTIMIZED:-$([ "$ram_gb" -ge 24 ] && echo 1 || echo 0)}"
export PATH="$HOME/.local/bin:$PATH"
[ "$os" = Darwin ] && export PATH="/opt/homebrew/bin:$(python3 -c 'import site; print(site.USER_BASE)' 2>/dev/null)/bin:$PATH"

out="$repo/results/m3_$(hostname -s)_$(date +%Y%m%d_%H%M).txt"
mkdir -p "$repo/results"

gz="$data/$day.NASDAQ_ITCH50.gz"
raw="$data/$day.NASDAQ_ITCH50"
[ -f "$gz" ] || { echo "missing $gz (run scripts/fetch_data.sh $day)" >&2; exit 1; }
if [ ! -f "$raw" ]; then
  echo "decompressing $gz once (about 13 GB for 01302020)..."
  gzip -dc "$gz" > "$raw.tmp" && mv "$raw.tmp" "$raw"
fi

echo "building itchbook (release preset)..."
# Presets are looked up in the current directory, so build from the repo root.
(cd "$repo" && cmake --preset release > /dev/null && cmake --build --preset release --target itch-stats itch-bench > /dev/null)
bin="$repo/build/release"
max_ref="$("$bin/itch-stats" "$raw" 2>&1 >/dev/null | sed -n 's/.*max order reference \([0-9]*\).*/\1/p')"
echo "max order reference in $day: $max_ref"

echo "building CppTrader $ct_rev..."
[ -d "$ct/.git" ] || git clone -q https://github.com/chronoxor/CppTrader.git "$ct"
git -C "$ct" fetch -q origin
git -C "$ct" checkout -q -f "$ct_rev"
if ! command -v gil >/dev/null; then
  # Homebrew and Ubuntu 24.04 Pythons are "externally managed" (PEP 668) and
  # need --break-system-packages; older pips reject that flag. gil is a small
  # pure-Python tool, so a user-level install is safe either way.
  python3 -m pip install -q --user gil 2>/dev/null || python3 -m pip install -q --user --break-system-packages gil
  export PATH="$(python3 -c 'import site; print(site.USER_BASE)')/bin:$PATH"
fi
(cd "$ct" && gil update > /dev/null)
targets=(cpptrader-performance-itch_handler cpptrader-performance-market_manager)
size=0
if [ "$optimized" = 1 ]; then
  # Its order array is indexed by reference number with no bounds check and is
  # sized for 300,000,000; 01302020 reaches 494,466,152. Enlarge that one
  # constant and change nothing else.
  size=$(( (max_ref / 100000000 + 1) * 100000000 ))
  perl -pi -e "s/_orders\\.resize\\(\\d+\\);/_orders.resize($size);/" "$ct/performance/market_manager_optimized.cpp"
  grep -n "_orders.resize" "$ct/performance/market_manager_optimized.cpp"
  targets+=(cpptrader-performance-market_manager_optimized)
fi
cmake -S "$ct" -B "$ct/temp" -G Ninja -DCMAKE_BUILD_TYPE=Release > /dev/null
ninja -C "$ct/temp" "${targets[@]}" > /dev/null

{
  echo "# itchbook vs CppTrader, $(date '+%Y-%m-%dT%H:%M:%S%z')"
  echo "# host: $(hostname -s); $cpu_name; $ncpu logical CPUs; ${ram_gb} GB RAM; $os $(uname -r)"
  echo "# compiler: $(c++ --version | head -1)"
  echo "# itchbook $(git -C "$repo" rev-parse --short HEAD), release preset; CppTrader $ct_rev, Release"
  [ "$optimized" = 1 ] && echo "#   CppTrader optimized: orders array resized 300000000 -> $size (max order reference $max_ref)"
  [ "$optimized" = 1 ] || echo "#   CppTrader optimized: skipped (${ram_gb} GB RAM < 24 GB needed for its order array)"
  for m in "$ct/modules/"*/; do echo "#   CppTrader module $(basename "$m"): $(git -C "$m" rev-parse --short HEAD 2>/dev/null || echo '?')"; done
  echo "# file: $raw ($(fsize "$raw") bytes); ${pin[*]:-no CPU pinning}; $repeats interleaved rounds"
} | tee "$out"

run() {
  local label="$1"; shift
  echo "## $label" | tee -a "$out"
  ${pin[@]+"${pin[@]}"} "$@" 2>&1 | tee -a "$out" | grep -Ei "processing time|message latency|ns_per_message|busy_s|orphans" || true
}
for i in $(seq 1 "$repeats"); do
  echo "### round $i" | tee -a "$out"
  run "cpptrader itch_handler (parse only)" "$ct/temp/cpptrader-performance-itch_handler" -i "$raw"
  run "cpptrader market_manager (L3, time priority)" "$ct/temp/cpptrader-performance-market_manager" -i "$raw"
  if [ "$optimized" = 1 ]; then
    run "cpptrader market_manager_optimized (aggregated levels)" \
      "$ct/temp/cpptrader-performance-market_manager_optimized" -i "$raw"
  fi
  run "itchbook fast (L3, time priority)" "$bin/itch-bench" --inline-input --engine fast "$raw"
  run "itchbook fast, with software prefetch" "$bin/itch-bench" --inline-input --prefetch --engine fast "$raw"
done
run "itchbook fast, per-update latency" "$bin/itch-bench" --inline-input --latency --engine fast "$raw"
run "itchbook reference (M2)" "$bin/itch-bench" --inline-input --engine reference "$raw"
echo "results: $out"
