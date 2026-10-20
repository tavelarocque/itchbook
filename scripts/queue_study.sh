#!/usr/bin/env bash
# M4: the queue-position study on both sample days.
#
# Pass 1 runs each day with the fixed models and measures where cancellations
# come from (PREFIX.cancels.tsv). Pass 2 reruns each day with the calibrated
# model fitted on the OTHER day, so the calibrated model is always scored out
# of sample. Pass 1 checks the truth against the L3 book's queue walk at every
# evaluation point (--verify 1). Results go to $ITCH_DATA/m4/ (default ~/CODE/itch-data/m4).
#
#   scripts/queue_study.sh [SAMPLE]      track 1 in SAMPLE eligible orders (default 16)
set -euo pipefail
sample="${1:-16}"
repo="$(cd "$(dirname "$0")/.." && pwd)"
data="${ITCH_DATA:-$HOME/CODE/itch-data}"
out="$data/m4"
mkdir -p "$out"
(cd "$repo" && cmake --preset release -DITCHBOOK_BUILD_TESTS=OFF -DITCHBOOK_BUILD_BENCH=OFF > /dev/null &&
  cmake --build --preset release --target itch-queue > /dev/null)
bin="$repo/build/release/itch-queue"
days=(12302019 01302020)
log="$out/run_$(date +%Y%m%d_%H%M).log"
echo "# itch-queue, itchbook $(git -C "$repo" rev-parse --short HEAD), sample 1/$sample, $(date)" | tee "$log"
for d in "${days[@]}"; do
  echo "## pass 1: $d" | tee -a "$log"
  "$bin" --sample "$sample" --events 16 --verify 1 --out "$out/$d" "$data/$d.NASDAQ_ITCH50.gz" | tee -a "$log"
done
for d in "${days[@]}"; do
  other="${days[0]}"; [ "$d" = "${days[0]}" ] && other="${days[1]}"
  echo "## pass 2: $d calibrated on $other" | tee -a "$log"
  "$bin" --sample "$sample" --verify 1000 --calibration "$out/$other.cancels.tsv" \
    --out "$out/$d.cal_$other" "$data/$d.NASDAQ_ITCH50.gz" | tee -a "$log"
done
echo "DONE" | tee -a "$log"
