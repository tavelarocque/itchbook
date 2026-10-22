#!/usr/bin/env bash
# M5 data preparation, safe to rerun as days arrive.
#
#   1. On the first training day, rank symbols by trade prints and fix the
#      universe (study.json) before any other day is looked at.
#   2. Extract 250 ms book features for every study day that is downloaded
#      and not yet extracted (about a minute per day).
#   3. When every day is extracted, fit the ridge OFI baseline.
#
#   scripts/m5_prepare.sh        outputs in $ITCH_DATA/m5 (default ~/CODE/itch-data/m5)
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
data="${ITCH_DATA:-$HOME/CODE/itch-data}"
out="$data/m5"
mkdir -p "$out/features" "$out/preds"
export PATH="/opt/homebrew/bin:$HOME/.local/bin:$PATH"
(cd "$repo" && cmake --preset release -DITCHBOOK_BUILD_TESTS=OFF -DITCHBOOK_BUILD_BENCH=OFF > /dev/null &&
  cmake --build --preset release --target itch-lob > /dev/null)
bin="$repo/build/release/itch-lob"
first=01302019
if [ ! -f "$out/study.json" ]; then
  [ -f "$data/$first.NASDAQ_ITCH50.gz" ] || { echo "waiting for $first" >&2; exit 0; }
  "$bin" --rank --top 200 "$data/$first.NASDAQ_ITCH50.gz" > "$out/rank_$first.tsv"
  (cd "$repo" && uv run -q --group analysis python tools/m5/run.py universe --rank "$out/rank_$first.tsv" --out "$out/study.json")
fi
syms="$(python3 -c "import json,sys; print(','.join(json.load(open(sys.argv[1]))['universe']))" "$out/study.json")"
days="$(python3 -c "import json,sys; j=json.load(open(sys.argv[1])); print(' '.join(j['train']+j['val']+j['test']))" "$out/study.json")"
pending=0
for d in $days; do
  gz="$data/$d.NASDAQ_ITCH50.gz"
  if [ -f "$out/features/$d.meta.json" ]; then
    continue
  fi
  # Only a fully downloaded day: fetch_data.sh writes the whole file before verifying it.
  if [ ! -f "$gz" ] || [ "$(stat -f%z "$gz" 2>/dev/null || stat -c%s "$gz")" != "$(awk -F'\t' -v d="$d" '$1==d {print $3}' "$repo/data/manifest.tsv")" ]; then
    echo "$d: not downloaded yet"
    pending=1
    continue
  fi
  echo "== extracting $d"
  "$bin" --symbols "$syms" --out "$out/features/$d" "$gz"
done
if [ "$pending" = 0 ] && [ ! -f "$out/preds/ridge.meta.json" ]; then
  echo "== ridge baseline"
  (cd "$repo" && uv run -q --group analysis python tools/m5/run.py ridge --study "$out/study.json" \
    --data "$out/features" --preds "$out/preds")
fi
echo DONE
