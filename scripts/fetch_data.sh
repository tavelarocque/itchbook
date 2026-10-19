#!/usr/bin/env bash
# Download NASDAQ TotalView-ITCH 5.0 sample days from emi.nasdaq.com and
# verify them. Data never goes in git; it lives in $ITCH_DATA (default
# ~/CODE/itch-data).
#
#   scripts/fetch_data.sh 12302019 01302020
#
# Verification, in order:
#   1. byte size equals the size in the server listing (data/manifest.tsv);
#   2. `gzip -t` passes (checks every member's CRC32 and length);
#   3. SHA-256 equals the hash recorded in data/manifest.tsv, once one has been
#      recorded. The server's own .md5sum files return 404, so the first
#      verified download defines the hash and later downloads must match it.
set -euo pipefail

repo="$(cd "$(dirname "$0")/.." && pwd)"
manifest="$repo/data/manifest.tsv"
dest="${ITCH_DATA:-$HOME/CODE/itch-data}"
base='https://emi.nasdaq.com/ITCH/Nasdaq%20ITCH'
mkdir -p "$dest"

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
fsize() { if stat -f%z "$1" >/dev/null 2>&1; then stat -f%z "$1"; else stat -c%s "$1"; fi; }

[ $# -gt 0 ] || { echo "usage: $0 DAY [DAY...]   (DAY as in data/manifest.tsv, e.g. 01302020)" >&2; exit 1; }

for day in "$@"; do
  row="$(awk -F'\t' -v d="$day" '$1 == d' "$manifest")"
  [ -n "$row" ] || { echo "$day: not in $manifest" >&2; exit 1; }
  file="$(cut -f2 <<<"$row")"; size="$(cut -f3 <<<"$row")"; hash="$(cut -f4 <<<"$row")"
  out="$dest/$file"

  echo "== $file ($size bytes) -> $out"
  if [ ! -f "$out" ] || [ "$(fsize "$out")" != "$size" ]; then
    curl -fL --retry 5 --retry-delay 10 -C - -o "$out" "$base/$file"
  fi

  got="$(fsize "$out")"
  [ "$got" = "$size" ] || { echo "$file: size $got, expected $size (re-run to resume)" >&2; exit 1; }
  echo "   size ok; gzip -t (takes a minute or two)..."
  gzip -t "$out"
  echo "   gzip CRC ok; hashing..."
  h="$(sha256 "$out")"
  if [ "$hash" = "-" ]; then
    echo "   sha256 $h (no recorded hash yet: add it to data/manifest.tsv and commit)"
  elif [ "$h" = "$hash" ]; then
    echo "   sha256 matches manifest"
  else
    echo "$file: sha256 $h does not match manifest $hash" >&2; exit 1
  fi
done
