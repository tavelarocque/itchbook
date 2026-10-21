#!/usr/bin/env bash
# M5 overnight run: train DeepLOB (book input) and DeepOFI (order-flow input), each within a time
# budget, then score everything (ridge, both deep models, the oracle) on the validation and test
# days. Keeps the Mac awake for the duration. Log: $ITCH_DATA/m5/train.log
#
#   scripts/m5_train.sh [BUDGET_MIN_PER_MODEL]      default 150
set -uo pipefail
budget="${1:-150}"
repo="$(cd "$(dirname "$0")/.." && pwd)"
data="${ITCH_DATA:-$HOME/CODE/itch-data}"
m5="$data/m5"
export PATH="/opt/homebrew/bin:$HOME/.local/bin:$PATH"
cd "$repo"
run() { uv run -q --group analysis --group deep python tools/m5/run.py "$@" --study "$m5/study.json" --data "$m5/features" --preds "$m5/preds"; }
{
  echo "# m5_train $(date) itchbook $(git rev-parse --short HEAD) budget ${budget} min/model"
  sysctl -n vm.swapusage 2>/dev/null
  run deep --input lob --epochs 40 --steps 1000 --batch 256 --patience 3 --val-samples 50000 --budget-min "$budget" || echo "DEEPLOB_FAILED"
  run deep --input ofi --epochs 40 --steps 1000 --batch 256 --patience 3 --val-samples 50000 --budget-min "$budget" || echo "DEEPOFI_FAILED"
  uv run -q --group analysis python tools/m5/run.py evaluate --study "$m5/study.json" --data "$m5/features" \
    --preds "$m5/preds" --repo "$m5" || echo "EVALUATE_FAILED"
  echo "# finished $(date)"
} > "$m5/train.log" 2>&1
