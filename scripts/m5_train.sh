#!/usr/bin/env bash
# M5 overnight run: train DeepLOB (book input) and DeepOFI (order-flow input), each within a time
# budget, then score everything (ridge, both deep models, the oracle) on the validation and test
# days. Keeps the Mac awake for the duration. Log: $ITCH_DATA/m5/train.log
#
#   scripts/m5_train.sh [BUDGET_MIN_PER_MODEL] [MODELS]
#
# BUDGET defaults to 150; MODELS to "lob ofi". A model whose checkpoint already exists in
# preds/ is not retrained, only used for prediction (move preds/deeplob.pt away to retrain).
# PATIENCE (default 3): epochs without validation improvement before stopping.
set -uo pipefail
budget="${1:-150}"
patience="${PATIENCE:-3}"
models="${2:-lob ofi}"
repo="$(cd "$(dirname "$0")/.." && pwd)"
data="${ITCH_DATA:-$HOME/CODE/itch-data}"
m5="$data/m5"
export PATH="/opt/homebrew/bin:$HOME/.local/bin:$PATH"
cd "$repo"
run() { uv run -q --group analysis --group deep python tools/m5/run.py "$@" --study "$m5/study.json" --data "$m5/features" --preds "$m5/preds"; }
{
  echo "# m5_train $(date) itchbook $(git rev-parse --short HEAD) budget ${budget} min/model"
  sysctl -n vm.swapusage 2>/dev/null
  for m in $models; do
    name=$([ "$m" = lob ] && echo deeplob || echo deepofi)
    extra=()
    [ -f "$m5/preds/$name.pt" ] && extra=(--predict-only)
    run deep --input "$m" --epochs 40 --steps 1000 --batch 256 --patience "$patience" --val-samples 50000 \
      --budget-min "$budget" ${extra[@]+"${extra[@]}"} || echo "${name}_FAILED"
  done
  uv run -q --group analysis python tools/m5/run.py evaluate --study "$m5/study.json" --data "$m5/features" \
    --preds "$m5/preds" --repo "$m5" || echo "EVALUATE_FAILED"
  echo "# finished $(date)"
} >> "$m5/train.log" 2>&1
