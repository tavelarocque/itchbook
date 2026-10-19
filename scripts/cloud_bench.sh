#!/usr/bin/env bash
# One-shot benchmark on a fresh Ubuntu 24.04 cloud machine (x86-64, >= 32 GB
# RAM). Copy the repo over as a git bundle, then run this script from the clone:
#
#   # on your Mac
#   scp ~/CODE/itchbook/target/itchbook.bundle ubuntu@HOST:
#   ssh ubuntu@HOST 'git clone -q itchbook.bundle itchbook && itchbook/scripts/cloud_bench.sh'
#   scp ubuntu@HOST:itchbook/results/m3_*.txt ~/CODE/itchbook/results/
#
# It installs the toolchain, downloads and verifies 01302020 from NASDAQ,
# builds itchbook and CppTrader, and runs scripts/bench_cpptrader.sh including
# CppTrader's market_manager_optimized. About an hour end to end, most of it
# the download.
set -euo pipefail
repo="$(cd "$(dirname "$0")/.." && pwd)"
export ITCH_DATA="$HOME/itch-data" CPPTRADER_DIR="$HOME/CppTrader" OPTIMIZED=1

sudo apt-get update -qq
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y -qq build-essential cmake ninja-build zlib1g-dev git \
  python3-pip uuid-dev libssl-dev binutils-dev util-linux > /dev/null

echo "== machine"
lscpu | grep -E "Model name|^CPU\(s\)|Thread|MHz" || true
free -g | head -2

"$repo/scripts/fetch_data.sh" 01302020
# Pick a core away from 0, where the kernel handles most interrupts.
CPU="${CPU:-$(( $(nproc) > 2 ? 2 : 0 ))}" "$repo/scripts/bench_cpptrader.sh" 01302020 3
