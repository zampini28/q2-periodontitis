#!/usr/bin/env bash
# Runs the analysis of the processed cohorts.
#
#   bash analysis/analyse.sh <data_dir> <results_dir>
set -euo pipefail
DATA=$(readlink -f "$1")
RESULTS=$(readlink -f "$2")
HERE=$(cd "$(dirname "$0")" && pwd)
ENV=${Q2_ENV:-q2-amplicon-2024.10}
run() { micromamba run -n "$ENV" python -I "$@"; }
run "$HERE/collect.py" "$DATA" "$RESULTS"
run "$HERE/run_analysis.py" "$RESULTS"
run "$HERE/robustness.py" "$DATA" "$RESULTS"
run "$HERE/performance.py" "$RESULTS"
run "$HERE/simulation.py" "$RESULTS"
run "$HERE/figures.py" "$RESULTS"
run "$HERE/report.py" "$RESULTS"
