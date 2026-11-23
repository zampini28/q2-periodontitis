#!/usr/bin/env bash
# Runs the analysis of the processed cohorts.
#
#   bash analysis/analyse.sh <data_dir> <results_dir>
set -euo pipefail
DATA=$(readlink -f "$1")
RESULTS=$(readlink -f "$2")
HERE=$(cd "$(dirname "$0")" && pwd)
ENV=${Q2_ENV:-q2-amplicon-2024.10}
R_ENV=${R_ENV:-r-microbiome-1.6.0}
run() { micromamba run -n "$ENV" python -I "$@"; }
run "$HERE/collect.py" "$DATA" "$RESULTS"
run "$HERE/run_analysis.py" "$RESULTS"
run "$HERE/robustness.py" "$DATA" "$RESULTS"
run "$HERE/performance.py" "$RESULTS"
run "$HERE/simulation.py" "$RESULTS"
# A failed R cross-check still lets the figures and report be written
RSCRIPT=$(micromamba run -n "$R_ENV" which Rscript)
rcheck=0
run "$HERE/rcheck.py" "$DATA" "$RESULTS" --rscript "$RSCRIPT" || rcheck=$?
run "$HERE/parser_check.py" "$DATA/references" "$RESULTS" \
  --rscript "$RSCRIPT" || rcheck=$?
run "$HERE/figures.py" "$RESULTS"
run "$HERE/report.py" "$RESULTS"
exit "$rcheck"
