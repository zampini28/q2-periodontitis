#!/usr/bin/env bash
# One command, from a fresh clone on Linux x86-64 to every result of the
# validation (docs/ANALYSIS_PLAN.md):
#
#   bash validation/reproduce.sh [data_dir]      (default ~/q2-validation-data)
#
# 1. installs micromamba in ~/.local/bin if it is missing;
# 2. creates the QIIME 2 amplicon 2024.10 environment (Q2_ENV, default
#    q2-amplicon-2024.10), adds environment.yml and installs this plugin;
#    only the channels of each yml are used, not those of ~/.condarc;
# 3. creates the R environment with microbiome 1.6.0 (R_ENV, default
#    r-microbiome-1.6.0, from r-microbiome-1.6.0.yml);
# 4. runs the plugin and validation test suites, including the R
#    cross-check on simulated data, and stops if they fail;
# 5. runs run_all.sh: every cohort of cohorts.tsv, then the multi-cohort
#    analysis, robustness, performance, simulation, the R cross-check and
#    plugin-vs-R benchmark, figures and REPORT.md in validation/results.
#
# Every step is skipped when already done, so the command can be re-run after
# an interruption. Needs curl, bzip2 and about 60 GB free in data_dir; takes
# about a day on 4 cores. Run it inside tmux or with nohup. THREADS defaults
# to the number of cores.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(dirname "$HERE")
DATA=${1:-$HOME/q2-validation-data}
export Q2_ENV=${Q2_ENV:-q2-amplicon-2024.10}
export R_ENV=${R_ENV:-r-microbiome-1.6.0}
export THREADS=${THREADS:-$(nproc)}
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba}
export PATH=$HOME/.local/bin:$PATH
MICROMAMBA_VERSION=2.0.5
Q2_YML=https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

log() { echo "$(date '+%F %H:%M') $*"; }
mkdir -p "$DATA"

free_gb=$(df -BG --output=avail "$DATA" | tail -1 | tr -dc 0-9)
[ "$free_gb" -ge 60 ] || log "WARNING: only ${free_gb} GB free in $DATA (about 60 GB needed)"

if ! command -v micromamba > /dev/null; then
  log "installing micromamba $MICROMAMBA_VERSION"
  mkdir -p "$HOME/.local"
  curl -sfL "https://micro.mamba.pm/api/micromamba/linux-64/$MICROMAMBA_VERSION" \
    | tar -xj -C "$HOME/.local" bin/micromamba
fi

if ! micromamba env list | awk '{print $1}' | grep -qx "$Q2_ENV"; then
  log "creating $Q2_ENV"
  curl -sfL -o "$DATA/qiime2-amplicon-2024.10.yml" "$Q2_YML"
  micromamba create -y --override-channels -n "$Q2_ENV" -f "$DATA/qiime2-amplicon-2024.10.yml"
fi
# Outside the check above, so that an interrupted creation is completed
micromamba install -y --override-channels -n "$Q2_ENV" -f "$REPO/environment.yml"

if ! micromamba env list | awk '{print $1}' | grep -qx "$R_ENV"; then
  log "creating $R_ENV"
  micromamba create -y --override-channels -n "$R_ENV" -f "$HERE/r-microbiome-1.6.0.yml"
fi
RSCRIPT=$(micromamba run -n "$R_ENV" which Rscript)
export RSCRIPT

run() { micromamba run -n "$Q2_ENV" "$@"; }
log "installing q2-periodontitis"
run pip install -q -e "$REPO"
run qiime dev refresh-cache

log "tests"
(cd "$REPO" && run pytest -q q2_periodontitis/tests validation/analysis/tests)

log "validation (log per cohort in $DATA/<cohort>.log)"
bash "$HERE/run_all.sh" "$DATA"
log "done: $HERE/results/REPORT.md"
