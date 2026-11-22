#!/usr/bin/env bash
# Downloads the FASTQ files of one SRA/ENA study and checks their MD5.
#
#   bash download.sh <study accession> <data_dir> [jobs]
#
# Writes <data_dir>/<accession>.runs.tsv (ENA file report) and
# <data_dir>/<accession>/fastq/*.fastq.gz. Files whose MD5 already matches
# are skipped, so the script can be re-run after an interruption.
set -euo pipefail

ACC=${1:?usage: download.sh <accession> <data_dir> [jobs]}
DATA=${2:?usage: download.sh <accession> <data_dir> [jobs]}
JOBS=${3:-12}
FIELDS=run_accession,sample_accession,sample_alias,sample_title,library_name,experiment_title,instrument_model,library_layout,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes
mkdir -p "$DATA/$ACC/fastq"
REPORT="$DATA/$ACC.runs.tsv"

if [ ! -s "$REPORT" ]; then
  for attempt in 1 2 3 4 5; do
    if curl -sf -o "$REPORT.tmp" \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$ACC&result=read_run&fields=$FIELDS"; then
      mv "$REPORT.tmp" "$REPORT"
      break
    fi
    sleep $((attempt * 30))
  done
fi
[ -s "$REPORT" ] || { echo "no file report for $ACC" >&2; exit 1; }

col() { head -1 "$REPORT" | tr '\t' '\n' | grep -n "^$1\$" | cut -d: -f1; }
URL=$(col fastq_ftp)
MD5=$(col fastq_md5)

# One line per file: <md5> <url>
tail -n +2 "$REPORT" | awk -F'\t' -v u="$URL" -v m="$MD5" '{
  n = split($u, urls, ";"); split($m, sums, ";")
  for (i = 1; i <= n; i++) if (urls[i] != "") print sums[i], urls[i]
}' > "$DATA/$ACC/files.txt"

fetch() {
  local sum=$1 url=$2 out
  out="$DATA/$ACC/fastq/$(basename "$url")"
  if [ -s "$out" ] && [ "$(md5sum "$out" | cut -d' ' -f1)" = "$sum" ]; then
    return 0
  fi
  for attempt in 1 2 3; do
    curl -sf --retry 5 -o "$out" "https://$url" \
      && [ "$(md5sum "$out" | cut -d' ' -f1)" = "$sum" ] && return 0
    sleep 10
  done
  echo "FAILED $url" >&2
  return 1
}
export -f fetch
export DATA ACC
xargs -P "$JOBS" -L 1 bash -c 'fetch "$0" "$1"' < "$DATA/$ACC/files.txt"
echo "$(date +%H:%M) $ACC: $(ls "$DATA/$ACC/fastq" | wc -l) files," \
  "$(du -sh "$DATA/$ACC/fastq" | cut -f1)"
