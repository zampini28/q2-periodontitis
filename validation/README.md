# Validation workflow

This directory contains the code for the pre-specified multi-cohort
validation described in [docs/ANALYSIS_PLAN.md](../docs/ANALYSIS_PLAN.md).
The single-cohort SRP075100 results currently reported in
[docs/VALIDATION.md](../docs/VALIDATION.md) are the legacy outputs in
[results/](results/); the full multi-cohort analysis has not yet been
reported.

## Requirements

- Linux with `micromamba`.
- A QIIME 2 amplicon 2024.10 environment named `q2-amplicon-2024.10`
  (override with `Q2_ENV`).
- `q2-periodontitis` installed in that environment.
- matplotlib for the validation figures (installed by `environment.yml`).
- Network access to ENA and SRA/BioSample.
- About 60 GB of free disk for the raw reads, 4 cores and 15 GB of RAM.
- About one day to process all 16 cohorts, plus the rarefaction robustness
  runs.

## Quick start

Run the full pre-specified analysis for every cohort in `cohorts.tsv`:

```bash
bash validation/run_all.sh /path/to/workdir
```

Process and analyse one cohort:

```bash
bash validation/run_cohort.sh SRP075100 /path/to/workdir
bash validation/analysis/analyse.sh /path/to/workdir /path/to/results
```

Each cohort has a configuration in `configs/<cohort>.sh`. `cohorts.tsv`
records its accession, design, independence from the SMDI derivation, and
the analyses it enters.

## Pipeline

- `download.sh` fetches runs from ENA and checks their MD5 sums.
- `prepare.py` builds standardised metadata from SRA/BioSample records.
- `read_qc.py` records primer positions and quality quartiles, which are
  used to fix trimming before any index is computed.
- `process_dataset.sh` runs DADA2 or Deblur, fits the naive Bayes and
  VSEARCH classifiers, and runs `qiime periodontitis collapse-complexes`.
- `analysis/collect.py` joins the per-sample indices, complex fractions
  and metadata, and writes `qc.tsv`.
- `analysis/run_analysis.py` runs the pre-specified Q1-Q6 analyses and
  writes `tables/`.
- `analysis/robustness.py` runs the pseudocount and rarefaction analyses.
- `analysis/performance.py` measures runtime and peak memory on synthetic
  sparse tables.
- `analysis/figures.py` writes `figures/`.
- `analysis/report.py` writes `REPORT.md`.

`run_all.sh` passes `validation/results` to `analyse.sh`; when analysing a
single cohort, pass another output directory if you do not want to add
files to the repository.

The files currently in `results/` are the legacy single-cohort outputs
described in [docs/VALIDATION.md](../docs/VALIDATION.md). The current
pipeline writes `samples/`, `qc.tsv`, `tables/` (including
`performance.tsv`), `performance_environment.json`, `figures/` and
`REPORT.md` to the results directory passed to `analyse.sh`. The
performance table is a synthetic engineering benchmark, not a biological
result.

## Reproducibility

- The random seed is fixed at `20261022` in `analysis/common.py`.
- Per-cohort software versions are written to
  `<cohort>/out/versions.txt` by `process_dataset.sh`.
- `.qza` and `.qzv` artifacts are not tracked in git. Archive them with the
  analysis data (for example on Zenodo) if provenance must be preserved.

## Current status

The multi-cohort pipeline is implemented, but its full 16-cohort results
are not committed. The validation reported in `docs/VALIDATION.md` is the
SRP075100 case study. See [docs/VALIDATION.md](../docs/VALIDATION.md) and
[docs/ANALYSIS_PLAN.md](../docs/ANALYSIS_PLAN.md) for the scientific scope
and the planned analyses.
