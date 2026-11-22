# Validation workflow

This directory contains the code for the pre-specified multi-cohort
validation described in [docs/ANALYSIS_PLAN.md](../docs/ANALYSIS_PLAN.md).
Its results, from the run of 21 November 2026, are in [results/](results/)
and are reported in [docs/VALIDATION.md](../docs/VALIDATION.md).

## Requirements

- Linux with `micromamba`.
- A QIIME 2 amplicon 2024.10 environment named `q2-amplicon-2024.10`
  (override with `Q2_ENV`).
- `q2-periodontitis` installed in that environment.
- An R environment with microbiome 1.6.0 named `r-microbiome-1.6.0`
  (override with `R_ENV`), from `r-microbiome-1.6.0.yml`.
- matplotlib for the validation figures (installed by `environment.yml`).
- Network access to ENA and SRA/BioSample.
- About 60 GB of free disk for the raw reads plus about 12 GB for the
  environments, 4 cores and 15 GB of RAM.
- About a day on 4 cores to process all 16 cohorts and run the analysis
  (12 hours on 48 cores in the reported run, mostly downloads).

## Quick start

From a fresh clone on a Linux machine (only `curl` and `bzip2` needed),
one command installs micromamba, the QIIME 2 environment and the R
microbiome 1.6.0 environment, runs the tests, processes every cohort and
writes all results (multi-cohort validation, robustness, simulation,
performance, R cross-check and plugin-vs-R benchmark), figures and
`results/REPORT.md`:

```bash
git clone https://github.com/LaBiOmicS/q2-periodontitis.git
cd q2-periodontitis
tmux new -s validation 'bash validation/reproduce.sh ~/q2-validation-data 2>&1 | tee ~/reproduce.log'
```

It can be re-run after an interruption; finished steps are skipped.

The environments and package cache take about 12 GB in
`MAMBA_ROOT_PREFIX` (default `~/micromamba`); set it, and the data
directory, to a larger disk if home is small, for example
`MAMBA_ROOT_PREFIX=/data/micromamba bash validation/reproduce.sh /data/q2`.
The tests, including the R cross-check, run before any download and stop
the command if they fail. On a shared cluster, submit the command as a
batch job rather than running it for a day on a login node; downloads from
conda-forge, bioconda, ENA and SRA must be reachable (set `https_proxy`
if needed).

If the environment is already installed, run the full pre-specified
analysis for every cohort in `cohorts.tsv`:

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
- `analysis/simulation.py` runs the simulated known-truth validation.
- `analysis/rcheck.py` cross-checks the SMDI of the plugin against
  `analysis/smdi_microbiome.R` (R microbiome 1.6.0, the implementation used
  to define the index) on the simulated, synthetic and cohort tables, and
  benchmarks the two implementations (`r_crosscheck.tsv`,
  `r_benchmark.tsv`, `fig7_r_benchmark`).
- `analysis/figures.py` writes `figures/`.
- `analysis/report.py` writes `REPORT.md`.

`run_all.sh` passes `validation/results` to `analyse.sh`; when analysing a
single cohort, pass another output directory if you do not want to add
files to the repository.

The pipeline writes `samples/`, `qc.tsv`, `tables/` (including
`performance.tsv` and `simulation.tsv`), `performance_environment.json`,
`figures/` and `REPORT.md` to the results directory passed to
`analyse.sh`. The performance and simulation tables are synthetic
engineering and validation results, not biological findings.
`results/provenance/` holds the per-cohort software versions and the note
on the one file taken from NCBI SRA.

## Reproducibility

- The random seed is fixed at `20261022` in `analysis/common.py`; the
  simulation uses `20261104` in `analysis/simulation.py`.
- Per-cohort software versions are written to
  `<cohort>/out/versions.txt` by `process_dataset.sh`.
- `.qza` and `.qzv` artifacts are not tracked in git. Archive them with the
  analysis data (for example on Zenodo) if provenance must be preserved.

## Current status

The full run of all 16 cohorts finished on 21 November 2026 (12 hours on a
48-core server). Its outputs are in `results/` and are reported in
[docs/VALIDATION.md](../docs/VALIDATION.md), including the deviations from
[docs/ANALYSIS_PLAN.md](../docs/ANALYSIS_PLAN.md). The R cross-check and
benchmark were rerun on the same day after a fix to the memory
measurement. On 22 November the analysis steps (`collect.py` to
`report.py`, except the R cross-check, performance and simulation) were
rerun after PRJNA1294584's controls were relabelled as gingivitis.
