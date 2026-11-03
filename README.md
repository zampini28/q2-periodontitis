# `q2-periodontitis`

**QIIME 2 plugin that summarizes subgingival 16S rRNA profiles as Socransky microbial complexes and dysbiosis indices**

[![CI](https://github.com/LaBiOmicS/q2-periodontitis/actions/workflows/ci.yml/badge.svg)](https://github.com/LaBiOmicS/q2-periodontitis/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![QIIME 2](https://img.shields.io/badge/QIIME%202-2024.10%2B-purple.svg)](https://qiime2.org)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB.svg)](https://www.python.org)

`q2-periodontitis` takes a feature table and its taxonomy and reports, for every sample:

* the relative abundance of the **Socransky complexes**. Only species documented by Socransky et al. (1998) and Socransky & Haffajee (2002) are assigned; every entry of the reference cites its source.
* the **SMDI** (Subgingival Microbial Dysbiosis Index; Chen, Marsh & Al-Hebshi, 2022): the full 49-species index and the simplified three-genus index, both computed as in the original publication.
* **PC2014_moderate**, the fraction of reads from the 17 taxa with moderate evidence of association with periodontitis (Pérez-Chaparro et al., 2014).
* **quality-control measures**: sequencing depth, the fraction of reads assigned to a complex, and the fraction of reads named to species.

The scientific basis, formulas and limits are in **[docs/THEORY.md](docs/THEORY.md)**. **[docs/VALIDATION.md](docs/VALIDATION.md)** reports a case study on one independent public cohort (SRP075100). On its 69 subjects, the 49-species SMDI separated periodontitis from health with AUC 0.90, and the published cut-off of 0.5 gave sensitivity 0.85 and specificity 0.90. The pre-specified multi-cohort validation is described in **[docs/ANALYSIS_PLAN.md](docs/ANALYSIS_PLAN.md)**; its results are not yet reported. See also [docs/TUTORIAL.md](docs/TUTORIAL.md), [docs/SPEC.md](docs/SPEC.md) and the [CHANGELOG](CHANGELOG.md).

> **Research use only.** The indices summarize microbial composition. They are not diagnostic tests and do not replace clinical periodontal examination.

---

## 1. Installation

Developed and tested with the official **QIIME 2 amplicon 2024.10** distribution:

```bash
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml
conda env create -n q2-amplicon-2024.10 --file qiime2-amplicon-2024.10-py310-linux-conda.yml
conda activate q2-amplicon-2024.10
pip install git+https://github.com/LaBiOmicS/q2-periodontitis.git
qiime dev refresh-cache
```

For development:

```bash
git clone https://github.com/LaBiOmicS/q2-periodontitis.git
cd q2-periodontitis
conda env update -n q2-amplicon-2024.10 --file environment.yml   # adds pytest
pip install -e .
pytest
```

## 2. Recommended input

* A feature table of **raw counts** (do not rarefy or normalize).
* A taxonomy classified to **species** against an **oral** reference, ideally eHOMD. [validation/README.md](validation/README.md) documents the validation workflow, including how to train a naive Bayes classifier on eHOMD for your amplicon region. With generic databases, fewer reads reach species level, and the complex profiles then cover a smaller part of the community. The `qc_*` columns show how much.

## 3. Usage

### 3.1 `collapse-complexes`

```bash
qiime periodontitis collapse-complexes \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --o-complex-table complexes.qza \
  --o-indices indices.qza \
  --o-mapping mapping.qza
```

| Option | Type | Description |
| :--- | :--- | :--- |
| `--i-table` | `FeatureTable[Frequency]` | Raw counts. |
| `--i-taxonomy` | `FeatureData[Taxonomy]` | `Taxon` and, optionally, `Confidence`. |
| `--i-custom-mapping` | `SocranskyReference` (optional) | Extra or overriding taxon→complex entries (`taxon`, `complex`, optional `level` = `species`/`genus`). |
| `--i-smdi-signature` | `DysbiosisSignature` (optional) | An additional signature (`taxon`, `direction` = `dysbiotic`/`normobiotic`, `level`), computed with the SMDI method as `SMDI_custom`. |
| `--p-min-confidence` | `Float` 0–1 (default 0) | Features below this classifier confidence go to `Other` (`low-confidence`) and are not counted as SMDI or PC2014 taxa. |
| `--p-smdi-pseudocount` | `Float` in (0, 1) (default: not set) | Fixed pseudocount for the SMDI CLR. When it is not set, the published definition is used: half the smallest non-zero relative abundance of the table. A fixed value makes each sample's SMDI independent of the other samples of the table. |

**Outputs**

* `complexes` (`SocranskyComplexes`): Red, Orange, Yellow, Green, Purple, Blue and Other per sample; rows sum to 1.
* `indices` (`PeriodontalIndices`): `SMDI`, `SMDI_simplified`, `SMDI_custom` (with a signature), `PC2014_moderate`, `qc_depth`, `qc_assigned_fraction`, `qc_species_fraction`. Undefined values (e.g. empty samples) are NaN. A warning lists the SMDI taxa absent from the table, since they enter the index as zeros.
* `mapping` (`SocranskyMapping`): for each feature, the parsed genus and species, the complex, the resolution (`species`, `genus`, `ambiguous`, `unassigned`, `low-confidence`) and whether it belongs to the Pérez-Chaparro list.

**Custom signatures.** To compute another published or in-house signature with the SMDI method:

```text
taxon	direction	level
Treponema denticola	dysbiotic	species
Streptococcus sanguinis	normobiotic	species
```

```bash
qiime tools import --type DysbiosisSignature \
  --input-path my-signature.tsv --output-path my-signature.qza
qiime periodontitis collapse-complexes ... --i-smdi-signature my-signature.qza
```

### 3.2 `summarize`

```bash
qiime periodontitis summarize \
  --i-complexes complexes.qza --i-indices indices.qza --i-mapping mapping.qza \
  --m-metadata-file metadata.tsv \
  --p-group-column status \
  --p-subject-column subject \
  --o-visualization summary.qzv
```

The report shows the complex profiles, the distribution of each index by group, and the quality-control summary. Tests run only for `--p-group-column`: Mann-Whitney U with AUC for two groups, Kruskal-Wallis with epsilon-squared for more, and Benjamini-Hochberg q-values. Set `--p-subject-column` whenever a subject contributed several samples (sites, visits). Their values are then averaged, so each subject counts once. The report works offline (Vega libraries are bundled).

### 3.3 `delta-indices`

```bash
qiime periodontitis delta-indices \
  --i-indices indices.qza --m-metadata-file metadata.tsv \
  --p-subject-col subject --p-group-col visit \
  --p-baseline-group baseline --p-target-group follow-up \
  --o-delta-table deltas.qza
```

Gives follow-up minus baseline for every index, one row per subject. For site-level change, pass a column that combines subject and site as `--p-subject-col`. Several samples of one unit in a group are averaged, and units without both visits are skipped.

## 4. Semantic types

| Type | Content | File |
| :--- | :--- | :--- |
| `SocranskyComplexes` | Complex relative abundances per sample | `complexes.tsv` |
| `PeriodontalIndices` | Indices and QC columns per sample | `indices.tsv` |
| `SampleData[PeriodontalIndicesDelta]` | Per-subject deltas (`subject-id`) | `deltas.tsv` |
| `SocranskyMapping` | Per-feature assignment audit | `mapping.tsv` |
| `SocranskyReference` | Custom taxon→complex entries | `reference.tsv` |
| `DysbiosisSignature` | SMDI signature taxa | `signature.tsv` |

Complexes, indices, deltas and mapping can be used as QIIME 2 Metadata, e.g. with `qiime metadata tabulate` or `qiime longitudinal linear-mixed-effects --p-metric SMDI_simplified`.

## 5. Citations

Run `qiime periodontitis --citations` for BibTeX. The key references are:

1. Socransky SS, Haffajee AD, Cugini MA, Smith C, Kent RL Jr (1998). Microbial complexes in subgingival plaque. *J Clin Periodontol* 25:134–144. https://doi.org/10.1111/j.1600-051X.1998.tb02419.x
2. Socransky SS, Haffajee AD (2002). Dental biofilms: difficult therapeutic targets. *Periodontol 2000* 28:12–55. https://doi.org/10.1034/j.1600-0757.2002.280102.x
3. Chen T, Marsh PD, Al-Hebshi NN (2022). SMDI: an index for measuring subgingival microbial dysbiosis. *J Dent Res* 101:331–338. https://doi.org/10.1177/00220345211035775
4. Pérez-Chaparro PJ, Gonçalves C, Figueiredo LC, et al. (2014). Newly identified pathogens associated with periodontitis: a systematic review. *J Dent Res* 93:846–858. https://doi.org/10.1177/0022034514542468

Validation data: Chen C, Hemme C, Beleno J, et al. (2018). Oral microbiota of periodontal health and disease and their changes after nonsurgical periodontal therapy. *ISME J* 12:1210–1224. https://doi.org/10.1038/s41396-017-0037-1

## 6. License

MIT ([LICENSE](LICENSE)). The bundled Vega, Vega-Lite and Vega-Embed libraries are BSD-3-Clause (`q2_periodontitis/assets/web/js/`).
