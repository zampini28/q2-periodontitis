# Validation on public clinical data

This page reports a case study on one independent public cohort
(SRP075100). The pre-specified multi-cohort validation is described in
[ANALYSIS_PLAN.md](ANALYSIS_PLAN.md); its results are not reported here
yet.

To reprocess SRP075100 with the current pipeline:

```bash
bash validation/run_cohort.sh SRP075100 /path/to/workdir
bash validation/analysis/analyse.sh /path/to/workdir /path/to/results
```

The full multi-cohort workflow is
`bash validation/run_all.sh /path/to/workdir`. See
[validation/README.md](../validation/README.md) for requirements and
outputs.

The per-subject values and the result tables described on this page are
the legacy single-cohort outputs in
[`validation/results/`](../validation/results/). The current pipeline
writes `samples/`, `qc.tsv`, `tables/`, `figures/` and `REPORT.md` to the
results directory passed to `analyse.sh`.

## Dataset

Chen et al. (2018, *ISME J* 12:1210–1224), SRA **SRP075100**. The data are
Illumina MiSeq reads of the 16S rRNA V4 region (515F/806R) from
subgingival plaque collected with paper points at two contralateral
maxillary posterior teeth per subject:

| Group | Samples | Subjects |
|---|---|---|
| Healthy (HP) | 42 | 21 |
| Periodontitis, before therapy (D1P) | 96 | 48 |
| Periodontitis, ≥ 4 weeks after non-surgical therapy (D2P) | 19 | 19 |

Groups, subjects and teeth are taken from the SRA sample titles
(`GROUP_SUBJECT_TOOTH…`). Clinical scores per sample are not in the SRA
records, so the analysis uses the diagnostic groups.

The dataset is **independent of the SMDI's derivation**. Chen, Marsh &
Al-Hebshi (2022) trained it on Abusleme 2013, Griffen 2012 and Kirst 2015,
and tested it on Pei 2020, Bizzarro 2016 and HMP.

## Processing

QIIME 2 amplicon 2024.10. The reads are already merged and primer-free
(about 253 nt). DADA2 `denoise-single` (`--p-trunc-len 250`) was used. A
naive Bayes classifier was trained on eHOMD 16S RefSeq V15.23 trimmed to
515F/806R. `collapse-complexes` ran with default parameters.

## Quality control

| Measure (157 samples) | Median | IQR |
|---|---|---|
| Reads per sample | 28,600 | 24,000–33,900 |
| Reads in a Socransky complex | 0.29 | 0.20–0.40 |
| Reads named to species | 0.83 | 0.75–0.89 |

Most reads are named to species, but the Socransky panel covers only
about 30% of the community. Species resolution in V4 is uneven across
genera:

| Genus | Reads named to species |
|---|---|
| *Campylobacter*, *Porphyromonas*, *Prevotella* | ≥ 0.94 |
| *Treponema*, *Veillonella*, *Capnocytophaga*, *Fusobacterium* | 0.76–0.87 |
| ***Streptococcus***, ***Actinomyces***, *Fretibacterium* | **0.27–0.29** |

The *Streptococcus* species of the **Yellow complex (mitis/oralis/sanguinis
group) are largely unresolvable in V4**, so the Yellow complex is
underestimated here (median 0 in periodontitis subjects). The Blue group is
less affected, because it is assigned at genus level: all of its features
here were assigned through the *Actinomyces* genus entry. This is a
property of the amplicon, not of the plugin. Prefer longer or more
variable regions (e.g. V1–V3) when the early-colonizer complexes matter.

## 1. Health vs periodontitis

Subjects are the unit of analysis: the two sites of each subject are
averaged. The 95% CI comes from a stratified bootstrap over subjects
(2000 resamples). AUC is P(periodontitis > health); values below 0.5 mean
higher in health.

| Measure | Median health | Median periodontitis | AUC (95% CI) | Mann-Whitney p |
|---|---|---|---|---|
| **SMDI** (49 species) | −0.07 | 2.43 | **0.90** (0.81–0.97) | 1.6 × 10⁻⁷ |
| SMDI_simplified | −3.68 | 1.26 | 0.86 (0.76–0.95) | 1.8 × 10⁻⁶ |
| PC2014_moderate | 0.015 | 0.095 | 0.84 (0.71–0.94) | 1.1 × 10⁻⁵ |
| Red complex | 0.002 | 0.099 | 0.86 (0.75–0.95) | 2.3 × 10⁻⁶ |
| Orange complex | 0.121 | 0.162 | 0.71 (0.56–0.84) | 0.0065 |
| Yellow complex | 0.003 | 0 | 0.18 (0.08–0.31) | 1.7 × 10⁻⁵ |
| Green complex | 0.008 | 0.005 | 0.38 (0.25–0.52) | 0.13 |
| Purple complex | 0.012 | 0.002 | 0.21 (0.09–0.35) | 1.5 × 10⁻⁴ |
| Blue (*Actinomyces*) | 0.047 | 0.009 | 0.16 (0.07–0.25) | 6.6 × 10⁻⁶ |

The complexes behave as the Socransky model predicts. The late
colonizers (Red, Orange) are higher in periodontitis. The early colonizers
(Yellow, Purple, Blue) are lower. Green is lower too, but not
significantly.

**The full SMDI reproduces its published accuracy.** Chen et al. (2022)
reported AUC 0.92 in their test sets, against 0.90 here. Its published
cut-offs give:

| Cut-off | Unit | Sensitivity | Specificity | Accuracy |
|---|---|---|---|---|
| 0 | subjects | 0.90 | 0.52 | 0.78 |
| 0.4 (best in training) | subjects | 0.88 | 0.76 | 0.84 |
| **0.5 (best in test)** | **subjects** | **0.85** | **0.90** | **0.87** |
| 0.5 | samples | 0.79 | 0.74 | 0.78 |

The cut-off that maximizes Youden's J in this dataset is 0.54, which
practically coincides with the published 0.5. This is an independent
confirmation of the threshold.

For SMDI_simplified at its published cut-off of 0: subjects 0.60 / 0.90 /
0.70 (sensitivity / specificity / accuracy), samples 0.62 / 0.88 / 0.70.
That is specific but less sensitive. This is consistent with the authors'
note that V4 data place the SMDI about one unit lower than V1–V3. The
in-sample Youden cut-off is −1.46 (0.85 / 0.81). It was chosen on the
same data, so it is descriptive only.

**16 of the 49 SMDI species are absent from this V4 table**, including the
top health-associated *S. sanguinis*, *A. naeslundii* and *S. mitis*, and
*Fretibacterium* HMT 359/360/362. The reason is the poor species
resolution of these genera in V4 (see Quality control). Absent taxa enter
at the pseudocount, as the method prescribes. The index still performed
best, but on other data it can behave differently. The plugin lists the
absent taxa in a warning.

**The indices largely agree.** Spearman ρ, n = 69 subjects: SMDI vs
SMDI_simplified 0.91, SMDI vs PC2014_moderate 0.73.

## 2. Response to non-surgical therapy

Same tooth before (D1P) and after (D2P) therapy, one site per subject, n =
19. Change = after − before.

| Index | Median change | Sites that decreased | Wilcoxon signed-rank p |
|---|---|---|---|
| SMDI | −0.25 | 13/19 | 0.096 |
| SMDI_simplified | −0.57 | 12/19 | 0.29 |
| PC2014_moderate | −0.019 | 13/19 | 0.13 |

All indices tend to decrease after therapy, but no change is significant
with 19 sites. The original study also described heterogeneous site-level
responses. These data neither confirm nor rule out responsiveness to
treatment. A larger longitudinal dataset is needed.

## Conclusions and limits

* **The 49-species SMDI is the best-performing and best-supported index.**
  In independent V4 data it reached AUC 0.90, and its published cut-off of
  0.5 gave sensitivity 0.85 and specificity 0.90 at subject level. Use it
  as the primary measure.
* SMDI_simplified, PC2014_moderate and the Red complex also discriminate
  (AUC 0.84–0.86). SMDI_simplified needs only genus names, which makes it
  the fallback when species resolution is poor.
* The early-colonizer complexes and several SMDI species depend on
  species resolution, which V4 does not provide for oral streptococci and
  *Actinomyces*. Always report the `qc_*` columns and the SMDI
  absent-taxa warning.
* One dataset, one region, one country. Diagnostic accuracy should be
  re-estimated in each study population, and neither index is a clinical
  test.
