# Tutorial

This tutorial runs `q2-periodontitis` on a subgingival 16S rRNA dataset in
QIIME 2 (amplicon 2024.10). [validation/README.md](../validation/README.md)
documents the validation workflow, and [validation/run_all.sh](../validation/run_all.sh)
runs it on public data from raw reads to results.

## 1. Inputs

You need:

1. **A feature table of raw counts** (`FeatureTable[Frequency]`, e.g. from
   DADA2). Do not rarefy it: the indices are log-ratios of counts.
2. **A species-level taxonomy** (`FeatureData[Taxonomy]`). Use a classifier
   trained on the eHOMD 16S RefSeq, trimmed to your primers:

   ```bash
   qiime feature-classifier extract-reads --i-sequences ehomd-refseq.qza \
     --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer GGACTACNVGGGTWTCTAAT \
     --o-reads ehomd-v4.qza                       # 515F/806R (V4) example
   qiime feature-classifier fit-classifier-naive-bayes \
     --i-reference-reads ehomd-v4.qza --i-reference-taxonomy ehomd-taxonomy.qza \
     --o-classifier ehomd-v4-classifier.qza
   qiime feature-classifier classify-sklearn --i-classifier ehomd-v4-classifier.qza \
     --i-reads rep-seqs.qza --o-classification taxonomy.qza
   ```

   SILVA, GTDB, Greengenes2 and NCBI strings are also parsed. Database
   suffixes (`Prevotella_7`, `Fusobacterium_A`) and placeholders
   (`uncultured_bacterium`) are handled. These databases usually name fewer
   oral species, though.
3. **Sample metadata** with the grouping you want to compare and, if
   subjects gave more than one sample, a subject column:

   ```text
   sample-id	subject	site	status	visit
   #q2:types	categorical	categorical	categorical	categorical
   P1-16-BL	P1	P1-16	periodontitis	baseline
   P1-16-FU	P1	P1-16	periodontitis	follow-up
   H1-3-BL	H1	H1-3	health	baseline
   ```

## 2. Collapse and compute the indices

```bash
qiime periodontitis collapse-complexes \
  --i-table table.qza --i-taxonomy taxonomy.qza \
  --o-complex-table complexes.qza --o-indices indices.qza --o-mapping mapping.qza
```

Check the quality-control columns first:

```bash
qiime metadata tabulate --m-input-file indices.qza --o-visualization indices.qzv
```

* `qc_species_fraction`: fraction of reads named to a species. When it is
  low, complex abundances and species-level indices cover little of the
  community.
* `qc_assigned_fraction`: fraction of reads in any Socransky complex.
* `qc_depth`: reads per sample. Very shallow samples give noisy
  log-ratios.

`mapping.qza` lists, for each feature, the genus and species the plugin
read from the taxonomy, its complex and the resolution of the assignment.
Use it to audit unexpected assignments.

`collapse-complexes` warns about SMDI signature taxa absent from the
table. A few absences are normal. If the top health-associated species
(*S. sanguinis*, *A. naeslundii*) are missing, your region or classifier
does not resolve them, and `SMDI_simplified` (genus level) is the more
robust index.

### Optional: another signature

```bash
qiime tools import --type DysbiosisSignature \
  --input-path my-signature.tsv --output-path my-signature.qza
qiime periodontitis collapse-complexes \
  --i-table table.qza --i-taxonomy taxonomy.qza --i-smdi-signature my-signature.qza \
  --o-complex-table complexes.qza --o-indices indices.qza --o-mapping mapping.qza
```

`my-signature.tsv` has the columns `taxon`, `direction` and `level`. The
result is reported as `SMDI_custom`.

### Optional: custom complex assignments

```bash
printf 'taxon\tcomplex\tlevel\nSelenomonas noxia\tGreen\tspecies\n' > extra.tsv
qiime tools import --type SocranskyReference \
  --input-path extra.tsv --output-path extra.qza
# then: --i-custom-mapping extra.qza
```

Custom entries are your own hypotheses, not part of the Socransky model.
Report them as such.

## 3. Compare groups

```bash
qiime periodontitis summarize \
  --i-complexes complexes.qza --i-indices indices.qza --i-mapping mapping.qza \
  --m-metadata-file metadata.tsv \
  --p-group-column status --p-subject-column subject \
  --o-visualization summary.qzv
```

Restrict the metadata to the samples you want to compare, e.g. baseline
only. Otherwise a subject's baseline and follow-up samples would be pooled
into one average. For more complex designs (covariates, repeated measures),
export the indices and model them, for example with `qiime longitudinal
linear-mixed-effects --p-metric SMDI_simplified`.

## 4. Treatment response

```bash
qiime periodontitis delta-indices \
  --i-indices indices.qza --m-metadata-file metadata.tsv \
  --p-subject-col site --p-group-col visit \
  --p-baseline-group baseline --p-target-group follow-up \
  --o-delta-table deltas.qza
qiime metadata tabulate --m-input-file deltas.qza --o-visualization deltas.qzv
```

Negative deltas mean the index decreased (less dysbiotic). Test whether the
median change differs from zero with a paired test, such as the Wilcoxon
signed-rank test in `validation/analysis/run_analysis.py` (Q4). Use one
row per subject, or per site if sites are the unit you analyse.

## 5. Reporting

State the amplicon region, reference database and version, classifier,
the plugin version, and the median `qc_species_fraction`. Report SMDI
values as computed within one run (see [THEORY.md](THEORY.md)). If you used
a custom reference or signature, include it as supplementary material.
