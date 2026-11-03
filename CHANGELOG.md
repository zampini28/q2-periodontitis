# Changelog

All notable changes to `q2-periodontitis` are documented in this file.

## [Unreleased] — documentation and CI

- Fix stale validation references and add `validation/README.md`.
- Clarify that the reported 69-subject validation is the SRP075100 case
  study, not the multi-cohort analysis.
- Restore the canonical MIT license text; the previous file had a
  duplicated clause and GitHub detected it as "Other".
- Run CI on every push and pull request, not only `main` and `master`.
- Add separate pytest-cov coverage reporting for the plugin and the
  validation analysis, with GitHub-native coverage badges.

## [Unreleased] — scientific revision

Every reported quantity is now either a published index or a documented
derivation from published complexes. The changes were checked against
public clinical data (see `docs/VALIDATION.md`).

### ⚠️ Breaking changes
- **Removed `SWL`, `KRS` and `pMDI_ilr`.** The SWL weights and the KRS
  threshold and sigmoid had no published basis. `pMDI_ilr` was dominated by
  the pseudocount when complexes were absent. `--p-keystone-threshold` and
  `--p-keystone-steepness` were removed.
- **Removed the "Extended" complex.** It was not part of the Socransky
  model. The moderate-evidence taxa of Pérez-Chaparro et al. (2014) are now
  reported separately as `PC2014_moderate` and flagged in the mapping
  (`PC2014` column).
- **Species-only reference** (`data/socransky_complexes.tsv` replaces
  `ehomd_map.tsv`). Only the species of Socransky et al. (1998) and
  Socransky & Haffajee (2002, Fig. 1 and Fig. 11) are assigned, and each
  row cites its figure. The only genus-level entry is *Actinomyces*
  (Blue), which the 2002 figure labels as a genus. Removed: every other
  genus entry (*Streptococcus*, *Fusobacterium*, *Veillonella*,
  *Capnocytophaga*, …), the unreferenced species (*S. mutans*, *Rothia*,
  *Neisseria*, *Corynebacterium*, …) and *A. actinomycetemcomitans*, whose
  complex depends on a serotype that 16S cannot resolve. Added:
  *F. periodonticum* and *F. nucleatum* subsp. *vincentii* (Orange,
  2002 Fig. 11).
- **`--p-assignment-mode` removed.** Assignment always uses the species,
  and falls back to the genus only for explicit genus-level entries.
- **Removed `pMDI` and `--p-pseudocount`.** The index was never
  published or validated, and in the validation data it was almost
  redundant with the SMDI (Spearman ρ 0.86–0.90). The SMDI is now the only
  dysbiosis index. The complex table is unchanged.
- **`delta-indices`** averages several samples of a unit within a group
  instead of raising an error, and skips `qc_*` columns.
- **`summarize`** no longer tests every categorical metadata column. It
  tests only `--p-group-column`, optionally after averaging per
  `--p-subject-column`.

### Added
- `--p-smdi-pseudocount` for `collapse-complexes`: an optional fixed
  pseudocount for the SMDI CLR, which makes each sample's SMDI
  independent of the rest of the table. When it is not set, the
  published, table-dependent pseudocount is used, as before.
- **SMDI** (Chen, Marsh & Al-Hebshi 2022), computed exactly as in the
  original implementation (R *microbiome* 1.6.0 CLR). Both are built in
  and always computed: `SMDI`, the 49 species with MDG ≥ 0.4 from the
  authors' Data Set 1, and `SMDI_simplified`, *Fretibacterium* and
  *Treponema* vs *Actinomyces*. Other signatures can be supplied as the new
  `DysbiosisSignature` type (`--i-smdi-signature`, column
  `SMDI_custom`).
- `PC2014_moderate`: fraction of reads from the 17 moderate-evidence taxa,
  matched by name or HMT number.
- Quality-control columns `qc_depth`, `qc_assigned_fraction` and
  `qc_species_fraction`.
- Mapping columns `Genus`, `Species` and `PC2014`, and the `ambiguous`
  resolution for multi-species labels that span complexes.
- Effect sizes in `summarize` (AUC for two groups, epsilon-squared for
  more).
- Parser support for eHOMD labels: `[Eubacterium] nodatum`,
  `Peptostreptococcaceae [G-6]`, `sp. HMT 360`, subspecies, and 16S
  clades (`infantis_clade_431`), which are merged into their species.
  Phylotypes are matched by HMT number.
- `validation/`: a reproducible re-analysis of SRA SRP075100 (Chen et al.
  2018) with `docs/VALIDATION.md`.
- Citations: Chen et al. (2022), Chen et al. (2018).

## [Unreleased] — review corrections

### ⚠️ Breaking changes
- **`pMDI_clr` removed.** It was exactly `pMDI_ilr / sqrt(6/5)`. Use `pMDI_ilr`, the canonical compositional balance.
- **Pseudocount.** The per-sample dynamic pseudocount (half of the smallest non-zero complex abundance) was replaced by the `--p-pseudocount` parameter (default `1e-5`), applied identically to all samples. `pMDI` and `pMDI_ilr` values change.
- **KRS formula.** Now `(1 + tanh(alpha * (Red - theta) / theta)) * (Red + Extended + 0.5 * Orange)`, so `theta` (`--p-keystone-threshold`) is a real threshold (factor = 1 at `Red == theta`). KRS values change. `--p-keystone-threshold` must be > 0.
- **Confidence handling.** Abundances are no longer multiplied by the taxonomy `Confidence` (which lost mass, so complex profiles did not sum to 1). Use the new `--p-min-confidence` to send low-confidence features to `Other`.
- **`--p-assignment-mode probabilistic` removed** (it was never implemented and behaved like `genus-majority`).
- **`--i-custom-mapping` now takes a `SocranskyReference` artifact** (TSV with a mandatory `taxon`/`complex` header) instead of `FeatureData[SocranskyMapping]`. Invalid complex names raise a clear error.
- **`delta-indices`:**
  - `--p-sample-id-col` renamed to `--p-subject-col`.
  - Output type is now `SampleData[PeriodontalIndicesDelta]` (indexed by subject ID, column `subject-id`) instead of `PeriodontalIndices`.
  - Deltas are computed for every index in the input (pMDI, pMDI_ilr, SWL, KRS).
  - Subjects with more than one sample per group, or no paired subjects at all, now raise an error instead of silently using the first sample or returning an empty table.
- **Reference map changes** (`ehomd_map.tsv` / `synonyms_map.tsv`), which change complex assignments:
  - *Actinomyces odontolyticus* / *Schaalia odontolytica*: Blue → **Purple** (Socransky et al., 1998).
  - Genus-level entries for *Porphyromonas* and *Prevotella* removed (species entries kept); unmatched species of these genera are now `Other`.
  - Wrong synonyms removed (*Campylobacter sputorum* → *C. rectus*, *Eubacterium brachy* → *Mogibacterium timidum*, *Streptococcus mitis oralis*, phylum names); added *Wolinella recta* → *C. rectus* and *Eubacterium timidum* → *M. timidum*.
- **Taxonomy parsing.** SILVA/GTDB/Greengenes2 suffixes (`Prevotella_7`, `Fusobacterium_A`, `nucleatum_D`) are stripped and placeholder names (`uncultured_bacterium`, `metagenome`, `unidentified`, …) ignored, so more features are assigned than before.
- **Packaging.** `setup.py` removed (everything is in `pyproject.toml`); the `scikit-bio` dependency was dropped.

### Added
- `--p-min-confidence` and `--p-pseudocount` parameters for `collapse-complexes`.
- `SocranskyReference` and `PeriodontalIndicesDelta` semantic types with formats and transformers.
- `source` column in `ehomd_map.tsv` (`Socransky1998`, a `citations.bib` key, or `unverified`).
- Citations: Socransky & Haffajee (2002); Pérez-Chaparro et al. (2014).
- Benjamini-Hochberg *q*-values in `summarize`.
- Format validation for complexes, indices, deltas and references.
- Test suite based on `qiime2.plugin.testing.TestPluginBase` with hand-computed expectations.
- `CHANGELOG.md`.

### Fixed
- `summarize` never displayed statistics: `numpy` was not imported and the resulting `NameError` was swallowed by a bare `except`. Failures now emit a warning with the reason.
- `summarize` only tests categorical metadata columns.
- `summarize` works offline: Vega 5.30.0, Vega-Lite 5.21.0 and Vega-Embed 6.26.0 are bundled instead of loaded from a CDN.
- Sample IDs containing `</script>` can no longer break the report.
- `delta-indices`: metadata columns named like an index no longer collide; numeric group columns work.
- Incorrect metadata of the Curtis et al. (2020), Chen et al. (2010) and Olsen & Yilmaz (2019) entries in `citations.bib`.
- Collapsing no longer densifies the feature table; indices are computed vectorized.

### Changed
- CI and `environment.yml` use the official QIIME 2 amplicon 2024.10 environment.
- `espec.md` moved to `docs/SPEC.md`.
- Author e-mail updated.
