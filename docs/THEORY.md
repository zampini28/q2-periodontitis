# Scientific basis of q2-periodontitis

This document states what each output measures, where every number and
taxon list comes from, and what the plugin cannot do. Every reference
named here is in `q2_periodontitis/citations.bib`.

## 1. Socransky complexes

### 1.1 What the complexes are

Socransky et al. (1998) measured 40 cultivable subgingival species by
checkerboard DNA-DNA hybridization in 13,261 plaque samples from 185
adults. They grouped the species by how often they co-occurred, using
cluster analysis and ordination. Five complexes were consistently
recovered. The Red complex, *Porphyromonas gingivalis*, *Tannerella
forsythia* and *Treponema denticola*, was strongly associated with pocket
depth and bleeding on probing. Socransky & Haffajee (2002, Fig. 1) later
drew the complexes as a pyramid with the *Actinomyces* species (blue) as a
sixth group. In that drawing the *Actinomyces*, Yellow, Green and Purple
groups are early colonizers, Orange is the bridge, and Red is the
late-stage group.

The complexes therefore describe a **fixed panel of 40 species** and how
they co-occur. They are not a classification of the whole oral microbiome.

### 1.2 How the plugin applies them to amplicon data

`q2_periodontitis/data/socransky_complexes.tsv` lists every reference
entry, with its source and the figure or sentence that supports it. The
rules are:

* **Only species from the published panel are assigned.** Membership comes
  from the 1998 abstract and from Socransky & Haffajee (2002), Fig. 1 (the
  complexes) and Fig. 11 (the 40-species panel ordered by complex). Current
  names are used, and historical names are resolved through
  `synonyms_map.tsv` (*Bacteroides forsythus* → *Tannerella forsythia*,
  *Peptostreptococcus micros* → *Parvimonas micra*, *Actinomyces
  odontolyticus* → *Schaalia odontolytica*, *Streptococcus sanguis* → *S.
  sanguinis*, …).
* **Genus-level extrapolation is not used**, with one exception. The 2002
  figure labels the blue group "*Actinomyces* species", without naming
  species. That is the only genus-level entry: unlisted *Actinomyces*
  species are assigned to Blue with resolution `genus`. *Schaalia
  odontolytica* (formerly *A. odontolyticus*) is Purple, following the
  original assignment.
* **Taxa the panel cannot represent are left in `Other`:**
  * *Aggregatibacter actinomycetemcomitans*. Serotype a is in the Green
    complex and serotype b is an outlier (1998). 16S rRNA genes do not
    distinguish serotypes.
  * *Fusobacterium nucleatum* subsp. *animalis*, which was not on the
    panel. *F. nucleatum* without a subspecies is assigned to Orange,
    because every subspecies on the panel is Orange.
  * The unnamed-species probes "*Streptococcus* sp." and "*Capnocytophaga*
    sp." of Fig. 11. Each probe targets one organism, not a whole genus,
    and the organism cannot be identified from the figure.
* **Multi-species labels** such as "*Streptococcus mitis/oralis*" are
  assigned only when every candidate falls in the same complex. Otherwise
  the feature is `ambiguous` and goes to `Other`.

### 1.3 Limits you must keep in mind

* **16S resolution.** Short amplicons often cannot separate panel species
  from close relatives. Examples are *S. mitis* / *S. oralis* /
  *S. pneumoniae* and several *Campylobacter* and *Capnocytophaga* species.
  A complex abundance is only as good as the species call behind it. The
  indices table reports `qc_species_fraction` (reads named to species) and
  `qc_assigned_fraction` (reads in any complex) for each sample, and the
  visualization shows them. Read them before interpreting complex
  profiles.
* **Different denominators.** Checkerboard data measure 40 species.
  Amplicon data measure everything. Complex proportions computed by this
  plugin are fractions of **all** reads, so they are not comparable in
  absolute value with the checkerboard percentages of the original papers.
* **Reference database.** Use a classifier trained on an oral reference
  (eHOMD) whenever possible. Generic databases name fewer oral species and
  assign more reads to placeholder taxa.

## 2. Indices

All indices are computed from **raw read counts**. Do not rarefy, and do
not supply relative abundances.

### 2.1 SMDI (Chen, Marsh & Al-Hebshi 2022)

The Subgingival Microbial Dysbiosis Index was derived by random-forest
selection of species that discriminate periodontitis from health in
published 16S datasets. It was validated on independent test sets (49
species: AUC 0.96 training, 0.92 test). It is defined as:

$$\mathrm{SMDI} = \overline{\mathrm{clr}}(\text{dysbiotic taxa}) - \overline{\mathrm{clr}}(\text{normobiotic taxa})$$

The paper's Methods state the software used for the CLR: "the 'transform'
function in R package 'microbiome' (version 1.6.0)". The plugin follows
the behaviour of that function, read from its source, and
`validation/analysis/rcheck.py` checks the plugin against microbiome 1.6.0
itself (see [VALIDATION.md](VALIDATION.md#cross-check-against-r-microbiome-160)):

1. Counts are aggregated to the signature's taxonomic level (species or
   genus). Features without a name at that level keep their reads as
   separate units, so sample totals do not change.
2. Counts are closed to relative abundances.
3. **If the table contains a zero**, one pseudocount, equal to **half the
   smallest non-zero relative abundance of the whole table**, is added to
   **every** entry. A table without zeros gets no pseudocount. In
   microbiome 1.6.0:
   `if (any(xt == 0)) { minval <- min(v[v > 0])/2; xt <- xt + minval }`.
4. `clr(x) = ln(x) - mean(ln(x))` is computed over all taxa of the sample.

**Why the code and not the paper's wording.** The paper's text, and the
guide of the authors' web calculator
(<https://bioinformatics.forsyth.org/smdi/>), describe the pseudocount as
added to zero entries only. The plugin follows the microbiome 1.6.0
behaviour instead, because that is the computation the paper's Methods
report, and so the one behind the published index values and cut-offs.
The two conventions are not equivalent: on sparse tables they differ by
up to about 0.2 SMDI units (simulated tables of `validation/analysis`),
which matters next to the published cut-offs of 0–0.5. Which convention
the web calculator uses is tested in `validation/calculator/`.

Two consequences follow from the original definition.

* SMDI values depend on the table they are computed in, through the
  pseudocount: it is taken from the smallest non-zero relative abundance
  of the whole table, so it depends on the other samples and features.
  Compare SMDI values computed in the same run. The CLR centre itself
  cancels out: since the SMDI is a difference of two CLR means of the same
  sample, it equals `mean ln(p + ε)` over the dysbiotic taxa minus
  `mean ln(p + ε)` over the normobiotic taxa. The other taxa of the table
  matter only through ε.
* **Fixed pseudocount (optional).** `--p-smdi-pseudocount` replaces the
  table-dependent ε by a fixed value added to every relative abundance.
  The SMDI of a sample then depends on that sample only, and values from
  different runs can be compared. The values are not identical to the
  published definition, so the published cut-offs may not transfer
  exactly. `docs/VALIDATION.md` reports how much this matters. The default
  keeps the original definition.
* Signature taxa absent from the whole table enter as all-zero columns,
  that is, at the pseudocount. They would appear that way in a table
  merged with samples that contain them. The plugin warns and lists them.

**Built-in: `SMDI_simplified`.** The paper also reports a three-genus
index: *Fretibacterium* and *Treponema* (dysbiotic) against *Actinomyces*
(normobiotic). Its accuracy was close to the full index (AUC 0.937
training, 0.930 test), with 0 as the cut-off. It needs only genus-level
names, so it is the more robust choice when species resolution is poor.

**Built-in: `SMDI` (49 species).** The paper's Appendix Data Set 1 lists
the 201 discriminating species with their mean decrease in Gini (MDG). It
is published by the authors on the eHOMD server
(`homd.org/ftp/publication_data/20210324/Supplementary file 1.xlsx`). The
SMDI uses the species with MDG ≥ 0.4: exactly 49, of which 31 are
dysbiotic and 18 normobiotic. They are in `smdi_chen2022_species.tsv`,
with each original label and MDG value. Labels from the 2021 eHOMD were
updated to current names:

* `Bacteroidaceae [G-1] sp. oral taxon 272` → *Phocaeicola abscessus*
* *Propionibacterium acnes* → *Cutibacterium acnes*
* `Veillonella parvula_group` → *V. parvula*
* `[Eubacterium]` species keep their epithet

Unnamed phylotypes are matched by **HMT number**, so `Fretibacterium
sp. oral taxon 360` matches `Fretibacterium sp. HMT 360` whatever
placeholder genus a database version uses. eHOMD sub-species "clades"
(`Streptococcus infantis clade 431`) are merged into their species.

The authors report a cut-off around 0, with 0.4 the best in training and
0.5 the best in test. They also report that V4 data shift the index about
one unit lower than V1–V3.

**Custom signatures.** Any other list (`taxon`, `direction`, `level`) can
be imported as a `DysbiosisSignature` and passed with
`--i-smdi-signature`. It is computed with the same method and reported
as `SMDI_custom`.

### 2.2 PC2014_moderate

The fraction of reads from the 17 species and phylotypes with **moderate
evidence** of association with periodontitis in the systematic review of
Pérez-Chaparro et al. (2014, Table 3). The list is in
`perezchaparro2014_moderate.tsv`. Taxa are matched by name, or by Human
Microbial Taxon (HMT / HOT / "oral taxon") number for unnamed phylotypes
such as *Fretibacterium* sp. HMT 360 or Bacteroidales [G-2] HMT 274. The
review also places the phylum Saccharibacteria (TM7) and the domain
Archaea in this category. They are not species, so they are not counted.
This is a descriptive measure, not a diagnostic index.

### 2.3 Removed indices

Earlier versions reported `pMDI_ilr`, the Socransky Weighted Load (SWL)
and the Keystone Risk Score (KRS). The SWL weights and the KRS threshold
and sigmoid were not taken from any publication, and the keystone-pathogen
hypothesis (Hajishengallis et al. 2012) does not imply a threshold on Red
complex abundance. `pMDI_ilr` took geometric means over complexes that are
often zero, so its value was driven by the pseudocount. All three were
removed.

An intermediate version also had a **pMDI**: the log10 ratio of late
(Red + Orange) to early (Yellow + Green + Purple + Blue) colonizer reads.
The grouping comes from Socransky & Haffajee (2002), but the index itself
was never published or validated. In the validation data it was almost
redundant with the SMDI (Spearman ρ 0.86 with SMDI, 0.90 with
SMDI_simplified). It was removed, so that the only dysbiosis index is the
published one. The comparison of early and late colonizers is still
available from the complex table. Removing it does not change the
complexes in any way.

## 3. Statistics in `summarize`

* Tests run only on the metadata column you name (`--p-group-column`).
  Testing every metadata column against every index is a form of data
  dredging.
* With `--p-subject-column`, the samples of each subject within a group
  are averaged before testing. Several sites or visits from one subject
  are not independent observations, and treating them as such makes the
  p-values too small.
* Two groups: Mann-Whitney U, with AUC as the effect size. Three or more
  groups: Kruskal-Wallis, with epsilon-squared. q-values are
  Benjamini-Hochberg adjusted across the indices tested.
* Quality-control columns (`qc_*`) are never tested.

## 4. Deltas

`delta-indices` subtracts baseline from follow-up for every index (not for
`qc_*` columns) within the unit given by `--p-subject-col`. Use a subject
column for subject-level change. For site-level change, use a column that
combines subject and site. Several samples of one unit within a group are
averaged.
