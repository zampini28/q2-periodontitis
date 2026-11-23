# Check against the SMDI web calculator

The authors of the SMDI (Chen, Marsh & Al-Hebshi 2022) provide a web
calculator, <https://bioinformatics.forsyth.org/smdi/> (Calculator V1.2,
DS/DG List V1.0). It is the only implementation from the index's authors.
Its guide describes the pseudocount as the paper's text does,
"min(relative abundance)/2 to exact zero relative abundance entries",
whereas microbiome 1.6.0, which the paper's Methods name and the plugin
follows, adds that value to every entry when the table has a zero. This
check tells which convention the calculator uses, and whether it gives the
plugin's values.

The calculator asks for a typed verification code, so it is run by hand.

## Files

`make_tables.py` writes three species-level count tables in the
calculator's input format (`#SampleID` header row, one row per species,
raw counts, 8 samples, no metadata rows) and `expected.tsv`:

| Table | Content | What it tests |
|---|---|---|
| `dense.txt` | every species in every sample | the definition itself; no pseudocount under either convention, so both columns of `expected.tsv` are equal |
| `zeros.txt` | dysbiotic species at 0–2 reads, others at 50–499 with ~30% zeros | the pseudocount convention: the two columns differ by 0.1–0.25 |
| `absent.txt` | `zeros.txt` without 5 SMDI species | absent species entering as zeros (zero-fill) |

Rows are the 49 SMDI species, named with the labels of the paper's
Appendix Data Set 1 (`original_label` in
`q2_periodontitis/data/smdi_chen2022_species.tsv`, the naming of the
calculator's sample input), and 30 other oral species.

To regenerate (QIIME 2 environment):

```bash
python -I validation/calculator/make_tables.py validation/calculator
```

## Running the calculator

For each table:

1. Open <https://bioinformatics.forsyth.org/smdi/>.
2. Taxonomy level: **Species**. Zero-fill missing DS/DG: **Yes**. Count
   transformation: **None (raw count)**.
3. Upload the table, type the verification code, submit.
4. Check that the report lists no SMDI species as missing for `dense.txt`
   and `zeros.txt`, and exactly the 5 removed ones for `absent.txt`
   (*Treponema medium*, *Rothia aeria*, *Fretibacterium* sp. HMT 359,
   *Mycoplasma faucium*, *Prevotella loescheii*). If others are listed,
   the calculator did not recognise their names; note them.
5. Copy each sample's SMDI into the `calculator` column of `expected.tsv`,
   and save the result page (PDF or HTML) next to it as
   `calculator_<table>.pdf`, with the date.

## Reading the result

* `dense.txt` equal to `plugin_microbiome_1.6.0`: the calculator computes
  the index as defined (CLR over all rows, mean of dysbiotic minus mean of
  normobiotic species).
* `zeros.txt` and `absent.txt` equal to `plugin_microbiome_1.6.0`: the
  calculator follows the microbiome 1.6.0 behaviour, like the plugin.
* Equal to `zeros_only_pseudocount` instead: the calculator adds the
  pseudocount to zeros only. Its values then differ from the plugin's (and
  from microbiome 1.6.0's) by up to a few tenths on sparse tables; report
  this in the article.
* Equal to neither: compare the `dense.txt` results first; a difference
  there points to the species names or the definition, not the
  pseudocount.
