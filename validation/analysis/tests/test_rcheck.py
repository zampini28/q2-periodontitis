"""Tests of rcheck.py. The R cross-check runs only when RSCRIPT points to
the Rscript of the microbiome 1.6.0 environment (reproduce.sh sets it)."""
import os
import sys

import biom
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import rcheck  # noqa: E402

RSCRIPT = os.environ.get("RSCRIPT")


def test_compare_matches_nan_and_tolerance():
    plugin = pd.Series([1.0, np.nan, 2.0], index=["a", "b", "c"])
    assert rcheck.compare(plugin, plugin + 1e-12)["pass"]
    assert not rcheck.compare(plugin, plugin + 1e-6)["pass"]
    r = pd.Series([1.0, 0.5, 2.0], index=["a", "b", "c"])
    result = rcheck.compare(plugin, r)
    assert result["nan_mismatch"] == 1 and not result["pass"]


def test_child_peak_excludes_parent_memory():
    # ru_maxrss from os.wait4 would include this 300 MB block
    block = np.ones(300 * 1024 * 1024 // 8)
    seconds, peak = rcheck.timing(rcheck.run_child([
        sys.executable, "-I", rcheck.__file__, "--plugin-seconds", "100",
        "10"]))
    assert seconds > 0 and 0 < peak < rcheck.own_peak_mb() - 250
    del block


# Labels from q2_periodontitis/tests/test_methods.py (TestParser) plus edge
# cases: every feature must get the same unit from both parsers.
LABELS = [
    "d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; "
    "f__Porphyromonadaceae; g__Porphyromonas; s__Porphyromonas_gingivalis",
    "d__Bacteria; f__Prevotellaceae; g__Prevotella_7; "
    "s__Prevotella_melaninogenica",
    "d__Bacteria; f__Spirochaetaceae; g__Treponema_2; s__uncultured_bacterium",
    "d__Bacteria; g__Fusobacterium_A; s__Fusobacterium_A nucleatum_D",
    "d__Bacteria; g__Prevotella_D; s__Prevotella_D intermedia",
    "k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; "
    "f__Porphyromonadaceae; g__Porphyromonas; s__gingivalis",
    "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Porphyromonadaceae;"
    "Porphyromonas;gingivalis",
    "d__Bacteria; g__Fretibacterium; s__Fretibacterium_human_oral_taxon_361",
    "d__Bacteria; g__Porphyromonas; s__human_oral_metagenome",
    "d__Bacteria; g__Porphyromonas; s__Ambiguous_taxa",
    "d__Bacteria; g__uncultured; s__metagenome",
    "Unassigned",
    "",
    "d__Bacteria; g__Tannerella; s__Bacteroides_forsythus",
    "d__Bacteria; g__Wolinella; s__Wolinella_recta",
    "d__Bacteria; g__Peptostreptococcus; s__Peptostreptococcus_micros",
    "d__Bacteria; g__Actinomyces; s__odontolyticus",
    "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;"
    "f__Peptostreptococcaceae;g__Peptostreptococcaceae_[G-6];"
    "s__[Eubacterium]_nodatum",
    "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;"
    "f__Bacteroidales_[F-2];g__Bacteroidales_[G-2];s__bacterium_HMT_274",
    "k__Bacteria;g__Fusobacterium;s__nucleatum_subsp._polymorphum",
    "k__Bacteria;g__Fusobacterium;s__nucleatum_subsp._vincentii",
    "k__Bacteria;g__Fretibacterium;s__sp._HMT_360",
    "k__Bacteria;g__Desulfobulbus;s__sp._HMT_041",
    "d__Bacteria; g__Fretibacterium; s__Fretibacterium sp. oral taxon 359",
    "d__Bacteria; g__Streptococcus; s__Streptococcus_mitis/oralis",
    "d__Bacteria; g__Streptococcus; s__mitis/sanguis",
    "k__Bacteria;g__Streptococcus;s__infantis_clade_431",
    "k__Bacteria;g__Streptococcus;s__oralis_subsp._dentisani_clade_058",
    "k__Bacteria;g__Treponema;s__Treponema Denticola",
    "k__Bacteria;g__Treponema",
    "k__Bacteria;g__Treponema;s__",
    "k__Bacteria;g__Actinomyces;s__sp.",
    "k__Bacteria;g__Streptococcus;s__sp. str. 12",
]


@pytest.mark.skipif(not RSCRIPT, reason="RSCRIPT not set")
@pytest.mark.parametrize("level", ["species", "genus"])
def test_r_parser_gives_plugin_units(tmp_path, level):
    ids = ["F%d" % i for i in range(len(LABELS))]
    taxonomy = pd.DataFrame({"Taxon": LABELS}, index=pd.Index(ids))
    table = biom.Table(np.ones((len(ids), 1)), ids, ["S1"])
    path, out = tmp_path / "taxonomy.tsv", tmp_path / "units.tsv"
    rcheck.write_taxonomy(taxonomy, path)
    rcheck.run_child([RSCRIPT, rcheck.R_SCRIPT, "--units", str(path),
                      rcheck.SYNONYMS, level, str(out)])
    r_units = rcheck.read_units(out)
    plugin = rcheck.feature_units(table, taxonomy, rcheck.load_synonyms(),
                                  level)
    assert list(r_units.reindex(ids)) == plugin


@pytest.mark.skipif(not RSCRIPT, reason="RSCRIPT (microbiome 1.6.0) not set")
def test_plugin_matches_r_microbiome():
    check = rcheck.crosscheck(RSCRIPT, "", quick=True)
    assert len(check) == 3 * len(rcheck.SETTINGS) * len(rcheck.INDICES)
    assert check["pass"].all(), check.to_string()
    # both branches of microbiome's pseudocount rule are covered
    assert check["zeros"].any() and not check["zeros"].all()
    # the min_confidence setting must actually exclude features
    low = check[check["setting"] == "min_confidence=0.7"]
    assert (low["excluded"] > 0).all()


@pytest.mark.skipif(not RSCRIPT, reason="RSCRIPT (microbiome 1.6.0) not set")
def test_empty_sample_is_nan_in_both():
    _, table, taxonomy = next(rcheck.simulated_tables(quick=True))
    counts = table.matrix_data.toarray()[:, :5]
    counts = np.hstack([counts, np.zeros((counts.shape[0], 1))])
    empty = biom.Table(counts, table.ids(axis="observation"),
                       list(table.ids()[:5]) + ["empty"])
    rows = rcheck.check_table(RSCRIPT, "empty", "", empty, taxonomy,
                              rcheck.load_synonyms())
    assert all(r["pass"] for r in rows), rows
