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


def test_signature_units_are_feature_unit_keys():
    synonyms = rcheck.load_synonyms()
    level, sig = rcheck.signature_units(rcheck.SMDI_SPECIES, synonyms)
    assert level == "species" and len(sig) == 49
    level, sig = rcheck.signature_units(rcheck.SMDI_SIMPLIFIED, synonyms)
    assert level == "genus" and set(sig["unit"]) == set(
        sig["unit"].str.lower())


@pytest.mark.skipif(not RSCRIPT, reason="RSCRIPT (microbiome 1.6.0) not set")
def test_plugin_matches_r_microbiome():
    check = rcheck.crosscheck(RSCRIPT, "", quick=True)
    assert len(check) == 4
    assert check["pass"].all(), check.to_string()


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
