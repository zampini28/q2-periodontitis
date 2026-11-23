"""Tests of parser_check.py."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import parser_check  # noqa: E402


def test_labels_are_every_rank_truncation_with_both_separators():
    out = parser_check.labels(["k__B;g__Porphyromonas;s__gingivalis"])
    assert out == sorted({
        "Unassigned", "k__B", "k__B;g__Porphyromonas",
        "k__B; g__Porphyromonas", "k__B;g__Porphyromonas;s__gingivalis",
        "k__B; g__Porphyromonas; s__gingivalis"})


def test_labels_strip_spaces_around_ranks():
    assert "d__B; g__T" in parser_check.labels(["d__B; g__T"])
    assert "d__B;g__T" in parser_check.labels(["d__B; g__T"])
