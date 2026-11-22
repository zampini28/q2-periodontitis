import os
import re
import unittest
import warnings

import numpy as np
import pandas as pd
import qiime2
from qiime2.plugin.testing import TestPluginBase

from q2_periodontitis._methods import COMPLEXES
from q2_periodontitis._visualizer import (
    summarize, benjamini_hochberg, analysis_units, group_comparisons)


class TestBenjaminiHochberg(unittest.TestCase):

    def test_hand_computed(self):
        # p sorted: 0.01, 0.02, 0.03, 0.04 (n=4)
        # -> p*n/rank = 0.04, 0.04, 0.04, 0.04
        np.testing.assert_allclose(
            benjamini_hochberg([0.04, 0.01, 0.03, 0.02]),
            [0.04, 0.04, 0.04, 0.04])
        # 0.01*3/1=0.03, 0.04*3/2=0.06, 0.5*3/3=0.5
        np.testing.assert_allclose(benjamini_hochberg([0.5, 0.01, 0.04]),
                                   [0.5, 0.03, 0.06])
        # monotonicity: 0.03*2/1=0.06 > 0.04*2/2=0.04 -> both 0.04
        np.testing.assert_allclose(benjamini_hochberg([0.03, 0.04]),
                                   [0.04, 0.04])
        self.assertEqual(len(benjamini_hochberg([])), 0)


class TestSummarize(TestPluginBase):
    package = "q2_periodontitis.tests"

    def setUp(self):
        super().setUp()
        rng = np.random.default_rng(42)
        # One sample ID tries to close the <script> tag
        self.ids = ["S0</script><b>x"] + [f"S{i}" for i in range(1, 12)]
        n = len(self.ids)
        self.indices = pd.DataFrame(
            rng.random((n, 3)), index=self.ids,
            columns=["SMDI", "SMDI_simplified", "PC2014_moderate"])
        self.indices.iloc[6:, 0] += 10  # strong group effect on SMDI
        self.indices["qc_depth"] = 1000.0
        complexes = rng.random((n, len(COMPLEXES)))
        self.complexes = pd.DataFrame(
            complexes / complexes.sum(axis=1, keepdims=True),
            index=self.ids, columns=COMPLEXES)
        self.mapping = pd.DataFrame({"Complex": ["Red", "Other"],
                                     "Resolution": ["species", "unassigned"]},
                                    index=pd.Index(["F1", "F2"],
                                                   name="Feature ID"))
        self.metadata = qiime2.Metadata(pd.DataFrame({
            "status": ["health"] * 6 + ["periodontitis"] * 6,
            "subject": ["a", "a", "b", "b", "c", "c",
                        "d", "d", "e", "e", "f", "f"],
            "site": ["x", "y", "z"] * 4,
            "age": np.arange(n, dtype=float),
        }, index=pd.Index(self.ids, name="id")))

    def render(self, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            summarize(self.temp_dir.name, self.complexes, self.indices,
                      kwargs.pop("mapping", None), **kwargs)
        with open(os.path.join(self.temp_dir.name, "index.html")) as fh:
            return fh.read()

    def stats_table(self, html):
        match = re.search(r"Group comparisons(.*?)</table>", html, re.S)
        return match.group(1) if match else None

    def test_generates_index_html_without_metadata(self):
        html = self.render(mapping=self.mapping)
        self.assertIsNone(self.stats_table(html))
        for f in ["index.html", "css/style.css", "data/complexes.tsv",
                  "data/indices.tsv"]:
            self.assertTrue(
                os.path.exists(os.path.join(self.temp_dir.name, f)), f)
        self.assertIn("Reads per sample", html)

    def test_metadata_alone_runs_no_test(self):
        # No blanket testing of every metadata column
        self.assertIsNone(
            self.stats_table(self.render(metadata=self.metadata)))

    def test_two_group_comparison(self):
        html = self.render(metadata=self.metadata, group_column="status")
        table = self.stats_table(html)
        self.assertIn("Mann-Whitney U", table)
        self.assertIn("q-value", table)
        # Three indices tested, QC columns never
        self.assertEqual(table.count("<td><strong>"), 3)
        self.assertNotIn("qc_depth", table)
        # SMDI separates the groups perfectly: AUC = 1, significant after FDR
        self.assertIn("AUC (periodontitis &gt; health) = 1.000", table)
        self.assertIn("#27ae60", table)
        self.assertIn("Unit of analysis: <strong>samples</strong>", html)

    def test_subject_aggregation(self):
        units, metrics = analysis_units(self.indices, self.metadata,
                                        "status", "subject")
        self.assertEqual(len(units), 6)
        self.assertNotIn("qc_depth", metrics)
        expected = self.indices.loc[self.ids[:2], "SMDI"].mean()
        got = units.loc[units["__subject__"] == "a", "SMDI"].iloc[0]
        self.assertAlmostEqual(got, expected)
        html = self.render(metadata=self.metadata, group_column="status",
                           subject_column="subject")
        self.assertIn("n=3", self.stats_table(html))

    def test_auc_hand_computed(self):
        units = pd.DataFrame({"__group__": ["a", "a", "b", "b"],
                              "x": [1.0, 3.0, 2.0, 4.0]})
        result, = group_comparisons(units, ["x"])
        # pairs (b > a): 2>1, 4>1, 4>3 -> 3 of 4
        self.assertAlmostEqual(result["effect"], 0.75)
        self.assertEqual(result["effect_name"], "AUC (b > a)")

    def test_kruskal_for_more_groups(self):
        html = self.render(metadata=self.metadata, group_column="site")
        table = self.stats_table(html)
        self.assertIn("Kruskal-Wallis", table)
        self.assertIn("epsilon-squared", table)

    def test_invalid_columns(self):
        with self.assertRaisesRegex(ValueError, "categorical"):
            self.render(metadata=self.metadata, group_column="age")
        with self.assertRaisesRegex(ValueError, "not found"):
            self.render(metadata=self.metadata, group_column="nope")
        with self.assertRaisesRegex(ValueError, "require metadata"):
            self.render(group_column="status")
        with self.assertRaisesRegex(ValueError, "requires group_column"):
            self.render(metadata=self.metadata, subject_column="subject")

    def test_skipped_test_warns(self):
        self.indices["PC2014_moderate"] = 1.0  # identical values
        with self.assertWarnsRegex(UserWarning,
                                   "Skipping test of 'PC2014_moderate'"):
            summarize(self.temp_dir.name, self.complexes, self.indices, None,
                      self.metadata, "site")

    def test_nan_values_are_skipped(self):
        self.indices.iloc[0, 0] = np.nan
        html = self.render(metadata=self.metadata, group_column="status")
        self.assertIn("n=5", self.stats_table(html))
        self.assertNotIn("NaN,", html[html.rindex("<script>"):])

    def test_javascript_is_local(self):
        html = self.render(metadata=self.metadata, group_column="status")
        self.assertNotIn("cdn.jsdelivr.net", html)
        self.assertNotRegex(html, r'<script[^>]+src="https?://')
        for js in ["vega.min.js", "vega-lite.min.js", "vega-embed.min.js"]:
            self.assertIn(f'src="js/{js}"', html)
            self.assertTrue(os.path.getsize(
                os.path.join(self.temp_dir.name, "js", js)) > 0)

    def test_script_json_is_escaped(self):
        html = self.render()
        script = html[html.rindex("<script>"):]
        self.assertNotIn("S0</script>", script)
        self.assertIn("S0<\\/script>", script)
        self.assertIn("S0&lt;/script&gt;&lt;b&gt;x", html)
        self.assertEqual(html.count("</script>"), 4)

    def test_color_map_from_python(self):
        html = self.render()
        self.assertIn('"Red": "#E41A1C"', html)
        self.assertNotIn("Extended", html)


if __name__ == "__main__":
    unittest.main()
