"""Tests of validation/analysis/stats.py.

Reference values come from metafor 4.8.0 and pROC 1.19.0.1 (R 4.3.3), see
reference_values.R; the others are computed by hand in the comments.

    python -m pytest validation/analysis/tests
"""
import json
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import stats as st  # noqa: E402

with open(os.path.join(HERE, "reference_values.json")) as fh:
    REF = json.load(fh)


class TestAuc(unittest.TestCase):

    def test_separation_and_ties(self):
        self.assertEqual(st.auc([3, 4], [1, 2]), 1.0)
        self.assertEqual(st.auc([1, 2], [3, 4]), 0.0)
        # one tie of four pairs: (1 + 1 + 0.5 + 1) / 4
        self.assertEqual(st.auc([2, 3], [1, 2]), 0.875)

    def test_hanley_mcneil(self):
        # A = 0.8, 50 vs 50: Q1 = 2/3, Q2 = 1.28/1.8,
        # var = (0.16 + 49 (Q1 - 0.64) + 49 (Q2 - 0.64)) / 2500
        q1, q2 = 0.8 / 1.2, 1.28 / 1.8
        var = (0.16 + 49 * (q1 - 0.64) + 49 * (q2 - 0.64)) / 2500
        self.assertAlmostEqual(st.hanley_mcneil_se(0.8, 50, 50),
                               np.sqrt(var), places=12)

    def test_bootstrap_is_reproducible_and_brackets(self):
        rng = np.random.default_rng(1)
        pos, neg = rng.normal(1, 1, 40), rng.normal(0, 1, 30)
        a = st.bootstrap_auc(pos, neg, 500, seed=3)
        self.assertEqual(a, st.bootstrap_auc(pos, neg, 500, seed=3))
        self.assertLess(a[0], st.auc(pos, neg))
        self.assertGreater(a[1], st.auc(pos, neg))


class TestDeLong(unittest.TestCase):

    def test_against_proc(self):
        aucs, cov = st.delong(np.vstack([REF["cases_a"], REF["cases_b"]]),
                              np.vstack([REF["ctrl_a"], REF["ctrl_b"]]))
        self.assertAlmostEqual(aucs[0], REF["auc_a"], places=10)
        self.assertAlmostEqual(aucs[1], REF["auc_b"], places=10)
        self.assertAlmostEqual(cov[0, 0], REF["var_a"], places=10)
        self.assertAlmostEqual(aucs[0], st.auc(REF["cases_a"], REF["ctrl_a"]))
        diff, se, z, p = st.delong_paired(REF["cases_a"], REF["ctrl_a"],
                                          REF["cases_b"], REF["ctrl_b"])
        self.assertAlmostEqual(z, REF["delong_z"], places=8)
        self.assertAlmostEqual(p, REF["delong_p"], places=8)


class TestMeta(unittest.TestCase):

    def test_reml_hksj_against_metafor(self):
        res = st.random_effects(REF["meta_y"], REF["meta_v"], "REML", True)
        self.assertAlmostEqual(res["tau2"], REF["reml_tau2"], places=7)
        self.assertAlmostEqual(res["estimate"], REF["reml_est"], places=7)
        self.assertAlmostEqual(res["se"], REF["reml_se"], places=7)
        self.assertAlmostEqual(res["ci_low"], REF["reml_ci"][0], places=6)
        self.assertAlmostEqual(res["ci_high"], REF["reml_ci"][1], places=6)

    def test_dl_against_metafor(self):
        res = st.random_effects(REF["meta_y"], REF["meta_v"], "DL", False)
        self.assertAlmostEqual(res["tau2"], REF["dl_tau2"], places=10)
        self.assertAlmostEqual(res["estimate"], REF["dl_est"], places=10)
        self.assertAlmostEqual(res["se"], REF["dl_se"], places=10)
        self.assertAlmostEqual(res["ci_low"], REF["dl_ci"][0], places=8)
        self.assertAlmostEqual(res["i2"], REF["dl_i2"], places=8)
        self.assertAlmostEqual(res["Q"], REF["Q"], places=8)

    def test_homogeneous_studies_give_zero_tau2(self):
        res = st.random_effects([0.8, 0.8, 0.8], [0.01, 0.02, 0.01])
        self.assertEqual(res["tau2"], 0.0)
        self.assertEqual(res["i2"], 0.0)
        self.assertAlmostEqual(res["estimate"], 0.8)


class TestMetaAuc(unittest.TestCase):

    def test_logit_pooling_is_back_transformed(self):
        a = np.array([0.90, 0.95, 0.99, 0.92])
        se = np.array([0.04, 0.02, 0.01, 0.03])
        res = st.random_effects_auc(a, se)
        y = np.log(a / (1 - a))
        ref = st.random_effects(y, (se / (a * (1 - a))) ** 2)
        inv = lambda x: 1 / (1 + np.exp(-x))  # noqa: E731
        self.assertAlmostEqual(res["estimate"], inv(ref["estimate"]))
        self.assertAlmostEqual(res["pi_high"], inv(ref["pi_high"]))
        self.assertAlmostEqual(res["tau2"], ref["tau2"])
        for key in ("ci_low", "ci_high", "pi_low", "pi_high"):
            self.assertTrue(0 < res[key] < 1)

    def test_identical_aucs_pool_to_themselves(self):
        res = st.random_effects_auc([0.8, 0.8, 0.8], [0.05, 0.04, 0.05])
        self.assertAlmostEqual(res["estimate"], 0.8)
        self.assertEqual(res["tau2"], 0.0)

    def test_auc_of_one_raises(self):
        with self.assertRaises(ValueError):
            st.random_effects_auc([1.0, 0.9, 0.9], [0.01, 0.05, 0.05])


class TestProportions(unittest.TestCase):

    def test_wilson(self):
        lo, hi = st.wilson(5, 10)
        self.assertAlmostEqual(lo, 0.2365931, places=6)
        self.assertAlmostEqual(hi, 0.7634069, places=6)

    def test_pooled_proportion_identical_studies(self):
        res = st.pooled_proportion([8, 8, 8], [10, 10, 10])
        self.assertAlmostEqual(res["estimate"], 0.8, places=10)
        # 100 % studies get the continuity correction instead of failing
        res = st.pooled_proportion([10, 9, 10], [10, 10, 10])
        self.assertTrue(0.8 < res["estimate"] < 1.0)


class TestPaired(unittest.TestCase):

    def test_rank_biserial(self):
        self.assertEqual(st.rank_biserial_paired([-1, -2, -3]), -1.0)
        # ranks 1, 2, 3; positive sum 3 (rank 3), negative sum 3
        self.assertEqual(st.rank_biserial_paired([-1, -2, 3, 0]), 0.0)

    def test_bh(self):
        q = st.bh([0.01, 0.04, 0.03, np.nan, 0.5])
        np.testing.assert_allclose(q[[0, 1, 2, 4]],
                                   [0.04, 0.16 / 3, 0.16 / 3, 0.5])
        self.assertTrue(np.isnan(q[3]))


if __name__ == "__main__":
    unittest.main()
