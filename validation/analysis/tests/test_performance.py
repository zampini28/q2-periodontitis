"""Tests of the synthetic benchmark data generator."""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import performance as perf  # noqa: E402


class TestMakeTable(unittest.TestCase):

    def test_shape_and_nonzero_samples(self):
        table, taxonomy = perf.make_table(50, 6, density=0.1, seed=7)
        self.assertEqual(table.shape, (50, 6))
        self.assertEqual(list(taxonomy.index),
                         ["F%d" % i for i in range(50)])
        self.assertEqual(list(table.ids(axis="observation"))[:2],
                         ["F0", "F1"])
        self.assertEqual(list(table.ids(axis="sample"))[:2], ["S0", "S1"])
        self.assertTrue((np.asarray(table.sum(axis="sample")) > 0).all())

    def test_seed_is_reproducible(self):
        first, _ = perf.make_table(20, 4, density=0.2, seed=11)
        second, _ = perf.make_table(20, 4, density=0.2, seed=11)
        np.testing.assert_array_equal(
            first.matrix_data.toarray(), second.matrix_data.toarray())

    def test_positive_arguments(self):
        with self.assertRaises(ValueError):
            perf.make_table(0, 1)


if __name__ == "__main__":
    unittest.main()
