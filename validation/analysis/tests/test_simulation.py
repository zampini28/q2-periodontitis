"""Tests of the simulated known-truth validation."""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import simulation as sim  # noqa: E402


class TestSimulation(unittest.TestCase):

    def test_taxa_are_unique_and_complete(self):
        taxa, direction = sim.build_taxa()
        self.assertEqual(taxa["key"].nunique(), len(taxa))
        self.assertEqual(len(direction), 49)
        self.assertTrue(set(direction).issubset(set(taxa["key"])))
        self.assertTrue(set(taxa["complex"]).issubset(set(sim.COMPLEXES)))

    def test_communities_and_counts(self):
        taxa, direction = sim.build_taxa()
        relative, groups = sim.generate_communities(
            taxa, direction, n_health=4, n_periodontitis=4, seed=3)
        self.assertEqual(relative.shape, (8, len(taxa)))
        np.testing.assert_allclose(relative.sum(axis=1), 1.0)
        self.assertEqual(list(groups).count("health"), 4)
        self.assertEqual(list(groups).count("periodontitis"), 4)
        counts = sim.simulate_counts(relative, 1000, seed=3)
        self.assertEqual(counts.shape, relative.shape)
        np.testing.assert_array_equal(counts.sum(axis=1), 1000)

    def test_end_to_end(self):
        frame, detail = sim.run_simulation()
        self.assertEqual(len(frame), 6)
        self.assertEqual(set(frame["scenario"]), {"clean", "unassigned"})
        self.assertEqual(set(frame["depth"]), set(sim.DEPTHS))
        self.assertTrue((frame["smdi_pearson"] > 0.9).all())
        self.assertTrue((frame["assignment_accuracy"] > 0.8).all())
        self.assertIsNotNone(detail)
        self.assertEqual(len(detail["true_smdi"]), 100)


if __name__ == "__main__":
    unittest.main()
