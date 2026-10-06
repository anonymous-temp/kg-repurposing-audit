import unittest
import numpy as np
from sklearn.metrics import average_precision_score
from kg_audit.comparison import paired_interval


class ComparisonTests(unittest.TestCase):
    def test_point_is_observed_paired_difference(self):
        y = np.array([1, 0, 1, 0, 0, 1, 0, 0])
        a = np.array([[0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2], [0.8, 0.9, 0.4, 0.3, 0.2, 0.7, 0.5, 0.6]])
        b = np.tile(np.arange(8), (2, 1))
        c = np.array(["a", "a", "b", "b", "c", "c", "d", "d"])
        d = np.tile(["x", "y"], 4)
        result = paired_interval(y, a, b, c, d, n_bootstrap=100)
        expected = np.mean([average_precision_score(y, x) - average_precision_score(y, z) for x, z in zip(a, b)])
        self.assertAlmostEqual(result["difference"], expected)

    def test_mismatched_shapes_rejected(self):
        with self.assertRaises(ValueError):
            paired_interval([1, 0], [0.8, 0.1], [0.8], ["a", "b"], ["x", "y"], n_bootstrap=100)

    def test_non_binary_labels_rejected(self):
        with self.assertRaises(ValueError):
            paired_interval([2, 0], [0.8, 0.1], [0.8, 0.1], ["a", "b"], ["x", "y"], n_bootstrap=100)
