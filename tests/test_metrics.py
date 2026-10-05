"""Behavioural checks for the revised evaluation contract."""
import sys
import unittest
from pathlib import Path

import numpy as np
from sklearn.metrics import average_precision_score

from kg_audit.metrics import ap_function, degree_bin, sample_candidates, assert_pair_contract, fit_degree_reference


class EvaluationTests(unittest.TestCase):
    def test_weighted_ap_with_ties_and_zero_weights(self):
        y = np.array([1, 0, 1, 0, 1, 0])
        score = np.array([.9, .9, .3, .2, .1, .1])
        weights = np.array([2., 0., 1., 3., 1., 2.])
        self.assertAlmostEqual(ap_function(y, score)(weights),
                               average_precision_score(y, score, sample_weight=weights))

    def test_zero_degree_has_its_own_bin(self):
        self.assertNotEqual(degree_bin(0), degree_bin(1))
        self.assertEqual(degree_bin(2), degree_bin(3))

    def test_cold_start_candidates_stay_in_heldout_drugs(self):
        positives = [("new", "a")]
        known = {("new", "a"), ("old", "b")}
        draws = sample_candidates(positives, ["new"], ["a", "b", "c"],
                                  known, {"old": 3}, {"a": 1, "b": 1, "c": 1},
                                  "matched", ratio=20, seed=8)
        self.assertEqual(len(draws), 20)
        self.assertTrue(all(c == "new" for c, d in draws))
        self.assertFalse(set(draws) & known)

    def test_insufficient_matched_support_fails_explicitly(self):
        with self.assertRaises(ValueError):
            sample_candidates([("c", "d")], ["c"], ["d"], {("c", "d")},
                              {}, {}, "matched", ratio=2, seed=1)

    def test_baseline_train_eval_overlap_is_rejected(self):
        with self.assertRaises(ValueError):
            assert_pair_contract([("a", "b")], [1], [("a", "b")])

    def test_duplicate_pair_conflicting_labels_is_rejected(self):
        with self.assertRaises(ValueError):
            assert_pair_contract([("a", "b"), ("a", "b")], [1, 0], [])

    def test_baseline_is_finite_and_only_training_labels_contribute_degree(self):
        train = [("a", "x"), ("b", "y"), ("a", "y")]
        neg = [("a", "z"), ("b", "z")]
        score, cdeg, ddeg = fit_degree_reference(train, neg)
        self.assertEqual(cdeg["new"], 0)
        self.assertEqual(ddeg["z"], 0)
        self.assertTrue(np.isfinite(score([("new", "x"), ("new", "z")])).all())

    def test_matching_preserves_joint_bins_and_repeated_draw_counts(self):
        from collections import Counter
        positives = [("a", "x"), ("b", "y")]
        cd = {"a": 1, "b": 3, "c": 1, "d": 2}
        dd = {"x": 2, "y": 0, "z": 3, "w": 0}
        draws = sample_candidates(positives,list(cd),list(dd),set(positives),cd,dd,'matched',20,12)
        bins = lambda pairs: Counter((degree_bin(cd[c]),degree_bin(dd[d])) for c,d in pairs)
        self.assertEqual(bins(draws),Counter({k:20*v for k,v in bins(positives).items()}))

    def test_uniform_is_reproducible_and_excludes_recorded_positives(self):
        args=([('a','x')],['a','b'],['x','y'],{('a','x')},{},{},'uniform')
        first=sample_candidates(*args,ratio=10,seed=8)
        self.assertEqual(first,sample_candidates(*args,ratio=10,seed=8))
        self.assertNotIn(('a','x'),first)

    def test_weighted_ap_equals_explicit_cluster_multiplicity(self):
        y=np.array([1,0,0,1,0])
        scores=np.array([.8,.8,.5,.4,.3])
        weights=np.array([2,2,0,3,1])
        expanded=np.repeat(np.arange(len(y)),weights)
        self.assertAlmostEqual(ap_function(y,scores)(weights),average_precision_score(y[expanded],scores[expanded]))


if __name__ == "__main__":
    unittest.main()
