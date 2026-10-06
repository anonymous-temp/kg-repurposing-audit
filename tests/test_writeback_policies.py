"""Write-back policies act on training labels as specified; metrics match their definitions."""

import importlib.util
import unittest

import numpy as np
from sklearn.metrics import average_precision_score

HAS_TORCH = importlib.util.find_spec("torch") is not None


@unittest.skipUnless(HAS_TORCH, "optional torch dependency not installed")
class PolicyTests(unittest.TestCase):
    def setUp(self):
        from kg_audit import writeback as W
        self.W = W
        self.ev = {"wb_all": {("c1", "d1"), ("c2", "d1"), ("c3", "d2"), ("c1", "d2")},
                   "wb_scientific": {("c1", "d1"), ("c3", "d2")}, "wb_scientific_scoped": {("c1", "d1")},
                   "wb_other": {("c2", "d1"), ("c1", "d2")}}
        self.comps, self.dises = ["c1", "c2", "c3", "c4"], ["d1", "d2", "d3"]

    def test_named_policies(self):
        P = lambda n: self.W.policy_sets(self.ev, self.W.NAMED[n])
        self.assertEqual(set(P("no_writeback").values()), {"ignore"})
        self.assertEqual(set(P("flat_negative").values()), {"negate"})
        self.assertEqual(set(P("mask_all").values()), {"mask"})
        typed = P("typed")
        self.assertEqual(typed[("c1", "d1")], "negate"); self.assertEqual(typed[("c2", "d1")], "mask")
        scoped = P("typed_scoped")
        self.assertEqual(scoped[("c1", "d1")], "negate"); self.assertEqual(scoped[("c3", "d2")], "mask")

    def test_training_matrix_weights(self):
        acts = self.W.policy_sets(self.ev, self.W.NAMED["flat_negative"])
        fit = [("c1", "d1"), ("c4", "d3")]                      # (c1, d1) is both recorded and stopped
        Y, Wt = self.W.training_matrix(self.comps, self.dises, fit, {("c4", "d1")}, acts, w_neg=10.0, masked_rows={"c3"})
        self.assertEqual(Y[0, 0], 1.0)                            # recorded positive takes precedence
        self.assertEqual(Wt[1, 0], 10.0)                          # explicit negative
        self.assertEqual(Wt[3, 0], 0.0)                           # CpD excluded
        self.assertEqual(Wt[2, 1], 10.0)                          # write-back negative on a held-out compound row
        self.assertEqual(Wt[2, 2], 0.0)                           # otherwise held-out rows carry no weight
        self.assertAlmostEqual(Wt[Y == 1].sum(), Wt[Y == 0].sum(), places=4)
        acts = self.W.policy_sets(self.ev, self.W.NAMED["mask_all"])
        Y, Wt = self.W.training_matrix(self.comps, self.dises, fit, set(), acts)
        self.assertEqual(Wt[1, 0], 0.0)

    def test_scorer_runs_and_returns_requested_epochs(self):
        rng = np.random.default_rng(0)
        Y = (rng.random((6, 5)) < 0.3).astype(np.float32); Y[0, 0] = 1
        Wt = np.ones_like(Y); Wt[Y == 1] = (Y == 0).sum() / (Y == 1).sum()
        out = self.W.fit_scorer(Y, Wt, mf_dim=2, epochs=(5, 10))
        self.assertEqual(sorted(out), [5, 10]); self.assertEqual(out[10].shape, Y.shape)


class MetricTests(unittest.TestCase):
    @unittest.skipUnless(HAS_TORCH, "optional torch dependency not installed")
    def test_weighted_ap_matches_sklearn(self):
        from kg_audit.writeback import weighted_ap
        rng = np.random.default_rng(1)
        y = (rng.random(200) < 0.1).astype(float); y[0] = 1
        s = rng.normal(size=200)
        self.assertAlmostEqual(weighted_ap(y, s), average_precision_score(y, s), places=10)
        s2 = np.round(s, 1)                                       # ties
        self.assertAlmostEqual(weighted_ap(y, s2), average_precision_score(y, s2), places=10)


if __name__ == "__main__":
    unittest.main()
