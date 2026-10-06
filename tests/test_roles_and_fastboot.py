"""Drug-role rules for stopped trials, role-restricted write-back policies and the inline bootstrap."""

import importlib.util
import unittest

import numpy as np

HAS_TORCH = importlib.util.find_spec("torch") is not None


class RoleRuleTests(unittest.TestCase):
    def setUp(self):
        from kg_audit import roles
        self.R = roles
        self.pat = {"doc": roles.name_patterns(["Docetaxel", "Taxotere"]), "gef": roles.name_patterns(["Gefitinib"]),
                    "pac": roles.name_patterns(["Paclitaxel"]), "mtx": roles.name_patterns(["Methotrexate"])}

    def role(self, c, study):
        return self.R.role_of(self.pat[c], study)[0]

    def test_comparator_backbone_investigational(self):
        vaccine = {"armGroups": [{"label": "1", "type": "EXPERIMENTAL", "interventionNames": ["Biological: vaccine"]},
                                 {"label": "2", "type": "OTHER", "interventionNames": ["Drug: Chemotherapy (Taxotere and prednisone)"]}],
                   "interventions": [{"type": "BIOLOGICAL", "name": "vaccine", "armGroupLabels": ["1"]},
                                     {"type": "DRUG", "name": "Chemotherapy (Taxotere and prednisone)", "armGroupLabels": ["2"]}]}
        self.assertEqual(self.role("doc", vaccine), "comparator")
        addon = {"armGroups": [{"label": "A", "type": "ACTIVE_COMPARATOR"}, {"label": "B", "type": "EXPERIMENTAL"}],
                 "interventions": [{"type": "DRUG", "name": "docetaxel", "armGroupLabels": ["A", "B"]},
                                   {"type": "DRUG", "name": "gefitinib", "armGroupLabels": ["B"]}]}
        self.assertEqual(self.role("doc", addon), "backbone")
        self.assertEqual(self.role("gef", addon), "investigational")

    def test_placebo_combination_arm_counts(self):
        study = {"armGroups": [{"label": "E", "type": "EXPERIMENTAL"}, {"label": "P", "type": "PLACEBO_COMPARATOR"}],
                 "interventions": [{"type": "DRUG", "name": "Sorafenib + Paclitaxel + Carboplatin", "armGroupLabels": ["E"]},
                                   {"type": "DRUG", "name": "Placebo + Paclitaxel + Carboplatin", "armGroupLabels": ["P"]}]}
        self.assertEqual(self.role("pac", study), "backbone")

    def test_described_only_is_not_investigational(self):
        study = {"armGroups": [{"label": "1", "type": "EXPERIMENTAL"}],
                 "interventions": [{"type": "DRUG", "name": "Rituximab", "description": "given with stable methotrexate",
                                    "armGroupLabels": ["1"]}]}
        self.assertEqual(self.role("mtx", study), "described_only")


@unittest.skipUnless(HAS_TORCH, "optional torch dependency not installed")
class RolePolicyTests(unittest.TestCase):
    def test_role_policies_negate_only_investigational_stops(self):
        from kg_audit import writeback as W
        ev = {"wb_all": {("c1", "d1"), ("c2", "d1"), ("c3", "d2")}, "wb_scientific": {("c1", "d1"), ("c2", "d1")},
              "wb_scientific_scoped": {("c1", "d1"), ("c2", "d1")}, "wb_other": {("c3", "d2")},
              "wb_scientific_inv": {("c1", "d1")}, "wb_scientific_inv_scoped": set()}
        role = W.policy_sets(ev, W.NAMED["typed_role"])
        self.assertEqual(role[("c1", "d1")], "negate"); self.assertEqual(role[("c2", "d1")], "mask"); self.assertEqual(role[("c3", "d2")], "mask")
        self.assertEqual(set(W.policy_sets(ev, W.NAMED["typed_role_scoped"]).values()), {"mask"})


@unittest.skipUnless(HAS_TORCH, "optional torch dependency not installed")
class FastBootTests(unittest.TestCase):
    def test_inline_bootstrap_matches_weighted_ap(self):
        from kg_audit import fastboot as F
        from kg_audit.writeback import weighted_ap
        rng = np.random.default_rng(0)
        n, nd = 5000, 30
        d = rng.integers(nd, size=n); y = (rng.random(n) < 0.05).astype(float)
        for s in (rng.random(n) + 0.3 * y, np.round(rng.random(n) * 10) + y):      # continuous and tied scores
            W = np.stack([np.bincount(rng.integers(nd, size=nd), minlength=nd).astype(float) for _ in range(4)])
            fast = F.pooled_ap_draws(y, s, d, W)
            slow = [weighted_ap(y, s, W[b][d]) for b in range(4)]
            np.testing.assert_allclose(fast, slow, rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
