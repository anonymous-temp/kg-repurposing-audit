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

    def test_remaining_role_branches(self):
        R = self.R
        pat = R.name_patterns(["Sorafenib"])
        role = lambda s: R.role_of(pat, s)[0]
        self.assertEqual(role({"interventions": [{"type": "DRUG", "name": "Imatinib"}]}), "unmatched")
        self.assertEqual(role({"interventions": [{"type": "DRUG", "name": "sorafenib"}]}), "investigational")              # no arms, one drug
        self.assertEqual(role({"interventions": [{"type": "DRUG", "name": "sorafenib"}, {"type": "DRUG", "name": "erlotinib"}]}), "no_arm_data")
        self.assertEqual(role({"armGroups": [{"label": "A", "type": "EXPERIMENTAL"}],
                               "interventions": [{"type": "DRUG", "name": "sorafenib", "armGroupLabels": ["A"]},
                                                 {"type": "DRUG", "name": "erlotinib", "armGroupLabels": ["A"]}]}), "single_arm_combination")
        self.assertEqual(role({"armGroups": [{"label": "A", "type": "ACTIVE_COMPARATOR"}, {"label": "B", "type": "ACTIVE_COMPARATOR"}],
                               "interventions": [{"type": "DRUG", "name": "sorafenib", "armGroupLabels": ["A"]},
                                                 {"type": "DRUG", "name": "sunitinib", "armGroupLabels": ["B"]}]}), "head_to_head")
        self.assertEqual(role({"armGroups": [{"label": "E", "type": "EXPERIMENTAL"}, {"label": "C1", "type": "ACTIVE_COMPARATOR"},
                                             {"label": "C2", "type": "PLACEBO_COMPARATOR"}],
                               "interventions": [{"type": "DRUG", "name": "sorafenib", "armGroupLabels": ["E", "C1"]}]}), "mixed")
        # intervention without arm labels: arms found through their intervention names
        self.assertEqual(role({"armGroups": [{"label": "E", "type": "EXPERIMENTAL", "interventionNames": ["Drug: sorafenib"]},
                                             {"label": "P", "type": "PLACEBO_COMPARATOR", "interventionNames": ["Drug: placebo"]}],
                               "interventions": [{"type": "DRUG", "name": "sorafenib"}]}), "investigational")
        self.assertEqual(role({"armGroups": [{"label": "E", "type": "EXPERIMENTAL"}, {"label": "P", "type": "PLACEBO_COMPARATOR"}],
                               "interventions": [{"type": "DRUG", "name": "sorafenib"}]}), "unmatched")
        self.assertFalse(R.usable({"name": "Placebo for sorafenib"}))
        self.assertTrue(R.usable({"name": "Placebo + sorafenib"}))
        self.assertEqual(R.name_patterns(["Na", "sodium", "12345", "Docetaxel"]), [" docetaxel "])

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

    def test_per_disease_draws(self):
        from kg_audit import fastboot as F
        from kg_audit.writeback import weighted_ap
        rng = np.random.default_rng(3)
        n, nd = 3000, 12
        d = rng.integers(nd, size=n); y = (rng.random(n) < 0.08).astype(float); s = rng.random(n) + 0.4 * y
        vals = F.per_disease_ap(y, s, d, nd, weighted_ap)
        for j in range(nd):
            m = d == j
            if 0 < y[m].sum() < m.sum():
                self.assertAlmostEqual(vals[j], weighted_ap(y[m], s[m]))
        W = np.ones((2, nd)); W[1, 0] = 3.0
        ok = ~np.isnan(vals)
        np.testing.assert_allclose(F.macro_draws(vals, W)[0], vals[ok].mean())
        np.testing.assert_allclose(F.macro_draws(vals, W)[1], (vals[ok] * W[1, ok]).sum() / W[1, ok].sum())


if __name__ == "__main__":
    unittest.main()
