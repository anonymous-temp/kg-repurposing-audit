"""Failure typing, write-back rule, schema consistency and the encoded clinical examples."""

import json
import unittest
from pathlib import Path

import jsonschema
import yaml

from kg_audit.evidence import assess_strategy, validate_record
from kg_audit.failures import (CATEGORY_TO_TYPE, FAILURE_TYPES, IMPLICATIONS, SCOPE_MATCHES, failure_type,
                               implication, pair_implication, scope_match, training_action)

ROOT = Path(__file__).resolve().parents[1]


class FailureRuleTests(unittest.TestCase):
    def test_failure_type_precedence_and_missing_reason(self):
        self.assertEqual(failure_type(["Business_Administrative", "Negative"]), "efficacy")
        self.assertEqual(failure_type(["Insufficient_Enrollment", "Safety_Sideeffects"]), "safety")
        self.assertEqual(failure_type(["Study_Design", "Covid19"]), "design")
        self.assertEqual(failure_type([]), "not_reported")
        self.assertEqual(failure_type(None), "not_reported")
        with self.assertRaises(ValueError):
            failure_type(["Made_Up"])

    def test_scope_match_from_depth_gap(self):
        self.assertEqual(scope_match(0), "same_concept")
        self.assertEqual(scope_match(3), "narrower_concept")
        self.assertEqual(scope_match(-1), "broader_concept")
        self.assertEqual(scope_match(None), "unmapped")

    def test_write_back_rule_table(self):
        for ft in FAILURE_TYPES:
            for sm in SCOPE_MATCHES:
                impl = implication(ft, sm)
                if ft in ("efficacy", "safety"):
                    self.assertEqual(impl, "negate" if sm == "same_concept" else "qualify")
                else:
                    self.assertEqual(impl, "defer")
                self.assertEqual(implication(ft, sm, recorded_indication=True), "retain")

    def test_pair_implication_order_and_training_action(self):
        self.assertEqual(pair_implication(["defer", "qualify", "negate"]), "negate")
        self.assertEqual(pair_implication(["negate", "retain"]), "retain")
        self.assertEqual(pair_implication(["defer"]), "defer")
        self.assertEqual([training_action(x) for x in IMPLICATIONS], ["negative", "mask", "positive", "mask"])


class SchemaConsistencyTests(unittest.TestCase):
    def setUp(self):
        self.schema = yaml.safe_load((ROOT / "schemas/handoff.linkml.yaml").read_text())

    def values(self, enum):
        return set(self.schema["enums"][enum]["permissible_values"])

    def test_enums_match_code(self):
        self.assertEqual(self.values("StopReasonCategory"), set(CATEGORY_TO_TYPE))
        self.assertEqual(self.values("FailureType"), set(FAILURE_TYPES))
        self.assertEqual(self.values("ScopeMatch"), set(SCOPE_MATCHES))
        self.assertEqual(self.values("Implication"), set(IMPLICATIONS))

    def test_json_schema_has_generated_enums(self):
        js = json.loads((ROOT / "schemas/handoff.schema.json").read_text())
        self.assertEqual(set(js["$defs"]["FailureType"]["enum"]), set(FAILURE_TYPES))
        self.assertIn("failure_type", js["$defs"]["FailureAnnotation"]["required"])


class ExampleTests(unittest.TestCase):
    def test_examples_validate_and_checker_output_is_reproducible(self):
        js = json.loads((ROOT / "schemas/handoff.schema.json").read_text())
        stored = json.loads((ROOT / "examples/checker_output.json").read_text())
        for name in ("baricitinib_covid19", "pimozide_als", "evacetrapib_vascular", "plazomicin_cuti"):
            rec = json.loads((ROOT / f"examples/{name}.json").read_text())
            jsonschema.validate(rec, js)
            self.assertEqual(validate_record(rec), [])
            self.assertEqual(assess_strategy(rec, "2026-10-06"), stored[name]["documentation_status"])
        self.assertEqual(stored["evacetrapib_vascular"]["documentation_status"]["status"], "scoped_counterevidence")
        self.assertEqual(stored["plazomicin_cuti"]["failure_actions"], ["defer"])


if __name__ == "__main__":
    unittest.main()
