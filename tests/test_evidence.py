"""Declared evidence-handling invariants, independent of clinical classification."""

import unittest
from kg_audit.evidence import assess_strategy, validate_record, writeback_action, DOMAINS


def complete_record():
    return {
        "id": "synthetic:strategy1",
        "drug": "synthetic:drug1",
        "indication": "synthetic:disease1",
        "population": "synthetic population",
        "comparator": "synthetic comparator",
        "regimen": "synthetic regimen",
        "endpoint": "synthetic endpoint",
        "required_exposure": "specified for this synthetic fixture",
        "assessment_date": "2026-01-01",
        "evidence": [
            {"id": "synthetic:source1", "source": "https://example.org/synthetic", "source_date": "2025-12-01"}
        ],
        "constraints": [{"domain": d, "status": "pass", "evidence_ids": ["synthetic:source1"]} for d in DOMAINS],
        "failures": [],
    }


def scoped_failure(record, day, refs):
    return {
        "event_type": "stopped_trial",
        "trial_id": "NCT00000000",
        "registry_status": "TERMINATED",
        "stop_reason_categories": ["Negative"],
        "failure_type": "efficacy",
        "scope_match": "same_concept",
        "implication": "negate",
        "source_date": day,
        "evidence_ids": refs,
        "scope": {k: record[k] for k in ("drug", "indication", "population", "regimen", "endpoint", "comparator")},
    }


class EvidenceTests(unittest.TestCase):
    def test_complete_fixture_can_pass_documentation_gate(self):
        self.assertEqual(assess_strategy(complete_record(), "2026-01-01")["status"], "ready_for_evidence_review")

    def test_missing_every_domain_never_passes(self):
        for domain in DOMAINS:
            record = complete_record()
            record["constraints"] = [c for c in record["constraints"] if c["domain"] != domain]
            with self.subTest(domain=domain):
                self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "incomplete")

    def test_empty_or_unknown_critical_values_never_pass(self):
        for key in ["population", "comparator", "regimen", "endpoint", "required_exposure"]:
            for value in [None, "", "unknown", "not reported"]:
                record = complete_record()
                record[key] = value
                self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "incomplete")

    def test_future_evidence_cannot_support_current_readiness(self):
        record = complete_record()
        record["evidence"][0]["source_date"] = "2026-02-01"
        self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "incomplete")

    def test_future_failure_does_not_supersede_past_assessment(self):
        record = complete_record()
        record["failures"] = [scoped_failure(record, "2026-02-01", ["synthetic:source2"])]
        record["evidence"].append(
            {"id": "synthetic:source2", "source": "https://example.org/future", "source_date": "2026-02-01"}
        )
        self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "ready_for_evidence_review")
        self.assertEqual(assess_strategy(record, "2026-03-01")["status"], "scoped_counterevidence")

    def test_missing_timestamp_is_not_silently_current(self):
        record = complete_record()
        record["evidence"][0].pop("source_date")
        self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "incomplete")

    def test_unknown_exposure_with_explanatory_suffix_stays_unknown(self):
        record = complete_record()
        record["required_exposure"] = "unknown for the proposed effect"
        self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "incomplete")

    def test_future_assessment_is_not_applied_retroactively(self):
        record = complete_record()
        self.assertEqual(assess_strategy(record, "2025-12-15")["status"], "not_yet_assessed")

    def test_percentile_interval_may_exclude_observed_point(self):
        record = complete_record()
        record["evidence"][0]["effect"] = {
            "measure": "bootstrap difference",
            "estimate": 0.1,
            "lower": 0.11,
            "upper": 0.2,
        }
        self.assertEqual(validate_record(record), [])

    def test_orphan_evidence_is_validation_error(self):
        record = complete_record()
        record["constraints"][0]["evidence_ids"] = ["missing"]
        self.assertTrue(validate_record(record))

    def test_contradictory_domain_records_do_not_overwrite(self):
        record = complete_record()
        record["constraints"].append({"domain": DOMAINS[0], "status": "fail", "evidence_ids": ["synthetic:source1"]})
        self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "conflicting_records")

    def test_negate_without_scientific_type_or_same_concept_is_not_counterevidence(self):
        base = {"event_type": "stopped_trial", "trial_id": "NCT0", "registry_status": "TERMINATED",
                "implication": "negate", "source_date": "2025-01-01", "evidence_ids": ["x"]}
        self.assertEqual(writeback_action({**base, "failure_type": "efficacy", "scope_match": "same_concept"}),
                         "store_scoped_counterevidence")
        self.assertEqual(writeback_action({**base, "failure_type": "efficacy", "scope_match": "narrower_concept"}),
                         "store_context_qualification")
        self.assertEqual(writeback_action({**base, "failure_type": "operational", "scope_match": "same_concept"}),
                         "store_context_qualification")
        self.assertEqual(writeback_action({**base, "evidence_ids": [], "failure_type": "efficacy",
                                           "scope_match": "same_concept"}), "defer")

    def test_registry_status_alone_cannot_make_negative(self):
        for status in ["TERMINATED", "WITHDRAWN", "SUSPENDED", "UNKNOWN", "COMPLETED"]:
            self.assertEqual(writeback_action({"registry_status": status}), "defer")

    def test_implication_must_follow_rule(self):
        record = complete_record()
        bad = scoped_failure(record, "2025-12-01", ["synthetic:source1"])
        bad["failure_type"] = "operational"          # operational stop cannot be stated as 'negate'
        record["failures"] = [bad]
        self.assertTrue(any("write-back rule" in e for e in validate_record(record)))

    def test_recorded_indication_is_retained(self):
        record = complete_record()
        f = scoped_failure(record, "2025-12-01", ["synthetic:source1"])
        f["recorded_indication"] = True
        f["implication"] = "retain"
        record["failures"] = [f]
        self.assertEqual(validate_record(record), [])
        self.assertEqual(writeback_action(f), "retain_existing_evidence")

    def test_scope_comparison_ignores_case_and_spacing(self):
        record = complete_record()
        f = scoped_failure(record, "2025-12-01", ["synthetic:source1"])
        f["scope"]["population"] = "  SYNTHETIC   Population "
        record["failures"] = [f]
        self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "scoped_counterevidence")

    def test_different_regimen_or_comparator_does_not_supersede(self):
        for field in ("regimen", "comparator"):
            record = complete_record()
            f = scoped_failure(record, "2025-12-01", ["synthetic:source1"])
            f["scope"][field] = "a different context"
            record["failures"] = [f]
            self.assertEqual(assess_strategy(record, "2026-01-01")["status"], "ready_for_evidence_review")


if __name__ == "__main__":
    unittest.main()
