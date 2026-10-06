import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import pandas as pd
import jsonschema
from kg_audit import cli
from kg_audit.bridge import attach_clinical_evidence, export_ranked, prediction_record


class IOTests(unittest.TestCase):
    def test_cli_demo_and_paired_comparison(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            cli.demo(tmp)
            frame = pd.read_csv(Path(tmp) / "synthetic_scores.csv")
            self.assertFalse(frame.duplicated(["compound", "disease"]).any())
            output = Path(tmp) / "comparison.json"
            with patch(
                "sys.argv",
                [
                    "kg-audit",
                    "compare",
                    "--scores",
                    str(Path(tmp) / "synthetic_scores.csv"),
                    "--model",
                    "model",
                    "--reference",
                    "reference",
                    "--bootstrap",
                    "100",
                    "--output",
                    str(output),
                ],
            ):
                cli.main()
            self.assertIn("input_sha256", json.loads(output.read_text()))

    def test_cli_evidence_preserves_missingness(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            record = Path(tmp) / "record.json"
            record.write_text(json.dumps({"id": "synthetic:1", "drug": "synthetic:d", "indication": "synthetic:i"}))
            output = Path(tmp) / "result.json"
            with patch(
                "sys.argv",
                ["kg-audit", "evidence", "--record", str(record), "--as-of", "2026-01-01", "--output", str(output)],
            ):
                cli.main()
            self.assertEqual(json.loads(output.read_text())["status"], "incomplete")

    def test_bridge_attaches_typed_trial_evidence_without_inferring_clinical_fields(self):
        prov = {"model": "synthetic", "graph": "synthetic", "task": "all", "candidate_universe": "grid",
                "write_back_policy": "no_writeback", "seed": 1, "rank_scope": "global", "source_sha256": "0" * 64}
        reports = [
            {"report_id": "nct00000001", "origin": "CLINICAL_TRIAL", "stage": "PHASE_3", "source": "ClinicalTrials.gov",
             "phase": "PHASE3", "status": "TERMINATED", "stop_categories": "Negative", "why_stopped": "Futility",
             "start_date": "2012-01-01", "url": "https://clinicaltrials.gov/study/NCT00000001", "gap": 0},
            {"report_id": "nct00000002", "origin": "CLINICAL_TRIAL", "stage": "PHASE_2", "source": "ClinicalTrials.gov",
             "phase": "PHASE2", "status": "WITHDRAWN", "stop_categories": "Business_Administrative", "why_stopped": "Funding",
             "start_date": "2013-01-01", "url": "https://clinicaltrials.gov/study/NCT00000002", "gap": 0},
            {"report_id": "label1", "origin": "DRUG_LABEL", "stage": "APPROVAL", "source": "DailyMed", "phase": None,
             "status": None, "stop_categories": None, "why_stopped": None, "start_date": None, "url": None, "gap": 0},
        ]
        rec = prediction_record("DB00001", "DOID:1", 0.5, 3, 10, prov, "2026-01-01")
        impl = attach_clinical_evidence(rec, reports)
        self.assertEqual(impl, "negate")
        self.assertEqual([f["implication"] for f in rec["failures"]], ["negate", "defer"])
        self.assertTrue(any(e.get("clinical_approval_status") for e in rec["evidence"]))
        self.assertTrue(all(k not in rec for k in ("regimen", "population", "comparator", "endpoint")))
        rec2 = prediction_record("DB00001", "DOID:1", 0.5, 3, 10, prov, "2026-01-01")
        self.assertEqual(attach_clinical_evidence(rec2, reports, recorded_indication=True), "retain")
        narrower = [dict(reports[0], gap=2)]
        rec3 = prediction_record("DB00001", "DOID:1", 0.5, 3, 10, prov, "2026-01-01")
        self.assertEqual(attach_clinical_evidence(rec3, narrower), "qualify")
        schema = json.loads((Path(__file__).resolve().parents[1] / "schemas/handoff.schema.json").read_text())
        for r in (rec, rec2, rec3):
            jsonschema.validate(r, schema)
        with tempfile.TemporaryDirectory() as tmp:
            summary = export_ranked([("DB00001", "DOID:1", 0.5, 1)], {("DB00001", "DOID:1"): reports}, set(),
                                    Path(tmp) / "r.jsonl", {**prov, "universe_size": 10}, "2026-01-01")
            self.assertEqual(summary["pair_implication"], {"negate": 1})
            self.assertEqual(summary["biomedical_fields_inferred_from_scores"], 0)
