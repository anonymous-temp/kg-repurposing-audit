"""Small synthetic integration checks; no source database is downloaded."""

import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from kg_audit.writeback import run as writeback_run
from kg_audit.observational import run as observational_run, FULL_CONT, FULL_BIN


class OptionalWorkflowTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("torch"), "optional torch dependency not installed")
    def test_four_policy_workflow_with_synthetic_pairs(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            pairs = [(f"synthetic:c{i}", f"synthetic:d{j}") for i in range(10) for j in range(10)]
            positives = [(f"synthetic:c{i}", f"synthetic:d{i}") for i in range(10)] + [
                ("synthetic:c0", "synthetic:d1"),
                ("synthetic:c2", "synthetic:d3"),
            ]
            test = [("synthetic:c1", "synthetic:d2"), ("synthetic:c3", "synthetic:d4")]
            neg = [p for p in pairs if p not in positives + test]
            path = Path(tmp) / "splits.json"
            path.write_text(
                json.dumps(
                    {
                        "random": {
                            "train_pos": positives,
                            "train_neg": neg,
                            "test_pos": test,
                            "neg_random": neg,
                            "neg_degmatch": list(reversed(neg)),
                        }
                    }
                )
            )
            result = writeback_run(path, Path(tmp) / "out", seeds=[1, 2], epochs=2, threads=1)
            self.assertEqual(set(result.policy), {"no_writeback", "flat_negative", "uncertainty_mask", "typed"})
            zero = result[result.selected_fraction == 0]
            for _, cell in zero.groupby("sampler"):
                self.assertLess(cell.AP_mean.max() - cell.AP_mean.min(), 1e-12)
            self.assertTrue((Path(tmp) / "out/protocol.json").exists())

    @unittest.skipUnless(importlib.util.find_spec("duckdb"), "optional duckdb dependency not installed")
    def test_observational_workflow_exports_only_aggregates(self):
        import duckdb

        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            n = 80
            rng = np.random.default_rng(7)
            frame = pd.DataFrame(
                {
                    "subject_id": np.arange(n),
                    "stay_id": np.arange(n),
                    "t0": ["2020-01-01"] * n,
                    "sepsis3": [1] * n,
                    "race_grp": np.tile(["White", "Black", "Asian", "Other"], 20),
                    "treat": [0] * n,
                    "prevalent_user": [0] * n,
                    "alive_48h": [1] * n,
                    "event28": [0] * n,
                }
            )
            for col in FULL_CONT:
                frame[col] = rng.uniform(1, 4, n)
            for col in FULL_BIN:
                frame[col] = rng.integers(0, 2, n)
            frame["age"] = 30 + np.arange(n) % 40
            frame.to_csv(Path(tmp) / "cohort.csv", index=False)
            timing = pd.DataFrame(
                {
                    "stay_id": np.arange(n),
                    "t0": pd.to_datetime(["2020-01-01"] * n),
                    "pre_t0_statin": pd.to_datetime([None] * n),
                    "first_post_statin": pd.to_datetime(["2020-01-02" if i % 2 else "2020-01-05" for i in range(n)]),
                    "deathtime": pd.to_datetime(["2020-01-10" if i % 7 == 0 else None for i in range(n)]),
                    "dod": pd.to_datetime([None] * n),
                    "dischtime": pd.to_datetime(["2020-02-10"] * n),
                    "anchor_age": frame.age,
                    "anchor_year": [2020] * n,
                }
            )
            con = duckdb.connect(str(Path(tmp) / "data.duckdb"))
            con.register("fixture", timing)
            con.execute("CREATE TABLE cohort AS SELECT * FROM fixture")
            con.execute(
                "CREATE TABLE statin_exposure AS SELECT stay_id,t0,pre_t0_statin,first_post_statin FROM fixture"
            )
            con.close()
            out = Path(tmp) / "result.json"
            result = observational_run(
                Path(tmp) / "cohort.csv", None, out, bootstrap=5, database=Path(tmp) / "data.duckdb"
            )
            self.assertEqual(result["n"], n)
            self.assertIn("risk_ratio_CI95", result["models"]["demographic"])
            self.assertNotIn("subject_id", out.read_text())
            self.assertNotIn("stay_id", out.read_text())
