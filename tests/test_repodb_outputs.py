"""Released repoDB audit outputs use DrugCentral indication rows only (contraindications excluded)."""
import json
from pathlib import Path

import pandas as pd

R = Path(__file__).resolve().parents[1] / "paper_results" / "repodb"


def test_drugcentral_indication_rows_only():
    audit = json.loads((R / "repodb_audit.json").read_text())
    assert audit["inputs"]["drugcentral_indication_rows"] == 12231
    assert audit["mapping"]["mapped_pairs"] == 1497
    assert audit["mapped_pairs"]["same"]["k"] == 148
    assert audit["mapped_pairs"]["same_or_broader"]["k"] == 237


def test_label_audit_matches_audit():
    lab = json.loads((R / "repodb_label_audit.json").read_text())
    assert lab["contamination"]["n"] == 1554
    assert lab["popularity_null"]["observed"] == lab["contamination"]["k"] == 153


def test_precision_sample_contains_only_current_matches():
    chk = pd.read_csv(R / "precision_check.tsv", sep="\t")
    pairs = pd.read_csv(R / "repodb_pairs.tsv", sep="\t")
    col = {"same concept": "same", "stage or subtype of an indication": "broader"}
    for r in chk.itertuples():
        m = pairs[(pairs.drug_id == r.drug_id) & (pairs.ind_name == r.ind_name)]
        assert m[col[r.level]].astype(bool).any(), r.check_id
