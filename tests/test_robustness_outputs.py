"""Consistency of the released robustness outputs (pipelines/26_robustness_checks.py, release 0.7.0)."""
import csv
import json
from pathlib import Path

import pytest

ROB = Path(__file__).resolve().parents[1] / "paper_results" / "robustness"
pytestmark = pytest.mark.skipif(not ROB.exists(), reason="robustness outputs not present")


def load(name):
    return json.loads((ROB / name).read_text())


def test_mantel_haenszel_matches_strata():
    m = load("mh_or_hetionet.json")
    num = den = 0.0
    n_pairs = 0
    for s in m["strata"]:
        a, b, c, d = s["stop_indication"], s["stop_other"], s["nostop_indication"], s["nostop_other"]
        n = a + b + c + d
        num += a * d / n
        den += b * c / n
        n_pairs += n
    assert n_pairs == m["pairs"]
    assert abs(num / den - m["mh_or"]) < 1e-9
    assert m["mh_ci"][0] < m["mh_or"] < m["mh_ci"][1]


def test_temporal_shares_and_rows():
    t = load("temporal_order_hetionet.json")
    for group in ("indications", "recorded", "external", "other_unlabelled"):
        for subset in ("all_stops", "scientific", "strict"):
            x = t[group][subset]
            a = x["after_first_approval"]
            assert a["n"] == 0 or abs(a["share"] - a["k"] / a["n"]) < 1e-4
            d = x["dated"]
            assert d["strictly_after"]["k"] <= d["in_or_after"]["k"] <= d["in_or_after"]["n"]
    with open(ROB / "temporal_rows_hetionet.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    after = sum(float(r["start_year"]) >= float(r["dated_approval"]) for r in rows)
    dated = t["indications"]["all_stops"]["dated"]["in_or_after"]
    assert len(rows) == dated["n"] and after == dated["k"]


def test_case_review_and_literature_sensitivity():
    c = load("case_review_sensitivity.json")
    assert c["final"]["yes"] <= c["any_coding_yes"]["cases"] <= c["any_coding_yes_or_unclear"]["cases"] <= 80
    assert len(c["any_coding_yes_or_unclear"]["case_ids"]) == c["any_coding_yes_or_unclear"]["cases"]
    lit = load("literature_second_screen.json")
    s = lit["changed_by_second_screen"]
    assert lit["first_decision_no"] == lit["sample"] == s["n"] and 0 <= s["k"] <= s["n"]


def test_approval_source_rate_ratios():
    a = load("approval_sources_hetionet.json")
    for key in ("all_non_trial_sources", "regulatory_labels", "us_labels"):
        x = a[key]
        assert abs(x["rate_ratio"] - x["share_stopped"] / x["share_other"]) < 1e-6
        assert x["ci"][0] < x["rate_ratio"] < x["ci"][1]
