"""Export ranked candidates as evidence-review records with their clinical evidence attached.

A score never fills a clinical field. Clinical fields are filled only from dated, identified
sources (here, Open Targets clinical reports): approval-stage records become evidence
assertions, trials become 'in clinical trials for' assertions, and stopped trials become
failure annotations whose implication follows kg_audit.failures.implication.
"""

import hashlib
import json
from collections import Counter
from pathlib import Path

from .evidence import DOMAINS, assess_strategy, writeback_action
from .failures import STOP_STATUSES, failure_type, implication, pair_implication, scope_match

PHASE = {"PHASE1": "clinical_trial_phase_1", "PHASE1/PHASE2": "clinical_trial_phase_1_to_2", "PHASE2": "clinical_trial_phase_2",
         "PHASE2/PHASE3": "clinical_trial_phase_2_to_3", "PHASE3": "clinical_trial_phase_3", "PHASE4": "clinical_trial_phase_4",
         "EARLY_PHASE1": "clinical_trial_phase_1"}


def prediction_record(compound, disease, score, rank, candidate_count, provenance, as_of):
    """A candidate record with prediction provenance and every clinical domain unknown."""
    key = json.dumps([provenance["graph"], provenance["task"], provenance["write_back_policy"], provenance["seed"],
                      compound, disease], separators=(",", ":"))
    return {
        "id": "prediction:" + hashlib.sha256(key.encode()).hexdigest()[:24],
        "drug": "DRUGBANK:" + compound,
        "indication": disease,
        "assessment_date": as_of,
        "knowledge_level": "prediction",
        "agent_type": "computational_model",
        "prediction": {**provenance, "score": float(score), "rank": int(rank), "candidate_count": int(candidate_count),
                       "score_is_probability": False},
        "evidence": [],
        "constraints": [{"domain": d, "status": "unknown", "evidence_ids": []} for d in DOMAINS],
        "failures": [],
    }


def _date(value):
    if value is None or value != value:      # None or NaN
        return None
    text = str(value)[:10]
    return text if len(text) == 10 else None


def attach_clinical_evidence(record, reports, recorded_indication=False):
    """Attach Open Targets clinical reports (rows for this drug-disease pair) to a record.

    reports: iterable of dicts with keys report_id, origin, stage, source, phase, status,
    stop_categories ('|'-joined), why_stopped, start_date, url, gap (ontology depth gap).
    """
    seen = set()
    trial_impl = []
    for r in reports:
        rid = str(r["report_id"])
        if rid in seen:
            continue
        seen.add(rid)
        if r["origin"] != "CLINICAL_TRIAL":
            if r["stage"] == "APPROVAL":
                record["evidence"].append({"id": rid, "source": f"{r['source']} via Open Targets", "predicate": "treats",
                                           "clinical_approval_status": "approved_for_condition"})
            continue
        ev = {"id": rid, "source": r.get("url") or rid, "predicate": "in_clinical_trials_for",
              "research_phase": PHASE.get(str(r.get("phase")), "not_provided")}
        if _date(r.get("start_date")):
            ev["source_date"] = _date(r.get("start_date"))
        record["evidence"].append(ev)
        if r.get("status") in STOP_STATUSES:
            cats = [c for c in str(r.get("stop_categories") or "").split("|") if c and c != "nan"]
            ft = failure_type(cats)
            gap = r.get("gap")
            sm = scope_match(None if gap is None or gap < 0 else int(gap))
            impl = implication(ft, sm, recorded_indication=recorded_indication)
            trial_impl.append(impl)
            fa = {"event_type": "stopped_trial", "trial_id": rid.upper(), "registry_status": r["status"],
                  "stop_reason_categories": cats, "category_source": "Open Targets Platform 26.09 stop-reason classifier",
                  "failure_type": ft, "scope_match": sm, "recorded_indication": recorded_indication,
                  "implication": impl, "evidence_ids": [rid]}
            if r.get("why_stopped") and r["why_stopped"] == r["why_stopped"]:
                fa["stop_reason_text"] = str(r["why_stopped"])
            if _date(r.get("start_date")):
                fa["trial_start_date"] = _date(r.get("start_date"))
                fa["source_date"] = _date(r.get("start_date"))
            if ev["research_phase"] != "not_provided":
                fa["trial_phase"] = ev["research_phase"]
            record["failures"].append(fa)
    return pair_implication(trial_impl) if trial_impl else None


def export_ranked(candidates, reports_by_pair, recorded, output, provenance, as_of):
    """Write ranked candidates with attached evidence; return a summary of what reviewers would see.

    candidates: list of (compound, disease, score, rank) sorted by rank
    """
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    status, pair_impl, actions, kinds = Counter(), Counter(), Counter(), Counter()
    n = len(candidates)
    with out.open("w") as fh:
        for compound, disease, score, rank in candidates:
            rec = prediction_record(compound, disease, score, rank, provenance["universe_size"],
                                    {k: v for k, v in provenance.items() if k != "universe_size"}, as_of)
            reps = reports_by_pair.get((compound, disease), [])
            impl = attach_clinical_evidence(rec, reps, recorded_indication=(compound, disease) in recorded)
            pair_impl[impl or "no_stopped_trial"] += 1
            has_appr = any(e.get("clinical_approval_status") for e in rec["evidence"])
            has_trial = any(e.get("predicate") == "in_clinical_trials_for" for e in rec["evidence"])
            kinds["approved_elsewhere" if has_appr else ("trialled_not_approved" if has_trial else "no_clinical_record")] += 1
            for f in rec["failures"]:
                actions[writeback_action(f)] += 1
            status[assess_strategy(rec, as_of)["status"]] += 1
            fh.write(json.dumps(rec, separators=(",", ":")) + "\n")
    summary = {"records": n, "clinical_record_kind": dict(kinds), "pair_implication": dict(pair_impl),
               "failure_actions": dict(actions), "documentation_status": dict(status),
               "biomedical_fields_inferred_from_scores": 0}
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2))
    return summary
