#!/usr/bin/env python
"""Encode the four worked clinical examples as handoff records and run the checker on them.

Every value below comes from the cited primary source or registry record (accessed
6 October 2026). Writes examples/<name>.json and examples/checker_output.json.
"""
import json
from pathlib import Path

from kg_audit.evidence import assess_strategy, validate_record, writeback_action
from kg_audit.failures import failure_type, implication

HERE = Path(__file__).resolve().parent
AS_OF = "2026-10-06"


def unknown_constraints(known=None):
    known = known or {}
    domains = ["mechanism", "systemic_exposure", "tissue_exposure", "regimen", "population",
               "comparator_endpoint", "safety", "evidence_path"]
    return [{"domain": d, "status": known.get(d, ("unknown", []))[0], "evidence_ids": known.get(d, ("unknown", []))[1]}
            for d in domains]


def failure(event_type, cats, scope_m, recorded=False, **kw):
    ft = failure_type(cats)
    impl = implication(ft, scope_m, recorded_indication=recorded)
    return {"event_type": event_type, "stop_reason_categories": cats, "failure_type": ft, "scope_match": scope_m,
            "recorded_indication": recorded, "implication": impl, **kw}


examples = {}

# 1. Baricitinib in hospitalised COVID-19: a graph-derived hypothesis followed by randomised evidence.
examples["baricitinib_covid19"] = {
    "id": "example:baricitinib-covid19",
    "drug": "DRUGBANK:DB11817",
    "indication": "MONDO:0100096",
    "population": "adults hospitalised with COVID-19",
    "regimen": "baricitinib 4 mg orally once daily for up to 14 days or until discharge",
    "comparator": "usual care",
    "endpoint": "28-day all-cause mortality",
    "required_exposure": "unknown",
    "assessment_date": AS_OF,
    "knowledge_level": "knowledge_assertion",
    "agent_type": "manual_agent",
    "evidence": [
        {"id": "doi:10.1016/S0140-6736(20)30304-4", "source": "Richardson et al. Lancet 2020;395:e30-e31",
         "source_date": "2020-02-04", "predicate": "studied_to_treat",
         "proposition": "Knowledge-graph-assisted hypothesis that baricitinib could reduce viral entry and inflammation in COVID-19"},
        {"id": "doi:10.1056/NEJMoa2031994", "source": "ACTT-2, Kalil et al. N Engl J Med 2021;384:795-807",
         "source_date": "2020-12-11", "predicate": "treats", "research_phase": "clinical_trial_phase_3",
         "proposition": "Baricitinib plus remdesivir versus remdesivir alone shortened time to recovery",
         "effect": {"measure": "rate ratio for recovery", "estimate": 1.16, "lower": 1.01, "upper": 1.32}},
        {"id": "doi:10.1016/S2213-2600(21)00331-3", "source": "COV-BARRIER, Marconi et al. Lancet Respir Med 2021;9:1407-1418",
         "source_date": "2021-09-01", "predicate": "treats", "research_phase": "clinical_trial_phase_3",
         "proposition": "Composite primary endpoint not met; 28-day mortality lower with baricitinib plus standard of care",
         "effect": {"measure": "hazard ratio for 28-day mortality", "estimate": 0.57, "lower": 0.41, "upper": 0.78}},
        {"id": "doi:10.1016/S0140-6736(22)01109-6", "source": "RECOVERY Collaborative Group. Lancet 2022;400:359-368",
         "source_date": "2022-07-30", "predicate": "treats", "research_phase": "clinical_trial_phase_3",
         "proposition": "Allocation to baricitinib reduced 28-day mortality versus usual care",
         "effect": {"measure": "age-adjusted rate ratio for 28-day mortality", "estimate": 0.87, "lower": 0.77, "upper": 0.99}},
        {"id": "url:fda-olumiant-covid19-2022", "source": "US FDA approval of OLUMIANT for certain hospitalised adults with COVID-19 (announced 11 May 2022)",
         "source_date": "2022-05-11", "predicate": "treats", "clinical_approval_status": "fda_approved_for_condition"},
    ],
    "constraints": unknown_constraints({
        "mechanism": ("concern", ["doi:10.1016/S0140-6736(20)30304-4"]),
        "population": ("pass", ["doi:10.1016/S0140-6736(22)01109-6"]),
        "comparator_endpoint": ("pass", ["doi:10.1016/S0140-6736(22)01109-6"]),
        "regimen": ("pass", ["doi:10.1016/S0140-6736(22)01109-6"]),
        "evidence_path": ("pass", ["url:fda-olumiant-covid19-2022"]),
    }),
    "failures": [],
}

# 2. Pimozide in ALS: a registry record with unknown status and no results; nothing to write back.
examples["pimozide_als"] = {
    "id": "example:pimozide-als",
    "drug": "DRUGBANK:DB01100",
    "indication": "MONDO:0004976",
    "population": "unknown",
    "regimen": "unknown",
    "comparator": "unknown",
    "endpoint": "unknown",
    "required_exposure": "unknown",
    "assessment_date": AS_OF,
    "knowledge_level": "not_provided",
    "agent_type": "manual_agent",
    "evidence": [
        {"id": "NCT03272503", "source": "ClinicalTrials.gov NCT03272503 (status UNKNOWN, last update 2020-05-28, no posted results)",
         "source_date": "2020-05-28", "predicate": "in_clinical_trials_for", "research_phase": "clinical_trial_phase_2"},
    ],
    "constraints": unknown_constraints(),
    "failures": [],
}

# 3. Evacetrapib: an efficacy stop. At graph level the registered/mapped condition differs from the
#    graph disease, so the label is qualified; at review level the same trial is scoped counterevidence.
acc_scope = {"drug": "DRUGBANK:DB11655", "indication": "high-risk vascular disease",
             "population": "adults with high-risk vascular disease on standard therapy",
             "regimen": "evacetrapib 130 mg orally once daily", "comparator": "placebo",
             "endpoint": "first major adverse cardiovascular event"}
examples["evacetrapib_vascular"] = {
    "id": "example:evacetrapib-high-risk-vascular-disease",
    **acc_scope,
    "required_exposure": "plasma exposure reported in trial",
    "assessment_date": AS_OF,
    "knowledge_level": "knowledge_assertion",
    "agent_type": "manual_agent",
    "evidence": [
        {"id": "NCT01687998", "source": "ClinicalTrials.gov NCT01687998 (TERMINATED: 'Study termination due to insufficient efficacy.')",
         "source_date": "2016-07-31", "predicate": "in_clinical_trials_for", "research_phase": "clinical_trial_phase_3"},
        {"id": "doi:10.1056/nejmoa1609581", "source": "ACCELERATE, Lincoff et al. N Engl J Med 2017;376:1933-1942",
         "source_date": "2017-05-18", "predicate": "treats", "negated": True, "research_phase": "clinical_trial_phase_3",
         "effect": {"measure": "hazard ratio for primary composite endpoint", "estimate": 1.01, "lower": 0.93, "upper": 1.08}},
        {"id": "doi:10.1056/NEJMoa1706444", "source": "REVEAL (anacetrapib), N Engl J Med 2017;377:1217-1227",
         "source_date": "2017-09-28", "predicate": "treats", "proposition": "A different CETP inhibitor reduced major coronary events; not evidence about evacetrapib",
         "effect": {"measure": "rate ratio for major coronary events", "estimate": 0.91, "lower": 0.85, "upper": 0.97}},
    ],
    "constraints": unknown_constraints({"comparator_endpoint": ("fail", ["doi:10.1056/nejmoa1609581"]),
                                        "population": ("pass", ["doi:10.1056/nejmoa1609581"]),
                                        "regimen": ("pass", ["doi:10.1056/nejmoa1609581"])}),
    "failures": [
        failure("stopped_trial", ["Negative"], "same_concept", trial_id="NCT01687998", registry_status="TERMINATED",
                trial_phase="clinical_trial_phase_3", trial_start_date="2012-10-31",
                stop_reason_text="Study termination due to insufficient efficacy.",
                category_source="Open Targets Platform 26.09 stop-reason classifier",
                trial_condition="high-risk vascular disease", source_date="2017-05-18",
                evidence_ids=["NCT01687998", "doi:10.1056/nejmoa1609581"], scope=acc_scope),
    ],
}
# graph-level view of the same trial: Open Targets maps its condition to endpoint events
# (myocardial infarction, cerebrovascular accident, unstable angina). Matched on identifiers,
# the rule would treat (evacetrapib, myocardial infarction) as a same-concept efficacy failure
# and negate it, although myocardial infarction was an outcome event rather than the trial
# population. Identifier-level scope matching is necessary but not sufficient.
graph_level_evacetrapib = failure("stopped_trial", ["Negative"], "same_concept", trial_id="NCT01687998",
                                  registry_status="TERMINATED",
                                  trial_condition="MONDO:0005068 myocardial infarction (as mapped; an endpoint event)")

# 4. Plazomicin: approved for complicated UTI; the sponsor's later insolvency is an operational event.
examples["plazomicin_cuti"] = {
    "id": "example:plazomicin-cuti",
    "drug": "DRUGBANK:DB12615",
    "indication": "MONDO:0100338",
    "population": "adults with complicated urinary tract infection including pyelonephritis",
    "regimen": "plazomicin 15 mg/kg intravenously once daily",
    "comparator": "meropenem",
    "endpoint": "composite cure at day 5 and test-of-cure visit",
    "required_exposure": "therapeutic drug monitoring recommended in renal impairment",
    "assessment_date": AS_OF,
    "knowledge_level": "knowledge_assertion",
    "agent_type": "manual_agent",
    "evidence": [
        {"id": "doi:10.1056/nejmoa1801467", "source": "EPIC, Wagenlehner et al. N Engl J Med 2019;380:729-740",
         "source_date": "2019-02-21", "predicate": "treats", "research_phase": "clinical_trial_phase_3"},
        {"id": "url:fda-zemdri-2018", "source": "US FDA approval of ZEMDRI (plazomicin), 25 June 2018",
         "source_date": "2018-06-25", "predicate": "treats", "clinical_approval_status": "fda_approved_for_condition"},
        {"id": "url:sec-achaogen-8k-2019", "source": "Achaogen Inc. Form 8-K (Chapter 11 petition, 15 April 2019)",
         "source_date": "2019-04-15", "proposition": "Sponsor insolvency"},
    ],
    "constraints": unknown_constraints({"evidence_path": ("pass", ["url:fda-zemdri-2018"])}),
    "failures": [
        failure("corporate_event", ["Business_Administrative"], "same_concept", source_date="2019-04-15",
                evidence_ids=["url:sec-achaogen-8k-2019"], category_source="manual",
                stop_reason_text="Sponsor filed for Chapter 11 protection"),
    ],
}

report = {}
for name, rec in examples.items():
    (HERE / f"{name}.json").write_text(json.dumps(rec, indent=2) + "\n")
    report[name] = {"structural_errors": validate_record(rec),
                    "documentation_status": assess_strategy(rec, AS_OF),
                    "failure_actions": [writeback_action(f) for f in rec["failures"]]}
report["evacetrapib_graph_level_trial"] = {"failure_type": graph_level_evacetrapib["failure_type"],
                                           "scope_match": graph_level_evacetrapib["scope_match"],
                                           "implication": graph_level_evacetrapib["implication"]}
(HERE / "checker_output.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
