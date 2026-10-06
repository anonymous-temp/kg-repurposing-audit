"""Conservative documentation rules, not clinical eligibility decisions.

A registry status alone never creates a negative indication. A stopped trial can support
scoped counterevidence only when its failure type is scientific (efficacy or safety) and its
condition is the same concept as the reviewed indication (kg_audit.failures.implication).
Dates constrain the information available at a requested review time.
"""

from datetime import date
import math
import re

from .failures import FAILURE_TYPES, IMPLICATIONS, SCOPE_MATCHES, implication as rule_implication

DOMAINS = (
    "mechanism",
    "systemic_exposure",
    "tissue_exposure",
    "regimen",
    "population",
    "comparator_endpoint",
    "safety",
    "evidence_path",
)
STATUSES = {"pass", "concern", "fail", "unknown"}
UNKNOWN = {"", "unknown", "not reported", "not available", "unspecified", "n/a", "none"}


def known(value):
    text = str(value).strip().lower()
    return (
        value is not None
        and text not in UNKNOWN
        and not text.startswith(("unknown ", "not reported ", "not available "))
    )


def iso_date(value):
    return date.fromisoformat(value) if isinstance(value, str) and value else None


def normalise(value):
    """Case- and whitespace-insensitive comparison key for scope fields."""
    return re.sub(r"\s+", " ", str(value).strip()).casefold() if value is not None else None


EVENT_TYPES = {"stopped_trial", "development_discontinuation", "market_withdrawal", "corporate_event"}


def validate_record(record):
    """Return structural errors. Missing scientific evidence is not fabricated."""
    errors = []
    if not isinstance(record, dict):
        return ["record must be an object"]
    for key in ("id", "drug", "indication"):
        if not known(record.get(key)):
            errors.append(f"missing {key}")
    for key in ("evidence", "constraints", "failures"):
        if not isinstance(record.get(key, []), list):
            errors.append(f"{key} must be a list")
    if errors:
        return errors
    evidence = record.get("evidence", [])
    if any(not isinstance(x, dict) for x in evidence):
        return ["evidence entries must be objects"]
    ids = [x.get("id") for x in evidence]
    if None in ids or len(set(ids)) != len(ids):
        errors.append("evidence identifiers must be present and unique")
    for entry in evidence:
        if not known(entry.get("source")):
            errors.append("evidence source missing")
        if entry.get("source_date"):
            try:
                iso_date(entry["source_date"])
            except ValueError:
                errors.append("invalid evidence source_date")
        effect = entry.get("effect")
        if effect is not None:
            if not isinstance(effect, dict):
                errors.append("effect must be an object")
                continue
            values = [effect.get(k) for k in ("estimate", "lower", "upper")]
            if any(not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v) for v in values):
                errors.append("effect and interval must be finite numbers")
            elif values[1] > values[2]:
                errors.append("effect interval bounds are reversed")
    for group in ("constraints", "failures"):
        for item in record.get(group, []):
            if not isinstance(item, dict):
                errors.append(f"{group} entry must be an object")
                continue
            refs = item.get("evidence_ids", [])
            if not isinstance(refs, list) or any(x not in ids for x in refs):
                errors.append("orphan evidence reference")
            if group == "constraints":
                if item.get("domain") not in DOMAINS:
                    errors.append("unknown constraint domain")
                if item.get("status") not in STATUSES:
                    errors.append("invalid constraint status")
            else:
                errors.extend(_failure_errors(item))
            if item.get("source_date"):
                try:
                    iso_date(item["source_date"])
                except ValueError:
                    errors.append("invalid failure source_date")
    return errors


def _failure_errors(item):
    """A failure annotation must state its type and scope, and its implication must follow the rule."""
    errs = []
    if item.get("event_type") not in EVENT_TYPES:
        errs.append("unknown failure event_type")
    if item.get("event_type") == "stopped_trial" and not (known(item.get("trial_id")) and item.get("registry_status")):
        errs.append("stopped trial needs trial_id and registry_status")
    ftype, scope, impl = item.get("failure_type"), item.get("scope_match"), item.get("implication")
    if ftype not in FAILURE_TYPES:
        errs.append("unknown failure_type")
    if scope not in SCOPE_MATCHES:
        errs.append("unknown scope_match")
    if impl not in IMPLICATIONS:
        errs.append("unknown implication")
    if not errs:
        expected = rule_implication(ftype, scope, recorded_indication=bool(item.get("recorded_indication")))
        if impl != expected:
            errs.append(f"implication '{impl}' does not follow the write-back rule (expected '{expected}')")
    return errs


def writeback_action(failure):
    """Return the evidence action for one failure annotation.

    Only a same-concept efficacy or safety failure with a dated source becomes scoped
    counterevidence; scientific failures in another concept are stored as qualifications;
    operational, design, uninformative and unreported reasons are deferred.
    """
    impl = failure.get("implication", "defer")
    refs = failure.get("evidence_ids")
    if impl == "retain":
        return "retain_existing_evidence" if refs else "defer"
    if impl == "qualify":
        return "store_context_qualification" if refs else "defer"
    if impl == "negate":
        if (
            failure.get("failure_type") in ("efficacy", "safety")
            and failure.get("scope_match") == "same_concept"
            and refs
            and failure.get("source_date")
        ):
            return "store_scoped_counterevidence"
        return "store_context_qualification" if refs else "defer"
    return "defer"


def assess_strategy(record, as_of):
    """Assess completeness and provenance as of an ISO date.

    A 'ready' result means a complete record for evidence review, not an actionable
    prescription, therapeutic efficacy, regulatory approval or expert consensus.
    """
    cutoff = iso_date(as_of)
    if cutoff is None:
        raise ValueError("as_of must be an ISO date")
    errors = validate_record(record)
    if errors:
        return {"status": "invalid_record", "reasons": errors, "as_of": as_of}
    if record.get("assessment_date"):
        try:
            assessment_date = iso_date(record["assessment_date"])
        except ValueError:
            return {"status": "invalid_record", "reasons": ["invalid assessment_date"], "as_of": as_of}
        if assessment_date > cutoff:
            return {
                "status": "not_yet_assessed",
                "reasons": ["assessment occurs after requested review date"],
                "as_of": as_of,
            }
    evidence = {x["id"]: x for x in record.get("evidence", [])}

    def available(refs):
        return bool(refs) and all(
            iso_date(evidence[i].get("source_date")) is not None and iso_date(evidence[i]["source_date"]) <= cutoff
            for i in refs
        )

    current = [x for x in record.get("constraints", []) if available(x.get("evidence_ids", []))]
    for domain in DOMAINS:
        values = {x["status"] for x in current if x["domain"] == domain}
        if len(values) > 1:
            return {"status": "conflicting_records", "reasons": [domain], "as_of": as_of}
    for failure in record.get("failures", []):
        day = iso_date(failure.get("source_date"))
        if day is not None and day <= cutoff and available(failure.get("evidence_ids", [])):
            action = writeback_action(failure)
            if action == "store_scoped_counterevidence":
                scope = failure.get("scope") or {}
                if all(
                    known(scope.get(k)) and normalise(scope.get(k)) == normalise(record.get(k))
                    for k in ("drug", "indication", "population", "regimen", "comparator", "endpoint")
                ):
                    return {
                        "status": "scoped_counterevidence",
                        "reasons": ["counterevidence retained with its regimen and population scope"],
                        "as_of": as_of,
                        "scope": scope,
                    }
    missing = [
        k for k in ("population", "comparator", "regimen", "endpoint", "required_exposure") if not known(record.get(k))
    ]
    observed = {x["domain"]: x["status"] for x in current}
    missing += [domain for domain in DOMAINS if domain not in observed or observed[domain] == "unknown"]
    if missing:
        return {"status": "incomplete", "reasons": sorted(set(missing)), "as_of": as_of}
    pending = [domain for domain, status in observed.items() if status != "pass"]
    if pending:
        return {"status": "requires_evidence_review", "reasons": pending, "as_of": as_of}
    return {"status": "ready_for_evidence_review", "reasons": [], "as_of": as_of}
