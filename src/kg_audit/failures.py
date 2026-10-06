"""Stop-reason typing, scope matching and the write-back rule for stopped trials.

Stop-reason categories follow the 17-class scheme of Razuvayevskaya et al. (Nat Genet 2024),
as distributed by Open Targets; Open Targets also emits 'Uncategorised'. Categories are
grouped into failure types. A stopped trial changes a drug-disease label only when its
failure type is scientific (efficacy or safety) and its condition is the same concept as
the graph disease; otherwise the evidence is kept as a qualification or deferred.
"""

STOP_STATUSES = ("TERMINATED", "WITHDRAWN", "SUSPENDED")

CATEGORY_TO_TYPE = {
    "Negative": "efficacy",
    "Safety_Sideeffects": "safety",
    "Business_Administrative": "operational",
    "Insufficient_Enrollment": "operational",
    "Logistics_Resources": "operational",
    "Study_Staff_Moved": "operational",
    "Covid19": "operational",
    "Another_Study": "operational",
    "Regulatory": "operational",
    "Study_Design": "design",
    "Interim_Analysis": "design",
    "Insufficient_Data": "uninformative",
    "Invalid_Reason": "uninformative",
    "No_Context": "uninformative",
    "Uncategorised": "uninformative",
    "Ethical_Reason": "operational",
    "Endpoint_Met": "success",
    "Success": "success",
}
FAILURE_TYPES = ("efficacy", "safety", "operational", "design", "uninformative", "success", "not_reported")
SCIENTIFIC = ("efficacy", "safety")
# precedence when one trial carries several categories: scientific evidence is never hidden
PRECEDENCE = ("efficacy", "safety", "success", "design", "operational", "uninformative", "not_reported")
SCOPE_MATCHES = ("same_concept", "narrower_concept", "broader_concept", "unmapped")
IMPLICATIONS = ("negate", "qualify", "retain", "defer")


def failure_type(categories):
    """Map one trial's stop-reason categories to a single failure type."""
    cats = [c for c in (categories or []) if c]
    if not cats:
        return "not_reported"
    unknown = [c for c in cats if c not in CATEGORY_TO_TYPE]
    if unknown:
        raise ValueError(f"Unknown stop-reason category: {unknown}")
    types = {CATEGORY_TO_TYPE[c] for c in cats}
    return next(t for t in PRECEDENCE if t in types)


def scope_match(depth_gap):
    """Scope of the trial condition relative to the graph disease.

    depth_gap is the number of ontology levels between the trial condition and the graph
    disease it was mapped to (0 = same concept, > 0 = trial condition is narrower).
    """
    if depth_gap is None:
        return "unmapped"
    if depth_gap == 0:
        return "same_concept"
    return "narrower_concept" if depth_gap > 0 else "broader_concept"


def implication(ftype, scope, recorded_indication=False):
    """Write-back implication for one stopped trial.

    retain  : the pair is already a recorded indication; a stopped trial does not remove it
    negate  : scientific failure (efficacy/safety) in the same disease concept
    qualify : scientific failure in a narrower or broader concept; stored as counterevidence
              for that context only
    defer   : operational, design, uninformative or unreported reasons; no label change
    """
    if ftype not in FAILURE_TYPES:
        raise ValueError(f"Unknown failure type: {ftype}")
    if scope not in SCOPE_MATCHES:
        raise ValueError(f"Unknown scope match: {scope}")
    if recorded_indication:
        return "retain"
    if ftype in SCIENTIFIC:
        return "negate" if scope == "same_concept" else "qualify"
    return "defer"


def pair_implication(trial_implications):
    """Combine trial-level implications for one drug-disease pair.

    One same-concept scientific failure is enough to negate the pair label; a recorded
    indication is always retained.
    """
    s = set(trial_implications)
    for x in ("retain", "negate", "qualify"):
        if x in s:
            return x
    return "defer"


def training_action(pair_impl):
    """Translate a pair implication into a training-label action."""
    return {"retain": "positive", "negate": "negative", "qualify": "mask", "defer": "mask"}[pair_impl]
