"""Role of a compound in a clinical trial, from ClinicalTrials.gov arm groups and interventions.

Roles: investigational | comparator | backbone | mixed | head_to_head | single_arm_combination |
described_only | unmatched | no_arm_data (Supplementary Section S16 of the manuscript).
"""
import re

CONTROL = {"ACTIVE_COMPARATOR", "PLACEBO_COMPARATOR", "SHAM_COMPARATOR", "NO_INTERVENTION"}
DRUG_TYPES = ("DRUG", "BIOLOGICAL", "COMBINATION_PRODUCT", "GENETIC")
GENERIC = {"acid", "sodium", "potassium", "calcium", "chloride", "hydrochloride", "vitamin", "oral", "tablet", "tablets",
           "injection", "standard", "control", "placebo", "therapy", "treatment", "drug", "dose", "solution", "cream",
           "water", "saline", "oxygen", "iron", "zinc", "magnesium", "glucose", "alcohol", "ethanol", "sugar"}
COMBO = re.compile(r"\+|\bplus\b|\band\b|\bwith\b|,|;|/")


def norm(s):
    return " " + re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip() + " "


def name_patterns(names):
    """Whole-phrase patterns for a compound's names (at least four characters, generic words removed)."""
    toks = sorted({norm(n).strip() for n in names if n} - GENERIC, key=len, reverse=True)
    return [" " + t + " " for t in toks if len(t) >= 4 and not t.isdigit()]


def mentions(pat, text):
    return any(t in text for t in pat)


def usable(intervention):
    """Placebo interventions count only when they combine placebo with active drugs ('Placebo + Paclitaxel + Carboplatin')."""
    name = intervention.get("name", "").lower()
    return "placebo" not in name or bool(COMBO.search(name))


def role_of(pat, study):
    """Return (role, number of arms) for a compound with name patterns `pat` in a study's armsInterventionsModule."""
    arms = study.get("armGroups", []); ivs = study.get("interventions", [])
    drug_ivs = [i for i in ivs if i.get("type") in DRUG_TYPES and "placebo" not in i.get("name", "").lower()]
    hit = [i for i in ivs if usable(i) and mentions(pat, norm(i.get("name", "") + " " + " ".join(i.get("otherNames", []))))]
    if not hit:
        described = any(mentions(pat, norm(i.get("description", ""))) for i in ivs) or \
            any(mentions(pat, norm(g.get("description", ""))) for g in arms)
        return ("described_only" if described else "unmatched"), len(arms)
    if not arms:
        return ("investigational" if len(drug_ivs) <= 1 else "no_arm_data"), 0
    labels = {g.get("label") for g in arms}
    c_arms = {l for i in hit for l in i.get("armGroupLabels", [])} & labels
    if not c_arms:
        c_arms = {g.get("label") for g in arms if mentions(pat, norm(" ".join(g.get("interventionNames", []))))}
    if not c_arms:
        return "unmatched", len(arms)
    if len(arms) == 1:
        n_drug = len({i.get("name") for i in drug_ivs})
        return ("investigational" if n_drug <= 1 else "single_arm_combination"), 1
    types = {g.get("label"): g.get("type", "NA") for g in arms}
    has_exp = any(t == "EXPERIMENTAL" for t in types.values())
    if has_exp:
        ctrl = {l for l, t in types.items() if t != "EXPERIMENTAL"}
    else:   # no arm labelled experimental: only placebo, sham and no-intervention arms are controls
        ctrl = {l for l, t in types.items() if t in CONTROL - {"ACTIVE_COMPARATOR"}}
        if len(labels - ctrl) >= 2 and not (c_arms >= (labels - ctrl)):
            return "head_to_head", len(arms)
    if c_arms >= labels:
        return "backbone", len(arms)
    if c_arms <= ctrl:
        return "comparator", len(arms)
    if c_arms & ctrl:
        return "mixed", len(arms)
    return "investigational", len(arms)
