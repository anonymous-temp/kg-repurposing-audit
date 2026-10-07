"""Shared helpers for the supplementary figures: markdown-table parser, label normalisation, paths.

All numbers in the supplementary figures are parsed from the markdown tables of the manuscript
(read only); nothing is typed by hand."""
import math, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # code/ (figstyle.py)
import figstyle as S                                # noqa: E402

MS = os.environ.get("SUPP_MD", "paper_results/supp_tables_md")   # generated supplementary tables (markdown)
F_S3S6, F_S16S18, F_S22, F_S7S14 = (f"{MS}/{n}" for n in ("supp_S3_S6.md", "supp_S16_S18.md", "supp_S22.md", "supp_S7_S14.md"))
OUT = os.environ.get("FIG_OUT", "work/figures_supp")
os.makedirs(OUT, exist_ok=True)

NAN = float("nan")
MINUS = "−"

# canonical policy keys (display order of the manuscript) and labels
POL_ORDER = ["no_writeback", "flat_negative", "typed", "typed_scoped", "typed_role", "typed_role_scoped", "typed_lowint2", "mask_all"]
PL = {"no_writeback": "No write-back", "flat_negative": "Flat negative", "typed": "Typed", "typed_scoped": "Typed + scope",
      "typed_role": "Typed + role", "typed_role_scoped": "Typed + role + scope", "typed_lowint2": "Typed + ≤2 trials",
      "mask_all": "Mask all", "typed_lowint1": "Typed + ≤1 trial", "typed_lowint4": "Typed + ≤4 trials"}
OUTCOME_KEYS = ["random", "compound", "e2"]


def _canon(label):
    return re.sub(r"[\s+_\-]", "", label.lower())


_POL = {"nowriteback": "no_writeback", "flatnegative": "flat_negative", "typed": "typed", "typedscope": "typed_scoped",
        "typedscoped": "typed_scoped", "typedrole": "typed_role", "typedrolescope": "typed_role_scoped",
        "typedrolescoped": "typed_role_scoped", "maskall": "mask_all", "typed≤2trials": "typed_lowint2",
        "typedlowint2": "typed_lowint2", "typed≤1trial": "typed_lowint1", "typed≤4trials": "typed_lowint4",
        # grid names used in some tables
        "sciignore/otherignore": "no_writeback", "scinegate/othernegate": "flat_negative",
        "scinegate/othermask": "typed", "scimask/othermask": "mask_all"}
_SCORER = {"degreereference": "degree", "degree": "degree", "mf": "mf", "labelonlymf": "mf", "labelonly": "mf",
           "graphhead": "graph", "graph": "graph", "hybrid": "hybrid", "diseasedegree": "disease_degree"}


def policy_key(label):
    k = _canon(label)
    if k not in _POL:
        raise KeyError(f"unknown policy label {label!r}")
    return _POL[k]


def scorer_key(label):
    k = _canon(label)
    if k not in _SCORER:
        raise KeyError(f"unknown scorer label {label!r}")
    return _SCORER[k]


def outcome_key(label):
    s = label.lower()
    if "random" in s:
        return "random"
    if "compound" in s:
        return "compound"
    if "e2" in s or "external" in s:
        return "e2"
    raise KeyError(f"unknown outcome label {label!r}")


# ------------------------------------------------------------------------------------------ parsing
def read_table(path, marker):
    """Return (header, rows) of the first markdown table after the line that starts with `marker`
    (ignoring leading '**'), e.g. marker='Table S4a.'."""
    lines = open(path, encoding="utf-8").read().splitlines()
    start = None
    for i, ln in enumerate(lines):
        if ln.lstrip("*").startswith(marker):
            start = i
            break
    if start is None:
        raise ValueError(f"{marker!r} not found in {path}")
    j = start + 1
    while j < len(lines) and not lines[j].startswith("|"):
        j += 1
    block = []
    while j < len(lines) and lines[j].startswith("|"):
        block.append([c.strip() for c in lines[j].strip().strip("|").split("|")])
        j += 1
    header, body = block[0], [r for r in block[2:]]
    assert all(len(r) == len(header) for r in body), f"ragged table {marker}"
    return header, body


def table_dicts(path, marker):
    header, rows = read_table(path, marker)
    return header, [dict(zip(header, r)) for r in rows]


def num(s):
    """'+0.049', '−0.095', '1,752', '' -> float (nan if empty / not numeric)."""
    s = s.strip().replace(MINUS, "-").replace(",", "").replace("+", "")
    if s in ("", "-", "–"):
        return NAN
    try:
        return float(s)
    except ValueError:
        return NAN


def mean_sd(s):
    """'0.336 (0.019)' -> (0.336, 0.019); '0.256' -> (0.256, nan)."""
    m = re.match(r"^\s*([+\-−]?[\d.]+)\s*(?:\(\s*([\d.]+)\s*\))?\s*$", s)
    if not m:
        return NAN, NAN
    return num(m.group(1)), (num(m.group(2)) if m.group(2) else NAN)


_NUM = r"[+\-−]?[\d.]+"
_RX_CI = re.compile(rf"^\s*({_NUM})\s*\(\s*({_NUM})\s*,\s*({_NUM})\s*\)\s*$")


def diff_ci(s):
    """'-0.095 (-0.125, -0.069)' (ASCII or U+2212 minus) -> (diff, lo, hi)."""
    m = _RX_CI.match(s)
    if not m:
        raise ValueError(f"cannot parse diff/CI from {s!r}")
    return num(m.group(1)), num(m.group(2)), num(m.group(3))


def isnan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


# ------------------------------------------------------------------------------------------ contrasts
_RX_VS = re.compile(r"^(?P<pol>.+?) vs (?P<base>.+?) \((?P<sc>[^)]+)\)$")
_RX_BR = re.compile(r"^(?P<pol>.+?) - (?P<base>\S+) \[(?P<sc>[^\]]+)\]$")


def parse_contrast(label):
    """Both label styles of the supplement -> (policy, baseline, scorer), all canonical keys.
    'Flat negative vs no write-back (graph head)'  and  'flat_negative - no_writeback [graph]'."""
    m = _RX_VS.match(label) or _RX_BR.match(label)
    if not m:
        raise ValueError(f"cannot parse contrast {label!r}")
    pol, base = m.group("pol"), m.group("base")
    sc = m.group("sc")
    return policy_key(pol), policy_key(base), scorer_key(sc)


# ------------------------------------------------------------------------------------------ absolute results (S4, S16a)
# column names in the tables -> measure keys
MEAS = {"Pooled AP": "pooled", "Per-disease AP": "perdis", "MRR": "mrr", "Hits@10": "hits10", "AUROC": "auroc"}


def load_s4(marker):
    """Hetionet Table S4a/b/c -> {(policy, scorer): {measure: (mean, sd)}} for the named policies at negative weight 10.
    The 'Sci X / other Y' grid rows are skipped except where the name is one of the named policies."""
    _, rows = table_dicts(F_S3S6, marker)
    out = {}
    for r in rows:
        lab = r["Policy"]
        if lab.lower().startswith("sci "):
            continue                                   # grid rows, not a named policy
        if num(r["Negative weight"]) != 10.0:
            continue
        pol, sc = policy_key(lab), scorer_key(r["Scorer"])
        key = (pol, sc)
        assert key not in out, f"duplicate {key} in {marker}"
        out[key] = {MEAS[c]: mean_sd(r[c]) for c in MEAS if c in r}
    return out


S16_SET = {"E1 random edge": "random", "E1 compound disjoint": "compound", "E2 external approvals": "e2"}


def load_s16a():
    """PrimeKG Table S16a -> {outcome: {(policy, scorer): {measure: (mean, sd)}}}; E3 AUROC under key 'e3auroc'."""
    _, rows = table_dicts(F_S16S18, "Table S16a.")
    out = {k: {} for k in S16_SET.values()}
    for r in rows:
        o = S16_SET[r["Set"]]
        d = {MEAS[c]: mean_sd(r[c]) for c in MEAS if c in r}
        d["e3auroc"] = mean_sd(r["E3 AUROC"])
        key = (policy_key(r["Policy"]), scorer_key(r["Scorer"]))
        assert key not in out[o]
        out[o][key] = d
    return out


# ------------------------------------------------------------------------------------------ contrasts vs no write-back
def _collect(rows, setcol, out, outcome_of, want_base="no_writeback"):
    for r in rows:
        try:
            pol, base, sc = parse_contrast(r["Contrast"])
        except KeyError:                                # between-scorer contrast ('graph - mf [no_writeback]')
            continue
        if base != want_base or pol == "no_writeback":
            continue
        o = outcome_of(r[setcol])
        key = (o, pol, sc)
        if "Per-disease AP Δ (95% CI)" in r:
            vals = {"pooled": diff_ci(r["Pooled AP Δ (95% CI)"]), "perdis": diff_ci(r["Per-disease AP Δ (95% CI)"])}
        else:                                           # Table S5 layout: separate diff and CI columns
            def pair(dcol, ccol):
                lo, hi = (num(x) for x in r[ccol].strip("()").split(","))
                return num(r[dcol]), lo, hi
            vals = {"pooled": pair("Pooled AP difference", "Pooled AP 95% CI"),
                    "perdis": pair("Per-disease AP difference", "Per-disease AP 95% CI")}
        assert key not in out, f"duplicate contrast {key}"
        out[key] = vals


def contrasts_hetionet():
    """{(outcome, policy, scorer): {'pooled'|'perdis': (diff, lo, hi)}}: Table S5 (flat, typed, typed+scope, mask all),
    Table S18b (role policies) and Table S22b (typed + <=1/2/4 trials)."""
    out = {}
    _, rows = table_dicts(F_S3S6, "Table S5.")
    _collect(rows, "Set", out, outcome_key)
    _, rows = table_dicts(F_S16S18, "Table S18b.")
    _collect(rows, "Set", out, outcome_key)
    _add_intensity(out, "Hetionet")
    return out


def contrasts_primekg():
    out = {}
    _, rows = table_dicts(F_S16S18, "Table S16b.")
    _collect(rows, "Set", out, outcome_key)
    _add_intensity(out, "PrimeKG")
    return out


def _add_intensity(out, graph):
    """Typed + <=k trials: per-disease AP change vs no write-back from Table S22b."""
    for r in intensity_rows():
        if r["graph"] != graph:
            continue
        key = (r["outcome"], r["policy"], r["scorer"])
        if key in out:
            continue
        out[key] = {"perdis": r["d_nowb"]}


# ------------------------------------------------------------------------------------------ S22
def placebo_rows():
    """Table S22a -> list of dicts (graph, outcome, scorer, set, n, real, degree, uniform) with (diff, lo, hi) triples."""
    _, rows = table_dicts(F_S22, "Table S22a.")
    out = []
    for r in rows:
        out.append(dict(graph=r["Graph"], outcome=outcome_key(r["Outcome"]), scorer=scorer_key(r["Scorer"]),
                        set=policy_key(r["Negated set"]), n=num(r["Negatives added"]),
                        real=diff_ci(r["Stopped pairs: Δ vs masking (95% CI)"]),
                        degree=diff_ci(r["Placebo, same drug and disease counts: Δ vs masking (95% CI)"]),
                        uniform=diff_ci(r["Placebo, uniform: Δ vs masking (95% CI)"])))
    return out


def intensity_rows():
    """Table S22b -> list of dicts."""
    _, rows = table_dicts(F_S22, "Table S22b.")
    out = []
    for r in rows:
        m = re.match(r"^\s*([\d,]+)\s*\((\d+)%\)\s*$", r["Negated pairs (recorded or approved)"])
        assert m, r["Negated pairs (recorded or approved)"]
        out.append(dict(graph=r["Graph"], policy=policy_key(r["Policy"]), n=num(m.group(1)), pct=num(m.group(2)),
                        outcome=outcome_key(r["Outcome"]), scorer=scorer_key(r["Scorer"]),
                        ap=num(r["Per-disease AP"]), d_nowb=diff_ci(r["Δ vs no write-back (95% CI)"]),
                        d_mask=diff_ci(r["Δ vs masking (95% CI)"])))
    return out
