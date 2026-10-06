"""Write-back policies, training labels, scorers and metrics for the failure write-back experiment.

Unit of analysis: a Hetionet (compound, disease) pair on the full 1,552 x 137 grid.
Training is positive-unlabelled: recorded fitting indications are positives; every other
pair is an unlabelled pair with weight 1 unless a write-back policy changes its weight.
"""
import gzip, os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

HET = os.environ.get("HET", "work/data/hetionet")
DATA = os.environ.get("DATA", "work/mapped")
STOP = {"TERMINATED", "WITHDRAWN", "SUSPENDED"}
SCIENTIFIC = {"Negative", "Safety_Sideeffects"}
CUTOFF_WRITEBACK = os.environ.get("WB_CUTOFF", "2015-01-01")   # stopped trials starting before this date form the write-back set
WB_PHASE13 = os.environ.get("WB_PHASE13", "0") == "1"           # sensitivity: exclude phase 4 trials from the write-back set
CUTOFF_FUTURE = "2017-01-01"      # trials starting on/after this date are later evidence
ACTIONS = ("ignore", "negate", "mask")


# ----------------------------------------------------------------------------- data
def load_graph():
    nodes = pd.read_csv(f"{HET}/hetionet-v1.0-nodes.tsv", sep="\t")
    nodes["kind"] = nodes.kind.str.strip()
    comps = sorted(nodes[nodes.kind == "Compound"].id.str.replace("Compound::", "", regex=False))
    dises = sorted(nodes[nodes.kind == "Disease"].id.str.replace("Disease::", "", regex=False))
    names = dict(zip(nodes.id.str.split("::").str[-1], nodes.name))
    ctd, cpd = [], set()
    with gzip.open(f"{HET}/hetionet-v1.0-edges.sif.gz", "rt") as fh:
        fh.readline()
        for line in fh:
            s, m, t = line.rstrip("\n").split("\t")
            if m == "CtD":
                ctd.append((s[10:], t[9:]))
            elif m == "CpD":
                cpd.add((s[10:], t[9:]))
    return comps, dises, names, ctd, cpd


def load_evidence():
    """Pair-level evidence from Open Targets clinical reports mapped to Hetionet."""
    tp = pd.read_csv(f"{DATA}/trial_pairs.tsv", sep="\t", low_memory=False)
    dm = pd.read_csv(f"{DATA}/disease_map.tsv", sep="\t")
    gap = dict(zip(zip(dm.ot_disease, dm.disease), dm.depth_gap))
    tp["gap"] = [gap.get((o, d), -1) for o, d in zip(tp.ot_disease, tp.disease)]
    tp["pair"] = list(zip(tp.compound, tp.disease))
    trials = tp[tp.origin == "CLINICAL_TRIAL"].copy()
    trials["start"] = pd.to_datetime(trials.start_date, errors="coerce")
    st = trials[trials.status.isin(STOP)].copy()
    cats = st.stop_categories.fillna("")
    st["scientific"] = cats.apply(lambda s: bool(set(s.split("|")) & SCIENTIFIC) if s else False)
    wb = st[st.start < CUTOFF_WRITEBACK]
    if WB_PHASE13:
        wb = wb[~wb.phase.fillna("").str.contains("PHASE4")]
    later = st[st.start >= CUTOFF_FUTURE]
    approved = set(tp[(tp.stage == "APPROVAL") & (tp.origin != "CLINICAL_TRIAL")].pair)
    ev = {
        "wb_all": set(wb.pair),
        "wb_scientific": set(wb[wb.scientific].pair),
        "wb_scientific_scoped": set(wb[wb.scientific & (wb.gap == 0)].pair),
        "later_scientific": set(later[later.scientific].pair),
        "stopped_any_date": set(st.pair),
        "approved": approved,
        "trials": trials,
        "stopped": st,
    }
    ev["wb_other"] = ev["wb_all"] - ev["wb_scientific"]
    rf = f"{DATA}/trial_roles.tsv"           # role of the compound in each stopped trial (ClinicalTrials.gov arm groups)
    if os.path.exists(rf):
        roles = pd.read_csv(rf, sep="\t")
        wr = wb.merge(roles, on=["report_id", "compound"], how="left")
        inv = wr[wr.scientific & (wr.role == "investigational")]
        ev["wb_scientific_inv"] = set(inv.pair)
        ev["wb_scientific_inv_scoped"] = set(inv[inv.gap == 0].pair)
    return ev


def policy_sets(ev, policy):
    """Return {pair: action} for one policy name.

    Grid policies are named 'sci=<action>,other=<action>'. 'typed_scoped' negates
    scientific stops only when a stopped trial's condition maps to the Hetionet disease
    itself (not to a narrower subtype); all other stopped pairs are masked. 'typed_role'
    negates only scientific stops in which the compound was the investigational agent
    (ClinicalTrials.gov arm groups), and 'typed_role_scoped' adds the scope condition.
    """
    out = {}
    only = {"typed_scoped": "wb_scientific_scoped", "typed_role": "wb_scientific_inv",
            "typed_role_scoped": "wb_scientific_inv_scoped"}
    if policy in only:          # negate the named subset of scientific stops; mask every other stopped pair
        for p in ev["wb_all"]:
            out[p] = "negate" if p in ev[only[policy]] else "mask"
        return out
    sci_a, oth_a = [x.split("=")[1] for x in policy.split(",")]
    for p in ev["wb_scientific"]:
        out[p] = sci_a
    for p in ev["wb_other"]:
        out[p] = oth_a
    return out


NAMED = {
    "no_writeback": "sci=ignore,other=ignore",
    "flat_negative": "sci=negate,other=negate",
    "mask_all": "sci=mask,other=mask",
    "typed": "sci=negate,other=mask",
    "typed_scoped": "typed_scoped",
    "typed_role": "typed_role",                 # scientific stops of the investigational drug only
    "typed_role_scoped": "typed_role_scoped",   # ... and only when the trial condition is the disease concept itself
}
GRID = [f"sci={a},other={b}" for a in ACTIONS for b in ACTIONS]


def training_matrix(comps, dises, fit_pos, cpd, actions, w_neg=10.0, masked_rows=()):
    """Labels and per-pair weights on the full grid.

    Recorded fitting positives always remain positives (they take precedence over a
    stopped trial). CpD (symptomatic) pairs are excluded. Rows of compounds whose labels
    are held out (compound-disjoint task) carry no unlabelled weight, but a write-back
    negative on such a compound is still applied because it is the only evidence there.
    """
    ci = {c: i for i, c in enumerate(comps)}; di = {d: j for j, d in enumerate(dises)}
    Y = np.zeros((len(comps), len(dises)), np.float32)
    Wt = np.ones_like(Y)
    for c in masked_rows:
        Wt[ci[c], :] = 0.0
    for c, d in cpd:
        Wt[ci[c], di[d]] = 0.0
    for (c, d), a in actions.items():
        if c not in ci or d not in di:
            continue
        if a == "negate":
            Wt[ci[c], di[d]] = w_neg
        elif a == "mask":
            Wt[ci[c], di[d]] = 0.0
    for c, d in fit_pos:
        Y[ci[c], di[d]] = 1.0
    pos = Y == 1
    Wt[pos] = 1.0
    neg_w = Wt[~pos].sum()
    Wt[pos] = neg_w / max(pos.sum(), 1)       # balance the two classes
    return Y, Wt


# ----------------------------------------------------------------------------- models
class Scorer(nn.Module):
    """Score matrix = MF term + bilinear graph term + compound/disease biases."""

    def __init__(self, nC, nD, mf_dim=0, emb_c=None, emb_d=None):
        super().__init__()
        self.bc = nn.Parameter(torch.zeros(nC)); self.bd = nn.Parameter(torch.zeros(nD)); self.b0 = nn.Parameter(torch.zeros(()))
        self.mf_dim = mf_dim
        if mf_dim:
            self.U = nn.Parameter(torch.randn(nC, mf_dim) * 0.1); self.V = nn.Parameter(torch.randn(nD, mf_dim) * 0.1)
        self.graph = emb_c is not None
        if self.graph:
            self.register_buffer("Ec", emb_c); self.register_buffer("Ed", emb_d)
            self.M = nn.Parameter(torch.zeros(emb_c.shape[1], emb_d.shape[1]))

    def forward(self):
        s = self.bc[:, None] + self.bd[None, :] + self.b0
        if self.mf_dim:
            s = s + self.U @ self.V.T
        if self.graph:
            s = s + (self.Ec @ self.M) @ self.Ed.T
        return s


def fit_scorer(Y, Wt, mf_dim=0, emb=None, wd=1e-4, lr=0.02, epochs=(200,), seed=1, threads=1):
    torch.manual_seed(seed)
    torch.set_num_threads(threads)
    nC, nD = Y.shape
    ec, ed = (None, None) if emb is None else emb
    model = Scorer(nC, nD, mf_dim, ec, ed)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    y = torch.tensor(Y); w = torch.tensor(Wt)
    wsum = w.sum()
    out = {}
    for ep in range(1, max(epochs) + 1):
        s = model()
        loss = (nn.functional.binary_cross_entropy_with_logits(s, y, reduction="none") * w).sum() / wsum
        opt.zero_grad(); loss.backward(); opt.step()
        if ep in epochs:
            with torch.no_grad():
                out[ep] = model().numpy().copy()
    return out


def degree_reference(comps, dises, Y, Wt):
    """Balanced logistic regression on log(1 + fitting degree) of each endpoint."""
    cdeg = Y.sum(1); ddeg = Y.sum(0)
    X = np.stack(np.broadcast_arrays(np.log1p(cdeg)[:, None], np.log1p(ddeg)[None, :]), -1).reshape(-1, 2)
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Xs = (X - mu) / sd
    y = Y.reshape(-1); w = Wt.reshape(-1)
    keep = w > 0
    clf = LogisticRegression(solver="liblinear", C=1.0, max_iter=1000).fit(Xs[keep], y[keep], sample_weight=w[keep])
    return (Xs @ clf.coef_[0] + clf.intercept_[0]).reshape(Y.shape)


def disease_degree(Y):
    return np.broadcast_to(np.log1p(Y.sum(0))[None, :], Y.shape).copy()


def load_embeddings(path, comps, dises):
    st = torch.load(path)
    idx = {e: i for i, e in enumerate(st["entities"])}
    E = st["er"].float()
    if "ei" in st:
        E = torch.cat([E, st["ei"].float()], 1)
    E = (E - E.mean(0)) / (E.std(0) + 1e-6)
    ec = E[[idx["Compound::" + c] for c in comps]]
    ed = E[[idx["Disease::" + d] for d in dises]]
    return ec, ed


# ----------------------------------------------------------------------------- metrics
def weighted_ap(y, s, w=None):
    """Non-interpolated AP with exact tie handling; optional row weights (bootstrap)."""
    o = np.argsort(-s, kind="stable")
    yy = y[o]; ss = s[o]
    ww = np.ones(len(y)) if w is None else w[o]
    ends = np.r_[np.flatnonzero(np.diff(ss)), len(ss) - 1]
    tp = np.cumsum(ww * yy)[ends]; tot = np.cumsum(ww)[ends]
    if tp[-1] <= 0:
        return float("nan")
    prec = tp / np.maximum(tot, 1e-12)
    return float(np.sum(np.diff(np.r_[0.0, tp]) * prec) / tp[-1])


def grid_metrics(pairs, labels, scores, comps_of_pairs=None):
    y = np.asarray(labels, float); s = np.asarray(scores, float)
    res = {"AP": weighted_ap(y, s), "AUROC": float(roc_auc_score(y, s)), "n": int(len(y)), "n_pos": int(y.sum())}
    df = pd.DataFrame({"c": [p[0] for p in pairs], "d": [p[1] for p in pairs], "y": y, "s": s})
    aps = [weighted_ap(g.y.values, g.s.values) for _, g in df.groupby("d") if g.y.sum() > 0 and g.y.sum() < len(g)]
    res["macroAP_disease"] = float(np.mean(aps))
    rr, h10 = [], []
    for c, g in df.groupby("c"):
        if g.y.sum() == 0:
            continue
        neg = g.s.values[g.y.values == 0]
        for sp in g.s.values[g.y.values == 1]:
            rank = 1 + int((neg > sp).sum()) + 0.5 * int((neg == sp).sum())   # filtered: other positives removed
            rr.append(1.0 / rank); h10.append(rank <= 10)
    res["MRR"] = float(np.mean(rr)); res["Hits@10"] = float(np.mean(h10))
    return res
