#!/usr/bin/env python
"""Experiment B - KGE training (ComplEx / DistMult / RotatE) in pure PyTorch.
For each split scheme: build training graph = all Hetionet edges minus held-out test
CtD edges and a held-out validation fraction (+ inverse relations), train with negative
sampling (logistic/BCE loss -> calibratable scores) and VALIDATION-BASED EARLY STOPPING
(model selection on validation MRR), then (a) score every needed (compound,CtD,disease)
pair for the AUPRC analysis and (b) compute full-candidate filtered rank metrics
(MRR, Hits@1/3/10) on the test positives.
Outputs: scores_<model>_<split><suffix>.tsv  and  rank_<model>_<split><suffix>.json
"""

import gzip, json, os, time, sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

SEED = int(os.environ.get("SEED", 7))
torch.manual_seed(SEED)
np.random.seed(SEED)
ROOT = os.environ["KG_AUDIT_WORKDIR"]
DATA = f"{ROOT}/data/hetionet"
OUT = f"{ROOT}/outputs/hetionet/raw"
os.makedirs(OUT, exist_ok=True)
# CPU is far faster than MPS for these small embedding ops (MPS per-op overhead dominates)
DEV = torch.device(os.environ.get("DEV", "cpu"))
torch.set_num_threads(int(os.environ.get("THREADS", max(1, os.cpu_count() - 1))))
DIM = 100
BATCH = int(os.environ.get("BATCH", 32768))
NNEG = 4
LR = 3e-3
L2 = 1e-6
# validation-based early stopping (equal model-selection budget across all models)
MAX_EPOCHS = int(os.environ.get("EPOCHS", os.environ.get("MAX_EPOCHS", 30)))
MIN_EPOCHS = 8
EVAL_EVERY = 2
PATIENCE = 3
VAL_FRAC = 0.10
ONLY = os.environ.get("ONLY", "")  # e.g. "random:complex" to run a single cell for smoke test
SUFFIX = os.environ.get("SUFFIX", "")  # appended to output filename (multi-seed)
RUN_SPLITS = os.environ.get("SPLITS", "random,coldstart,scaffold").split(",")
RUN_MODELS = os.environ.get("MODELS", "complex,distmult").split(",")

# ---------- load graph ----------
edges = []
with gzip.open(f"{DATA}/hetionet-edges.sif.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        s, m, t = line.rstrip("\n").split("\t")
        edges.append((s, m, t))
edges = np.array(edges, dtype=object)
# entity universe = ALL nodes (some disease/compound nodes are isolated -> still need embeddings)
_nodes = pd.read_csv(f"{DATA}/hetionet-v1.0-nodes.tsv", sep="\t")
ents = sorted(set(_nodes.id) | set(edges[:, 0]) | set(edges[:, 2]))
rels = sorted(set(edges[:, 1]))
e2i = {e: i for i, e in enumerate(ents)}
r2i = {r: i for i, r in enumerate(rels)}  # base relations
CtD = "CtD"
nE, nR = len(ents), len(rels)
print(f"entities={nE} relations={nR} dev={DEV} max_epochs={MAX_EPOCHS}", flush=True)

H = np.array([e2i[s] for s in edges[:, 0]], dtype=np.int64)
R = np.array([r2i[m] for m in edges[:, 1]], dtype=np.int64)
T = np.array([e2i[t] for t in edges[:, 2]], dtype=np.int64)
ctd_mask = edges[:, 1] == CtD
# full-candidate disease universe (all diseases that appear as a CtD target)
CAND = np.array(sorted({e2i[edges[i, 2]] for i in np.where(ctd_mask)[0]}), dtype=np.int64)

splits = json.load(open(f"{OUT}/splits.json"))


def to_idx_pairs(pairs):
    return [(e2i[c], e2i[d]) for c, d in pairs]


# ---------- models ----------
class KGE(nn.Module):
    def __init__(self, nE, nR, dim, kind="complex"):
        super().__init__()
        self.kind = kind
        self.er = nn.Embedding(nE, dim)
        self.rr = nn.Embedding(nR, dim)
        nn.init.xavier_uniform_(self.er.weight)
        nn.init.xavier_uniform_(self.rr.weight)
        if kind in ("complex", "rotate"):
            self.ei = nn.Embedding(nE, dim)
            nn.init.xavier_uniform_(self.ei.weight)
        if kind == "complex":
            self.ri = nn.Embedding(nR, dim)
            nn.init.xavier_uniform_(self.ri.weight)
        if kind == "rotate":
            # relation phase in (-pi, pi]; constrains |r|=1 (a rotation), per Sun et al. 2019
            nn.init.uniform_(self.rr.weight, -3.14159265, 3.14159265)

    def score(self, h, r, t):
        hr, rr_, tr = self.er(h), self.rr(r), self.er(t)
        if self.kind == "distmult":
            return (hr * rr_ * tr).sum(-1)
        if self.kind == "rotate":
            hi, ti = self.ei(h), self.ei(t)
            r_re, r_im = torch.cos(rr_), torch.sin(rr_)  # unit-modulus rotation
            rot_re = hr * r_re - hi * r_im  # complex Hadamard (h ∘ r)
            rot_im = hr * r_im + hi * r_re
            dist = torch.sqrt((rot_re - tr) ** 2 + (rot_im - ti) ** 2 + 1e-9).sum(-1)
            return 6.0 - dist  # gamma - distance -> higher = better
        hi, ri_, ti = self.ei(h), self.ri(r), self.ei(t)
        return (hr * rr_ * tr + hr * ri_ * ti + hi * rr_ * ti - hi * ri_ * tr).sum(-1)


@torch.no_grad()
def rank_metrics(model, queries, known):
    """Full-candidate FILTERED ranking of each test positive's true disease against CAND.
    queries: list of (c_idx, d_idx); known: dict c_idx -> set of training d_idx to filter out."""
    if not queries:
        return {"MRR": float("nan"), "Hits@1": float("nan"), "Hits@3": float("nan"), "Hits@10": float("nan"), "n": 0}
    ct = torch.tensor(CAND, device=DEV)
    r = torch.full((len(CAND),), r2i[CtD], dtype=torch.long, device=DEV)
    pos = {int(d): j for j, d in enumerate(CAND)}
    rr = []
    hits = {1: 0, 3: 0, 10: 0}
    for c, dt in queries:
        h = torch.full((len(CAND),), c, dtype=torch.long, device=DEV)
        s = model.score(h, r, ct).cpu().numpy().copy()
        for d in known.get(c, ()):  # filtered setting
            if d != dt and d in pos:
                s[pos[d]] = -1e9
        jt = pos[dt]
        rank = 1 + int((s > s[jt]).sum())
        rr.append(1.0 / rank)
        for k in hits:
            if rank <= k:
                hits[k] += 1
    n = len(queries)
    return {"MRR": float(np.mean(rr)), **{f"Hits@{k}": hits[k] / n for k in hits}, "n": n}


def train_model(kind, drop_idx, val_queries, known):
    """drop_idx: edges to EXCLUDE from training (held-out test CtD + validation CtD).
    Early-stops on validation MRR (model selection); restores the best epoch's weights."""
    keep = ~drop_idx
    h = torch.tensor(H[keep])
    r = torch.tensor(R[keep])
    t = torch.tensor(T[keep])
    h2 = torch.cat([h, t])
    t2 = torch.cat([t, h])
    r2 = torch.cat([r, r + nR])  # inverse relations
    h2, r2, t2 = h2.to(DEV), r2.to(DEV), t2.to(DEV)
    nTr = h2.shape[0]
    model = KGE(nE, 2 * nR, DIM, kind).to(DEV)
    opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=L2)
    bce = nn.BCEWithLogitsLoss()
    best_mrr, best_state, wait, best_ep = -1.0, None, 0, 0
    for ep in range(MAX_EPOCHS):
        perm = torch.randperm(nTr, device=DEV)
        tot = 0.0
        for b in range(0, nTr, BATCH):
            idx = perm[b : b + BATCH]
            ph, pr, pt = h2[idx], r2[idx], t2[idx]
            bs = ph.shape[0]
            neg_t = torch.randint(0, nE, (bs, NNEG), device=DEV)
            ph_e = ph.unsqueeze(1).expand(bs, NNEG)
            pr_e = pr.unsqueeze(1).expand(bs, NNEG)
            pos_s = model.score(ph, pr, pt)
            neg_s = model.score(ph_e.reshape(-1), pr_e.reshape(-1), neg_t.reshape(-1)).reshape(bs, NNEG)
            scores = torch.cat([pos_s.unsqueeze(1), neg_s], dim=1)
            labels = torch.zeros_like(scores)
            labels[:, 0] = 1.0
            loss = bce(scores, labels)
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item() * bs
        if ep + 1 >= MIN_EPOCHS and (ep % EVAL_EVERY == 0 or ep == MAX_EPOCHS - 1):
            vm = rank_metrics(model, val_queries, known)["MRR"]
            if vm > best_mrr + 1e-4:
                best_mrr, best_state, wait, best_ep = (
                    vm,
                    {k: v.detach().clone() for k, v in model.state_dict().items()},
                    0,
                    ep,
                )
            else:
                wait += 1
            if wait >= PATIENCE:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    print(f"   {kind}: stopped ep{ep} best_ep{best_ep} val_MRR={best_mrr:.3f}", flush=True)
    return model


@torch.no_grad()
def score_pairs(model, pairs):
    if not pairs:
        return np.array([])
    cd = np.array(pairs, dtype=np.int64)
    h = torch.tensor(cd[:, 0]).to(DEV)
    t = torch.tensor(cd[:, 1]).to(DEV)
    r = torch.full((len(pairs),), r2i[CtD], dtype=torch.long, device=DEV)
    return model.score(h, r, t).cpu().numpy()


for split in RUN_SPLITS:
    sset = splits[split]
    test_set = set(map(tuple, sset["test_pos"]))
    # all CtD (c,d) pairs and the training pool (CtD not held out for test)
    all_ctd = [(edges[i, 0], edges[i, 2]) for i in np.where(ctd_mask)[0]]
    train_ctd = [p for p in all_ctd if p not in test_set]
    # validation hold-out for model selection (seeded, per split)
    vrng = np.random.default_rng(SEED * 100 + len(split))
    nval = max(10, int(VAL_FRAC * len(train_ctd)))
    vsel = set(map(int, vrng.choice(len(train_ctd), size=nval, replace=False)))
    val_pairs = [train_ctd[i] for i in vsel]
    graph_ctd = [train_ctd[i] for i in range(len(train_ctd)) if i not in vsel]  # actually in graph
    # edges to drop from training = test CtD edges + validation CtD edges
    hold = test_set | set(val_pairs)
    drop = np.zeros(len(edges), dtype=bool)
    for idx in np.where(ctd_mask)[0]:
        if (edges[idx, 0], edges[idx, 2]) in hold:
            drop[idx] = True
    print(
        f"[{split}] drop {drop.sum()} CtD edges (test {len(test_set)} + val {len(val_pairs)}) from training", flush=True
    )
    # filtered-ranking known set: training-graph CtD targets per compound
    known = {}
    for c, d in graph_ctd:
        known.setdefault(e2i[c], set()).add(e2i[d])
    val_q = [(e2i[c], e2i[d]) for c, d in val_pairs]
    test_q = [(e2i[c], e2i[d]) for c, d in sset["test_pos"]]
    pack = {
        k: to_idx_pairs([tuple(x) for x in sset[k]])
        for k in ["test_pos", "neg_random", "neg_degmatch", "train_pos", "train_neg"]
    }
    for kind in RUN_MODELS:
        if ONLY and ONLY != f"{split}:{kind}":
            continue
        t0 = time.time()
        model = train_model(kind, drop, val_q, known)
        rows = []
        for grp, prs in pack.items():
            sc = score_pairs(model, prs)
            for (c, d), s in zip(sset[grp], sc):
                rows.append((grp, c, d, float(s)))
        pd.DataFrame(rows, columns=["group", "compound", "disease", "score"]).to_csv(
            f"{OUT}/scores_{kind}_{split}{SUFFIX}.tsv", sep="\t", index=False
        )
        rm = rank_metrics(model, test_q, known)  # full-candidate filtered rank metrics
        json.dump(
            {"model": kind, "split": split, "seed": SEED, **rm},
            open(f"{OUT}/rank_{kind}_{split}{SUFFIX}.json", "w"),
            indent=2,
        )
        print(
            f"[{split}/{kind}] {time.time() - t0:.0f}s -> MRR {rm['MRR']:.3f} Hits@10 {rm['Hits@10']:.3f}", flush=True
        )
print("KGE_DONE")
