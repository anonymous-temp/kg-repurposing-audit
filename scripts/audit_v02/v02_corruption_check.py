#!/usr/bin/env python
"""Check of the v0.2.0 tail-corruption scheme (Supplementary S8): all-entity vs type-constrained corruption.

Re-implements the training loop of pipelines/hetionet_train_kge.py (kg-repurposing-audit
v0.2.0, commit fea88c7) with the same hyper-parameters (dim 100, BCE, Adam 3e-3, wd 1e-6,
4 corrupted tails, batch 32768, validation-MRR checkpoint selection, <=20 epochs) and adds
ONE switch:
  NEG=uniform : corrupted tails drawn from all 47,031 entities (as in the authors' code)
  NEG=typed   : corrupted tails drawn from entities of the same node type as the true tail
Evaluation uses the reconstructed v0.2.0 split and candidates (reconstruct_v02_candidates.py).
Idempotent: a finished (split, model, neg, seed) cell is skipped.
"""
import gzip, json, os, sys, time
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from kg_audit.metrics import fit_degree_reference

DATA = os.environ.get("HET", "work/data/hetionet")
REPRO = os.environ.get("REPRO", "work/audit_v02")
OUT = os.environ.get("OUT", "work/audit_v02/kge_check")
os.makedirs(OUT, exist_ok=True)
SPLIT = os.environ.get("SPLIT", "random")
KIND = os.environ.get("MODEL", "rotate")
NEG = os.environ.get("NEG", "uniform")
SEED = int(os.environ.get("SEED", 1))
MAX_EPOCHS = int(os.environ.get("MAX_EPOCHS", 20))
DIM, BATCH, NNEG, LR, L2 = 100, 32768, 4, 3e-3, 1e-6
MIN_EPOCHS, EVAL_EVERY, PATIENCE, VAL_FRAC = 8, 2, 3, 0.10
tag = f"{SPLIT}_{KIND}_{NEG}_s{SEED}_e{MAX_EPOCHS}"
if os.path.exists(f"{OUT}/{tag}.json"):
    print("skip (done):", tag)
    sys.exit(0)
torch.manual_seed(SEED)
np.random.seed(SEED)
torch.set_num_threads(int(os.environ.get("THREADS", 3)))

nodes = pd.read_csv(f"{DATA}/hetionet-v1.0-nodes.tsv", sep="\t")
ents = sorted(set(nodes.id))
e2i = {e: i for i, e in enumerate(ents)}
kind_of = dict(zip(nodes.id, nodes.kind))
kinds = sorted(set(nodes.kind))
k2i = {k: i for i, k in enumerate(kinds)}
etype = np.array([k2i[kind_of[e]] for e in ents])
H, R, T, rels = [], [], [], {}
with gzip.open(f"{DATA}/hetionet-v1.0-edges.sif.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        s, m, t = line.rstrip("\n").split("\t")
        H.append(e2i[s]); T.append(e2i[t]); R.append(rels.setdefault(m, len(rels)))
# the authors index relations in sorted order; the order does not affect training
H, R, T = np.array(H), np.array(R), np.array(T)
nE, nR = len(ents), len(rels)
CtD = rels["CtD"]
ctd_mask = R == CtD
i2e = ents
all_ctd = [(i2e[H[i]], i2e[T[i]]) for i in np.where(ctd_mask)[0]]
CAND = np.array(sorted({T[i] for i in np.where(ctd_mask)[0]}), dtype=np.int64)

spec = json.load(open(f"{REPRO}/splits.json"))[SPLIT]
test = [tuple(x) for x in spec["test_pos"]]
test_set = set(test)
train_ctd = [p for p in all_ctd if p not in test_set]
assert train_ctd == [tuple(x) for x in spec["train_pos"]]
vrng = np.random.default_rng(SEED * 100 + len(SPLIT))
nval = max(10, int(VAL_FRAC * len(train_ctd)))
vsel = set(map(int, vrng.choice(len(train_ctd), size=nval, replace=False)))
val_pairs = [train_ctd[i] for i in vsel]
graph_ctd = [train_ctd[i] for i in range(len(train_ctd)) if i not in vsel]
hold = test_set | set(val_pairs)
drop = np.zeros(len(H), dtype=bool)
for idx in np.where(ctd_mask)[0]:
    if (i2e[H[idx]], i2e[T[idx]]) in hold:
        drop[idx] = True
known = {}
for c, d in graph_ctd:
    known.setdefault(e2i[c], set()).add(e2i[d])
val_q = [(e2i[c], e2i[d]) for c, d in val_pairs]
test_q = [(e2i[c], e2i[d]) for c, d in test]
print(f"[{tag}] entities={nE} relations={nR} drop={drop.sum()} (test {len(test)} + val {len(val_pairs)})", flush=True)

# typed sampling tables: entities grouped by node type
order = np.argsort(etype, kind="stable")
counts = np.bincount(etype, minlength=len(kinds))
starts = np.r_[0, np.cumsum(counts)[:-1]]
order_t, counts_t, starts_t, etype_t = map(torch.tensor, (order, counts, starts, etype))


class KGE(nn.Module):
    def __init__(self, nE, nR, dim, kind):
        super().__init__()
        self.kind = kind
        self.er = nn.Embedding(nE, dim); self.rr = nn.Embedding(nR, dim)
        nn.init.xavier_uniform_(self.er.weight); nn.init.xavier_uniform_(self.rr.weight)
        if kind in ("complex", "rotate"):
            self.ei = nn.Embedding(nE, dim); nn.init.xavier_uniform_(self.ei.weight)
        if kind == "complex":
            self.ri = nn.Embedding(nR, dim); nn.init.xavier_uniform_(self.ri.weight)
        if kind == "rotate":
            nn.init.uniform_(self.rr.weight, -3.14159265, 3.14159265)

    def score(self, h, r, t):
        hr, rr_, tr = self.er(h), self.rr(r), self.er(t)
        if self.kind == "distmult":
            return (hr * rr_ * tr).sum(-1)
        if self.kind == "rotate":
            hi, ti = self.ei(h), self.ei(t)
            r_re, r_im = torch.cos(rr_), torch.sin(rr_)
            rot_re = hr * r_re - hi * r_im
            rot_im = hr * r_im + hi * r_re
            return 6.0 - torch.sqrt((rot_re - tr) ** 2 + (rot_im - ti) ** 2 + 1e-9).sum(-1)
        hi, ri_, ti = self.ei(h), self.ri(r), self.ei(t)
        return (hr * rr_ * tr + hr * ri_ * ti + hi * rr_ * ti - hi * ri_ * tr).sum(-1)


@torch.no_grad()
def rank_metrics(model, queries):
    ct = torch.tensor(CAND); r = torch.full((len(CAND),), CtD, dtype=torch.long)
    pos = {int(d): j for j, d in enumerate(CAND)}
    rr, hits = [], {1: 0, 3: 0, 10: 0}
    for c, dt in queries:
        s = model.score(torch.full((len(CAND),), c, dtype=torch.long), r, ct).numpy().copy()
        for d in known.get(c, ()):
            if d != dt and d in pos:
                s[pos[d]] = -1e9
        rank = 1 + int((s > s[pos[dt]]).sum())
        rr.append(1.0 / rank)
        for k in hits:
            hits[k] += rank <= k
    return {"MRR": float(np.mean(rr)), **{f"Hits@{k}": hits[k] / len(queries) for k in hits}}


keep = ~drop
h, r, t = torch.tensor(H[keep]), torch.tensor(R[keep]), torch.tensor(T[keep])
h2, t2, r2 = torch.cat([h, t]), torch.cat([t, h]), torch.cat([r, r + nR])
nTr = h2.shape[0]
model = KGE(nE, 2 * nR, DIM, KIND)
opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=L2)
bce = nn.BCEWithLogitsLoss()
best_mrr, best_state, wait, best_ep, hist = -1.0, None, 0, 0, []
t0 = time.time()
disease_type = k2i["Disease"]
n_ctd_neg = n_ctd_neg_disease = 0
for ep in range(MAX_EPOCHS):
    perm = torch.randperm(nTr)
    for b in range(0, nTr, BATCH):
        idx = perm[b : b + BATCH]
        ph, pr, pt = h2[idx], r2[idx], t2[idx]
        bs = ph.shape[0]
        if NEG == "uniform":
            neg_t = torch.randint(0, nE, (bs, NNEG))
        else:
            tt = etype_t[pt]
            u = torch.rand(bs, NNEG)
            neg_t = order_t[starts_t[tt].unsqueeze(1) + (u * counts_t[tt].unsqueeze(1)).long().clamp_(max=nE - 1)]
        m = pr == CtD
        if m.any():
            n_ctd_neg += int(m.sum()) * NNEG
            n_ctd_neg_disease += int((etype_t[neg_t[m]] == disease_type).sum())
        pos_s = model.score(ph, pr, pt)
        neg_s = model.score(ph.unsqueeze(1).expand(bs, NNEG).reshape(-1), pr.unsqueeze(1).expand(bs, NNEG).reshape(-1),
                            neg_t.reshape(-1)).reshape(bs, NNEG)
        scores = torch.cat([pos_s.unsqueeze(1), neg_s], dim=1)
        labels = torch.zeros_like(scores); labels[:, 0] = 1.0
        loss = bce(scores, labels)
        opt.zero_grad(); loss.backward(); opt.step()
    line = f"  ep{ep} loss={loss.item():.4f} elapsed={time.time() - t0:.0f}s"
    if ep + 1 >= MIN_EPOCHS and (ep % EVAL_EVERY == 0 or ep == MAX_EPOCHS - 1):
        vm = rank_metrics(model, val_q)["MRR"]
        hist.append({"epoch": ep, "val_MRR": vm})
        line += f" val_MRR={vm:.4f}"
        if vm > best_mrr + 1e-4:
            best_mrr, best_state, wait, best_ep = vm, {k: v.detach().clone() for k, v in model.state_dict().items()}, 0, ep
        else:
            wait += 1
        print(line, flush=True)
        if wait >= PATIENCE:
            break
    else:
        print(line, flush=True)
if best_state is not None:
    model.load_state_dict(best_state)


@torch.no_grad()
def score_pairs(pairs):
    cd = np.array([(e2i[c], e2i[d]) for c, d in pairs])
    return model.score(torch.tensor(cd[:, 0]), torch.full((len(cd),), CtD, dtype=torch.long), torch.tensor(cd[:, 1])).numpy()


res = {"tag": tag, "split": SPLIT, "model": KIND, "neg": NEG, "seed": SEED, "max_epochs": MAX_EPOCHS, "best_epoch": best_ep,
       "stopped_epoch": ep, "best_val_MRR": best_mrr, "history": hist, "train_seconds": time.time() - t0,
       "ctd_forward_corruptions": n_ctd_neg, "ctd_forward_corruptions_that_are_diseases": n_ctd_neg_disease,
       "test_rank": rank_metrics(model, test_q), "cells": {}}
train_neg = [tuple(p) for p in spec["train_neg"]]
for j, nk in enumerate(["neg_random", "neg_degmatch"]):
    pool = [tuple(p) for p in spec[nk]]
    sel = np.random.default_rng(20261003 + len(SPLIT) + j).choice(len(pool), 20 * len(test), replace=False)
    ng = [pool[i] for i in sel]
    pairs = test + ng
    y = np.r_[np.ones(len(test)), np.zeros(len(ng))]
    s = score_pairs(pairs)
    pred, cdeg, ddeg = fit_degree_reference(graph_ctd, [p for p in train_neg if p not in set(pairs)], SEED)
    dref = pred(pairs)
    donly = np.array([np.log1p(ddeg[d]) for c, d in pairs])
    res["cells"][nk] = {"model_AP": float(average_precision_score(y, s)), "model_AUROC": float(roc_auc_score(y, s)),
                        "degree_AP": float(average_precision_score(y, dref)), "degree_AUROC": float(roc_auc_score(y, dref)),
                        "disease_only_AP": float(average_precision_score(y, donly))}
json.dump(res, open(f"{OUT}/{tag}.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "history"}, indent=1), flush=True)
