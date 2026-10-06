#!/usr/bin/env python
"""Knowledge-graph embedding training on Hetionet v1.0 for the revised analysis.

MODE=pretrain : every Compound-treats-Disease (CtD) and Compound-palliates-Disease (CpD)
                edge is removed, so the embeddings never see an indication label. Entity
                embeddings are saved (float16) at fixed epochs for the supervised heads.
MODE=e2e      : end-to-end ablation in the style of the original pipeline. Test and
                validation CtD edges of one partition are removed and the model is scored on
                the full evaluation grid. Used only to show sensitivity to tail corruption.

NEG=typed   : corrupted tails are drawn from entities of the true tail's node type
NEG=uniform : corrupted tails are drawn from all 47,031 entities (original recipe)

Hyper-parameters follow the original pipeline (dim 100, BCE, Adam 3e-3, weight decay 1e-6,
4 corruptions per triple, batch 32768, forward and inverse relations).
Idempotent: finished outputs are skipped.
"""
import gzip, json, os, sys, time
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

HET = os.environ.get("HET", "work/data/hetionet")
OUT = os.environ.get("OUT", "work/kge")
MODE = os.environ.get("MODE", "pretrain")
KIND = os.environ.get("MODEL", "distmult")
NEG = os.environ.get("NEG", "typed")
SEED = int(os.environ.get("SEED", 1))
EPOCHS = int(os.environ.get("EPOCHS", 20))
SAVE_AT = [int(x) for x in os.environ.get("SAVE_AT", "5,10,15,20").split(",")]
PARTITION = os.environ.get("PARTITION", "")          # e2e only: path to a partition json
DIM, BATCH, NNEG, LR, L2 = 100, 32768, 4, 3e-3, 1e-6
os.makedirs(OUT, exist_ok=True)
tag = f"{MODE}_{KIND}_{NEG}_s{SEED}" + (f"_{os.path.basename(PARTITION).replace('.json', '')}" if MODE == "e2e" else "")
done_flag = f"{OUT}/{tag}.done.json"
if os.path.exists(done_flag):
    print("skip (done):", tag)
    sys.exit(0)
torch.manual_seed(SEED)
np.random.seed(SEED)
torch.set_num_threads(int(os.environ.get("THREADS", 2)))

nodes = pd.read_csv(f"{HET}/hetionet-v1.0-nodes.tsv", sep="\t")
nodes["kind"] = nodes.kind.str.strip()
ents = sorted(set(nodes.id))
e2i = {e: i for i, e in enumerate(ents)}
kind_of = dict(zip(nodes.id, nodes.kind))
kinds = sorted(set(nodes.kind))
k2i = {k: i for i, k in enumerate(kinds)}
etype = np.array([k2i[kind_of[e]] for e in ents])
H, R, T, rels = [], [], [], {}
with gzip.open(f"{HET}/hetionet-v1.0-edges.sif.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        s, m, t = line.rstrip("\n").split("\t")
        H.append(e2i[s]); T.append(e2i[t]); R.append(rels.setdefault(m, len(rels)))
H, R, T = np.array(H), np.array(R), np.array(T)
nE, nR = len(ents), len(rels)
CtD, CpD = rels["CtD"], rels["CpD"]

if MODE == "pretrain":
    drop = (R == CtD) | (R == CpD)
else:
    part = json.load(open(PARTITION))
    hold = {(e2i["Compound::" + c], e2i["Disease::" + d]) for c, d in part["test_pos"] + part["val_pos"]}
    drop = np.array([(R[i] == CtD) and ((H[i], T[i]) in hold) for i in range(len(H))])
keep = ~drop
print(f"[{tag}] entities={nE} relations={nR} edges={len(H)} dropped={int(drop.sum())}", flush=True)

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


h, r, t = torch.tensor(H[keep]), torch.tensor(R[keep]), torch.tensor(T[keep])
h2, t2, r2 = torch.cat([h, t]), torch.cat([t, h]), torch.cat([r, r + nR])
nTr = h2.shape[0]
model = KGE(nE, 2 * nR, DIM, KIND)
opt = torch.optim.Adam(model.parameters(), lr=LR, weight_decay=L2)
bce = nn.BCEWithLogitsLoss()
comp_idx = np.array([e2i[e] for e in ents if kind_of[e] == "Compound"])
dis_idx = np.array([e2i[e] for e in ents if kind_of[e] == "Disease"])
t0 = time.time()
log = []


def save_embeddings(ep):
    state = {"entities": [ents[i] for i in np.r_[comp_idx, dis_idx]], "kind": KIND, "neg": NEG, "epoch": ep, "seed": SEED}
    with torch.no_grad():
        idx = torch.tensor(np.r_[comp_idx, dis_idx])
        state["er"] = model.er.weight[idx].half()
        if KIND in ("complex", "rotate"):
            state["ei"] = model.ei.weight[idx].half()
        state["rel_CtD"] = model.rr.weight[CtD].clone()
        if KIND == "complex":
            state["rel_CtD_im"] = model.ri.weight[CtD].clone()
    torch.save(state, f"{OUT}/{tag}_ep{ep}.pt")


for ep in range(1, EPOCHS + 1):
    perm = torch.randperm(nTr)
    tot = 0.0
    for b in range(0, nTr, BATCH):
        idx = perm[b: b + BATCH]
        ph, pr, pt = h2[idx], r2[idx], t2[idx]
        bs = ph.shape[0]
        if NEG == "uniform":
            neg_t = torch.randint(0, nE, (bs, NNEG))
        else:
            tt = etype_t[pt]
            u = torch.rand(bs, NNEG)
            neg_t = order_t[starts_t[tt].unsqueeze(1) + (u * counts_t[tt].unsqueeze(1)).long().clamp_(max=nE - 1)]
        pos_s = model.score(ph, pr, pt)
        neg_s = model.score(ph.unsqueeze(1).expand(bs, NNEG).reshape(-1), pr.unsqueeze(1).expand(bs, NNEG).reshape(-1),
                            neg_t.reshape(-1)).reshape(bs, NNEG)
        scores = torch.cat([pos_s.unsqueeze(1), neg_s], dim=1)
        labels = torch.zeros_like(scores); labels[:, 0] = 1.0
        loss = bce(scores, labels)
        opt.zero_grad(); loss.backward(); opt.step()
        tot += float(loss) * bs
    log.append({"epoch": ep, "loss": tot / nTr, "elapsed_s": time.time() - t0})
    print(f"  ep{ep} loss={tot / nTr:.4f} elapsed={time.time() - t0:.0f}s", flush=True)
    if MODE == "pretrain" and ep in SAVE_AT:
        save_embeddings(ep)
    if MODE == "e2e" and ep in SAVE_AT:
        # score the partition's evaluation grid with the CtD relation
        with torch.no_grad():
            grid = part["eval_grid"]
            cc = torch.tensor([e2i["Compound::" + c] for c, d in grid]); dd = torch.tensor([e2i["Disease::" + d] for c, d in grid])
            s = model.score(cc, torch.full((len(grid),), CtD, dtype=torch.long), dd).numpy()
            vg = part["val_grid"]
            vc = torch.tensor([e2i["Compound::" + c] for c, d in vg]); vd = torch.tensor([e2i["Disease::" + d] for c, d in vg])
            vs = model.score(vc, torch.full((len(vg),), CtD, dtype=torch.long), vd).numpy()
        np.savez_compressed(f"{OUT}/{tag}_ep{ep}_scores.npz", eval=s.astype(np.float32), val=vs.astype(np.float32))
json.dump({"tag": tag, "mode": MODE, "model": KIND, "neg": NEG, "seed": SEED, "epochs": EPOCHS, "save_at": SAVE_AT,
           "dropped_edges": int(drop.sum()), "train_triples_with_inverse": int(nTr), "log": log,
           "train_seconds": time.time() - t0}, open(done_flag, "w"), indent=1)
print("done", tag, flush=True)
