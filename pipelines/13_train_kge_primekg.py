#!/usr/bin/env python
"""DistMult embeddings of PrimeKG without any drug-disease edge (indication, off-label use,
contraindication), with type-constrained tail corruption. Same recipe as the Hetionet
pretraining (dim 100, BCE, Adam 3e-3, weight decay 1e-6, 4 corruptions, batch 32768,
forward and inverse relations). Each undirected PrimeKG edge is used once in each direction.
Idempotent: skipped when the done flag exists."""
import json, os, sys, time
import numpy as np, pandas as pd, torch, torch.nn as nn
PK, OUT = "work/data/primekg", "work/kge_primekg"
SEED, EPOCHS = 1, int(os.environ.get("EPOCHS", 10))
SAVE_AT = [int(v) for v in os.environ.get("SAVE_AT", "5,10").split(",")]
DIM, BATCH, NNEG, LR, L2 = 100, 32768, 4, 3e-3, 1e-6
tag = "pretrain_distmult_typed_s1"; flag = f"{OUT}/{tag}.done.json"
if os.path.exists(flag):
    sys.exit("skip (done)")
torch.manual_seed(SEED); np.random.seed(SEED); torch.set_num_threads(int(os.environ.get("THREADS", 2)))
n = pd.read_csv(f"{PK}/nodes.csv", dtype={"node_id": str})
E = np.load(f"{PK}/edges.npz"); rels = list(E["rels"])
x, y, r = E["x"].astype(np.int64), E["y"].astype(np.int64), E["r"].astype(np.int64)
drop_rel = [rels.index(k) for k in ("indication", "off-label use", "contraindication")]
keep = (x < y) & ~np.isin(r, drop_rel)          # one direction per undirected edge; no drug-disease edges
H, T, R = torch.tensor(x[keep]), torch.tensor(y[keep]), torch.tensor(r[keep])
nE, nR = len(n), len(rels)
kinds = sorted(n.node_type.unique()); etype = n.node_type.map({k: i for i, k in enumerate(kinds)}).values
order = np.argsort(etype, kind="stable"); counts = np.bincount(etype, minlength=len(kinds)); starts = np.r_[0, np.cumsum(counts)[:-1]]
order_t, counts_t, starts_t, etype_t = map(torch.tensor, (order, counts, starts, etype))
print(f"[{tag}] entities={nE} relations={nR} triples={int(keep.sum())} dropped_drug_disease={int(np.isin(r, drop_rel).sum() // 2)}", flush=True)
h2, t2, r2 = torch.cat([H, T]), torch.cat([T, H]), torch.cat([R, R + nR]); del H, T, R, x, y, r
nTr = h2.shape[0]
er = nn.Embedding(nE, DIM); rr = nn.Embedding(2 * nR, DIM); nn.init.xavier_uniform_(er.weight); nn.init.xavier_uniform_(rr.weight)
opt = torch.optim.Adam(list(er.parameters()) + list(rr.parameters()), lr=LR, weight_decay=L2); bce = nn.BCEWithLogitsLoss()
drug_idx = n.index[n.node_type == "drug"].values; dis_idx = n.index[n.node_type == "disease"].values
ent_names = ["Compound::" + v for v in n.node_id.values[drug_idx]] + ["Disease::" + v for v in n.node_id.values[dis_idx]]
t0 = time.time(); log = []
for ep in range(1, EPOCHS + 1):
    perm = torch.randperm(nTr); tot = 0.0
    for b in range(0, nTr, BATCH):
        idx = perm[b:b + BATCH]; ph, pr, pt = h2[idx], r2[idx], t2[idx]; bs = ph.shape[0]
        tt = etype_t[pt]; u = torch.rand(bs, NNEG)
        neg_t = order_t[starts_t[tt].unsqueeze(1) + (u * counts_t[tt].unsqueeze(1)).long().clamp_(max=nE - 1)]
        hv, rv = er(ph), rr(pr)
        pos_s = (hv * rv * er(pt)).sum(-1)
        neg_s = ((hv * rv).unsqueeze(1) * er(neg_t)).sum(-1)
        sc = torch.cat([pos_s.unsqueeze(1), neg_s], 1); lab = torch.zeros_like(sc); lab[:, 0] = 1
        loss = bce(sc, lab); opt.zero_grad(); loss.backward(); opt.step(); tot += float(loss) * bs
    log.append({"epoch": ep, "loss": tot / nTr, "elapsed_s": time.time() - t0})
    print(f"  ep{ep} loss={tot / nTr:.4f} elapsed={time.time() - t0:.0f}s", flush=True)
    if ep in SAVE_AT:
        with torch.no_grad():
            st = {"entities": ent_names, "er": er.weight[torch.tensor(np.r_[drug_idx, dis_idx])].half(), "epoch": ep, "seed": SEED}
        torch.save(st, f"{OUT}/{tag}_ep{ep}.pt")
json.dump({"tag": tag, "epochs": EPOCHS, "triples": int(nTr // 2), "log": log, "train_seconds": time.time() - t0}, open(flag, "w"), indent=1)
print("done", flush=True)
