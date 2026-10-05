#!/usr/bin/env python
"""Two-layer relational graph encoder with a DistMult treatment decoder.
Reads local, licensed Hetionet inputs; writes per-pair scores and validation-selected ranks.
Auxiliary graph identities remain available under treatment-label holdout.
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
DEV = torch.device(os.environ.get("DEV", "cpu"))
torch.set_num_threads(int(os.environ.get("THREADS", max(1, os.cpu_count() - 1))))

# encoder/decoder hyper-params (DIM matches KGE; RGCN-specific knobs kept modest for CPU budget)
DIM = int(os.environ.get("DIM", 100))
NBASES = int(os.environ.get("NBASES", 10))
NLAYERS = int(os.environ.get("NLAYERS", 2))
NNEG = 4
LR = float(os.environ.get("LR", 3e-3))
L2 = 1e-6
DROPOUT = 0.0
# scoring is full-batch over the held-out (compound,disease) pairs; training samples positive
# CtD edges in mini-batches but the encoder runs full-graph once per forward pass.
BATCH = int(os.environ.get("BATCH", 4096))
# validation-based early stopping (same selection budget shape as expB_kge.py; fewer max epochs
# because each R-GCN epoch does a full-graph propagation and is far more expensive than a KGE step)
MAX_EPOCHS = int(os.environ.get("EPOCHS", os.environ.get("MAX_EPOCHS", 30)))
MIN_EPOCHS = int(os.environ.get("MIN_EPOCHS", 4))
EVAL_EVERY = 2
PATIENCE = int(os.environ.get("PATIENCE", 3))
VAL_FRAC = 0.10
SUFFIX = os.environ.get("SUFFIX", "")
RUN_SPLITS = os.environ.get("SPLITS", "random,coldstart,scaffold").split(",")

# ---------- load graph (identical to expB_kge.py) ----------
edges = []
with gzip.open(f"{DATA}/hetionet-edges.sif.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        s, m, t = line.rstrip("\n").split("\t")
        edges.append((s, m, t))
edges = np.array(edges, dtype=object)
_nodes = pd.read_csv(f"{DATA}/hetionet-v1.0-nodes.tsv", sep="\t")
ents = sorted(set(_nodes.id) | set(edges[:, 0]) | set(edges[:, 2]))
rels = sorted(set(edges[:, 1]))
e2i = {e: i for i, e in enumerate(ents)}
r2i = {r: i for i, r in enumerate(rels)}
CtD = "CtD"
nE, nR = len(ents), len(rels)
print(
    f"entities={nE} relations={nR} dev={DEV} dim={DIM} bases={NBASES} layers={NLAYERS} max_epochs={MAX_EPOCHS}",
    flush=True,
)

H = np.array([e2i[s] for s in edges[:, 0]], dtype=np.int64)
R = np.array([r2i[m] for m in edges[:, 1]], dtype=np.int64)
T = np.array([e2i[t] for t in edges[:, 2]], dtype=np.int64)
ctd_mask = edges[:, 1] == CtD
CAND = np.array(sorted({e2i[edges[i, 2]] for i in np.where(ctd_mask)[0]}), dtype=np.int64)
CtD_IDX = r2i[CtD]

splits = json.load(open(f"{OUT}/splits.json"))


def to_idx_pairs(pairs):
    return [(e2i[c], e2i[d]) for c, d in pairs]


# ---------- R-GCN encoder ----------
# Build, ONCE per split, the per-relation normalised sparse adjacency used for message passing.
# We add inverse relations (a separate relation block) and a self-loop, mirroring the KGE setup
# (which also adds inverse relations). Total relation blocks = 2*nR (fwd+inv) + 1 (self-loop).
NREL_MP = 2 * nR + 1  # message-passing relation blocks
SELF_LOOP = 2 * nR  # index of the self-loop block


def build_adjacency(keep_mask):
    """Return a list of (indices, values, size) sparse-COO pieces, one per message-passing
    relation block, with symmetric per-relation degree normalisation 1/sqrt(d_i d_j).
    keep_mask: boolean over `edges` selecting which edges enter message passing."""
    h = H[keep_mask]
    r = R[keep_mask]
    t = T[keep_mask]
    blocks = []
    # forward block rb=r ; inverse block rb=nR+r
    for direction, (src, dst, base) in enumerate([(h, t, 0), (t, h, nR)]):
        for rb in range(nR):
            sel = r == rb
            if not sel.any():
                blocks.append(None)
                continue
            i = src[sel]
            j = dst[sel]  # message i -> j for this relation
            deg = np.zeros(nE, dtype=np.float64)
            np.add.at(deg, j, 1.0)  # in-degree per relation (normaliser)
            deg[deg == 0] = 1.0
            val = 1.0 / deg[j]  # mean-aggregation over relation-r neighbours
            idx = torch.tensor(np.stack([j, i]), dtype=torch.long)  # row=target j, col=source i
            v = torch.tensor(val, dtype=torch.float32)
            a = torch.sparse_coo_tensor(idx, v, (nE, nE)).coalesce()
            blocks.append(a)
    # self-loop (identity)
    diag = torch.arange(nE)
    a = torch.sparse_coo_tensor(torch.stack([diag, diag]), torch.ones(nE), (nE, nE)).coalesce()
    blocks.append(a)
    return blocks


class RGCNLayer(nn.Module):
    """Basis-decomposed R-GCN layer: W_rel = sum_b a_{rel,b} V_b  (Schlichtkrull eq. 3)."""

    def __init__(self, in_dim, out_dim, nrel, nbases):
        super().__init__()
        self.nrel, self.nbases = nrel, nbases
        self.V = nn.Parameter(torch.empty(nbases, in_dim, out_dim))
        self.comp = nn.Parameter(torch.empty(nrel, nbases))  # mixing coefficients
        nn.init.xavier_uniform_(self.V)
        nn.init.xavier_uniform_(self.comp)

    def forward(self, x, blocks):
        # W per relation = comp @ V  ->  (nrel, in, out)
        W = torch.einsum("rb,bio->rio", self.comp, self.V)
        out = torch.zeros(x.shape[0], W.shape[2], device=x.device)
        for rb, A in enumerate(blocks):
            if A is None:
                continue
            msg = torch.sparse.mm(A, x)  # aggregate neighbours under relation rb
            out = out + msg @ W[rb]  # then relation-specific transform
        return out


class RGCN(nn.Module):
    def __init__(self, nE, nrel, dim, nbases, nlayers):
        super().__init__()
        self.emb = nn.Embedding(nE, dim)
        nn.init.xavier_uniform_(self.emb.weight)
        self.layers = nn.ModuleList([RGCNLayer(dim, dim, nrel, nbases) for _ in range(nlayers)])
        self.act = nn.ReLU()
        self.drop = nn.Dropout(DROPOUT)
        # DistMult decoder relation vector for CtD (the only relation we score)
        self.rel_ctd = nn.Parameter(torch.empty(dim))
        nn.init.xavier_uniform_(self.rel_ctd.view(1, -1))

    def encode(self, blocks):
        x = self.emb.weight
        for i, layer in enumerate(self.layers):
            x = layer(x, blocks)
            if i < len(self.layers) - 1:
                x = self.drop(self.act(x))
        return x

    def decode(self, x, h_idx, t_idx):
        return (x[h_idx] * self.rel_ctd * x[t_idx]).sum(-1)  # DistMult on CtD


@torch.no_grad()
def rank_metrics(emb, queries, known):
    """Full-candidate FILTERED ranking, identical semantics to expB_kge.py."""
    if not queries:
        return {"MRR": float("nan"), "Hits@1": float("nan"), "Hits@3": float("nan"), "Hits@10": float("nan"), "n": 0}
    ct = torch.tensor(CAND, device=DEV)
    cand_t = emb[ct]  # (|CAND|, dim)
    pos = {int(d): j for j, d in enumerate(CAND)}
    rr = []
    hits = {1: 0, 3: 0, 10: 0}
    for c, dt in queries:
        s = (emb[c] * model.rel_ctd * cand_t).sum(-1).cpu().numpy().copy()
        for d in known.get(c, ()):
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


def train_model(blocks, pos_pairs, val_queries, known):
    """pos_pairs: list of (c_idx,d_idx) training CtD positives that ARE in the graph.
    Encoder runs full-graph each forward; decoder scored on sampled CtD pos/neg.
    Early-stops on validation MRR (same selection rule as expB_kge.py)."""
    model_ = RGCN(nE, NREL_MP, DIM, NBASES, NLAYERS).to(DEV)
    opt = torch.optim.Adam(model_.parameters(), lr=LR, weight_decay=L2)
    bce = nn.BCEWithLogitsLoss()
    P = torch.tensor(np.array(pos_pairs, dtype=np.int64), device=DEV)
    nP = P.shape[0]
    best_mrr, best_state, wait, best_ep = -1.0, None, 0, 0
    global model
    for ep in range(MAX_EPOCHS):
        model_.train()
        perm = torch.randperm(nP, device=DEV)
        tot = 0.0
        for b in range(0, nP, BATCH):
            emb = model_.encode(blocks)  # fresh differentiable encoder state for each optimizer step
            idx = perm[b : b + BATCH]
            ph = P[idx, 0]
            pt = P[idx, 1]
            bs = ph.shape[0]
            neg_t = torch.randint(0, nE, (bs, NNEG), device=DEV)  # corrupt the tail (disease side)
            pos_s = model_.decode(emb, ph, pt)
            ph_e = ph.unsqueeze(1).expand(bs, NNEG).reshape(-1)
            neg_s = model_.decode(emb, ph_e, neg_t.reshape(-1)).reshape(bs, NNEG)
            scores = torch.cat([pos_s.unsqueeze(1), neg_s], dim=1)
            labels = torch.zeros_like(scores)
            labels[:, 0] = 1.0
            loss = bce(scores, labels)
            opt.zero_grad()
            loss.backward()
            opt.step()
            tot += loss.item() * bs
        if ep + 1 >= MIN_EPOCHS and (ep % EVAL_EVERY == 0 or ep == MAX_EPOCHS - 1):
            model_.eval()
            with torch.no_grad():
                emb_eval = model_.encode(blocks)
            model = model_
            vm = rank_metrics(emb_eval, val_queries, known)["MRR"]
            if vm > best_mrr + 1e-4:
                best_mrr, best_state, wait, best_ep = (
                    vm,
                    {k: v.detach().clone() for k, v in model_.state_dict().items()},
                    0,
                    ep,
                )
            else:
                wait += 1
            if wait >= PATIENCE:
                break
    if best_state is not None:
        model_.load_state_dict(best_state)
    model = model_
    print(f"   rgcn: stopped ep{ep} best_ep{best_ep} val_MRR={best_mrr:.3f}", flush=True)
    with torch.no_grad():
        return model_.encode(blocks)


@torch.no_grad()
def score_pairs(emb, pairs):
    if not pairs:
        return np.array([])
    cd = np.array(pairs, dtype=np.int64)
    h = torch.tensor(cd[:, 0])
    t = torch.tensor(cd[:, 1])
    return model.decode(emb, h, t).cpu().numpy()


model = None
for split in RUN_SPLITS:
    sset = splits[split]
    test_set = set(map(tuple, sset["test_pos"]))
    all_ctd = [(edges[i, 0], edges[i, 2]) for i in np.where(ctd_mask)[0]]
    train_ctd = [p for p in all_ctd if p not in test_set]
    # validation hold-out for model selection (SAME seeded scheme as expB_kge.py)
    vrng = np.random.default_rng(SEED * 100 + len(split))
    nval = max(10, int(VAL_FRAC * len(train_ctd)))
    vsel = set(map(int, vrng.choice(len(train_ctd), size=nval, replace=False)))
    val_pairs = [train_ctd[i] for i in vsel]
    graph_ctd = [train_ctd[i] for i in range(len(train_ctd)) if i not in vsel]
    # edges to DROP from the message-passing graph = test CtD edges + validation CtD edges
    hold = test_set | set(val_pairs)
    keep = np.ones(len(edges), dtype=bool)
    for idx in np.where(ctd_mask)[0]:
        if (edges[idx, 0], edges[idx, 2]) in hold:
            keep[idx] = False
    print(
        f"[{split}] drop {(~keep).sum()} CtD edges (test {len(test_set)} + val {len(val_pairs)}) from message passing",
        flush=True,
    )
    blocks = build_adjacency(keep)
    # training positives for the decoder = CtD edges actually in the graph
    pos_pairs = [(e2i[c], e2i[d]) for c, d in graph_ctd]
    known = {}
    for c, d in graph_ctd:
        known.setdefault(e2i[c], set()).add(e2i[d])
    val_q = [(e2i[c], e2i[d]) for c, d in val_pairs]
    test_q = [(e2i[c], e2i[d]) for c, d in sset["test_pos"]]
    pack = {
        k: to_idx_pairs([tuple(x) for x in sset[k]])
        for k in ["test_pos", "neg_random", "neg_degmatch", "train_pos", "train_neg"]
    }
    t0 = time.time()
    emb = train_model(blocks, pos_pairs, val_q, known)
    rows = []
    for grp, prs in pack.items():
        sc = score_pairs(emb, prs)
        for (c, d), s in zip(sset[grp], sc):
            rows.append((grp, c, d, float(s)))
    pd.DataFrame(rows, columns=["group", "compound", "disease", "score"]).to_csv(
        f"{OUT}/scores_rgcn_{split}{SUFFIX}.tsv", sep="\t", index=False
    )
    rm = rank_metrics(emb, test_q, known)
    json.dump(
        {"model": "rgcn", "split": split, "seed": SEED, **rm},
        open(f"{OUT}/rank_rgcn_{split}{SUFFIX}.json", "w"),
        indent=2,
    )
    print(f"[{split}/rgcn] {time.time() - t0:.0f}s -> MRR {rm['MRR']:.3f} Hits@10 {rm['Hits@10']:.3f}", flush=True)
print("RGCN_DONE")
