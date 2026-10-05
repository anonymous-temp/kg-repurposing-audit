#!/usr/bin/env python
"""Experiment B - data preparation.
Builds the Hetionet graph, the 755 CtD positives, Bemis-Murcko scaffolds (RDKit),
the candidate universe, three split schemes (naive-random / cold-start / scaffold)
and two negative-sampling schemes (random / degree-matched). Saves to experiments/expB/.
"""
import gzip, json, os, sys
import numpy as np
import pandas as pd
from collections import defaultdict

np.random.seed(42)
ROOT = os.environ["KG_AUDIT_WORKDIR"]
DATA = f"{ROOT}/data/hetionet"
OUT  = f"{ROOT}/outputs/hetionet/raw"
os.makedirs(OUT, exist_ok=True)

# ---- load nodes / edges ----
nodes = pd.read_csv(f"{DATA}/hetionet-v1.0-nodes.tsv", sep="\t")
compounds = nodes[nodes.kind == "Compound"].id.tolist()           # 1552
diseases  = nodes[nodes.kind == "Disease"].id.tolist()            # 137
print(f"compounds={len(compounds)} diseases={len(diseases)}")

edges = []
with gzip.open(f"{DATA}/hetionet-edges.sif.gz", "rt") as fh:
    header = fh.readline()
    for line in fh:
        s, m, t = line.rstrip("\n").split("\t")
        edges.append((s, m, t))
print(f"edges={len(edges)}")

ctd = [(s, t) for (s, m, t) in edges if m == "CtD"]
print(f"CtD positives={len(ctd)}")
pos_df = pd.DataFrame(ctd, columns=["compound", "disease"])
pos_df.to_csv(f"{OUT}/positives.tsv", sep="\t", index=False)

# ---- compound -> SMILES -> Bemis-Murcko scaffold ----
annot = pd.read_csv(f"{DATA}/drugcentral_compound_annot.tsv", sep="\t", dtype=str)
db2smiles = dict(zip(annot.drugbank_id, annot.smiles))

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

def scaffold_of(dbid):
    dbid_short = dbid.split("::")[-1]
    smi = db2smiles.get(dbid_short)
    if not smi or pd.isna(smi):
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
        s = Chem.MolToSmiles(scaf)
        return s if s else "ACYCLIC"
    except Exception:
        return None

comp_scaffold = {c: scaffold_of(c) for c in compounds}
n_mapped = sum(v is not None for v in comp_scaffold.values())
print(f"compounds with scaffold: {n_mapped}/{len(compounds)}")
# coverage among positive compounds
pos_comps = pos_df.compound.unique()
n_pos_mapped = sum(comp_scaffold.get(c) is not None for c in pos_comps)
print(f"positive compounds with scaffold: {n_pos_mapped}/{len(pos_comps)}")

pd.DataFrame([(c, comp_scaffold[c]) for c in compounds],
             columns=["compound", "scaffold"]).to_csv(
    f"{OUT}/compound_scaffold.tsv", sep="\t", index=False)

# ---- candidate universe & degrees ----
pos_set = set(map(tuple, ctd))
comp_deg = pos_df.compound.value_counts().to_dict()   # CtD-degree among positives
dis_deg  = pos_df.disease.value_counts().to_dict()

# ---- SPLIT SCHEMES (test fraction ~0.20 of positives) ----
rng = np.random.RandomState(42)
TEST_FRAC = 0.20

def split_random(pos):
    idx = rng.permutation(len(pos))
    n_test = int(round(TEST_FRAC * len(pos)))
    test_i = set(idx[:n_test].tolist())
    train = [pos[i] for i in range(len(pos)) if i not in test_i]
    test  = [pos[i] for i in range(len(pos)) if i in test_i]
    return train, test

def split_coldstart(pos):
    comps = sorted(set(c for c, d in pos))
    rng.shuffle(comps)
    n_test_c = int(round(TEST_FRAC * len(comps)))
    test_c = set(comps[:n_test_c])
    train = [(c, d) for c, d in pos if c not in test_c]
    test  = [(c, d) for c, d in pos if c in test_c]
    return train, test, test_c

def split_scaffold(pos):
    # group positive compounds by scaffold; hold out whole scaffolds
    scafs = defaultdict(list)
    for c in set(c for c, d in pos):
        scafs[comp_scaffold.get(c) or f"__none__{c}"].append(c)
    scaf_keys = sorted(scafs.keys())
    rng.shuffle(scaf_keys)
    target = TEST_FRAC * len(set(c for c, d in pos))
    test_c, cum = set(), 0
    for k in scaf_keys:
        if cum >= target:
            break
        for c in scafs[k]:
            test_c.add(c); cum += 1
    train = [(c, d) for c, d in pos if c not in test_c]
    test  = [(c, d) for c, d in pos if c in test_c]
    return train, test, test_c

splits = {}
tr, te = split_random(ctd)
splits["random"] = {"train": tr, "test": te, "test_comps": sorted(set(c for c, d in te))}
tr, te, tc = split_coldstart(ctd)
splits["coldstart"] = {"train": tr, "test": te, "test_comps": sorted(tc)}
tr, te, tc = split_scaffold(ctd)
splits["scaffold"] = {"train": tr, "test": te, "test_comps": sorted(tc)}

for k, v in splits.items():
    print(f"{k}: train_pos={len(v['train'])} test_pos={len(v['test'])} "
          f"test_comps={len(v['test_comps'])}")

# ---- NEGATIVE SAMPLING ----
# For each scheme, build eval negatives for the TEST positives.
# random: uniform over non-CtD pairs whose compound is a test compound (so the
#         classification task concerns the same compounds being scored).
# degree-matched: negatives whose disease degree distribution matches the test positives'.
all_pairs = set((c, d) for c in compounds for d in diseases)
neg_pool = list(all_pairs - pos_set)
neg_pool_arr = np.array(neg_pool, dtype=object)

# Generation pool: 30 negatives per positive (this larger pool also serves the closed-world
# ranking analysis). The naive-vs-leakage AUPRC contrast is then evaluated at a COMMON 20:1
# prevalence (chance AUPRC 0.0476) by subsampling these negatives in aggregate_benchmark.py, so
# Hetionet, PrimeKG (ratio=20) and DRKG (20x) are all scored at one prevalence and one chance floor.
NEG_RATIO = 30

def sample_random_negs(test_pos, ratio=NEG_RATIO):
    test_comps = set(c for c, d in test_pos)
    pool = [(c, d) for (c, d) in neg_pool if c in test_comps]
    k = min(len(pool), ratio * len(test_pos))
    sel = rng.choice(len(pool), size=k, replace=False)
    return [pool[i] for i in sel]

def sample_degree_matched_negs(test_pos, ratio=NEG_RATIO):
    # match disease frequency: sample negatives with disease drawn from the test
    # positives' disease multiset, compound from test compounds -> removes "rare disease" shortcut
    test_comps = sorted(set(c for c, d in test_pos))
    test_dis   = [d for c, d in test_pos]
    negs, tries = [], 0
    want = ratio * len(test_pos)
    seen = set()
    while len(negs) < want and tries < want * 50:
        d = test_dis[rng.randint(len(test_dis))]
        c = test_comps[rng.randint(len(test_comps))]
        tries += 1
        if (c, d) in pos_set or (c, d) in seen:
            continue
        seen.add((c, d)); negs.append((c, d))
    return negs

eval_sets = {}
for scheme, v in splits.items():
    tp = [tuple(x) for x in v["test"]]
    eval_sets[scheme] = {
        "test_pos": tp,
        "neg_random": sample_random_negs(tp),
        "neg_degmatch": sample_degree_matched_negs(tp),
        "train_pos": [tuple(x) for x in v["train"]],
    }
    print(f"{scheme}: eval neg_random={len(eval_sets[scheme]['neg_random'])} "
          f"neg_degmatch={len(eval_sets[scheme]['neg_degmatch'])}")

# training negatives (for supervised baseline): random non-CtD among all compounds
def sample_train_negs(train_pos, ratio=NEG_RATIO):
    train_comps = set(c for c, d in train_pos)
    pool = [(c, d) for (c, d) in neg_pool if c in train_comps]
    k = min(len(pool), ratio * len(train_pos))
    sel = rng.choice(len(pool), size=k, replace=False)
    return [pool[i] for i in sel]

for scheme in splits:
    eval_sets[scheme]["train_neg"] = sample_train_negs(eval_sets[scheme]["train_pos"])

with open(f"{OUT}/splits.json", "w") as fh:
    json.dump({k: {kk: [list(x) for x in vv] if isinstance(vv, list) else vv
                   for kk, vv in v.items()} for k, v in eval_sets.items()}, fh)
print("saved splits.json")
print("PREP_DONE")
