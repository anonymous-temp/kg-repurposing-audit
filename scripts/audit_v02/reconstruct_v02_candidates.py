#!/usr/bin/env python
"""Reconstruction of the v0.2.0 Hetionet splits and candidate samplers (Supplementary S8).

Mirrors pipelines/hetionet_prepare.py of kg-repurposing-audit v0.2.0 (commit fea88c7)
for the random-edge and compound-disjoint splits (the scaffold split needs the
authors' DrugCentral mapping and is skipped). The random/cold-start *splits* use the
same RandomState(42) call order as the authors and should be identical; candidate
pools are redrawn here (the authors' RNG state after the scaffold shuffle is unknown),
so pool-level numbers are distributionally comparable, not row-identical.
"""
import gzip, json, sys, os
from collections import Counter
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from kg_audit.metrics import fit_degree_reference

DATA = os.environ.get("HET", "work/data/hetionet")
OUT = os.environ.get("OUT", "work/audit_v02")
os.makedirs(OUT, exist_ok=True)

nodes = pd.read_csv(f"{DATA}/hetionet-v1.0-nodes.tsv", sep="\t")
compounds = nodes[nodes.kind == "Compound"].id.tolist()
diseases = nodes[nodes.kind == "Disease"].id.tolist()
kind_counts = nodes.kind.value_counts().to_dict()
print(f"nodes={len(nodes)} compounds={len(compounds)} diseases={len(diseases)}")

ctd, cpd, meta = [], [], Counter()
with gzip.open(f"{DATA}/hetionet-v1.0-edges.sif.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        s, m, t = line.rstrip("\n").split("\t")
        meta[m] += 1
        if m == "CtD":
            ctd.append((s, t))
        elif m == "CpD":
            cpd.append((s, t))
print(f"edges={sum(meta.values())} metaedges={len(meta)} CtD={len(ctd)} CpD={len(cpd)}")
pos_set, cpd_set = set(ctd), set(cpd)
print("CtD compounds/diseases:", len({c for c, d in ctd}), len({d for c, d in ctd}))
print("CpD pairs that are also CtD:", len(pos_set & cpd_set))
print("diseases share of all nodes: %.4f ; compounds share: %.4f" % (len(diseases) / len(nodes), len(compounds) / len(nodes)))

# ---- splits: identical call order to the authors' script ----
rng = np.random.RandomState(42)
TEST_FRAC = 0.20


def split_random(pos):
    idx = rng.permutation(len(pos))
    n_test = int(round(TEST_FRAC * len(pos)))
    test_i = set(idx[:n_test].tolist())
    return [pos[i] for i in range(len(pos)) if i not in test_i], [pos[i] for i in range(len(pos)) if i in test_i]


def split_coldstart(pos):
    comps = sorted(set(c for c, d in pos))
    rng.shuffle(comps)
    n_test_c = int(round(TEST_FRAC * len(comps)))
    test_c = set(comps[:n_test_c])
    return [(c, d) for c, d in pos if c not in test_c], [(c, d) for c, d in pos if c in test_c], test_c


splits = {}
tr, te = split_random(ctd)
splits["random"] = {"train": tr, "test": te}
tr, te, tc = split_coldstart(ctd)
splits["coldstart"] = {"train": tr, "test": te}
for k, v in splits.items():
    print(f"[split {k}] train_pos={len(v['train'])} test_pos={len(v['test'])} "
          f"test_compounds={len({c for c, d in v['test']})} test_diseases={len({d for c, d in v['test']})}")

# ---- candidate samplers (same logic as v0.2.0; RNG state after the scaffold shuffle is unknown, so pools are redrawn) ----
all_pairs = set((c, d) for c in compounds for d in diseases)
neg_pool = sorted(all_pairs - pos_set)
NEG_RATIO = 30


def sample_random_negs(test_pos, r):
    test_comps = set(c for c, d in test_pos)
    pool = [(c, d) for (c, d) in neg_pool if c in test_comps]
    k = min(len(pool), NEG_RATIO * len(test_pos))
    sel = r.choice(len(pool), size=k, replace=False)
    return [pool[i] for i in sel], len(pool)


def sample_degree_matched_negs(test_pos, r):
    test_comps = sorted(set(c for c, d in test_pos))
    test_dis = [d for c, d in test_pos]
    negs, tries, want, seen = [], 0, NEG_RATIO * len(test_pos), set()
    while len(negs) < want and tries < want * 50:
        d = test_dis[r.randint(len(test_dis))]
        c = test_comps[r.randint(len(test_comps))]
        tries += 1
        if (c, d) in pos_set or (c, d) in seen:
            continue
        seen.add((c, d))
        negs.append((c, d))
    return negs, tries


def sample_train_negs(train_pos, r):
    train_comps = set(c for c, d in train_pos)
    pool = [(c, d) for (c, d) in neg_pool if c in train_comps]
    k = min(len(pool), NEG_RATIO * len(train_pos))
    sel = r.choice(len(pool), size=k, replace=False)
    return [pool[i] for i in sel]


full_ddeg = Counter(d for c, d in ctd)
report = {}
export = {}
for scheme, v in splits.items():
    test, train = v["test"], v["train"]
    r = np.random.RandomState(1000 + len(scheme))
    neg_u, pool_size = sample_random_negs(test, r)
    neg_d, tries = sample_degree_matched_negs(test, r)
    trneg = sample_train_negs(train, r)
    export[scheme] = {"test_pos": test, "neg_random": neg_u, "neg_degmatch": neg_d, "train_pos": train, "train_neg": trneg}
    test_comps = sorted({c for c, d in test})
    tdis = Counter(d for c, d in test)
    print(f"\n===== {scheme}: test_pos={len(test)} | uniform pool universe={pool_size} pairs "
          f"({len(test_comps)} test compounds x {len(diseases)} diseases) | "
          f"pools: uniform={len(neg_u)} disease-frequency={len(neg_d)} (target {NEG_RATIO * len(test)}; tries={tries})")
    # saturation: available unlabelled pairs per disease vs demand at 20:1 and 30:1
    sat20 = sum(1 for d, k in tdis.items() if 20 * k > len(test_comps) - sum(1 for c in test_comps if (c, d) in pos_set))
    sat30 = sum(1 for d, k in tdis.items() if 30 * k > len(test_comps) - sum(1 for c in test_comps if (c, d) in pos_set))
    pos_in_sat20 = sum(k for d, k in tdis.items() if 20 * k > len(test_comps) - sum(1 for c in test_comps if (c, d) in pos_set))
    print(f"  diseases among test positives: {len(tdis)}; diseases whose 20:1 demand exceeds ALL available unlabelled "
          f"(test compound, disease) pairs: {sat20} (30:1: {sat30}); these hold {pos_in_sat20}/{len(test)} test positives")
    print("  top test-positive diseases (count):", tdis.most_common(6))
    res = {}
    for nk, pool in [("neg_random", neg_u), ("neg_degmatch", neg_d)]:
        sel = np.random.default_rng(20261003 + len(scheme) + (nk == "neg_degmatch")).choice(len(pool), 20 * len(test), replace=False)
        ng = [pool[i] for i in sel]
        pairs = test + ng
        y = np.r_[np.ones(len(test)), np.zeros(len(ng))]
        ndis = Counter(d for c, d in ng)
        zero_deg_share = np.mean([full_ddeg[d] == 0 for c, d in ng])
        cpd_in_neg = sum(p in cpd_set for p in ng)
        cpd_in_pool = sum(p in cpd_set for p in pool)
        # distribution mismatch between positive and unlabelled disease identity
        alld = sorted(set(tdis) | set(ndis))
        pp = np.array([tdis[d] / len(test) for d in alld])
        qq = np.array([ndis[d] / len(ng) for d in alld])
        tv = 0.5 * np.abs(pp - qq).sum()
        top = [d for d, _ in tdis.most_common(5)]
        aps, aus, daps, daus = [], [], [], []
        for seed in (1, 2, 3):
            vrng = np.random.default_rng(seed * 100 + len(scheme))
            vsel = set(map(int, vrng.choice(len(train), max(10, int(0.1 * len(train))), replace=False)))
            fit = [p for i, p in enumerate(train) if i not in vsel]
            evalset = set(pairs)
            tn = [p for p in trneg if p not in evalset]
            pred, cdeg, ddeg = fit_degree_reference(fit, tn, seed)
            s = pred(pairs)
            aps.append(average_precision_score(y, s)); aus.append(roc_auc_score(y, s))
            ds = np.array([np.log1p(ddeg[d]) for c, d in pairs])
            daps.append(average_precision_score(y, ds)); daus.append(roc_auc_score(y, ds))
        res[nk] = dict(n_neg=len(ng), zero_treatment_degree_disease_share=float(zero_deg_share), cpd_in_candidates=cpd_in_neg,
                       cpd_in_pool=cpd_in_pool, tv_distance_disease_identity=float(tv),
                       degree_AP=float(np.mean(aps)), degree_AUROC=float(np.mean(aus)),
                       disease_only_AP=float(np.mean(daps)), disease_only_AUROC=float(np.mean(daus)))
        print(f"  [{nk}] negatives={len(ng)} | share with a disease that has NO treatment label anywhere: {zero_deg_share:.3f} | "
              f"CpD (palliates) pairs scored as unlabelled: {cpd_in_neg} (pool {cpd_in_pool})")
        print(f"      disease-identity total-variation distance positives vs unlabelled: {tv:.3f}; "
              f"top-5 positive diseases hold {sum(tdis[d] for d in top) / len(test):.3f} of positives vs {sum(ndis[d] for d in top) / len(ng):.3f} of unlabelled")
        print(f"      degree reference AP={np.mean(aps):.3f} AUROC={np.mean(aus):.3f} | disease-degree only AP={np.mean(daps):.3f} AUROC={np.mean(daus):.3f}")
    report[scheme] = res

json.dump({k: {kk: [list(x) for x in vv] for kk, vv in v.items()} for k, v in export.items()}, open(f"{OUT}/splits.json", "w"))
json.dump(report, open(f"{OUT}/sampler_report.json", "w"), indent=1)
print("\nsaved", OUT)
