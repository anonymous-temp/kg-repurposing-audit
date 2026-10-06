#!/usr/bin/env python
"""Hetionet treatment-label partitions for the revised analysis.

Two tasks x five partition seeds. Seed 42 reproduces the partitions of the original
pipeline exactly (one RandomState(42) used first for the random-edge split and then for
the compound-disjoint split); seeds 1-4 use a fresh RandomState each.

Each partition file holds:
  fit_pos   training positives used for fitting
  val_pos   validation positives (10% of non-test labels; compound-grouped for the
            compound-disjoint task, so validation is task-consistent)
  test_pos  held-out positives
  eval_grid every (test compound, labelled disease) pair except fitting/validation
            positives and Compound-palliates-Disease pairs, i.e. the full candidate grid
            rather than a sample of it
  val_grid  the same construction for validation compounds
"""
import gzip, json, os, sys
import numpy as np
import pandas as pd

HET = sys.argv[1] if len(sys.argv) > 1 else "work/data/hetionet"
OUT = sys.argv[2] if len(sys.argv) > 2 else "work/partitions"
os.makedirs(OUT, exist_ok=True)
TEST_FRAC, VAL_FRAC = 0.20, 0.10

ctd, cpd = [], set()
with gzip.open(f"{HET}/hetionet-v1.0-edges.sif.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        s, m, t = line.rstrip("\n").split("\t")
        if m == "CtD":
            ctd.append((s.replace("Compound::", ""), t.replace("Disease::", "")))
        elif m == "CpD":
            cpd.add((s.replace("Compound::", ""), t.replace("Disease::", "")))
labelled = sorted({d for c, d in ctd})


def split_random(pos, rng):
    idx = rng.permutation(len(pos))
    n_test = int(round(TEST_FRAC * len(pos)))
    test_i = set(idx[:n_test].tolist())
    return [pos[i] for i in range(len(pos)) if i not in test_i], [pos[i] for i in range(len(pos)) if i in test_i]


def split_compound(pos, rng):
    comps = sorted(set(c for c, d in pos))
    rng.shuffle(comps)
    test_c = set(comps[: int(round(TEST_FRAC * len(comps)))])
    return [(c, d) for c, d in pos if c not in test_c], [(c, d) for c, d in pos if c in test_c]


def grid(compounds, exclude):
    return [(c, d) for c in sorted(compounds) for d in labelled if (c, d) not in exclude]


summary = []
for seed in [42, 1, 2, 3, 4]:
    rng = np.random.RandomState(seed)
    parts = {}
    parts["random"] = split_random(ctd, rng)
    parts["compound"] = split_compound(ctd, rng)   # for seed 42 this follows split_random, as in the original code
    for task, (train, test) in parts.items():
        vr = np.random.default_rng(1000 + seed)
        if task == "random":
            vi = set(map(int, vr.choice(len(train), size=int(round(VAL_FRAC * len(train))), replace=False)))
            val = [train[i] for i in sorted(vi)]
        else:
            tc = sorted({c for c, d in train})
            vc = set(vr.choice(tc, size=int(round(VAL_FRAC * len(tc))), replace=False).tolist())
            val = [(c, d) for c, d in train if c in vc]
        vset = set(val)
        fit = [p for p in train if p not in vset]
        excl = set(fit) | vset | cpd
        eval_grid = grid({c for c, d in test}, excl)
        val_grid = grid({c for c, d in val}, set(fit) | cpd)
        assert set(test) <= set(eval_grid) and set(val) <= set(val_grid)
        rec = {"task": task, "seed": seed, "fit_pos": fit, "val_pos": val, "test_pos": test,
               "eval_grid": eval_grid, "val_grid": val_grid}
        json.dump(rec, open(f"{OUT}/{task}_p{seed}.json", "w"))
        summary.append({"task": task, "seed": seed, "fit": len(fit), "val": len(val), "test": len(test),
                        "test_compounds": len({c for c, d in test}), "eval_grid": len(eval_grid),
                        "test_compounds_without_fit_label": len({c for c, d in test} - {c for c, d in fit}),
                        "test_pos_from_compounds_without_fit_label": sum(c not in {x for x, _ in fit} for c, d in test)})
pd.DataFrame(summary).to_csv(f"{OUT}/partition_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(summary).to_string(index=False))
