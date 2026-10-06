"""PrimeKG indication partitions: same construction as the Hetionet partitions (20% test,
10% validation; random-edge and compound-disjoint), three partition seeds. Grid compounds are
drugs with at least one indication; diseases are the 1,363 diseases with an indication;
off-label pairs play the role of Hetionet's palliative (CpD) pairs and are excluded."""
import json, os, numpy as np, pandas as pd
G = json.load(open("work/mapped_primekg/graph.json")); OUT = "work/partitions_primekg"; os.makedirs(OUT, exist_ok=True)
ctd = [tuple(p) for p in G["indication"]]; cpd = {tuple(p) for p in G["off_label"]} - set(ctd)   # 125 pairs are both; the indication label wins
labelled = G["diseases"]
def grid(compounds, exclude):
    return [(c, d) for c in sorted(compounds) for d in labelled if (c, d) not in exclude]
summ = []
for seed in [42, 1, 2]:
    rng = np.random.RandomState(seed)
    idx = rng.permutation(len(ctd)); nt = int(round(0.2 * len(ctd))); ti = set(idx[:nt].tolist())
    parts = {"random": ([ctd[i] for i in range(len(ctd)) if i not in ti], [ctd[i] for i in range(len(ctd)) if i in ti])}
    comps = sorted({c for c, d in ctd}); rng.shuffle(comps); tc = set(comps[: int(round(0.2 * len(comps)))])
    parts["compound"] = ([p for p in ctd if p[0] not in tc], [p for p in ctd if p[0] in tc])
    for task, (train, test) in parts.items():
        vr = np.random.default_rng(1000 + seed)
        if task == "random":
            vi = set(map(int, vr.choice(len(train), size=int(round(0.1 * len(train))), replace=False))); val = [train[i] for i in sorted(vi)]
        else:
            tcs = sorted({c for c, d in train}); vc = set(vr.choice(tcs, size=int(round(0.1 * len(tcs))), replace=False).tolist())
            val = [p for p in train if p[0] in vc]
        vs = set(val); fit = [p for p in train if p not in vs]
        eg = grid({c for c, d in test}, set(fit) | vs | cpd); vg = grid({c for c, d in val}, set(fit) | cpd)
        assert set(test) <= set(eg)
        ci = {c: i for i, c in enumerate(sorted({c for c, d in ctd}))}; di = {d: j for j, d in enumerate(labelled)}
        arr = {f"{k}_{a}": np.array([ci[c] if a == "c" else di[d] for c, d in v], np.int32)
               for k, v in [("fit_pos", fit), ("val_pos", val), ("test_pos", test), ("eval_grid", eg), ("val_grid", vg)] for a in "cd"}
        np.savez_compressed(f"{OUT}/{task}_p{seed}.npz", task=task, seed=seed, **arr)   # index form (compounds sorted, diseases as in graph.json)
        summ.append({"task": task, "seed": seed, "fit": len(fit), "val": len(val), "test": len(test), "test_compounds": len({c for c, d in test}),
                     "eval_grid": len(eg), "val_grid": len(vg)})
pd.DataFrame(summ).to_csv(f"{OUT}/partition_summary.tsv", sep="\t", index=False); print(pd.DataFrame(summ).to_string(index=False))
