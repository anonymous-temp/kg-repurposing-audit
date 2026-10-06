"""Merge the main Hetionet run (results/) with the role-policy run (results_role2/) for the paired bootstrap.
Scores of the five main policies are bit-identical between runs (checked on E2 earlier), so the merged
file is equivalent to a single run of all seven policies."""
import os, shutil, numpy as np, pandas as pd
A, B, O = "work/results", "work/results_role_only", "work/results_role"
os.makedirs(O, exist_ok=True)
shutil.copy(f"{A}/selection.json", O)
for stem in ["e1_random", "e1_compound", "e2"]:
    za, zb = np.load(f"{A}/{stem}_scores.npz"), np.load(f"{B}/{stem}_scores.npz")
    keep = {k: za[k] for k in za.files if "__w10__" in k or k.endswith("labels")}
    for k in zb.files:
        if k.endswith("labels"):
            assert np.array_equal(zb[k], keep[k]), k
        else:
            keep[k] = zb[k]
    np.savez_compressed(f"{O}/{stem}_scores.npz", **keep)
    ma, mb = pd.read_csv(f"{A}/{stem}_metrics.tsv", sep="\t"), pd.read_csv(f"{B}/{stem}_metrics.tsv", sep="\t")
    pd.concat([ma[ma.w_neg == 10], mb]).to_csv(f"{O}/{stem}_metrics.tsv", sep="\t", index=False)
    print(stem, len(keep))
