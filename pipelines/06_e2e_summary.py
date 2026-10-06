#!/usr/bin/env python
"""Summarise the end-to-end KGE ablation (all-entity vs type-constrained tail corruption).

For each run, the checkpoint (saved every two epochs) with the highest validation AP on the
partition's validation grid is selected; test metrics are computed on the full evaluation grid.
"""
import glob, json, os, re, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from kg_audit import writeback as L

KGE = os.environ.get("KGE", "work/kge")
PART = os.environ.get("PART", "work/partitions")
OUT = os.environ.get("OUT", "work/results")
rows = []
for done in sorted(glob.glob(f"{KGE}/e2e_*.done.json")):
    meta = json.load(open(done))
    tag = meta["tag"]
    pname = re.search(r"_(random|compound)_p(\d+)$", tag)
    p = json.load(open(f"{PART}/{pname.group(1)}_p{pname.group(2)}.json"))
    test = set(map(tuple, p["test_pos"])); val = set(map(tuple, p["val_pos"]))
    grid = [tuple(x) for x in p["eval_grid"]]; vgrid = [tuple(x) for x in p["val_grid"]]
    y = np.array([q in test for q in grid], float); vy = np.array([q in val for q in vgrid], float)
    hist = []
    for f in sorted(glob.glob(f"{KGE}/{tag}_ep*_scores.npz"), key=lambda s: int(re.search(r"_ep(\d+)_", s).group(1))):
        ep = int(re.search(r"_ep(\d+)_", f).group(1))
        z = np.load(f)
        hist.append({"epoch": ep, "val_AP": L.weighted_ap(vy, z["val"].astype(float)),
                     "test": L.grid_metrics(grid, y, z["eval"].astype(float))})
    best = max(hist, key=lambda h: h["val_AP"])
    rows.append({"tag": tag, "neg": meta["neg"], "task": pname.group(1), "pseed": int(pname.group(2)),
                 "selected_epoch": best["epoch"], "val_AP": best["val_AP"], **{f"test_{k}": v for k, v in best["test"].items()},
                 "history": [{"epoch": h["epoch"], "val_AP": h["val_AP"], "test_AP": h["test"]["AP"]} for h in hist],
                 "train_seconds": meta["train_seconds"]})
json.dump(rows, open(f"{OUT}/e2e_ablation.json", "w"), indent=1)
for r in rows:
    print(r["tag"], "ep", r["selected_epoch"], "valAP %.3f" % r["val_AP"], "testAP %.3f macro %.3f MRR %.3f H10 %.3f AUROC %.3f" % (
        r["test_AP"], r["test_macroAP_disease"], r["test_MRR"], r["test_Hits@10"], r["test_AUROC"]))
