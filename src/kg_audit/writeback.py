"""Four-policy semi-synthetic write-back stress test on a supplied split.

Selected recorded positive training edges receive hypothetical retain/defer
annotations. These are not observed trial failures or expert annotations.
"""

import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


def policy_edges(positive, negative, selected, retained, policy):
    selected = set(selected)
    retained = set(retained)
    if not retained <= selected:
        raise ValueError("Retained indices must be selected")
    if policy == "no_writeback":
        return list(positive), list(negative)
    if policy == "flat_negative":
        return [p for i, p in enumerate(positive) if i not in selected], list(negative) + [
            positive[i] for i in sorted(selected)
        ]
    if policy == "uncertainty_mask":
        return [p for i, p in enumerate(positive) if i not in selected], list(negative)
    if policy == "typed":
        return [p for i, p in enumerate(positive) if i not in selected or i in retained], list(negative)
    raise ValueError("Unknown write-back policy")


def run(splits_path, output, seeds=range(1, 11), epochs=200, threads=2):
    import torch
    from torch import nn

    torch.set_num_threads(threads)
    spec = json.loads(Path(splits_path).read_text())["random"]
    positive = list(map(tuple, spec["train_pos"]))
    testpos = list(map(tuple, spec["test_pos"]))
    rng = np.random.default_rng(20261005)
    eval_sets = {}
    for name in ["neg_random", "neg_degmatch"]:
        pool = list(map(tuple, spec[name]))
        idx = rng.choice(len(pool), 20 * len(testpos), replace=False)
        eval_sets[name] = testpos + [pool[i] for i in idx]
    excluded = set().union(*map(set, eval_sets.values()))
    negative = [tuple(p) for p in spec["train_neg"] if tuple(p) not in excluded]
    allpairs = positive + negative + list(excluded)
    c2i = {x: i for i, x in enumerate(sorted({c for c, d in allpairs}))}
    d2i = {x: i for i, x in enumerate(sorted({d for c, d in allpairs}))}

    def indexes(pairs):
        return torch.tensor([[c2i[c], d2i[d]] for c, d in pairs], dtype=torch.long)

    eval_index = {k: indexes(pairs) for k, pairs in eval_sets.items()}
    labels = {k: np.r_[np.ones(len(testpos)), np.zeros(len(pairs) - len(testpos))] for k, pairs in eval_sets.items()}

    class MF(nn.Module):
        def __init__(self):
            super().__init__()
            self.u = nn.Embedding(len(c2i), 16)
            self.v = nn.Embedding(len(d2i), 16)
            self.bc = nn.Embedding(len(c2i), 1)
            self.bd = nn.Embedding(len(d2i), 1)
            self.intercept = nn.Parameter(torch.zeros(()))
            nn.init.normal_(self.u.weight, 0, 0.1)
            nn.init.normal_(self.v.weight, 0, 0.1)
            nn.init.zeros_(self.bc.weight)
            nn.init.zeros_(self.bd.weight)

        def forward(self, ix):
            c, d = ix[:, 0], ix[:, 1]
            return (self.u(c) * self.v(d)).sum(1) + self.bc(c).squeeze(1) + self.bd(d).squeeze(1) + self.intercept

    rows = []
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    all_scores = {}
    start = time.time()
    fit_count = 0
    for seed in seeds:
        cache = {}
        prng = np.random.default_rng(seed + 50000)
        order = prng.permutation(len(positive))
        for fraction in [0.0, 0.1, 0.25, 0.5]:
            selected = order[: int(round(fraction * len(positive)))].tolist()
            for retain_fraction in [0.0, 0.5, 1.0]:
                retained = selected[: int(round(retain_fraction * len(selected)))]
                for policy in ["no_writeback", "flat_negative", "uncertainty_mask", "typed"]:
                    p, n = policy_edges(positive, negative, selected, retained, policy)
                    key = hashlib.sha256(json.dumps([p, n], separators=(",", ":")).encode()).hexdigest()
                    if key not in cache:
                        torch.manual_seed(seed)
                        model = MF()
                        opt = torch.optim.Adam(model.parameters(), lr=0.02, weight_decay=1e-4)
                        ix = indexes(p + n)
                        y = torch.tensor([1.0] * len(p) + [0.0] * len(n))
                        lossfn = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(len(n) / len(p)))
                        for _ in range(epochs):
                            loss = lossfn(model(ix), y)
                            opt.zero_grad()
                            loss.backward()
                            opt.step()
                        with torch.no_grad():
                            scores = {name: model(vals).numpy() for name, vals in eval_index.items()}
                        cache[key] = scores
                        fit_count += 1
                    for sampler, scores in cache[key].items():
                        row = {
                            "seed": seed,
                            "selected_fraction": fraction,
                            "n_selected": len(selected),
                            "retain_fraction": retain_fraction,
                            "policy": policy,
                            "sampler": sampler,
                            "AP": float(average_precision_score(labels[sampler], scores)),
                            "AUROC": float(roc_auc_score(labels[sampler], scores)),
                            "fit_id": key,
                            "n_train_positive": len(p),
                        }
                        rows.append(row)
                        all_scores[f"s{seed}_f{fraction}_r{retain_fraction}_{policy}_{sampler}"] = scores
        pd.DataFrame(rows).to_csv(out / "replicate_results.csv", index=False)
        print(f"seed {seed} complete; unique fits={fit_count}; elapsed={time.time() - start:.1f}s", flush=True)
    data = pd.DataFrame(rows)
    summary = data.groupby(["selected_fraction", "retain_fraction", "policy", "sampler"], as_index=False).agg(
        AP_mean=("AP", "mean"), AP_sd=("AP", "std"), AUROC_mean=("AUROC", "mean"), n_replicates=("seed", "nunique")
    )
    summary.to_csv(out / "summary.csv", index=False)
    np.savez_compressed(out / "scores.npz", **all_scores)
    for name, pairs in eval_sets.items():
        frame = pd.DataFrame(pairs, columns=["compound", "disease"])
        frame["label"] = labels[name]
        frame.to_csv(out / f"{name}_candidates.tsv", sep="\t", index=False)
    protocol = {
        "input_sha256": hashlib.sha256(Path(splits_path).read_bytes()).hexdigest(),
        "seeds": list(seeds),
        "epochs": epochs,
        "latent_dimension": 16,
        "learning_rate": 0.02,
        "weight_decay": 0.0001,
        "sampler_seed": 20261005,
        "unique_fits": fit_count,
        "n_training_positive": len(positive),
        "n_training_unlabelled": len(negative),
        "n_test_positive": len(testpos),
        "design": "Semi-synthetic annotations on recorded positive training edges, not observed failures. Typed retention is supplied correctly by construction; no real-world annotation reliability is estimated.",
    }
    (out / "protocol.json").write_text(json.dumps(protocol, indent=2))
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--splits", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--epochs", type=int, default=200)
    p.add_argument("--seeds", type=int, nargs="+", default=list(range(1, 11)))
    p.add_argument("--threads", type=int, default=2)
    a = p.parse_args()
    run(a.splits, a.output, a.seeds, a.epochs, a.threads)


if __name__ == "__main__":
    main()
