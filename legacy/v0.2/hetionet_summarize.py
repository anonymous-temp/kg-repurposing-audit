"""Auditable reanalysis of the frozen Hetionet prediction files.

Outputs include aligned candidate-level scores, reconstructed training and
validation labels, paired cluster intervals and original-number reconciliation.
"""

from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import os

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score, roc_auc_score

from kg_audit.metrics import ap_function, assert_pair_contract, fit_degree_reference

ROOT = Path(os.environ["KG_AUDIT_WORKDIR"])
OUT = ROOT / "outputs/hetionet/analysis"
OUT.mkdir(exist_ok=True, parents=True)
EXPB = ROOT / "outputs/hetionet/raw"
SEEDS = [1, 2, 3]
MODELS = ["complex", "distmult", "rotate", "rgcn"]
SAMPLING_SEED = 20261003
B = int(os.environ.get("BOOTSTRAPS", "5000"))


def permutation_probability(train, pairs, seed, samples=500):
    """Bipartite double-edge-swap baseline inspired by Zietz et al.

    Both endpoint degrees and simple-graph constraints are preserved. This is
    an explicit local implementation, not an execution of the xswap package.
    Burn-in 20E proposals; 5E proposals between 500 recorded graphs.
    """
    rng = np.random.default_rng(seed)
    edges, present = list(train), set(train)
    n = len(edges)
    query = set(pairs)
    counts = Counter()
    accepted = 0
    for step in range((20 + 5 * samples) * n):
        i, j = rng.integers(n, size=2)
        a, b = edges[i]
        c, d = edges[j]
        if a == c or b == d or (a, d) in present or (c, b) in present:
            pass
        else:
            present.remove((a, b))
            present.remove((c, d))
            present.add((a, d))
            present.add((c, b))
            edges[i], edges[j] = (a, d), (c, b)
            accepted += 1
        if step >= 20 * n and (step - 20 * n + 1) % (5 * n) == 0:
            counts.update(present & query)
    assert Counter(c for c, d in edges) == Counter(c for c, d in train)
    assert Counter(d for c, d in edges) == Counter(d for c, d in train)
    return np.array([counts[p] / samples for p in pairs]), accepted


def cluster_intervals(frame, scorecols, contrasts, scheme, seed):
    """Percentile CI of the mean seed-specific AP difference, on common draws."""
    y = frame.label.to_numpy()
    functions = {col: ap_function(y, frame[col].to_numpy()) for col in scorecols}
    cids, _ = pd.factorize(frame.compound, sort=True)
    dids, _ = pd.factorize(frame.disease, sort=True)
    nc, nd = cids.max() + 1, dids.max() + 1
    rng = np.random.default_rng(seed)
    distributions = {name: [] for name in contrasts}
    invalid = 0
    for _ in range(B):
        cw = rng.multinomial(nc, np.full(nc, 1 / nc)) if scheme in {"compound", "two_way"} else np.ones(nc)
        dw = rng.multinomial(nd, np.full(nd, 1 / nd)) if scheme in {"disease", "two_way"} else np.ones(nd)
        weights = cw[cids] * dw[dids]
        if weights[y == 1].sum() == 0 or weights[y == 0].sum() == 0:
            invalid += 1
            continue
        aps = {col: fn(weights) for col, fn in functions.items()}
        for name, (left, right) in contrasts.items():
            distributions[name].append(float(np.mean([aps[f"{left}_s{s}"] - aps[f"{right}_s{s}"] for s in SEEDS])))
    out = {}
    observed = {col: fn() for col, fn in functions.items()}
    for name, (left, right) in contrasts.items():
        dist = distributions[name]
        out[name] = {
            "difference": float(np.mean([observed[f"{left}_s{s}"] - observed[f"{right}_s{s}"] for s in SEEDS])),
            "CI95": np.quantile(dist, [0.025, 0.975]).tolist(),
            "valid_resamples": len(dist),
            "invalid_resamples": invalid,
        }
    return out


def main():
    splits = json.loads((EXPB / "splits.json").read_text())
    with gzip.open(ROOT / "data/hetionet/hetionet-edges.sif.gz", "rt") as fh:
        next(fh)
        all_ctd = [(a, b) for a, r, b in (line.strip().split("\t") for line in fh) if r == "CtD"]
    audit, results, legacy, hashes = {}, {}, {}, {}
    for split in os.environ.get("TASKS", "random,scaffold,coldstart").split(","):
        spec = splits[split]
        test = list(map(tuple, spec["test_pos"]))
        train_pool = [p for p in all_ctd if p not in set(test)]
        assert train_pool == list(map(tuple, spec["train_pos"]))
        candidates = {}
        for j, negkey in enumerate(["neg_random", "neg_degmatch"]):
            pool = list(map(tuple, spec[negkey]))
            # Frozen once across all models and training seeds; names encode actual sampling.
            sel = np.random.default_rng(SAMPLING_SEED + j + len(split)).choice(len(pool), 20 * len(test), replace=False)
            ng = [pool[i] for i in sel]
            f = pd.DataFrame(test + ng, columns=["compound", "disease"])
            f["label"] = np.r_[np.ones(len(test), dtype=int), np.zeros(len(ng), dtype=int)]
            candidates[negkey] = f
        eval_pairs = set().union(*(set(zip(f.compound, f.disease)) for f in candidates.values()))
        audit[split] = {
            "n_test_positive": len(test),
            "n_train_pool": len(train_pool),
            "n_test_drug": len({c for c, d in test}),
            "n_test_disease": len({d for c, d in test}),
            "seeds": {},
        }
        for seed in SEEDS:
            vrng = np.random.default_rng(seed * 100 + len(split))
            vsel = set(map(int, vrng.choice(len(train_pool), max(10, int(0.1 * len(train_pool))), replace=False)))
            train = [p for i, p in enumerate(train_pool) if i not in vsel]
            val = [p for i, p in enumerate(train_pool) if i in vsel]
            trneg = [tuple(p) for p in spec["train_neg"] if tuple(p) not in eval_pairs]
            assert_pair_contract(list(eval_pairs), [int(p in set(test)) for p in eval_pairs], train + trneg)
            degree_predict, cdeg, ddeg = fit_degree_reference(train, trneg, seed)
            union_pairs = sorted(eval_pairs)
            perm, accepted = permutation_probability(train, union_pairs, seed + 1000)
            permmap = dict(zip(union_pairs, perm))
            for nk, frame in candidates.items():
                pairs = list(zip(frame.compound, frame.disease))
                frame[f"degree_s{seed}"] = degree_predict(pairs)
                frame[f"disease_s{seed}"] = [np.log1p(ddeg[d]) for c, d in pairs]
                frame[f"permutation_s{seed}"] = [permmap[p] for p in pairs]
                frame[f"compound_degree_s{seed}"] = [cdeg[c] for c, d in pairs]
                frame[f"disease_degree_s{seed}"] = [ddeg[d] for c, d in pairs]
            for model in MODELS:
                path = EXPB / f"scores_{model}_{split}_seed{seed}.tsv"
                if split == "coldstart" and model != "rgcn":
                    path = EXPB / path.name
                hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
                source = pd.read_csv(path, sep="\t")
                for nk, frame in candidates.items():
                    selected = source[source.group.isin(["test_pos", nk])]
                    values = {(c, d): score for c, d, score in zip(selected.compound, selected.disease, selected.score)}
                    assert len(values) == len(selected)
                    frame[f"{model}_s{seed}"] = [values[p] for p in zip(frame.compound, frame.disease)]
                # Reconstruct the historical table and the incompatible significance sample.
                if split == "coldstart":
                    continue  # These are new runs, not historical table inputs.
                nk = "neg_random" if split == "random" else "neg_degmatch"
                positive = source[source.group == "test_pos"].score.to_numpy()
                negative = source[source.group == nk].score.to_numpy()
                old = np.random.default_rng(seed).choice(negative, 20 * len(positive), replace=False)
                yy = np.r_[np.ones(len(positive)), np.zeros(len(old))]
                legacy.setdefault(f"{split}/{model}", []).append(
                    {
                        "seed": seed,
                        "table_AP": float(average_precision_score(yy, np.r_[positive, old])),
                        "significance_AP": float(average_precision_score(yy, np.r_[positive, negative[: len(old)]])),
                    }
                )
            audit[split]["seeds"][seed] = {
                "actual_train_positive": len(train),
                "validation_positive": len(val),
                "baseline_train_negative": len(trneg),
                "removed_train_eval_negative_overlap": len(spec["train_neg"]) - len(trneg),
                "permutation_accepted_swaps": accepted,
            }
            (OUT / f"{split}_seed{seed}_training.json").write_text(
                json.dumps({"graph_positive": train, "validation_positive": val, "baseline_negative": trneg})
            )
            print(f"Hetionet {split} seed {seed} aligned", flush=True)
        for nk, frame in candidates.items():
            key = f"{split}/{nk}"
            frame.to_csv(OUT / f"{split}_{nk}_scores.tsv.gz", sep="\t", index=False)
            methods = MODELS + ["degree", "disease", "permutation"]
            scores = {}
            for method in methods:
                aps = [float(average_precision_score(frame.label, frame[f"{method}_s{s}"])) for s in SEEDS]
                aus = [float(roc_auc_score(frame.label, frame[f"{method}_s{s}"])) for s in SEEDS]
                scores[method] = {
                    "AP_mean": float(np.mean(aps)),
                    "AP_sd": float(np.std(aps, ddof=1)),
                    "AP_seeds": aps,
                    "AUROC_mean": float(np.mean(aus)),
                    "AUROC_sd": float(np.std(aus, ddof=1)),
                }
            contrasts = {f"{m}_minus_degree": (m, "degree") for m in MODELS}
            contrasts.update(
                {"rotate_minus_complex": ("rotate", "complex"), "rotate_minus_distmult": ("rotate", "distmult")}
            )
            used = [f"{m}_s{s}" for m in MODELS + ["degree"] for s in SEEDS]
            intervals = {
                scheme: cluster_intervals(frame, used, contrasts, scheme, SAMPLING_SEED)
                for scheme in ["disease", "compound", "two_way"]
            }
            results[key] = {
                "n_positive": len(test),
                "n_negative": 20 * len(test),
                "scores": scores,
                "intervals": intervals,
            }
            print(key, {m: round(v["AP_mean"], 4) for m, v in scores.items()}, flush=True)
        (OUT / "results.json").write_text(json.dumps(results, indent=2))
    hashes["outputs/hetionet/raw/splits.json"] = hashlib.sha256((EXPB / "splits.json").read_bytes()).hexdigest()
    source = ROOT / "data/hetionet/hetionet-edges.sif.gz"
    hashes["data/hetionet/hetionet-edges.sif.gz"] = hashlib.sha256(source.read_bytes()).hexdigest()
    (OUT / "source_hashes.json").write_text(json.dumps(hashes, indent=2))
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2))
    (OUT / "historical_reconciliation.json").write_text(json.dumps(legacy, indent=2))


if __name__ == "__main__":
    main()
