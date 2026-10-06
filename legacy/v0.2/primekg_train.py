"""Fresh PrimeKG benchmark with task-compatible candidates and disjoint fitting.

Inputs are supplied locally under KG_AUDIT_WORKDIR. Checkpoint selection uses
task-consistent validation labels; test scores are computed after selection.
"""

from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score, roc_auc_score

from kg_audit.metrics import assert_pair_contract, sample_candidates, degree_bin, fit_degree_reference
from kg_audit.splits import validation_partition

ROOT = Path(os.environ["KG_AUDIT_WORKDIR"])
OUT = ROOT / "outputs/primekg"
OUT.mkdir(exist_ok=True, parents=True)
DIM, EPOCHS, BATCH, NNEG = 100, int(os.environ.get("MAX_EPOCHS", "40")), 16384, 4
torch.set_num_threads(int(os.environ.get("THREADS", "4")))
FEATURE_RELS = {"drug_protein", "disease_protein", "protein_protein"}
MODELS = ["complex", "distmult", "rotate"]
SEEDS = [int(x) for x in os.environ.get("SEEDS", "1 2 3").split()]


class KGE(torch.nn.Module):
    """Same scoring functions and optimiser budget as the historical script."""

    def __init__(self, nentity, nrelation, kind):
        super().__init__()
        self.kind = kind
        self.er = torch.nn.Embedding(nentity, DIM)
        self.rr = torch.nn.Embedding(nrelation, DIM)
        torch.nn.init.xavier_uniform_(self.er.weight)
        torch.nn.init.xavier_uniform_(self.rr.weight)
        if kind != "distmult":
            self.ei = torch.nn.Embedding(nentity, DIM)
            torch.nn.init.xavier_uniform_(self.ei.weight)
        if kind == "complex":
            self.ri = torch.nn.Embedding(nrelation, DIM)
            torch.nn.init.xavier_uniform_(self.ri.weight)
        if kind == "rotate":
            torch.nn.init.uniform_(self.rr.weight, -np.pi, np.pi)

    def score(self, h, r, t):
        hr, rr, tr = self.er(h), self.rr(r), self.er(t)
        if self.kind == "distmult":
            return (hr * rr * tr).sum(-1)
        hi, ti = self.ei(h), self.ei(t)
        if self.kind == "rotate":
            ar = hr * torch.cos(rr) - hi * torch.sin(rr)
            ai = hr * torch.sin(rr) + hi * torch.cos(rr)
            return 6.0 - torch.sqrt((ar - tr) ** 2 + (ai - ti) ** 2 + 1e-9).sum(-1)
        ri = self.ri(r)
        return (hr * rr * tr + hr * ri * ti + hi * rr * ti - hi * ri * tr).sum(-1)


def train(kind, triples, nentity, nrelation, seed, forbidden, indication_id, validation, e2i, destination):
    torch.manual_seed(seed)
    model = KGE(nentity, 2 * nrelation, kind)
    opt = torch.optim.Adam(model.parameters(), lr=0.003, weight_decay=1e-6)
    lossfn = torch.nn.BCEWithLogitsLoss()
    h, r, t = (torch.tensor(triples[:, j]) for j in range(3))
    hh, rr, tt = torch.cat([h, t]), torch.cat([r, r + nrelation]), torch.cat([t, h])
    # Known indication and evaluation pairs are never sampled as training negatives.
    forbidden_forward = {}
    forbidden_reverse = {}
    for c, d in forbidden:
        forbidden_forward.setdefault(c, set()).add(d)
        forbidden_reverse.setdefault(d, set()).add(c)
    start = time.time()
    best = -np.inf
    best_state = None
    wait = 0
    best_epoch = 0
    history = []
    vpairs, vy = validation
    for epoch in range(EPOCHS):
        perm = torch.randperm(len(hh))
        for begin in range(0, len(hh), BATCH):
            ix = perm[begin : begin + BATCH]
            ph, pr, pt = hh[ix], rr[ix], tt[ix]
            bs = len(ix)
            nt = torch.randint(nentity, (bs, NNEG))
            for relation, lookup in [
                (indication_id, forbidden_forward),
                (indication_id + nrelation, forbidden_reverse),
            ]:
                rows = torch.where(pr == relation)[0].tolist()
                for j in rows:
                    blocked = lookup.get(int(ph[j]), set())
                    for k in range(NNEG):
                        while int(nt[j, k]) in blocked:
                            nt[j, k] = torch.randint(nentity, ())
            p = model.score(ph, pr, pt)
            n = model.score(
                ph[:, None].expand(-1, NNEG).reshape(-1), pr[:, None].expand(-1, NNEG).reshape(-1), nt.reshape(-1)
            ).reshape(bs, NNEG)
            scores = torch.cat([p[:, None], n], 1)
            labels = torch.zeros_like(scores)
            labels[:, 0] = 1
            loss = lossfn(scores, labels)
            opt.zero_grad()
            loss.backward()
            opt.step()
        if epoch % 4 == 3:
            print(f"{kind} seed={seed} epoch={epoch + 1}/{EPOCHS} elapsed={time.time() - start:.0f}s", flush=True)
        if epoch + 1 >= 8 and (epoch + 1) % 2 == 0:
            score = predict(model, vpairs, e2i, indication_id)
            ap = float(average_precision_score(vy, score))
            history.append({"epoch": epoch + 1, "validation_AP": ap})
            if ap > best + 1e-4:
                best = ap
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
                wait = 0
                best_epoch = epoch + 1
            else:
                wait += 1
            if wait >= 3:
                break
    if best_state is None:
        raise ValueError("Training did not reach a validation checkpoint")
    model.load_state_dict(best_state)
    torch.save(best_state, destination / f"{kind}_state.pt")
    (destination / f"{kind}_training.json").write_text(
        json.dumps(
            {
                "seed": seed,
                "best_epoch": best_epoch,
                "stopped_epoch": epoch + 1,
                "best_validation_AP": best,
                "history": history,
                "torch_version": torch.__version__,
            },
            indent=2,
        )
    )
    return model.eval()


@torch.no_grad()
def predict(model, pairs, e2i, rid):
    out = []
    for start in range(0, len(pairs), 8192):
        block = pairs[start : start + 8192]
        h = torch.tensor([e2i[c] for c, d in block])
        t = torch.tensor([e2i[d] for c, d in block])
        out.append(model.score(h, torch.full_like(h, rid), t).numpy())
    return np.concatenate(out)


def main():
    features, positives, proxy = [], set(), set()
    for chunk in pd.read_csv(
        ROOT / "data/primekg/kg.csv",
        chunksize=300000,
        dtype=str,
        usecols=["relation", "x_type", "x_id", "y_type", "y_id"],
    ):
        for row in chunk.itertuples(index=False):
            # Access column names because the source column order is not a contract.
            rel, xt, xi, yt, yi = row.relation, row.x_type, row.x_id, row.y_type, row.y_id
            if rel in FEATURE_RELS:
                features.append((f"{xt}:{xi}", rel, f"{yt}:{yi}"))
            elif rel in {"indication", "contraindication", "off-label use"}:
                pair = (f"drug:{xi}", f"disease:{yi}") if xt == "drug" else (f"drug:{yi}", f"disease:{xi}")
                (positives if rel == "indication" else proxy).add(pair)
    positives = sorted(positives)
    allpos = set(positives)
    drugs = sorted({c for c, d in positives})
    diseases = sorted({d for c, d in positives})
    entities = sorted({x for h, r, t in features for x in [h, t]} | set(drugs) | set(diseases))
    e2i = {x: i for i, x in enumerate(entities)}
    rels = sorted(FEATURE_RELS) + ["indication"]
    r2i = {x: i for i, x in enumerate(rels)}
    ft = np.array([(e2i[h], r2i[r], e2i[t]) for h, r, t in features], dtype=np.int64)
    del features
    print(f"PrimeKG: {len(entities)} nodes, {len(ft)} feature edges, {len(positives)} indications", flush=True)
    digest = hashlib.sha256()
    with (ROOT / "data/primekg/kg.csv").open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    meta = {
        "max_epochs": EPOCHS,
        "dimension": DIM,
        "feature_edges": len(ft),
        "entities": len(entities),
        "positives": len(positives),
        "excluded_proxy_pairs": len(proxy),
        "seeds": SEEDS,
        "sampling": "20 draws per positive, with replacement; both endpoints matched on training-degree bins",
        "training": "Validation-AP checkpoint selection; minimum 8 epochs, evaluation every 2 epochs, patience 3. Ten percent edge validation for random task; ten percent compound-label validation for compound-disjoint task.",
        "input_sha256": digest.hexdigest(),
        "training_code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    signature = hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest()
    (OUT / "protocol.json").write_text(json.dumps(meta, indent=2))
    for seed in SEEDS:
        rng = np.random.default_rng(seed)
        ix = set(rng.permutation(len(positives))[: int(0.3 * len(positives))])
        held = set(rng.choice(drugs, int(0.3 * len(drugs)), replace=False))
        for split in ["random", "compound_disjoint"]:
            destination = OUT / f"seed{seed}_{split}"
            destination.mkdir(exist_ok=True)
            if (destination / "recipe.json").exists() and json.loads((destination / "recipe.json").read_text()).get(
                "signature"
            ) != signature:
                raise ValueError(
                    "Output cache belongs to a different dataset or training recipe; select a new work directory"
                )
            (destination / "recipe.json").write_text(json.dumps({"signature": signature, "protocol": meta}, indent=2))
            if (destination / "complete.json").exists():
                print(f"Reuse completed {destination.name}", flush=True)
                continue
            test = [p for i, p in enumerate(positives) if (i in ix if split == "random" else p[0] in held)]
            testset = set(test)
            non_test = [p for p in positives if p not in testset]
            trainpos, valpos = validation_partition(non_test, split, 700000 + seed + len(split))
            cdeg, ddeg = Counter(c for c, d in trainpos), Counter(d for c, d in trainpos)
            candidate_drugs = drugs if split == "random" else sorted(held)
            evalneg = {
                scheme: sample_candidates(
                    test,
                    candidate_drugs,
                    diseases,
                    allpos | proxy,
                    cdeg,
                    ddeg,
                    scheme,
                    20,
                    10000 * seed + 10 * len(split) + j,
                )
                for j, scheme in enumerate(["uniform", "matched"])
            }
            evalall = set(test) | set(evalneg["uniform"]) | set(evalneg["matched"])
            valdrugs = drugs if split == "random" else sorted({c for c, d in valpos})
            valneg = sample_candidates(
                valpos,
                valdrugs,
                diseases,
                allpos | proxy | evalall,
                cdeg,
                ddeg,
                "uniform",
                20,
                800000 + seed + len(split),
            )
            valall = set(valpos) | set(valneg)
            traindrugs = sorted({c for c, d in trainpos})
            trneg = sample_candidates(
                trainpos,
                traindrugs,
                diseases,
                allpos | proxy | evalall | valall,
                cdeg,
                ddeg,
                "uniform",
                1,
                900000 + seed + len(split),
            )
            trneg = sorted(set(trneg))
            # Fitting negatives are explicit and disjoint from both test candidate sets.
            assert_pair_contract(list(evalall), [int(p in testset) for p in evalall], trainpos + trneg)
            degree_predict, cdeg, ddeg = fit_degree_reference(trainpos, trneg, seed)
            manifests = {}
            for scheme, ng in evalneg.items():
                pairs = test + ng
                frame = pd.DataFrame(pairs, columns=["compound", "disease"])
                frame["label"] = np.r_[np.ones(len(test), dtype=int), np.zeros(len(ng), dtype=int)]
                frame["degree"] = degree_predict(pairs)
                frame["compound_degree"] = [cdeg[c] for c, d in pairs]
                frame["disease_degree"] = [ddeg[d] for c, d in pairs]
                frame["disease_only"] = [np.log1p(ddeg[d]) for c, d in pairs]
                manifests[scheme] = frame
                if split == "compound_disjoint":
                    assert frame.compound_degree.eq(0).all()
                frame.to_csv(destination / f"{scheme}.tsv.gz", sep="\t", index=False)
            np.savetxt(destination / "train_positive.tsv", np.asarray(trainpos), fmt="%s", delimiter="\t")
            np.savetxt(destination / "baseline_train_negative.tsv", np.asarray(trneg), fmt="%s", delimiter="\t")
            np.savetxt(destination / "validation_positive.tsv", np.asarray(valpos), fmt="%s", delimiter="\t")
            np.savetxt(destination / "validation_negative.tsv", np.asarray(valneg), fmt="%s", delimiter="\t")
            triples = np.concatenate([ft, np.array([(e2i[c], r2i["indication"], e2i[d]) for c, d in trainpos])])
            forbidden = {(e2i[c], e2i[d]) for c, d in allpos | proxy | evalall | valall if c in e2i and d in e2i}
            summary = {}
            for kind in MODELS:
                scorepath = destination / f"{kind}_scores.npz"
                if scorepath.exists():
                    arrays = np.load(scorepath)
                else:
                    model = train(
                        kind,
                        triples,
                        len(entities),
                        len(rels),
                        seed,
                        forbidden,
                        r2i["indication"],
                        (valpos + valneg, np.r_[np.ones(len(valpos)), np.zeros(len(valneg))]),
                        e2i,
                        destination,
                    )
                    arrays = {
                        scheme: predict(model, list(zip(f.compound, f.disease)), e2i, r2i["indication"])
                        for scheme, f in manifests.items()
                    }
                    np.savez_compressed(scorepath, **arrays)
                    del model
                for scheme, frame in manifests.items():
                    frame[kind] = arrays[scheme]
                    frame.to_csv(destination / f"{scheme}.tsv.gz", sep="\t", index=False)
            for scheme, frame in manifests.items():
                summary[scheme] = {
                    method: {
                        "AP": float(average_precision_score(frame.label, frame[method])),
                        "AUROC": float(roc_auc_score(frame.label, frame[method])),
                    }
                    for method in ["degree", "disease_only"] + MODELS
                }
                summary[scheme].update(
                    n_positive=len(test),
                    n_negative_draws=len(evalneg[scheme]),
                    n_unique_negative=len(set(evalneg[scheme])),
                )
            (destination / "complete.json").write_text(json.dumps(summary, indent=2))
            print(destination.name, summary, flush=True)


if __name__ == "__main__":
    main()
