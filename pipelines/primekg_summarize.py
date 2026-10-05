"""Re-fit the reference consistently and aggregate the completed PrimeKG runs."""
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from kg_audit.metrics import ap_function, fit_degree_reference, assert_pair_contract, degree_bin

BASE = Path(os.environ["KG_AUDIT_WORKDIR"]) / "outputs/primekg"
METHODS = ["degree", "disease_only", "complex", "distmult", "rotate"]


def seed_intervals(frame, n=2000):
    """Disease and compound cluster CIs; do not average endpoints between seeds."""
    y = frame.label.to_numpy()
    fs = {m: ap_function(y, frame[m]) for m in METHODS}
    result = {}
    for cluster in ["disease", "compound"]:
        ids, names = pd.factorize(frame[cluster], sort=True)
        rng = np.random.default_rng(20261003)
        dist = {m: [] for m in METHODS if m not in {"degree", "disease_only"}}
        for _ in range(n):
            w = rng.multinomial(len(names), np.full(len(names), 1/len(names)))[ids]
            base = fs['degree'](w)
            for m in dist:
                value = fs[m](w)-base
                if np.isfinite(value):
                    dist[m].append(value)
        result[cluster] = {m:{"difference":fs[m]()-fs['degree'](), "CI95":np.quantile(d,[.025,.975]).tolist(),
                              "resamples":len(d)} for m,d in dist.items()}
    return result


def main():
    all_results = {}
    for split in ["random", "compound_disjoint"]:
        for scheme in ["uniform", "matched"]:
            cells = []
            for seed in [1,2,3]:
                folder = BASE / f"seed{seed}_{split}"
                if not (folder / 'complete.json').exists():
                    raise RuntimeError(f"Missing completed run: {folder.name}")
                frame = pd.read_csv(folder / f'{scheme}.tsv.gz',sep='\t')
                tp = list(pd.read_csv(folder/'train_positive.tsv',sep='\t',header=None).itertuples(index=False,name=None))
                tn = list(pd.read_csv(folder/'baseline_train_negative.tsv',sep='\t',header=None).itertuples(index=False,name=None))
                pairs = list(zip(frame.compound,frame.disease))
                assert_pair_contract(pairs,frame.label.tolist(),tp+tn)
                predict,cdeg,ddeg = fit_degree_reference(tp,tn,seed)
                frame['degree'] = predict(pairs)
                if split == 'compound_disjoint':
                    assert all(cdeg[c]==0 for c,d in pairs)
                # Exact joint-bin counts, not just independently matched marginal counts.
                if scheme == 'matched':
                    posbins = {}
                    negbins = {}
                    for (c,d),label in zip(pairs,frame.label):
                        key = (degree_bin(cdeg[c]),degree_bin(ddeg[d]))
                        dest = posbins if label else negbins
                        dest[key] = dest.get(key,0)+1
                    assert {k:20*v for k,v in posbins.items()}==negbins
                frame.to_csv(folder/f'{scheme}.tsv.gz',sep='\t',index=False)
                scores = {m:{'AP':float(average_precision_score(frame.label,frame[m])),
                              'AUROC':float(roc_auc_score(frame.label,frame[m]))} for m in METHODS}
                intervals = seed_intervals(frame)
                cell={'seed':seed,'scores':scores,'intervals':intervals,
                      'n_positive':int(frame.label.sum()),'n_negative_draws':int((frame.label==0).sum()),
                      'n_unique_negative':len(set(zip(frame.loc[frame.label==0,'compound'],frame.loc[frame.label==0,'disease']))),
                      'n_train_positive':len(tp),'n_baseline_train_negative':len(tn)}
                completed=json.loads((folder/'complete.json').read_text())
                completed[scheme]={**scores,'n_positive':cell['n_positive'],
                               'n_negative_draws':cell['n_negative_draws'],
                               'n_unique_negative':cell['n_unique_negative']}
                completed['reference_solver']='LIBLINEAR; recalculated on stored fitting pairs'
                (folder/'complete.json').write_text(json.dumps(completed,indent=2))
                cells.append(cell)
                print(split,scheme,seed,{m:round(s['AP'],4) for m,s in scores.items()},flush=True)
            summary={m:{'AP_mean':float(np.mean([c['scores'][m]['AP'] for c in cells])),
                        'AP_sd':float(np.std([c['scores'][m]['AP'] for c in cells],ddof=1)),
                        'AP_seeds':[c['scores'][m]['AP'] for c in cells],
                        'AUROC_mean':float(np.mean([c['scores'][m]['AUROC'] for c in cells])),
                        'AUROC_sd':float(np.std([c['scores'][m]['AUROC'] for c in cells],ddof=1))} for m in METHODS}
            all_results[f'{split}/{scheme}']={'scores':summary,'seed_results':cells}
            (BASE/'results.json').write_text(json.dumps(all_results,indent=2))


if __name__=='__main__':
    main()
