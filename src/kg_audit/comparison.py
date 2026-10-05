"""Paired average-precision comparisons on identical candidate rows."""
import numpy as np
import pandas as pd
from .metrics import ap_function


def paired_interval(labels,model_scores,reference_scores,compounds,diseases,
                    clusters='disease',n_bootstrap=5000,seed=20261003):
    labels=np.asarray(labels)
    model=np.atleast_2d(np.asarray(model_scores,dtype=float))
    reference=np.atleast_2d(np.asarray(reference_scores,dtype=float))
    if model.shape!=reference.shape or model.shape[1]!=len(labels):
        raise ValueError('Model, reference, and labels must share the same rows and seed count')
    if not np.isin(labels,[0,1]).all() or len(np.unique(labels))!=2:
        raise ValueError('Binary labels with both classes are required')
    if len(compounds)!=len(labels) or len(diseases)!=len(labels):raise ValueError('Cluster labels are not aligned')
    if clusters not in {'compound','disease','two_way'}:raise ValueError('Unknown cluster scheme')
    if n_bootstrap<100:raise ValueError('At least 100 bootstrap draws are required')
    mf=[ap_function(labels,s) for s in model];rf=[ap_function(labels,s) for s in reference]
    point=float(np.mean([m()-r() for m,r in zip(mf,rf)]))
    cids,cs=pd.factorize(compounds,sort=True);dids,ds=pd.factorize(diseases,sort=True)
    if min(cids.min(),dids.min())<0:raise ValueError('Missing cluster identifiers')
    rng=np.random.default_rng(seed);draws=[]
    for _ in range(n_bootstrap):
        cw=rng.multinomial(len(cs),np.full(len(cs),1/len(cs))) if clusters in {'compound','two_way'} else np.ones(len(cs))
        dw=rng.multinomial(len(ds),np.full(len(ds),1/len(ds))) if clusters in {'disease','two_way'} else np.ones(len(ds))
        weight=cw[cids]*dw[dids]
        value=float(np.mean([m(weight)-r(weight) for m,r in zip(mf,rf)]))
        if np.isfinite(value):draws.append(value)
    if not draws:raise ValueError('No non-degenerate bootstrap draws')
    return {'difference':point,'CI95':np.quantile(draws,[.025,.975]).tolist(),
            'valid_resamples':len(draws),'requested_resamples':n_bootstrap,'clusters':clusters,
            'interpretation':'Exploratory percentile interval conditional on supplied fitted scores and candidates; unadjusted for multiplicity.'}
