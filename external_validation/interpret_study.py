"""Compute paired scientific contrasts and cohort-level association from real study outputs."""
from __future__ import annotations

import argparse
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from external_validation.run import dataset_dir, _write_json


def fit_rank_scale(x, weights):
    center = x - np.average(x,axis=0,weights=weights)
    s = np.sqrt(np.average(center**2,axis=0,weights=weights))
    good = s>1e-12
    return center[:,good]/s[good]*np.sqrt(weights[:,None])


def association(subjects, protein, metab):
    """Equal-patient rank association, with a separate permutation reference for each estimand."""
    ids = subjects.patient_id.to_numpy()
    members = [np.flatnonzero(ids==p) for p in np.unique(ids)]
    members = [i for i in members if len(i)>=2]
    ix = np.concatenate(members)
    offsets = np.cumsum([0]+[len(i) for i in members])
    blocks = [slice(a,b) for a,b in zip(offsets[:-1],offsets[1:])]
    matrices = []
    for array in (protein,metab):
        x = array[ix]
        # The association estimand uses the complete-feature subset of the
        # model's seed-42 panel, rather than counting an unobserved cell as zero.
        keep = np.isfinite(x).all(0) & (np.nanstd(x,axis=0)>0)
        matrices.append(rankdata(x[:,keep],axis=0))
    rng = np.random.default_rng(1341)
    n = len(blocks)
    def rms(x,y,weights,permutations):
        a,b = fit_rank_scale(x,weights),fit_rank_scale(y,weights)
        if not a.shape[1] or not b.shape[1]:
            return None
        ga,gb = a@a.T,b@b.T
        denom = a.shape[1]*b.shape[1]
        score = float(np.sqrt(max(0,(ga*gb).sum()/denom)))
        ref = [np.sqrt(max(0,(ga*gb[np.ix_(j,j)]).sum()/denom)) for j in permutations]
        return {'rms':score,'reference_median':float(np.median(ref)),
                'reference_95':np.quantile(ref,[.025,.975]).tolist(),
                'protein_features':a.shape[1],'metab_features':b.shape[1]}
    between = [np.stack([x[s].mean(0) for s in blocks]) for x in matrices]
    within = [np.concatenate([x[s]-x[s].mean(0) for s in blocks]) for x in matrices]
    w = np.concatenate([np.full(len(i),1/(n*len(i))) for i in members])
    result = {'patients':n,'visits':len(ix),'between':rms(*between,np.full(n,1/n),[rng.permutation(n) for _ in range(199)]),
              'within':rms(*within,w,[np.concatenate([np.roll(np.arange(s.start,s.stop),rng.integers(s.stop-s.start))
                        for s in blocks]) for _ in range(199)])}
    return result


def effects(file):
    data = np.load(file,allow_pickle=False)
    direction = file.parent.name
    cohort = file.parents[2].name
    seed = int(file.parents[1].name.removeprefix('seed-'))
    truth,valid = data['truth'],data['valid']
    patients = data['patient_id']
    arms = [k for k in data.files if k not in {'truth','valid','patient_id','sample_id','gap','history_n'}]
    residual = {k:np.abs(data[k]-truth)*valid for k in arms}
    errors = {k:a.sum(1)/valid.sum(1) for k,a in residual.items()}
    paired = [('history','current','all'),('history','no_alignment','all'),('history','persistence','all'),
              ('forecast','persistence','all'),('forecast','mean','all'),('current','mean','all'),
              ('fusion_off','history','all'),('past_source_only','history','all'),
              ('reverse_own','history','multiple_history'),('history','no_chronology','multiple_history')]
    masks = {'all':np.ones(len(patients),bool),'multiple_history':data['history_n']>=2}
    if cohort in ('ipop','hmp2'):
        lo,hi = (30,120) if cohort=='ipop' else (8,12)
        for label,keep in (('gap_short',data['gap']<=lo),('gap_middle',(data['gap']>lo)&(data['gap']<=hi)),
                           ('gap_long',data['gap']>hi)):
            masks[label] = keep
            paired.extend((('history','current',label),('forecast','persistence',label)))
    results = []
    rng = np.random.default_rng(772)
    for a,b,subset in paired:
        keep = masks[subset]
        if a not in errors or b not in errors or not keep.any():
            continue
        delta = pd.Series(errors[a][keep]-errors[b][keep]).groupby(patients[keep]).mean().to_numpy()
        ci = np.quantile(rng.choice(delta,(5000,len(delta)),replace=True).mean(1),[.025,.975])
        results.append({'cohort':cohort,'seed':seed,'direction':direction,'arm':a,'reference':b,'subset':subset,
                        'patients':len(delta),'visits':int(keep.sum()),'delta_zmae':float(delta.mean()),
                        'ci_low':float(ci[0]) if len(delta)>=2 else None,'ci_high':float(ci[1]) if len(delta)>=2 else None})
    feature_names = json.loads((file.parent/'features.json').read_text())[direction.split('_to_')[1]]
    sums = residual['history'].sum(0)
    worst = np.argsort(-sums)[:3]
    diagnostic = {'cohort':cohort,'seed':seed,'direction':direction,
                  'top_3_feature_fraction_of_cell_error':float(sums[worst].sum()/sums.sum()),
                  'max_absolute_observed_target_z':float(np.abs(truth[valid]).max()),
                  'history_zmae_median_per_feature':float(np.median(sums[valid.sum(0)>0]/valid.sum(0)[valid.sum(0)>0])),
                  'top_3_features':[feature_names[i] for i in worst]}
    return results,diagnostic


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,required=True)
    p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    records, diagnostics, datasets = [],[],{}
    for f in sorted({f for root in a.runs for f in root.rglob('heldout.npz')}):
        rows,d = effects(f)
        records.extend(rows)
        diagnostics.append(d)
        if d['seed']==42:
            datasets[d['cohort']] = f.parent/'features.json'
    if not records:
        raise ValueError('no completed held-out experiments found')
    pd.DataFrame(records).to_csv(a.output/'matched_effects.csv',index=False)
    _write_json(a.output/'feature_error_diagnostics.json',diagnostics)
    results={}
    for name, features in datasets.items():
        folder=dataset_dir(a.data_root,name)
        subjects=pd.read_csv(folder/'subjects.csv',dtype={'sample_id':str,'patient_id':str}).sort_values(['patient_id','visit_order'])
        selected=json.loads(features.read_text())
        data=[]
        for mod, filename in [('protein','proteomics.csv'),('metab','metabolomics.csv')]:
            d=pd.read_csv(folder/filename,usecols=['sample_id',*selected[mod]],dtype={'sample_id':str})
            data.append(d.set_index('sample_id').reindex(subjects.sample_id).to_numpy(dtype=float))
        paired=np.isfinite(data[0]).any(1)&np.isfinite(data[1]).any(1)
        results[name]=association(subjects.loc[paired],data[0][paired],data[1][paired])
    _write_json(a.output/'rank_association.json',results)
    print(a.output/'matched_effects.csv')


if __name__=='__main__':
    main()
