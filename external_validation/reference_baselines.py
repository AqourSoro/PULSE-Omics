"""Evaluate order-invariant references on the saved PULSE patient splits and panels."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from external_validation import run
from external_validation.real_study import metrics

PENALTIES = np.array([.01, .1, 1., 10., 100.])


def ridge_paths(train_x, test_x, train_y, valid):
    """Fit with one kernel eigendecomposition per target observation pattern."""
    gram = train_x @ train_x.T
    cross = test_x @ train_x.T
    out = np.zeros((len(PENALTIES), len(test_x), train_y.shape[1]))
    groups = {}
    for col in range(valid.shape[1]):
        groups.setdefault(valid[:, col].tobytes(), []).append(col)
    for cols in groups.values():
        rows = np.flatnonzero(valid[:, cols[0]])
        if not len(rows):
            continue
        y = train_y[np.ix_(rows, cols)]
        center = y.mean(0)
        k = gram[np.ix_(rows, rows)]
        col_mean = k.mean(0)
        eval_k = cross[:, rows]
        eval_k = eval_k - eval_k.mean(1, keepdims=True) - col_mean + col_mean.mean()
        k = k - col_mean[:, None] - col_mean + col_mean.mean()
        eigenvalues, vectors = np.linalg.eigh(k)
        projected = vectors.T @ (y - center)
        alpha = projected[None, :, :] / (np.maximum(eigenvalues, 0)[None, :, None] + PENALTIES[:, None, None])
        out[:, :, cols] = (eval_k @ vectors) @ alpha + center
    return out


def reference_risk(file, root):
    folder = file.parent
    cohort = file.parents[2].name
    seed = int(file.parents[1].name.removeprefix('seed-'))
    source, target = folder.name.split('_to_')
    inputs = run.dataset_dir(root, cohort)
    subjects = pd.read_csv(inputs/'subjects.csv',dtype={'sample_id':str,'patient_id':str})
    subjects = subjects.sort_values(['patient_id','visit_order'])
    panel = json.loads((folder/'features.json').read_text())
    raw = {}
    for mod, name in zip(run.MODALITIES, ('proteomics.csv', 'metabolomics.csv')):
        frame = pd.read_csv(inputs/name,usecols=['sample_id',*panel[mod]],dtype={'sample_id':str})
        raw[mod] = frame.set_index('sample_id').reindex(subjects.sample_id)[panel[mod]]
    split = json.loads((folder/'split.json').read_text())
    history_limit = json.loads((file.parents[3]/'config.json').read_text())['max_history']
    cases = run.eligible_cases(subjects, raw, source, target, True)
    samples = subjects.loc[subjects.patient_id.isin(split['train']), 'sample_id']
    norm = run.Normalizer(raw, samples.tolist())
    if norm.names != panel:
        raise ValueError(f'{folder}: selected feature panel differs')
    values, known = {}, {}
    for mod, frame in norm.values.items():
        matrix = frame.to_numpy()
        known[mod] = np.isfinite(matrix)
        values[mod] = np.where(known[mod], (matrix-norm.mean[mod])/norm.std[mod], 0)
    lookup = {str(s):i for i,s in enumerate(norm.values[source].index)}

    def indices(case):
        return lookup[case.sample], [lookup[s] for s in case.history]

    def usable(case):
        i, history = indices(case)
        return all(m[i].any() for m in known.values()) and any(
            all(m[h].any() for m in known.values()) for h in history)

    stored = np.load(file, allow_pickle=False)
    groups = {g:[c for c in cases if c.patient in ids and usable(c)] for g,ids in split.items()}
    if [c.sample for c in groups['test']] != stored['sample_id'].tolist():
        raise ValueError(f'{folder}: held-out case order differs')
    parts = {}
    for g, part in groups.items():
        rows, historical = zip(*(indices(c) for c in part))
        historical = [[i for i in h if any(mask[i].any() for mask in known.values())][-history_limit:]
                      for h in historical]
        averaged = {
            mod: np.stack([values[mod][h].sum(0) / known[mod][h].sum(0).clip(min=1) for h in historical])
            for mod in run.MODALITIES
        }
        current = values[source][list(rows)] / np.sqrt(values[source].shape[1])
        past = np.concatenate([v / np.sqrt(v.shape[1]) for v in averaged.values()], axis=1)
        parts[g] = {'current':current, 'history':np.concatenate([current,past],axis=1), 'forecast':past,
                    'history_mean':averaged[target], 'truth':values[target][list(rows)], 'known':known[target][list(rows)]}
    train, val, test = (parts[g] for g in ('train','val','test'))
    n_val = len(val['truth'])
    truth, valid = stored['truth'], stored['valid']
    if not np.array_equal(test['known'],valid) or not np.allclose(test['truth'],truth,rtol=1e-4,atol=1e-5):
        raise ValueError(f'{folder}: held-out target masks or standardized values differ')
    result = []
    rng = np.random.default_rng(772)
    for arm in ('current','history','forecast','history_mean'):
        if arm == 'history_mean':
            pred, penalty = test[arm], None
        else:
            x = np.concatenate([val[arm],test[arm]])
            all_predictions = ridge_paths(train[arm], x, train['truth'], train['known'])
            squared = (all_predictions[:, :n_val] - val['truth'])**2
            validation = ((squared * val['known']).sum(2) / val['known'].sum(1)).mean(1)
            best = int(validation.argmin())
            pred, penalty = all_predictions[best,n_val:], float(PENALTIES[best])
        base, base_patient = metrics(groups['test'], truth, pred, valid)
        pulse_arm = 'forecast' if arm == 'history_mean' else arm
        _, pulse_patient = metrics(groups['test'], truth, stored[pulse_arm], valid)
        delta = (pulse_patient-base_patient).to_numpy()
        ci = np.quantile(rng.choice(delta,(5000,len(delta)),replace=True).mean(1),[.025,.975])
        result.append({'cohort':cohort,'seed':seed,'direction':folder.name,'pulse_arm':pulse_arm,
                       'reference':'history_mean' if arm=='history_mean' else f'ridge_{arm}',
                       'ridge_penalty':penalty,**base,'pulse_minus_reference_zmae':float(delta.mean()),
                       'ci_low':float(ci[0]),'ci_high':float(ci[1])})
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,required=True)
    p.add_argument('--runs',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    rows = [row for file in sorted({f for root in a.runs for f in root.rglob('heldout.npz')})
            for row in reference_risk(file,a.data_root)]
    pd.DataFrame(rows).to_csv(a.output,index=False)
    print(a.output)
