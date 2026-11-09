"""Independent real-omics experiments using the unmodified PULSE model.

Inputs: verified canonical sample / assay CSVs. Outputs: F: drive only.
No synthetic data are used for scientific results. Current targets are always masked.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from pulse.model import PulseModel

DEFAULT_ROOT = Path('/mnt/f/datasets/pulse-omics')
DATASETS = ('su2020', 'ipop', 'hmp2', 'covaar', 'csiro', 'p100')
MODALITIES = ('protein', 'metab')


@dataclass(frozen=True)
class Case:
    patient: str
    sample: str
    visit: float
    history: tuple[str, ...]


def load_canonical(folder: Path):
    """Match only explicit sample IDs, never spreadsheet row positions."""
    required = (folder / 'subjects.csv', folder / 'proteomics.csv', folder / 'metabolomics.csv')
    missing = [p.name for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError(f'{folder}: missing canonical files: {missing}')
    subjects = pd.read_csv(required[0], dtype={'sample_id': str, 'patient_id': str})
    for col in ('sample_id', 'patient_id', 'visit_order'):
        if col not in subjects or subjects[col].isna().any():
            raise ValueError(f'{folder}: invalid subjects.{col}')
    subjects['visit_order'] = pd.to_numeric(subjects['visit_order'], errors='raise')
    if not np.isfinite(subjects['visit_order']).all():
        raise ValueError('visit_order must be finite')
    if subjects.sample_id.duplicated().any() or subjects.duplicated(['patient_id', 'visit_order']).any():
        raise ValueError('duplicate sample ID or patient/visit assignment')
    subjects = subjects.sort_values(['patient_id', 'visit_order']).reset_index(drop=True)
    raw = {}
    valid_ids = set(subjects.sample_id)
    for mod, path in zip(MODALITIES, required[1:]):
        table = pd.read_csv(path, dtype={'sample_id': str})
        if 'sample_id' not in table or table.sample_id.isna().any() or table.sample_id.duplicated().any():
            raise ValueError(f'{path}: missing/duplicate sample ID')
        if not set(table.sample_id).issubset(valid_ids):
            raise ValueError(f'{path}: assay samples absent from subjects table')
        fields = [col for col in table if col != 'sample_id']
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f'{path}: missing/duplicate features')
        values = table[fields].apply(pd.to_numeric, errors='raise').astype('float64')
        if np.isinf(values.to_numpy()).any():
            raise ValueError(f'{path}: infinite measurements')
        values.index = table.sample_id
        raw[mod] = values.reindex(subjects.sample_id)
    return subjects, raw


def eligible_cases(subjects, raw, source, target, history):
    cases = []
    for patient, group in subjects.groupby('patient_id', sort=True):
        preceding = []
        for rec in group.itertuples(index=False):
            sample = rec.sample_id
            source_ok = raw[source].loc[sample].notna().any()
            target_ok = raw[target].loc[sample].notna().any()
            paired_preceding = [s for s in preceding if
                                 raw[source].loc[s].notna().any() and raw[target].loc[s].notna().any()]
            if source_ok and target_ok and (not history or paired_preceding):
                cases.append(Case(str(patient), str(sample), float(rec.visit_order),
                                  tuple(preceding)))
            preceding.append(str(sample))
    return cases


def split_cases(cases, seed):
    patients = sorted({c.patient for c in cases})
    if len(patients) < 8:
        raise ValueError(f'need >=8 eligible patients; found {len(patients)}')
    rng = np.random.default_rng(seed)
    rng.shuffle(patients)
    n_test = max(2, round(len(patients) * .15))
    n_val = max(2, round(len(patients) * .15))
    n_train = len(patients) - n_test - n_val
    if n_train < 4:
        raise ValueError('not enough training patients after splitting')
    ids = {'train': sorted(patients[:n_train]),
           'val': sorted(patients[n_train:n_train+n_val]),
           'test': sorted(patients[n_train+n_val:])}
    groups = {k: [c for c in cases if c.patient in set(v)] for k, v in ids.items()}
    if any(not v for v in groups.values()):
        raise ValueError('empty split after patient allocation')
    return ids, groups


class Normalizer:
    """Fit feature selection and z-scaling strictly on training patients."""

    def __init__(self, raw, train_samples, minimum=5):
        self.names, self.mean, self.std, self.values = {}, {}, {}, {}
        for mod in MODALITIES:
            frame = raw[mod]
            train = frame.loc[train_samples].to_numpy(dtype=float)
            n = np.isfinite(train).sum(axis=0)
            var = np.nanvar(train, axis=0)
            use = (n >= minimum) & np.isfinite(var) & (var > 1e-12)
            if use.sum() == 0:
                raise ValueError(f'{mod}: no nonconstant features with >= {minimum} training observations')
            self.names[mod] = np.array(frame.columns)[use].tolist()
            avg = np.nanmean(train[:, use], axis=0)
            std = np.nanstd(train[:, use], axis=0)
            self.mean[mod], self.std[mod] = avg, std
            self.values[mod] = frame.loc[:, self.names[mod]]

    def raw(self, mod, sample):
        return self.values[mod].loc[sample].to_numpy(dtype=np.float32)

    def scaled(self, mod, sample):
        original = self.raw(mod, sample)
        normal = (original - self.mean[mod]) / self.std[mod]
        return np.where(np.isfinite(normal), normal, 0).astype(np.float32)


def model_input(case, norm, source, target, use_history, max_history, device):
    history = case.history[-max_history:] if use_history else ()
    inputs, masks = [], []
    for sample in history:
        data, mask = {}, {}
        for mod in MODALITIES:
            val = norm.raw(mod, sample)
            observed = bool(np.isfinite(val).any())
            data[mod] = torch.from_numpy(norm.scaled(mod, sample)).to(device).unsqueeze(0) if observed else None
            mask[mod] = 0 if observed else 2
        if any(value == 0 for value in mask.values()):
            inputs.append(data)
            masks.append(mask)
    current = {
        source: torch.from_numpy(norm.scaled(source, case.sample)).to(device).unsqueeze(0),
        target: None,
    }
    inputs.append(current)
    masks.append({source: 0, target: 1})
    truth = norm.raw(target, case.sample)
    valid = np.isfinite(truth)
    if not valid.any():
        raise ValueError('case has no observed target values')
    scaled = norm.scaled(target, case.sample)
    return inputs, masks, (torch.from_numpy(scaled).to(device).unsqueeze(0),
                           torch.from_numpy(valid).to(device).unsqueeze(0))


def predict_and_loss(model, case, norm, source, target, history, args, train):
    visits, masks, (truth, valid) = model_input(
        case, norm, source, target, history, args.max_history, args.device)
    out = model(visits, masks, kl_annealing_weight=1.0)
    pred = out['reconstructions'][-1][target]
    # Missing feature cells are never used as supervised labels.
    masked_error = ((pred - truth).square() * valid).sum() / valid.sum()
    if train:
        losses = out['losses']
        loss = (masked_error + args.context_weight * losses['recon_loss']
                + args.kl_weight * losses['kl_loss']
                + (0 if args.no_alignment else args.align_weight * losses['alignment_loss']))
        return loss, pred
    return masked_error, pred


def fit_arm(groups, norm, source, target, history, no_alignment, args, folder):
    torch.manual_seed(args.seed)
    random.seed(args.seed)
    np.random.seed(args.seed)
    model = PulseModel(
        modality_dims={m: len(norm.names[m]) for m in MODALITIES},
        latent_dim=args.latent_dim, hidden_dim=args.hidden_dim,
        lambda_kl=0, lambda_align=0, lambda_adv=0,
        use_adversarial_loss=False, temporal_model=args.temporal_model,
    ).to(args.device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    args.no_alignment = no_alignment
    best, best_state, no_improve = float('inf'), None, 0
    rows = []
    for epoch in range(args.epochs):
        model.train()
        order = np.random.default_rng(args.seed + epoch).permutation(len(groups['train']))
        total = 0.0
        for idx in order:
            case = groups['train'][int(idx)]
            optimizer.zero_grad(set_to_none=True)
            loss, _ = predict_and_loss(model, case, norm, source, target, history, args, True)
            if not torch.isfinite(loss):
                raise ValueError(f'nonfinite train loss for {case.sample}')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            total += float(loss.detach())
        model.eval()
        with torch.no_grad():
            val = np.mean([float(predict_and_loss(model, c, norm, source, target,
                                                  history, args, False)[0])
                           for c in groups['val']])
        rows.append({'epoch': epoch+1, 'train_loss': total/len(groups['train']), 'val_masked_mse': val})
        if val < best - 1e-7:
            best = val
            best_state = {k: v.detach().cpu().clone() for k,v in model.state_dict().items()}
            no_improve = 0
        else:
            no_improve += 1
        if no_improve >= args.patience:
            break
    if best_state is None:
        raise ValueError('no finite validation checkpoint')
    folder.mkdir(parents=True, exist_ok=True)
    torch.save({'state_dict': best_state, 'config': model.get_config()}, folder / 'best.pt')
    model.load_state_dict(best_state)
    model.eval()
    pd.DataFrame(rows).to_csv(folder / 'train_history.csv', index=False)
    return model


def _pearson(x, y):
    if len(x) < 3 or np.std(x) <= 1e-12 or np.std(y) <= 1e-12:
        return float('nan')
    return float(np.corrcoef(x, y)[0, 1])


def metric_tables(actual, predicted, names):
    actual = np.asarray(actual, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    assert actual.shape == predicted.shape and actual.shape[1] == len(names)
    records, valid_mae, valid_mse = [], [], []
    for col, name in enumerate(names):
        good = np.isfinite(actual[:,col]) & np.isfinite(predicted[:,col])
        v, p = actual[good,col], predicted[good,col]
        if len(v):
            mae = float(np.mean(abs(v-p)))
            mse = float(np.mean((v-p)**2))
            valid_mae.extend(abs(v-p))
            valid_mse.extend((v-p)**2)
        else:
            mae = mse = float('nan')
        records.append({'feature': name, 'n': int(good.sum()), 'pearson_r': _pearson(v,p),
                        'mae': mae, 'mse': mse})
    detail = pd.DataFrame.from_records(records)
    r = detail.pearson_r.dropna()
    summary = {'n_samples': int(actual.shape[0]), 'n_features': int(actual.shape[1]),
               'n_valid_cells': len(valid_mae), 'n_features_with_r': int(len(r)),
               'pearson_mean': float(r.mean()) if len(r) else None,
               'pearson_median': float(r.median()) if len(r) else None,
               'mae_micro': float(np.mean(valid_mae)) if valid_mae else None,
               'mse_micro': float(np.mean(valid_mse)) if valid_mse else None,
               'mae_macro': float(detail.mae.mean()) if valid_mae else None,
               'mse_macro': float(detail.mse.mean()) if valid_mse else None}
    return detail, summary


def evaluate_arm(model, cases, norm, source, target, use_history, args, folder):
    actual, predicted, identifiers = [], [], []
    delta_actual, delta_pred = [], []
    model.eval()
    with torch.no_grad():
        for c in cases:
            _, output = predict_and_loss(model,c,norm,source,target,use_history,args,False)
            scaled = output.squeeze(0).cpu().numpy()
            pred = scaled * norm.std[target] + norm.mean[target]
            truth = norm.raw(target, c.sample)
            actual.append(truth)
            predicted.append(pred)
            identifiers.append((c.patient,c.sample,c.visit))
            previous = [s for s in c.history if np.isfinite(norm.raw(target,s)).any()]
            if previous:
                previous_value = norm.raw(target,previous[-1])
                delta_actual.append(truth - previous_value)
                delta_pred.append(pred - previous_value)
            else:
                delta_actual.append(np.full_like(truth, np.nan))
                delta_pred.append(np.full_like(truth, np.nan))
    actual, predicted = np.stack(actual), np.stack(predicted)
    features, summary = metric_tables(actual, predicted, norm.names[target])
    features.to_csv(folder / 'feature_metrics.csv', index=False)
    changes, d_summary = metric_tables(np.stack(delta_actual), np.stack(delta_pred), norm.names[target])
    changes.to_csv(folder / 'change_feature_metrics.csv',index=False)
    summary.update({f'change_{k}': v for k,v in d_summary.items()})
    prediction_rows = []
    for idx, (p,s,t) in enumerate(identifiers):
        for j,feat in enumerate(norm.names[target]):
            if np.isfinite(actual[idx,j]):
                prediction_rows.append((p,s,t,feat,actual[idx,j],predicted[idx,j],
                                        delta_actual[idx][j],delta_pred[idx][j]))
    pd.DataFrame(prediction_rows, columns=[
        'patient_id','sample_id','visit_order','feature','true','predicted','true_change','predicted_change'
    ]).to_csv(folder / 'predictions.csv.gz', index=False,compression='gzip')
    (folder/'metrics.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return summary



def evaluate_persistence(cases, norm, target, folder):
    """Last actually measured target at a prior visit, with no future leakage."""
    truth, pred = [], []
    for c in cases:
        prev = [s for s in c.history if np.isfinite(norm.raw(target, s)).any()]
        if not prev:
            continue
        truth.append(norm.raw(target, c.sample))
        pred.append(norm.raw(target, prev[-1]))
    if not truth:
        return None
    folder.mkdir(parents=True, exist_ok=True)
    features, summary = metric_tables(np.stack(truth), np.stack(pred), norm.names[target])
    features.to_csv(folder/'feature_metrics.csv', index=False)
    # Persistence predicts no change; its change error can be computed without a model.
    change_true = np.stack(truth) - np.stack(pred)
    change_pred = np.zeros_like(change_true)
    change_features, change_summary = metric_tables(change_true, change_pred, norm.names[target])
    change_features.to_csv(folder/'change_feature_metrics.csv', index=False)
    summary.update({f'change_{key}': value for key, value in change_summary.items()})
    _write_json(folder/'metrics.json', summary)
    return summary


def _write_json(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def dataset_dir(root, dataset):
    return root / 'processed' / 'external_validation' / dataset


def run_dataset(dataset, args, root, output):
    folder = dataset_dir(root,dataset)
    if dataset == 'su2020' and args.prepare_su2020 and not (folder/'subjects.csv').is_file():
        from omics.su2020 import prepare
        prep_report=prepare(root/'raw'/'su2020'/'table_s1.xlsx', folder)
        if prep_report['status'] != 'canonical_raw_omics_exported':
            raise ValueError('Su2020: no verified patient-assay links; see preparation_report.json')
    subjects, raw = load_canonical(folder)
    results=[]
    for source,target in [('protein','metab'),('metab','protein')]:
        for task in ('current_completion','historical_completion'):
            history = task == 'historical_completion'
            cases = eligible_cases(subjects, raw, source, target, history)
            if len({c.patient for c in cases}) < 8:
                results.append({'dataset':dataset,'source':source,'target':target,
                                'task':task,'arm':'all','status':'skipped_insufficient_paired_patients',
                                'patients':len({c.patient for c in cases})})
                continue
            ids, groups=split_cases(cases,args.seed)
            train_samples=subjects.loc[subjects.patient_id.isin(ids['train']),'sample_id'].tolist()
            norm=Normalizer(raw,train_samples,minimum=args.min_feature_observations)
            taskfolder=output/dataset/f'{source}_to_{target}'/task
            _write_json(taskfolder/'split.json',ids)
            _write_json(taskfolder/'features.json',norm.names)
            arms=(('pulse_full',True,False),('no_history_input',False,False),('no_alignment',True,True)) if history else (('pulse_current',False,False),)
            for arm,use_history,no_align in arms:
                armfolder=taskfolder/arm
                model=fit_arm(groups,norm,source,target,use_history,no_align,args,armfolder)
                metrics=evaluate_arm(model,groups['test'],norm,source,target,use_history,args,armfolder)
                results.append({'dataset':dataset,'source':source,'target':target,
                                'task':task,'arm':arm,'status':'completed',
                                'train_patients':len(ids['train']),'test_patients':len(ids['test']),
                                **metrics})
                print(f'{dataset} {source}->{target} {task} {arm}: '
                      f'r(median)={metrics["pearson_median"]} MAE={metrics["mae_micro"]:.4g}',flush=True)
            if history:
                persistence = evaluate_persistence(groups['test'], norm, target,
                                                   taskfolder/'last_observation')
                if persistence is not None:
                    results.append({'dataset':dataset,'source':source,'target':target,
                                    'task':task,'arm':'last_observation','status':'completed',
                                    'train_patients':len(ids['train']),'test_patients':len(ids['test']),
                                    **persistence})
    return results


def main(argv=None):
    parser=argparse.ArgumentParser(description='Isolated external validation of the original PULSE model')
    parser.add_argument('--data-root',type=Path,default=Path(os.environ.get('PULSE_OMICS_DATA_ROOT',DEFAULT_ROOT)))
    parser.add_argument('--datasets',nargs='+',choices=DATASETS,default=list(DATASETS))
    parser.add_argument('--prepare-su2020',action='store_true')
    parser.add_argument('--inspect-only',action='store_true')
    parser.add_argument('--device',default='cuda' if torch.cuda.is_available() else 'cpu')
    parser.add_argument('--epochs',type=int,default=30)
    parser.add_argument('--patience',type=int,default=8)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--lr',type=float,default=1e-3)
    parser.add_argument('--latent-dim',type=int,default=16)
    parser.add_argument('--hidden-dim',type=int,default=64)
    parser.add_argument('--temporal-model',choices=('recurrent','transformer'),default='recurrent')
    parser.add_argument('--max-history',type=int,default=4)
    parser.add_argument('--min-feature-observations',type=int,default=5)
    parser.add_argument('--context-weight',type=float,default=.25)
    parser.add_argument('--kl-weight',type=float,default=.001)
    parser.add_argument('--align-weight',type=float,default=.01)
    args=parser.parse_args(argv)
    if args.epochs<1 or args.patience<1 or args.max_history<1 or args.min_feature_observations<2:
        parser.error('positive epochs/patience/history and min feature observations >= 2 required')
    if args.lr <= 0 or args.context_weight < 0 or args.kl_weight < 0 or args.align_weight < 0:
        parser.error('learning rate must be positive; loss weights must be nonnegative')
    if args.device.startswith('cuda') and not torch.cuda.is_available():
        parser.error('CUDA requested but not available')
    root=args.data_root.expanduser()
    runid=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output=root/'results'/'pulse_validation'/runid
    if not args.inspect_only:
        output.mkdir(parents=True,exist_ok=False)
        _write_json(output/'config.json', {k: str(v) if isinstance(v,Path) else v
                                          for k,v in vars(args).items()})
    rows=[]
    for dataset in dict.fromkeys(args.datasets):
        try:
            if args.inspect_only:
                folder=dataset_dir(root,dataset)
                subjects,raw=load_canonical(folder)
                info={'dataset':dataset,'subjects':subjects.patient_id.nunique(),
                      'samples':len(subjects),
                      'paired_current':len(eligible_cases(subjects,raw,'protein','metab',False)),
                      'eligible_historical':len(eligible_cases(subjects,raw,'protein','metab',True))}
                print(json.dumps(info),flush=True)
                continue
            rows.extend(run_dataset(dataset,args,root,output))
        except (FileNotFoundError, ValueError) as exc:
            msg=f'{dataset}: {type(exc).__name__}: {exc}'
            print('SKIPPED '+msg,file=sys.stderr,flush=True)
            if not args.inspect_only:
                rows.append({'dataset':dataset,'task':'all','arm':'all',
                             'status':'skipped','reason':msg})
    if args.inspect_only:
        return 0
    table=pd.DataFrame(rows)
    table.to_csv(output/'summary.csv',index=False)
    finished=bool(len(table) and (table.status=='completed').any())
    _write_json(output/'status.json',{'completed_arms':int((table.status=='completed').sum()) if len(table) else 0,
                                     'run_utc':runid,'note':'independent external validation, not same-cohort paper reproduction'})
    print(f'Summary: {output}/summary.csv',flush=True)
    return 0 if finished else 2


if __name__ == '__main__':
    raise SystemExit(main())