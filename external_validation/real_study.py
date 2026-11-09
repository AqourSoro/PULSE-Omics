"""Matched-patient experiments on actual paired omics with the original PulseModel."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from external_validation import run as legacy
from pulse.model import PulseModel

COHORTS = ('su2020', 'csiro', 'ipop', 'hmp2')
ARMS = ('current', 'history', 'no_alignment', 'no_chronology', 'forecast')
MODS = legacy.MODALITIES
DEFAULT_ROOT = Path('/mnt/f/datasets/pulse-omics')
SU_TABLE_S1 = ('https://data.mendeley.com/public-files/datasets/tzydswhhb5/files/'
               'ead2d262-fb88-4bc7-b697-9b0f2bc7d4db/file_downloaded')
SU_TABLE_S1_SHA256 = '1082603fe92429d2d115c6ee9e2ed8223f8ce0699a95669347fca8307dbc7744'


def fetch_missing(root, dataset):
    from external_validation.scripts.download_real_omics import Downloader, HMP2, IPOP_SEASONAL_URL, checksum

    d = Downloader(root,False,.5,1.)
    raw = root/'raw'/dataset
    if dataset == 'su2020':
        file = raw/'table_s1.xlsx'
        if not file.is_file() and not (raw/'pmc/mmc1.xlsx').is_file():
            d.https(dataset,file.name,SU_TABLE_S1)
            if file.is_file() and checksum(file) != SU_TABLE_S1_SHA256:
                raise ValueError('Su2020 Table S1 differs from the Mendeley v5 checksum')
    elif dataset == 'hmp2':
        for name,url in HMP2.items():
            if not (raw/name).is_file():
                d.https(dataset,name,url)
    elif dataset == 'ipop':
        name = 'seasonal/Multi_Omics_Seasonal.RData'
        if not (raw/name).is_file():
            d.https(dataset,name,IPOP_SEASONAL_URL)
    elif dataset == 'csiro':
        if not all(any(raw.glob(f'{n:03d}_*.xlsx')) for n in range(2,7)):
            d.csiro()
    elif dataset == 'covaar':
        from external_validation.scripts.download_real_omics import COVAAR
        for name, url in COVAAR.items():
            if not (raw/'geo'/name).is_file():
                if 'RFUs' in name:
                    url = f'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE325178&format=file&file={name}'
                d.https(dataset,f'geo/{name}',url)
        figshare = {'Metabolite_Meta_Data.csv': 62724259, 'Metabolomics_Data.csv': 62724262,
                    'Sample_Meta_Data.csv': 62724265}
        for name, file_id in figshare.items():
            if not (raw/'figshare'/name).is_file():
                d.https(dataset,f'figshare/{name}',f'https://ndownloader.figshare.com/files/{file_id}')
    if dataset in d.report:
        d.flush(dataset)
    if d.failed:
        raise RuntimeError(f'{dataset}: download failed, see raw/{dataset}/download_manifest.json')


def prepare(root, dataset, download=False):
    directory = legacy.dataset_dir(root, dataset)
    if all((directory / f).is_file() for f in ('subjects.csv','proteomics.csv','metabolomics.csv')):
        return
    if download:
        fetch_missing(root,dataset)
    if dataset == 'su2020':
        from external_validation.omics.su2020 import prepare as extract
        source = root / 'raw/su2020/table_s1.xlsx'
        if not source.is_file():
            source = root / 'raw/su2020/pmc/mmc1.xlsx'
        extract(source, directory)
    elif dataset == 'csiro':
        from external_validation.prepare_csiro import prepare as extract
        extract(root)
    elif dataset == 'ipop':
        from external_validation.prepare_ipop import prepare as extract
        extract(root, permit_release_row_order=True)
    elif dataset == 'covaar':
        from external_validation.prepare_covaar import prepare as extract
        extract(root)
    elif dataset == 'hmp2':
        from external_validation.prepare_hmp2 import prepare as extract
        extract(root)
    if not all((directory / f).is_file() for f in ('subjects.csv','proteomics.csv','metabolomics.csv')):
        raise RuntimeError(f'{dataset}: source preparation did not produce all three linked tables')


class Observations:
    def __init__(self, raw, norm, device):
        self.row = {str(sample): i for i, sample in enumerate(raw['protein'].index)}
        self.x, self.valid, self.available = {}, {}, {}
        for mod in MODS:
            values = norm.values[mod].to_numpy(dtype=np.float32)
            valid = np.isfinite(values)
            scaled = (values - norm.mean[mod]) / norm.std[mod]
            self.available[mod] = valid.any(axis=1)
            self.valid[mod] = torch.from_numpy(valid).to(device)
            self.x[mod] = torch.from_numpy(np.where(valid, scaled, 0).astype(np.float32)).to(device)

    def usable(self, case):
        i = self.row[case.sample]
        return all(self.available[m][i] for m in MODS) and any(
            all(self.available[m][self.row[s]] for m in MODS) for s in case.history
        )

    def batches(self, cases, source, arm, batch_size, max_history, rng=None, *, reverse=False,
                past_source_only=False):
        target = 'metab' if source == 'protein' else 'protein'
        groups = defaultdict(list)
        order = rng.permutation(len(cases)) if rng is not None else range(len(cases))
        for idx in order:
            case = cases[idx]
            history = [self.row[s] for s in case.history if any(self.available[m][self.row[s]] for m in MODS)]
            history = history[-max_history:] if arm != 'current' else []
            if arm == 'no_chronology':
                # Independent orders at each epoch. Evaluation uses a fixed RNG supplied by caller.
                history = rng.permutation(history).tolist()
            if reverse:
                history = history[::-1]
            availability = tuple(tuple(bool(self.available[m][i] and (m == source or not past_source_only))
                                       for m in MODS) for i in history)
            # Ignore historical visits containing no remaining input assay.
            history = [i for i, ok in zip(history, availability) if any(ok)]
            availability = tuple(ok for ok in availability if any(ok))
            groups[availability].append((int(idx), history + [self.row[case.sample]]))
        chunks = []
        for pattern, entries in groups.items():
            chunks.extend((pattern, entries[i:i+batch_size]) for i in range(0,len(entries),batch_size))
        if rng is not None:
            rng.shuffle(chunks)
        for pattern, entries in chunks:
            ids, sequences = zip(*entries)
            indices = np.asarray(sequences, dtype=int).T
            visits, masks, context = [], [], []
            current = tuple(m == source and arm != 'forecast' for m in MODS)
            for row_ids, observed in zip(indices, (*pattern, current)):
                data, mask, valid = {}, {}, {}
                for mod, known in zip(MODS, observed):
                    data[mod] = self.x[mod][row_ids] if known else None
                    mask[mod] = 0 if known else 2
                    valid[mod] = self.valid[mod][row_ids] if known else None
                visits.append(data)
                masks.append(mask)
                context.append(valid)
            masks[-1][target] = 1
            yield np.asarray(ids), visits, masks, context, self.x[target][indices[-1]], self.valid[target][indices[-1]]


def risk(output, truth, valid):
    error = output - truth
    count = valid.sum(1).clamp_min(1)
    return (error.abs()*valid).sum(1)/count, (error.square()*valid).sum(1)/count


def objective(output, visits, masks, context, target, truth, valid, align):
    _, sq = risk(output['reconstructions'][-1][target], truth, valid)
    context_loss = truth.new_zeros(())
    for mod in MODS:
        terms = []
        for rec, visit, mask, observed in zip(output['reconstructions'],visits,masks,context):
            if mask[mod] == 0:
                _, mse = risk(rec[mod],visit[mod],observed[mod])
                terms.append(mse.mean())
        if terms:
            context_loss = context_loss + torch.stack(terms).mean()
    losses = output['losses']
    return sq.mean() + .25*context_loss + .001*losses['kl_loss'] + align*losses['alignment_loss']


def predict(model, obs, cases, source, arm, args, seed, **controls):
    target = 'metab' if source == 'protein' else 'protein'
    pred = np.empty((len(cases), obs.x[target].shape[1]), dtype=np.float32)
    rng = np.random.default_rng(seed) if arm == 'no_chronology' else None
    model.eval()
    with torch.no_grad():
        for ids, visits, masks, _, _, _ in obs.batches(cases,source,arm,args.batch_size,args.max_history,rng,
                                                       reverse=controls.get('reverse',False),
                                                       past_source_only=controls.get('past_source_only',False)):
            out = model(visits,masks,use_history_in_fusion=controls.get('fusion',True))
            pred[ids] = out['reconstructions'][-1][target].cpu().numpy()
    return pred


def train(obs, groups, source, arm, args, seed, output):
    torch.manual_seed(seed)
    model = PulseModel(modality_dims={m:obs.x[m].shape[1] for m in MODS}, latent_dim=16, hidden_dim=64,
                       lambda_kl=0, lambda_align=0, lambda_adv=0, use_adversarial_loss=False,
                       temporal_model='recurrent').to(args.device)
    if args.resume and (output/'best.pt').is_file() and (output/'epochs.csv').is_file():
        state = torch.load(output/'best.pt',map_location=args.device,weights_only=True)
        model.load_state_dict(state['state_dict'])
        return model
    target = 'metab' if source == 'protein' else 'protein'
    optimizer = torch.optim.AdamW(model.parameters(), lr=.001, weight_decay=1e-4)
    val_keys = [obs.row[c.sample] for c in groups['val']]
    truth = obs.x[target][val_keys].cpu().numpy()
    valid = obs.valid[target][val_keys].cpu().numpy()
    best, stagnant, checkpoint = float('inf'), 0, None
    history = []
    for epoch in range(args.epochs):
        model.train()
        rng = np.random.default_rng(seed+epoch)
        losses = []
        for _, visits, masks, context, y, keep in obs.batches(groups['train'],source,arm,args.batch_size,args.max_history,rng):
            optimizer.zero_grad(set_to_none=True)
            result = model(visits,masks)
            loss = objective(result,visits,masks,context,target,y,keep,0 if arm=='no_alignment' else .01)
            if not torch.isfinite(loss):
                raise FloatingPointError(f'{source} {arm}: nonfinite training objective')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),5.)
            optimizer.step()
            losses.append(float(loss.detach()))
        pred = predict(model,obs,groups['val'],source,arm,args,seed)
        mse = float((((pred-truth)**2*valid).sum(1)/valid.sum(1)).mean())
        history.append({'epoch':epoch+1,'train_objective':float(np.mean(losses)),'val_zmse':mse})
        if mse < best-1e-7:
            best, stagnant = mse, 0
            checkpoint = {k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
        else:
            stagnant += 1
        if stagnant >= args.patience:
            break
    output.mkdir(parents=True,exist_ok=True)
    staged = output/'best.pt.part'
    torch.save({'state_dict':checkpoint,'config':model.get_config()},staged)
    staged.replace(output/'best.pt')
    pd.DataFrame(history).to_csv(output/'epochs.csv',index=False)
    model.load_state_dict(checkpoint)
    return model


def median_correlation(x,y,valid):
    n = valid.sum(0)
    x = np.where(valid,x,0).astype(np.float64)
    y = np.where(valid,y,0).astype(np.float64)
    x = np.where(valid,x-x.sum(0)/n.clip(1),0)
    y = np.where(valid,y-y.sum(0)/n.clip(1),0)
    scale = np.sqrt((x*x).sum(0)*(y*y).sum(0))
    good = (n>=3) & (scale>1e-12)
    return float(np.median((x*y).sum(0)[good]/scale[good])) if good.any() else None


def metrics(cases, truth, pred, valid):
    patients = pd.Series([c.patient for c in cases])
    abs_err = np.abs(pred-truth)
    visit_loss = (abs_err*valid).sum(1)/valid.sum(1)
    per_patient = pd.Series(visit_loss).groupby(patients).mean()
    a,b = [pd.DataFrame(np.where(valid,x,np.nan)) for x in (truth,pred)]
    a,b = a-a.groupby(patients).transform('mean'), b-b.groupby(patients).transform('mean')
    return {'zmae':float(per_patient.mean()), 'median_feature_r':median_correlation(truth,pred,valid),
            'median_within_subject_r':median_correlation(a.to_numpy(),b.to_numpy(),valid),
            'patients':len(per_patient),'visits':len(cases),'observed_cells':int(valid.sum())}, per_patient


def record_scores(cases, truth, valid, values, output):
    rows, risks = [], {}
    for name, pred in values.items():
        result, risks[name] = metrics(cases,truth,pred,valid)
        rows.append({'arm':name,**result})
    rng = np.random.default_rng(772)
    contrasts = []
    for name, ref in (('history','current'),('history','no_alignment'),('history','no_chronology'),
                     ('history','persistence'),('forecast','persistence'),('forecast','mean'),('current','mean')):
        if name in risks and ref in risks:
            delta = (risks[name]-risks[ref]).to_numpy()
            ci = np.quantile(rng.choice(delta,(5000,len(delta)),replace=True).mean(1),[.025,.975])
            contrasts.append({'test':name,'reference':ref,'delta_zmae':float(delta.mean()),'ci_low':ci[0],'ci_high':ci[1]})
    pd.DataFrame(rows).to_csv(output/'metrics.csv',index=False)
    pd.DataFrame(contrasts).to_csv(output/'contrasts.csv',index=False)
    return rows,contrasts


def study_dataset(root,dataset,args,output):
    prepare(root,dataset,args.fetch_missing)
    source_dir = legacy.dataset_dir(root,dataset)
    signatures = {}
    for name in ('subjects.csv','proteomics.csv','metabolomics.csv'):
        path = source_dir/name
        with path.open('rb') as handle:
            signatures[name] = {'bytes':path.stat().st_size,'sha256':hashlib.file_digest(handle,'sha256').hexdigest()}
    signature_file = output/dataset/'input_files.json'
    if signature_file.exists() and json.loads(signature_file.read_text()) != signatures:
        raise ValueError(f'{dataset}: prepared input changed after the run was created')
    legacy._write_json(signature_file,signatures)
    subjects,raw = legacy.load_canonical(source_dir)
    all_cases = legacy.eligible_cases(subjects,raw,'protein','metab',True)
    visit = subjects.set_index('sample_id').visit_order.to_dict()
    outputs = []
    for seed in args.seeds:
        ids,groups = legacy.split_cases(all_cases,seed)
        samples = subjects.loc[subjects.patient_id.isin(ids['train']),'sample_id'].tolist()
        norm = legacy.Normalizer(raw,samples,**legacy.feature_policy(dataset))
        obs = Observations(raw,norm,args.device)
        groups = {name:[case for case in part if obs.usable(case)] for name,part in groups.items()}
        if not all(groups.values()):
            raise ValueError(f'{dataset}, seed {seed}: no eligible case in a patient split after feature selection')
        cases = groups['test']
        target_keys = [obs.row[c.sample] for c in cases]
        for source,target in (('protein','metab'),('metab','protein')):
            folder = output / dataset / f'seed-{seed}' / f'{source}_to_{target}'
            folder.mkdir(parents=True,exist_ok=True)
            if args.resume and (folder/'heldout.npz').is_file() and (folder/'metrics.csv').is_file():
                scores = pd.read_csv(folder/'metrics.csv').to_dict('records')
                outputs.extend({'dataset':dataset,'seed':seed,'source':source,'target':target,**row} for row in scores)
                continue
            legacy._write_json(folder/'split.json',ids)
            legacy._write_json(folder/'features.json',norm.names)
            legacy._write_json(folder/'feature_selection.json',norm.selection)
            truth,valid = obs.x[target][target_keys].cpu().numpy(),obs.valid[target][target_keys].cpu().numpy()
            # Persistence uses the last available measurement for each feature, with train mean as fallback.
            last = np.zeros_like(truth)
            for i,c in enumerate(cases):
                for sample in c.history:
                    row = obs.row[sample]
                    known = obs.valid[target][row].cpu().numpy()
                    last[i,known] = obs.x[target][row].cpu().numpy()[known]
            values = {'mean':np.zeros_like(truth),'persistence':last}
            for arm in ARMS:
                if arm == 'no_chronology' and not any(len(c.history)>=2 for c in groups['train']):
                    continue
                model = train(obs,groups,source,arm,args,seed,folder/arm)
                values[arm] = predict(model,obs,cases,source,arm,args,seed)
                if arm=='history':
                    values['fusion_off'] = predict(model,obs,cases,source,arm,args,seed,fusion=False)
                    values['past_source_only'] = predict(model,obs,cases,source,arm,args,seed,past_source_only=True)
                    values['reverse_own'] = predict(model,obs,cases,source,arm,args,seed,reverse=True)
                print(dataset,seed,source,target,arm,'complete',flush=True)
            scores, contrasts = record_scores(cases,truth,valid,values,folder)
            staged = folder/'heldout.part.npz'
            np.savez_compressed(staged, sample_id=[c.sample for c in cases],
                                patient_id=[c.patient for c in cases],
                                gap=[c.visit-visit[c.history[-1]] for c in cases],
                                history_n=[min(len(c.history),args.max_history) for c in cases],
                                truth=truth,valid=valid,**values)
            staged.replace(folder/'heldout.npz')
            outputs.extend({'dataset':dataset,'seed':seed,'source':source,'target':target,**row} for row in scores)
    return outputs


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=Path(os.environ.get('PULSE_OMICS_DATA_ROOT',DEFAULT_ROOT)))
    p.add_argument('--datasets',nargs='+',choices=(*COHORTS, 'covaar'),default=list(COHORTS))
    p.add_argument('--seeds',nargs='+',type=int,default=[42,137,509])
    p.add_argument('--device',default='cuda' if torch.cuda.is_available() else 'cpu')
    p.add_argument('--epochs',type=int,default=50)
    p.add_argument('--patience',type=int,default=10)
    p.add_argument('--batch-size',type=int,default=8)
    p.add_argument('--max-history',type=int,default=4)
    p.add_argument('--output',type=Path)
    p.add_argument('--fetch-missing',action='store_true',help='Fetch public sources only when needed for preparation')
    p.add_argument('--resume',action='store_true',help='Reuse completed arms from the exact saved run configuration')
    a=p.parse_args(argv)
    if min(a.batch_size,a.epochs,a.patience,a.max_history) < 1:
        p.error('batch size, epochs, patience and history length must be positive')
    torch.set_num_threads(1)
    root=a.data_root.expanduser().resolve()
    output=a.output or root/'results/real_study'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    config = {k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items() if k not in ('resume','fetch_missing')}
    config_path = output/'config.json'
    if a.resume:
        saved = json.loads(config_path.read_text())
        if config != {k:v for k,v in saved.items() if k not in ('resume','fetch_missing')}:
            p.error('resume requires the exact saved experiment arguments')
    else:
        output.mkdir(parents=True,exist_ok=False)
        legacy._write_json(config_path,config)
    rows=[]
    for d in dict.fromkeys(a.datasets):
        rows.extend(study_dataset(root,d,a,output))
        pd.DataFrame(rows).to_csv(output/'summary.csv',index=False)
    print(output/'summary.csv',flush=True)


if __name__=='__main__':
    main()
