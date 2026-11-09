"""Data-readiness report: actual matched patients, observability and held-out feature sets.

No training, checkpoint creation, or changes to original publisher data.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np

ROOT=Path('/mnt/f/datasets/pulse-omics')
COHORTS=('su2020','hmp2','ipop','covaar','csiro')


def preflight(root: Path, dataset: str, seed=42, min_feature_observations=5,
              hmp2_min_nonzero_fraction=.10, hmp2_max_features=4096):
    from external_validation.run import load_canonical, eligible_cases, split_cases, Normalizer, feature_policy
    folder=root/'processed'/'external_validation'/dataset
    subjects, raw=load_canonical(folder)
    report={'dataset':dataset,
            'patients':int(subjects.patient_id.nunique()),
            'visits':len(subjects),
            'modalities':{}, 'directions':{}}
    for modality, matrix in raw.items():
        a=matrix.to_numpy(dtype=float)
        observed=np.isfinite(a)
        finite_values = a[observed]
        nonzero = observed & (a != 0)
        report['modalities'][modality]={
            'features':int(a.shape[1]),
            'samples_with_measurements':int(observed.any(axis=1).sum()),
            'measured_cell_fraction':round(float(observed.mean()),5),
            'all_missing_features':int((~observed.any(axis=0)).sum()),
            'zero_fraction_among_finite':round(float((finite_values == 0).mean()),5)
                if finite_values.size else None,
            'all_zero_features':int((observed.any(axis=0) & ~nonzero.any(axis=0)).sum()),
            'features_with_fewer_than_five_nonzero_samples':
                int((nonzero.sum(axis=0) < 5).sum()),
            'sparse_zero_semantics':('HMP2 BIOM sparse implicit zeros are zeros after export;'
                 ' zero is not evidence of a measured, quantified abundance')
                 if dataset=='hmp2' and modality=='metab' else None,
        }
    for source,target in [('protein','metab'),('metab','protein')]:
        key=source+'_to_'+target
        report['directions'][key]={}
        for task,history in [('current_completion',False),('historical_completion',True)]:
            cases=eligible_cases(subjects,raw,source,target,history)
            patient_count=len({c.patient for c in cases})
            rec={'eligible_visits':len(cases),'eligible_patients':patient_count}
            if patient_count>=8:
                ids,groups=split_cases(cases,seed)
                train_samples=subjects.loc[subjects.patient_id.isin(ids['train']),'sample_id'].tolist()
                policy=feature_policy(dataset,hmp2_min_nonzero_fraction=hmp2_min_nonzero_fraction,
                                      hmp2_max_features=hmp2_max_features)
                norm=Normalizer(raw,train_samples,minimum=min_feature_observations,**policy)
                rec['feature_selection']=norm.selection
                rec['split_patients']={k:len(v) for k,v in ids.items()}
                rec['split_visits']={k:len(v) for k,v in groups.items()}
                rec['selected_training_features']={m:len(features) for m,features in norm.names.items()}
                rec['status']='ready_for_training'
            else:
                rec['status']='insufficient_paired_patients'
            report['directions'][key][task]=rec
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,default=Path(os.environ.get('PULSE_OMICS_DATA_ROOT',ROOT)))
    parser.add_argument('--datasets',nargs='+',choices=COHORTS,default=['su2020','hmp2'])
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--min-feature-observations',type=int,default=5)
    parser.add_argument('--hmp2-min-nonzero-fraction',type=float,default=.10)
    parser.add_argument('--hmp2-max-features',type=int,default=4096)
    parser.add_argument('--save',action='store_true',help='Write preparation QC JSON to F drive')
    args=parser.parse_args(argv)
    for ds in dict.fromkeys(args.datasets):
        try:
            report=preflight(args.data_root,ds,args.seed,args.min_feature_observations,
                             args.hmp2_min_nonzero_fraction,args.hmp2_max_features)
        except (ValueError,FileNotFoundError) as e:
            print(f'{ds}: not ready: {e}')
            continue
        print(f'{ds}: {report["patients"]} patients, {report["visits"]} visits')
        for m,meta in report['modalities'].items():
            zero_rate = meta['zero_fraction_among_finite']
            zero_label = f'{zero_rate:.1%}' if zero_rate is not None else 'NA'
            print(f'  {m}: {meta["features"]} features; {meta["samples_with_measurements"]} measured samples; '
                  f'{meta["measured_cell_fraction"]:.1%} finite cells; '
                  f'{zero_label} zeros; '
                  f'{meta["all_zero_features"]} all-zero features; '
                  f'{meta["features_with_fewer_than_five_nonzero_samples"]} features with <5 nonzero samples')
        for direction,tasks in report['directions'].items():
            for name,info in tasks.items():
                print(f'  {direction}/{name}: {info["eligible_visits"]} visits, '
                      f'{info["eligible_patients"]} patients, {info["status"]}')
                if 'selected_training_features' in info:
                    print(f'    train-only selected features: {info["selected_training_features"]}')
        if args.save:
            dest=args.data_root/'processed/external_validation'/ds/'preflight.json'
            dest.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
            print(f'  report: {dest}')
    return 0


if __name__=='__main__':
    raise SystemExit(main())