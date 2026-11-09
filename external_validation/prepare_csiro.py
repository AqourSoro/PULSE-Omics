"""Prepare exact source-linked CHPV proteomics/metabolomics longitudinal visits.

Uses the published 005 protein sample-identity workbook plus 006 CCM
metabolomics identifiers and the three 002-004 proteomics workbooks.
No row-position joins; ambiguous identities never become model input.
"""
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

from external_validation.cohort_join_audit import (_find_excel, _excel_sheet,
                                                   _vac_id, _canonical_text,
                                                   audit_csiro)

ROOT=Path('/mnt/f/datasets/pulse-omics')


# Ordinal study stages. These are NOT elapsed weeks or patient-specific dates.
# The published 006 "Times" field contains these four named stages.
CSIRO_VISIT_ORDER = {
    'baseline': 0.,
    '2wk 2nd dose': 1.,
    '6mo 2nd dose': 2.,
    '2wk 3rd dose': 3.,
}


def _visits(values: pd.Series) -> pd.Series:
    def decode(value):
        label = _canonical_text(value)
        if label is None:
            return np.nan
        label = ' '.join(label.lower().split())
        if label in CSIRO_VISIT_ORDER:
            return CSIRO_VISIT_ORDER[label]
        m = re.fullmatch(r'T(?:ime)?\s*([0-9]+)', label, re.I)
        if m:
            return float(m.group(1))
        try:
            return float(label)
        except ValueError:
            return np.nan

    order = values.map(decode)
    if order.isna().any() or not np.isfinite(order).all():
        labels = sorted(set(str(v) for v in values[order.isna()].dropna()))
        raise ValueError(f'Unknown native CCM visit stage(s): {labels[:12]}')
    return order.astype(float)


def _proteins(folder: Path) -> pd.DataFrame:
    by_group=[]
    for prefix in ('002_Task','003_Task','004_Task'):
        frame=_excel_sheet(_find_excel(folder,prefix),'Proteins')
        if 'Accession' not in frame:
            raise ValueError(f'{prefix}: missing native protein accession')
        names=frame['Accession'].map(_canonical_text)
        if names.isna().any() or names.duplicated().any():
            raise ValueError(f'{prefix}: ambiguous duplicate or empty protein accession')
        columns=[c for c in frame.columns if str(c).startswith('Abundance: Vac_Plasma_')]
        if not columns: raise ValueError(f'{prefix}: no quantitative Vac_Plasma abundance columns')
        vals=frame.loc[:,columns].apply(pd.to_numeric,errors='coerce')
        for col in columns:
            old=frame[col].astype('string').str.strip()
            missing=old.isna() | old.isin(['','NA','N/A','NaN','-'])
            if (vals[col].isna() & ~missing).any():
                raise ValueError(f'{prefix}: nonnumeric abundance {col}')
        vals.index=names
        vals=vals.T
        vals.index=[_vac_id(s) for s in vals.index]
        if any(s is None for s in vals.index) or vals.index.duplicated().any():
            raise ValueError(f'{prefix}: invalid or duplicate sample names')
        by_group.append(vals)
    joined=pd.concat(by_group,axis=0,join='outer',sort=True)
    if joined.index.duplicated().any():
        raise ValueError('Vac_Plasma sample ID reused between separate published protein files')
    return joined


def prepare(root: Path, *, inspect_only: bool=False) -> dict:
    folder=root/'raw/csiro'
    audit=audit_csiro(root)
    link=_excel_sheet(_find_excel(folder,'005_Task'),'Sheet1')
    metab=_excel_sheet(_find_excel(folder,'006_Task'),'COVID HPV CCM')
    link=link.dropna(subset=['Sample IDs','proteomic#'])
    barcode='Sample IDS (Barcode)'
    other='IDs & Treatments'
    overlap={name:audit['crosswalk_sample_id_to_metabolomics'][name]['matching_distinct']
             for name in (barcode,other)}
    if overlap[barcode] > 0:
        identity=barcode
    elif overlap[other] > 0:
        identity=other
    else:
        return {'dataset':'csiro','status':'needs_verified_crosswalk',
                'reason':'No exact published 005 Sample IDs to 006 sample identity join',
                'overlap':overlap}
    link=link.assign(protein_id=link['proteomic#'].map(_vac_id))
    if link['protein_id'].isna().any():
        return {'dataset':'csiro','status':'needs_verified_crosswalk',
                'reason':'Some proteomic# values cannot be matched to native Vac_Plasma_ IDs'}
    for key in ['Sample IDs','protein_id']:
        if link[key].duplicated().any():
            raise ValueError(f'005 native crosswalk has duplicate {key}')
    link['_join_sample'] = link['Sample IDs'].map(_canonical_text)
    if link['_join_sample'].isna().any():
        raise ValueError('005 native sample ID is empty')
    metab = metab.copy()
    metab['_join_sample'] = metab[identity].map(_canonical_text)
    # 006 may include QC or blank rows outside the mapped 005 specimens.
    # Discard unmatched source records, never assume a row number is an ID.
    blank_006 = int(metab['_join_sample'].isna().sum())
    candidates = metab.loc[metab['_join_sample'].isin(link['_join_sample'])].copy()
    duplicate_ids = set(candidates.loc[candidates['_join_sample'].duplicated(keep=False), '_join_sample'])
    ambiguous_006 = int(candidates['_join_sample'].isin(duplicate_ids).sum())
    candidates = candidates.loc[~candidates['_join_sample'].isin(duplicate_ids)]
    joined = link.merge(candidates, on='_join_sample', how='inner',
                        validate='one_to_one', suffixes=('_mapping',''))
    protein=_proteins(folder)
    joined=joined.loc[joined.protein_id.isin(protein.index)].copy()
    if joined.empty:
        return {'dataset':'csiro','status':'needs_verified_crosswalk',
                'reason':'No exact 005 proteomic# names match 002-004 abundance headers'}
    joined['patient_id'] = joined['Subject ID'].map(_canonical_text)
    missing_patient = int(joined.patient_id.isna().sum())
    joined = joined.loc[joined.patient_id.notna()].copy()
    joined['visit_order'] = _visits(joined['Times'])
    # Multiple specimens at one patient/stage are not unambiguously a single
    # target observation. Retain neither unless the source provides a rule.
    duplicated_visits = joined.duplicated(['patient_id', 'visit_order'], keep=False)
    excluded_visits = int(duplicated_visits.sum())
    joined = joined.loc[~duplicated_visits].copy()
    if joined.empty:
        return {'dataset':'csiro', 'status':'needs_verified_crosswalk',
                'reason':'No unique patient-stage paired specimens after identity checks',
                'empty_native_006_ids':blank_006,
                'ambiguous_matched_006_rows':ambiguous_006,
                'excluded_duplicate_patient_stage_rows':excluded_visits}
    sample_ids=['csiro_'+str(s) for s in joined['_join_sample']]
    if len(sample_ids)!=len(set(sample_ids)):raise ValueError('non-unique canonical IDs')
    join_counts=joined.patient_id.value_counts()
    report={'dataset':'csiro','status':'inspected',
            'source':'CSIRO 002-006 CHPV workbooks; identity by 005 Sample IDs and 006 '+identity,
            'crosswalk_entries':len(link),'verified_paired_samples':len(joined),
            'empty_native_006_ids':blank_006,
            'ambiguous_matched_006_rows':ambiguous_006,
            'missing_matched_subject_ids':missing_patient,
            'excluded_duplicate_patient_stage_rows':excluded_visits,
            'visit_encoding':'published stages: baseline=0, 2wk 2nd dose=1, 6mo 2nd dose=2, 2wk 3rd dose=3 (ordinal)',
            'paired_subjects':int(joined.patient_id.nunique()),
            'subjects_with_two_or_more_pairs':int((join_counts>=2).sum()),
            'eligible_historical':int((join_counts-1).clip(lower=0).sum()),
            'protein_features':int(protein.shape[1]),
            'metab_features':int(len(metab.columns)-9),
            'join_key':identity}
    if inspect_only:return report
    if report['paired_subjects']<8:
        return {**report,'status':'insufficient_paired_subjects'}
    excluded={'Data File','Subject ID','Sample IDS (Barcode)','Times','IDs & Treatments',
              'COR','Year of Birth','Vaccine Received','_join_sample'}
    fields=[c for c in metab.columns if c not in excluded]
    target=joined.loc[:,fields].apply(pd.to_numeric, errors='coerce')
    for field in fields:
        orig=joined[field].astype('string').str.strip()
        acceptable=orig.isna() | orig.isin(['','NA','N/A','NaN','-'])
        if (target[field].isna() & ~acceptable).any():
            raise ValueError(f'006 CCM metabolomics nonnumeric feature {field}')
    if not np.isfinite(target.to_numpy(dtype=float)[~target.isna().to_numpy()]).all():
        raise ValueError('Nonfinite CCM values')
    prot=protein.loc[joined.protein_id].reset_index(drop=True)
    prot.columns=['protein:'+str(c) for c in prot.columns]
    prot.insert(0,'sample_id',sample_ids)
    target.insert(0,'sample_id',sample_ids)
    subjects=pd.DataFrame({'sample_id':sample_ids,'patient_id':joined.patient_id.values,
                           'visit_order':joined.visit_order.values})
    if not prot.drop(columns='sample_id').notna().any(axis=1).all() or not target.drop(columns='sample_id').notna().any(axis=1).all():
        raise ValueError('One or more crosswalked visits contains no observed assay values')
    out=root/'processed/external_validation/csiro';out.mkdir(parents=True,exist_ok=True)
    staged={name:out/(name+'.building') for name in ('subjects.csv','proteomics.csv','metabolomics.csv')}
    try:
        subjects.to_csv(staged['subjects.csv'],index=False)
        prot.to_csv(staged['proteomics.csv'],index=False)
        target.to_csv(staged['metabolomics.csv'],index=False)
        for name,path in staged.items():path.replace(out/name)
        report['status']='canonical_raw_omics_exported'
        (out/'preparation_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    finally:
        for path in staged.values():path.unlink(missing_ok=True)
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=Path(os.getenv('PULSE_OMICS_DATA_ROOT',ROOT)))
    p.add_argument('--inspect-only',action='store_true')
    args=p.parse_args(argv)
    try: report=prepare(args.data_root,inspect_only=args.inspect_only)
    except (OSError,ValueError,KeyError) as exc:
        print(f'CSIRO preparation error: {type(exc).__name__}: {exc}')
        return 2
    print(json.dumps(report,ensure_ascii=False,indent=2))
    return 0 if report['status'] in ('inspected','canonical_raw_omics_exported') else 2


if __name__=='__main__':raise SystemExit(main())