"""Prepare iPOP plasma omics from its published paired RData tables.

Sample-to-annotation correspondence in the published release is positional
when native row labels differ, as used by netOmics data-raw/hmp_T2D_raw.R.
That provenance is recorded, not silently promoted to native specimen IDs.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from external_validation.scripts.inspect_restored_omics import read_ipop

DEFAULT_ROOT = Path('/mnt/f/datasets/pulse-omics')
SOURCES = {
    'protein': ('Proteomics_df_Data', 'Proteomics_annotation_colData', 'proteomics.csv'),
    'metab': ('Metabolomics_df_Data', 'Metabolomics_annotation_colData', 'metabolomics.csv'),
}


def _date(value):
    """Decode R POSIXct; fix year 00xx only when documented two-digit-year pathology applies."""
    if pd.isna(value):
        return None, None
    if isinstance(value, (int, float, np.integer, np.floating)):
        if not np.isfinite(value):
            return None, None
        d = datetime(1970, 1, 1) + timedelta(seconds=float(value))
    else:
        d = pd.Timestamp(value).to_pydatetime()
    repaired = False
    if 10 <= d.year <= 30:
        d = d.replace(year=d.year + 2000)
        repaired = True
    if not 2010 <= d.year <= 2023:
        raise ValueError(f'iPOP Date decodes to implausible year {d.year}; no unverified calendar reconstruction')
    return d, repaired


def _assay(objects, name, *, permit_release_row_order):
    matrix_name, metadata_name, _ = SOURCES[name]
    matrix = objects[matrix_name]
    metadata = objects[metadata_name]
    if len(matrix) != len(metadata) or not matrix.index.is_unique or not metadata.index.is_unique:
        raise ValueError(f'{name}: assay and annotation lengths/indices are ambiguous')
    for key in ('SubjectID', 'Date', 'Time'):
        if key not in metadata:
            raise ValueError(f'{name}: missing native annotation {key}')
    if not matrix.columns.is_unique or matrix.columns.isna().any():
        raise ValueError(f'{name}: duplicate or empty feature names')
    if not pd.Index(matrix.index.astype(str)).is_unique or not pd.Index(metadata.index.astype(str)).is_unique:
        raise ValueError(f'{name}: row identifiers collide after string conversion')
    if matrix.index.astype(str).equals(metadata.index.astype(str)):
        ordered = metadata.reset_index(drop=True)
        join_mode = 'native_row_index'
    else:
        if not permit_release_row_order:
            raise ValueError(f'{name}: native assay row labels differ from annotation labels; '
                             'source-release parallel row order must be explicitly enabled')
        ordered = metadata.reset_index(drop=True)
        join_mode = 'published_parallel_row_order'
    dates = [_date(x) for x in ordered['Date']]
    repaired = sum(bool(r) for d, r in dates if d is not None)
    day_disagreement = 0
    for (d, _), day in zip(dates, ordered['Time']):
        if d is not None and pd.notna(day):
            day_disagreement += int(int(float(day)) != d.timetuple().tm_yday)
    if day_disagreement:
        raise ValueError(f'{name}: {day_disagreement} native day-of-year / calendar date conflicts')
    vals = matrix.apply(pd.to_numeric, errors='raise').astype('float64').reset_index(drop=True)
    if np.isinf(vals.to_numpy()).any():
        raise ValueError(f'{name}: infinite assay values')
    keys = pd.DataFrame({
        'patient_id': ordered['SubjectID'].astype('string').str.strip(),
        'date': [d.date().isoformat() if d is not None else None for d, _ in dates],
    })
    good = (keys.patient_id.notna() & keys.patient_id.ne('') & keys.date.notna()
            & vals.notna().any(axis=1))
    keys = keys.loc[good].copy()
    vals = vals.loc[good].copy()
    # Duplicate same-subject/same-day records cannot be aligned across modalities.
    duplicated = keys.duplicated(['patient_id','date'],keep=False)
    excluded = int(duplicated.sum())
    keys = keys.loc[~duplicated].copy()
    vals = vals.loc[keys.index].copy()
    keys['row'] = np.arange(len(keys))
    vals.index = keys['row']
    keys.index = keys['row']
    return keys, vals, {
        'raw_rows': len(matrix), 'candidate_rows': int(good.sum()),
        'excluded_duplicate_patient_date_rows': excluded,
        'rows_ready': len(keys), 'features': vals.shape[1],
        'join_mode': join_mode, 'two_digit_year_repairs': repaired,
        'nan_fraction': round(float(vals.isna().to_numpy().mean()),6),
        'source_processing': 'published processed data; metabolite values may already be KNN-imputed',
    }


def prepare(root: Path, *, inspect_only=False, permit_release_row_order=False, reader=None):
    source = root/'raw/ipop/seasonal/Multi_Omics_Seasonal.RData'
    if not source.is_file():
        raise FileNotFoundError(source)
    objects = (reader or read_ipop)(source)
    pkeys,pvals,preport = _assay(objects,'protein',permit_release_row_order=permit_release_row_order)
    mkeys,mvals,mreport = _assay(objects,'metab',permit_release_row_order=permit_release_row_order)
    paired = pkeys.merge(mkeys,on=['patient_id','date'],validate='one_to_one',suffixes=('_p','_m'))
    if paired.empty:
        raise ValueError('iPOP: no shared verified participant/calendar-day records')
    paired = paired.sort_values(['patient_id','date']).reset_index(drop=True)
    counts = paired.groupby('patient_id').size()
    paired['visit_order'] = pd.to_datetime(paired['date']).sub(
        pd.to_datetime(paired.groupby('patient_id')['date'].transform('min'))).dt.days.astype(float)
    paired['sample_id'] = [hashlib.sha256(f'{pid}|{date}'.encode()).hexdigest()[:20]
                           for pid,date in zip(paired.patient_id,paired.date)]
    report = {
        'dataset':'ipop','source':str(source),
        'status':'inspected' if inspect_only else 'canonical_raw_omics_exported',
        'protein':preport,'metab':mreport,
        'paired_samples':len(paired), 'paired_subjects':len(counts),
        'subjects_with_two_or_more_pairs':int((counts>=2).sum()),
        'eligible_historical':int((counts-1).clip(lower=0).sum()),
        'time_units':'real elapsed calendar days per subject, not seasonal day-of-year',
        'identity_limitation':('Rows paired by publication-conventional parallel order, '
                               'not an independently supplied specimen crosswalk')
                               if any(r['join_mode']=='published_parallel_row_order'
                                      for r in (preport,mreport)) else None,
    }
    if inspect_only:
        return report
    if len(counts)<8 or (counts>=2).sum()<8:
        raise ValueError('iPOP: insufficient longitudinal paired participants')
    out = root/'processed/external_validation/ipop'
    out.mkdir(parents=True,exist_ok=True)
    subjects = paired[['sample_id','patient_id','visit_order','date']]
    protein = pvals.loc[paired.row_p].reset_index(drop=True)
    metab = mvals.loc[paired.row_m].reset_index(drop=True)
    protein.insert(0,'sample_id',paired.sample_id.tolist())
    metab.insert(0,'sample_id',paired.sample_id.tolist())
    staged=[]
    try:
        for filename,frame in [('subjects.csv',subjects),('proteomics.csv',protein),('metabolomics.csv',metab)]:
            temp = out/(filename+'.building')
            frame.to_csv(temp,index=False,na_rep='')
            staged.append((temp,out/filename))
        for temp,target in staged:
            temp.replace(target)
        (out/'preparation_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    finally:
        for temp,_ in staged:
            temp.unlink(missing_ok=True)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=DEFAULT_ROOT)
    p.add_argument('--inspect-only',action='store_true')
    p.add_argument('--permit-release-row-order',action='store_true',
                   help='Use documented parallel RData row order when native indexes differ')
    a=p.parse_args()
    print(json.dumps(prepare(a.data_root,inspect_only=a.inspect_only,
                             permit_release_row_order=a.permit_release_row_order),
                     indent=2,ensure_ascii=False))

if __name__=='__main__': main()
