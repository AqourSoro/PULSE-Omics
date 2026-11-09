"""Prepare a strictly matched HMP2 stool proteomics/metabolomics cohort.

Use only the three published merged files and explicit IBDMDB sample metadata.
The source matrices and the original PULSE implementation remain unchanged.

  python -m external_validation.prepare_hmp2 --inspect-only
  python -m external_validation.prepare_hmp2

Requires `ijson` for the large metabolomics BIOM stream (pip install ijson).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path('/mnt/f/datasets/pulse-omics')
PROTEIN = 'HMP2_proteomics_1pep-pro_1p_FDR_5ppm.tsv.gz'
METAB = 'HMP2_metabolomics_w_metadata.biom.gz'
METADATA = 'hmp2_metadata_2018-08-20.csv'
ID_FIELDS = ('External ID', 'Participant ID', 'week_num')


def _text(s: object) -> str:
    return str(s).strip() if pd.notna(s) else ''


def metadata_map(path: Path) -> tuple[pd.DataFrame, dict]:
    data = pd.read_csv(path, dtype=str, low_memory=False, encoding='utf-8-sig')
    absent = sorted(set(ID_FIELDS) - set(data.columns))
    if absent:
        raise ValueError(f'HMP2 metadata lacks required columns {absent}; actual={list(data.columns)[:30]}')
    data = data.loc[:, list(ID_FIELDS)].copy()
    for col in ID_FIELDS:
        data[col] = data[col].map(_text)
    data = data.loc[(data['External ID'] != '') & (data['Participant ID'] != '')].copy()
    data['week_num'] = pd.to_numeric(data['week_num'], errors='coerce')
    data = data.loc[np.isfinite(data['week_num'])].copy()
    # Metadata can list the same sample more than once, e.g. separate assay pipelines.
    # It must, however, always identify exactly one patient and one visit.
    collisions = data.groupby('External ID')[['Participant ID', 'week_num']].nunique()
    bad = collisions.index[(collisions > 1).any(axis=1)].tolist()
    if bad:
        raise ValueError(f'{len(bad)} sample IDs conflict across patient or week, e.g. {bad[:5]}')
    unique = data.drop_duplicates('External ID').set_index('External ID')
    return unique, {'metadata_rows_valid': int(len(data)),
                    'metadata_sample_ids_unique': int(len(unique)),
                    'metadata_duplicate_rows_consistent': int(len(data)-len(unique))}


def protein_ids(path: Path) -> tuple[str, list[str], int]:
    skipped = 0
    with gzip.open(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        for header in csv.reader(handle, delimiter='\t'):
            if len(header) >= 3:
                break
            skipped += 1
        else:
            raise ValueError('HMP2 protein TSV lacks sample columns')
    first, ids = header[0], [x.strip() for x in header[1:]]
    if '' in ids or len(set(ids)) != len(ids):
        raise ValueError('HMP2 protein TSV has missing or repeated sample column IDs')
    return first, ids, skipped


def biom_metadata(path: Path, small_json_limit: int = 32*1024*1024) -> tuple[list[str], list[str], str, tuple[int, int]]:
    """Read BIOM identifiers with a streaming parser; never load full real matrix."""
    try:
        import ijson
    except ImportError:
        if path.stat().st_size > small_json_limit:
            raise RuntimeError('Run `python -m pip install ijson` to parse the large HMP2 BIOM file') from None
        with gzip.open(path, 'rt', encoding='utf-8') as fp:
            obj = json.load(fp)
        ids = [str(x['id']).strip() for x in obj['columns']]
        features = [str(x['id']).strip() for x in obj['rows']]
        matrix_type = str(obj['matrix_type'])
        shape = tuple(int(v) for v in obj['shape'])
    else:
        ids, features, shape = [], [], []
        matrix_type = None
        with gzip.open(path, 'rb') as fp:
            for prefix, event, value in ijson.parse(fp):
                if event in ('string', 'number') and prefix == 'columns.item.id':
                    ids.append(str(value).strip())
                elif event in ('string', 'number') and prefix == 'rows.item.id':
                    features.append(str(value).strip())
                elif event in ('string', 'number') and prefix == 'shape.item':
                    shape.append(int(value))
                elif event == 'string' and prefix == 'matrix_type':
                    matrix_type = value
                # Stop as soon as all published identifiers are available.
                # There is no need to parse millions of sparse matrix elements.
                if len(shape)==2 and matrix_type is not None and \
                   len(features)==shape[0] and len(ids)==shape[1]:
                    break
        shape = tuple(shape)
    if len(shape) != 2 or shape != (len(features), len(ids)):
        raise ValueError(f'BIOM shape mismatch: metadata rows={len(features)}, cols={len(ids)}, shape={shape}')
    if len(set(ids)) != len(ids) or not all(ids):
        raise ValueError('BIOM column sample IDs are missing or repeated')
    if len(set(features)) != len(features) or not all(features):
        raise ValueError('BIOM feature IDs are missing or repeated')
    if matrix_type not in ('sparse', 'dense'):
        raise ValueError(f'Unsupported BIOM matrix_type {matrix_type}')
    return ids, features, matrix_type, shape


def matched_visits(metadata: pd.DataFrame, protein: list[str], metab: list[str]) -> tuple[pd.DataFrame, dict]:
    p, m = set(protein), set(metab)
    shared = sorted(p & m)
    known = [sample for sample in shared if sample in metadata.index]
    result = metadata.loc[known].rename_axis('External ID').reset_index().rename(columns={
        'External ID':'sample_id', 'Participant ID':'patient_id', 'week_num':'visit_order'})
    duplicates = result.duplicated(['patient_id', 'visit_order'], keep=False)
    conflict_count = int(duplicates.sum())
    if conflict_count:
        # Preserve the mappings in the report; do not select an arbitrary specimen.
        result = result.loc[~duplicates].copy()
    result = result.sort_values(['patient_id','visit_order','sample_id']).reset_index(drop=True)
    counts = result.groupby('patient_id').size()
    report = {
        'protein_sample_ids':len(p), 'metabolite_sample_ids':len(m),
        'protein_ids_with_metadata':len(p & set(metadata.index)),
        'metabolite_ids_with_metadata':len(m & set(metadata.index)),
        'exact_protein_metabolite_shared_ids':len(shared),
        'exact_shared_ids_with_patient_and_week':len(known),
        'rows_excluded_due_to_duplicate_patient_week':conflict_count,
        'canonical_paired_visits':int(len(result)),
        'canonical_subjects':int(result.patient_id.nunique()),
        'subjects_with_two_or_more_paired_visits':int((counts >= 2).sum()),
        'join_basis':'identical sample IDs in protein and BIOM; identical ID mapped to one patient/week in clinical metadata',
    }
    return result, report


def stream_biom_matrix(path: Path, target_ids: list[str], all_ids: list[str], nfeatures: int,
                       matrix_type: str, dest: Path) -> np.memmap:
    """Create temporary sample x feature float32 matrix for paired specimens only."""
    try:
        import ijson
    except ImportError:
        ijson = None
    index = {i:j for j,i in enumerate(all_ids)}
    col_out = {index[s]:pos for pos,s in enumerate(target_ids)}
    x = np.memmap(dest, dtype='float32', mode='w+', shape=(len(target_ids), nfeatures))
    x[:] = 0.0  # BIOM sparse implicit values are explicitly stored as zero.
    if ijson is None:
        with gzip.open(path, 'rt', encoding='utf-8') as stream:
            values = json.load(stream)['data']
    else:
        with gzip.open(path, 'rb') as stream:
            values = ijson.items(stream, 'data.item')
            _fill_biom(x, values, col_out, nfeatures, matrix_type)
            x.flush()
            return x
    _fill_biom(x, values, col_out, nfeatures, matrix_type)
    x.flush()
    return x


def _fill_biom(x, rows, col_out, nfeatures, matrix_type):
    if matrix_type == 'sparse':
        for triple in rows:
            row, col, val = triple
            row, col = int(row), int(col)
            if not 0 <= row < nfeatures:
                raise ValueError(f'BIOM sparse row index out of bounds: {row}')
            if col in col_out:
                x[col_out[col], row] = float(val)
    else:
        row_count=0
        for i, vals in enumerate(rows):
            if i >= nfeatures:
                raise ValueError('BIOM dense array has too many rows')
            for j, k in col_out.items():
                x[k, i] = float(vals[j])
            row_count+=1
        if row_count != nfeatures:
            raise ValueError(f'BIOM dense data rows {row_count} != features {nfeatures}')


def export_proteomics(path: Path, sample_ids: list[str], header: list[str], skiprows: int, out: Path) -> int:
    """Read a feature-by-sample TSV by chunks and write sample-by-feature CSV."""
    columns = pd.read_csv(path, compression='gzip', sep='\t', nrows=0, skiprows=skiprows).columns.to_list()
    if len(columns) != 1+len(header):
        raise ValueError('Protein column count differs from inspected header')
    # Header names must be unchanged; else matching would be ambiguous.
    if [x.strip() for x in columns[1:]] != header:
        raise ValueError('Protein sample column names changed during parsing')
    protein = pd.read_csv(path, sep='\t', compression='gzip', dtype=str, low_memory=False, skiprows=skiprows)
    features = protein.iloc[:,0].map(_text).tolist()
    if len(set(features)) != len(features) or not all(features):
        raise ValueError('Protein feature names contain blanks or duplicates')
    subset = protein.loc[:, sample_ids].apply(pd.to_numeric, errors='coerce')
    if np.isinf(subset.to_numpy(dtype=float)).any():
        raise ValueError('Protein matrix contains infinities')
    with out.open('w', encoding='utf-8', newline='') as fh:
        writer=csv.writer(fh)
        writer.writerow(['sample_id', *features])
        for sample in sample_ids:
            writer.writerow([sample, *subset[sample].where(subset[sample].notna(), '').tolist()])
    return len(features)


def export_metabolomics(matrix: np.memmap, sample_ids: list[str], features: list[str], out: Path):
    with out.open('w', encoding='utf-8', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['sample_id', *features])
        for i, sample in enumerate(sample_ids):
            values = matrix[i]
            writer.writerow([sample, *[float(v) if np.isfinite(v) else '' for v in values]])


def prepare(root: Path, inspect_only=False) -> dict:
    raw = root / 'raw' / 'hmp2'
    output = root / 'processed' / 'external_validation' / 'hmp2'
    for p in (METADATA, PROTEIN, METAB):
        if not (raw/p).is_file():
            raise FileNotFoundError(f'Missing HMP2 primary source: {raw/p}')
    metadata, meta_report = metadata_map(raw/METADATA)
    _, p_ids, p_header_lines = protein_ids(raw/PROTEIN)
    m_ids, features, matrix_type, shape = biom_metadata(raw/METAB)
    subjects, join_report = matched_visits(metadata, p_ids, m_ids)
    report = {'source': {key: str(raw / name) for key,name in (
        ('metadata',METADATA),('protein',PROTEIN),('metabolomics',METAB))},
        'metadata':meta_report, 'biom':{'matrix_type':matrix_type, 'shape':shape},
        'join':join_report, 'status':'inspected'}
    if inspect_only:
        return report
    output.mkdir(parents=True,exist_ok=True)
    report_file = output/'preparation_report.json'
    if len(subjects) < 8:
        report['status'] = 'insufficient_exact_matched_visits'
        report_file.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
        return report
    names = subjects.sample_id.tolist()
    fd, mempath = tempfile.mkstemp(prefix='hmp2_',suffix='.memmap',dir=output)
    os.close(fd)
    temporary = {name: output/f'{name}.building' for name in ('subjects.csv','proteomics.csv','metabolomics.csv')}
    try:
        subjects.to_csv(temporary['subjects.csv'],index=False,columns=['sample_id','patient_id','visit_order'])
        n_protein=export_proteomics(raw/PROTEIN,names,p_ids,p_header_lines,temporary['proteomics.csv'])
        x=stream_biom_matrix(raw/METAB,names,m_ids,len(features),matrix_type,Path(mempath))
        export_metabolomics(x,names,features,temporary['metabolomics.csv'])
        del x
        for name,p in temporary.items():
            p.replace(output/name)
        report['status']='canonical_raw_omics_exported'
        report['features']={'protein':n_protein,'metabolomics':len(features)}
        report['biom']['sparse_missing_cell_convention']='implicit BIOM entries exported as zero'
    finally:
        Path(mempath).unlink(missing_ok=True)
        for p in temporary.values():
            p.unlink(missing_ok=True)
    report_file.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=Path(os.getenv('PULSE_OMICS_DATA_ROOT',ROOT)))
    p.add_argument('--inspect-only',action='store_true',help='Only report exact available assay matches; do not export')
    args=p.parse_args(argv)
    result=prepare(args.data_root,inspect_only=args.inspect_only)
    print('\nHMP2 两种组学实际配对')
    for key,val in result['join'].items():
        if key!='join_basis':print(f'  {key}: {val}')
    print('状态：', result['status'])
    return 0 if result['status'] in ('inspected','canonical_raw_omics_exported') else 2


if __name__=='__main__':
    raise SystemExit(main())