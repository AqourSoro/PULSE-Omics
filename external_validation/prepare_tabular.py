"""Explicit, provenance-preserving adapter for published longitudinal assay tables.

A cohort-specific JSON mapping supplies every identifier column. No guessing IDs,
visit order, orientation, or equivalence between specimens. Data stays on F:.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path('/mnt/f/datasets/pulse-omics')
KINDS = ('protein', 'metab')


def _file(root: Path, record: dict) -> Path:
    rel = Path(record['path'])
    if rel.is_absolute() or '..' in rel.parts:
        raise ValueError('Mapping paths must be relative to the data root')
    resolved = (root / rel).resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as exc:
        raise ValueError('Source path escapes the data root') from exc
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    return resolved


def _read(root: Path, record: dict) -> pd.DataFrame:
    path = _file(root, record)
    lower = path.name.lower()
    if lower.endswith(('.csv', '.csv.gz')):
        frame = pd.read_csv(path, dtype=str)
    elif lower.endswith(('.tsv', '.tsv.gz', '.txt', '.txt.gz')):
        frame = pd.read_csv(path, sep='\t', dtype=str)
    elif lower.endswith(('.xlsx', '.xlsm')):
        if not record.get('sheet'):
            raise ValueError(f'{path.name}: Excel sheet must be specified explicitly')
        frame = pd.read_excel(path, sheet_name=record['sheet'], dtype=str)
    else:
        raise ValueError(f'Unsupported tabular file: {path.name}')
    if frame.empty or frame.columns.duplicated().any():
        raise ValueError(f'{path.name}: empty or duplicate headers')
    return frame


def _ids(frame: pd.DataFrame, column: str, name: str) -> pd.Series:
    if column not in frame:
        raise ValueError(f'{name}: missing identifier column {column!r}')
    ids = frame[column].astype('string').str.strip()
    if ids.isna().any() or ids.eq('').any() or ids.duplicated().any():
        raise ValueError(f'{name}: identifiers must be present and unique')
    return ids.astype(str)


def _assay(root: Path, config: dict, name: str) -> pd.DataFrame:
    frame = _read(root, config)
    orientation = config.get('orientation', 'samples_as_rows')
    identifier = config['sample_id_column']
    if orientation == 'samples_as_rows':
        ids = _ids(frame, identifier, name)
        excluded = set(config.get('exclude_columns', [])) | {identifier}
        if not excluded.issubset(frame.columns):
            raise ValueError(f'{name}: unknown excluded columns')
        features = config.get('feature_columns')
        if features is None:
            raise ValueError(f'{name}: feature_columns must be explicit')
        if len(features) != len(set(features)) or not features or any(c not in frame or c in excluded for c in features):
            raise ValueError(f'{name}: invalid or duplicate assay features')
        numeric = frame.loc[:, features].apply(pd.to_numeric, errors='coerce')
        # Nonnumeric codes must not silently become missing measurements.
        for feature in features:
            originally = frame[feature].astype('string').str.strip()
            allowed_missing = originally.isna() | originally.isin(config.get('missing_tokens', ['', 'NA', 'N/A', 'NaN']))
            if (numeric[feature].isna() & ~allowed_missing).any():
                raise ValueError(f'{name}: nonnumeric measurement in {feature!r}')
        numeric.index = ids
    elif orientation == 'samples_as_columns':
        if identifier not in frame:
            raise ValueError(f'{name}: feature identifier column {identifier!r} absent')
        names = _ids(frame, identifier, f'{name} feature names')
        columns = config.get('sample_columns')
        if not columns or len(columns) != len(set(columns)) or any(c not in frame or c == identifier for c in columns):
            raise ValueError(f'{name}: provide explicit sample_columns')
        numeric = frame[columns].apply(pd.to_numeric, errors='coerce')
        for col in columns:
            orig = frame[col].astype('string').str.strip()
            miss = orig.isna() | orig.isin(config.get('missing_tokens', ['', 'NA', 'N/A', 'NaN']))
            if (numeric[col].isna() & ~miss).any():
                raise ValueError(f'{name}: nonnumeric assay value in {col!r}')
        numeric.index = names
        numeric = numeric.T
        numeric.index.name = 'source_sample_id'
    else:
        raise ValueError('orientation must be samples_as_rows or samples_as_columns')
    if not np.isfinite(numeric.to_numpy(dtype=float)[~numeric.isna().to_numpy()]).all():
        raise ValueError(f'{name}: infinite assay measurement')
    if numeric.shape[1] == 0:
        raise ValueError(f'{name}: no assay features')
    return numeric.astype('float64')


def _visit_order(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        raise ValueError(f'Missing visit column {column!r}')
    values = pd.to_numeric(frame[column], errors='coerce')
    if values.isna().any() or not np.isfinite(values).all():
        raise ValueError('visit_order must be an explicit finite numeric measurement')
    return values


def prepare(root: Path, config: dict, *, inspect_only: bool = False) -> dict:
    dataset = config['dataset']
    if dataset not in {'ipop', 'covaar', 'csiro'}:
        raise ValueError('This adapter handles only iPOP, COVAAR and CSIRO')
    protein = _assay(root, config['protein'], 'protein')
    metab = _assay(root, config['metab'], 'metab')
    bridge_cfg = config['sample_crosswalk']
    link = _read(root, bridge_cfg)
    columns = [bridge_cfg[k] for k in ('protein_sample_id', 'metab_sample_id', 'patient_id', 'visit_order')]
    if len(columns) != len(set(columns)) or any(x not in link for x in columns):
        raise ValueError('Crosswalk must contain four distinct, named columns')
    bridge = link.loc[:, columns].copy()
    bridge.columns = ['protein_id', 'metab_id', 'patient_id', 'visit_order']
    for col in ('protein_id', 'metab_id', 'patient_id'):
        bridge[col] = bridge[col].astype('string').str.strip()
        if bridge[col].isna().any() or bridge[col].eq('').any():
            raise ValueError(f'Crosswalk: missing {col}')
    bridge['visit_order'] = _visit_order(bridge, 'visit_order')
    # Repeated source specimen in two different patient/visit records is never safe.
    for col in ('protein_id', 'metab_id'):
        if bridge[col].duplicated().any():
            raise ValueError(f'Crosswalk: duplicate assay specimen {col}')
    if bridge.duplicated(['patient_id', 'visit_order']).any():
        raise ValueError('Crosswalk: duplicate patient/visit')
    pair = bridge.loc[bridge.protein_id.isin(protein.index) & bridge.metab_id.isin(metab.index)].copy()
    pair = pair.sort_values(['patient_id', 'visit_order']).reset_index(drop=True)
    pair.insert(0, 'sample_id', [f'{dataset}_{i:06d}' for i in range(len(pair))])
    counts = pair.patient_id.value_counts()
    report = {
        'dataset': dataset,
        'protein_source_samples': int(len(protein)),
        'metab_source_samples': int(len(metab)),
        'crosswalk_entries': int(len(bridge)),
        'verified_paired_samples': int(len(pair)),
        'paired_subjects': int(pair.patient_id.nunique()),
        'subjects_with_two_or_more_pairs': int((counts >= 2).sum()),
        'eligible_historical': int((counts - 1).clip(lower=0).sum()),
        'join': 'source-specific explicit crosswalk with distinct specimen IDs and patient/time',
        'status': 'inspected',
    }
    if inspect_only:
        return report
    if pair.patient_id.nunique() < 8:
        return {**report, 'status': 'insufficient_paired_subjects'}
    output = root / 'processed' / 'external_validation' / dataset
    output.mkdir(parents=True, exist_ok=True)
    subj = pair.loc[:, ['sample_id', 'patient_id', 'visit_order']]
    prot = protein.loc[pair.protein_id].reset_index(drop=True)
    m = metab.loc[pair.metab_id].reset_index(drop=True)
    prot.insert(0, 'sample_id', pair.sample_id)
    m.insert(0, 'sample_id', pair.sample_id)
    pending = [output / (n + '.building') for n in ('subjects.csv', 'proteomics.csv', 'metabolomics.csv')]
    try:
        for frame, path in zip((subj, prot, m), pending):
            frame.to_csv(path, index=False)
        for path, name in zip(pending, ('subjects.csv', 'proteomics.csv', 'metabolomics.csv')):
            path.replace(output/name)
        report['status'] = 'canonical_raw_omics_exported'
        (output/'preparation_report.json').write_text(json.dumps(report, indent=2)+'\n',encoding='utf-8')
    finally:
        for path in pending:
            path.unlink(missing_ok=True)
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,default=Path(os.environ.get('PULSE_OMICS_DATA_ROOT', ROOT)))
    parser.add_argument('--mapping', type=Path, required=True, help='Explicit observed source schema and sample crosswalk, JSON')
    parser.add_argument('--inspect-only', action='store_true')
    args = parser.parse_args(argv)
    report = prepare(args.data_root.expanduser(), json.loads(args.mapping.read_text(encoding='utf-8')),
                     inspect_only=args.inspect_only)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report['status'] in {'inspected','canonical_raw_omics_exported'} else 2


if __name__ == '__main__':
    raise SystemExit(main())