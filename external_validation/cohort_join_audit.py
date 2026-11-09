"""Source-specific, identity-preserving joint-omics readiness audit.

Reads already-downloaded COVAAR and CSIRO releases. Writes only aggregate
statistics to F:. An exact published crosswalk is required before exporting
joint assay measurements. No numerical sample IDs are inferred by row order.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
from collections import Counter
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

ROOT = Path('/mnt/f/datasets/pulse-omics')


def _canonical_text(value) -> str | None:
    if pd.isna(value):
        return None
    value = str(value).strip()
    if not value or value.casefold() in {'na', 'none', 'nan'}:
        return None
    return value


def _column_overlap(left, right) -> dict:
    a = [_canonical_text(v) for v in left]
    b = [_canonical_text(v) for v in right]
    sa, sb = set(a)-{None}, set(b)-{None}
    return {'left_distinct': len(sa), 'right_distinct': len(sb),
            'matching_distinct': len(sa & sb),
            'left_duplicate_values': sum(n > 1 for key,n in Counter(a).items() if key is not None),
            'right_duplicate_values': sum(n > 1 for key,n in Counter(b).items() if key is not None)}


def audit_covaar(root: Path) -> dict:
    folder = root / 'raw/covaar'
    meta = pd.read_csv(folder / 'figshare/Sample_Meta_Data.csv', dtype=str)
    for col in ('Sample_ID', 'Patient', 'Time'):
        if col not in meta.columns:
            raise ValueError(f'COVAAR metadata missing {col}')
    with (folder/'figshare/Metabolomics_Data.csv').open(encoding='utf-8-sig') as f:
        metabolism_cols = pd.read_csv(f, nrows=0).columns.tolist()
    with gzip.open(folder/'geo/GSE325178_SomaScan_RFUs_allSamples.txt.gz', 'rt', encoding='utf-8-sig') as f:
        protein_cols = pd.read_csv(f, sep='\t', nrows=0).columns.tolist()
    meta_ids = meta['Sample_ID'].map(_canonical_text)
    met_ids = [s for s in metabolism_cols if s != 'Metabolite']
    protein_ids = protein_cols[5:]  # first five are assay-feature annotation fields
    # GEO sample metadata are column-wise '!Sample_*' entries, not the assay data.
    series = folder / 'geo/GSE325178_series_matrix.txt.gz'
    available = {'!Sample_title': [], '!Sample_geo_accession': []}
    characteristics = []
    geo_attributes = {}
    protein_id_set = set(protein_ids)
    # A GEO title/accession need not be the assay specimen ID. Search every
    # published per-sample metadata field for exact VA identifiers, including
    # the value following a named "sample id:" characteristic.
    geo_id_locations = {}
    with gzip.open(series, 'rt', encoding='utf-8-sig', errors='replace') as f:
        for line in f:
            if not line.startswith('!Sample_'):
                if line.startswith('!series_matrix_table_begin'):
                    break
                continue
            cells = next(__import__('csv').reader([line.rstrip('\n')], delimiter='\t'))
            if not cells: continue
            name = cells[0]
            values = cells[1:]
            for index, value in enumerate(values):
                direct = value.strip().strip('"')
                after_label = direct.rsplit(':', 1)[-1].strip()
                for candidate in {direct, after_label} & protein_id_set:
                    geo_id_locations.setdefault(candidate, set()).add(index)
            if name in available: available[name].extend(values)
            if name.startswith('!Sample_characteristics_ch1'):
                for index,value in enumerate(values):
                    if ':' not in value: continue
                    field,raw_value=value.split(':',1)
                    field=field.strip().casefold().replace(' ','').replace('_','')
                    if field in {'patient','subject','participant','subjectid','participantid','time','timepoint','visittime'}:
                        geo_attributes.setdefault(index,{})[field]=raw_value.strip()
            if name.startswith('!Sample_characteristics_ch1') or name in ('!Sample_description', '!Sample_source_name_ch1'):
                characteristics.append({'field': name,
                    'n_values':len(values),
                    'examples':[v[:100] for v in values[:2]]})
    patient_time = meta[['Patient','Time']].copy()
    patient_time = patient_time.dropna()
    time_counts = patient_time.groupby(['Patient','Time']).size()
    mapping_status = ('cross_assay_ids_match' if set(protein_ids) & set(met_ids)
                      else 'cross_assay_ids_differ_need_verified_crosswalk')
    # A joint patient/time lookup is informative only when those exact
    # identity fields occur in official GEO per-sample characteristics.
    unique_geo_locations = {
        sample: next(iter(indices)) for sample,indices in geo_id_locations.items()
        if len(indices) == 1
    }
    ambiguous_geo_locations = sum(len(indices) != 1 for indices in geo_id_locations.values())
    geo_patient_time = set()
    for sample, index in unique_geo_locations.items():
        attr = geo_attributes.get(index, {})
        patient = next((attr.get(k) for k in ('patient','subject','participant','subjectid','participantid') if attr.get(k)),None)
        visit = next((attr.get(k) for k in ('time','timepoint','visittime') if attr.get(k)),None)
        if patient and visit:geo_patient_time.add((patient.casefold(),visit.casefold()))
    figshare_patient_time = set((str(x).strip().casefold(),str(y).strip().casefold())
                               for x,y in zip(meta['Patient'],meta['Time'])
                               if pd.notna(x) and pd.notna(y))
    tuple_matches = len(geo_patient_time & figshare_patient_time)
    return {'dataset':'covaar',
            'metabolomics_metadata_rows':len(meta), 'metabolomics_columns':len(met_ids),
            'proteomics_columns':len(protein_ids),
            'metabolite_to_metadata':_column_overlap(met_ids,meta_ids),
            'protein_to_metabolite':_column_overlap(protein_ids,met_ids),
            'protein_to_geo_title':_column_overlap(protein_ids,available['!Sample_title']),
            'protein_to_geo_accession':_column_overlap(protein_ids,available['!Sample_geo_accession']),
            'protein_ids_found_in_any_geo_sample_field':len(unique_geo_locations),
            'ambiguous_geo_sample_id_locations':ambiguous_geo_locations,
            'metadata_patients':int(meta.Patient.nunique()),
            'patient_time_unique':bool((time_counts <= 1).all()),
            'patient_time_duplicate_groups':int((time_counts > 1).sum()),
            'geo_characteristic_fields':characteristics[:12],
            'geo_patient_time_pairs':len(geo_patient_time),
            'patient_time_tuples_matching_figshare':tuple_matches,
            'status':mapping_status}


def _find_excel(folder: Path, prefix: str) -> Path:
    matching = sorted(p for p in folder.glob('*.xlsx') if p.name.startswith(prefix) and '%' not in p.name)
    if len(matching) != 1:
        raise ValueError(f'CSIRO expected one decoded filename with prefix {prefix}, got {len(matching)}')
    return matching[0]


def _vac_id(value) -> str | None:
    value = _canonical_text(value)
    if value is None: return None
    m = re.fullmatch(r'(?:Abundance:\s*)?(?:Vac[_ ]?Plasma[_ ]?)?(\d{1,5})(?:\.0)?', value, flags=re.I)
    return f'Vac_Plasma_{int(m.group(1)):03d}' if m else None


def _excel_sheet(path: Path, sheet: str) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet, dtype=str)


def _protein_header(path: Path, sheet: str) -> tuple[list[str], int]:
    wb = load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb[sheet]
        return ([str(x) for x in next(ws.iter_rows(min_row=1, max_row=1, values_only=True)) if x is not None], ws.max_row)
    finally:
        wb.close()


def audit_csiro(root: Path) -> dict:
    folder = root/'raw/csiro'
    mapping = _excel_sheet(_find_excel(folder, '005_Task'), 'Sheet1')
    metab = _excel_sheet(_find_excel(folder, '006_Task'), 'COVID HPV CCM')
    required_map = ['Sample IDs', 'Treatments', 'proteomic#', 'Vac']
    required_met = ['Subject ID', 'Sample IDS (Barcode)', 'Times', 'IDs & Treatments']
    for c in required_map:
        if c not in mapping: raise ValueError(f'CSIRO mapping missing column {c}')
    for c in required_met:
        if c not in metab: raise ValueError(f'CSIRO metabolomics missing column {c}')
    protein = {}
    protein_columns = []
    for code in ('AZ','Moderna','Pfizer'):
        path = _find_excel(folder, '002_Task' if code=='AZ' else '003_Task' if code=='Moderna' else '004_Task')
        names, row_count = _protein_header(path,'Proteins')
        p = [v for v in names if v.startswith('Abundance:')]
        protein[code]={'features_rows':int(row_count - 1),
                       'abundance_columns':len(p)}
        protein_columns.extend(p)
    vac_ids = [_vac_id(x) for x in protein_columns]
    map_vac = [_vac_id(x) for x in mapping['proteomic#']]
    joins={field:_column_overlap(mapping['Sample IDs'],metab[field])
           for field in ('Sample IDS (Barcode)','IDs & Treatments')}
    # Numeric IDs are only normalized for the explicitly documented 'proteomic#'
    # field and the actual Vac_Plasma_N column suffix, never for patient identities.
    valid_ids = _column_overlap(map_vac,vac_ids)
    metadata_keys=_column_overlap(mapping['Sample IDs'],metab['Subject ID'])
    return {'dataset':'csiro', 'crosswalk_rows':len(mapping),
            'metabolomics_rows':len(metab), 'metabolomics_features':int(len(metab.columns)-8),
            'protein_sources':protein,
            'mapped_proteomics_to_actual_columns':valid_ids,
            'crosswalk_sample_id_to_metabolomics':joins,
            'mapping_id_to_metabolomics_subject':metadata_keys,
            'metabolomics_subjects':int(metab['Subject ID'].nunique()),
            'metabolomics_times':sorted(str(v) for v in metab['Times'].dropna().unique())[:15],
            'status':('candidate_direct_crosswalk' if valid_ids['matching_distinct']>0 and
                     max(x['matching_distinct'] for x in joins.values())>0
                     else 'crosswalk_not_established')}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root', type=Path, default=Path(os.getenv('PULSE_OMICS_DATA_ROOT', ROOT)))
    p.add_argument('--datasets', nargs='+', choices=['covaar','csiro'], default=['covaar','csiro'])
    a=p.parse_args(argv)
    for name in a.datasets:
        try:
            report = audit_covaar(a.data_root) if name=='covaar' else audit_csiro(a.data_root)
            out = a.data_root/'processed/external_validation'/name/'join_audit.json'
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
            print(f'{name.upper()}: {report["status"]}')
            if name=='covaar':
                for k in ('metabolite_to_metadata','protein_to_metabolite','protein_to_geo_title','protein_to_geo_accession'):
                    print(f'  {k}: {report[k]["matching_distinct"]} exact IDs')
                print(f'  patients={report["metadata_patients"]}, unique patient/time={report["patient_time_unique"]}')
                print('  Exact protein IDs appearing in any GEO sample metadata field:',
                      report['protein_ids_found_in_any_geo_sample_field'])
                print('  GEO characteristics:', ', '.join(v['field'] for v in report['geo_characteristic_fields']))
                print('  GEO patient/time pairs matching Figshare:',report['patient_time_tuples_matching_figshare'])
            else:
                print('  mapped proteomics IDs:',report['mapped_proteomics_to_actual_columns']['matching_distinct'])
                for k,v in report['crosswalk_sample_id_to_metabolomics'].items():
                    print(f'  Sample IDs to metabolomics {k}: {v["matching_distinct"]} exact IDs')
                print(f'  patients={report["metabolomics_subjects"]}, times={report["metabolomics_times"]}')
            print('  saved:', out)
        except (OSError, ValueError, KeyError, ImportError) as exc:
            print(f'{name.upper()}: {type(exc).__name__}: {exc}')
            return 2
    return 0


if __name__=='__main__':
    raise SystemExit(main())