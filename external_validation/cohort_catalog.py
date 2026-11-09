"""Inventory publisher tables and every measured feature, keeping participant rows local."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook

from external_validation import run as validation
from external_validation.prepare_hmp2 import METAB, PROTEIN, METADATA

UNITS = {'su2020':'T1/T2 order', 'csiro':'study stage order',
         'ipop':'elapsed day', 'hmp2':'study week number', 'covaar':'trial-event order'}
XLSX_LABELS = {'S1.3 Plasma Proteomic Data':2,'S1.4 Plasma Metabolomic Data':2,'COVID HPV CCM':8,'HPC QAQC':1}


def workbook(path, feature_file):
    wb = load_workbook(path,read_only=True,data_only=True)
    sheets = []
    try:
        for ws in wb:
            rows = list(ws.values)
            nonempty = [i for i,row in enumerate(rows,1) if any(x is not None for x in row)]
            width = ws.max_column
            headers = [3] if path.name.startswith('000_') else [1]
            if ws.title.startswith('S1.5 '):
                headers = [1,2]
            if ws.title.startswith('S1.6 '):
                headers = [1,2,3]
            shape = {'sheet':ws.title,'worksheet_rows':len(rows),'columns':width,
                     'nonempty_rows':len(nonempty), 'last_nonempty_row':max(nonempty,default=0)}
            if 'QToF' in ws.title and 'data' in ws.title:
                sample_header = next(i for i,row in enumerate(rows) if row[0] == 'Data file')
                labels = [(i+1,row[8]) for i,row in enumerate(rows[:sample_header]) if isinstance(row[8],str)]
                shape['sample_header_row'] = sample_header+1
                shape['sample_metadata_fields'] = list(rows[sample_header][:9])
                shape['sample_rows_by_label'] = dict(Counter(r[7] for r in rows[sample_header+1:]))
                shape['feature_axis'] = 'columns from J, annotated by feature-property rows'
                shape['native_feature_property_rows'] = labels
                compound_cols = [j for j in range(9,width) if isinstance(rows[0][j],(int,float))]
                shape['feature_count'] = len(compound_cols)
                for j in compound_cols:
                    feature_file.write(json.dumps({'file':path.name,'sheet':ws.title,'column':j+1,
                        'raw_feature_id':rows[0][j],'properties':{f'{i}:{label}':rows[i-1][j] for i,label in labels}},
                        ensure_ascii=False,default=str)+'\n')
            else:
                shape['header_rows'] = {str(h):rows[h-1] for h in headers}
                start = headers[-1]
                shape['fields'] = []
                for col in range(width):
                    types = Counter(type(row[col]).__name__ for row in rows[start:])
                    empty = sum(row[col] is None or (isinstance(row[col],str) and not row[col].strip()) for row in rows[start:])
                    shape['fields'].append({'column':col+1,'name':rows[start-1][col],
                                            'empty_cells':empty,'native_types':dict(types)})
                if ws.title in XLSX_LABELS:
                    start_col = XLSX_LABELS[ws.title]
                    shape['feature_axis'] = 'columns'
                    shape['feature_count'] = width-start_col
                    for j in range(start_col,width):
                        feature_file.write(json.dumps({'file':path.name,'sheet':ws.title,'column':j+1,
                            'raw_feature_id':rows[start-1][j]},ensure_ascii=False,default=str)+'\n')
                if ws.title.startswith('S1.6 '):
                    for j in (0,11,22):
                        for i,row in enumerate(rows[3:],4):
                            if row[j] is not None:
                                feature_file.write(json.dumps({'file':path.name,'sheet':ws.title,
                                    'row':i,'column':j+1,'raw_feature_id':row[j],
                                    'comparison':rows[1][j],'annotation':dict(zip(rows[2][j+8:j+11],row[j+8:j+11]))},
                                    ensure_ascii=False,default=str)+'\n')
                is_feature_table = (ws.title == 'Proteins' or 'RAW DATA' in ws.title.upper()
                                    or ws.title.endswith('library name') or 'lipid id stats' in ws.title.lower()
                                    or 'Surface Protein Panel' in ws.title)
                if is_feature_table:
                    header = rows[0]
                    key = header.index('Accession') if ws.title=='Proteins' else 0
                    fields = [j for j,c in enumerate(header) if c is not None and not str(c).startswith('Abundance:')]
                    if 'RAW DATA' in ws.title.upper():
                        fields = [0]
                    shape['feature_axis'] = 'rows'
                    shape['feature_count'] = sum(row[key] is not None for row in rows[1:])
                    for i,row in enumerate(rows[1:],2):
                        if row[key] is not None:
                            feature_file.write(json.dumps({'file':path.name,'sheet':ws.title,'row':i,
                                'raw_feature_id':row[key],'annotation':{str(header[j]):row[j] for j in fields}},
                                ensure_ascii=False,default=str)+'\n')
            sheets.append(shape)
    finally:
        wb.close()
    return sheets


def hmp2_native(folder, out):
    import ijson
    p = folder/METAB
    with gzip.open(p,'rb') as f:
        prefix = f.read(512).decode()
    biom_header = json.loads(prefix.split(',"data":')[0]+'}')
    features = 0
    annotation_fields = Counter()
    with gzip.open(p,'rb') as source, gzip.open(out/'hmp2_biom_features.jsonl.gz','wt',encoding='utf-8') as target:
        for rec in ijson.items(source,'rows.item'):
            target.write(json.dumps(rec,ensure_ascii=False,default=str)+'\n')
            annotation_fields.update((rec.get('metadata') or {}).keys())
            features += 1
    sample_fields = Counter()
    with gzip.open(p,'rb') as source:
        for rec in ijson.items(source,'columns.item'):
            sample_fields.update((rec.get('metadata') or {}).keys())
    assert features == biom_header['shape'][0]
    with gzip.open(folder/PROTEIN,'rt',encoding='utf-8-sig',newline='') as src, gzip.open(out/'hmp2_protein_native_features.csv.gz','wt',encoding='utf-8',newline='') as dest:
        reader = csv.reader(src,delimiter='\t')
        header = next(reader)
        writer = csv.writer(dest)
        writer.writerow(['source_row','protein_identifier','taxon_label_prefix'])
        for n,row in enumerate(reader,2):
            if len(row) != len(header):
                raise ValueError(f'HMP2 ragged protein TSV row {n}')
            writer.writerow([n,row[0],row[0].split(':')[0] if ':' in row[0] else ''])
    metadata = pd.read_csv(folder/METADATA,low_memory=False)
    return {'biom':biom_header,'biom_feature_annotation_fields':dict(annotation_fields),
            'biom_sample_annotation_fields':dict(sample_fields),
            'protein':{'rows_features':n-1,'columns_samples':len(header)-1,'first_header':header[0]},
            'metadata':{'rows':len(metadata),'fields':[{ 'name':c,'dtype':str(metadata[c].dtype),
                        'missing':int(metadata[c].isna().sum())} for c in metadata]}}


def paired_catalog(root, dataset, out):
    folder = validation.dataset_dir(root,dataset)
    subjects = pd.read_csv(folder/'subjects.csv',dtype={'sample_id':str,'patient_id':str})
    raw = {}
    for mod, file in zip(validation.MODALITIES,('proteomics.csv','metabolomics.csv')):
        frame = pd.read_csv(folder/file,dtype={'sample_id':str}).set_index('sample_id')
        raw[mod] = frame.reindex(subjects.sample_id)
    p, m = [v.notna().any(axis=1).to_numpy() for v in raw.values()]
    both = p&m
    paired = subjects.loc[both].sort_values(['patient_id','visit_order'])
    counts = paired.groupby('patient_id').size()
    gaps = paired.groupby('patient_id').visit_order.diff().dropna().to_numpy()
    shapes = {'patients_in_linked_table':int(subjects.patient_id.nunique()),'visits_in_linked_table':len(subjects),
              'paired_visits':int(both.sum()),'paired_patients':len(counts),'paired_patients_ge2':int((counts>=2).sum()),
              'paired_patients_ge3':int((counts>=3).sum()),'paired_patients_ge4':int((counts>=4).sum()),
              'paired_visits_per_patient':{str(k):int(v) for k,v in counts.value_counts().sort_index().items()},
              'protein_only_visits':int((p&~m).sum()),'metabolite_only_visits':int((m&~p).sum()),
              'time_basis':UNITS[dataset],'gap_n':len(gaps),
              'gap_min_q25_median_q75_max':np.quantile(gaps,[0,.25,.5,.75,1]).tolist(), 'assays':{}}
    for mod,frame in raw.items():
        matrix = frame.to_numpy(dtype=float)
        observed = np.isfinite(matrix)
        available = observed.any(1)
        matrix = matrix[available]
        observed = observed[available]
        measured = observed.sum(0)
        stats = pd.DataFrame({'feature':frame.columns,'position_1based':np.arange(1,len(frame.columns)+1),
                              'finite':measured,'missing_in_assayed_rows':len(matrix)-measured,
                              'zero':(matrix==0).sum(0),'negative':(matrix<0).sum(0),
                              'minimum':np.min(np.where(observed,matrix,np.inf),axis=0),
                              'maximum':np.max(np.where(observed,matrix,-np.inf),axis=0)})
        stats.loc[measured==0,['minimum','maximum']] = np.nan
        stats.to_csv(out/f'{dataset}_{mod}_features.csv.gz',index=False,compression='gzip')
        shapes['assays'][mod] = {'assayed_rows':len(matrix),'features':len(stats),'finite_cells':int(measured.sum()),
            'missing_cells_in_assayed_rows':int(matrix.size-measured.sum()),'zero_cells':int(stats.zero.sum()),
            'negative_cells':int(stats.negative.sum()),'min':float(stats.minimum.min()),'max':float(stats.maximum.max()),
            'all_zero_observed_features':int(((stats.finite>0)&(stats.zero==stats.finite)).sum())}
    return shapes


def covaar_native(folder, out):
    result = {'figshare':{},'geo':{}}
    for path in sorted((folder/'figshare').glob('*.csv')):
        df = pd.read_csv(path)
        result['figshare'][path.name] = {'rows':len(df),'fields':[{
            'name':c,'dtype':str(df[c].dtype),'missing':int(df[c].isna().sum())} for c in df]}
        if path.name == 'Metabolite_Meta_Data.csv':
            df.to_csv(out/'covaar_metabolite_annotations.csv',index=False)
    for path in sorted((folder/'geo').glob('*.txt.gz')):
        if 'RFUs' in path.name:
            df = pd.read_csv(path,sep='\t')
            samples = [c for c in df if c.startswith('VA') and c[2:].isdigit()]
            annotations = df.drop(columns=samples).copy()
            a = df[samples].to_numpy(dtype=float)
            annotations['finite'] = np.isfinite(a).sum(1)
            annotations['zero'] = (a==0).sum(1)
            annotations['minimum'] = np.nanmin(a,axis=1)
            annotations['maximum'] = np.nanmax(a,axis=1)
            annotations.to_csv(out/'covaar_somascan_features.csv.gz',index=False,compression='gzip')
            result['geo'][path.name] = {'feature_rows':len(df),'sample_columns':a.shape[1],
                                       'annotation_fields':[c for c in df if c not in samples],
                                       'feature_id_header':str(df.columns[0]),'measurement_scale':'published SomaScan RFU'}
        else:
            headers, matrix_rows, in_table = [], 0, False
            with gzip.open(path,'rt',encoding='utf-8-sig',newline='') as f:
                for row in csv.reader(f,delimiter='\t'):
                    if not row:
                        continue
                    key = row[0]
                    if key == '!series_matrix_table_begin':
                        in_table = True
                    elif key == '!series_matrix_table_end':
                        in_table = False
                    elif in_table:
                        matrix_rows += 1
                    elif key.startswith('!Sample_'):
                        headers.append({'name':key,'sample_columns':len(row)-1,
                                        'characteristic':row[1].split(':',1)[0] if key.startswith('!Sample_characteristics') and len(row)>1 else None})
            result['geo'][path.name] = {'metadata_rows':headers,'matrix_rows_including_header':matrix_rows}
    return result


def run_catalog(root, out):
    out.mkdir(parents=True,exist_ok=True)
    report = {'canonical':{},'source_workbooks':{},'files':[]}
    for d in UNITS:
        if d == 'covaar' and not (validation.dataset_dir(root,d)/'subjects.csv').is_file():
            continue
        report['canonical'][d] = paired_catalog(root,d,out)
    for cohort in ('su2020','csiro'):
        with gzip.open(out/f'{cohort}_native_features.jsonl.gz','wt',encoding='utf-8') as feature_file:
            sources = sorted((root/'raw'/cohort).glob('*.xlsx'))
            report['source_workbooks'][cohort] = {p.name:workbook(p,feature_file) for p in sources}
    report['hmp2_native'] = hmp2_native(root/'raw/hmp2',out)
    if (root/'raw/covaar').is_dir():
        report['covaar_native'] = covaar_native(root/'raw/covaar',out)
    for cohort in UNITS:
        for file in sorted((root/'raw'/cohort).rglob('*')):
            if not file.is_file() or not file.name.lower().endswith(('.csv','.gz','.xlsx','.rdata','.txt')):
                continue
            with file.open('rb') as f:
                sha = hashlib.file_digest(f,'sha256').hexdigest()
            report['files'].append({'cohort':cohort,'source_file':str(file.relative_to(root/'raw')),
                                    'bytes':file.stat().st_size,'sha256':sha})
    validation._write_json(out/'structure.json',report)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=validation.DEFAULT_ROOT)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    out=a.output or a.data_root/'processed/real_study_catalog'
    run_catalog(a.data_root,out)
    print(out/'structure.json')


if __name__=='__main__':
    main()
