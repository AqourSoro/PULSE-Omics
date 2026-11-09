"""Inspect the real local iPOP/COVAAR/CSIRO release structures, without training.

Retains original files, does not infer specimen/patient/visit correspondence.
Report contains field names and counts, not participant/sample identifiers.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

ROOT=Path('/mnt/f/datasets/pulse-omics')
COHORTS=('ipop','covaar','csiro')
CANDIDATE=('subject','patient','participant','person','time','week','date','visit',
           'sample','specimen','id','hours','hour','day','phase')


def _text_columns(stream, name: str):
    # GEO series matrices use annotation lines before the actual expression table.
    delim = '\t' if '.tsv' in name or '.txt' in name else ','
    first, line_count = [], 0
    for line in stream:
        line_count += 1
        if line_count > 200:
            break
        if not line.strip():
            continue
        row=next(csv.reader([line], delimiter=delim))
        if len(row)>1:
            first.append({'line':line_count,'width':len(row),'columns':[str(x)[:110] for x in row[:20]]})
            # Preserve only the first tabular header; never serialize sample rows.
            break
    return first


def _file_info(path: Path, relative: Path) -> dict:
    n=path.name.lower()
    out={'path':relative.as_posix(), 'bytes':path.stat().st_size}
    try:
        if n.endswith(('.csv','.tsv','.txt','.csv.gz','.tsv.gz','.txt.gz')):
            op=gzip.open if n.endswith('.gz') else open
            with op(path,'rt',encoding='utf-8-sig',errors='replace',newline='') as f:
                preview=_text_columns(f,n)
            out.update(kind='delimited', headers=preview)
        elif n.endswith(('.xlsx','.xlsm')):
            from openpyxl import load_workbook
            book=load_workbook(path,read_only=True,data_only=True)
            try:
                tabs=[]
                for sh in book:
                    headers=[]
                    for row in sh.iter_rows(min_row=1,max_row=min(sh.max_row,8),max_col=min(sh.max_column,32),values_only=True):
                        values=[str(v)[:110] if v is not None else '' for v in row]
                        if sum(bool(v) for v in values)>=2:
                            headers.append(values)
                        # Stop at the first candidate header; do not record participant data.
                        break
                    tabs.append({'name':sh.title,'rows':sh.max_row,'columns':sh.max_column,'first_rows':headers})
                out.update(kind='excel',sheets=tabs)
            finally:
                book.close()
        elif n.endswith('.zip'):
            with ZipFile(path) as z:
                names=[{'name':i.filename,'bytes':i.file_size} for i in z.infolist() if not i.is_dir()]
                out.update(kind='zip',n_members=len(names),members=names[:80])
        elif n.endswith(('.rdata','.rda','.rds')):
            out.update(kind='r_object',note='iPOP source objects inspected separately when rdata is installed')
        elif n.endswith(('.pdf','.html')):
            out.update(kind='document')
        else:
            out.update(kind='binary_or_other')
    except Exception as exc:
        out.update(kind='unreadable',error=f'{type(exc).__name__}: {exc}')
    return out


def _ipop_rdata_structure(path: Path) -> dict:
    try:
        from external_validation.scripts.inspect_restored_omics import read_ipop
    except ImportError as e:
        return {'status':'rdata_dependency_unavailable','error':str(e)}
    try:
        obj=read_ipop(path)
    except Exception as e:
        return {'status':'cannot_read_rdata','error':f'{type(e).__name__}: {e}'}
    items={}
    for name, value in obj.items():
        if isinstance(value,pd.DataFrame):
            keys=value.index.astype(str)
            items[name]={'rows':len(value),'columns':value.shape[1],
                         'index_unique':bool(value.index.is_unique),
                         'columns_candidate_ids':[str(c) for c in value.columns if any(w in str(c).lower() for w in CANDIDATE)][:40],
                         'index_matches':{}}
            for other, table in obj.items():
                if other==name or not isinstance(table,pd.DataFrame):continue
                items[name]['index_matches'][other]=len(set(keys)&set(table.index.astype(str)))
        else:
            items[name]={'type':type(value).__name__}
    return {'status':'inspected','objects':items,
            'note':'Index intersection alone is not an established biological sample join.'}


def probe(root: Path, dataset: str, *, max_files: int = 300) -> dict:
    if dataset not in COHORTS:raise ValueError('Unsupported cohort')
    folder=root/'raw'/dataset
    if not folder.is_dir():raise FileNotFoundError(folder)
    paths=sorted(p for p in folder.rglob('*') if p.is_file() and not p.name.endswith('.part')
                 and p.name not in {'download_manifest.json','download_events.jsonl','data_inventory.json'})
    result={'dataset':dataset,'files':len(paths),'inspected_files':min(len(paths),max_files),
            'source_files':[_file_info(p,p.relative_to(folder)) for p in paths[:max_files]],
            'canonical_status':'not_verified',
            'required_for_export':['explicit protein and metabolite sample identities',
                                   'verified cross-assay specimen crosswalk',
                                   'patient identity and numeric visit order']}
    if dataset=='ipop':
        rdata=folder/'seasonal'/'Multi_Omics_Seasonal.RData'
        if rdata.is_file():result['rdata']=_ipop_rdata_structure(rdata)
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,default=Path(os.environ.get('PULSE_OMICS_DATA_ROOT',ROOT)))
    parser.add_argument('--datasets',nargs='+',choices=COHORTS,default=list(COHORTS))
    parser.add_argument('--max-files',type=int,default=300)
    args=parser.parse_args(argv)
    if args.max_files<1:parser.error('max-files must be >=1')
    for ds in dict.fromkeys(args.datasets):
        try:
            report=probe(args.data_root,ds,max_files=args.max_files)
        except (FileNotFoundError,ValueError) as e:
            print(f'{ds}: {e}')
            continue
        out=args.data_root/'processed'/'external_validation'/ds/'source_profile.json'
        out.parent.mkdir(parents=True,exist_ok=True)
        out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print(f'{ds}: source files {report["files"]}, scanned {report["inspected_files"]}, '
              f'schema={out}')
        for item in report['source_files']:
            kind = item['kind']
            if kind not in {'excel','delimited','zip','r_object','unreadable'}:
                continue
            print(f'  {item["path"]}: {kind}')
            if kind == 'delimited':
                for info in item.get('headers', [])[:1]:
                    print('    header:', ', '.join(info['columns'][:14]))
            elif kind == 'excel':
                for sheet in item.get('sheets', [])[:12]:
                    print(f'    {sheet["name"]}: {sheet["rows"]} rows x {sheet["columns"]} cols')
                    if sheet.get('first_rows'):
                        print('      header:', ', '.join(sheet['first_rows'][0][:12]))
            elif kind == 'zip':
                print('    members:', ', '.join(v['name'] for v in item.get('members', [])[:8]))
        if 'rdata' in report:
            details = report['rdata']
            print(f'  iPOP RData: {details["status"]}, {len(details.get("objects",{}))} objects')
            for name, obj in details.get('objects', {}).items():
                if 'rows' in obj:
                    overlap = sorted(obj['index_matches'].items(), key=lambda x: -x[1])[:2]
                    print(f'    {name}: {obj["rows"]} x {obj["columns"]}; '
                          f'row-index overlaps={overlap}')
    return 0

if __name__=='__main__':
    raise SystemExit(main())