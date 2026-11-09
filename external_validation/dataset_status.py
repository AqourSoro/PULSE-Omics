"""Check LOCAL omics download and canonical preparation readiness, without network calls.

Download-manifest state records coverage of the source discovery/download pass.
It must not be mistaken for a verified assay/visit join or training readiness.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

DATASETS = ('su2020','ipop','hmp2','covaar','csiro','p100')
ROOT = Path('/mnt/f/datasets/pulse-omics')
CANONICAL = ('subjects.csv','proteomics.csv','metabolomics.csv')


def status(root: Path) -> pd.DataFrame:
    rows=[]
    for name in DATASETS:
        raw=root/'raw'/name
        prepared=root/'processed'/'external_validation'/name
        manifest=raw/'download_manifest.json'
        raw_files = 0
        partial_files = 0
        if raw.is_dir():
            for file in raw.rglob('*'):
                if file.is_file():
                    partial_files += file.name.endswith('.part')
                    raw_files += not file.name.endswith('.part') and file.name not in ('download_manifest.json', 'download_events.jsonl')
        entry={'dataset':name, 'raw_present':raw.exists(),
               'raw_files_present':raw_files, 'partial_files':partial_files,
               'manifest_present':manifest.is_file(),
               'download_state':'no_manifest', 'files_checked_last_pass':None,
               'last_pass_issues':None, 'listed_sources':None,
               'files_downloaded_last_pass':None, 'files_existing_last_pass':None,
               'canonical_files_present':sum((prepared/f).is_file() for f in CANONICAL),
               'canonical_ready':all((prepared/f).is_file() for f in CANONICAL)}
        if manifest.is_file():
            try:
                body=json.loads(manifest.read_text(encoding='utf-8'))
                if body.get('dataset') not in (name,None):
                    raise ValueError('manifest dataset mismatch')
                records=body.get('files',[])
                entry.update(download_state=body.get('status','unreported'),
                             files_checked_last_pass=body.get('files_checked_this_run',len(records)),
                             last_pass_issues=body.get('issues'),
                             listed_sources=len(body.get('source_inventories',[])),
                             files_downloaded_last_pass=sum(x.get('status')=='downloaded' for x in records),
                             files_existing_last_pass=sum(x.get('status')=='existing' for x in records))
            except (OSError, ValueError, TypeError) as exc:
                entry['download_state']='invalid_manifest: '+str(exc)
        if name=='p100':
            entry['access_note']='public supplements only; individual-level omics controlled'
        else:entry['access_note']=''
        rows.append(entry)
    return pd.DataFrame(rows)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=ROOT)
    p.add_argument('--output',type=Path,help='Optional CSV output; no output file by default')
    a=p.parse_args(argv)
    table=status(a.data_root.expanduser())
    print('\nPULSE-Omics 数据状态')
    print('─' * 35)
    ready_count = 0
    for _, row in table.iterrows():
        name = str(row['dataset']).upper() if row['dataset'] == 'hmp2' else str(row['dataset'])
        ready = bool(row['canonical_ready'])
        ready_count += int(ready)
        prepared = int(row['canonical_files_present'])
        experimental = '可运行' if ready else '尚未准备好'
        if row['dataset'] == 'p100' and not ready:
            experimental = '受限（缺少授权的个体级组学数据）'
        print(f'{name}  |  实验：{experimental}（标准化输入 {prepared}/3）')
        state = str(row['download_state'])
        issues = row['last_pass_issues']
        if state == 'incomplete':
            note = (f'最近记录有 {int(issues)} 项异常'
                    if pd.notna(issues) else '最近记录显示未完成')
        elif state == 'enumerated_in_this_run':
            note = '最近一次检查无记录异常'
        elif state == 'unreported':
            note = '旧记录未报告状态'
        elif state == 'no_manifest':
            note = '没有下载记录'
        elif state.startswith('invalid_manifest'):
            note = '下载记录读取失败'
        else:
            note = f'下载记录状态：{state}'
        print(f"  下载：本地 {int(row['raw_files_present']):,} 个文件；"
              f"未完成片段 {int(row['partial_files'])} 个；{note}")
        print()
    print(f'实验输入齐全：{ready_count}/{len(table)} 个数据集')
    print('说明：本地文件数不代表下载完整；实验输入齐全仅表示所需文件存在。')
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        table.to_csv(a.output,index=False)
        print('Saved:',a.output)
    return 0


if __name__=='__main__':
    raise SystemExit(main())