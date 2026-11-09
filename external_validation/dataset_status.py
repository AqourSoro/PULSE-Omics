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
        entry={'dataset':name, 'raw_present':raw.exists(),
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
    print(table.to_string(index=False,na_rep='—'))
    print('\nNOTE: enumerated_in_this_run means the last source enumeration/download pass had no recorded issues; it does not prove that raw data are fully parsed, clinically paired, or ready for model training.')
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        table.to_csv(a.output,index=False)
        print('Saved:',a.output)
    return 0


if __name__=='__main__':
    raise SystemExit(main())