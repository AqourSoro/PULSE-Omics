"""Audit the CSIRO group coverage and score sensitivity to the published CCM cutoff."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from external_validation.cohort_join_audit import _find_excel
from external_validation.run import _write_json, dataset_dir


def audit(root, runs):
    folder = dataset_dir(root,'csiro')
    protein = pd.read_csv(folder/'proteomics.csv').set_index('sample_id')
    metabolite = pd.read_csv(folder/'metabolomics.csv').set_index('sample_id')
    visits = pd.read_csv(folder/'subjects.csv').set_index('sample_id')
    source = _find_excel(root/'raw/csiro','006_Task')
    meta = pd.read_excel(source,sheet_name='COVID HPV CCM')
    barcode = 'Sample IDS (Barcode)'
    label = meta.drop_duplicates(barcode,keep=False).set_index(barcode)['Vaccine Received']
    vaccine = label.loc[visits.index.str.removeprefix('csiro_')].set_axis(visits.index)
    by_vaccine = {}
    for name, ids in vaccine.groupby(vaccine).groups.items():
        p,m = protein.loc[ids],metabolite.loc[ids]
        by_vaccine[name] = {
            'visits':len(ids), 'participants':int(visits.loc[ids,'patient_id'].nunique()),
            'by_visit_ordinal':visits.loc[ids,'visit_order'].value_counts().sort_index().to_dict(),
            'protein_features_with_measurement':int(p.notna().any().sum()),
            'protein_finite_cells':int(p.notna().sum().sum()), 'metabolite_finite_cells':int(m.notna().sum().sum())}
    scores = []
    for file in sorted(runs.rglob('protein_to_metab/heldout.npz')):
        if file.parents[2].name != 'csiro':
            continue
        features = json.loads((file.parent/'features.json').read_text())['metab']
        with np.load(file, allow_pickle=False) as d:
            raw = metabolite.loc[d['sample_id'],features].to_numpy()
            valid = d['valid']
            above_cutoff = valid & (raw>=300)
            metrics = {}
            for arm in ('mean','current','history','forecast','persistence'):
                error = np.abs(d[arm]-d['truth'])
                metrics[arm] = {}
                for name, keep in (('all_numeric_zmae',valid),('raw_ge_300_zmae',above_cutoff)):
                    n = keep.sum(1)
                    per_visit = (error[n>0]*keep[n>0]).sum(1)/n[n>0]
                    metrics[arm][name] = float(pd.Series(per_visit).groupby(d['patient_id'][n>0]).mean().mean())
            scores.append({'seed':int(file.parents[1].name.removeprefix('seed-')),
                           'cells_below_300':int((valid&~above_cutoff).sum()),'metrics':metrics})
    return {'ccm_numeric_threshold_on_release_values':{
                'rule':'published CCM threshold 300, applied here only to numeric release values',
                'finite_cells':int(metabolite.notna().sum().sum()),
                'finite_cells_below_300':int(metabolite.lt(300).sum().sum()),
                'features_with_a_finite_value_below_300':int(metabolite.lt(300).any().sum())},
            'by_vaccine':by_vaccine, 'score_sensitivity':scores}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--runs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    _write_json(args.output,audit(args.data_root,args.runs))
