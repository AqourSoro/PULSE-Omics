"""Link the COVAAR SomaScan and Figshare releases at documented trial events."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
from pathlib import Path
import re

import numpy as np
import pandas as pd

from external_validation.run import _write_json, dataset_dir


ARM = {'Vaccine': 'V', 'Placebo': 'P', 'Booster': 'V_B'}
TIME = {'Baseline': 0, '35 min': .5, '2 hour': 2, '6 hour': 6, '24 hour': 24}
STAGE = {'Day 1': 0, 'Day 2': 1, 'Month 5': 2}


def read_geo(path):
    metadata = {}
    with gzip.open(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        reader = csv.reader(stream, delimiter='\t')
        for row in reader:
            if not row:
                continue
            if row[0] == '!series_matrix_table_begin':
                break
            if row[0] == '!Sample_geo_accession':
                accessions = row[1:]
            if row[0] == '!Sample_characteristics_ch1':
                label, _, _ = row[1].partition(':')
                metadata[label] = [value.split(':', 1)[1].strip() for value in row[1:]]
        header = next(reader)
        if header != ['ID_REF', *accessions]:
            raise ValueError('COVAAR GEO matrix header differs from metadata accessions')
        records = []
        for row in reader:
            if not row:
                continue
            if row[0] == '!series_matrix_table_end':
                break
            records.append(row)
    features = pd.DataFrame.from_records(records, columns=header).set_index('ID_REF').astype(float)
    events = pd.DataFrame(metadata, index=accessions)
    return events, features


def exact_rfu_columns(geo, rfu):
    if not geo.index.is_unique or not rfu.index.is_unique or set(geo.index) != set(rfu.index):
        raise ValueError('COVAAR RFU and GEO probe identifiers are not a common unique set')
    rfu = rfu.loc[geo.index]
    rfu_vectors = {vec.tobytes(): sample for sample, vec in zip(rfu.columns, rfu.to_numpy().T)}
    geo_vectors = [vec.tobytes() for vec in geo.to_numpy().T]
    if (len(rfu_vectors) != rfu.shape[1] or len(set(geo_vectors)) != len(geo_vectors)
            or set(rfu_vectors) != set(geo_vectors)):
        raise ValueError('COVAAR GEO and RFU columns do not correspond by exact probe vectors')
    return dict(zip(geo.columns, (rfu_vectors[key] for key in geo_vectors)))


def description_times(description):
    label = description.lower().replace('_', ' ')
    match = re.search(r'(\d+)\s*(min|hour|hr)', label)
    times = {0} if 'baseline' in label else set()
    if match:
        times.add(.5 if match.groups() == ('35', 'min') else int(match[1]))
    return times


def prepare(root, *, inspect_only=False):
    raw = root / 'raw/covaar'
    geo = raw / 'geo'
    fig = raw / 'figshare'
    events, protein = read_geo(geo / 'GSE325178_series_matrix.txt.gz')
    rfu = pd.read_csv(geo / 'GSE325178_SomaScan_RFUs_allSamples.txt.gz', sep='\t').set_index('AptName')
    rfu = rfu.loc[:, rfu.columns.str.fullmatch(r'VA\d+')].astype(float)
    rfu_ids = exact_rfu_columns(protein, rfu)
    metabolite = pd.read_csv(fig / 'Metabolomics_Data.csv').set_index('Metabolite').astype(float)
    annotations = pd.read_csv(fig / 'Metabolite_Meta_Data.csv')
    meta = pd.read_csv(fig / 'Sample_Meta_Data.csv', keep_default_na=False)
    keys = ['Patient', 'Vax', 'Time']
    if (meta.duplicated(keys).any() or meta.Sample_ID.duplicated().any()
            or not metabolite.index.is_unique or annotations.Metabolite.duplicated().any()
            or set(meta.Sample_ID) != set(metabolite.columns)
            or set(annotations.Metabolite) != set(metabolite.index)):
        raise ValueError('COVAAR Figshare identifiers do not resolve uniquely')
    trials = meta.set_index(keys)
    people = events.individual.str.extract(r'^ID(\d+)$', expand=False).astype(int)
    if set(people) != set(meta.Patient):
        raise ValueError('COVAAR GEO and Figshare participant labels differ after removing the ID prefix')

    pairs, excluded, counts, used_metabolite = [], [], Counter(), set()
    for accession, event in events.iterrows():
        patient = people.loc[accession]
        if event.day == 'Day 1/Day 2':
            if event.timepoint != '24 hour/Baseline':
                raise ValueError(f'COVAAR unexpected combined event for {accession}')
            first, second = event.treatment.split('/')
            matched = [(patient, ARM[first], 24), (patient, ARM[second], 0)]
            stage, within_stage = 1, 0
        else:
            matched = [(patient, ARM[event.treatment], TIME[event.timepoint])]
            stage, within_stage = STAGE[event.day], TIME[event.timepoint]
        ids = [trials.loc[key, 'Sample_ID'] for key in matched if key in trials.index]
        if len(ids) != len(matched):
            reason = 'no_complete_figshare_event_key'
        elif any(key[2] not in description_times(trials.loc[key, 'Description']) for key in matched):
            reason = 'figshare_description_time_disagreement'
        elif len(ids) == 2 and not np.array_equal(metabolite[ids[0]], metabolite[ids[1]], equal_nan=True):
            reason = 'dual_event_metabolite_vectors_differ'
        else:
            reason = None
        if reason:
            counts[reason] += 1
            excluded.append({'sample_id':accession, 'protein_assay_id':rfu_ids[accession],
                             'figshare_assay_ids':'|'.join(ids), 'reason':reason})
            continue
        if used_metabolite.intersection(ids):
            raise ValueError('COVAAR Figshare column maps to multiple GEO events')
        used_metabolite.update(ids)
        counts['duplicated_dual_event' if len(ids) == 2 else 'single_event'] += 1
        # Visit order is ordinal. Month 5 has no participant-specific elapsed time.
        pairs.append({'sample_id':accession,'patient_id':f'ID{patient}',
                      'visit_order':stage*5 + list(TIME.values()).index(within_stage),
                      'study_stage':event.day, 'intervention':event.treatment,
                      'within_stage_hours_nominal':within_stage,
                      'geo_timepoint':event.timepoint, 'protein_assay_id':rfu_ids[accession],
                      'metabolite_assay_id':ids[0],
                      'metabolite_dual_event_assay_id':ids[1] if len(ids) == 2 else ''})
    paired = pd.DataFrame(pairs).sort_values(['patient_id','visit_order'])
    if paired.duplicated(['patient_id','visit_order']).any():
        raise ValueError('COVAAR duplicate linked physical visits')
    visits = paired.groupby('patient_id').size()
    p = protein[paired.sample_id].T
    m = metabolite[paired.metabolite_assay_id].T.set_axis(p.index)
    report = {'dataset':'covaar', 'geo_protein_probes':len(protein),
              'geo_protein_samples':protein.shape[1], 'geo_rfu_exact_unique_vectors':len(rfu_ids),
              'figshare_metabolite_features':len(metabolite),
              'figshare_metabolite_samples':metabolite.shape[1],
              'participant_label_intersection':len(set(people)), 'paired_visits':len(paired),
              'paired_participants':len(visits), 'paired_participants_ge2':int((visits>=2).sum()),
              'paired_participants_ge3':int((visits>=3).sum()), 'paired_participants_ge4':int((visits>=4).sum()),
              'paired_visits_per_participant':dict(sorted(Counter(map(int,visits)).items())),
              'event_resolution_counts':dict(counts),
              'figshare_columns_unused':metabolite.shape[1]-len(used_metabolite),
              'paired_metabolite_missing':int(m.isna().sum().sum()),
              'paired_metabolite_negative':int(m.lt(0).sum().sum()),
              'paired_protein_missing':int(p.isna().sum().sum()),
              'pairing_basis':'same participant and unique trial event; specimen-level identity unavailable',
              'time_basis':'GEO study stage and native timepoint, with Figshare Time 0.5 as nominal 35 min bin'}
    out = dataset_dir(root, 'covaar')
    _write_json(out/'event_link_report.json', report)
    if not inspect_only:
        paired.to_csv(out/'subjects.csv', index=False)
        p.to_csv(out/'proteomics.csv', index_label='sample_id')
        m.to_csv(out/'metabolomics.csv', index_label='sample_id')
        pd.DataFrame(excluded).to_csv(out/'unmatched_geo_events.csv', index=False)
        meta.loc[~meta.Sample_ID.isin(used_metabolite)].to_csv(out/'unmatched_figshare_events.csv', index=False)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, default=Path('/mnt/f/datasets/pulse-omics'))
    parser.add_argument('--inspect-only', action='store_true')
    args = parser.parse_args()
    print(prepare(args.data_root, inspect_only=args.inspect_only))
