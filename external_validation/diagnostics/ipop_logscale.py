"""Seed-42 sensitivity of protein-to-metabolite reconstruction to intensity scale."""
import argparse
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from external_validation import run as validation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_root', type=Path)
    args = parser.parse_args()
    root = args.data_root
    subjects, raw = validation.load_canonical(validation.dataset_dir(root, 'ipop'))
    if not (raw['metab'].to_numpy() > 0).all():
        raise ValueError('iPOP log10 sensitivity requires the observed positive intensity matrix')
    source, target = 'protein', 'metab'
    raw['metab'] = np.log10(raw['metab'])
    cases = validation.eligible_cases(subjects, raw, source, target, True)
    ids, groups = validation.split_cases(cases, 42)
    training_samples = subjects.loc[subjects.patient_id.isin(ids['train']), 'sample_id']
    norm = validation.Normalizer(raw, training_samples.tolist())
    config = SimpleNamespace(seed=42, device='cpu', max_history=4, epochs=30, patience=8,
                             latent_dim=16, hidden_dim=64, temporal_model='recurrent', lr=.001,
                             context_weight=.25, kl_weight=.001, align_weight=.01)
    out = root / 'scale-sensitivity' / 'seed-42'
    validation._write_json(out / 'split.json', ids)
    validation._write_json(out / 'features.json', norm.names)
    for arm, history in (('pulse_full', True), ('no_history_input', False)):
        folder = out / arm
        model = validation.fit_arm(groups, norm, source, target, history, False, config, folder)
        metrics = validation.evaluate_arm(model, groups['test'], norm, source, target, history, config, folder)
        print('log10 metabolite intensity', arm, metrics['pearson_median'], flush=True)


if __name__ == '__main__':
    main()
