"""Next observed visit using the existing PULSE and historical split.

Both current assays are withheld. Fit and validation use the same masking.
"""
import argparse
from pathlib import Path
from types import SimpleNamespace

from external_validation import run as validation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--source', choices=validation.MODALITIES, required=True)
    args = parser.parse_args()
    source = args.source
    target = 'metab' if source == 'protein' else 'protein'
    data_dir = validation.dataset_dir(args.data_root, 'ipop')
    subjects, raw = validation.load_canonical(data_dir)
    cases = validation.eligible_cases(subjects, raw, source, target, True)
    patient_ids, groups = validation.split_cases(cases, args.seed)
    training_samples = subjects.loc[subjects.patient_id.isin(patient_ids['train']), 'sample_id']
    norm = validation.Normalizer(raw, training_samples.tolist())
    config = SimpleNamespace(seed=args.seed, device='cpu', max_history=4, epochs=30, patience=8,
                             latent_dim=16, hidden_dim=64, temporal_model='recurrent', lr=.001,
                             context_weight=.25, kl_weight=.001, align_weight=.01)
    existing_input = validation.model_input

    def both_hidden(case, normalizer, source_mod, target_mod, use_history, max_history, device):
        visits, masks, truth = existing_input(case, normalizer, source_mod, target_mod,
                                               use_history, max_history, device)
        visits[-1][source_mod] = None
        masks[-1][source_mod] = 1
        return visits, masks, truth

    validation.model_input = both_hidden
    out = args.data_root / 'forecast' / f'seed-{args.seed}' / f'{source}_to_{target}'
    out.mkdir(parents=True, exist_ok=True)
    validation._write_json(out / 'split.json', patient_ids)
    validation._write_json(out / 'config.json', vars(config))
    model = validation.fit_arm(groups, norm, source, target, True, False, config, out)
    score = validation.evaluate_arm(model, groups['test'], norm, source, target, True, config, out)
    print(args.seed, source, 'next-visit', score['n_samples'], score['pearson_median'], flush=True)


if __name__ == '__main__':
    main()
