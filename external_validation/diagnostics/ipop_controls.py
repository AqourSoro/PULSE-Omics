"""Report scientific contrasts for completed iPOP experiments without exporting participant rows."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.stats import rankdata
import torch

from external_validation import run as validation


def scaled_rank(data):
    return rankdata(data, axis=0)


def whiten(x, weights):
    center = x - np.average(x, axis=0, weights=weights)
    spread = np.sqrt(np.average(center**2, axis=0, weights=weights))
    return np.divide(center, spread, where=spread > 0, out=np.zeros_like(center)) * np.sqrt(weights[:, None])


def association(subjects, raw):
    counts = subjects.groupby('patient_id').size()
    ids = subjects.patient_id.to_numpy()
    groups = [np.flatnonzero(ids == pid) for pid in counts[counts >= 2].index]
    indices = np.concatenate(groups)
    x, y = [scaled_rank(raw[k].to_numpy()[indices]) for k in validation.MODALITIES]
    offsets = np.cumsum([0] + [len(g) for g in groups])
    slices = [slice(start, end) for start, end in zip(offsets[:-1], offsets[1:])]
    n = len(slices)
    w = np.concatenate([np.full(len(g), 1 / (n * len(g))) for g in groups])
    random = np.random.default_rng(1341)

    def strength(a, b, weights, permutations):
        a, b = whiten(a, weights), whiten(b, weights)
        observed = float(np.sqrt(np.mean((a.T @ b)**2)))
        gram_a, gram_b = a @ a.T, b @ b.T
        p, q = a.shape[1], b.shape[1]
        null = [np.sqrt(max(0., (gram_a * gram_b[np.ix_(idx, idx)]).sum() / (p * q)))
                for idx in permutations]
        return {'rms_rank_cross_correlation': observed,
                'permuted_median': float(np.median(null)),
                'permuted_2p5_97p5': np.quantile(null, [.025, .975]).tolist()}

    subject_means = [np.stack([a[s].mean(axis=0) for s in slices]) for a in (x, y)]
    between = strength(*subject_means, np.full(n, 1/n), [random.permutation(n) for _ in range(199)])
    within_data = [np.concatenate([a[s] - a[s].mean(axis=0) for s in slices]) for a in (x, y)]
    shifts = [np.concatenate([np.roll(np.arange(s.start, s.stop), random.integers(s.stop - s.start))
                              for s in slices]) for _ in range(199)]
    within = strength(*within_data, w, shifts)
    return {'patients': n, 'visits': len(indices), 'between_patient_means': between,
            'within_patient_deviations': within,
            'reference_draws': 199,
            'weights': 'equal patient weights; rank transform by feature',
            'within_reference': 'independent circular shifts inside each patient; descriptive reference'}


def predictions(file, cases, features):
    rows = pd.read_csv(file, dtype={'sample_id': str})
    matrix = rows.pivot(index='sample_id', columns='feature', values='predicted')
    return matrix.loc[[c.sample for c in cases], features].to_numpy()


def errors(pred, truth, std):
    return np.abs((pred - truth) / std).mean(axis=1)


def patient_mean(value, patients):
    return pd.Series(value).groupby(patients).mean().to_numpy()


def paired_contrast(test, reference, patients):
    effect = patient_mean(test - reference, patients)
    rng = np.random.default_rng(772)
    draws = rng.choice(effect, (5000, len(effect)), replace=True).mean(axis=1)
    return {'delta_zmae_test_minus_reference': float(effect.mean()),
            'patient_bootstrap_95_low': float(np.quantile(draws, .025)),
            'patient_bootstrap_95_high': float(np.quantile(draws, .975)),
            'patients': len(effect), 'visits': len(test)}


def within_r(actual, predicted, patients):
    x, y = (pd.DataFrame(a) for a in (actual, predicted))
    p = pd.Series(patients)
    eligible = p.map(p.value_counts()) > 1
    x, y, p = x.loc[eligible], y.loc[eligible], p.loc[eligible]
    x = x - x.groupby(p).transform('mean')
    y = y - y.groupby(p).transform('mean')
    out = x.corrwith(y).dropna()
    return float(out.median()) if len(out) else None, int(p.nunique())


def load_model(path):
    checkpoint = torch.load(path, map_location='cpu', weights_only=True)
    model = validation.PulseModel(**checkpoint['config'])
    model.load_state_dict(checkpoint['state_dict'])
    return model.eval()


def counterfactual_input(model, case, norm, source, target, *, history=None,
                         fusion=True, source_only=False, no_current=False):
    if history is not None:
        case = replace(case, history=tuple(history))
    visits, masks, _ = validation.model_input(case, norm, source, target, True, 4, 'cpu')
    if source_only:
        for visit, mask in zip(visits[:-1], masks[:-1]):
            visit[target] = None
            mask[target] = 2
    if no_current:
        visits[-1][source] = None
        masks[-1][source] = 1
    out = model(visits, masks, use_history_in_fusion=fusion)
    z = out['reconstructions'][-1][target].detach().numpy()[0]
    return z * norm.std[target] + norm.mean[target]


def donor_histories(case, test_visits, dates):
    needed = len(case.history[-4:])
    time = dates[case.sample]
    real_lags = np.array([time - dates[s] for s in case.history[-4:]])
    donor_windows = []
    for pid, samples in test_visits.items():
        if pid == case.patient:
            continue
        available = [s for s in samples if dates[s] < time]
        if len(available) < needed:
            continue
        candidates = [tuple(available[j:j+needed]) for j in range(len(available)-needed+1)]
        scores = []
        for window in candidates:
            lags = np.array([time - dates[s] for s in window])
            score = np.abs(np.log1p(lags) - np.log1p(real_lags)).mean()
            scores.append((float(score), window, int(abs(lags[-1] - real_lags[-1]))))
        donor_windows.append(min(scores))
    return sorted(donor_windows)[:3]


def latent_jeffrey(model, a, b, norm, source, target):
    arrays = []
    for mod, rows in ((source, a), (target, b)):
        z = torch.from_numpy(np.stack([norm.scaled(mod, s) for s in rows]))
        mean, logvar = model.encoders[mod](z)
        arrays.append((mean.double(), logvar.double().exp().clamp_min(1e-8)))
    (p, v), (q, w) = arrays
    return float((.5 * (v/w + w/v + (p-q).square()*(1/v + 1/w) - 2)).sum(1).median())


def evaluate(root, subjects, raw, seed, source, rows, contrasts, gap_rows, latent_rows):
    target = 'metab' if source == 'protein' else 'protein'
    dirs = list((root/'runs'/f'seed-{seed}'/'results'/'pulse_validation').glob('*/summary.csv'))
    if len(dirs) != 1:
        raise ValueError(f'seed {seed}: exactly one completed original run is required')
    task = dirs[0].parent/'ipop'/f'{source}_to_{target}'/'historical_completion'
    all_cases = validation.eligible_cases(subjects, raw, source, target, True)
    ids, groups = validation.split_cases(all_cases, seed)
    if ids != json.loads((task/'split.json').read_text()):
        raise ValueError('split differs from the original result')
    samples = subjects.loc[subjects.patient_id.isin(ids['train']), 'sample_id'].tolist()
    norm = validation.Normalizer(raw, samples)
    if norm.names != json.loads((task/'features.json').read_text()):
        raise ValueError('feature set differs from the original result')
    cases = groups['test']
    patients = np.asarray([c.patient for c in cases])
    truth = np.stack([norm.raw(target, c.sample) for c in cases])
    last = np.stack([norm.raw(target, c.history[-1]) for c in cases])
    std = norm.std[target]
    means = np.broadcast_to(norm.mean[target], truth.shape)
    values = {'training_mean': means, 'last_observation': last}
    metrics = {}
    for arm in ('pulse_full', 'no_history_input', 'no_alignment'):
        values[arm] = predictions(task/arm/'predictions.csv.gz', cases, norm.names[target])
        metrics[arm] = json.loads((task/arm/'metrics.json').read_text())
    fpath = root/'forecast'/f'seed-{seed}'/f'{source}_to_{target}'
    if json.loads((fpath/'split.json').read_text()) != ids:
        raise ValueError('forecast split differs from historical-completion split')
    values['forecast'] = predictions(fpath/'predictions.csv.gz', cases, norm.names[target])
    metrics['forecast'] = json.loads((fpath/'metrics.json').read_text())
    errors_by_arm = {}
    for arm, pred in values.items():
        e = errors(pred, truth, std)
        errors_by_arm[arm] = e
        w, n = within_r(truth, pred, patients)
        r = metrics.get(arm, {}).get('pearson_median')
        if arm == 'last_observation':
            r = json.loads((task/arm/'metrics.json').read_text())['pearson_median']
        rows.append({'seed': seed, 'direction': f'{source}_to_{target}', 'arm': arm,
                     'patients': len(set(patients)), 'visits': len(cases), 'target_features': len(std),
                     'zmae_equal_patient': float(patient_mean(e, patients).mean()),
                     'zrmse_equal_patient': float(np.sqrt(patient_mean((((pred-truth)/std)**2).mean(1), patients).mean())),
                     'median_r_across_visits': r,
                     'median_r_within_patient': w, 'repeated_test_patients': n})
    for arm, ref in (('pulse_full','no_history_input'), ('pulse_full','no_alignment'),
                     ('pulse_full','last_observation'), ('no_history_input','training_mean'),
                     ('forecast','last_observation'), ('forecast','training_mean')):
        contrasts.append({'seed': seed, 'direction': f'{source}_to_{target}', 'test': arm, 'reference': ref,
                          **paired_contrast(errors_by_arm[arm], errors_by_arm[ref], patients)})
    date_series = pd.to_datetime(subjects.date).astype('int64') // 86400000000000
    dates = dict(zip(subjects.sample_id, date_series))
    test_visits = {pid: g.sample_id.tolist() for pid, g in subjects[subjects.patient_id.isin(ids['test'])].groupby('patient_id')}
    gaps = np.array([dates[c.sample] - dates[c.history[-1]] for c in cases])
    for label, keep in (('1–30 d', gaps<=30), ('31–120 d', (gaps>30)&(gaps<=120)), ('>120 d',gaps>120)):
        if keep.any():
            for arm, ref in (('pulse_full','no_history_input'),('forecast','last_observation')):
                gap_rows.append({'seed': seed, 'direction': f'{source}_to_{target}', 'gap': label,
                                 'test': arm, 'reference': ref,
                                 **paired_contrast(errors_by_arm[arm][keep], errors_by_arm[ref][keep], patients[keep])})
    model = load_model(task/'pulse_full'/'best.pt')
    without_alignment = load_model(task/'no_alignment'/'best.pt')
    forecast_model = load_model(fpath/'best.pt')
    # Prediction effects for trained weights and identical targets.
    controls = {k: [] for k in ('no_fusion', 'past_source_only', 'reverse_order', 'one_of_four',
                               'two_of_four', 'wrong_patient', 'no_current_inference',
                               'forecast_reverse_order', 'forecast_wrong_patient')}
    age_mismatch = []
    with torch.no_grad():
        for i, c in enumerate(cases):
            def append(name, **kwargs):
                p = counterfactual_input(model, c, norm, source, target, **kwargs)
                controls[name].append((i, float(errors(p[None], truth[i:i+1], std)[0])))
            append('no_fusion', fusion=False)
            append('past_source_only', source_only=True)
            append('no_current_inference', no_current=True)
            if len(c.history) >= 2:
                append('reverse_order', history=c.history[-4:][::-1])
                p = counterfactual_input(forecast_model, c, norm, source, target,
                                         history=c.history[-4:][::-1], no_current=True)
                controls['forecast_reverse_order'].append((i, float(errors(p[None], truth[i:i+1], std)[0])))
            if len(c.history) >= 4:
                append('one_of_four', history=c.history[-1:])
                append('two_of_four', history=c.history[-2:])
            donors = donor_histories(c, test_visits, dates)
            if donors:
                risks, forecast_risks = [], []
                for _, history, lag_diff in donors:
                    p = counterfactual_input(model, c, norm, source, target, history=history)
                    risks.append(float(errors(p[None], truth[i:i+1], std)[0]))
                    p = counterfactual_input(forecast_model, c, norm, source, target,
                                             history=history, no_current=True)
                    forecast_risks.append(float(errors(p[None], truth[i:i+1], std)[0]))
                    age_mismatch.append(lag_diff)
                controls['wrong_patient'].append((i, np.mean(risks)))
                controls['forecast_wrong_patient'].append((i, np.mean(forecast_risks)))
        current_keys = [c.sample for c in cases]
        previous_keys = [c.history[-1] for c in cases]
        for arm, fitted in (('pulse_full', model), ('no_alignment', without_alignment)):
            latent_rows.append({'seed': seed, 'direction': f'{source}_to_{target}', 'arm': arm,
                                'same_date_jeffreys_median': latent_jeffrey(fitted,current_keys,current_keys,norm,source,target),
                                'previous_same_subject_jeffreys_median': latent_jeffrey(fitted,current_keys,previous_keys,norm,source,target)})
    for control, pairs in controls.items():
        if not pairs:
            continue
        indices, risk = zip(*pairs)
        ii = np.asarray(indices)
        reference = 'forecast' if control.startswith('forecast_') else 'pulse_full'
        result = {'seed': seed, 'direction': f'{source}_to_{target}', 'test': control, 'reference': reference,
                  **paired_contrast(np.asarray(risk), errors_by_arm[reference][ii], patients[ii])}
        if control.endswith('wrong_patient'):
            result['median_donor_last_gap_difference_days'] = float(np.median(age_mismatch))
        contrasts.append(result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('data_root', type=Path)
    parser.add_argument('--association-only', action='store_true')
    args = parser.parse_args()
    torch.set_num_threads(1)
    root = args.data_root
    out = root/'aggregate'
    out.mkdir(exist_ok=True)
    subjects, raw = validation.load_canonical(validation.dataset_dir(root,'ipop'))
    result = association(subjects, raw)
    validation._write_json(out/'association.json', result)
    print(json.dumps(result, indent=2), flush=True)
    if args.association_only:
        return
    rows, contrasts, gaps, latents = [], [], [], []
    for seed in (42, 137, 509):
        for source in validation.MODALITIES:
            evaluate(root, subjects, raw, seed, source, rows, contrasts, gaps, latents)
            print('Evaluated', seed, source, flush=True)
    for name, table in (('scores',rows), ('paired_contrasts',contrasts), ('by_gap',gaps), ('alignment',latents)):
        pd.DataFrame(table).to_csv(out/f'{name}.csv', index=False)


if __name__ == '__main__':
    main()
