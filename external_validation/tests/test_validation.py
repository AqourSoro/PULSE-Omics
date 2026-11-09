"""Deterministic tests with a small model double; no clinical data or GPU needed."""
from __future__ import annotations

import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import torch


class SmallModel(torch.nn.Module):
    def __init__(self, modality_dims, **kwargs):
        super().__init__()
        self.dims = modality_dims
        self.weights = torch.nn.ParameterDict({
            dest: torch.nn.Parameter(torch.ones(max(modality_dims.values()), modality_dims[dest]) * .02)
            for dest in modality_dims
        })

    def get_config(self):
        return {'dummy': True}

    def forward(self, visits_data, missing_masks, **kwargs):
        outputs = []
        for visit in visits_data:
            x = next(v for v in visit.values() if v is not None)
            length = x.shape[1]
            output = {}
            for dest, weight in self.weights.items():
                padded = torch.nn.functional.pad(x, (0,weight.shape[0]-length))
                output[dest] = padded @ weight
            outputs.append(output)
        base = torch.sum(self.weights['protein'] ** 2) * 0.0
        return {'reconstructions': outputs, 'losses': {
            'recon_loss': base, 'kl_loss': base, 'alignment_loss': base}}


pulse = types.ModuleType('pulse')
model_module = types.ModuleType('pulse.model')
model_module.PulseModel = SmallModel
pulse.model = model_module
sys.modules.setdefault('pulse', pulse)
sys.modules.setdefault('pulse.model', model_module)
from external_validation import run as experiment


def make_canonical(folder, n=24):
    folder.mkdir(parents=True, exist_ok=True)
    subjects, protein, metab = [], [], []
    for i in range(n):
        for visit in (1, 2):
            sid=f'p{i:03d}-v{visit}'
            subjects.append({'sample_id':sid,'patient_id':f'p{i:03d}','visit_order':visit})
            protein.append({'sample_id':sid,'p1':float(i+visit),'p2':float(i*i/10+visit),'p3':float(i+3*visit)})
            metab.append({'sample_id':sid,'m1':float(i+visit*2),'m2':float(i*i/11+visit),'m3':float(i+4*visit)})
    pd.DataFrame(subjects).to_csv(folder/'subjects.csv',index=False)
    pd.DataFrame(protein).to_csv(folder/'proteomics.csv',index=False)
    pd.DataFrame(metab).to_csv(folder/'metabolomics.csv',index=False)


class ExternalValidationTests(unittest.TestCase):
    def test_split_grouping_and_missing_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            make_canonical(folder)
            subjects,raw=experiment.load_canonical(folder)
            cases=experiment.eligible_cases(subjects,raw,'protein','metab',True)
            self.assertEqual(len(cases),24)
            splits, groups=experiment.split_cases(cases,42)
            self.assertFalse(set(splits['train']) & set(splits['test']))
            train_samples=subjects.loc[subjects.patient_id.isin(splits['train']),'sample_id'].tolist()
            norm=experiment.Normalizer(raw,train_samples)
            input_visits,masks,_=experiment.model_input(cases[0],norm,'protein','metab',True,4,'cpu')
            self.assertIsNone(input_visits[-1]['metab'])
            self.assertEqual(masks[-1]['metab'],1)
            self.assertEqual(masks[0]['metab'],0)
            self.assertEqual(len(input_visits),2)

    def test_metrics_include_per_feature_r_and_change(self):
        y=np.array([[0.,2.],[1.,3.],[2.,4.],[3.,5.]])
        detail,summary=experiment.metric_tables(y,y+1,['a','b'])
        self.assertAlmostEqual(summary['pearson_mean'],1)
        self.assertAlmostEqual(summary['mae_micro'],1)
        self.assertEqual(detail.n.tolist(),[4,4])

    def test_hmp2_sparse_feature_selection_is_train_only_and_bounded(self):
        ids=[f's{i}' for i in range(20)]
        train_ids=ids[:10]
        x=np.stack([
            [float(10+i), float(i%3+1), float(i%2), float(i+1) if i>=10 else 0.]
            for i in range(20)
        ])
        frame=pd.DataFrame(x, index=ids, columns=['common','moderate','rare','test_only'])
        raw={'protein':frame,'metab':frame.copy()}
        policy=experiment.feature_policy('hmp2',hmp2_min_nonzero_fraction=.5,
                                         hmp2_max_features=2)
        norm=experiment.Normalizer(raw,train_ids,minimum=3,**policy)
        self.assertEqual(len(norm.names['protein']),2)
        self.assertNotIn('test_only',norm.names['protein'])
        self.assertEqual(norm.selection['protein']['maximum_features'],2)
        self.assertEqual(norm.selection['protein']['minimum_nonzero_training_samples'],5)
        self.assertEqual(experiment.feature_policy('su2020')['max_features'],0)

    def test_run_end_to_end_and_isolation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            folder=experiment.dataset_dir(root,'su2020')
            make_canonical(folder)
            code=experiment.main(['--data-root',str(root),'--datasets','su2020',
                                  '--epochs','1','--patience','1','--device','cpu',
                                  '--min-feature-observations','3'])
            self.assertEqual(code,0)
            runs=list((root/'results'/'pulse_validation').iterdir())
            self.assertEqual(len(runs),1)
            results=pd.read_csv(runs[0]/'summary.csv')
            self.assertEqual((results.status=='completed').sum(),12)
            self.assertTrue((runs[0]/'report.md').is_file())
            self.assertTrue((runs[0]/'ablation_comparison.csv').is_file())
            arm=runs[0]/'su2020'/'protein_to_metab'/'historical_completion'/'pulse_full'
            self.assertTrue((arm/'metrics.json').is_file())
            self.assertTrue((arm/'predictions.csv.gz').is_file())
            self.assertTrue((arm/'best.pt').is_file())
            ridge=runs[0]/'su2020'/'protein_to_metab'/'historical_completion'/'ridge_current'
            self.assertTrue((ridge/'feature_metrics.csv').is_file())
            self.assertTrue((ridge/'baseline_config.json').is_file())
            self.assertTrue((folder/'subjects.csv').is_file())
            self.assertFalse((root/'raw').exists())


if __name__=='__main__':
    unittest.main()