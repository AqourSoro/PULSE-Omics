"""Strict source mapping tests; only synthetic fixtures, never scientific data."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from external_validation.prepare_tabular import prepare
from external_validation.source_probe import probe


class PrepareTabularTests(unittest.TestCase):
    def setup_data(self, root):
        folder=root/'raw'/'csiro'
        folder.mkdir(parents=True)
        patients=[f'S{i}' for i in range(12)]
        pairs=[(p,v,f'protein-{p}-{v}',f'metab-{p}-{v}') for p in patients for v in (1,2)]
        pd.DataFrame([{'subject':p,'time':v,'protein_ref':pid,'metab_ref':mid}
                      for p,v,pid,mid in pairs]).to_csv(folder/'mapping.csv',index=False)
        pd.DataFrame([{'specimen':pid,'protA':v+1,'protB':float(i)}
                      for i,(p,v,pid,mid) in enumerate(pairs)]).to_csv(folder/'protein.csv',index=False)
        pd.DataFrame([{'specimen':mid,'metabA':v+2,'metabB':float(i**2)}
                      for i,(p,v,pid,mid) in enumerate(pairs)]).to_csv(folder/'metab.csv',index=False)
        conf={'dataset':'csiro',
              'protein':{'path':'raw/csiro/protein.csv','sample_id_column':'specimen',
                         'feature_columns':['protA','protB']},
              'metab':{'path':'raw/csiro/metab.csv','sample_id_column':'specimen',
                       'feature_columns':['metabA','metabB']},
              'sample_crosswalk':{'path':'raw/csiro/mapping.csv','protein_sample_id':'protein_ref',
                                  'metab_sample_id':'metab_ref','patient_id':'subject','visit_order':'time'}}
        return folder,conf

    def test_strict_crosswalk_and_longitudinal_export(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            folder,config=self.setup_data(root)
            summary=prepare(root,config,inspect_only=True)
            self.assertEqual(summary['verified_paired_samples'],24)
            self.assertEqual(summary['eligible_historical'],12)
            self.assertFalse((root/'processed').exists())
            report=prepare(root,config)
            self.assertEqual(report['status'],'canonical_raw_omics_exported')
            out=root/'processed/external_validation/csiro'
            subjects=pd.read_csv(out/'subjects.csv')
            protein=pd.read_csv(out/'proteomics.csv')
            metab=pd.read_csv(out/'metabolomics.csv')
            self.assertEqual(len(subjects),24)
            self.assertEqual(protein.shape,(24,3))
            self.assertEqual(metab.shape,(24,3))
            self.assertEqual(subjects.groupby('patient_id').size().nunique(),1)
            self.assertTrue((out/'preparation_report.json').is_file())
            self.assertEqual(len(probe(root,'csiro')['source_files']),3)

    def test_ambiguous_crosswalk_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); folder,config=self.setup_data(root)
            table=pd.read_csv(folder/'mapping.csv')
            table.loc[1,'protein_ref']=table.loc[0,'protein_ref']
            table.to_csv(folder/'mapping.csv',index=False)
            with self.assertRaisesRegex(ValueError,'duplicate assay specimen'):
                prepare(root,config)

    def test_same_patient_visit_duplicate_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); folder,config=self.setup_data(root)
            table=pd.read_csv(folder/'mapping.csv')
            table.loc[1,'time']=table.loc[0,'time']
            table.to_csv(folder/'mapping.csv',index=False)
            with self.assertRaisesRegex(ValueError,'duplicate patient/visit'):
                prepare(root,config)

    def test_crosswalk_must_resolve_real_assay_ids(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,config=self.setup_data(root)
            table=pd.read_csv(folder/'mapping.csv')
            table['metab_ref']=[f'unknown-{i}' for i in range(len(table))]
            table.to_csv(folder/'mapping.csv',index=False)
            report=prepare(root,config,inspect_only=True)
            self.assertEqual(report['verified_paired_samples'],0)
            self.assertEqual(prepare(root,config)['status'],'insufficient_paired_subjects')
            self.assertFalse((root/'processed').exists())

    def test_unmeasured_assay_row_not_counted_as_a_pair(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d); folder,config=self.setup_data(root)
            t=pd.read_csv(folder/'metab.csv')
            t.loc[0,['metabA','metabB']]=float('nan')
            t.to_csv(folder/'metab.csv',index=False)
            report=prepare(root,config,inspect_only=True)
            self.assertEqual(report['verified_paired_samples'],23)
            self.assertEqual(report['listed_pairs_without_both_measurements'],1)

    def test_nonnumeric_measurements_not_silently_imputed(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,config=self.setup_data(root)
            table=pd.read_csv(folder/'protein.csv').astype({'protA':str});table.loc[0,'protA']='text'
            table.to_csv(folder/'protein.csv',index=False)
            with self.assertRaisesRegex(ValueError,'nonnumeric'):
                prepare(root,config)

    def test_columns_oriented_assay(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);folder,config=self.setup_data(root)
            t=pd.read_csv(folder/'protein.csv')
            trans=t.set_index('specimen').T.reset_index().rename(columns={'index':'feature'})
            trans.to_csv(folder/'protein2.csv',index=False)
            config['protein']={'path':'raw/csiro/protein2.csv','sample_id_column':'feature',
                               'orientation':'samples_as_columns',
                               'sample_columns':t.specimen.tolist()}
            report=prepare(root,config)
            self.assertEqual(report['status'],'canonical_raw_omics_exported')
            p=pd.read_csv(root/'processed/external_validation/csiro/proteomics.csv')
            self.assertEqual(p.shape,(24,3))


if __name__=='__main__':unittest.main()