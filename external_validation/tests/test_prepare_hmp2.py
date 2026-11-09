"""Focused offline checks for real HMP2 source parsing and patient linkage."""
import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from external_validation.prepare_hmp2 import metadata_map, matched_visits, prepare


class Hmp2PreparationTests(unittest.TestCase):
    def write_fixture(self, base):
        raw=base/'raw'/'hmp2'
        raw.mkdir(parents=True)
        ids=[f'S{i:02d}_{visit}' for i in range(12) for visit in (0,2)]
        pd.DataFrame([
            {'External ID':sid,'Participant ID':sid.split('_')[0], 'week_num':sid.split('_')[1]}
            for sid in ids
        ]).to_csv(raw/'hmp2_metadata_2018-08-20.csv',index=False)
        with gzip.open(raw/'HMP2_proteomics_1pep-pro_1p_FDR_5ppm.tsv.gz','wt') as fh:
            w=csv.writer(fh,delimiter='\t');w.writerow(['feature',*ids])
            for feat in range(5):
                w.writerow([f'p{feat}',*[1+feat+i*0.1 for i in range(len(ids))]])
        obj={
          'format':'Biological Observation Matrix 1.0.0','matrix_type':'sparse',
          'shape':[3,len(ids)],
          'rows':[{'id':f'm{i}','metadata':None} for i in range(3)],
          'columns':[{'id':sid,'metadata':None} for sid in ids],
          'data':[[i,j, float(i+j+1)] for i in range(3) for j in range(len(ids))]
        }
        with gzip.open(raw/'HMP2_metabolomics_w_metadata.biom.gz','wt') as fh:json.dump(obj,fh)
        return raw,ids

    def test_complete_pairing_and_nonleaking_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);raw,ids=self.write_fixture(root)
            inspection=prepare(root,inspect_only=True)
            self.assertEqual(inspection['join']['exact_protein_metabolite_shared_ids'],24)
            self.assertEqual(inspection['join']['subjects_with_two_or_more_paired_visits'],12)
            self.assertFalse((root/'processed').exists())
            outcome=prepare(root)
            self.assertEqual(outcome['status'],'canonical_raw_omics_exported')
            output=root/'processed/external_validation/hmp2'
            s=pd.read_csv(output/'subjects.csv',dtype={'sample_id':str,'patient_id':str})
            p=pd.read_csv(output/'proteomics.csv');m=pd.read_csv(output/'metabolomics.csv')
            self.assertEqual(len(s),24)
            self.assertEqual(p.shape,(24,6))
            self.assertEqual(m.shape,(24,4))
            self.assertEqual(s.visit_order.loc[0],0)
            self.assertEqual(m.loc[0,'m0'],1)
            self.assertTrue((raw/'HMP2_metabolomics_w_metadata.biom.gz').exists())

    def test_metadata_identical_duplicates_are_allowed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'metadata.csv'
            pd.DataFrame([
                {'External ID':'s1','Participant ID':'P1','week_num':1},
                {'External ID':'s1','Participant ID':'P1','week_num':1}
            ]).to_csv(p,index=False)
            m,info=metadata_map(p)
            self.assertEqual(len(m),1)
            self.assertEqual(info['metadata_duplicate_rows_consistent'],1)

    def test_metadata_conflicting_patient_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'metadata.csv'
            pd.DataFrame([
                {'External ID':'s1','Participant ID':'P1','week_num':1},
                {'External ID':'s1','Participant ID':'P2','week_num':1}
            ]).to_csv(p,index=False)
            with self.assertRaisesRegex(ValueError,'conflict'):
                metadata_map(p)

    def test_duplicate_patient_week_excluded_instead_of_arbitrary_join(self):
        m=pd.DataFrame({'Participant ID':['P1','P1','P1'], 'week_num':[1,1,2]},index=['S1','S2','S3'])
        patients, info = matched_visits(m,['S1','S2','S3'],['S1','S2','S3'])
        self.assertEqual(info['rows_excluded_due_to_duplicate_patient_week'],2)
        self.assertEqual(len(patients),1)

    def test_no_fake_match_on_similar_sample_ids(self):
        m=pd.DataFrame({'Participant ID':['P1','P1'],'week_num':[1,2]},index=['S1','S2'])
        patients,info=matched_visits(m,['S1'],['S2'])
        self.assertEqual(len(patients),0)
        self.assertEqual(info['exact_protein_metabolite_shared_ids'],0)


if __name__=='__main__':unittest.main()