"""Check cohort readiness without running PULSE training or touching original data."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from external_validation.preflight import preflight


class PreflightTests(unittest.TestCase):
    def test_reports_true_patient_held_out_sets(self):
        with tempfile.TemporaryDirectory() as root:
            folder=Path(root)/'processed/external_validation/hmp2'
            folder.mkdir(parents=True)
            subjects=[];protein=[];metab=[]
            for i in range(14):
                for v in (1,2):
                    sid=f'p{i}-v{v}'
                    subjects.append({'sample_id':sid,'patient_id':f'p{i}','visit_order':v})
                    protein.append({'sample_id':sid,'p1':float(i+v),'p2':float(i*i+v)})
                    metab.append({'sample_id':sid,'m1':float(i+2*v),'m2':float(i*i+2*v)})
            pd.DataFrame(subjects).to_csv(folder/'subjects.csv',index=False)
            pd.DataFrame(protein).to_csv(folder/'proteomics.csv',index=False)
            pd.DataFrame(metab).to_csv(folder/'metabolomics.csv',index=False)
            report=preflight(Path(root),'hmp2')
            self.assertEqual(report['patients'],14)
            self.assertEqual(report['visits'],28)
            h=report['directions']['protein_to_metab']['historical_completion']
            self.assertEqual(h['eligible_visits'],14)
            self.assertEqual(h['eligible_patients'],14)
            self.assertEqual(sum(h['split_patients'].values()),14)
            self.assertEqual(h['status'],'ready_for_training')

if __name__=='__main__':unittest.main()