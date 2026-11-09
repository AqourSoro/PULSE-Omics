"""Source-specific CSIRO canonical preparation tests on synthetic fixtures."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from external_validation.prepare_csiro import prepare


def fixture(root: Path):
    base=root/'raw/csiro';base.mkdir(parents=True)
    rows=[];meta=[]
    for j in range(16):
        sid=f'B{j+1:03d}'
        rows.append({'Sample IDs':sid,'Treatments':'vax','proteomic#':str(j+1),'Vac':'x'})
        meta.append({'Data File':f'S{j+1}', 'Subject ID':f'patient{j//2}',
                     'Sample IDS (Barcode)':sid,'Times':'T'+str(j%2+1),
                     'IDs & Treatments':sid+' vax', 'COR':'none','Year of Birth':'1990',
                     'Vaccine Received':'AZ', 'Lactate':[2.3,3.2][j%2]})
    pd.DataFrame(rows).to_excel(base/'005_Task mapping.xlsx',sheet_name='Sheet1',index=False)
    pd.DataFrame(meta).to_excel(base/'006_Task metabolomics.xlsx',sheet_name='COVID HPV CCM',index=False)
    for i,(prefix,ids) in enumerate([('002_Task',range(1,7)),('003_Task',range(7,12)),('004_Task',range(12,17))]):
        record={'Accession':['ProteinA','ProteinB']}
        for n in ids:record[f'Abundance: Vac_Plasma_{n:03d}']=[n*3.1,n*0.3]
        pd.DataFrame(record).to_excel(base/f'{prefix} result.xlsx',sheet_name='Proteins',index=False)


class CsiroPreparerTests(unittest.TestCase):
    def test_16_matched_visits_and_8_subjects(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixture(root)
            pre=prepare(root, inspect_only=True)
            self.assertEqual((pre['verified_paired_samples'],pre['paired_subjects'],pre['eligible_historical']),(16,8,8))
            self.assertFalse((root/'processed').exists())
            report=prepare(root)
            self.assertEqual(report['status'],'canonical_raw_omics_exported')
            out=root/'processed/external_validation/csiro'
            prot=pd.read_csv(out/'proteomics.csv')
            metab=pd.read_csv(out/'metabolomics.csv')
            subs=pd.read_csv(out/'subjects.csv')
            self.assertEqual(len(subs),16)
            self.assertIn('protein:ProteinA',prot)
            self.assertIn('Lactate',metab)
            self.assertEqual(subs.groupby('patient_id').size().max(),2)
            self.assertEqual(prot.loc[0,'sample_id'],metab.loc[0,'sample_id'])

    def test_no_implicit_mapping_by_row_number(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);fixture(root)
            path=root/'raw/csiro/005_Task mapping.xlsx'
            tab=pd.read_excel(path)
            tab['Sample IDs']=['UNKNOWN'+str(i) for i in range(len(tab))]
            tab.to_excel(path,index=False,sheet_name='Sheet1')
            report=prepare(root,inspect_only=True)
            self.assertEqual(report['status'],'needs_verified_crosswalk')


if __name__=='__main__':unittest.main()