"""Tiny fixture checks for actual header conventions in COVAAR and CSIRO."""
import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from external_validation.cohort_join_audit import audit_covaar, audit_csiro, _vac_id


class CohortJoinAuditTests(unittest.TestCase):
    def test_covaar_exact_metadata_vs_geo(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); base=root/'raw/covaar'
            (base/'figshare').mkdir(parents=True)
            (base/'geo').mkdir(parents=True)
            pd.DataFrame([{'Sample_ID':'Sample_001','Patient':'P1','Time':'T1'},
                          {'Sample_ID':'Sample_002','Patient':'P1','Time':'T2'}]).to_csv(
                              base/'figshare/Sample_Meta_Data.csv',index=False)
            pd.DataFrame({'Metabolite':['X'],'Sample_001':[1.0],'Sample_002':[2.0]}).to_csv(
                base/'figshare/Metabolomics_Data.csv',index=False)
            with gzip.open(base/'geo/GSE325178_SomaScan_RFUs_allSamples.txt.gz','wt') as f:
                f.write('AptName\tSeqId\tTarget\tEntrezGeneSymbol\tUniProt\tVA0001\tVA0002\n')
                f.write('X\tY\tZ\tQ\tR\t2\t3\n')
            with gzip.open(base/'geo/GSE325178_series_matrix.txt.gz','wt') as f:
                f.write('!Sample_title\t"VA0001"\t"VA0002"\n')
                f.write('!Sample_characteristics_ch1\t"patient: P1"\t"patient: P1"\n')
            report=audit_covaar(root)
            self.assertEqual(report['metabolite_to_metadata']['matching_distinct'],2)
            self.assertEqual(report['protein_to_geo_title']['matching_distinct'],2)
            self.assertEqual(report['protein_to_metabolite']['matching_distinct'],0)
            self.assertEqual(report['metadata_patients'],1)

    def test_covaar_sample_ids_can_live_in_geo_characteristics(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=root/'raw/covaar'
            (base/'figshare').mkdir(parents=True)
            (base/'geo').mkdir(parents=True)
            pd.DataFrame([
                {'Sample_ID':'Sample_001','Patient':'P1','Time':'T1'},
                {'Sample_ID':'Sample_002','Patient':'P1','Time':'T2'},
            ]).to_csv(base/'figshare/Sample_Meta_Data.csv',index=False)
            pd.DataFrame({'Metabolite':['X'],'Sample_001':[1], 'Sample_002':[2]}).to_csv(
                base/'figshare/Metabolomics_Data.csv',index=False)
            with gzip.open(base/'geo/GSE325178_SomaScan_RFUs_allSamples.txt.gz','wt') as stream:
                stream.write('AptName\\tSeqId\\tTarget\\tEntrezGeneSymbol\\tUniProt\\tVA0001\\tVA0002\\n')
                stream.write('X\\tY\\tZ\\tQ\\tR\\t2\\t3\\n')
            with gzip.open(base/'geo/GSE325178_series_matrix.txt.gz','wt') as stream:
                stream.write('!Sample_title\\tGEO-title-1\\tGEO-title-2\\n')
                stream.write('!Sample_characteristics_ch1\\t"sample id: VA0001"\\t"sample id: VA0002"\\n')
                stream.write('!Sample_characteristics_ch1\\t"patient: P1"\\t"patient: P1"\\n')
                stream.write('!Sample_characteristics_ch1\\t"time: T1"\\t"time: T2"\\n')
            result=audit_covaar(root)
            self.assertEqual(result['protein_to_geo_title']['matching_distinct'],0)
            self.assertEqual(result['protein_ids_found_in_any_geo_sample_field'],2)
            self.assertEqual(result['patient_time_tuples_matching_figshare'],2)

    def test_csiro_published_crosswalk_headers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);base=root/'raw/csiro';base.mkdir(parents=True)
            pd.DataFrame({'Sample IDs':['BC01','BC02'], 'Treatments':['A','A'],
                          'proteomic#':['1','2'],'Vac':['AZ','AZ']}).to_excel(
                            base/'005_Task mapping.xlsx',index=False,sheet_name='Sheet1')
            pd.DataFrame({'Subject ID':['Patient1','Patient1'],
                          'Sample IDS (Barcode)':['BC01','BC02'],'Times':['T1','T2'],
                          'IDs & Treatments':['BC01 A','BC02 A'], 'Met1':[1.,2.]}).to_excel(
                            base/'006_Task metabolomics.xlsx',index=False,sheet_name='COVID HPV CCM')
            for prefix,code,cols in [('002_Task','AZ',['Abundance: Vac_Plasma_001','Abundance: Vac_Plasma_002']),
                                ('003_Task','Moderna',['Abundance: Vac_Plasma_015']),
                                ('004_Task','Pfizer',['Abundance: Vac_Plasma_063'])]:
                pd.DataFrame({'Accession':['P001'],cols[0]:[1.],**({cols[1]:[2.]} if len(cols)>1 else {})}).to_excel(
                    base/(prefix+' result.xlsx'),index=False,sheet_name='Proteins')
            report=audit_csiro(root)
            self.assertEqual(report['mapped_proteomics_to_actual_columns']['matching_distinct'],2)
            self.assertEqual(report['crosswalk_sample_id_to_metabolomics']['Sample IDS (Barcode)']['matching_distinct'],2)
            self.assertEqual(report['status'],'candidate_direct_crosswalk')

    def test_vac_id_rejects_ambiguous_values(self):
        self.assertEqual(_vac_id('Abundance: Vac_Plasma_001'), 'Vac_Plasma_001')
        self.assertEqual(_vac_id('1'), 'Vac_Plasma_001')
        self.assertIsNone(_vac_id('Patient1 T1'))


if __name__=='__main__':unittest.main()