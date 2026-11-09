"""Offline invariants for source-based iPOP pairing."""
import tempfile
import unittest
from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from external_validation.prepare_ipop import prepare, _date


def records(n=12):
    labels=[f'p{i}' for i in range(n) for _ in range(2)]
    dates=[datetime(2015,1,15) if j==0 else datetime(2015,3,15)
           for _ in range(n) for j in range(2)]
    annotation=pd.DataFrame({
        'SubjectID': labels,
        'Date':[(d-datetime(1970,1,1)).total_seconds() for d in dates],
        'Time':[d.timetuple().tm_yday for d in dates],
    },index=[f'annotations_{i}' for i in range(n*2)])
    protein=pd.DataFrame({'p1':np.arange(n*2)*1.,'p2':np.arange(n*2)**2+1.},
                         index=[f'protein_{i}' for i in range(n*2)])
    metab=pd.DataFrame({'m1':np.arange(n*2)**2+1.,'m2':np.arange(n*2)*2.+3.},
                       index=[f'metabolite_{i}' for i in range(n*2)])
    return {'Proteomics_annotation_colData':annotation,
            'Metabolomics_annotation_colData':annotation.copy(),
            'Proteomics_df_Data':protein,
            'Metabolomics_df_Data':metab}


class IpopPreparationTests(unittest.TestCase):
    def test_release_row_order_requires_explicit_provenance(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            f=root/'raw/ipop/seasonal/Multi_Omics_Seasonal.RData'
            f.parent.mkdir(parents=True);f.write_bytes(b'fixture')
            with patch('external_validation.prepare_ipop.read_ipop',return_value=records()):
                with self.assertRaisesRegex(ValueError,'row labels differ'):
                    prepare(root,inspect_only=True)
                report=prepare(root,inspect_only=True,permit_release_row_order=True)
                self.assertEqual(report['paired_samples'],24)
                self.assertEqual(report['eligible_historical'],12)
                self.assertFalse((root/'processed').exists())
                prepare(root,permit_release_row_order=True)
            out=root/'processed/external_validation/ipop'
            a=pd.read_csv(out/'subjects.csv')
            p=pd.read_csv(out/'proteomics.csv')
            m=pd.read_csv(out/'metabolomics.csv')
            self.assertEqual(a.shape,(24,4))
            self.assertEqual(p.shape,(24,3))
            self.assertEqual(m.shape,(24,3))
            self.assertEqual(set(a.visit_order),{0.,59.})
            self.assertTrue(a.sample_id.is_unique)

    def test_two_digit_r_date_requires_documented_offset(self):
        seconds=(datetime(15,3,1)-datetime(1970,1,1)).total_seconds()
        d,repaired=_date(seconds)
        self.assertEqual(d.year,2015)
        self.assertTrue(repaired)
        with self.assertRaisesRegex(ValueError,'implausible year'):
            _date(0)

if __name__=='__main__': unittest.main()
