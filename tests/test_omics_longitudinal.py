"""Focused checks for sample identity and modality missingness."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from omics.longitudinal import load_sequences


class LongitudinalOmicsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        pd.DataFrame([
            {"sample_id": "s2", "patient_id": "p1", "visit_order": 2},
            {"sample_id": "s3", "patient_id": "p2", "visit_order": 1},
            {"sample_id": "s1", "patient_id": "p1", "visit_order": 1},
        ]).to_csv(self.root / "subjects.csv", index=False)
        pd.DataFrame([
            {"sample_id": "s1", "pA": 1.0, "pB": 2.0},
            {"sample_id": "s2", "pA": 3.0, "pB": 4.0},
            {"sample_id": "s3", "pA": 5.0, "pB": 6.0},
        ]).to_csv(self.root / "proteomics.csv", index=False)
        pd.DataFrame([
            {"sample_id": "s1", "mA": 7.0},
            {"sample_id": "s3", "mA": 8.0},
        ]).to_csv(self.root / "metabolomics.csv", index=False)

    def tearDown(self):
        self.tmp.cleanup()

    def load(self):
        return load_sequences(
            self.root / "subjects.csv",
            self.root / "proteomics.csv",
            self.root / "metabolomics.csv",
        )

    def test_patient_visit_order_and_absent_whole_modality(self):
        ids, visits, masks, features = self.load()
        self.assertEqual(ids, ["p1", "p2"])
        self.assertEqual([len(v) for v in visits], [2, 1])
        self.assertEqual(visits[0][0]["protein"].tolist(), [1.0, 2.0])
        self.assertEqual(visits[0][1]["protein"].tolist(), [3.0, 4.0])
        self.assertIsNone(visits[0][1]["metab"])
        self.assertEqual(masks[0][1], {"protein": 0, "metab": 2})
        self.assertEqual(features["protein"], ["pA", "pB"])

    def test_partial_missingness_is_not_marked_observed(self):
        path = self.root / "proteomics.csv"
        data = pd.read_csv(path)
        data.loc[0, "pB"] = float("nan")
        data.to_csv(path, index=False)
        with self.assertRaisesRegex(ValueError, "partially missing"):
            self.load()


if __name__ == "__main__":
    unittest.main()
