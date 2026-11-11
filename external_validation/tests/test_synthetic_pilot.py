"""Data-free unit test for the independent patient split and generation contract."""
import unittest

import numpy as np
import pandas as pd

from external_validation.synthetic_pilot import FactorCopula, correlation_gap, split_patients


class FastPilotTest(unittest.TestCase):
    def test_synthetic_patient_pairs_and_covariance(self):
        rng = np.random.default_rng(3)
        patients = 50
        samples = patients * 2
        meta = pd.DataFrame({
            "patient_id": np.repeat([f"p{i}" for i in range(patients)], 2),
            "visit_order": np.tile([1, 2], patients),
            "sample_id": [f"s{i}" for i in range(samples)],
        })
        f = np.repeat(rng.standard_normal((patients, 3)), 2, axis=0)
        f += rng.standard_normal((samples, 3)) * .25
        x = f @ rng.standard_normal((3, 12)) + rng.standard_normal((samples, 12)) * .2
        split = split_patients(meta, 42)
        self.assertFalse(split["train"] & split["val"])
        self.assertFalse(split["train"] & split["test"])
        model = FactorCopula(n_factors=3, seed=5).fit(x, meta)
        synthetic_metadata, synthetic = model.sample(70, rng)
        self.assertEqual(synthetic.shape, (140, 12))
        self.assertEqual(synthetic_metadata.patient_id.nunique(), 70)
        self.assertTrue(np.isfinite(synthetic).all())
        self.assertTrue(np.isfinite(list(correlation_gap(x, synthetic, 6).values())).all())
        for pid, records in synthetic_metadata.groupby("patient_id"):
            self.assertEqual(records.visit_order.tolist(), [1, 2], pid)


if __name__ == "__main__":
    unittest.main()