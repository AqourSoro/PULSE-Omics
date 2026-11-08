"""Check that the Su 2020 schema profiler measures genuine ID overlap only."""
import tempfile
import unittest
from pathlib import Path

from openpyxl import Workbook

from scripts.su2020_schema import build_report


class SchemaProfilerTests(unittest.TestCase):
    def test_sample_overlap_and_repeated_patients(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table_s1.xlsx"
            wb = Workbook()
            clinical = wb.active
            clinical.title = "S1.1 Patient Clinical Data"
            clinical.append(["sample_id", "patient_id", "visit_order"])
            protein = wb.create_sheet("S1.3 Plasma Proteomic Data")
            protein.append(["sample_id", "P1", "P2"])
            metab = wb.create_sheet("S1.4 Plasma Metabolomic Data")
            metab.append(["sample_id", "M1", "M2"])
            for i in range(36):
                sample_id = f"s{i:03d}"
                clinical.append([sample_id, f"p{i // 2:03d}", 1 + (i % 2)])
                protein.append([sample_id, float(i), float(i + 10)])
                if i < 32:
                    metab.append([sample_id, float(i * 2), float(i * 3)])
            wb.save(path)

            report = build_report(path)
            self.assertEqual(
                report["cross_sheet_overlap"]["protein_to_metab"][0]["common_ids"], 32
            )
            self.assertEqual(
                report["cross_sheet_overlap"]["clinical_to_protein"][0]["common_ids"], 36
            )
            self.assertTrue(
                any(x["ids_with_multiple_rows"] >= 18
                    for x in report["clinical_repeated_identifier_candidates"])
            )
            self.assertEqual(report["sheets"]["protein"]["dimensions"], [37, 3])
            self.assertIn("not verified", report["note"])


if __name__ == "__main__":
    unittest.main()
