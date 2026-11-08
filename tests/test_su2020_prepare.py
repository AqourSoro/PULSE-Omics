"""Su2020 extraction checks using a tiny workbook with known IDs."""
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from openpyxl import Workbook

from omics.su2020 import prepare


def fixture(path: Path, bridge: bool, direct_labels: bool = False) -> None:
    wb = Workbook()
    clinical = wb.active
    clinical.title = "S1.1 Patient Clinical Data"
    headings = ["Study Subject ID", "Sample ID", "Blood draw time point"]
    if bridge:
        headings += ["Plasma protein sample", "Plasma metab sample"]
    clinical.append(headings)
    proteins = wb.create_sheet("S1.3 Plasma Proteomic Data")
    proteins.append(["Healthy donor sample or COVID19 sample", "sample_id", "protein A", "protein B"])
    metabolites = wb.create_sheet("S1.4 Plasma Metabolomic Data")
    metabolites.append(["Healthy donor sample or COVID19 sample", "sample_id", "metabolite A", "metabolite B"])
    for i in range(25):
        for visit in (1, 2):
            p = f"p{i:03d}"
            blood = f"{p}-{visit}"
            protein_id = f"{10000 + 2 * i + visit}"
            metab_id = f"{20000 + 2 * i + visit}"
            row = [p, blood, f"T{visit}"]
            if bridge:
                row += [protein_id, metab_id]
            clinical.append(row)
            label = blood if direct_labels else "COVID19 sample"
            proteins.append([label, protein_id, float(i), float(visit)])
            metabolites.append([label.lower() if direct_labels else label, metab_id, float(i + visit), None if i == 0 else 3.0])
    for i in range(4):
        proteins.append(["Healthy Donor", str(30000 + i), 10.0, 11.0])
        metabolites.append(["Healthy Donor", str(30000 + i), 3.0, 4.0])
    wb.save(path)


class Su2020PreparationTests(unittest.TestCase):
    def test_clinical_assay_bridge_and_raw_nan(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root / "table_s1.xlsx"
            fixture(path, bridge=True)
            report = prepare(path, root / "processed")
            self.assertEqual(report["clinical_patients"], 25)
            self.assertEqual(report["clinical_patients_with_two_draws"], 25)
            self.assertEqual(report["assay_shared_sample_ids"], 4)
            self.assertEqual(report["linked_clinical_draws"], 50)
            self.assertEqual(report["linked_paired_draws"], 50)
            self.assertEqual(report["patients_with_two_paired_omics_visits"], 25)
            self.assertEqual(report["status"], "canonical_raw_omics_exported")
            subjects = pd.read_csv(root / "processed" / "subjects.csv")
            protein = pd.read_csv(root / "processed" / "proteomics.csv")
            metab = pd.read_csv(root / "processed" / "metabolomics.csv")
            self.assertEqual(len(subjects), 50)
            self.assertEqual(len(protein), 50)
            self.assertEqual(len(metab), 50)
            self.assertEqual(set(subjects.visit_order), {1, 2})
            self.assertEqual(int(metab["metabolite B"].isna().sum()), 2)

    def test_explicit_covid_labels_connect_without_numeric_bridge(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root / "table_s1.xlsx"
            fixture(path, bridge=False, direct_labels=True)
            report = prepare(path, root / "processed")
            self.assertEqual(report["exact_clinical_label_matches"], {"protein": 50, "metab": 50})
            self.assertEqual(report["linked_clinical_draws"], 50)
            self.assertEqual(report["linked_paired_draws"], 50)
            self.assertEqual(report["linked_patients_with_two_draws"], 25)
            mapping = pd.read_csv(root / "processed" / "assay_mapping.csv")
            self.assertEqual(len(mapping), 50)
            self.assertEqual(set(mapping["Study Subject ID"].unique()), {f"p{i:03d}" for i in range(25)})

    def test_incov_clinical_and_bl_ac_omics_ids(self):
        """Scientific source uses INCOV001-1/2 and INCOV001-BL/AC."""
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root / "su2020.xlsx"
            wb = Workbook()
            clinical = wb.active
            clinical.title = "S1.1 Patient Clinical Data"
            clinical.append(["Study Subject ID", "Sample ID", "Blood draw time point"])
            prot = wb.create_sheet("S1.3 Plasma Proteomic Data")
            metab = wb.create_sheet("S1.4 Plasma Metabolomic Data")
            for sheet, feature in ((prot, "protein A"), (metab, "metabolite A")):
                sheet.append(["Healthy donor sample or COVID19 sample", "sample_id", feature])
                sheet.append(["Healthy Donor", 1065558, 5.0])
            for p in range(1, 5):
                for visit in (1, 2):
                    clinical.append([f"INCOV{p:03d}", f"INCOV{p:03d}-{visit}", f"T{visit}"])
                    suffix = "BL" if visit == 1 else "AC"
                    prot.append(["COVID19 sample", f"INCOV{p:03d}-{suffix}", 3.0 * p])
                    metab.append(["COVID19 sample", f"COVID_{p}_T{visit}", 4.0 * p])
            wb.save(path)
            report = prepare(path, root / "processed")
            self.assertEqual(report["clinical_sample_matches"], {"protein": 8, "metab": 8})
            self.assertEqual(report["exact_clinical_label_matches"], {"protein": 0, "metab": 0})
            self.assertEqual(report["linked_clinical_draws"], 8)
            self.assertEqual(report["linked_paired_draws"], 8)
            self.assertEqual(report["linked_patients_with_two_draws"], 4)
            mapping = pd.read_csv(root / "processed" / "assay_mapping.csv")
            self.assertEqual(set(mapping["protein_assay_id"]),
                             {f"INCOV{p:03d}-{s}" for p in range(1, 5) for s in ("BL", "AC")})

    def test_numeric_donor_ids_are_not_converted_into_patient_ids(self):
        self.assertIsNone(
            __import__("omics.su2020", fromlist=["normalized_patient_sample"]).normalized_patient_sample(1065558)
        )

    def test_no_bridge_never_invents_patient_links(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root / "table_s1.xlsx"
            fixture(path, bridge=False)
            report = prepare(path, root / "processed")
            self.assertEqual(report["status"], "raw_assays_exported_no_verified_patient_join")
            self.assertFalse((root / "processed" / "subjects.csv").exists())
            self.assertTrue((root / "processed" / "protein_raw_by_assay_id.csv").exists())
            self.assertTrue((root / "processed" / "metab_raw_by_assay_id.csv").exists())


if __name__ == "__main__":
    unittest.main()
