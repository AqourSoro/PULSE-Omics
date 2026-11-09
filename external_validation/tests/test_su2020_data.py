"""Local-only validation of the Su 2020 Table S1 acquisition code."""
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from openpyxl import Workbook

from external_validation.scripts.su2020_data import install_workbook, obtain_workbook, validate_workbook


class Su2020DataTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.root = Path(self.folder.name)
        wb = Workbook()
        wb.active.title = "Sample metadata"
        wb.active.append(["sample_id", "subject_id", "visit"])
        wb.active.append(["s1", "p1", 1])
        self.workbook = self.root / "mmc1.xlsx"
        wb.save(self.workbook)

    def tearDown(self):
        self.folder.cleanup()

    def test_direct_xlsx_import(self):
        output = self.root / "raw" / "table_s1.xlsx"
        install_workbook(self.workbook, output, "test data", overwrite=False)
        validate_workbook(output)
        self.assertTrue((output.parent / "download_manifest.json").is_file())

    def test_supplementary_zip_import(self):
        archive = self.root / "supplements.zip"
        with ZipFile(archive, "w") as out:
            out.write(self.workbook, "supplement/mmc1.xlsx")
            out.writestr("supplement/mmc2.xlsx", b"ignore")
        result = self.root / "extracted.xlsx"
        obtain_workbook(archive, result)
        validate_workbook(result)

    def test_html_page_rejected(self):
        html = self.root / "response.html"
        html.write_text("<html>Preparing to download</html>", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Neither XLSX"):
            obtain_workbook(html, self.root / "invalid.xlsx")


if __name__ == "__main__":
    unittest.main()
