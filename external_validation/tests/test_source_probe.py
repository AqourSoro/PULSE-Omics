"""Offline checks for source parsing diagnostics without real participant data."""
import unittest
from pathlib import Path
from unittest.mock import patch

from external_validation.source_probe import _ipop_rdata_structure


class SourceProbeTests(unittest.TestCase):
    def test_missing_rdata_dependency_is_reported_as_setup_issue(self):
        error = ModuleNotFoundError("No module named 'rdata'", name="rdata")
        with patch(
            "external_validation.scripts.inspect_restored_omics.read_ipop",
            side_effect=error,
        ):
            result = _ipop_rdata_structure(Path("/not/read.rdata"))
        self.assertEqual(result["status"], "dependency_missing")
        self.assertIn("requirements-structure.txt", result["install"])

    def test_real_decoder_error_is_distinguished_from_missing_dependency(self):
        with patch(
            "external_validation.scripts.inspect_restored_omics.read_ipop",
            side_effect=ValueError("R serialization was not decoded"),
        ):
            result = _ipop_rdata_structure(Path("/not/read.rdata"))
        self.assertEqual(result["status"], "cannot_read_rdata")
        self.assertIn("R serialization", result["error"])


if __name__ == "__main__":
    unittest.main()
