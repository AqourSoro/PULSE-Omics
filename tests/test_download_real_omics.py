"""Offline, focused tests for real-data download safety and resumable transfers."""
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.download_real_omics import Downloader, check_format, safe_path


class FakeResponse:
    def __init__(self, payload: bytes, status=200):
        self.payload = io.BytesIO(payload)
        self.status = status
        self.headers = {"Content-Length": str(len(payload))}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.payload.close()

    def read(self, size):
        return self.payload.read(size)


class DownloaderTests(unittest.TestCase):
    def test_rejects_path_traversal(self):
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            safe_path(Path("/tmp/omics"), "../../private.txt")

    def test_html_is_not_a_zip(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "supplement.zip"
            path.write_bytes(b"<html>please sign in</html>")
            with self.assertRaises(ValueError):
                check_format(path)

    def test_https_transfer_is_recorded_and_idempotent(self):
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), False, 1.0, 1.0)
            payload = b"PK\x03\x04" + b"test-supplement-payload"
            with patch(
                "scripts.download_real_omics.urlopen",
                return_value=FakeResponse(payload),
            ) as request:
                runner.https("p100", "test.zip", "https://example.org/test.zip")
                runner.https("p100", "test.zip", "https://example.org/test.zip")
            runner.flush("p100")
            self.assertEqual(request.call_count, 1)
            self.assertEqual(runner.downloaded, 1)
            self.assertEqual(runner.existing, 1)
            self.assertTrue((Path(root) / "raw/p100/download_manifest.json").exists())

    def test_ipop_uses_published_processed_source_when_s3_is_unavailable(self):
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), False, 1.0, 1.0)
            with patch.dict(sys.modules, {"boto3": None}), patch.object(runner, "https") as download:
                runner.ipop()
            download.assert_called_once_with(
                "ipop", "seasonal/Multi_Omics_Seasonal.RData",
                "https://ndownloader.figshare.com/files/22807673",
            )
            self.assertEqual(runner.failed, 0)

    def test_complete_hmp2_indexes_all_links(self):
        """The downloader must enumerate source links, not just hard-coded matrices."""
        html = b'<html><a href="https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MPX/first.biom.gz">first</a><a href="https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MPX/second.tsv.gz">second</a></html>'
        with tempfile.TemporaryDirectory() as root:
            with patch("scripts.download_real_omics.urlopen", return_value=FakeResponse(html)):
                from scripts.download_real_omics import list_http_files
                found = list_http_files(
                    "https://www.ibdmdb.org/downloads/html/products_HMP2_Proteomics.html",
                    {"g-227ca.190ebd.75bc.data.globus.org"},
                )
            self.assertEqual(len(found), 2)
            self.assertTrue(any(item.endswith("first.biom.gz") for item in found))

    def test_zero_download_caps_are_unlimited(self):
        runner = Downloader(Path("/tmp"), True, 0, 0)
        self.assertEqual(runner.max_file, float("inf"))
        self.assertEqual(runner.max_total, float("inf"))

    def test_http_range_resumes_incomplete_file(self):
        with tempfile.TemporaryDirectory() as root:
            dataset = Path(root) / "raw" / "covaar"
            dataset.mkdir(parents=True)
            (dataset / "sample.txt.gz.part").write_bytes(b"\x1f\x8b")
            runner = Downloader(Path(root), False, 1.0, 1.0)
            with patch(
                "scripts.download_real_omics.urlopen",
                return_value=FakeResponse(b"\x08rest-of-payload", status=206),
            ) as request:
                runner.https(
                    "covaar", "sample.txt.gz",
                    "https://example.org/sample.txt.gz",
                )
            req = request.call_args.args[0]
            self.assertEqual(req.headers["Range"], "bytes=2-")
            self.assertEqual(runner.downloaded, 1)
            self.assertFalse((dataset / "sample.txt.gz.part").exists())


if __name__ == "__main__":
    unittest.main()

