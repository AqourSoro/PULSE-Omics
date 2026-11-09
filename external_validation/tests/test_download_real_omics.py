"""Offline, focused tests for real-data download safety and resumable transfers."""
import io
from contextlib import redirect_stdout
from http.client import IncompleteRead
import sys
import tempfile
import unittest
from pathlib import Path
from urllib.error import URLError
from unittest.mock import patch

from external_validation.scripts.download_real_omics import Downloader, check_format, safe_path


class FakeResponse:
    def __init__(self, payload: bytes, status=200, range_start=0):
        self.payload = io.BytesIO(payload)
        self.status = status
        self.headers = {"Content-Length": str(len(payload))}
        if status == 206:
            self.headers["Content-Range"] = (
                f"bytes {range_start}-{range_start+len(payload)-1}/*"
            )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.payload.close()

    def read(self, size=-1):
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
                "external_validation.scripts.download_real_omics.urlopen",
                return_value=FakeResponse(payload),
            ) as request:
                runner.https("p100", "test.zip", "https://example.org/test.zip")
                runner.https("p100", "test.zip", "https://example.org/test.zip")
            runner.flush("p100")
            self.assertEqual(request.call_count, 1)
            self.assertEqual(runner.downloaded, 1)
            self.assertEqual(runner.existing, 1)
            self.assertTrue((Path(root) / "raw/p100/download_manifest.json").exists())

    def test_csiro_link_list_accepts_direct_and_downloadfile_urls(self):
        from external_validation.scripts.download_real_omics import csiro_download_targets
        raw = (
            "https://data.csiro.au/dap/ws/v2/collections/77992/downloadfile?"
            "fileName=Proteomics%2Fintensities.csv&amp;format=raw\n"
            "https://data.csiro.au/dap/files/metadata.xlsx\n"
            "https://evilcsiro.au/dap/files/not-trusted.csv\n"
            "https://data.csiro.au/dap/collection/landing-page\n"
        )
        found = csiro_download_targets(raw)
        self.assertEqual(len(found), 2)
        self.assertEqual(set(name for _, name in found),
                         {"intensities.csv", "metadata.xlsx"})
        self.assertTrue(any("fileName=Proteomics%2F" in url
                            for url, _ in found))

    def test_ipop_downloads_only_the_published_seasonal_release(self):
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), False, 1.0, 1.0)
            with patch.dict(sys.modules, {"boto3": None}), patch.object(
                runner, "figshare_article"
            ) as figshare, patch.object(runner, "https") as download:
                runner.ipop()
            figshare.assert_called_once_with("ipop", 12376508, "seasonal")
            download.assert_called_once_with(
                "ipop", "seasonal/Multi_Omics_Seasonal.RData",
                "https://ndownloader.figshare.com/files/22807673",
            )
            self.assertEqual(runner.failed, 0)
            self.assertNotIn("ipop", runner.incomplete)

    def test_figshare_uses_publisher_file_size_and_md5(self):
        import json
        metadata = {"files": [{
            "name": "assays.xlsx", "download_url": "https://ndownloader.figshare.com/files/12345",
            "size": 12, "computed_md5": "abcdef",
        }]}
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), False, 0, 0)
            with patch(
                "external_validation.scripts.download_real_omics.urlopen",
                return_value=FakeResponse(json.dumps(metadata).encode()),
            ), patch.object(runner, "https") as download:
                runner.figshare_article("covaar", 31685023, "figshare")
            download.assert_called_once_with(
                "covaar", "figshare/assays.xlsx",
                "https://ndownloader.figshare.com/files/12345",
                expected_size=12, expected_md5="abcdef",
            )
            self.assertEqual(runner.inventory["covaar"][0]["files"][0]["bytes"], 12)

    def test_figshare_bad_existing_file_is_replaced_using_md5(self):
        import hashlib
        payload = b"PK\x03\x04" + b"correct-data"
        expected = hashlib.md5(payload).hexdigest()
        with tempfile.TemporaryDirectory() as root:
            folder = Path(root) / "raw" / "covaar"
            folder.mkdir(parents=True)
            (folder / "sample.xlsx").write_bytes(b"PK\x03\x04" + b"old-wrong-value")
            runner = Downloader(Path(root), False, 0, 0)
            with patch(
                "external_validation.scripts.download_real_omics.urlopen",
                return_value=FakeResponse(payload),
            ) as request:
                runner.https(
                    "covaar", "sample.xlsx", "https://example.org/sample.xlsx",
                    expected_size=len(payload), expected_md5=expected,
                )
            self.assertEqual(request.call_count, 1)
            self.assertEqual((folder / "sample.xlsx").read_bytes(), payload)
            self.assertEqual(runner.downloaded, 1)

    def test_figshare_fully_downloaded_part_is_finalized_without_http(self):
        import hashlib
        payload = b"PK\x03\x04" + b"complete-file"
        with tempfile.TemporaryDirectory() as root:
            folder = Path(root) / "raw" / "covaar"
            folder.mkdir(parents=True)
            (folder / "sample.xlsx.part").write_bytes(payload)
            runner = Downloader(Path(root), False, 0, 0)
            with patch("external_validation.scripts.download_real_omics.urlopen") as request:
                runner.https(
                    "covaar", "sample.xlsx", "https://example.org/sample.xlsx",
                    expected_size=len(payload), expected_md5=hashlib.md5(payload).hexdigest(),
                )
            request.assert_not_called()
            self.assertEqual((folder / "sample.xlsx").read_bytes(), payload)
            self.assertFalse((folder / "sample.xlsx.part").exists())

    def test_complete_hmp2_indexes_all_links(self):
        """The downloader must enumerate source links, not just hard-coded matrices."""
        html = b'<html><a href="https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MPX/first.biom.gz">first</a><a href="https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb/products/HMP2/MPX/second.tsv.gz">second</a></html>'
        with tempfile.TemporaryDirectory() as root:
            with patch("external_validation.scripts.download_real_omics.urlopen", return_value=FakeResponse(html)):
                from external_validation.scripts.download_real_omics import list_http_files
                found = list_http_files(
                    "https://www.ibdmdb.org/downloads/html/products_HMP2_Proteomics.html",
                    {"g-227ca.190ebd.75bc.data.globus.org"},
                )
            self.assertEqual(len(found), 2)
            self.assertTrue(any(item.endswith("first.biom.gz") for item in found))

    def test_hmp2_double_slash_catalog_path_is_relative(self):
        """Regresses crash for /ibdmdb//products/HMP2/MGX/... paths."""
        href = (
            "https://g-227ca.190ebd.75bc.data.globus.org/"
            "ibdmdb//products/HMP2/MGX/2018-05-04/contigs/"
            "CSM5FZ4M_contigs.fna.gz"
        )
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), True, 0, 0)
            with patch(
                "external_validation.scripts.download_real_omics.HMP2_PRODUCT_PAGES",
                ("products_MGX_2017-08-12.html",),
            ), patch(
                "external_validation.scripts.download_real_omics.list_http_files", return_value=[href],
            ), patch.object(runner, "https") as downloader:
                runner.hmp2_all_products()
            downloader.assert_any_call(
                "hmp2",
                "products/HMP2/MGX/2018-05-04/contigs/CSM5FZ4M_contigs.fna.gz",
                href,
            )
            self.assertNotIn("hmp2", runner.incomplete)

    def test_hmp2_rejects_traversal_even_after_url_decoding(self):
        href = (
            "https://g-227ca.190ebd.75bc.data.globus.org/"
            "ibdmdb//products/HMP2/%2e%2e/%2e%2e/escape.txt"
        )
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), True, 0, 0)
            with patch(
                "external_validation.scripts.download_real_omics.HMP2_PRODUCT_PAGES", ("one.html",),
            ), patch(
                "external_validation.scripts.download_real_omics.list_http_files", return_value=[href],
            ), patch.object(runner, "https") as downloader:
                runner.hmp2_all_products()
            # The independent metadata CSV is still considered, but the
            # unsafe source object must not reach the downloader.
            self.assertEqual(downloader.call_count, 1)
            self.assertEqual(downloader.call_args.args[1], "hmp2_metadata_2018-08-20.csv")
            self.assertIn("hmp2", runner.incomplete)

    def test_su2020_attempts_all_six_pmc_supplements(self):
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), True, 0, 0)
            with patch.object(runner, "https") as download:
                runner.su2020_full_release()
            self.assertEqual(download.call_count, 8)
            self.assertEqual(
                [args.args[1] for args in download.call_args_list[-6:]],
                [f"pmc/mmc{i}.xlsx" for i in range(1, 7)],
            )

    def test_zero_download_caps_are_unlimited(self):
        runner = Downloader(Path("/tmp"), True, 0, 0)
        self.assertEqual(runner.max_file, float("inf"))
        self.assertEqual(runner.max_total, float("inf"))

    def test_chunked_incomplete_read_resumes_and_keeps_exact_bytes(self):
        content = b"\x1f\x8b" + b"abcdefghijklm"
        class BrokenStream(FakeResponse):
            def __init__(self):
                super().__init__(content)

            def read(self, size):
                if self.payload.tell() < 7:
                    return self.payload.read(7)
                raise IncompleteRead(content[7:10], len(content) - 7)

        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), False, 0, 0)
            with patch("external_validation.scripts.download_real_omics.urlopen", side_effect=[
                BrokenStream(),
                FakeResponse(content[7:], status=206, range_start=7),
            ]) as urlopen, patch("external_validation.scripts.download_real_omics.time.sleep"):
                runner.https("hmp2", "sample.gz", "https://example.org/sample.gz")
            self.assertEqual(urlopen.call_count, 2)
            self.assertEqual(
                urlopen.call_args_list[1].args[0].headers["Range"], "bytes=7-",
            )
            self.assertEqual(
                (Path(root) / "raw/hmp2/sample.gz").read_bytes(), content,
            )
            self.assertEqual(runner.downloaded, 1)
            self.assertEqual(runner.failed, 0)

    def test_range_ignored_restarts_instead_of_corrupting(self):
        content = b"\x1f\x8b" + b"entire-content"
        with tempfile.TemporaryDirectory() as root:
            folder = Path(root) / "raw" / "hmp2"
            folder.mkdir(parents=True)
            (folder / "sample.gz.part").write_bytes(content[:5])
            runner = Downloader(Path(root), False, 0, 0)
            with patch("external_validation.scripts.download_real_omics.urlopen",
                       return_value=FakeResponse(content)):
                runner.https("hmp2", "sample.gz", "https://example.org/sample.gz")
            self.assertEqual((folder / "sample.gz").read_bytes(), content)
            self.assertEqual(runner.failed, 0)

    def test_persistent_error_does_not_abort_subsequent_file(self):
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), False, 0, 0)
            with patch("external_validation.scripts.download_real_omics.urlopen",
                       side_effect=URLError("dropped")), \
                 patch("external_validation.scripts.download_real_omics.time.sleep"):
                runner.https("hmp2", "first.tsv", "https://example.org/first.tsv")
            self.assertEqual(runner.failed, 1)
            self.assertEqual(len(runner.report["hmp2"]), 1)
            with patch("external_validation.scripts.download_real_omics.urlopen",
                       return_value=FakeResponse(b"sample\tfeature\n1\t2\n")):
                runner.https("hmp2", "second.tsv", "https://example.org/second.tsv")
            self.assertEqual(runner.downloaded, 1)
            self.assertTrue((Path(root) / "raw/hmp2/second.tsv").exists())

    def test_noninteractive_download_status_only_prints_summary(self):
        output = io.StringIO()
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), True, 0, 0)
            runner.planned["hmp2"] = 200
            with redirect_stdout(output):
                for i in range(200):
                    runner.record("hmp2", {
                        "status": "candidate",
                        "file": f"private/individual_sample_{i}.tsv",
                    })
                runner.finish_dataset("hmp2")
        self.assertEqual(output.getvalue().count("\n"), 1)
        self.assertIn("checked 200/200", output.getvalue())
        self.assertNotIn("individual_sample_", output.getvalue())

    def test_interactive_progress_rewrites_one_line(self):
        class TTYBuffer(io.StringIO):
            def isatty(self):
                return True

        output = TTYBuffer()
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), True, 0, 0)
            runner.planned["hmp2"] = 2
            with redirect_stdout(output):
                runner.record("hmp2", {"status": "candidate", "file": "one.biom"})
                runner.record("hmp2", {"status": "candidate", "file": "two.biom"})
                runner.progress("hmp2", force=True)
                runner.finish_dataset("hmp2")
        self.assertEqual(output.getvalue().count("\n"), 1)
        self.assertIn("\r", output.getvalue())
        self.assertIn("checked 2/2", output.getvalue())
        self.assertNotIn("one.biom", output.getvalue())

    def test_eta_uses_recent_file_completion_rate(self):
        class TTYBuffer(io.StringIO):
            def isatty(self):
                return True

        output = TTYBuffer()
        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), True, 0, 0)
            runner.planned["hmp2"] = 20
            with patch("external_validation.scripts.download_real_omics.time.monotonic",
                       side_effect=[1.0, 11.0]), redirect_stdout(output):
                runner.record("hmp2", {"status": "candidate", "file": "a.tsv"})
                runner.record("hmp2", {"status": "candidate", "file": "b.tsv"})
                runner.finish_dataset("hmp2")
        self.assertIn("ETA ~3m", output.getvalue())
        self.assertIn("0.1 files/s", output.getvalue())
        self.assertEqual(output.getvalue().count("\n"), 1)

    def test_hmp2_uses_bounded_parallel_workers(self):
        from threading import Lock
        from time import sleep
        urls = [
            f"https://g-227ca.190ebd.75bc.data.globus.org/"
            f"ibdmdb/products/HMP2/MPX/data_{i}.tsv"
            for i in range(8)
        ]
        active = 0
        peak = 0
        calls = []
        lock = Lock()

        def fake_download(dataset, relative, url):
            nonlocal active, peak
            with lock:
                calls.append(relative)
                active += 1
                peak = max(peak, active)
            sleep(0.02)
            with lock:
                active -= 1

        with tempfile.TemporaryDirectory() as root:
            runner = Downloader(Path(root), True, 0, 0, workers=4)
            with patch("external_validation.scripts.download_real_omics.HMP2_PRODUCT_PAGES",
                       ("one.html",)), patch(
                "external_validation.scripts.download_real_omics.list_http_files", return_value=urls
            ), patch.object(runner, "https", side_effect=fake_download):
                runner.hmp2_all_products()
        self.assertEqual(len(calls), 9)  # eight products and metadata
        self.assertGreaterEqual(peak, 2)
        self.assertLessEqual(peak, 4)

    def test_http_range_resumes_incomplete_file(self):
        with tempfile.TemporaryDirectory() as root:
            dataset = Path(root) / "raw" / "covaar"
            dataset.mkdir(parents=True)
            (dataset / "sample.txt.gz.part").write_bytes(b"\x1f\x8b")
            runner = Downloader(Path(root), False, 1.0, 1.0)
            with patch(
                "external_validation.scripts.download_real_omics.urlopen",
                return_value=FakeResponse(b"\x08rest-of-payload", status=206, range_start=2),
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

