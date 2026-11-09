#!/usr/bin/env python3
"""Acquire public longitudinal omics resources to the WSL-mounted F: drive.

Enumerates files from each available publisher/source index rather than selecting
a few hand-picked matrices. No private credentials; resumable HTTPS transfers.
All data and manifests stay on the F: data drive, never in GitHub.
"""
from __future__ import annotations

import argparse
import errno
import hashlib
from http.client import HTTPException
import json
import os
import re
import shutil
import sys
import time
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse, urljoin, unquote
from urllib.request import Request, urlopen

ROOT = "/mnt/f/datasets/pulse-omics"
# The official IBDMDB results page links to all these HMP2 product pages.
# Follow their actual file hrefs; do not preselect one table per modality.
HMP2_PRODUCT_PAGES = (
    "products_16s_2018-01-08.html",
    "products_MBX_2017-03-01.html",
    "products_MGX_2017-08-12.html",
    "products_MTX_2017-12-14.html",
    "products_HMP2_HTX.html",
    "products_HMP2_Proteomics.html",
    "products_HMP2_vir.html",
    "products_HMP2_serology.html",
)
HMP2_RAW_PAGES = (
    "rawfiles_16s_2018-01-08.html",
    "rawfiles_MBX_2017-03-01.html",
    "rawfiles_MGX_2017-08-12.html",
    "rawfiles_MTX_HMP2.html",
    "rawfiles_proteomics_HMP2.html",
    "rawfiles_Viromics_HMP2.html",
)
HMP2_CATALOG_BASE = "https://www.ibdmdb.org/downloads/html/"
SU2020_RELEASE_ZIP = (
    "https://api.data.mendeley.com/datasets/"
    "tzydswhhb5/zip/file_downloaded?version=5"
)
SU2020_EUROPEPMC_ZIP = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/"
    "PMC7598382/supplementaryFiles"
)
# Cell's PMC version lists supplements mmc1.xlsx through mmc6.xlsx.
# mmc1.xlsx is already a known source in scripts/su2020_data.py.
SU2020_PMC_AWS_PREFIX = (
    "https://pmc-oa-opendata.s3.amazonaws.com/PMC7598382.1/"
)
GEO_SERIES = ("GSE325178",)
GEO_SECTIONS = ("suppl", "matrix", "soft", "miniml")


class AnchorLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            link = dict(attrs).get("href")
            if link:
                self.hrefs.append(link)


def list_http_files(url: str, allow_hosts: set[str]) -> list[str]:
    req = Request(url, headers={"User-Agent": "PULSE-Omics/1.0"})
    with urlopen(req, timeout=80) as response:
        raw = response.read(32 * 1024 * 1024)
    parser = AnchorLinks()
    parser.feed(raw.decode("utf-8", errors="replace"))
    candidates = [urljoin(url, href) for href in parser.hrefs]
    return sorted(set(link for link in candidates
                      if urlparse(link).scheme == "https"
                      and (urlparse(link).hostname or "") in allow_hosts
                      and not urlparse(link).path.endswith("/")))

HMP2_BASE = "https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb"
HMP2 = {
    "hmp2_metadata_2018-08-20.csv": (
        HMP2_BASE + "/metadata/hmp2_metadata_2018-08-20.csv"
    ),
    "HMP2_proteomics_1pep-pro_1p_FDR_5ppm.tsv.gz": (
        HMP2_BASE + "/products/HMP2/MPX/2017-03-20/"
        "HMP2_proteomics_1pep-pro_1p_FDR_5ppm.tsv.gz"
    ),
    "HMP2_metabolomics_w_metadata.biom.gz": (
        HMP2_BASE + "/products/HMP2/MBX/"
        "HMP2_metabolomics_w_metadata.biom.gz"
    ),
}
COVAAR = {
    "GSE325178_SomaScan_RFUs_allSamples.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/"
        "GSE325178/suppl/GSE325178_SomaScan_RFUs_allSamples.txt.gz"
    ),
    "GSE325178_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/"
        "GSE325178/matrix/GSE325178_series_matrix.txt.gz"
    ),
}
# Direct attachments from the publisher's actual supplementary-information list.
# Dataset 1/2 are Cytoscape sessions, not participant-level omics matrices.
P100_MEDIA = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.3870/"
    "MediaObjects/"
)
P100 = {
    f"41587_2017_BFnbt3870_MOESM{number}_ESM.{extension}":
        P100_MEDIA + f"41587_2017_BFnbt3870_MOESM{number}_ESM.{extension}"
    for number, extension in (
        (26, "pdf"), (27, "zip"), (28, "zip"), (29, "zip"),
        (30, "xlsx"), (31, "xlsx"), (32, "xlsx"), (33, "xlsx"),
        (34, "xlsx"), (35, "xlsx"), (36, "xlsx"),
    )
}
IPOP_SEASONAL_URL = "https://ndownloader.figshare.com/files/22807673"
CSIRO_LIST_API = (
    "https://data.csiro.au/dap/ws/v2/collections/77992/files/list"
)
FIGSHARE_ARTICLE = "https://api.figshare.com/v2/articles/31685023"
S3_PREFIXES = {
    "host_proteome": "ihmp/t2d/host/proteome/",
    "host_proteome_nonpride": "ihmp/t2d/host/proteome_nonpride/",
    "host_metabolome": "ihmp/t2d/host/metabolome/",
}
DATA_EXTENSIONS = (
    ".csv", ".tsv", ".txt", ".xlsx", ".xls", ".parquet",
    ".biom", ".biom.gz", ".tsv.gz", ".csv.gz", ".txt.gz",
    ".rds", ".h5", ".h5ad", ".h5mu",
)


def safe_path(base: Path, relative: str) -> Path:
    rel = PurePosixPath(relative)
    if rel.is_absolute() or ".." in rel.parts or not rel.parts:
        raise ValueError(f"Unsafe relative filename: {relative}")
    return base.joinpath(*rel.parts)


def check_format(path: Path) -> None:
    if path.stat().st_size == 0:
        raise ValueError("Empty downloaded file")
    with path.open("rb") as stream:
        header = stream.read(256)
    low = header.lower().lstrip()
    if low.startswith((b"<!doctype html", b"<html", b"<?xml", b"{\"error\"")):
        raise ValueError("Server returned an HTML or error response, not data")
    name = path.name.removesuffix(".part").lower()
    if name.endswith(".zip") and not header.startswith(b"PK"):
        raise ValueError("Expected a ZIP archive")
    if name.endswith((".gz", ".gzip")) and not header.startswith(b"\x1f\x8b"):
        raise ValueError("Expected gzip-compressed data")
    if name.endswith((".xlsx", ".xlsm")) and not header.startswith(b"PK"):
        raise ValueError("Expected an Excel workbook")


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for piece in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(piece)
    return digest.hexdigest()


class Downloader:
    def __init__(self, root: Path, dry_run: bool, max_file_gb: float, max_total_gb: float):
        self.root = root
        self.dry_run = dry_run
        # Zero means unlimited. A cap is opt-in, not a hidden file-selection rule.
        self.max_file = int(max_file_gb * 1024 ** 3) if max_file_gb else float("inf")
        self.max_total = int(max_total_gb * 1024 ** 3) if max_total_gb else float("inf")
        self.used = 0
        self.downloaded = 0
        self.existing = 0
        self.failed = 0
        self.report = {}
        self.inventory = {}
        self.incomplete = set()
        self.planned = {}
        self.stats = {}
        self._last_status = 0.0
        self._status_visible = False

    def progress(self, dataset: str, *, force: bool = False, note: str = "") -> None:
        """One in-place terminal status line; no per-file stdout messages."""
        if not sys.stdout.isatty():
            return
        now = time.monotonic()
        if not force and now - self._last_status < 0.5:
            return
        self._last_status = now
        stats = self.stats.get(dataset, {})
        done = sum(stats.get(k, 0) for k in ("downloaded", "existing", "failed", "candidate"))
        total = self.planned.get(dataset)
        count = f"{done}/{total}" if total is not None else str(done)
        size = self.used / (1024 ** 2)
        line = (
            f"{dataset.upper()} {count} | new {stats.get('downloaded', 0)} "
            f"| existing {stats.get('existing', 0)} | failed {stats.get('failed', 0)} "
            f"| transferred {size:.1f} MiB"
        )
        if note:
            line += f" | {note}"
        cols = shutil.get_terminal_size(fallback=(100, 20)).columns
        sys.stdout.write("\r\033[2K" + line[:max(1, cols - 1)])
        sys.stdout.flush()
        self._status_visible = True

    def finish_dataset(self, dataset: str) -> None:
        if self._status_visible:
            sys.stdout.write("\r\033[2K")
            sys.stdout.flush()
            self._status_visible = False
        stats = self.stats.get(dataset, {})
        done = sum(stats.get(k, 0) for k in ("downloaded", "existing", "failed", "candidate"))
        total = self.planned.get(dataset)
        checked = f"{done}/{total}" if total is not None else str(done)
        print(
            f"{dataset.upper()}: checked {checked}, downloaded {stats.get('downloaded', 0)}, "
            f"existing {stats.get('existing', 0)}, failed {stats.get('failed', 0)}"
        )

    def mark_discovery(self, dataset: str, source: str, found: int | None,
                       error: str | None = None) -> None:
        self.inventory.setdefault(dataset, []).append({
            "source": source, "listed_files": found,
            "status": "failed" if error else "listed", "error": error,
        })
        if error:
            self.incomplete.add(dataset)

    def record(self, dataset: str, entry: dict) -> None:
        self.report.setdefault(dataset, []).append(entry)
        status = entry.get("status")
        if status in ("downloaded", "existing", "failed", "candidate"):
            counts = self.stats.setdefault(dataset, {})
            counts[status] = counts.get(status, 0) + 1
        if not self.dry_run:
            folder = self.root / "raw" / dataset
            folder.mkdir(parents=True, exist_ok=True)
            with (folder / "download_events.jsonl").open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(entry, ensure_ascii=False) + "\n")
        self.progress(dataset)

    def flush(self, dataset: str) -> None:
        if self.dry_run:
            return
        out = self.root / "raw" / dataset / "download_manifest.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        entries = self.report.get(dataset, [])
        issues = [x for x in entries if x.get("status") in ("failed", "discovery_failed")]
        state = "incomplete" if issues or dataset in self.incomplete else "enumerated_in_this_run"
        payload = {
            "dataset": dataset,
            "updated_utc": datetime.now(timezone.utc).isoformat(),
            "status": state,
            "source_inventories": self.inventory.get(dataset, []),
            "files_checked_this_run": len(entries),
            "issues": len(issues),
            "files": entries,
        }
        temp = out.with_suffix(".json.part")
        temp.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        temp.replace(out)

    def https(self, dataset: str, relative: str, url: str) -> None:
        folder = self.root / "raw" / dataset
        dest = safe_path(folder, relative)
        if dest.exists():
            try:
                check_format(dest)
                sha = checksum(dest)
                self.existing += 1
                self.record(dataset, {
                    "file": relative, "source": url, "sha256": sha,
                    "bytes": dest.stat().st_size, "status": "existing",
                })
                return
            except (OSError, ValueError):
                pass
        if self.dry_run:
            self.record(dataset, {"status": "candidate", "file": relative, "source": url})
            return

        dest.parent.mkdir(parents=True, exist_ok=True)
        temp = dest.with_name(dest.name + ".part")
        max_attempts = 5
        for attempt in range(1, max_attempts + 1):
            try:
                already = temp.stat().st_size if temp.exists() else 0
                if already >= self.max_file:
                    raise ValueError("Incomplete transfer already exceeds file-size cap")

                headers = {"User-Agent": "PULSE-Omics/1.0 academic research"}
                if already:
                    headers["Range"] = f"bytes={already}-"
                request = Request(url, headers=headers)
                with urlopen(request, timeout=100) as response:
                    status = response.status
                    if status not in (200, 206):
                        raise ValueError(f"Unexpected HTTP status {status}")
                    if status == 206:
                        content_range = response.headers.get("Content-Range", "")
                        if not content_range.startswith(f"bytes {already}-"):
                            raise ValueError(
                                f"Invalid HTTP Content-Range for offset {already}: "
                                f"{content_range!r}"
                            )
                    # HTTP 200 means the remote endpoint ignored Range:
                    # restart this particular file, never append a full response.
                    append = status == 206 and already > 0
                    expected = int(response.headers.get("Content-Length") or 0)
                    if expected and expected + (already if append else 0) > self.max_file:
                        raise ValueError("File exceeds per-file size cap")

                    received = 0
                    with temp.open("ab" if append else "wb") as stream:
                        while block := response.read(1024 * 1024):
                            if stream.tell() + len(block) > self.max_file:
                                raise ValueError("File exceeds per-file size cap")
                            if self.used + len(block) > self.max_total:
                                raise ValueError("Session exceeds total download cap")
                            if shutil.disk_usage(dest.parent).free < len(block) + 10 * 1024 ** 3:
                                raise ValueError("Insufficient free disk space (10 GiB reserve)")
                            stream.write(block)
                            received += len(block)
                            self.used += len(block)
                            self.progress(dataset)
                    if expected and received != expected:
                        raise OSError(
                            f"HTTP body ended early: received {received} of {expected} bytes"
                        )

                check_format(temp)
                size, sha = temp.stat().st_size, checksum(temp)
                temp.replace(dest)
                self.downloaded += 1
                self.record(dataset, {
                    "file": relative, "source": url, "sha256": sha,
                    "bytes": size, "status": "downloaded",
                })
                return
            except (HTTPError, URLError, OSError, HTTPException, ValueError) as exc:
                # IncompleteRead is an HTTPException, not an OSError. The bytes
                # already written to .part are preserved and requested by Range.
                # Bytes in IncompleteRead.partial were not written, so the new
                # request must start at the actual on-disk file size.
                permanent_http = isinstance(exc, HTTPError) and exc.code not in {
                    408, 425, 429, 500, 502, 503, 504,
                }
                disk_error = isinstance(exc, OSError) and exc.errno in {
                    errno.ENOSPC, errno.EACCES, errno.EPERM, errno.EROFS,
                }
                transient = (
                    not permanent_http
                    and not disk_error
                    and not isinstance(exc, ValueError)
                )
                if transient and attempt < max_attempts:
                    partial_bytes = temp.stat().st_size if temp.exists() else 0
                    self.progress(dataset, force=True, note=f"retry {attempt}/{max_attempts}")
                    time.sleep(min(2 ** (attempt - 1), 8))
                    continue
                self.failed += 1
                self.incomplete.add(dataset)
                message = f"{type(exc).__name__}: {exc}"
                self.record(dataset, {
                    "file": relative, "source": url, "status": "failed",
                    "error": message, "attempts": attempt,
                    "partial_bytes": temp.stat().st_size if temp.exists() else 0,
                })
                return

    def hmp2_all_products(self, include_raw: bool = False) -> None:
        """Enumerate all HMP2 processed products, optionally source raw files."""
        listed = {}
        pages = HMP2_PRODUCT_PAGES + (HMP2_RAW_PAGES if include_raw else ())
        for page in pages:
            catalog_url = urljoin(HMP2_CATALOG_BASE, page)
            try:
                links = list_http_files(
                    catalog_url,
                    {"g-227ca.190ebd.75bc.data.globus.org"},
                )
                if not links:
                    raise ValueError("No downloadable Globus file links in source page")
                self.mark_discovery("hmp2", catalog_url, len(links))
                for url in links:
                    parsed = urlparse(url)
                    if not parsed.path.startswith("/ibdmdb/"):
                        continue
                    # Some official Globus hrefs contain "/ibdmdb//products/...".
                    # After removing the remote root they still begin with "/".
                    # Treat them as dataset-relative paths, never filesystem absolutes.
                    rel = unquote(parsed.path.removeprefix("/ibdmdb/")).lstrip("/")
                    if not rel:
                        continue
                    # Check every discovered path before starting any transfer.
                    safe_path(self.root / "raw" / "hmp2", rel)
                    # Reuse original three files already downloaded in earlier runs.
                    name = Path(rel).name
                    local = name if name in HMP2 else rel
                    listed[local] = url
                print(f"HMP2 {page}: {len(links)} source files enumerated")
            except Exception as exc:
                self.failed += 1
                self.mark_discovery("hmp2", catalog_url, None, str(exc))
                self.record("hmp2", {
                    "status": "discovery_failed", "source": catalog_url,
                    "error": f"{type(exc).__name__}: {exc}",
                })
        self.planned["hmp2"] = len(listed) + 1
        for name, url in sorted(listed.items()):
            self.https("hmp2", name, url)
        # Existing clinical metadata from the published HMP2 release.
        self.https("hmp2", "hmp2_metadata_2018-08-20.csv", HMP2["hmp2_metadata_2018-08-20.csv"])

    def covaar_all_geo(self) -> None:
        """Retrieve complete public series-level GEO files in every standard section."""
        for accession in GEO_SERIES:
            base = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
                    f"{accession[:-3]}nnn/{accession}/")
            for section in GEO_SECTIONS:
                page = base + section + "/"
                try:
                    urls = list_http_files(page, {"ftp.ncbi.nlm.nih.gov"})
                    # Standard FTP indexes have parent links; discard all items
                    # outside the exact section and directory.
                    urls = [u for u in urls if urlparse(u).path.startswith(
                        f"/geo/series/{accession[:-3]}nnn/{accession}/{section}/"
                    )]
                    self.mark_discovery("covaar", page, len(urls))
                    if not urls and section in ("suppl", "matrix"):
                        raise ValueError("No files listed for mandatory GEO section")
                    print(f"COVAAR {accession}/{section}: {len(urls)} source files")
                    for url in urls:
                        name = unquote(Path(urlparse(url).path).name)
                        relative = f"geo/{name}" if name in COVAAR else f"geo/{section}/{name}"
                        self.https("covaar", relative, url)
                except Exception as exc:
                    self.failed += 1
                    self.mark_discovery("covaar", page, None, str(exc))
                    self.record("covaar", {
                        "status": "discovery_failed", "source": page,
                        "error": f"{type(exc).__name__}: {exc}",
                    })
        # Retain earlier published filenames on disk rather than overwrite them.
        for name, url in COVAAR.items():
            self.https("covaar", f"geo/{name}", url)

    def su2020_full_release(self) -> None:
        """Acquire all published Mendeley v5 tables in one source archive."""
        self.https(
            "su2020", "mendeley_v5_complete_supplementary.zip",
            SU2020_RELEASE_ZIP,
        )
        self.mark_discovery(
            "su2020", "https://data.mendeley.com/datasets/tzydswhhb5/5",
            6,
        )
        # Preserve a second public source even if the Mendeley API requires
        # an authenticated token or redirects to a blocked object.
        self.https(
            "su2020", "europepmc_supplementary.zip",
            SU2020_EUROPEPMC_ZIP,
        )
        # Individual paper supplements avoid dependence on Mendeley's
        # OAuth-gated bulk archive and Europe PMC's HTML response.
        for number in range(1, 7):
            name = f"mmc{number}.xlsx"
            self.https("su2020", f"pmc/{name}", SU2020_PMC_AWS_PREFIX + name)

    def ipop(self) -> None:
        self.figshare_article("ipop", 12376508, "seasonal")
        try:
            import boto3
            from botocore import UNSIGNED
            from botocore.config import Config
            client = boto3.client(
                "s3", config=Config(
                    signature_version=UNSIGNED,
                    retries={"max_attempts": 3, "mode": "standard"},
                ),
            )
            paginator = client.get_paginator("list_objects_v2")
            entries = []
            for label, prefix in S3_PREFIXES.items():
                for page in paginator.paginate(Bucket="hmpdcc", Prefix=prefix):
                    for obj in page.get("Contents", []):
                        key = obj["Key"]
                        if key.endswith("/"):
                            continue
                        size = int(obj["Size"])
                        relative = key[len(prefix):]
                        if not relative:
                            continue
                        if size > self.max_file:
                            self.incomplete.add("ipop")
                                    continue
                        entries.append((size, label, relative, key))
            self.mark_discovery("ipop", "s3://hmpdcc/ihmp/t2d/host/", len(entries))
            print(f"iPOP public host objects enumerated: {len(entries)}")
            if not entries:
                raise ValueError("No candidate processed omics objects found")
            # Small processed tables first; do not copy arbitrary raw MS/sequencing data.
            for size, label, relative, key in sorted(entries):
                if self.used + size > self.max_total:
                    self.incomplete.add("ipop")
                    continue
                self.s3("ipop", f"{label}/{relative}", client, key, size)
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
            self.failed += 1
            self.mark_discovery("ipop", "s3://hmpdcc/ihmp/t2d/host/", None, message)
            print(f"iPOP S3 discovery unavailable: {message}; using the published "
                  "seasonal-study processed data", file=sys.stderr)
            self.record("ipop", {
                "status": "discovery_failed", "source": "s3://hmpdcc/ihmp/t2d/host/",
                "error": message,
            })
            # Sailani et al. 2020 explicitly links Figshare article 12376508 v1.
            # This is an annotated processed release, not a copy of all S3 data.
            self.https("ipop", "seasonal/Multi_Omics_Seasonal.RData", IPOP_SEASONAL_URL)

    def s3(self, dataset: str, relative: str, client, key: str, size: int) -> None:
        dest = safe_path(self.root / "raw" / dataset, relative)
        uri = "s3://hmpdcc/" + key
        if dest.exists():
            try:
                check_format(dest)
                if dest.stat().st_size == size:
                    self.existing += 1
                    self.record(dataset, {
                        "file": relative, "source": uri,
                        "bytes": size, "sha256": checksum(dest),
                        "status": "existing",
                    })
                    return
            except (ValueError, OSError):
                pass
        if self.dry_run:
            self.record(dataset, {"status": "candidate", "file": relative, "source": uri, "bytes": size})
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        temp = dest.with_name(dest.name + ".part")
        try:
            offset = temp.stat().st_size if temp.exists() else 0
            if offset > size:
                temp.unlink()
                offset = 0
            kwargs = {"Bucket": "hmpdcc", "Key": key}
            if offset:
                kwargs["Range"] = f"bytes={offset}-"
            response = client.get_object(**kwargs)
            body = response["Body"]
            try:
                with temp.open("ab" if offset else "wb") as stream:
                    while block := body.read(1024 * 1024):
                        if self.used + len(block) > self.max_total:
                            raise ValueError("Session exceeds total download cap")
                        stream.write(block)
                        self.used += len(block)
            finally:
                body.close()
            if temp.stat().st_size != size:
                raise ValueError(f"Object size mismatch: expected {size}")
            check_format(temp)
            sha = checksum(temp)
            temp.replace(dest)
            self.downloaded += 1
            self.record(dataset, {
                "file": relative, "source": uri,
                "bytes": size, "sha256": sha, "status": "downloaded",
            })
        except Exception as exc:
            self.failed += 1
            message = f"{type(exc).__name__}: {exc}"
            self.record(dataset, {
                "file": relative, "source": uri,
                "status": "failed", "error": message,
            })

    def figshare_article(self, dataset: str, article_id: int, subdir: str) -> None:
        api = f"https://api.figshare.com/v2/articles/{article_id}"
        try:
            req = Request(api, headers={"User-Agent": "PULSE-Omics/1.0"})
            with urlopen(req, timeout=40) as response:
                metadata = json.load(response)
            files = metadata.get("files", [])
            if not files:
                raise ValueError("Figshare record has no public downloadable files")
            self.mark_discovery(dataset, api, len(files))
            for item in files:
                name = item.get("name")
                url = item.get("download_url")
                if name and url and urlparse(url).scheme == "https":
                    self.https(dataset, f"{subdir}/{Path(name).name}", url)
                else:
                    raise ValueError("Figshare listed file without a usable public URL")
        except Exception as exc:
            self.failed += 1
            self.mark_discovery(dataset, api, None, str(exc))
            message = f"{type(exc).__name__}: {exc}"
            print(f"FAILED {dataset} Figshare discovery: {message}", file=sys.stderr)
            self.record(dataset, {
                "status": "discovery_failed", "source": api, "error": message,
            })

    def figshare(self) -> None:
        self.figshare_article("covaar", 31685023, "figshare")

    def csiro(self) -> None:
        # Discover first: collection numeric ID and public file list are API-specific.
        # Never treat HTML errors or arbitrary portal content as study data.
        try:
            req = Request(CSIRO_LIST_API, headers={
                "User-Agent": "PULSE-Omics/1.0",
                "Accept": "text/plain, application/json",
            })
            with urlopen(req, timeout=40) as response:
                body = response.read(5 * 1024 * 1024)
            text = body.decode("utf-8", errors="replace")
            urls = sorted(set(re.findall(r"https://[^\s<>\"']+", text)))
            downloads = [
                url.rstrip("),;") for url in urls
                if (urlparse(url).hostname or "").endswith(
                    ("csiro.au", "csirocdn.com.au")
                ) and Path(urlparse(url).path).suffix.lower() in {
                    ".zip", ".xlsx", ".xls", ".csv", ".tsv", ".txt",
                    ".pdf", ".gz", ".rdata", ".rds", ".h5", ".biom",
                    ".mzml", ".raw", ".7z", ".tar", ".parquet",
                }
            ]
            if not downloads:
                raise ValueError(
                    "CSIRO collection file list did not expose direct public file URLs"
                )
            self.mark_discovery("csiro", CSIRO_LIST_API, len(downloads))
            for index, url in enumerate(downloads):
                name = Path(urlparse(url).path).name
                self.https("csiro", f"{index:03d}_{name}", url)
        except Exception as exc:
            self.failed += 1
            message = f"{type(exc).__name__}: {exc}"
            self.mark_discovery("csiro", CSIRO_LIST_API, None, message)
            print(f"CSIRO collection discovery unresolved: {message}", file=sys.stderr)
            self.record("csiro", {
                "status": "discovery_failed", "source": CSIRO_LIST_API,
                "error": message,
            })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=os.environ.get("PULSE_OMICS_DATA_ROOT", ROOT))
    parser.add_argument("--datasets", nargs="+", choices=(
        "su2020", "ipop", "hmp2", "covaar", "p100", "csiro",
    ), default=["su2020", "ipop", "hmp2", "covaar", "p100", "csiro"])
    parser.add_argument("--dry-run", action="store_true", help="Discover only, do not download")
    parser.add_argument("--include-raw", action="store_true",
                        help="Also fetch publicly indexed HMP2 raw FASTQ/instrument files")
    parser.add_argument("--max-file-gb", type=float, default=0.0,
                        help="Optional maximum individual file size in GiB; 0 means unlimited")
    parser.add_argument("--max-total-gb", type=float, default=0.0,
                        help="Optional total transfer cap in GiB; 0 means unlimited")
    args = parser.parse_args()
    if args.max_file_gb < 0 or args.max_total_gb < 0:
        parser.error("Download size limits must be nonnegative")
    runner = Downloader(
        Path(args.data_root).expanduser(),
        args.dry_run,
        args.max_file_gb,
        args.max_total_gb,
    )
    for dataset in dict.fromkeys(args.datasets):
        print(f"\n=== {dataset.upper()} ===", flush=True)
        if dataset == "su2020":
            runner.su2020_full_release()
        elif dataset == "ipop":
            runner.ipop()
        elif dataset == "hmp2":
            runner.hmp2_all_products(include_raw=args.include_raw)
        elif dataset == "covaar":
            runner.covaar_all_geo()
            runner.figshare()
        elif dataset == "p100":
            runner.mark_discovery(
                "p100", "https://doi.org/10.1038/nbt.3870",
                len(P100),
            )
            for name, url in P100.items():
                runner.https("p100", name, url)
            # Public publisher supplements have been enumerated, but individual
            # participant-level P100 omics are a separate controlled dbGaP release.
            runner.incomplete.add("p100")
            runner.record("p100", {
                "status": "controlled_access_not_downloaded",
                "source": "https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001363.v1.p1",
                "detail": "Individual-level omics require dbGaP authorization",
            })
        elif dataset == "csiro":
            runner.csiro()
        runner.flush(dataset)
        runner.finish_dataset(dataset)
    print(
        f"\nDownloads: {runner.downloaded}, existing: {runner.existing}, "
        f"failures: {runner.failed}; data root: {runner.root / 'raw'}"
    )
    if runner.failed or runner.incomplete:
        print("Incomplete source coverage or failed downloads: "
              + ", ".join(sorted(runner.incomplete)))
        sys.exit(2)


if __name__ == "__main__":
    main()

