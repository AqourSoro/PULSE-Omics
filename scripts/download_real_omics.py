#!/usr/bin/env python3
"""Acquire public longitudinal omics resources to the WSL-mounted F: drive.

Data catalog and biological interpretation: docs/real_omics_data.md
No private credentials, anonymous AWS S3, resumable HTTPS transfers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen

ROOT = "/mnt/f/datasets/pulse-omics"
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
        self.max_file = int(max_file_gb * 1024 ** 3)
        self.max_total = int(max_total_gb * 1024 ** 3)
        self.used = 0
        self.downloaded = 0
        self.existing = 0
        self.failed = 0
        self.report = {}

    def record(self, dataset: str, entry: dict) -> None:
        self.report.setdefault(dataset, []).append(entry)
        manifest = self.root / "raw" / dataset / "download_manifest.json"
        if not self.dry_run:
            manifest.parent.mkdir(parents=True, exist_ok=True)
            manifest.write_text(
                json.dumps({
                    "dataset": dataset,
                    "updated_utc": datetime.now(timezone.utc).isoformat(),
                    "files": self.report[dataset],
                }, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )

    def https(self, dataset: str, relative: str, url: str) -> None:
        folder = self.root / "raw" / dataset
        dest = safe_path(folder, relative)
        if dest.exists():
            try:
                check_format(dest)
                sha = checksum(dest)
                self.existing += 1
                print(f"EXISTS {dataset}/{relative} [{dest.stat().st_size} bytes]")
                self.record(dataset, {
                    "file": relative, "source": url, "sha256": sha,
                    "bytes": dest.stat().st_size, "status": "existing"
                })
                return
            except (OSError, ValueError) as exc:
                print(f"INVALID {dataset}/{relative}: {exc}; retrying")
        if self.dry_run:
            print(f"CANDIDATE {dataset}/{relative} <= {url}")
            return
        dest.parent.mkdir(parents=True, exist_ok=True)
        temp = dest.with_name(dest.name + ".part")
        try:
            already = temp.stat().st_size if temp.exists() else 0
            if already >= self.max_file:
                raise ValueError("Incomplete transfer already exceeds file-size cap")
            request_headers = {"User-Agent": "PULSE-Omics/1.0 academic research"}
            if already:
                request_headers["Range"] = f"bytes={already}-"
            req = Request(url, headers=request_headers)
            with urlopen(req, timeout=100) as response:
                status = response.status
                if status not in (200, 206):
                    raise ValueError(f"Unexpected HTTP status {status}")
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
                        stream.write(block)
                        received += len(block)
                        self.used += len(block)
                if expected and received != expected:
                    raise ValueError(
                        f"Incomplete HTTP body: received {received} of {expected} bytes"
                    )
            check_format(temp)
            size, sha = temp.stat().st_size, checksum(temp)
            temp.replace(dest)
            self.downloaded += 1
            print(f"OK {dataset}/{relative} [{size} bytes] sha256={sha}")
            self.record(dataset, {
                "file": relative, "source": url, "sha256": sha,
                "bytes": size, "status": "downloaded",
            })
        except (HTTPError, URLError, OSError, ValueError) as exc:
            self.failed += 1
            message = f"{type(exc).__name__}: {exc}"
            print(f"FAILED {dataset}/{relative}: {message}", file=sys.stderr)
            self.record(dataset, {
                "file": relative, "source": url,
                "status": "failed", "error": message,
            })

    def ipop(self) -> None:
        processed = self.root / "raw/ipop/seasonal/Multi_Omics_Seasonal.RData"
        if processed.exists():
            self.https("ipop", "seasonal/Multi_Omics_Seasonal.RData", IPOP_SEASONAL_URL)
            return
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
                        if not relative or not key.lower().endswith(DATA_EXTENSIONS):
                            continue
                        if size > self.max_file:
                            continue
                        entries.append((size, label, relative, key))
            print(f"iPOP candidate processed omics objects: {len(entries)}")
            if not entries:
                raise ValueError("No candidate processed omics objects found")
            # Small processed tables first; do not copy arbitrary raw MS/sequencing data.
            for size, label, relative, key in sorted(entries):
                if self.used + size > self.max_total:
                    print(f"SKIP iPOP/{key}: total download budget")
                    continue
                self.s3("ipop", f"{label}/{relative}", client, key, size)
        except Exception as exc:
            message = f"{type(exc).__name__}: {exc}"
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
            print(f"CANDIDATE {dataset}/{relative} [{size} bytes] <= {uri}")
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
            print(f"OK {dataset}/{relative} [{size} bytes] sha256={sha}")
            self.record(dataset, {
                "file": relative, "source": uri,
                "bytes": size, "sha256": sha, "status": "downloaded",
            })
        except Exception as exc:
            self.failed += 1
            message = f"{type(exc).__name__}: {exc}"
            print(f"FAILED {dataset}/{relative}: {message}", file=sys.stderr)
            self.record(dataset, {
                "file": relative, "source": uri,
                "status": "failed", "error": message,
            })

    def figshare(self) -> None:
        try:
            req = Request(FIGSHARE_ARTICLE, headers={"User-Agent": "PULSE-Omics/1.0"})
            with urlopen(req, timeout=40) as response:
                metadata = json.load(response)
            files = metadata.get("files", [])
            if not files:
                raise ValueError("Figshare article has no downloadable files")
            for item in files:
                name = item.get("name")
                url = item.get("download_url")
                if name and url and urlparse(url).scheme == "https":
                    self.https("covaar", f"figshare/{Path(name).name}", url)
        except Exception as exc:
            self.failed += 1
            message = f"{type(exc).__name__}: {exc}"
            print(f"FAILED COVAAR Figshare discovery: {message}", file=sys.stderr)
            self.record("covaar", {"status": "discovery_failed", "error": message})

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
                ) and urlparse(url).path.lower().endswith(DATA_EXTENSIONS + (".zip",))
            ]
            if not downloads:
                raise ValueError(
                    "CSIRO collection file list did not expose direct public omics URLs"
                )
            for index, url in enumerate(downloads):
                name = Path(urlparse(url).path).name
                self.https("csiro", f"{index:03d}_{name}", url)
        except Exception as exc:
            self.failed += 1
            message = f"{type(exc).__name__}: {exc}"
            print(f"CSIRO collection discovery unresolved: {message}", file=sys.stderr)
            self.record("csiro", {
                "status": "discovery_failed", "source": CSIRO_LIST_API,
                "error": message,
            })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=os.environ.get("PULSE_OMICS_DATA_ROOT", ROOT))
    parser.add_argument("--datasets", nargs="+", choices=(
        "ipop", "hmp2", "covaar", "p100", "csiro",
    ), default=["ipop", "hmp2", "covaar", "p100", "csiro"])
    parser.add_argument("--dry-run", action="store_true", help="Discover only, do not download")
    parser.add_argument("--max-file-gb", type=float, default=2.0)
    parser.add_argument("--max-total-gb", type=float, default=12.0)
    args = parser.parse_args()
    if args.max_file_gb <= 0 or args.max_total_gb <= 0:
        parser.error("Download size limits must be positive")
    runner = Downloader(
        Path(args.data_root).expanduser(),
        args.dry_run,
        args.max_file_gb,
        args.max_total_gb,
    )
    for dataset in dict.fromkeys(args.datasets):
        print(f"\n=== {dataset.upper()} ===", flush=True)
        if dataset == "ipop":
            runner.ipop()
        elif dataset == "hmp2":
            for name, url in HMP2.items():
                runner.https("hmp2", name, url)
        elif dataset == "covaar":
            for name, url in COVAAR.items():
                runner.https("covaar", f"geo/{name}", url)
            runner.figshare()
        elif dataset == "p100":
            for name, url in P100.items():
                runner.https("p100", name, url)
        elif dataset == "csiro":
            runner.csiro()
    print(
        f"\nDownloads: {runner.downloaded}, existing: {runner.existing}, "
        f"failures: {runner.failed}; data root: {runner.root / 'raw'}"
    )
    if runner.failed:
        sys.exit(2)


if __name__ == "__main__":
    main()

