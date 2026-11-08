#!/usr/bin/env python3
"""Fetch and inspect Su et al. 2020 Table S1 without guessing its schema.

Raw data: /mnt/f/datasets/pulse-omics/raw/su2020/
Inspection: /mnt/f/datasets/pulse-omics/processed/su2020/
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen
from zipfile import BadZipFile, ZipFile

DATASET_DOI = "10.17632/tzydswhhb5.5"
DATASET_PAGE = "https://data.mendeley.com/datasets/tzydswhhb5/5"
PMC_S1_URL = "https://pmc.ncbi.nlm.nih.gov/articles/instance/7754857/bin/mmc1.xlsx"
MENDELEY_API = "https://api.data.mendeley.com/datasets/tzydswhhb5"
DEFAULT_ROOT = "/mnt/f/datasets/pulse-omics"


def data_root(value: str | None) -> Path:
    return Path(value or os.environ.get("PULSE_OMICS_DATA_ROOT", DEFAULT_ROOT)).expanduser()


def request(url: str):
    return Request(url, headers={"User-Agent": "PULSE-Omics/1.0 (academic dataset retrieval)", "Accept": "*/*"})


def mendeley_table_url() -> str:
    """Resolve a v5 Table S1 file UUID using the current public Mendeley API."""
    with urlopen(request(MENDELEY_API + "/files?version=5"), timeout=40) as response:
        metadata = json.load(response)
    if isinstance(metadata, dict):
        metadata = metadata.get("files", metadata.get("results", []))
    if not isinstance(metadata, list):
        raise RuntimeError("Unexpected Mendeley file metadata response")
    matching = []
    for record in metadata:
        name = str(record.get("filename") or record.get("file_name") or record.get("name") or "")
        if name.lower().startswith("table s1") and name.lower().endswith(".xlsx"):
            matching.append(record)
    if len(matching) != 1:
        raise RuntimeError(f"Expected one Table S1 workbook in Mendeley v5, found {len(matching)}")
    file_id = matching[0].get("id") or matching[0].get("file_id")
    if not file_id:
        raise RuntimeError("Mendeley metadata is missing the file UUID")
    return f"{MENDELEY_API}/files/{quote(str(file_id))}/file_downloaded?version=5"


def fetch(destination: Path, source: str, overwrite: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        print(f"Already present: {destination}. Use --overwrite to retrieve again.")
        validate_workbook(destination)
        return

    attempts = [("PMC", lambda: PMC_S1_URL), ("Mendeley", mendeley_table_url)]
    if source != "auto":
        attempts = [pair for pair in attempts if pair[0].lower() == source]
    errors = []
    for label, url_factory in attempts:
        temp = destination.with_suffix(destination.suffix + ".part")
        temp.unlink(missing_ok=True)
        try:
            url = url_factory()
            sha = hashlib.sha256()
            total = 0
            with urlopen(request(url), timeout=90) as response, temp.open("wb") as output:
                while block := response.read(1024 * 1024):
                    output.write(block)
                    sha.update(block)
                    total += len(block)
            validate_workbook(temp)
            temp.replace(destination)
            record = {
                "dataset": "Su et al. 2020 Table S1",
                "dataset_doi": DATASET_DOI,
                "dataset_page": DATASET_PAGE,
                "retrieval_source": label,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                "file": str(destination),
                "bytes": total,
                "sha256": sha.hexdigest(),
            }
            (destination.parent / "download_manifest.json").write_text(
                json.dumps(record, indent=2) + "\n", encoding="utf-8"
            )
            print(f"Downloaded: {destination} ({total} bytes)")
            print(f"SHA256: {sha.hexdigest()}")
            return
        except Exception as exc:
            temp.unlink(missing_ok=True)
            errors.append(f"{label}: {type(exc).__name__}: {exc}")
    raise RuntimeError("All download sources failed:\n" + "\n".join(errors))


def validate_workbook(path: Path) -> None:
    try:
        with ZipFile(path) as archive:
            names = set(archive.namelist())
            if "xl/workbook.xml" not in names or "[Content_Types].xml" not in names:
                raise ValueError("Downloaded file is not an Excel workbook")
    except BadZipFile as exc:
        raise ValueError(f"Downloaded file is not a valid XLSX: {path}") from exc


def compact_cell(value: object) -> str | None:
    if value is None:
        return None
    value = str(value).strip()
    return value[:120] if value else None


def inspect(workbook: Path, output: Path, preview_rows: int, preview_cols: int) -> None:
    try:
        from openpyxl import load_workbook
    except ImportError as exc:
        raise SystemExit("Install openpyxl: python -m pip install -r requirements-omics.txt") from exc
    validate_workbook(workbook)
    report = {
        "source_file": str(workbook),
        "source_sha256": hashlib.sha256(workbook.read_bytes()).hexdigest(),
        "dataset_doi": DATASET_DOI,
        "sheets": [],
    }
    with load_workbook(workbook, read_only=True, data_only=True) as excel:
        for sheet in excel:
            preview = [
                [compact_cell(value) for value in row]
                for row in sheet.iter_rows(
                    min_row=1, max_row=min(preview_rows, sheet.max_row or preview_rows),
                    min_col=1, max_col=min(preview_cols, sheet.max_column or preview_cols),
                    values_only=True,
                )
            ]
            info = {
                "name": sheet.title,
                "max_row": sheet.max_row,
                "max_column": sheet.max_column,
                "top_left_preview": preview,
            }
            report["sheets"].append(info)
            print(f"{sheet.title}: {sheet.max_row} rows x {sheet.max_column} cols")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Workbook structure saved to: {output}")
    print("No protein/metabolite matrices or patient visit matches have been inferred.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=None, help="Default: PULSE_OMICS_DATA_ROOT or /mnt/f/datasets/pulse-omics")
    commands = parser.add_subparsers(dest="command", required=True)
    fetch_cmd = commands.add_parser("fetch", help="Download Su 2020 v5 Table S1 to F: drive")
    fetch_cmd.add_argument("--source", choices=["auto", "pmc", "mendeley"], default="auto")
    fetch_cmd.add_argument("--overwrite", action="store_true")
    inspect_cmd = commands.add_parser("inspect", help="Report real workbook sheets and cell layout")
    inspect_cmd.add_argument("--preview-rows", type=int, default=8)
    inspect_cmd.add_argument("--preview-cols", type=int, default=12)
    args = parser.parse_args()
    root = data_root(args.data_root)
    workbook = root / "raw" / "su2020" / "table_s1.xlsx"
    if args.command == "fetch":
        fetch(workbook, args.source, args.overwrite)
    else:
        if not workbook.is_file():
            parser.error(f"Workbook not found: {workbook}; run 'fetch' first")
        if args.preview_rows < 1 or args.preview_cols < 1:
            parser.error("Preview dimensions must be positive")
        inspect(
            workbook,
            root / "processed" / "su2020" / "workbook_structure.json",
            args.preview_rows,
            args.preview_cols,
        )


if __name__ == "__main__":
    main()
