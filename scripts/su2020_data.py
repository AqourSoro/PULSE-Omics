#!/usr/bin/env python3
"""Acquire and inspect Su et al. (2020) Table S1 on the F: data drive."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from zipfile import BadZipFile, ZipFile

DATASET_DOI = "10.17632/tzydswhhb5.5"
DATASET_PAGE = "https://data.mendeley.com/datasets/tzydswhhb5/5"
# Europe PMC delivers an archive containing the paper's supplementary files.
EUROPE_PMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7598382/supplementaryFiles"
# PMC moved its article datasets to AWS in August 2026; availability depends on license.
PMC_AWS_URL = "https://pmc-oa-opendata.s3.amazonaws.com/PMC7598382.1/mmc1.xlsx"
DEFAULT_ROOT = "/mnt/f/datasets/pulse-omics"
MAX_ARCHIVE_BYTES = 128 * 1024 * 1024
MAX_WORKBOOK_BYTES = 64 * 1024 * 1024


def data_root(value: str | None) -> Path:
    return Path(value or os.environ.get("PULSE_OMICS_DATA_ROOT", DEFAULT_ROOT)).expanduser()


def request(url: str) -> Request:
    return Request(
        url,
        headers={
            "User-Agent": "PULSE-Omics/1.0 academic research dataset retrieval",
            "Accept": "application/zip, application/octet-stream, */*",
        },
    )


def validate_workbook(path: Path) -> None:
    try:
        with ZipFile(path) as archive:
            names = set(archive.namelist())
            if "xl/workbook.xml" not in names or "[Content_Types].xml" not in names:
                raise ValueError("File is a ZIP but not an XLSX workbook")
    except BadZipFile as exc:
        with path.open("rb") as handle:
            preview = handle.read(100)
        if b"<html" in preview.lower() or b"<!doctype" in preview.lower():
            raise ValueError("Server returned an HTML page instead of an XLSX file") from exc
        raise ValueError(f"Not a valid XLSX: {path}") from exc


def obtain_workbook(source_path: Path, target: Path) -> None:
    """Read an XLSX, or extract mmc1.xlsx / Table S1 from a ZIP without unzip-all."""
    try:
        validate_workbook(source_path)
        if source_path.stat().st_size > MAX_WORKBOOK_BYTES:
            raise ValueError("Workbook exceeds the 64 MiB safety limit")
        shutil.copyfile(source_path, target)
        return
    except ValueError as exc:
        workbook_error = str(exc)
    try:
        with ZipFile(source_path) as archive:
            choices = [
                entry for entry in archive.infolist()
                if not entry.is_dir()
                and (Path(entry.filename).name.lower() == "mmc1.xlsx"
                     or Path(entry.filename).name.lower().startswith("table s1"))
                and entry.filename.lower().endswith(".xlsx")
            ]
            if len(choices) != 1:
                raise ValueError(f"Expected exactly one Table S1 workbook in archive; found {len(choices)}")
            chosen = choices[0]
            if chosen.file_size > MAX_WORKBOOK_BYTES:
                raise ValueError("Workbook in archive exceeds the 64 MiB safety limit")
            with archive.open(chosen) as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
        validate_workbook(target)
    except BadZipFile as exc:
        raise ValueError(f"Neither XLSX nor usable supplementary ZIP ({workbook_error})") from exc


def install_workbook(source_file: Path, destination: Path, provenance: str, overwrite: bool) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        validate_workbook(destination)
        print(f"Already present and valid: {destination}")
        print("Use --overwrite if replacement is intended.")
        return
    tmp = destination.with_name(destination.name + ".tmp")
    tmp.unlink(missing_ok=True)
    try:
        obtain_workbook(source_file, tmp)
        with tmp.open("rb") as handle:
            sha = hashlib.file_digest(handle, "sha256").hexdigest()
        size = tmp.stat().st_size
        tmp.replace(destination)
        manifest = {
            "dataset": "Su et al. 2020 Table S1",
            "doi": DATASET_DOI,
            "dataset_page": DATASET_PAGE,
            "retrieval_source": provenance,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(),
            "path": str(destination),
            "bytes": size,
            "sha256": sha,
        }
        (destination.parent / "download_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        print(f"Workbook installed: {destination} ({size} bytes, SHA-256 {sha})")
    finally:
        tmp.unlink(missing_ok=True)


def fetch(destination: Path, source: str, overwrite: bool) -> None:
    if destination.exists() and not overwrite:
        validate_workbook(destination)
        print(f"Already present: {destination}")
        return
    sources = {
        "europepmc": EUROPE_PMC_URL,
        "pmc-aws": PMC_AWS_URL,
    }
    candidates = list(sources.items()) if source == "auto" else [(source, sources[source])]
    errors: list[str] = []
    destination.parent.mkdir(parents=True, exist_ok=True)
    for label, url in candidates:
        raw = destination.with_name(destination.name + ".download")
        raw.unlink(missing_ok=True)
        try:
            with urlopen(request(url), timeout=90) as response, raw.open("wb") as output:
                total = 0
                while block := response.read(1024 * 1024):
                    total += len(block)
                    if total > MAX_ARCHIVE_BYTES:
                        raise ValueError("Download exceeds the 128 MiB safety limit")
                    output.write(block)
            install_workbook(raw, destination, f"{label}: {url}", overwrite=True)
            return
        except (HTTPError, OSError, ValueError, BadZipFile) as exc:
            errors.append(f"{label}: {type(exc).__name__}: {exc}")
        finally:
            raw.unlink(missing_ok=True)
    raise RuntimeError(
        "No automated source supplied a valid Table S1 workbook:\n"
        + "\n".join(errors)
        + "\nIf publisher download succeeds in your browser, run:\n"
        + "  python scripts/su2020_data.py import --file /path/to/Table_S1.xlsx\n"
        + f"Dataset page: {DATASET_PAGE}"
    )


def compact_cell(value: object) -> str | None:
    if value is None:
        return None
    string = str(value).strip()
    return string[:120] if string else None


def inspect(workbook: Path, output: Path, preview_rows: int, preview_cols: int) -> None:
    try:
        from openpyxl import load_workbook
    except ImportError as exc:
        raise SystemExit("Install: python -m pip install -r requirements-omics.txt") from exc
    validate_workbook(workbook)
    with workbook.open("rb") as handle:
        sha = hashlib.file_digest(handle, "sha256").hexdigest()
    report = {
        "source_file": str(workbook),
        "source_sha256": sha,
        "dataset_doi": DATASET_DOI,
        "sheets": [],
    }
    excel = load_workbook(workbook, read_only=True, data_only=True)
    try:
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
    finally:
        excel.close()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Workbook structure saved to: {output}")
    print("Patient/visit and assay linkage are not inferred from spreadsheet position.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=None, help=f"Data root (default: {DEFAULT_ROOT})")
    commands = parser.add_subparsers(dest="command", required=True)
    fetch_cmd = commands.add_parser("fetch", help="Fetch published supplements without API credentials")
    fetch_cmd.add_argument("--source", choices=["auto", "europepmc", "pmc-aws"], default="auto")
    fetch_cmd.add_argument("--overwrite", action="store_true")
    import_cmd = commands.add_parser("import", help="Validate and place a browser-downloaded XLSX/ZIP on F:")
    import_cmd.add_argument("--file", type=Path, required=True)
    import_cmd.add_argument("--overwrite", action="store_true")
    inspect_cmd = commands.add_parser("inspect", help="Inspect real worksheet dimensions and sample identifier hints")
    inspect_cmd.add_argument("--preview-rows", type=int, default=8)
    inspect_cmd.add_argument("--preview-cols", type=int, default=12)
    args = parser.parse_args()
    root = data_root(args.data_root)
    workbook = root / "raw" / "su2020" / "table_s1.xlsx"
    if args.command == "fetch":
        fetch(workbook, args.source, args.overwrite)
    elif args.command == "import":
        source = args.file.expanduser().resolve()
        if not source.is_file():
            parser.error(f"Source file does not exist: {source}")
        if source == workbook.resolve():
            parser.error("Import source is already the destination; run 'inspect'")
        install_workbook(source, workbook, f"local import: {source}", args.overwrite)
    else:
        if not workbook.is_file():
            parser.error(f"Workbook not found: {workbook}; run 'fetch' or 'import' first")
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
