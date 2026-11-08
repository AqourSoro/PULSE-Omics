#!/usr/bin/env python3
"""Inventory public omics files on the F: drive without loading large matrices.

Reads a small preview, compression headers, and worksheet/archive metadata.
Never guesses participant/visit joins or changes source files.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
from zipfile import BadZipFile, ZipFile

DEFAULT_ROOT = "/mnt/f/datasets/pulse-omics"
IGNORE = {"download_manifest.json", "data_inventory.json"}
LIMIT = 32768


def preview_delimited(path: Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    delimiter = "\t" if ".tsv" in path.name or ".txt" in path.name else ","
    with opener(path, "rt", encoding="utf-8-sig", errors="replace", newline="") as source:
        rows = []
        for row in csv.reader(source, delimiter=delimiter):
            rows.append([field[:90] for field in row[:15]])
            if len(rows) == 5:
                break
    return {
        "delimiter": "\\t" if delimiter == "\t" else ",",
        "first_rows": rows,
        "header_columns_in_preview": len(rows[0]) if rows else 0,
    }


def preview_excel(path: Path) -> dict:
    from openpyxl import load_workbook
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        sheets = []
        for sheet in workbook.worksheets:
            top_left = list(sheet.iter_rows(
                min_row=1, max_row=min(3, sheet.max_row),
                min_col=1, max_col=min(8, sheet.max_column),
                values_only=True,
            ))
            sheets.append({
                "name": sheet.title,
                "rows": sheet.max_row,
                "columns": sheet.max_column,
                "preview": [
                    [None if cell is None else str(cell)[:100] for cell in row]
                    for row in top_left
                ],
            })
    finally:
        workbook.close()
    return {"sheets": sheets}


def preview_zip(path: Path) -> dict:
    with ZipFile(path) as archive:
        names = archive.namelist()
        return {
            "member_count": len(names),
            "first_members": names[:30],
            "total_uncompressed_bytes": sum(item.file_size for item in archive.infolist()),
        }


def preview_biom(path: Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        prefix = handle.read(LIMIT)
    content = prefix.decode("utf-8", errors="replace")
    markers = ["matrix_type", "shape", "rows", "columns", "data", "generated_by"]
    return {
        "prefix_excerpt": content[:500],
        "found_json_keys_in_prefix": [key for key in markers if f'"{key}"' in content],
        "note": "BIOM metadata/shape require proper BIOM JSON/HDF5 parsing; no dimensions inferred from prefix.",
    }


def describe(path: Path) -> dict:
    name = path.name.lower()
    try:
        if name.endswith((".tsv", ".tsv.gz", ".csv", ".csv.gz", ".txt", ".txt.gz")):
            return {"kind": "delimited_matrix_or_metadata", **preview_delimited(path)}
        if name.endswith((".xlsx", ".xlsm")):
            return {"kind": "excel", **preview_excel(path)}
        if name.endswith(".zip"):
            return {"kind": "archive", **preview_zip(path)}
        if name.endswith((".biom", ".biom.gz")):
            return {"kind": "biom", **preview_biom(path)}
        if name.endswith((".h5", ".h5mu", ".h5ad")):
            return {"kind": "hierarchical_binary", "note": "Use h5py/MuData for semantic dimensions."}
        if name.endswith(".rds"):
            return {"kind": "r_object", "note": "Read with R to inspect SummarizedExperiment or assay metadata."}
        return {"kind": "other"}
    except Exception as exc:
        return {"kind": "unreadable", "error": f"{type(exc).__name__}: {exc}"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=os.environ.get("PULSE_OMICS_DATA_ROOT", DEFAULT_ROOT))
    parser.add_argument("--datasets", nargs="+",
                        choices=("ipop", "hmp2", "covaar", "p100", "csiro"),
                        default=("ipop", "hmp2", "covaar", "p100", "csiro"))
    args = parser.parse_args()
    base = Path(args.data_root).expanduser() / "raw"
    for cohort in args.datasets:
        folder = base / cohort
        if not folder.exists():
            print(f"{cohort}: no files yet")
            continue
        files = []
        for path in sorted(folder.rglob("*")):
            if not path.is_file() or path.name in IGNORE or path.name.endswith(".part"):
                continue
            record = {
                "relative_path": path.relative_to(folder).as_posix(),
                "bytes": path.stat().st_size,
                **describe(path),
            }
            files.append(record)
        out = folder / "data_inventory.json"
        out.write_text(json.dumps({
            "dataset": cohort, "files": files,
            "semantic_join_status": "unverified_pending_cohort_specific_metadata_mapping",
        }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        errors = sum(record["kind"] == "unreadable" for record in files)
        print(f"{cohort}: {len(files)} files, {errors} unreadable, report={out}")


if __name__ == "__main__":
    main()
