#!/usr/bin/env python3
"""Inspect Su 2020 subject and assay IDs without making row-order joins.

Writes a local-only schema report under F: (never into GitHub).
"""
from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from pathlib import Path

from openpyxl import load_workbook

SHEETS = {
    "clinical": "S1.1 Patient Clinical Data",
    "protein": "S1.3 Plasma Proteomic Data",
    "metab": "S1.4 Plasma Metabolomic Data",
}
DATA_ROOT = "/mnt/f/datasets/pulse-omics"


def norm(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip().casefold()
    if not text or text in {"nan", "none", "na", "n/a", "-"}:
        return None
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return text


def inspect_sheet(sheet, max_preview_rows: int = 12, max_preview_cols: int = 24) -> dict:
    """Read only the first few columns of each sheet for ID inference."""
    id_scan_cols = sheet.max_column if "Patient Clinical Data" in sheet.title else min(sheet.max_column, 2)
    max_rows = min(sheet.max_row, 600)
    rows = list(sheet.iter_rows(
        min_row=1, max_row=max_rows, max_col=id_scan_cols, values_only=True
    ))
    candidates = []
    for header_row in range(min(1, len(rows) - 1)):
        for col in range(id_scan_cols):
            entries = [norm(row[col]) for row in rows[header_row + 1:]]
            nonempty = [value for value in entries if value is not None]
            distinct = set(nonempty)
            if len(nonempty) < 15 or len(distinct) < 15:
                continue
            is_assay_patient_label = (
                col == 0
                and str(rows[header_row][col]).strip() == "Healthy donor sample or COVID19 sample"
            )
            if len(distinct) / len(nonempty) < 0.85 and not is_assay_patient_label:
                continue
            header = rows[header_row][col]
            candidates.append({
                "header_row_1based": header_row + 1,
                "id_column_1based": col + 1,
                "label": "" if header is None else str(header)[:100],
                "nonempty": len(nonempty),
                "distinct": len(distinct),
                "_values": distinct,
            })
    preview = []
    for row in sheet.iter_rows(
        min_row=1, max_row=min(max_preview_rows, sheet.max_row),
        min_col=1, max_col=min(max_preview_cols, sheet.max_column),
        values_only=True,
    ):
        preview.append([None if value is None else str(value)[:100] for value in row])
    return {
        "sheet": sheet.title,
        "max_row": sheet.max_row,
        "max_column": sheet.max_column,
        "preview": preview,
        "candidates": candidates,
        "first_columns": rows,
    }


def overlap_candidates(left: dict, right: dict, limit: int = 12) -> list[dict]:
    scores = []
    for a in left["candidates"]:
        for b in right["candidates"]:
            overlap = len(a["_values"] & b["_values"])
            if overlap < 12:
                continue
            relative = overlap / min(len(a["_values"]), len(b["_values"]))
            scores.append({
                "left": {
                    "header_row": a["header_row_1based"],
                    "column": a["id_column_1based"],
                    "label": a["label"],
                },
                "right": {
                    "header_row": b["header_row_1based"],
                    "column": b["id_column_1based"],
                    "label": b["label"],
                },
                "common_ids": overlap,
                "min_side_coverage": round(relative, 4),
            })
    return sorted(
        scores,
        key=lambda entry: (entry["common_ids"], entry["min_side_coverage"]),
        reverse=True,
    )[:limit]


def clinical_group_candidates(info: dict) -> list[dict]:
    """Look for recurring patient identifiers, separately from unique sample IDs."""
    rows = info["first_columns"]
    choices = []
    for header_row in range(min(1, len(rows) - 1)):
        for col in range(len(rows[0])):
            label = "" if rows[header_row][col] is None else str(rows[header_row][col])
            if label.strip().casefold() not in {"study subject id", "patient_id"}:
                continue
            values = [norm(row[col]) for row in rows[header_row + 1:]]
            counts = Counter(value for value in values if value is not None)
            n_repeated = sum(number > 1 for number in counts.values())
            if len(counts) >= 20 and n_repeated >= 10:
                choices.append({
                    "header_row": header_row + 1,
                    "column": col + 1,
                    "label": label[:100],
                    "distinct_ids": len(counts),
                    "ids_with_multiple_rows": n_repeated,
                    "maximum_rows_per_id": max(counts.values()),
                })
    return sorted(
        choices, key=lambda c: (c["ids_with_multiple_rows"], -c["maximum_rows_per_id"]),
        reverse=True,
    )[:12]


def build_report(workbook: Path) -> dict:
    excel = load_workbook(workbook, read_only=True, data_only=True)
    try:
        missing = [name for name in SHEETS.values() if name not in excel.sheetnames]
        if missing:
            raise ValueError(f"Missing expected Su2020 sheets: {missing}")
        analyzed = {
            alias: inspect_sheet(excel[name])
            for alias, name in SHEETS.items()
        }
    finally:
        excel.close()
    result = {
        "source": str(workbook),
        "sheets": {},
        "cross_sheet_overlap": {},
        "clinical_repeated_identifier_candidates": clinical_group_candidates(analyzed["clinical"]),
        "note": "Candidate overlaps are diagnostics, not verified patient/visit links. No training data has been exported.",
    }
    for alias, info in analyzed.items():
        result["sheets"][alias] = {
            "sheet": info["sheet"],
            "dimensions": [info["max_row"], info["max_column"]],
            "preview": info["preview"],
            "candidate_id_columns": [
                {key: value for key, value in c.items() if key != "_values"}
                for c in info["candidates"]
            ],
        }
    for left, right in (("clinical", "protein"), ("clinical", "metab"), ("protein", "metab")):
        key = f"{left}_to_{right}"
        result["cross_sheet_overlap"][key] = overlap_candidates(analyzed[left], analyzed[right])
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=os.environ.get("PULSE_OMICS_DATA_ROOT", DATA_ROOT))
    args = parser.parse_args()
    root = Path(args.data_root).expanduser()
    workbook = root / "raw" / "su2020" / "table_s1.xlsx"
    if not workbook.is_file():
        parser.error(f"Workbook missing: {workbook}")
    report = build_report(workbook)
    out = root / "processed" / "su2020" / "schema_profile.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print("Su2020 patient/omics schema profile")
    for alias, info in report["sheets"].items():
        print(f"  {alias}: {info['dimensions'][0]} rows x {info['dimensions'][1]} cols")
    for pair, matches in report["cross_sheet_overlap"].items():
        if matches:
            best = matches[0]
            print(f"  {pair}: {best['common_ids']} shared candidate sample IDs "
                  f"({best['left']['label']!r} versus {best['right']['label']!r})")
        else:
            print(f"  {pair}: no unique-ID overlap established")
    print(f"  Repeated patient-ID candidates: {len(report['clinical_repeated_identifier_candidates'])}")
    print(f"Saved locally: {out}")
    print("No rows were joined and no longitudinal samples were inferred from row position.")


if __name__ == "__main__":
    main()
