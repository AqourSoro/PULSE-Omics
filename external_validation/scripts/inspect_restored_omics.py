#!/usr/bin/env python3
"""Read the restored iPOP seasonal RData and P100 publisher attachments.

Outputs native field/feature catalogues and aggregate statistics. Participant
keys are written only when --private-output-dir is explicitly supplied.
No assay matrices are exported and no identifiers are inferred from row order.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timedelta
import hashlib
import io
import json
from pathlib import Path
import warnings
from zipfile import ZipFile

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")


def digest(path: Path) -> dict:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(part)
    return {"file": path.name, "bytes": path.stat().st_size, "sha256": h.hexdigest()}


def read_ipop(path: Path) -> dict:
    import rdata

    def factor(obj, attrs):
        # rdata 1.1.0's default factor constructor mishandles masked NA codes
        # in this release. Preserve R factor levels and NA explicitly.
        levels = list(attrs["levels"])
        values = [None if np.ma.is_masked(i) or not np.isfinite(i) or i < 1
                  else levels[int(i) - 1] for i in obj]
        return pd.Categorical(values, categories=levels)

    constructors = dict(rdata.conversion.DEFAULT_CLASS_MAP)
    constructors["factor"] = factor
    # Keep serialized POSIXct seconds intact; do not silently repair the years.
    constructors["POSIXct"] = lambda obj, attrs: obj
    # R's serialization includes harmless attributes unsupported by rdata.
    # Do not turn these warnings into fatal parse errors; validate objects below.
    with warnings.catch_warnings(record=True) as received:
        warnings.simplefilter("always")
        objects = rdata.read_rda(path, constructor_dict=constructors)
    if not isinstance(objects, dict) or not objects:
        raise ValueError("No named objects decoded from published seasonal RData")
    return objects


def inspect_ipop(path: Path, output: Path, private: Path | None) -> dict:
    objects = read_ipop(path)
    summary = {**digest(path), "objects": [], "join_basis": ["SubjectID", "Date"]}
    fields = []
    for name, frame in objects.items():
        if not isinstance(frame, pd.DataFrame):
            raise TypeError(f"Unexpected native object: {name}: {type(frame)}")
        metadata = "annotation_colData" in name
        info = {"object": name, "rows": len(frame), "columns": len(frame.columns),
                "row_index_unique": bool(frame.index.is_unique),
                "row_index_examples": [str(i) for i in frame.index[:3]],
                "missing_cells": int(frame.isna().sum().sum())}
        if metadata:
            info.update({
                "subjects_with_nonmissing_ID": int(frame.SubjectID.nunique()),
                "missing_subject_rows": int(frame.SubjectID.isna().sum()),
                "duplicate_subject_date_rows": int(frame.duplicated(["SubjectID", "Date"]).sum()),
                "missing_by_field": {str(k): int(v) for k, v in frame.isna().sum().items()},
            })
        else:
            a = frame.to_numpy(dtype=float)
            info.update({"finite_cells": int(np.isfinite(a).sum()),
                         "constant_features": int((frame.nunique() <= 1).sum()),
                         "value_quantiles_0_1_50_99_100_percent": np.nanquantile(
                             a, [0, .01, .5, .99, 1]).tolist()})
        summary["objects"].append(info)
        for pos, (field, values) in enumerate(frame.items(), 1):
            record = {"file": path.name, "object": name, "column_1based": pos,
                      "raw_name": str(field), "role": "metadata" if metadata else "measurement",
                      "dtype": str(values.dtype), "rows": len(values),
                      "missing": int(values.isna().sum()), "unit": "",
                      "transform": "not specified in serialized object"}
            if metadata:
                if field == "Date":
                    record.update(unit="seconds since 1970-01-01; R POSIXct", transform="raw serialized value")
                elif field == "Time":
                    record.update(unit="day of year", transform="seasonal index; not elapsed time")
                elif field == "BMI":
                    record.update(unit="kg/m^2", transform="stored as R factor")
                elif field == "IRIS":
                    record["transform"] = "R factor; IR / IS; NA retained"
            fields.append(record)
    pd.DataFrame(fields).to_csv(output / "ipop_native_fields.csv.gz", index=False,
                                compression={"method": "gzip", "mtime": 0})

    p0 = objects["Proteomics_annotation_colData"].copy()
    m0 = objects["Metabolomics_annotation_colData"].copy()
    p0["P_metadata_row_key"] = p0.index.astype(str)
    m0["M_metadata_row_key"] = m0.index.astype(str)
    p = p0[p0.SubjectID.notna()].copy()
    m = m0[m0.SubjectID.notna()].copy()
    links = p.merge(m, on=["SubjectID", "Date"], how="outer", suffixes=("_P", "_M"),
                    validate="one_to_one", indicator=True)
    both = links[links["_merge"] == "both"]
    counts = both.groupby("SubjectID", observed=True).size()
    gaps = (both.sort_values(["SubjectID", "Date"])
            .groupby("SubjectID", observed=True).Date.diff().dropna() / 86400)
    summary["protein_metabolite_metadata_intersection"] = {
        "paired_subject_date_records": len(both), "paired_subjects": len(counts),
        "subjects_at_least_two_paired_dates": int((counts >= 2).sum()),
        "subjects_at_least_three_paired_dates": int((counts >= 3).sum()),
        "paired_date_count_distribution": {str(k): int(v) for k, v in counts.value_counts().sort_index().items()},
        "protein_only_valid_key_records": int((links["_merge"] == "left_only").sum()),
        "metabolite_only_valid_key_records": int((links["_merge"] == "right_only").sum()),
        "protein_missing_subject_records_excluded": int(p0.SubjectID.isna().sum()),
        "metabolite_missing_subject_records_excluded": int(m0.SubjectID.isna().sum()),
        "gap_count": len(gaps),
        "gap_days_min_q25_median_q75_max": np.quantile(gaps, [0, .25, .5, .75, 1]).tolist(),
        "specimen_identity": "subject/date metadata match; no native specimen ID in this release",
        "matrix_row_identity": "matrix row labels and metadata row labels differ; no positional measurement join exported",
    }
    summary["native_date_encoding"] = {}
    for name, frame in (("protein", p0), ("metabolite", m0)):
        dates = [datetime(1970, 1, 1) + timedelta(seconds=float(v)) for v in frame.Date]
        summary["native_date_encoding"][name] = {
            "raw_seconds_min": float(frame.Date.min()), "raw_seconds_max": float(frame.Date.max()),
            "decoded_UTC_min": str(min(dates)), "decoded_UTC_max": str(max(dates)),
            "time_field_disagrees_with_day_of_year": sum(d.timetuple().tm_yday != t for d, t in zip(dates, frame.Time)),
            "calendar_status": "years outside study period; retained without correction",
        }
    if private is not None:
        private.mkdir(parents=True, exist_ok=True)
        links.to_csv(private / "ipop_subject_date_link_keys.csv", index=False)
        p0[p0.SubjectID.isna()].to_csv(private / "ipop_protein_missing_subject.csv", index=True)
    dump(output / "ipop_structure.json", summary)
    return summary


def inspect_p100(folder: Path, output: Path) -> dict:
    summary = {"files": [], "worksheets": [], "zip_members": []}
    fields = []
    for path in sorted(folder.glob("41587_2017_BFnbt3870_MOESM*")):
        summary["files"].append(digest(path))
        if path.suffix == ".xlsx":
            workbook = load_workbook(path, read_only=True, data_only=True)
            try:
                for sheet in workbook:
                    rows = list(sheet.values)
                    nonempty = [i for i, row in enumerate(rows, 1) if any(v is not None for v in row)]
                    legend = sheet.title == "Legend"
                    header = list(rows[0])
                    info = {"file": path.name, "sheet": sheet.title,
                            "worksheet_rows": sheet.max_row, "worksheet_columns": sheet.max_column,
                            "nonempty_rows": len(nonempty), "first_nonempty_row": min(nonempty),
                            "last_nonempty_row": max(nonempty),
                            "header_row": None if legend else 1,
                            "data_rows": None if legend else len(nonempty) - 1,
                            "raw_headers": None if legend else header}
                    if legend:
                        # Legends contain field descriptions, not patient records.
                        info["nonempty_cells_by_row"] = [{"row": i, "cells": list(rows[i-1])} for i in nonempty]
                    summary["worksheets"].append(info)
                    if legend:
                        continue
                    for j, raw_name in enumerate(header):
                        values = [row[j] for row in rows[1:]]
                        fields.append({"file": path.name, "sheet": sheet.title,
                                       "column_1based": j+1, "excel_column": get_column_letter(j+1),
                                       "raw_header": raw_name, "header_row": 1,
                                       "data_rows": len(values), "empty_cells": values.count(None),
                                       "native_type_counts": json.dumps(dict(Counter(type(v).__name__ for v in values))),
                                       "literal_missing_tokens": json.dumps(dict(Counter(v for v in values if isinstance(v, str) and v in {"NA", "N/A", "-", "No response", ""})))})
                    if path.name.endswith("MOESM30_ESM.xlsx") and sheet.title in {
                        "Metabolomics", "Proteomics", "Clinical Lab Tests", "Genetic Traits", "Microbiome"
                    }:
                        # These are analyte annotations, not participant measurements.
                        name = "p100_" + sheet.title.lower().replace(" ", "_") + "_annotations.csv"
                        if sheet.title == "Microbiome":
                            pd.DataFrame(rows[1:], columns=header).to_csv(
                                output / (name + ".gz"), index=False,
                                compression={"method": "gzip", "mtime": 0})
                        else:
                            pd.DataFrame(rows[1:], columns=header).to_csv(output / name, index=False)
                    if path.name.endswith("MOESM30_ESM.xlsx") and sheet.title == "Demographics":
                        demo = pd.DataFrame(rows[1:], columns=header)
                        summary["availability_count_fields"] = {}
                        for col in demo.columns:
                            if col.startswith(("clinical.", "metabolomics.", "proteomics.", "microbiome.")):
                                values = pd.to_numeric(demo[col], errors="raise")
                                summary["availability_count_fields"][col] = {
                                    "nonmissing_subjects": int(values.notna().sum()),
                                    "positive_count_subjects": int((values > 0).sum()),
                                    "min_feature_count": float(values.min()),
                                    "max_feature_count": float(values.max()),
                                }
                        summary["paired_round_availability"] = {}
                        for platform in ("olink", "srm_mass_spec", "either"):
                            rounds = []
                            for timepoint in (1, 2, 3):
                                metab = pd.to_numeric(demo[f"metabolomics.metabolon.{timepoint}"], errors="raise") > 0
                                olink = pd.to_numeric(demo[f"proteomics.olink.{timepoint}"], errors="raise") > 0
                                srm = pd.to_numeric(demo[f"proteomics.srm_mass_spec.{timepoint}"], errors="raise") > 0
                                protein = olink if platform == "olink" else srm if platform == "srm_mass_spec" else olink | srm
                                rounds.append(metab & protein)
                            count = pd.concat(rounds, axis=1).sum(axis=1)
                            summary["paired_round_availability"][platform] = {
                                "round_1_2_3": [int(v.sum()) for v in rounds],
                                "subjects_at_least_two_rounds": int((count >= 2).sum()),
                                "subjects_three_rounds": int((count >= 3).sum()),
                                "basis": "positive per-participant analyte counts; no measurement matrix or specimen identifiers supplied",
                            }
            finally:
                workbook.close()
        elif path.suffix == ".zip":
            with ZipFile(path) as archive:
                for item in archive.infolist():
                    record = {"file": path.name, "member": item.filename, "bytes": item.file_size}
                    if item.filename.endswith(".cys"):
                        with ZipFile(io.BytesIO(archive.read(item))) as cys:
                            record["nested_members"] = cys.namelist()
                            record["content_type"] = "Cytoscape session / correlation network"
                    summary["zip_members"].append(record)
    pd.DataFrame(fields).to_csv(output / "p100_native_fields.csv", index=False)
    dump(output / "p100_structure.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    root = Path("/mnt/f/datasets/pulse-omics")
    parser.add_argument("--ipop", type=Path, default=root / "raw/ipop/seasonal/Multi_Omics_Seasonal.RData")
    parser.add_argument("--p100-dir", type=Path, default=root / "raw/p100")
    parser.add_argument("--output-dir", type=Path, default=root / "processed/restored_sources/structure")
    parser.add_argument("--private-output-dir", type=Path)
    args = parser.parse_args()
    if not args.ipop.is_file() or not list(args.p100_dir.glob("*MOESM*.xlsx")):
        parser.error("Restore iPOP and P100 files first; no missing-file counts are reported as zero")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ipop = inspect_ipop(args.ipop, args.output_dir, args.private_output_dir)
    p100 = inspect_p100(args.p100_dir, args.output_dir)
    print(f"iPOP: {len(ipop['objects'])} native objects; "
          f"{ipop['protein_metabolite_metadata_intersection']['paired_subject_date_records']} matched metadata keys")
    print(f"P100: {len(p100['files'])} attachments, {len(p100['worksheets'])} worksheets")
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    main()
