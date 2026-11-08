"""Extract Su 2020 plasma assays using verified worksheet headers and ID joins.

This module produces raw, unscaled assay matrices. Missing feature values remain NaN.
A numeric assay sample identifier must never be assumed to be a patient ID.
"""
from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

SHEETS = {
    "clinical": "S1.1 Patient Clinical Data",
    "protein": "S1.3 Plasma Proteomic Data",
    "metab": "S1.4 Plasma Metabolomic Data",
}
CLINICAL_COLUMNS = ("Study Subject ID", "Sample ID", "Blood draw time point")
ASSAY_COLUMNS = ("Healthy donor sample or COVID19 sample", "sample_id")


def as_id(value: object) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    if not text:
        return None
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return None if text.casefold() in ("na", "n/a", "none", "nan") else text


def _require_unique_headers(frame: pd.DataFrame, sheet: str) -> None:
    if frame.columns.duplicated().any():
        raise ValueError(f"{sheet}: duplicate worksheet headers")


def _clinical(frame: pd.DataFrame) -> pd.DataFrame:
    _require_unique_headers(frame, "clinical")
    missing = set(CLINICAL_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f"Clinical sheet is missing headers: {sorted(missing)}")
    df = frame.copy()
    for col in ("Study Subject ID", "Sample ID"):
        df[col] = df[col].map(as_id)
        if df[col].isna().any():
            raise ValueError(f"Clinical sheet has missing {col}")
    if df["Sample ID"].duplicated().any():
        raise ValueError("Clinical sheet has duplicate Sample ID")
    df["visit_order"] = df["Blood draw time point"].map(
        lambda x: {"T1": 1, "T2": 2}.get(str(x).strip().upper())
    )
    if df["visit_order"].isna().any():
        raise ValueError("Unexpected clinical timepoint (expected T1/T2)")
    if df.duplicated(["Study Subject ID", "visit_order"]).any():
        raise ValueError("A clinical patient appears more than once at one timepoint")
    return df


def _assay(frame: pd.DataFrame, name: str) -> tuple[pd.DataFrame, list[str]]:
    _require_unique_headers(frame, name)
    missing = set(ASSAY_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f"{name} sheet is missing headers: {sorted(missing)}")
    df = frame.copy()
    df["sample_id"] = df["sample_id"].map(as_id)
    if df["sample_id"].isna().any() or df["sample_id"].duplicated().any():
        raise ValueError(f"{name}: sample_id must be present and unique")
    features = [name for name in df.columns if name not in ASSAY_COLUMNS]
    if not features:
        raise ValueError(f"{name}: no measured features")
    numeric = df.loc[:, features].apply(pd.to_numeric, errors="raise")
    if np.isinf(numeric.to_numpy(dtype=np.float64)).any():
        raise ValueError(f"{name}: infinite assay abundance")
    # One concat prevents fragmentation in wide metabolomics tables.
    df = pd.concat([df.loc[:, list(ASSAY_COLUMNS)], numeric], axis=1)
    df = df.set_index("sample_id")
    return df, features


def _bridge_candidates(
    clinical: pd.DataFrame,
    assay: pd.DataFrame,
) -> list[dict]:
    assay_ids = set(assay.index)
    candidates = []
    for column in clinical.columns:
        if column in (*CLINICAL_COLUMNS, "visit_order"):
            continue
        ids = clinical[column].map(as_id)
        matched = ids.isin(assay_ids)
        n_matches = int(matched.sum())
        if n_matches < 10:
            continue
        vals = ids[matched].dropna()
        if vals.duplicated().any():
            continue
        candidates.append({
            "column": str(column),
            "matched_clinical_rows": n_matches,
            "unique_assay_ids": int(vals.nunique()),
            "fraction_of_clinical_rows": round(n_matches / len(clinical), 5),
        })
    return sorted(candidates, key=lambda x: x["matched_clinical_rows"], reverse=True)


def _choose_bridge(candidates: list[dict], modality: str) -> str | None:
    """Allow only explicitly named clinical-to-assay ID fields as a fallback."""
    token = "prot" if modality == "protein" else "metab"
    named = [
        candidate for candidate in candidates
        if token in candidate["column"].casefold()
        and "sample" in candidate["column"].casefold()
    ]
    if not named or named[0]["matched_clinical_rows"] < 10:
        return None
    ranked = sorted(named, key=lambda item: item["matched_clinical_rows"], reverse=True)
    if len(ranked) > 1 and ranked[0]["matched_clinical_rows"] == ranked[1]["matched_clinical_rows"]:
        return None
    return ranked[0]["column"]


def normalized_patient_sample(value: object) -> str | None:
    """Use the documented INCOV BL/AC -> T1/T2 naming convention.

    Reproduced from ByronSyun/DIVAS_COVID19_CaseStudy's process_bulk scripts.
    This function NEVER converts numeric assay/donor IDs to patient IDs.
    """
    value = as_id(value)
    if value is None:
        return None
    match = re.fullmatch(
        r"(?:INCOV|COVID[_-])0*([0-9]+)[_-](BL|AC|T1|T2|1|2)",
        value, flags=re.IGNORECASE,
    )
    if not match:
        return None
    subject = int(match.group(1))
    time = match.group(2).upper()
    visit = 1 if time in {"1", "BL", "T1"} else 2
    return f"INCOV{subject:03d}-{visit}"


def _patient_sample_labels(
    clinical: pd.DataFrame, assay: pd.DataFrame, name: str
) -> tuple[dict[str, str], dict[str, int]]:
    """Match clinical samples against BOTH assay label and assay sample_id.

    Earlier implementation checked only the assay group/label column, missing
    the actual sample_id field and the INCOVxxx-BL/AC suffix variants.
    """
    clinical_samples = {
        normalized_patient_sample(value) or str(value).casefold(): value
        for value in clinical["Sample ID"]
    }
    if len(clinical_samples) != len(clinical):
        raise ValueError("Ambiguous clinical Sample ID normalization")
    matched = {}
    sources = {"assay_sample_id": 0, "assay_label": 0}
    for assay_id, label in assay[ASSAY_COLUMNS[0]].items():
        hits = set()
        for candidate, source in (
            (assay_id, "assay_sample_id"),
            (label, "assay_label"),
        ):
            raw = as_id(candidate)
            key = normalized_patient_sample(candidate) or (
                raw.casefold() if raw is not None else None
            )
            if key in clinical_samples:
                hits.add((clinical_samples[key], source))
        samples = {sample for sample, _ in hits}
        if len(samples) > 1:
            raise ValueError(f"{name}: conflicting sample identities for assay row {assay_id!r}")
        if not samples:
            continue
        clinical_sample = samples.pop()
        if clinical_sample in matched:
            raise ValueError(f"{name}: two assay rows map to {clinical_sample!r}")
        matched[clinical_sample] = assay_id
        for _, source in hits:
            sources[source] += 1
    return matched, sources


def _assay_id_formats(assay: pd.DataFrame) -> dict[str, int]:
    counts = {"incov": 0, "covid": 0, "numeric": 0, "other": 0}
    for sample_id in assay.index:
        name = str(sample_id).upper()
        if name.startswith("INCOV"):
            counts["incov"] += 1
        elif name.startswith("COVID_"):
            counts["covid"] += 1
        elif name.isdecimal():
            counts["numeric"] += 1
        else:
            counts["other"] += 1
    return counts


def _atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    try:
        frame.to_csv(tmp, index=False, na_rep="")
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def prepare(workbook: Path, output_dir: Path) -> dict:
    if not workbook.is_file():
        raise FileNotFoundError(workbook)
    tables = pd.read_excel(workbook, sheet_name=list(SHEETS.values()), engine="openpyxl")
    clinical = _clinical(tables[SHEETS["clinical"]])
    protein, protein_features = _assay(tables[SHEETS["protein"]], "protein")
    metab, metab_features = _assay(tables[SHEETS["metab"]], "metab")

    # Always preserve both full published assay matrices independently, including
    # healthy donors and unmatched samples. No imputation or scaling is performed.
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, data, features in (
        ("protein", protein, protein_features),
        ("metab", metab, metab_features),
    ):
        original = data.reset_index().loc[:, ["sample_id", *features]]
        _atomic_csv(original, output_dir / f"{name}_raw_by_assay_id.csv")

    p_labels, p_sources = _patient_sample_labels(clinical, protein, "protein")
    m_labels, m_sources = _patient_sample_labels(clinical, metab, "metab")
    p_candidates = _bridge_candidates(clinical, protein)
    m_candidates = _bridge_candidates(clinical, metab)
    p_col = _choose_bridge(p_candidates, "protein")
    m_col = _choose_bridge(m_candidates, "metab")

    # The same published sample_id can occur in both assays; that's evidence for
    # a cross-sectional matched specimen, not evidence for a longitudinal patient.
    shared_assay_ids = set(protein.index) & set(metab.index)
    counts = Counter(clinical["Study Subject ID"])
    report = {
        "source": str(workbook),
        "clinical_rows": int(len(clinical)),
        "clinical_patients": int(len(counts)),
        "clinical_patients_with_two_draws": sum(n == 2 for n in counts.values()),
        "assay_rows": {"protein": int(len(protein)), "metab": int(len(metab))},
        "assay_id_formats": {
            "protein": _assay_id_formats(protein),
            "metab": _assay_id_formats(metab),
        },
        "assay_features": {"protein": len(protein_features), "metab": len(metab_features)},
        "assay_shared_sample_ids": len(shared_assay_ids),
        "assay_missing_feature_fractions": {
            "protein": round(float(protein[protein_features].isna().to_numpy().mean()), 5),
            "metab": round(float(metab[metab_features].isna().to_numpy().mean()), 5),
        },
        "clinical_sample_matches": {"protein": len(p_labels), "metab": len(m_labels)},
        "exact_clinical_label_matches": {
            "protein": p_sources["assay_label"],
            "metab": m_sources["assay_label"],
        },
        "match_sources": {"protein": p_sources, "metab": m_sources},
        "clinical_bridge_candidates": {"protein": p_candidates, "metab": m_candidates},
        "selected_bridge_columns": {"protein": p_col, "metab": m_col},
        "status": "raw_assays_exported_no_verified_patient_join",
    }

    if p_labels or m_labels or p_col or m_col:
        linked = clinical.loc[:, [*CLINICAL_COLUMNS, "visit_order"]].copy()
        linked["protein_assay_id"] = pd.Series(
            [p_labels.get(value) for value in clinical["Sample ID"]],
            index=clinical.index, dtype="object",
        )
        linked["metab_assay_id"] = pd.Series(
            [m_labels.get(value) for value in clinical["Sample ID"]],
            index=clinical.index, dtype="object",
        )
        for name, bridge_col in (("protein", p_col), ("metab", m_col)):
            if bridge_col:
                target = f"{name}_assay_id"
                fallback = clinical[bridge_col].map(as_id)
                conflicting = (
                    linked[target].notna()
                    & fallback.notna()
                    & (linked[target] != fallback)
                )
                if conflicting.any():
                    raise ValueError(
                        f"{name}: clinical bridge disagrees with explicit COVID19 label"
                    )
                unset = linked[target].isna() & fallback.notna()
                linked.loc[unset, target] = fallback.loc[unset].astype(object)
        for name, assay in (("protein", protein), ("metab", metab)):
            col = f"{name}_assay_id"
            linked[col] = linked[col].where(linked[col].isin(assay.index), None)

        linked = linked.loc[
            linked["protein_assay_id"].notna() | linked["metab_assay_id"].notna()
        ].copy()
        for col in ("protein_assay_id", "metab_assay_id"):
            ids = linked[col].dropna()
            if ids.duplicated().any():
                raise ValueError(f"Clinical mapping is not one-to-one: {col}")

        linked = linked.sort_values(["Study Subject ID", "visit_order"])
        subjects = pd.DataFrame({
            "sample_id": linked["Sample ID"],
            "patient_id": linked["Study Subject ID"],
            "visit_order": linked["visit_order"].astype(int),
        })
        if not subjects.empty:
            _atomic_csv(subjects, output_dir / "subjects.csv")
            for name, assay, features in (
                ("proteomics", protein, protein_features),
                ("metabolomics", metab, metab_features),
            ):
                col = "protein_assay_id" if name == "proteomics" else "metab_assay_id"
                selected = linked.loc[linked[col].notna(), ["Sample ID", col]]
                values = assay.loc[selected[col].tolist(), features].copy()
                values.insert(0, "sample_id", selected["Sample ID"].tolist())
                _atomic_csv(values, output_dir / f"{name}.csv")
            _atomic_csv(
                linked.loc[:, ["Sample ID", "Study Subject ID", "visit_order", "protein_assay_id", "metab_assay_id"]],
                output_dir / "assay_mapping.csv",
            )
            per_patient = Counter(linked["Study Subject ID"])
            paired = linked["protein_assay_id"].notna() & linked["metab_assay_id"].notna()
            report["linked_clinical_draws"] = len(linked)
            report["linked_paired_draws"] = int(paired.sum())
            report["linked_patients_with_two_draws"] = sum(n == 2 for n in per_patient.values())
            report["status"] = "canonical_raw_omics_exported"

    output = output_dir / "preparation_report.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report
