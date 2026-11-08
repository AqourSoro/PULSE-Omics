"""Validated omics CSVs to the existing PULSE LongitudinalDataset interface.

Expected input files:
- subjects.csv: sample_id, patient_id, visit_order (one row per sample)
- proteomics.csv: sample_id, protein_1, ... (sample-by-feature)
- metabolomics.csv: sample_id, metabolite_1, ... (sample-by-feature)

The source-specific Excel extraction must establish these identifiers first.
This adapter does not invent visits, fill partially missing features, or
normalize using held-out patients.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pulse.data.dataset import LongitudinalDataset


def _load_assay(csv_path: str | Path, name: str) -> tuple[pd.DataFrame, list[str]]:
    table = pd.read_csv(csv_path, dtype={"sample_id": str})
    if "sample_id" not in table.columns:
        raise ValueError(f"{name}: missing sample_id column")
    if table["sample_id"].isna().any() or table["sample_id"].duplicated().any():
        raise ValueError(f"{name}: sample_id values must be present and unique")
    feature_names = [col for col in table.columns if col != "sample_id"]
    if not feature_names or len(feature_names) != len(set(feature_names)):
        raise ValueError(f"{name}: empty or duplicate feature names")
    values = table[feature_names].apply(pd.to_numeric, errors="raise")
    values.index = table["sample_id"].astype(str)
    return values.astype(np.float32), feature_names


def load_sequences(
    subjects_path: str | Path,
    proteomics_path: str | Path,
    metabolomics_path: str | Path,
) -> tuple[list[str], list[list[dict]], list[list[dict]], dict[str, list[str]]]:
    """Return patient IDs, visit dictionaries, 0/2 masks, and feature names.

    0 = available; 2 = naturally missing. A sample missing in an assay file is
    treated as an absent *whole assay*, not as a zero vector.
    """
    subjects = pd.read_csv(
        subjects_path, dtype={"sample_id": str, "patient_id": str}
    )
    required = {"sample_id", "patient_id", "visit_order"}
    if not required.issubset(subjects.columns):
        raise ValueError(f"subjects: missing {sorted(required - set(subjects.columns))}")
    if subjects[list(required)].isna().any().any():
        raise ValueError("subjects: sample, patient and visit must be non-null")
    if subjects["sample_id"].duplicated().any():
        raise ValueError("subjects: duplicate sample_id")
    subjects["visit_order"] = pd.to_numeric(subjects["visit_order"], errors="raise")
    if subjects["visit_order"].isna().any():
        raise ValueError("subjects: invalid visit_order")
    if subjects.duplicated(["patient_id", "visit_order"]).any():
        raise ValueError("subjects: multiple samples assigned to a patient/visit")

    protein, protein_features = _load_assay(proteomics_path, "proteomics")
    metab, metab_features = _load_assay(metabolomics_path, "metabolomics")
    subject_ids = set(subjects["sample_id"])
    for name, assay in (("proteomics", protein), ("metabolomics", metab)):
        unknown = set(assay.index) - subject_ids
        if unknown:
            raise ValueError(f"{name}: {len(unknown)} sample IDs absent from subjects")
    patient_ids = []
    visits_data = []
    missing_masks = []
    for patient_id, group in subjects.sort_values(
        ["patient_id", "visit_order"]
    ).groupby("patient_id", sort=True):
        visits = []
        masks = []
        for sample_id in group["sample_id"]:
            data = {}
            mask = {}
            for name, assay in (("protein", protein), ("metab", metab)):
                if sample_id not in assay.index:
                    data[name], mask[name] = None, 2
                    continue
                values = assay.loc[sample_id].to_numpy(dtype=np.float32)
                if np.isnan(values).all():
                    data[name], mask[name] = None, 2
                elif not np.isfinite(values).all():
                    raise ValueError(
                        f"{name}: sample {sample_id} contains partially missing "
                        "or non-finite features; define preprocessing explicitly"
                    )
                else:
                    data[name], mask[name] = values, 0
            if not any(value == 0 for value in mask.values()):
                raise ValueError(f"Sample {sample_id} has neither modality")
            visits.append(data)
            masks.append(mask)
        patient_ids.append(str(patient_id))
        visits_data.append(visits)
        missing_masks.append(masks)
    features = {"protein": protein_features, "metab": metab_features}
    return patient_ids, visits_data, missing_masks, features


def to_pulse_dataset(
    subjects_path: str | Path,
    proteomics_path: str | Path,
    metabolomics_path: str | Path,
    train_mask_rate: float | None = None,
) -> tuple[LongitudinalDataset, dict[str, list[str]]]:
    patient_ids, visits, masks, features = load_sequences(
        subjects_path, proteomics_path, metabolomics_path
    )
    dataset = LongitudinalDataset(
        patient_ids=patient_ids,
        visits_data=visits,
        missing_masks=masks,
        train_mask_rate=train_mask_rate,
    )
    return dataset, features
