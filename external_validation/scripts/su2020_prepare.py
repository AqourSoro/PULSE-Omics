#!/usr/bin/env python3
"""Prepare Su2020 raw assay matrices and verified patient-level sample links."""
from __future__ import annotations

import argparse
import os
from pathlib import Path

from external_validation.omics.su2020 import prepare

DATA_ROOT = "/mnt/f/datasets/pulse-omics"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root", default=os.environ.get("PULSE_OMICS_DATA_ROOT", DATA_ROOT)
    )
    args = parser.parse_args()
    root = Path(args.data_root).expanduser()
    workbook = root / "raw" / "su2020" / "table_s1.xlsx"
    output = root / "processed" / "su2020"
    report = prepare(workbook, output)
    print("Su2020 preparation completed")
    print(f"  Clinical patients: {report['clinical_patients']}")
    print(f"  Clinical draws: {report['clinical_rows']}")
    print(f"  Patients with T1/T2: {report['clinical_patients_with_two_draws']}")
    print(f"  Assay-matched samples (no time inference): {report['assay_shared_sample_ids']}")
    for name, kinds in report["assay_id_formats"].items():
        print(f"  {name} sample_id formats: {kinds}")
    for name, count in report["clinical_sample_matches"].items():
        breakdown = report["match_sources"][name]
        print(f"  Clinical samples matched to {name}: {count} "
              f"(assay sample_id={breakdown['assay_sample_id']}, "
              f"label={breakdown['assay_label']})")
    for name, info in report["clinical_bridge_candidates"].items():
        print(f"  Clinical bridge candidates for {name}: {info[:4]}")
    if "linked_clinical_draws" in report:
        print(f"  Linked clinical draws: {report['linked_clinical_draws']}")
        print(f"  With both assays: {report['linked_paired_draws']}")
        print(f"  Linked patients with T1/T2: {report['linked_patients_with_two_draws']}")
        print(f"  Patients with BOTH assays at BOTH visits: "
              f"{report['patients_with_two_paired_omics_visits']}")
    print(f"  Result status: {report['status']}")
    print(f"  Report: {output / 'preparation_report.json'}")


if __name__ == "__main__":
    main()
