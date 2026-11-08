# Su et al. (2020) plasma proteomics/metabolomics data

Reference: Su Y. et al., *Multi-Omics Resolves a Sharp Disease-State Shift between Mild and Moderate COVID-19*, Cell (2020), DOI: 10.1016/j.cell.2020.10.037.

Source: [Mendeley Data v5](https://data.mendeley.com/datasets/tzydswhhb5/5), DOI: 10.17632/tzydswhhb5.5. Dataset license: CC BY 4.0. Table S1 comprises subject details and plasma proteomic/metabolomic data. The same supplementary workbook is linked from [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC7598382/).

## Storage contract

Keep code at `~/projects/PULSE-Omics` and all datasets on the Windows F: drive:

```
/mnt/f/datasets/pulse-omics/
  raw/su2020/table_s1.xlsx
  raw/su2020/download_manifest.json
  processed/su2020/workbook_structure.json
  synthetic/
  results/
  download_logs/
```

The downloader honors `PULSE_OMICS_DATA_ROOT` or `--data-root`. It never silently writes downloaded data under the Git repository.

## Download and inspect (WSL)

```bash
cd ~/projects/PULSE-Omics
conda activate pulse-omics
git pull --ff-only
python -m pip install -r requirements-omics.txt
python scripts/su2020_data.py fetch
python scripts/su2020_data.py inspect
```

The downloader tries [Europe PMC's supplementary ZIP endpoint](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7598382/supplementaryFiles) first and then PMC's [post-August-2026 AWS article object storage](https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/). These sources may not distribute the supplement depending on licensing or connectivity; they are attempted, not guaranteed. It selects only the Table S1 workbook from an archive, validates that the file is a real XLSX, and records SHA-256. In November 2026 the old PMC attachment path returned a browser interstitial and the anonymous Mendeley API returned HTTP 401; neither is used as an automated fallback.

If automated retrieval fails, visit the [original Mendeley Data v5 page](https://data.mendeley.com/datasets/tzydswhhb5/5), download the version 5 Table S1 workbook using your browser, then import the file from its WSL-visible Windows Downloads path:

```bash
python scripts/su2020_data.py import --file '/mnt/c/Users/YOUR_WINDOWS_USER/Downloads/Table S1.xlsx'
python scripts/su2020_data.py inspect
```

You can instead use a ZIP downloaded through the site's `Download All` button; the importer extracts only Table S1. Keep the downloaded file's actual name in the `--file` argument. No account token is embedded or logged. The importer checks the XLSX structure before placing it on F:. Downloads and imports never overwrite existing data unless `--overwrite` is passed.

Inspection saves actual worksheet names, dimensions, and a small top-left cell preview in `workbook_structure.json`.

**Important:** the Su study reports 139 COVID-19 patients and serial samples, but that does not establish that every patient has two usable *paired* proteomic/metabolomic timepoints. Actual sample IDs, assay availability, and visit times must be checked from Table S1 before constructing a patient-level longitudinal split. Do not infer linkage from row order or construct a synthetic time variable from spreadsheet positions.

After inspection, implement a schema-specific adapter matching the observed worksheet identifiers, retain assay missingness, and split data by patient before fitting any feature scaling. All raw clinical/omics measurements remain on F:, outside version control.

Dataset download and analysis on the user's WSL workstation are separate from modifications committed to GitHub.

## PULSE longitudinal data adapter

`omics/longitudinal.py` now accepts three **validated canonical CSV files**, once the workbook's actual identifiers have been resolved:

- `subjects.csv`: columns `sample_id,patient_id,visit_order` (one row per blood draw).
- `proteomics.csv`: one row per `sample_id` with protein abundance columns.
- `metabolomics.csv`: one row per `sample_id` with metabolite abundance columns.

```python
from omics.longitudinal import to_pulse_dataset
dataset, features = to_pulse_dataset(
    subjects_path="/mnt/f/datasets/pulse-omics/processed/su2020/subjects.csv",
    proteomics_path="/mnt/f/datasets/pulse-omics/processed/su2020/proteomics.csv",
    metabolomics_path="/mnt/f/datasets/pulse-omics/processed/su2020/metabolomics.csv",
    train_mask_rate=0.6,
)
print(len(dataset), {k: len(v) for k, v in features.items()})
```

The adapter rejects duplicate sample IDs, contradictory visit assignments, unmatched assay IDs, and partially missing/nonfinite feature vectors; an entire missing assay uses PULSE mask 2. It does not preprocess or normalize measurements. The three canonical CSVs are **not created by the workbook inspector**; that extraction step must be implemented after examining the real workbook layout.

## Identify the actual patient and assay join keys

After the Table S1 download and workbook inspection succeed, run:

```bash
python scripts/su2020_schema.py
```

This uses the three actual sheets `S1.1 Patient Clinical Data`, `S1.3 Plasma Proteomic Data`, and `S1.4 Plasma Metabolomic Data`. It reports unique-ID overlaps between tables and potential repeating patient IDs, then saves a **local-only** diagnostic:

`/mnt/f/datasets/pulse-omics/processed/su2020/schema_profile.json`

The diagnostics distinguish assay-row identities from patient identities; matching candidates are **not** automatically promoted into confirmed longitudinal links. The report contains top-left sheet previews to support explicit schema mapping, so review before sharing it outside the local machine. No training/validation split or omics matrix has been generated at this stage.
