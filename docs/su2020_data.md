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

The downloader first tries the paper's PMC supplementary Table S1. If it is blocked, it resolves the v5 file using Mendeley Data's **current** public API (`api.data.mendeley.com`). It validates the downloaded XLSX container and writes a SHA-256 manifest. If both sources fail, it reports their errors and does not create a fake workbook. It does not overwrite an existing file unless `--overwrite` is specified.

Inspection saves actual worksheet names, dimensions, and a small top-left cell preview in `workbook_structure.json`.

**Important:** the Su study reports 139 COVID-19 patients and serial samples, but that does not establish that every patient has two usable *paired* proteomic/metabolomic timepoints. Actual sample IDs, assay availability, and visit times must be checked from Table S1 before constructing a patient-level longitudinal split. Do not infer linkage from row order or construct a synthetic time variable from spreadsheet positions.

After inspection, implement a schema-specific adapter matching the observed worksheet identifiers, retain assay missingness, and split data by patient before fitting any feature scaling. All raw clinical/omics measurements remain on F:, outside version control.

Dataset download and analysis on the user's WSL workstation are separate from modifications committed to GitHub.
