# PULSE independent omics validation

## Strict project isolation

Only `external_validation/` differs from the original `main` branch.
Our source adapters, download tools, additional Python requirements,
tests, instructions and experimental code are all under this directory.
The original `pulse/` implementation, `scripts/`, `tests/`, `docs/`,
root `README.md`, `requirements.txt` and `setup.py` are unchanged.

Run our utilities as Python modules from the repository root:
`python -m external_validation.scripts.download_real_omics`,
`python -m external_validation.scripts.su2020_data`,
`python -m external_validation.run`, etc. The older
`python scripts/download_real_omics.py` entry point has been moved.

Activate `conda activate pulse-omics` first. Install independent
requirements with `python -m pip install -r external_validation/requirements-structure.txt`;
HMP2 additionally requires `ijson>=3.2,<4`. Plotting uses Seaborn and
Matplotlib. Original data, prepared matrices and trained checkpoints remain
on F: under `/mnt/f/datasets/pulse-omics/`.


This directory is a **separate external-validation experiment**. It imports the unmodified `pulse.model.PulseModel` and the source-specific preparation code for Su2020 and HMP2. It does not edit the paper's implementation, the existing `scripts/` runners, the downloader, or original raw data.

## What is implemented

On each validated cohort, in **both directions** (`protein -> metab` and `metab -> protein`):

- `current_completion/pulse_current`: reconstruct the hidden current target from the current observed source.
- `historical_completion/pulse_full`: reconstruct hidden current target from current source **plus previous measured visits**.
- `historical_completion/no_history_input`: train a separate model on the **same eligible current targets**, but with no previous visits as input; this is a history-input ablation, not just an inference-time switch.
- `historical_completion/no_alignment`: retrain the historical model with the Jeffrey-divergence alignment coefficient set to zero.
- `historical_completion/last_observation`: no-learning baseline using the last available target measurement.

All learned arms use the existing PULSE model architecture without modifying its source. **The experimental training objective is adapted**: observed-cell target MSE + `0.25 * context reconstruction` + `0.001 * KL` + (`0.01 * Jeffrey alignment` for non-ablated models). Adversarial regularization is disabled to keep the first external validation tractable. This is *not* a bitwise or protocol-identical reproduction of the paper's training procedure.

Reported metrics: **feature-wise Pearson r, MAE and MSE**, unweighted macro summary and cell-weighted micro errors; **longitudinal change** metrics computed from `target(t)-target(previous actual observation)` and `predicted(t)-target(previous actual observation)` for those with a real previous measurement. For features with fewer than 3 valid paired values or zero variance, Pearson is undefined, not zero. `change_pearson` for the persistence baseline is undefined because it predicts no change.

## Prepare and inspect

Activate the existing `pulse-omics` conda environment from the repository root:

```bash
cd ~/projects/PULSE-Omics
conda activate pulse-omics
git pull --ff-only
python -m py_compile external_validation/run.py
python -m unittest discover -s external_validation/tests -q
```

The original Su2020 Table S1 must exist at either `raw/su2020/table_s1.xlsx` or `raw/su2020/pmc/mmc1.xlsx` under the data root. Build canonical inputs **only under the isolated namespace**:

`--inspect-only` examines existing canonical files and **does not run the Su2020 preparer**. For first-time preparation without training, run the existing source adapter explicitly:

```bash
python -c 'from pathlib import Path; from external_validation.omics.su2020 import prepare; p=Path("/mnt/f/datasets/pulse-omics"); print(prepare(p/"raw/su2020/table_s1.xlsx", p/"processed/external_validation/su2020"))'
python -m external_validation.run --datasets su2020 --inspect-only
```

A real run (GPU if available):

```bash
python -m external_validation.run --datasets su2020 --prepare-su2020 --device cuda --epochs 30
```

For a CPU smoke experiment reduce `--epochs 1`; do **not** treat one-epoch numbers as scientific performance.

Once the other datasets have **verified** sample/visit/assay joins, use the *same* three files under `processed/external_validation/<dataset>/`:

- `subjects.csv`: unique `sample_id`, `patient_id`, numeric `visit_order`.
- `proteomics.csv`: unique `sample_id`, then numerical protein columns; missing cells are blank/NaN.
- `metabolomics.csv`: unique `sample_id`, then numerical metabolite columns; missing cells are blank/NaN.

Missing complete assay rows are allowed. **Do not** join on row number, or use group-level labels in place of sample identifiers. In particular, currently downloaded iPOP, HMP2, COVAAR and CSIRO source matrices are **not assumed** to have verified patient/visit mapping. They are reported as skipped until canonical joins have been demonstrated. P100 individual-level assay values remain access-restricted.

```bash
python -m external_validation.run --datasets su2020 ipop hmp2 covaar csiro --inspect-only
python -m external_validation.run --datasets su2020 ipop hmp2 covaar csiro --device cuda --epochs 30
```

## Data isolation and outputs

No data or numerical results are committed to Git. They remain under:

```
/mnt/f/datasets/pulse-omics/
  raw/                                      # untouched publisher originals
  processed/external_validation/<dataset>/   # verified experimental inputs
  results/pulse_validation/<UTC run ID>/
    config.json
    status.json
    summary.csv
    ablation_comparison.csv
    report.md
    <dataset>/<source>_to_<target>/<task>/
      split.json
      features.json
      <arm>/
        metrics.json
        feature_metrics.csv
        change_feature_metrics.csv
        predictions.csv.gz                # learned arms: real/predicted per sample/feature
        train_history.csv                 # learned arms
        best.pt                           # learned arms
```

Splits are by **patient**, not visit; feature selection, means and standard deviations are computed using **training-patient samples only**. Current target values are masked from model inputs, but retained for observed-cell reconstruction supervision and held-out evaluation. Historical inputs may include real observed measurements from earlier visits only. There is no direct next-visit forecasting task in this first implementation.

The previous paper's aggregate values cannot be directly ranked against a new cohort's raw MAE/MSE because assays, units, preprocessing, and populations differ. The aim here is **same metric definitions plus independently estimated improvement over history and alignment ablations**, not a numerical reproduction of UK Biobank. The original paper's exact aggregate reference tables should be linked in the eventual comparative report with figure/table IDs and their original preprocessing definitions; this runner does not invent those values.

## Readiness / error handling

If no valid longitudinal patient pairs exist, `historical_completion` is skipped while cross-sectional completion may still run. Fewer than eight eligible patients causes a task to be marked `skipped_insufficient_paired_patients`. An inaccessible or unverified cohort is logged `SKIPPED`; this does not mean its scientific validity has been established. A run with no completed arms returns status 2.

The unit suite uses **small synthetic test fixtures only** to check plumbing, leakage prevention, output structure, and statistics. No synthetic benchmark output is included in the scientific report. After pulling on WSL, run the local tests before the expensive GPU experiments.

## Plot existing results and compare with the published paper

Plotting is **separate from model training**. It reads already completed `summary.csv`
and `feature_metrics.csv` files; it does not rerun or modify the experiment.

```bash
conda activate pulse-omics
cd ~/projects/PULSE-Omics
conda activate pulse-omics
git pull --ff-only
python -m pip install "seaborn>=0.13,<0.14" "matplotlib>=3.7"

python -m external_validation.dataset_status

python -m external_validation.visualize \
  --run-dir /mnt/f/datasets/pulse-omics/results/pulse_validation/20261109T150325Z
```

To plot the latest completed result automatically, omit `--run-dir`:
`python -m external_validation.visualize`.

Figures and comparison tables go into the **selected experiment's** `figures/` directory:

- `01_model_overview.{png,svg}`: per-direction model performance (median per-feature Pearson r and micro-averaged MAE).
- `02_feature_correlations.{png,svg}`: Pearson-r distributions across measured assay features.
- `03_paired_feature_effects.{png,svg}`: feature-matched full PULSE versus no-history and no-alignment ablations.
- `04_paper_context.{png,svg}`: explicit two-panel comparison with the **six exact Pearson-r values printed in the original paper's Fig. 2c**.
- `paper_comparison.csv`, `paired_feature_effects.csv`, `comparison_report.md`: machine-readable and narrative comparison.

Reference values are in `external_validation/paper_reference.csv`, transcribed
from Wu et al., *Nature Computational Science* 6, 954–970 (2026),
Fig. 2c, p. 957, DOI `10.1038/s43588-026-01026-5`.
That figure evaluates **longitudinal change in one large-LDL lipid-ratio marker**.
Our independent evaluation uses median per-feature change correlation across
an assay with different inputs and biomarkers. The two panels are shown for
context, **not** as an identical-target numerical reproduction. The precise
Fig. 2a summary medians have **not** been extracted from Supplementary Table 2,
so they are not populated or invented.

**Important for the previously completed Su2020 run:** the current-only
and historical task arms used separately selected sample populations/splits.
Do not subtract their medians and call this the historical-information effect.
Within the historical task, the full/no-history/no-alignment arms share the
same held-out cases. Those are the intended within-cohort contrasts.

## Raw downloads versus experimental readiness

`python -m external_validation.dataset_status` summarizes, for each dataset:
local raw file count, incomplete `.part` downloads, last downloader
`download_manifest.json` state, issues from the latest pass, and whether
the three required canonical matrices exist. It makes **no network requests**.
The manifest status `enumerated_in_this_run` is a description of the last
source enumeration attempt, **not proof of complete public-source coverage**.
A dataset with downloaded raw files may still lack a verified patient/visit
mapping and therefore remain ineligible for experiments.

To resume HMP2 only, using the existing resumable downloader:

```bash
conda activate pulse-omics
cd ~/projects/PULSE-Omics
python -m external_validation.scripts.download_real_omics --datasets hmp2 --workers 4
```

The initial HMP2 clinical metadata, proteomics and metabolomics matrices
can be investigated before downloading every unrelated processed product.
The external validation script will not run HMP2 until its identifiers
and visit relationships are genuinely verified.


## HMP2 official matrices: read, verify, prepare

The downloader has already obtained the HMP2 processed product inventory. The
independent preparer reads **only three primary published source files**:

- `raw/hmp2/hmp2_metadata_2018-08-20.csv`: clinical patient and visit mapping
- `raw/hmp2/HMP2_proteomics_1pep-pro_1p_FDR_5ppm.tsv.gz`: protein matrix
- `raw/hmp2/HMP2_metabolomics_w_metadata.biom.gz`: metabolomics matrix

It requires an **identical sample ID** in both matrices, an unambiguous
`External ID -> (Participant ID, week_num)` in the metadata and a unique
patient/week assignment. Conflicting or repeated patient-week mappings are
excluded; row-position joins and fabricated visit identities are forbidden.

```bash
cd ~/projects/PULSE-Omics
conda activate pulse-omics
git pull --ff-only
python -m pip install "ijson>=3.2,<4"
python -m unittest discover -s external_validation/tests -p "test_prepare_hmp2.py" -q

# First: number of truly pairable specimens, with no output files written.
python -m external_validation.prepare_hmp2 --inspect-only

# If inspection yields sufficient paired longitudinal specimens:
python -m external_validation.prepare_hmp2
python -m external_validation.run --datasets hmp2 --inspect-only
```

Exports are written exclusively to `processed/external_validation/hmp2/`:
three canonical CSVs and `preparation_report.json`. No raw files, original
PULSE implementation, or earlier Su2020 experiment outputs are changed.

The BIOM file is parsed by `ijson` and a temporary disk-backed matrix to
avoid loading the entire metabolomics JSON into memory. **BIOM sparse implicit
entries are kept as zeros, with this convention recorded in the report.**
These values should be examined in the downstream assay missingness analysis
before interpreting quantitative accuracy or training on all features.

If no exactly identified dual-omics samples exist, the preparer reports
`insufficient_exact_matched_visits` and writes **no** canonical assay CSVs.
Do not force pairing by patient-week alone: that identifies a visit but does
not establish that two sample IDs represent the same biological specimen.
