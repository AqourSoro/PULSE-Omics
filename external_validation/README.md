# PULSE independent omics validation

This directory is a **separate external-validation experiment**. It imports the unmodified `pulse.model.PulseModel` and (only for Su2020 preparation) `omics.su2020.prepare`. It does not edit the paper's implementation, the existing `scripts/` runners, the downloader, or original raw data.

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
git pull --ff-only
python -m py_compile external_validation/run.py
python -m unittest discover -s external_validation/tests -q
```

The original Su2020 file must exist at `/mnt/f/datasets/pulse-omics/raw/su2020/table_s1.xlsx`. Build canonical inputs **only under the isolated namespace**:

`--inspect-only` examines existing canonical files and **does not run the Su2020 preparer**. For first-time preparation without training, run the existing source adapter explicitly:

```bash
python -c 'from pathlib import Path; from omics.su2020 import prepare; p=Path("/mnt/f/datasets/pulse-omics"); print(prepare(p/"raw/su2020/table_s1.xlsx", p/"processed/external_validation/su2020"))'
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