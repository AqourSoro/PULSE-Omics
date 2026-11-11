"""Su2020 fast pilot: paired low-rank longitudinal state + conditional Gaussian copula.

This is a deliberately *lightweight*, real-calibrated FEASIBILITY BASELINE inspired
by MEFISTO's shared factors and scDesign3's conditional marginals/copulas. It does
NOT fit either published software package; see the JSON output for model provenance.
T1/T2 are stage indicators, not measured elapsed times.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, wasserstein_distance
from sklearn.covariance import LedoitWolf
from sklearn.decomposition import FactorAnalysis
from sklearn.linear_model import Ridge

from external_validation.run import load_canonical

MODS = ("protein", "metab")


def split_patients(subjects: pd.DataFrame, seed: int = 42) -> dict[str, set[str]]:
    patients = np.array(sorted(subjects.patient_id.unique()), dtype=str)
    if len(patients) < 12:
        raise ValueError("Need >=12 patients for an independent test")
    rng = np.random.default_rng(seed)
    rng.shuffle(patients)
    nt = max(3, round(.15 * len(patients)))
    nv = max(3, round(.15 * len(patients)))
    return {
        "test": set(patients[:nt]),
        "val": set(patients[nt:nt + nv]),
        "train": set(patients[nt + nv:]),
    }


def select_features(raw: dict[str, pd.DataFrame], subjects: pd.DataFrame,
                    ids: set[str], nfeatures: int) -> dict[str, list[str]]:
    train_samples = subjects.loc[subjects.patient_id.isin(ids), "sample_id"]
    result = {}
    for mod in MODS:
        values = raw[mod].loc[train_samples]
        finite = np.isfinite(values.to_numpy(dtype=float))
        coverage = finite.mean(axis=0)
        spread = np.nanstd(values.to_numpy(dtype=float), axis=0)
        eligible = np.where((coverage >= .65) & (spread > 1e-8))[0]
        if eligible.size < nfeatures:
            raise ValueError(f"{mod}: only {eligible.size} features have >=65% observed values and nonzero variance")
        # Favor reliable measurements; use spread only as a tie breaker.
        chosen = sorted(eligible, key=lambda k: (-coverage[k], -spread[k], k))[:nfeatures]
        result[mod] = values.columns[chosen].tolist()
    return result


def assemble(subjects: pd.DataFrame, raw: dict[str, pd.DataFrame],
             features: dict[str, list[str]], ids: set[str]) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    metadata = subjects.loc[subjects.patient_id.isin(ids)].copy()
    metadata = metadata.sort_values(["patient_id", "visit_order", "sample_id"]).reset_index(drop=True)
    # Both modalities must be present at this visit (partial feature missingness is allowed).
    p = raw["protein"].reindex(metadata.sample_id)[features["protein"]].to_numpy(dtype=float)
    m = raw["metab"].reindex(metadata.sample_id)[features["metab"]].to_numpy(dtype=float)
    valid = np.isfinite(p).any(axis=1) & np.isfinite(m).any(axis=1)
    return metadata.loc[valid].reset_index(drop=True), p[valid], m[valid]


def estimate_scaling(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.nanmean(x, axis=0)
    sd = np.nanstd(x, axis=0)
    if not np.all(np.isfinite(mean)) or np.any(sd < 1e-8):
        raise ValueError("Insufficient nonmissing data for selected features")
    return mean, sd


def scale(x: np.ndarray, mean: np.ndarray, sd: np.ndarray) -> np.ndarray:
    out = (x - mean) / sd
    return np.where(np.isfinite(out), out, 0.0)


def covariance_draw(rng: np.random.Generator, values: np.ndarray, n: int,
                    *, minimum: float = 1e-4) -> np.ndarray:
    """New independent latent draws; shrink covariance for small longitudinal cohorts."""
    if values.shape[0] < 3:
        raise ValueError("Not enough patients to fit a latent distribution")
    if values.shape[1] == 1:
        variance = max(float(np.var(values, ddof=1)), minimum)
        return rng.normal(values.mean(axis=0), np.sqrt(variance), size=(n, 1))
    cov = LedoitWolf().fit(values).covariance_
    eig, vec = np.linalg.eigh(cov)
    transform = vec @ np.diag(np.sqrt(np.maximum(eig, minimum)))
    return rng.standard_normal((n, values.shape[1])) @ transform.T + values.mean(axis=0)


class FactorCopula:
    """Shared-factor conditional mean + Gaussian-copula empirical residuals."""

    def __init__(self, n_factors: int = 6, seed: int = 42):
        self.n_factors = n_factors
        self.seed = seed

    def fit(self, x: np.ndarray, metadata: pd.DataFrame) -> "FactorCopula":
        if not np.isfinite(x).all():
            raise ValueError("Inputs must be train-normalized with missing cells imputed")
        self.factor = FactorAnalysis(n_components=min(self.n_factors, x.shape[1] - 1),
                                     random_state=self.seed, max_iter=250)
        z = self.factor.fit_transform(x)
        self.dim = x.shape[1]
        self.training_scores = z
        self.predict_mean = self.factor.mean_
        self.loading = self.factor.components_.T
        # The factor-analysis posterior has shrinkage; fit the conditional mean
        # to the *inferred* scores, and model remaining dependence separately.
        conditional = np.c_[np.ones(len(z)), z]
        self.coef = np.linalg.lstsq(conditional, x, rcond=1e-6)[0]
        residual = x - conditional @ self.coef
        n = len(residual)
        q = np.column_stack([norm.ppf(np.clip(rankdata(residual[:, j]) / (n + 1), .001, .999))
                             for j in range(self.dim)])
        corr = LedoitWolf().fit(q).covariance_
        std = np.sqrt(np.diag(corr))
        corr = corr / np.outer(std, std)
        eig, vec = np.linalg.eigh(corr)
        corr = (vec * np.maximum(eig, 1e-6)) @ vec.T
        std = np.sqrt(np.diag(corr))
        self.corr = corr / np.outer(std, std)
        self.residual_quantiles = np.sort(residual, axis=0)
        self.metadata = metadata.reset_index(drop=True).copy()
        self.pairs = {}
        for pid, indices in self.metadata.groupby("patient_id").indices.items():
            stage = self.metadata.iloc[indices].visit_order.to_numpy(dtype=float)
            if len(indices) >= 2 and np.unique(stage).size >= 2:
                ix = np.argsort(stage)
                self.pairs[str(pid)] = (z[indices[ix[0]]], z[indices[ix[-1]]])
        if len(self.pairs) < 8:
            raise ValueError(f"Only {len(self.pairs)} usable longitudinal patient pairs; need at least 8")
        self.baseline = np.vstack([p[0] for p in self.pairs.values()])
        self.transitions = np.vstack([p[1] - p[0] for p in self.pairs.values()])
        return self

    def sample(self, n_patients: int, rng: np.random.Generator) -> tuple[pd.DataFrame, np.ndarray]:
        base = covariance_draw(rng, self.baseline, n_patients)
        delta = covariance_draw(rng, self.transitions, n_patients)
        z = np.stack([base, base + delta], axis=1).reshape(n_patients * 2, -1)
        mu = np.c_[np.ones(len(z)), z] @ self.coef
        # Gaussian copula draws with empirical conditional residual marginals.
        eig, vec = np.linalg.eigh(self.corr)
        noise = rng.standard_normal((len(z), self.dim)) @ (vec * np.sqrt(np.maximum(eig, 1e-8))).T
        uniform = norm.cdf(noise)
        empirical_positions = np.linspace(0.0, 1.0, self.residual_quantiles.shape[0])
        residual = np.column_stack([np.interp(uniform[:, j], empirical_positions,
                                               self.residual_quantiles[:, j]) for j in range(self.dim)])
        synthetic = mu + residual
        meta = pd.DataFrame({
            "patient_id": np.repeat([f"SYN_{i:05d}" for i in range(n_patients)], 2),
            "visit_order": np.tile([1, 2], n_patients),
            "sample_id": [f"SYN_{i:05d}_T{t}" for i in range(n_patients) for t in (1, 2)],
        })
        return meta, synthetic


def distributions(real: np.ndarray, fake: np.ndarray) -> float:
    return float(np.median([wasserstein_distance(real[:, j], fake[:, j])
                            for j in range(real.shape[1])]))


def correlation_gap(real: np.ndarray, fake: np.ndarray, n_protein: int) -> dict[str, float]:
    a = np.corrcoef(real, rowvar=False)
    b = np.corrcoef(fake, rowvar=False)
    p = n_protein
    def error(r0, r1):
        return float(np.mean(np.abs(r0-r1)))
    return {
        "protein_pearson_mae": error(a[:p, :p], b[:p, :p]),
        "metab_pearson_mae": error(a[p:, p:], b[p:, p:]),
        "cross_pearson_mae": error(a[:p, p:], b[:p, p:]),
    }


def evaluate_ridge(train: np.ndarray, test: np.ndarray, np_: int) -> dict[str, float]:
    source, target = train[:, :np_], train[:, np_:]
    predictor = Ridge(alpha=100.).fit(source, target)
    truth = test[:, np_:]
    pred = predictor.predict(test[:, :np_])
    return {"zmae": float(np.mean(np.abs(pred-truth))), "zmse": float(np.mean((pred-truth)**2))}


def evaluate_pulse(train: np.ndarray, validation: np.ndarray, test: np.ndarray,
                   np_: int, seed: int, epochs: int) -> dict[str, float]:
    import torch
    from pulse.model import PulseModel
    torch.manual_seed(seed)
    torch.set_num_threads(1)
    model = PulseModel(modality_dims={"protein": np_, "metab": train.shape[1]-np_},
                       latent_dim=8, hidden_dim=32, lambda_adv=0,
                       use_adversarial_loss=False)
    optim = torch.optim.AdamW(model.parameters(), lr=.001)
    best, best_state, stale = float("inf"), None, 0
    def predict(x):
        model.eval()
        with torch.no_grad():
            xx = torch.from_numpy(x[:, :np_].astype(np.float32))
            out = model([{"protein": xx, "metab": None}],
                        [{"protein": 0, "metab": 2}])
        return out["reconstructions"][0]["metab"].numpy()
    rng = np.random.default_rng(seed)
    for epoch in range(epochs):
        model.train()
        for subset in np.array_split(rng.permutation(len(train)), max(1, int(np.ceil(len(train) / 32)))):
            batch = torch.from_numpy(train[subset].astype(np.float32))
            source, target = batch[:, :np_], batch[:, np_:]
            result = model([{"protein": source, "metab": None}],
                           [{"protein": 0, "metab": 1}],
                           recon_targets=[{"protein": source, "metab": target}])
            decoded = result["reconstructions"][0]
            penalty = result["losses"]
            loss = (torch.mean((decoded["metab"] - target)**2)
                    + .25 * torch.mean((decoded["protein"] - source)**2)
                    + .001 * penalty["kl_loss"] + .01 * penalty["alignment_loss"])
            optim.zero_grad()
            loss.backward()
            optim.step()
        zval = np.mean((predict(validation) - validation[:, np_:])**2)
        if zval < best-1e-6:
            best, stale = zval, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            stale += 1
        if stale >= 8:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    estimate = predict(test)
    return {"zmae": float(np.mean(np.abs(estimate - test[:, np_:]))),
            "zmse": float(np.mean((estimate - test[:, np_:])**2)),
            "epochs_max": epochs}


def draw_figures(out: Path, true: np.ndarray, synth: np.ndarray, np_: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    axes[0].scatter(np.quantile(true, [.1, .25, .5, .75, .9], axis=0).flatten(),
                    np.quantile(synth, [.1, .25, .5, .75, .9], axis=0).flatten(),
                    alpha=.35, s=8)
    axes[0].set_xlabel("Held-out real quantiles (training z-score)")
    axes[0].set_ylabel("Synthetic quantiles")
    for ax, x, name in zip(axes[1:], (true, synth), ("Real test", "Synthetic")):
        corr = np.corrcoef(x, rowvar=False)
        ax.imshow(corr[:np_, np_:], vmin=-1, vmax=1, cmap="coolwarm", aspect="auto")
        ax.set_title(f"{name}: protein–metabolite correlations")
        ax.set_xlabel("Metabolites")
        ax.set_ylabel("Proteins")
    fig.tight_layout()
    fig.savefig(out / "diagnostics.png", dpi=160)
    plt.close(fig)


def run(args: argparse.Namespace) -> dict:
    folder = args.data_root / "processed" / "external_validation" / "su2020"
    if args.prepare:
        from external_validation.real_study import prepare
        prepare(args.data_root, "su2020", download=True)
    subjects, raw = load_canonical(folder)
    partitions = split_patients(subjects, args.seed)
    features = select_features(raw, subjects, partitions["train"], args.features)
    rows, assay = {}, {}
    for part in ("train", "val", "test"):
        metadata, p, m = assemble(subjects, raw, features, partitions[part])
        rows[part] = metadata
        assay[part] = np.hstack((p, m))
        if len(metadata) < 8:
            raise ValueError(f"Too few paired measurements in {part}: {len(metadata)}")
    center, sd = estimate_scaling(assay["train"])
    x = {key: scale(value, center, sd) for key, value in assay.items()}
    gen = FactorCopula(args.factors, args.seed).fit(x["train"], rows["train"])
    rng = np.random.default_rng(args.seed + 999)
    synth_meta, synth = gen.sample(args.patients, rng)
    # Matched training patient count for fair synthetic-to-real prediction.
    eligible = max(12, len(partitions["train"]))
    matched_synth = synth[:min(len(synth), 2 * eligible)]
    # Validation in the synthetic domain is independent of training data.
    _, extra = gen.sample(max(16, len(partitions["val"])), rng)
    # Negative control: preserve all univariate marginals but remove dependence.
    shuffled = np.column_stack([rng.permutation(synth[:, j]) for j in range(synth.shape[1])])
    np_ = args.features
    metrics = {
        "model": "shared FactorAnalysis + T1/T2 latent transitions + conditional empirical residual Gaussian copula",
        "model_provenance": "MEFISTO/scDesign3-inspired pilot, NOT the original packages",
        "time_resolution": "T1/T2 categorical visits; no elapsed-day or GP calibration",
        "seed": args.seed,
        "n_features_each": args.features,
        "n_factors": gen.factor.n_components,
        "patients": {key: len(ids) for key, ids in partitions.items()},
        "paired_rows": {key: len(rows[key]) for key in rows},
        "paired_training_patients": len(gen.pairs),
        "synthetic_patients": args.patients,
        "synthetic_rows": len(synth),
        "median_univariate_wasserstein": {
            "synthetic_vs_test": distributions(x["test"], synth),
            "training_vs_test": distributions(x["test"], x["train"]),
            "shuffled_control_vs_test": distributions(x["test"], shuffled),
        },
        "correlation_error": {
            "synthetic": correlation_gap(x["test"], synth, np_),
            "training": correlation_gap(x["test"], x["train"], np_),
            "shuffled_control": correlation_gap(x["test"], shuffled, np_),
        },
        "conditional_prediction_on_real_test": {
            "ridge_train_real": evaluate_ridge(x["train"], x["test"], np_),
            "ridge_train_synthetic": evaluate_ridge(matched_synth, x["test"], np_),
            "ridge_train_shuffled_control": evaluate_ridge(shuffled[:len(matched_synth)], x["test"], np_),
        },
    }
    if not args.skip_pulse:
        metrics["conditional_prediction_on_real_test"]["pulse_train_real"] = evaluate_pulse(
            x["train"], x["val"], x["test"], np_, args.seed, args.epochs)
        metrics["conditional_prediction_on_real_test"]["pulse_train_synthetic"] = evaluate_pulse(
            matched_synth, extra, x["test"], np_, args.seed, args.epochs)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "metrics.json").write_text(json.dumps(metrics, ensure_ascii=False, indent=2)+"\n")
    (args.output / "selected_features.json").write_text(json.dumps(features, ensure_ascii=False, indent=2)+"\n")
    (args.output / "patient_splits.json").write_text(json.dumps({k: sorted(v) for k, v in partitions.items()}, indent=2)+"\n")
    output_frame = pd.DataFrame(synth * sd + center, columns=features["protein"] + features["metab"])
    pd.concat([synth_meta, output_frame], axis=1).to_csv(args.output / "generated_complete.csv", index=False)
    draw_figures(args.output, x["test"], synth, np_)
    return metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, default=Path("/mnt/f/datasets/pulse-omics"))
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--prepare", action="store_true", help="Fetch and prepare Su2020 if canonical tables are missing")
    ap.add_argument("--features", type=int, default=50)
    ap.add_argument("--factors", type=int, default=6)
    ap.add_argument("--patients", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--epochs", type=int, default=35)
    ap.add_argument("--skip-pulse", action="store_true")
    args = ap.parse_args()
    if min(args.features, args.factors, args.patients, args.epochs) < 1:
        ap.error("All size parameters must be positive")
    report = run(args)
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()