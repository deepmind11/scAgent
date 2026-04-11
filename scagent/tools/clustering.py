"""Leiden clustering with resolution sweep and seed sensitivity check."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import warnings
import scanpy as sc

# Suppress Scanpy's Leiden backend deprecation warning
warnings.filterwarnings("ignore", message=".*default backend for leiden.*")


def run_leiden(
    adata: ad.AnnData,
    *,
    resolution: float = 1.0,
    key_added: str = "leiden",
    random_state: int = 0,
    n_iterations: int = -1,
    checkpoint_dir: str | None = None,
) -> dict:
    """Run Leiden clustering. Modifies *adata* in place.

    Adds cluster labels to ``adata.obs[key_added]``.

    Parameters
    ----------
    resolution
        Higher → more clusters. 0.3–0.5 for major lineages, 0.8–1.2
        standard, 1.5–3.0 for subtypes.
    key_added
        Key in ``adata.obs`` for cluster labels.
    random_state
        Random seed for reproducibility.
    n_iterations
        -1 means iterate until convergence.
    checkpoint_dir
        Save checkpoint after clustering.

    Returns
    -------
    Result dict with metrics, provenance.
    """
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=key_added,
        random_state=random_state,
        n_iterations=n_iterations,
    )

    labels = adata.obs[key_added]
    n_clusters = labels.nunique()
    cluster_sizes = labels.value_counts().sort_index()
    min_cluster_size = int(cluster_sizes.min())

    warnings = []
    if n_clusters < 2:
        warnings.append("Only 1 cluster found. Resolution may be too low.")
    if n_clusters > adata.n_obs * 0.1:
        warnings.append(
            f"{n_clusters} clusters for {adata.n_obs} cells — likely over-clustered. "
            "Lower the resolution."
        )
    if min_cluster_size < 10:
        small = cluster_sizes[cluster_sizes < 10]
        warnings.append(
            f"Clusters with <10 cells: {dict(small)}. These may be artifacts."
        )

    result = {
        "metrics": {
            "n_clusters": n_clusters,
            "resolution": resolution,
            "cluster_sizes": {str(k): int(v) for k, v in cluster_sizes.items()},
            "min_cluster_size": min_cluster_size,
            "max_cluster_size": int(cluster_sizes.max()),
        },
        "plots": [],
        "provenance": {
            "tool_id": "leiden_clustering",
            "parameters": {
                "resolution": resolution,
                "key_added": key_added,
                "random_state": random_state,
                "n_iterations": n_iterations,
            },
            "n_clusters": n_clusters,
        },
        "warnings": warnings,
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_leiden")

    return result


def sweep_resolution(
    adata: ad.AnnData,
    *,
    resolutions: list[float] | None = None,
    random_state: int = 0,
    n_seeds_stability: int = 5,
    plot_dir: str | None = None,
) -> dict:
    """Run Leiden at multiple resolutions with clustering quality metrics.

    For each resolution, computes:
    - **Number of clusters**
    - **Silhouette score** — how well-separated clusters are in PCA space
      (range -1 to 1; higher is better)
    - **Calinski-Harabasz index** — ratio of between- to within-cluster
      variance (higher is better; not comparable across datasets)
    - **Stability (mean ARI)** — consistency across random seeds
      (range 0 to 1; ≥0.9 is stable)

    Recommends the resolution that maximises silhouette score among
    resolutions with stability ≥ 0.9. If no resolution is stable, falls
    back to the highest silhouette score and warns.

    Each resolution gets its own key: ``adata.obs['leiden_{res}']``.

    Parameters
    ----------
    resolutions
        Resolutions to try. Defaults to
        ``[0.3, 0.5, 0.8, 1.0, 1.5, 2.0]``.
    random_state
        Base seed — used for the primary clustering at each resolution.
    n_seeds_stability
        Number of seeds for the stability check at each resolution.
        Set to 0 to skip stability analysis (faster).
    plot_dir
        Directory for the sweep summary plot.

    Returns
    -------
    Result dict with per-resolution metrics and a recommended resolution.
    """
    from sklearn.metrics import (
        silhouette_score,
        calinski_harabasz_score,
        adjusted_rand_score,
    )

    if resolutions is None:
        resolutions = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]

    # Get PCA embedding for metric computation
    if "X_pca" not in adata.obsm:
        raise ValueError("PCA not found. Run run_pca() before sweep_resolution().")
    X_pca = adata.obsm["X_pca"]

    sweep_rows = []

    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=key,
            random_state=random_state,
            n_iterations=-1,
        )
        labels = adata.obs[key]
        n_clusters = labels.nunique()

        # --- quality metrics (need ≥2 clusters) ---
        if n_clusters >= 2:
            labels_int = labels.astype(int).values
            sil = float(silhouette_score(X_pca, labels_int, sample_size=min(5000, len(labels_int)), random_state=random_state))
            ch = float(calinski_harabasz_score(X_pca, labels_int))
        else:
            sil = float("nan")
            ch = float("nan")

        # --- stability across seeds ---
        mean_ari = float("nan")
        if n_seeds_stability >= 2 and n_clusters >= 2:
            all_labels = []
            seeds = list(range(n_seeds_stability))
            for s in seeds:
                tmp_key = f"_sweep_stab_{res}_{s}"
                sc.tl.leiden(
                    adata, resolution=res, key_added=tmp_key,
                    random_state=s, n_iterations=-1,
                )
                all_labels.append(adata.obs[tmp_key].values.copy())
                del adata.obs[tmp_key]

            ari_vals = []
            for i in range(len(seeds)):
                for j in range(i + 1, len(seeds)):
                    ari_vals.append(adjusted_rand_score(all_labels[i], all_labels[j]))
            mean_ari = float(np.mean(ari_vals))

        sweep_rows.append({
            "resolution": res,
            "n_clusters": n_clusters,
            "silhouette": round(sil, 4) if not np.isnan(sil) else None,
            "calinski_harabasz": round(ch, 1) if not np.isnan(ch) else None,
            "stability_ari": round(mean_ari, 4) if not np.isnan(mean_ari) else None,
        })

    # --- recommend ---
    # Prefer: stable (ARI ≥ 0.9) AND highest silhouette
    stable_rows = [r for r in sweep_rows if r["stability_ari"] is not None and r["stability_ari"] >= 0.9 and r["silhouette"] is not None]
    if stable_rows:
        best = max(stable_rows, key=lambda r: r["silhouette"])
        recommendation_note = (
            f"Resolution {best['resolution']} recommended: highest silhouette "
            f"({best['silhouette']:.3f}) among stable resolutions (ARI ≥ 0.9)."
        )
    else:
        # Fallback: highest silhouette regardless of stability
        valid_rows = [r for r in sweep_rows if r["silhouette"] is not None]
        if valid_rows:
            best = max(valid_rows, key=lambda r: r["silhouette"])
            recommendation_note = (
                f"Resolution {best['resolution']} has the highest silhouette "
                f"({best['silhouette']:.3f}), but NO resolution achieved ARI ≥ 0.9. "
                "Clustering may be fragile — inspect carefully."
            )
        else:
            best = sweep_rows[0]
            recommendation_note = "Could not compute quality metrics."

    recommended = best["resolution"]

    plots = []
    if plot_dir is not None:
        plots = _plot_sweep_metrics(sweep_rows, recommended, plot_dir)

    warnings = []
    if all(r.get("stability_ari") is not None and r["stability_ari"] < 0.9 for r in sweep_rows):
        warnings.append(
            "No resolution achieved stable clustering (ARI ≥ 0.9). "
            "Consider adjusting n_neighbors or inspecting the neighbor graph."
        )

    return {
        "metrics": {
            "sweep": sweep_rows,
            "recommended_resolution": recommended,
            "recommendation_note": recommendation_note,
        },
        "plots": plots,
        "provenance": {
            "tool_id": "leiden_clustering",
            "parameters": {
                "resolutions": resolutions,
                "random_state": random_state,
                "n_seeds_stability": n_seeds_stability,
                "mode": "sweep",
            },
        },
        "warnings": warnings,
    }


def check_seed_sensitivity(
    adata: ad.AnnData,
    *,
    resolution: float = 1.0,
    seeds: list[int] | None = None,
) -> dict:
    """Run Leiden at the same resolution with different seeds.

    Computes Adjusted Rand Index (ARI) between each pair of runs.
    If ARI < 0.9, the clustering at this resolution is fragile.

    Parameters
    ----------
    resolution
        Resolution to test.
    seeds
        List of random seeds. Defaults to ``[0, 1, 2, 3, 4]``.

    Returns
    -------
    Result dict with pairwise ARI values and stability assessment.
    """
    from sklearn.metrics import adjusted_rand_score

    if seeds is None:
        seeds = [0, 1, 2, 3, 4]

    # Run Leiden with each seed
    all_labels = {}
    for seed in seeds:
        key = f"_seed_test_{seed}"
        sc.tl.leiden(
            adata, resolution=resolution, key_added=key,
            random_state=seed, n_iterations=-1,
        )
        all_labels[seed] = adata.obs[key].values.copy()
        # Clean up temporary column
        del adata.obs[key]

    # Compute pairwise ARI
    ari_values = []
    for i, s1 in enumerate(seeds):
        for s2 in seeds[i + 1:]:
            ari = adjusted_rand_score(all_labels[s1], all_labels[s2])
            ari_values.append({"seed_a": s1, "seed_b": s2, "ari": round(ari, 4)})

    mean_ari = float(np.mean([x["ari"] for x in ari_values]))
    min_ari = float(np.min([x["ari"] for x in ari_values]))
    n_clusters_per_seed = {
        seed: len(set(labels)) for seed, labels in all_labels.items()
    }

    stable = min_ari >= 0.9

    warnings = []
    if not stable:
        warnings.append(
            f"Clustering is NOT stable at resolution={resolution}. "
            f"Min ARI={min_ari:.3f} (threshold: 0.9). "
            f"Cluster counts per seed: {n_clusters_per_seed}. "
            "Consider a different resolution or inspecting the neighbor graph."
        )

    return {
        "metrics": {
            "resolution": resolution,
            "seeds": seeds,
            "mean_ari": round(mean_ari, 4),
            "min_ari": round(min_ari, 4),
            "stable": stable,
            "n_clusters_per_seed": n_clusters_per_seed,
            "pairwise_ari": ari_values,
        },
        "plots": [],
        "provenance": {
            "tool_id": "leiden_clustering",
            "parameters": {
                "resolution": resolution,
                "seeds": seeds,
                "mode": "seed_sensitivity",
            },
        },
        "warnings": warnings,
    }


def _plot_sweep_metrics(
    sweep_rows: list[dict], recommended: float, plot_dir: str
) -> list[str]:
    d = Path(plot_dir)
    d.mkdir(parents=True, exist_ok=True)

    resolutions = [r["resolution"] for r in sweep_rows]
    n_clusters = [r["n_clusters"] for r in sweep_rows]
    silhouettes = [r["silhouette"] for r in sweep_rows]
    stabilities = [r["stability_ari"] for r in sweep_rows]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Panel 1: n_clusters
    ax = axes[0]
    ax.plot(resolutions, n_clusters, "o-", markersize=6, color="steelblue")
    for r, n in zip(resolutions, n_clusters):
        ax.annotate(str(n), (r, n), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=9)
    ax.axvline(recommended, color="red", linestyle="--", alpha=0.5)
    ax.set_xlabel("Resolution")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Cluster Count")

    # Panel 2: silhouette
    ax = axes[1]
    sil_valid = [(r, s) for r, s in zip(resolutions, silhouettes) if s is not None]
    if sil_valid:
        ax.plot([x[0] for x in sil_valid], [x[1] for x in sil_valid],
                "o-", markersize=6, color="forestgreen")
    ax.axvline(recommended, color="red", linestyle="--", alpha=0.5)
    ax.set_xlabel("Resolution")
    ax.set_ylabel("Silhouette Score")
    ax.set_title("Cluster Separation")

    # Panel 3: stability
    ax = axes[2]
    stab_valid = [(r, s) for r, s in zip(resolutions, stabilities) if s is not None]
    if stab_valid:
        ax.plot([x[0] for x in stab_valid], [x[1] for x in stab_valid],
                "o-", markersize=6, color="darkorange")
    ax.axhline(0.9, color="grey", linestyle=":", alpha=0.7, label="ARI = 0.9 threshold")
    ax.axvline(recommended, color="red", linestyle="--", alpha=0.5, label=f"recommended ({recommended})")
    ax.set_xlabel("Resolution")
    ax.set_ylabel("Mean ARI (stability)")
    ax.set_title("Seed Stability")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=8)

    fig.suptitle("Leiden Resolution Sweep", fontsize=14)
    fig.tight_layout()

    path = str(d / "leiden_resolution_sweep.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
