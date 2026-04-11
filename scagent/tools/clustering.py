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
    plot_dir: str | None = None,
) -> dict:
    """Run Leiden at multiple resolutions and report cluster counts.

    Each resolution gets its own key: ``adata.obs['leiden_{res}']``.
    Does NOT pick a resolution — that is the user's decision.

    Parameters
    ----------
    resolutions
        List of resolutions to try. Defaults to
        ``[0.3, 0.5, 0.8, 1.0, 1.5, 2.0]``.
    random_state
        Same seed for every resolution — makes comparison fair.
    plot_dir
        Directory for the resolution-vs-clusters plot.

    Returns
    -------
    Result dict with a table of resolution → n_clusters.
    """
    if resolutions is None:
        resolutions = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]

    sweep_results = {}
    for res in resolutions:
        key = f"leiden_{res}"
        sc.tl.leiden(
            adata,
            resolution=res,
            key_added=key,
            random_state=random_state,
            n_iterations=-1,
        )
        n_clusters = adata.obs[key].nunique()
        sweep_results[res] = n_clusters

    plots = []
    if plot_dir is not None:
        plots = _plot_sweep(sweep_results, plot_dir)

    return {
        "metrics": {
            "resolution_to_clusters": {str(k): v for k, v in sweep_results.items()},
        },
        "plots": plots,
        "provenance": {
            "tool_id": "leiden_clustering",
            "parameters": {
                "resolutions": resolutions,
                "random_state": random_state,
                "mode": "sweep",
            },
        },
        "warnings": [],
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


def _plot_sweep(sweep_results: dict, plot_dir: str) -> list[str]:
    d = Path(plot_dir)
    d.mkdir(parents=True, exist_ok=True)

    resolutions = sorted(sweep_results.keys())
    n_clusters = [sweep_results[r] for r in resolutions]

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(resolutions, n_clusters, "o-", markersize=6)
    for r, n in zip(resolutions, n_clusters):
        ax.annotate(str(n), (r, n), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=9)
    ax.set_xlabel("Resolution")
    ax.set_ylabel("Number of clusters")
    ax.set_title("Leiden Resolution Sweep")
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
