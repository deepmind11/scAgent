"""Principal component analysis."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def run_pca(
    adata: ad.AnnData,
    *,
    n_comps: int = 50,
    use_highly_variable: bool = True,
    scale: bool = True,
    max_value: float | None = 10.0,
    svd_solver: str = "arpack",
    random_state: int = 0,
    plot_dir: str | None = None,
    checkpoint_dir: str | None = None,
) -> dict:
    """Z-score scale and run PCA. Modifies *adata* in place.

    The standard Scanpy pipeline is::

        normalize_total → log1p → [set adata.raw] → HVG → **scale** → PCA

    This function handles the last two steps: z-score scaling (each gene
    to mean 0, variance 1) followed by PCA. Scaling is required before
    PCA so that highly-expressed genes don't dominate the components.

    After scaling, ``adata.X`` contains z-scores (not log-counts).
    ``adata.raw`` (set during normalization) still holds the log-normalized
    values for all genes — used for plotting and marker detection.

    Adds ``adata.obsm['X_pca']`` and ``adata.uns['pca']``.

    Parameters
    ----------
    n_comps
        Number of principal components. 30–50 is standard.
    use_highly_variable
        Restrict PCA to highly variable genes. Should be *True* after
        HVG selection.
    scale
        Z-score scale genes before PCA (mean=0, var=1). Standard and
        recommended. Set *False* only if data is already scaled.
    max_value
        Clip scaled values to this maximum (avoids outlier genes
        dominating PCA). 10 is the Scanpy default. *None* to disable.
    svd_solver
        SVD solver. ``'arpack'`` is default and deterministic.
        ``'randomized'`` is faster for large datasets (>50K cells)
        but depends on ``random_state``.
    random_state
        Random seed. Relevant for ``svd_solver='randomized'``.
    plot_dir
        Directory for the elbow plot (variance ratio).
    checkpoint_dir
        Save checkpoint after PCA.

    Returns
    -------
    Result dict with metrics, plots, provenance.
    """
    # Step 1: Z-score scaling (standard pre-PCA step)
    if scale:
        sc.pp.scale(adata, max_value=max_value)

    # Step 2: PCA
    sc.pp.pca(
        adata,
        n_comps=n_comps,
        use_highly_variable=use_highly_variable,
        svd_solver=svd_solver,
        random_state=random_state,
    )

    variance_ratio = adata.uns["pca"]["variance_ratio"]
    total_var = float(np.sum(variance_ratio))
    pc1_var = float(variance_ratio[0])

    plots = []
    if plot_dir is not None:
        plots = _plot_elbow(variance_ratio, plot_dir)

    warnings = []
    if pc1_var > 0.5:
        warnings.append(
            f"PC1 explains {pc1_var:.1%} of variance — unusually high. "
            "This often indicates a dominant technical artifact (e.g., total "
            "counts, cell cycle). Inspect PC1 loadings to identify the source."
        )

    result = {
        "metrics": {
            "n_comps": n_comps,
            "total_variance_explained": round(total_var, 4),
            "pc1_variance": round(pc1_var, 4),
            "variance_ratio_top10": [round(float(v), 4) for v in variance_ratio[:10]],
        },
        "plots": plots,
        "provenance": {
            "tool_id": "pca",
            "parameters": {
                "n_comps": n_comps,
                "use_highly_variable": use_highly_variable,
                "scale": scale,
                "max_value": max_value,
                "svd_solver": svd_solver,
                "random_state": random_state,
            },
            "total_variance_explained": round(total_var, 4),
        },
        "warnings": warnings,
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_pca")

    return result


def _plot_elbow(variance_ratio: np.ndarray, plot_dir: str) -> list[str]:
    d = Path(plot_dir)
    d.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 4))
    n = len(variance_ratio)
    ax.plot(range(1, n + 1), variance_ratio, "o-", markersize=3)
    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Variance Ratio")
    ax.set_title("PCA Elbow Plot")

    # Mark cumulative 90% variance
    cumvar = np.cumsum(variance_ratio)
    idx_90 = int(np.searchsorted(cumvar, 0.9)) + 1
    if idx_90 <= n:
        ax.axvline(idx_90, color="red", linestyle="--", alpha=0.7,
                   label=f"90% variance at PC {idx_90}")
        ax.legend()

    fig.tight_layout()
    path = str(d / "pca_elbow.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
