"""Highly variable gene selection."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def _get_raw_counts(adata: ad.AnnData) -> np.ndarray | None:
    """Find raw counts in the AnnData, checking common locations."""
    import scipy.sparse as sp

    def _is_counts(x) -> bool:
        """Check if a matrix looks like raw integer counts."""
        if sp.issparse(x):
            data = x.data[:1000] if len(x.data) > 1000 else x.data
        else:
            data = np.asarray(x).flat[:1000]
        return np.allclose(data, np.round(data)) and np.all(data >= 0)

    # Check named layers first
    for key in ("counts", "raw_counts"):
        if key in adata.layers and _is_counts(adata.layers[key]):
            return adata.layers[key]

    # Check .X itself
    if _is_counts(adata.X):
        return adata.X

    # Check .raw
    if adata.raw is not None and _is_counts(adata.raw.X):
        return adata.raw.X

    return None


def select_hvg(
    adata: ad.AnnData,
    *,
    n_top_genes: int = 2000,
    flavor: str = "seurat_v3",
    batch_key: str | None = None,
    min_mean: float = 0.0125,
    max_mean: float = 3.0,
    min_disp: float = 0.5,
    plot_dir: str | None = None,
    checkpoint_dir: str | None = None,
) -> dict:
    """Select highly variable genes. Modifies *adata* in place.

    Marks genes in ``adata.var['highly_variable']`` but does NOT subset
    ``adata.X``. PCA uses ``use_highly_variable=True`` to restrict
    automatically, keeping all genes available in ``adata.raw``.

    By default uses ``seurat_v3`` on raw counts, which produces more
    biologically meaningful HVGs that recover known cell-type markers
    (Hafemeister & Satija 2019). Falls back to ``seurat`` on log-normalized
    data if raw counts are not available.

    Parameters
    ----------
    n_top_genes
        Number of HVGs to select. 2000 is standard; use 3000–5000 for
        heterogeneous tissues with many cell types.
    flavor
        HVG selection method. ``'seurat_v3'`` (default, recommended) uses
        variance-stabilizing transformation on raw counts. ``'seurat'``
        works on log-normalized data but is biased toward lowly-expressed
        genes.
    batch_key
        If set, select HVGs per batch and take the union. Use when
        integrating multi-batch data.
    plot_dir
        Directory for the mean-vs-dispersion plot.
    checkpoint_dir
        Save checkpoint after selection.

    Returns
    -------
    Result dict with metrics, plots, provenance.
    """
    actual_flavor = flavor
    used_raw_from = None

    if flavor == "seurat_v3":
        raw = _get_raw_counts(adata)
        if raw is not None:
            # Temporarily swap .X to raw counts for seurat_v3
            original_X = adata.X.copy()
            adata.X = raw
            used_raw_from = next(
                (k for k in ("counts", "raw_counts") if k in adata.layers),
                ".X" if raw is adata.X else ".raw",
            )
            try:
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=n_top_genes,
                    flavor="seurat_v3",
                    batch_key=batch_key,
                )
            finally:
                # Restore original .X (log-normalized) so downstream PCA works
                adata.X = original_X
        else:
            # No raw counts found — fall back to seurat
            actual_flavor = "seurat"
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                flavor="seurat",
                batch_key=batch_key,
                min_mean=min_mean,
                max_mean=max_mean,
                min_disp=min_disp,
            )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor=flavor,
            batch_key=batch_key,
            min_mean=min_mean,
            max_mean=max_mean,
            min_disp=min_disp,
        )

    n_hvgs = int(adata.var["highly_variable"].sum())

    plots = []
    if plot_dir is not None:
        plots = _plot_hvg(adata, plot_dir)

    warnings = []
    if n_hvgs < 500:
        warnings.append(
            f"Only {n_hvgs} HVGs selected — very low. Consider increasing "
            "n_top_genes or lowering min_disp."
        )
    if n_hvgs > 8000:
        warnings.append(
            f"{n_hvgs} HVGs selected — unusually high. Results may include "
            "noisy genes. Consider lowering n_top_genes."
        )

    result = {
        "metrics": {
            "n_hvgs": n_hvgs,
            "n_total_genes": adata.n_vars,
            "pct_hvg": round(n_hvgs / adata.n_vars * 100, 1),
        },
        "plots": plots,
        "provenance": {
            "tool_id": "highly_variable_genes",
            "parameters": {
                "n_top_genes": n_top_genes,
                "flavor": actual_flavor,
                "flavor_requested": flavor,
                "raw_counts_source": used_raw_from,
                "batch_key": batch_key,
                "min_mean": min_mean if actual_flavor != "seurat_v3" else None,
                "max_mean": max_mean if actual_flavor != "seurat_v3" else None,
                "min_disp": min_disp if actual_flavor != "seurat_v3" else None,
            },
            "n_hvgs_selected": n_hvgs,
        },
        "warnings": warnings,
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_hvg")

    return result


def _plot_hvg(adata: ad.AnnData, plot_dir: str) -> list[str]:
    d = Path(plot_dir)
    d.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    hvg = adata.var["highly_variable"]
    means = adata.var["means"]
    dispersions_norm = adata.var["dispersions_norm"]

    ax.scatter(means[~hvg], dispersions_norm[~hvg], s=1, alpha=0.3, c="grey", label="Other")
    ax.scatter(means[hvg], dispersions_norm[hvg], s=1, alpha=0.5, c="black", label="HVG")
    ax.set_xlabel("Mean expression")
    ax.set_ylabel("Normalized dispersion")
    ax.set_title(f"Highly Variable Genes ({hvg.sum()} selected)")
    ax.legend(markerscale=5)
    fig.tight_layout()

    path = str(d / "hvg_selection.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
