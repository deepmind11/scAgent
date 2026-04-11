"""Highly variable gene selection."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def select_hvg(
    adata: ad.AnnData,
    *,
    n_top_genes: int = 2000,
    flavor: str = "seurat",
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

    Parameters
    ----------
    n_top_genes
        Number of HVGs to select. 2000 is standard; use 3000–5000 for
        heterogeneous tissues with many cell types.
    flavor
        HVG selection method. ``'seurat'`` works on log-normalized data
        (our default pipeline). ``'seurat_v3'`` requires raw counts.
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
                "flavor": flavor,
                "batch_key": batch_key,
                "min_mean": min_mean,
                "max_mean": max_mean,
                "min_disp": min_disp,
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
