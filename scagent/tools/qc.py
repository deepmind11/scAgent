"""Quality control: compute metrics, generate plots, filter cells and genes."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")  # non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def calculate_qc_metrics(
    adata: ad.AnnData,
    *,
    species: str | None = None,
    plot_dir: str | None = None,
) -> dict:
    """Compute QC metrics and generate diagnostic plots. Modifies *adata* in place.

    Adds to ``adata.obs``: n_genes_by_counts, total_counts, pct_counts_mt.
    Adds to ``adata.var``: mt (bool), n_cells_by_counts.
    Stores detected species in ``adata.uns['species']`` if not already set.

    Parameters
    ----------
    species
        ``'human'`` or ``'mouse'``. Auto-detected from gene names if *None*.
    plot_dir
        Directory to save diagnostic plots. No plots if *None*.

    Returns
    -------
    Result dict with metrics (medians, MAD-based threshold recommendations),
    plot paths, and warnings.
    """
    # Detect / confirm species
    if species is None:
        species = adata.uns.get("species", _detect_species_from_genes(adata))
    adata.uns["species"] = species

    # Annotate mitochondrial genes
    mt_prefix = "MT-" if species == "human" else "mt-"
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    n_mt = adata.var["mt"].sum()

    # Compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # Compute summary statistics
    obs = adata.obs
    metrics = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_mito_genes": int(n_mt),
        "species": species,
        "median_genes_per_cell": float(np.median(obs["n_genes_by_counts"])),
        "median_counts_per_cell": float(np.median(obs["total_counts"])),
        "median_pct_mito": float(np.median(obs["pct_counts_mt"])),
        "mean_pct_mito": float(np.mean(obs["pct_counts_mt"])),
    }

    # MAD-based threshold recommendations
    thresholds = _compute_mad_thresholds(obs)
    metrics["recommended_thresholds"] = thresholds

    # Generate plots
    plots = []
    if plot_dir is not None:
        plots = _generate_qc_plots(adata, plot_dir)

    warnings = []
    if n_mt == 0:
        warnings.append(
            f"No mitochondrial genes found with prefix '{mt_prefix}'. "
            "Check if the species is correct or if gene names use a different convention."
        )

    return {
        "metrics": metrics,
        "plots": plots,
        "provenance": {
            "tool_id": "calculate_qc_metrics",
            "parameters": {"species": species, "mt_prefix": mt_prefix},
        },
        "warnings": warnings,
    }


def filter_cells(
    adata: ad.AnnData,
    *,
    min_genes: int = 200,
    max_genes: int = 5000,
    max_pct_mito: float = 10.0,
    min_counts: int | None = None,
    checkpoint_dir: str | None = None,
) -> tuple[ad.AnnData, dict]:
    """Filter low-quality cells based on QC metrics.

    Requires :func:`calculate_qc_metrics` to have been run first.

    Returns
    -------
    (filtered_adata, result) — a new AnnData object (view materialised to copy).
    """
    if "n_genes_by_counts" not in adata.obs.columns:
        raise ValueError(
            "QC metrics not found in adata.obs. Run calculate_qc_metrics() first."
        )

    cells_before = adata.n_obs

    # Build boolean mask
    keep = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["n_genes_by_counts"] <= max_genes)
        & (adata.obs["pct_counts_mt"] <= max_pct_mito)
    )
    if min_counts is not None:
        keep = keep & (adata.obs["total_counts"] >= min_counts)

    adata_filtered = adata[keep, :].copy()
    cells_after = adata_filtered.n_obs
    cells_removed = cells_before - cells_after

    result = {
        "metrics": {
            "cells_before": cells_before,
            "cells_after": cells_after,
            "cells_removed": cells_removed,
            "pct_removed": round(cells_removed / cells_before * 100, 1),
        },
        "plots": [],
        "provenance": {
            "tool_id": "filter_cells",
            "parameters": {
                "min_genes": min_genes,
                "max_genes": max_genes,
                "max_pct_mito": max_pct_mito,
                "min_counts": min_counts,
            },
            "cells_before": cells_before,
            "cells_after": cells_after,
        },
        "warnings": [],
    }

    if cells_after / cells_before < 0.5:
        result["warnings"].append(
            f"More than 50% of cells removed ({result['metrics']['pct_removed']}%). "
            "Consider relaxing thresholds."
        )

    if checkpoint_dir is not None:
        _save_checkpoint(adata_filtered, checkpoint_dir, "after_filter_cells")

    return adata_filtered, result


def filter_genes(
    adata: ad.AnnData,
    *,
    min_cells: int = 3,
    checkpoint_dir: str | None = None,
) -> tuple[ad.AnnData, dict]:
    """Filter rarely-expressed genes.

    Returns
    -------
    (filtered_adata, result) — a new AnnData object.
    """
    genes_before = adata.n_vars

    # sc.pp.filter_genes modifies in place and returns None,
    # but we want to return a clean new object for consistency.
    sc.pp.filter_genes(adata, min_cells=min_cells)
    genes_after = adata.n_vars

    result = {
        "metrics": {
            "genes_before": genes_before,
            "genes_after": genes_after,
            "genes_removed": genes_before - genes_after,
        },
        "plots": [],
        "provenance": {
            "tool_id": "filter_genes",
            "parameters": {"min_cells": min_cells},
            "genes_before": genes_before,
            "genes_after": genes_after,
        },
        "warnings": [],
    }

    if genes_after / genes_before < 0.5:
        result["warnings"].append(
            f"More than 50% of genes removed. Check if min_cells={min_cells} is too aggressive."
        )

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_filter_genes")

    return adata, result


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _detect_species_from_genes(adata: ad.AnnData) -> str:
    sample = adata.var_names[:200].tolist()
    n_upper = sum(1 for g in sample if g == g.upper())
    return "human" if n_upper / len(sample) > 0.7 else "mouse"


def _compute_mad_thresholds(obs) -> dict:
    """Compute MAD-based outlier thresholds for QC metrics."""
    def _mad_bounds(series, n_mads: float = 5.0, log: bool = False):
        values = series.values.copy().astype(float)
        if log:
            values = np.log1p(values)
        med = np.median(values)
        mad = np.median(np.abs(values - med))
        lower = med - n_mads * mad
        upper = med + n_mads * mad
        if log:
            lower = np.expm1(lower)
            upper = np.expm1(upper)
        return max(0, float(lower)), float(upper)

    min_genes, max_genes = _mad_bounds(obs["n_genes_by_counts"], n_mads=5.0, log=True)
    min_counts, max_counts = _mad_bounds(obs["total_counts"], n_mads=5.0, log=True)

    # For mito, only upper bound matters
    mito_values = obs["pct_counts_mt"].values.astype(float)
    mito_med = np.median(mito_values)
    mito_mad = np.median(np.abs(mito_values - mito_med))
    max_pct_mito_mad = float(mito_med + 5.0 * mito_mad)

    return {
        "min_genes": int(round(min_genes)),
        "max_genes": int(round(max_genes)),
        "min_counts": int(round(min_counts)),
        "max_counts": int(round(max_counts)),
        "max_pct_mito_mad": round(max_pct_mito_mad, 1),
    }


def _generate_qc_plots(adata: ad.AnnData, plot_dir: str) -> list[str]:
    """Generate and save QC diagnostic plots."""
    d = Path(plot_dir)
    d.mkdir(parents=True, exist_ok=True)
    plots = []

    # 1. Violin plots — n_genes, total_counts, pct_mito
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, metric, title in zip(
        axes,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        ["Genes per cell", "UMI counts per cell", "% mitochondrial"],
    ):
        ax.violinplot(adata.obs[metric].values, showmedians=True)
        ax.set_title(title)
        ax.set_ylabel(metric)
        ax.set_xticks([])
    fig.suptitle("QC Metrics Distribution", fontsize=14)
    fig.tight_layout()
    path = str(d / "qc_violin.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    plots.append(path)

    # 2. Scatter — total_counts vs n_genes, colored by pct_mito
    fig, ax = plt.subplots(figsize=(7, 5))
    scatter = ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        cmap="RdYlGn_r",
        s=1,
        alpha=0.5,
    )
    ax.set_xlabel("Total UMI counts")
    ax.set_ylabel("Genes detected")
    ax.set_title("Counts vs Genes (color = % mito)")
    fig.colorbar(scatter, ax=ax, label="% mitochondrial")
    fig.tight_layout()
    path = str(d / "qc_scatter.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    plots.append(path)

    # 3. Mito histogram
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(adata.obs["pct_counts_mt"], bins=100, edgecolor="none")
    ax.set_xlabel("% mitochondrial counts")
    ax.set_ylabel("Number of cells")
    ax.set_title("Mitochondrial Content Distribution")
    # Add median line
    med = np.median(adata.obs["pct_counts_mt"])
    ax.axvline(med, color="blue", linestyle="--", label=f"median = {med:.1f}%")
    ax.legend()
    fig.tight_layout()
    path = str(d / "qc_mito_hist.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    plots.append(path)

    return plots


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
