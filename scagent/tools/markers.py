"""Marker gene detection for cluster characterization."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def find_marker_genes(
    adata: ad.AnnData,
    *,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 100,
    use_raw: bool = True,
    corr_method: str = "benjamini-hochberg",
    plot_dir: str | None = None,
    checkpoint_dir: str | None = None,
) -> dict:
    """Find marker genes per cluster.

    Uses ``sc.tl.rank_genes_groups`` with the specified method.

    Parameters
    ----------
    groupby
        obs column with cluster or cell type labels.
    method
        Statistical test. ``'wilcoxon'`` rank-sum is the standard for
        scRNA-seq marker detection.
    n_genes
        Number of top genes to report per group.
    use_raw
        Use ``adata.raw`` for testing. Should be *True* to test on all
        genes, not just HVGs.
    corr_method
        Multiple testing correction.
    plot_dir
        Directory for marker gene plots (dotplot, heatmap).
    checkpoint_dir
        Save checkpoint after marker detection.

    Returns
    -------
    Result dict with top markers per cluster, plots, provenance.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(
            f"Column '{groupby}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        n_genes=n_genes,
        use_raw=use_raw,
        corr_method=corr_method,
    )

    # Extract top markers per cluster
    groups = adata.obs[groupby].cat.categories.tolist()
    top_markers = {}
    n_top = 5  # For the summary
    for group in groups:
        names = adata.uns["rank_genes_groups"]["names"][group][:n_top]
        scores = adata.uns["rank_genes_groups"]["scores"][group][:n_top]
        logfc = adata.uns["rank_genes_groups"]["logfoldchanges"][group][:n_top]
        pvals = adata.uns["rank_genes_groups"]["pvals_adj"][group][:n_top]
        top_markers[str(group)] = [
            {
                "gene": str(names[i]),
                "score": round(float(scores[i]), 2),
                "logfc": round(float(logfc[i]), 2),
                "pval_adj": float(pvals[i]),
            }
            for i in range(len(names))
        ]

    # Check quality: clusters with weak markers
    weak_clusters = []
    for group in groups:
        pvals = adata.uns["rank_genes_groups"]["pvals_adj"][group][:10]
        logfc = adata.uns["rank_genes_groups"]["logfoldchanges"][group][:10]
        n_sig = sum(1 for p, l in zip(pvals, logfc) if p < 0.05 and l > 1.0)
        if n_sig < 3:
            weak_clusters.append(str(group))

    plots = []
    if plot_dir is not None:
        plots = _plot_markers(adata, groupby, n_top, plot_dir)

    warnings = []
    if weak_clusters:
        warnings.append(
            f"Clusters with <3 significant markers (padj<0.05, logFC>1): "
            f"{weak_clusters}. These clusters may be poorly defined — "
            "consider merging them or adjusting resolution."
        )

    result = {
        "metrics": {
            "n_clusters": len(groups),
            "groupby": groupby,
            "method": method,
            "top_markers": top_markers,
            "weak_clusters": weak_clusters,
        },
        "plots": plots,
        "provenance": {
            "tool_id": "wilcoxon_markers",
            "parameters": {
                "groupby": groupby,
                "method": method,
                "n_genes": n_genes,
                "use_raw": use_raw,
                "corr_method": corr_method,
            },
            "n_clusters_tested": len(groups),
        },
        "warnings": warnings,
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_markers")

    return result


def _plot_markers(
    adata: ad.AnnData, groupby: str, n_top: int, plot_dir: str
) -> list[str]:
    d = Path(plot_dir)
    d.mkdir(parents=True, exist_ok=True)
    plots = []

    # Dotplot of top markers
    try:
        fig = sc.pl.rank_genes_groups_dotplot(
            adata, n_genes=n_top, groupby=groupby, show=False, return_fig=True,
        )
        path = str(d / "marker_dotplot.png")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        plots.append(path)
    except Exception:
        # Fall back if dotplot fails (some scanpy versions)
        pass

    return plots


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
