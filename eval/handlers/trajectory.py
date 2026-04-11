"""Trajectory eval: infer trajectory, find terminal markers."""

from __future__ import annotations
import warnings
import numpy as np
import anndata as ad
import scanpy as sc

warnings.filterwarnings("ignore", message=".*default backend for leiden.*")


def handle_trajectory(adata: ad.AnnData, task_prompt: str) -> dict:
    """Infer trajectory and find markers distinguishing the two terminal groups.

    Strategy:
    1. Standard preprocessing
    2. Diffusion map for pseudotime
    3. Identify cells at the two extremes of the trajectory
    4. DE between the two terminal groups
    5. Return top marker genes
    """
    # Preprocess
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=30, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=0)

    # Diffusion map for pseudotime
    sc.tl.diffmap(adata)

    # Use first diffusion component as trajectory axis
    dc1 = adata.obsm["X_diffmap"][:, 0]

    # Identify terminal groups: bottom 15% and top 15%
    low_thresh = np.percentile(dc1, 15)
    high_thresh = np.percentile(dc1, 85)

    adata.obs["terminal_group"] = "middle"
    adata.obs.loc[dc1 <= low_thresh, "terminal_group"] = "terminal_A"
    adata.obs.loc[dc1 >= high_thresh, "terminal_group"] = "terminal_B"
    adata.obs["terminal_group"] = adata.obs["terminal_group"].astype("category")

    # DE between the two terminal groups
    terminal_mask = adata.obs["terminal_group"].isin(["terminal_A", "terminal_B"])
    adata_terminal = adata[terminal_mask, :].copy()

    sc.tl.rank_genes_groups(
        adata_terminal,
        groupby="terminal_group",
        method="wilcoxon",
        use_raw=True,
        n_genes=50,
    )

    # Collect top genes from both terminal groups
    top_genes = []
    for group in ["terminal_A", "terminal_B"]:
        try:
            genes = list(adata_terminal.uns["rank_genes_groups"]["names"][group][:25])
            top_genes.extend(genes)
        except (KeyError, IndexError):
            pass

    # Deduplicate while preserving order
    seen = set()
    unique_genes = []
    for g in top_genes:
        if g not in seen:
            seen.add(g)
            unique_genes.append(g)

    return {"top_marker_genes": unique_genes}
