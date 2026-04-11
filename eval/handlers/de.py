"""DE eval: sub-cluster CAFs, find contractile subcluster markers."""

from __future__ import annotations
import warnings
import numpy as np
import anndata as ad
import scanpy as sc

warnings.filterwarnings("ignore", message=".*default backend for leiden.*")


def handle_de(adata: ad.AnnData, task_prompt: str) -> dict:
    """Cluster CAFs into 6 subclusters, find contractile markers.

    The data has cell_type labels. We:
    1. Subset to CAFs
    2. Preprocess the CAF subset
    3. Cluster into 6 subclusters
    4. Identify the contractile subcluster (highest expression of
       contractile genes like Acta2, Tagln, Myl9)
    5. Run DE: contractile vs. all other CAFs
    6. Return top 50 marker genes
    """
    cell_type_col = "cell_type"
    if cell_type_col not in adata.obs.columns:
        raise ValueError("cell_type column not found")

    # Subset to CAFs
    caf_mask = adata.obs[cell_type_col] == "CAFs"
    adata_caf = adata[caf_mask, :].copy()

    # Use raw counts if available
    if "raw_counts" in adata_caf.layers:
        adata_caf.X = adata_caf.layers["raw_counts"].copy()

    # Preprocess
    sc.pp.normalize_total(adata_caf, target_sum=1e4)
    sc.pp.log1p(adata_caf)
    adata_caf.raw = adata_caf.copy()

    sc.pp.highly_variable_genes(adata_caf, n_top_genes=2000, flavor="seurat")
    sc.pp.scale(adata_caf, max_value=10)
    sc.pp.pca(adata_caf, n_comps=30, random_state=0)
    sc.pp.neighbors(adata_caf, n_neighbors=15, random_state=0)

    # Find resolution that gives ~6 clusters
    resolution = _find_resolution_for_k(adata_caf, target_k=6, random_state=0)
    sc.tl.leiden(adata_caf, resolution=resolution, key_added="leiden",
                 random_state=0, n_iterations=-1)

    n_clusters = adata_caf.obs["leiden"].nunique()

    # Identify contractile subcluster by scoring contractile markers
    contractile_markers = ["Acta2", "Tagln", "Myl9", "Tpm1", "Tpm2", "Cnn2"]
    available = [g for g in contractile_markers if g in adata_caf.raw.var_names]

    if available:
        sc.tl.score_genes(adata_caf, gene_list=available, score_name="contractile_score")
        # Find cluster with highest mean contractile score
        cluster_scores = adata_caf.obs.groupby("leiden")["contractile_score"].mean()
        contractile_cluster = cluster_scores.idxmax()
    else:
        # Fallback: just pick cluster 0
        contractile_cluster = "0"

    # DE: contractile vs all other CAFs
    adata_caf.obs["is_contractile"] = (
        adata_caf.obs["leiden"] == contractile_cluster
    ).astype(str).astype("category")

    sc.tl.rank_genes_groups(
        adata_caf,
        groupby="is_contractile",
        groups=["True"],
        reference="False",
        method="wilcoxon",
        use_raw=True,
        n_genes=50,
    )

    # Extract top 50 genes
    top_genes = list(adata_caf.uns["rank_genes_groups"]["names"]["True"][:50])

    return {"top_marker_genes": top_genes}


def _find_resolution_for_k(adata, target_k: int, random_state: int = 0) -> float:
    """Binary search for a Leiden resolution giving ~target_k clusters."""
    lo, hi = 0.1, 3.0
    best_res = 1.0
    best_diff = float("inf")

    for _ in range(15):
        mid = (lo + hi) / 2
        sc.tl.leiden(adata, resolution=mid, key_added="_search",
                     random_state=random_state, n_iterations=-1)
        k = adata.obs["_search"].nunique()
        del adata.obs["_search"]

        diff = abs(k - target_k)
        if diff < best_diff:
            best_diff = diff
            best_res = mid

        if k < target_k:
            lo = mid
        elif k > target_k:
            hi = mid
        else:
            return mid

    return best_res
