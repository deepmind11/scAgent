"""Clustering eval: pericyte-CAF adjacency in PCA KNN space."""

from __future__ import annotations
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse import issparse


def handle_clustering(adata: ad.AnnData, task_prompt: str) -> dict:
    """Compute relative PCA-KNN adjacency of pericytes to CAFs.

    The data already has cell_type labels. We:
    1. Preprocess (normalize, HVG, scale, PCA)
    2. Build KNN in PCA space
    3. For pericyte cells, count KNN edges to each cell type
    4. Compute: CAF fraction / max(other non-pericyte fractions)
    """
    cell_type_col = "cell_type"
    if cell_type_col not in adata.obs.columns:
        raise ValueError("cell_type column not found in adata.obs")

    # Preprocessing
    if "raw_counts" in adata.layers:
        adata.X = adata.layers["raw_counts"].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50, random_state=0)

    # Build KNN graph in PCA space
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50, random_state=0)

    # Get the adjacency matrix
    conn = adata.obsp["connectivities"]
    if issparse(conn):
        conn = conn.toarray()

    cell_types = adata.obs[cell_type_col].values
    unique_types = sorted(set(cell_types))

    # Find pericyte cells
    pericyte_mask = cell_types == "Pericytes"
    pericyte_indices = np.where(pericyte_mask)[0]

    if len(pericyte_indices) == 0:
        raise ValueError("No Pericytes found in cell_type labels")

    # For each pericyte, count neighbor types (excluding self-type)
    type_fractions = {}
    for ct in unique_types:
        if ct == "Pericytes":
            continue
        ct_mask = cell_types == ct
        ct_indices = np.where(ct_mask)[0]
        # Sum of adjacency weights from pericytes to this cell type
        adjacency_sum = conn[np.ix_(pericyte_indices, ct_indices)].sum()
        type_fractions[ct] = float(adjacency_sum)

    # Normalize to fractions
    total = sum(type_fractions.values())
    if total > 0:
        type_fractions = {k: v / total for k, v in type_fractions.items()}

    # CAF fraction
    caf_fraction = type_fractions.get("CAFs", 0.0)

    # Max non-CAF, non-pericyte fraction
    other_fractions = {k: v for k, v in type_fractions.items() if k != "CAFs"}
    max_other = max(other_fractions.values()) if other_fractions else 0.0

    # Relative adjacency score
    ratio = caf_fraction / max_other if max_other > 0 else float("inf")

    return {"pericyte_caf_relative_pca_knn_adjacency": round(ratio, 4)}
