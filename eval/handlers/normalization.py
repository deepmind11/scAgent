"""Normalization eval: log-normalize + z-score, report specific gene value."""

from __future__ import annotations
import numpy as np
import anndata as ad
import scanpy as sc


def handle_normalization(adata: ad.AnnData, task_prompt: str) -> dict:
    """Log-normalize, z-score scale, and report a specific gene's value.

    Task: find cell with highest raw count for 'Mrc1', report its
    final normalized+scaled value.
    """
    # Step 1: Find cell with highest raw count for Mrc1 BEFORE normalizing
    gene = "Mrc1"
    if gene not in adata.var_names:
        raise ValueError(f"Gene '{gene}' not found in data")

    gene_idx = list(adata.var_names).index(gene)
    raw_expression = np.asarray(adata.X[:, gene_idx]).flatten()
    max_cell_idx = int(np.argmax(raw_expression))

    # Step 2: Normalize (CP10K + log1p)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Step 3: Z-score scale
    sc.pp.scale(adata)

    # Step 4: Get the final value for that cell
    final_value = float(adata.X[max_cell_idx, gene_idx])

    return {"gene_value": round(final_value, 2)}
