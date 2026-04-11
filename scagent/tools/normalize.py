"""Normalization: CP10K + log1p with raw count preservation."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc
from scipy.sparse import issparse


def log_normalize(
    adata: ad.AnnData,
    *,
    target_sum: float = 1e4,
    exclude_highly_expressed: bool = False,
    checkpoint_dir: str | None = None,
) -> dict:
    """Log-normalize counts (CP10K + log1p). Modifies *adata* in place.

    Order of operations:
    1. Copy raw counts to ``adata.layers['counts']``
    2. ``sc.pp.normalize_total(adata, target_sum=target_sum)``
    3. ``sc.pp.log1p(adata)``
    4. ``adata.raw = adata.copy()`` (freezes full gene set, normalized)

    This ensures:
    - ``adata.layers['counts']`` = raw integer UMIs (for DESeq2/edgeR pseudobulk)
    - ``adata.raw.X`` = log-normalized, all genes (for Wilcoxon markers, plotting)
    - ``adata.X`` = log-normalized (will be subset to HVGs later for PCA)

    Parameters
    ----------
    target_sum
        Normalize each cell to this total count before log-transform.
        1e4 (CP10K) is the standard for scRNA-seq.
    exclude_highly_expressed
        If *True*, exclude very highly expressed genes from the
        normalization factor computation.
    checkpoint_dir
        If provided, save checkpoint after normalization.

    Returns
    -------
    Result dict with metrics, provenance, and warnings.
    """
    # 1. Preserve raw counts
    adata.layers["counts"] = adata.X.copy()

    # 2. Normalize
    sc.pp.normalize_total(
        adata,
        target_sum=target_sum,
        exclude_highly_expressed=exclude_highly_expressed,
    )

    # 3. Log-transform
    sc.pp.log1p(adata)

    # 4. Freeze as raw (full gene set, normalized)
    adata.raw = adata.copy()

    # Validation
    X = adata.X
    if issparse(X):
        max_val = float(X.max())
    else:
        max_val = float(np.max(X))

    has_nan = False
    if issparse(X):
        has_nan = bool(np.isnan(X.data).any())
    else:
        has_nan = bool(np.isnan(X).any())

    # Check that counts layer has integers
    counts_layer = adata.layers["counts"]
    if issparse(counts_layer):
        sample = counts_layer.data[:1000]
    else:
        sample = counts_layer.flat[:1000]
    counts_are_integers = np.allclose(sample, np.round(sample))

    warnings_list = []
    if max_val >= 15:
        warnings_list.append(
            f"Max normalized value is {max_val:.1f} (expected <15). "
            "Data may have been pre-normalized or contain extreme outliers."
        )
    if has_nan:
        warnings_list.append("NaN values detected after normalization.")
    if not counts_are_integers:
        warnings_list.append(
            "Raw counts layer does not contain integers. "
            "Data may have been pre-normalized before this step."
        )

    result = {
        "metrics": {
            "max_normalized_value": round(max_val, 2),
            "has_nan": has_nan,
            "raw_counts_preserved": "counts" in adata.layers,
            "raw_frozen": adata.raw is not None,
            "counts_are_integers": counts_are_integers,
        },
        "plots": [],
        "provenance": {
            "tool_id": "log_normalize",
            "parameters": {
                "target_sum": target_sum,
                "exclude_highly_expressed": exclude_highly_expressed,
            },
        },
        "warnings": warnings_list,
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_normalize")

    return result


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
