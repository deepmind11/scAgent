"""Batch correction / integration methods."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import scanpy as sc


def run_harmony(
    adata: ad.AnnData,
    *,
    key: str,
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
    max_iter_harmony: int = 10,
    random_state: int = 0,
    checkpoint_dir: str | None = None,
) -> dict:
    """Run Harmony batch correction. Modifies *adata* in place.

    Adds ``adata.obsm['X_pca_harmony']`` (or *adjusted_basis*).

    Only needed for multi-sample or multi-batch experiments. For
    single-sample data, skip and use ``X_pca`` directly for neighbors.

    Parameters
    ----------
    key
        obs column containing batch labels (e.g., ``'sample'``,
        ``'donor'``, ``'batch'``).
    basis
        Embedding to correct. Almost always ``'X_pca'``.
    adjusted_basis
        Key for the corrected embedding.
    max_iter_harmony
        Maximum Harmony iterations. 10 is usually sufficient.
    random_state
        Random seed for Harmony's iterative optimization.
    checkpoint_dir
        Save checkpoint after integration.

    Returns
    -------
    Result dict with metrics, provenance.
    """
    if key not in adata.obs.columns:
        raise ValueError(
            f"Batch key '{key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    if basis not in adata.obsm:
        raise ValueError(
            f"Basis '{basis}' not found in adata.obsm. Run PCA first."
        )

    n_batches = adata.obs[key].nunique()
    batch_sizes = adata.obs[key].value_counts()

    sc.external.pp.harmony_integrate(
        adata,
        key=key,
        basis=basis,
        adjusted_basis=adjusted_basis,
        max_iter_harmony=max_iter_harmony,
        random_state=random_state,
    )

    result = {
        "metrics": {
            "n_batches": n_batches,
            "batch_sizes": {str(k): int(v) for k, v in batch_sizes.items()},
            "adjusted_basis": adjusted_basis,
        },
        "plots": [],
        "provenance": {
            "tool_id": "harmony",
            "parameters": {
                "key": key,
                "basis": basis,
                "adjusted_basis": adjusted_basis,
                "max_iter_harmony": max_iter_harmony,
                "random_state": random_state,
            },
            "n_batches": n_batches,
        },
        "warnings": [],
    }

    if n_batches < 2:
        result["warnings"].append(
            f"Only {n_batches} batch(es) found in '{key}'. "
            "Integration is unnecessary for single-batch data."
        )

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_integration")

    return result


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
