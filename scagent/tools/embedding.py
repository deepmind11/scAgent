"""UMAP embedding for visualization."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc


def run_umap(
    adata: ad.AnnData,
    *,
    min_dist: float = 0.5,
    spread: float = 1.0,
    random_state: int = 0,
    checkpoint_dir: str | None = None,
) -> dict:
    """Compute UMAP embedding for visualization. Modifies *adata* in place.

    Adds ``adata.obsm['X_umap']``.

    **UMAP is for VISUALIZATION ONLY.** Never cluster on UMAP coordinates.
    Never interpret distances between clusters as biologically meaningful.

    Parameters
    ----------
    min_dist
        Controls how tightly points pack together. Lower → tighter clusters
        visually. 0.1–0.3 often looks better than the default 0.5.
    spread
        Scale of the embedding. Usually left at 1.0.
    random_state
        Random seed. UMAP is stochastic — different seeds give different
        layouts, but topology should be stable.
    checkpoint_dir
        Save checkpoint after UMAP.

    Returns
    -------
    Result dict with metrics, provenance.
    """
    sc.tl.umap(
        adata,
        min_dist=min_dist,
        spread=spread,
        random_state=random_state,
    )

    embedding = adata.obsm["X_umap"]
    has_nan = bool(np.isnan(embedding).any())

    warnings = []
    if has_nan:
        warnings.append("NaN values in UMAP embedding. Check neighbor graph.")

    result = {
        "metrics": {
            "has_nan": has_nan,
            "x_range": [float(embedding[:, 0].min()), float(embedding[:, 0].max())],
            "y_range": [float(embedding[:, 1].min()), float(embedding[:, 1].max())],
        },
        "plots": [],
        "provenance": {
            "tool_id": "umap",
            "parameters": {
                "min_dist": min_dist,
                "spread": spread,
                "random_state": random_state,
            },
        },
        "warnings": warnings,
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_umap")

    return result


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
