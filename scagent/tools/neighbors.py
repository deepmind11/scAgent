"""K-nearest neighbor graph construction."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import scanpy as sc


def compute_neighbors(
    adata: ad.AnnData,
    *,
    n_neighbors: int = 15,
    n_pcs: int = 50,
    metric: str = "cosine",
    use_rep: str = "X_pca",
    method: str = "umap",
    random_state: int = 0,
    checkpoint_dir: str | None = None,
) -> dict:
    """Compute k-nearest neighbor graph. Modifies *adata* in place.

    Adds ``adata.obsp['connectivities']``, ``adata.obsp['distances']``,
    and ``adata.uns['neighbors']``.

    Parameters
    ----------
    n_neighbors
        Number of nearest neighbors. 10–15 for fine-grained, 20–30 for
        coarser structure. Sciaraffa et al. 2025 found k=10 optimal with
        high-resolution clustering.
    n_pcs
        Number of PCs to use. Should be ≤ n_comps from PCA.
    metric
        Distance metric. ``'cosine'`` preferred (Sciaraffa et al. 2025).
    use_rep
        Embedding to use. ``'X_pca'`` for standard, ``'X_pca_harmony'``
        after Harmony, ``'X_scVI'`` after scVI.
    method
        Method for computing connectivities. ``'umap'`` preferred
        (Sciaraffa et al. 2025).
    random_state
        Random seed for the approximate nearest neighbor search.
    checkpoint_dir
        Save checkpoint after computing neighbors.

    Returns
    -------
    Result dict with metrics, provenance.
    """
    if use_rep not in adata.obsm:
        available = list(adata.obsm.keys())
        raise ValueError(
            f"Representation '{use_rep}' not found in adata.obsm. "
            f"Available: {available}"
        )

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        metric=metric,
        use_rep=use_rep,
        method=method,
        random_state=random_state,
    )

    result = {
        "metrics": {
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "metric": metric,
            "use_rep": use_rep,
        },
        "plots": [],
        "provenance": {
            "tool_id": "neighbor_graph",
            "parameters": {
                "n_neighbors": n_neighbors,
                "n_pcs": n_pcs,
                "metric": metric,
                "use_rep": use_rep,
                "method": method,
                "random_state": random_state,
            },
        },
        "warnings": [],
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_neighbors")

    return result


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
