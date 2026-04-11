"""Cell type annotation — CellTypist and marker-based manual."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def annotate_celltypist(
    adata: ad.AnnData,
    *,
    model: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
    over_clustering: str | None = None,
    checkpoint_dir: str | None = None,
) -> dict:
    """Annotate cell types using CellTypist.

    Downloads the model on first use (~50 MB).

    Parameters
    ----------
    model
        Pre-trained model name. Common choices:
        - ``'Immune_All_Low.pkl'`` — broad immune subtypes (PBMCs)
        - ``'Immune_All_High.pkl'`` — fine-grained immune subtypes
        Run ``celltypist.models.models_description()`` for full list.
    majority_voting
        If *True*, refine predictions by majority vote within clusters.
        More robust than per-cell predictions.
    over_clustering
        obs key for an over-clustering to use with majority voting.
        If *None*, CellTypist generates its own internally.
    checkpoint_dir
        Save checkpoint after annotation.

    Returns
    -------
    Result dict with cell type counts, provenance.
    """
    import celltypist
    from celltypist import models as ct_models

    # Download model if needed
    ct_models.download_models(model=model)
    ct_model = ct_models.Model.load(model=model)

    # CellTypist expects log-normalized data (which we have)
    predictions = celltypist.annotate(
        adata,
        model=ct_model,
        majority_voting=majority_voting,
        over_clustering=over_clustering,
    )

    # Transfer predictions to adata
    adata.obs["celltypist_prediction"] = predictions.predicted_labels[
        "predicted_labels"
    ]
    if majority_voting:
        adata.obs["celltypist_majority_voting"] = predictions.predicted_labels[
            "majority_voting"
        ]
        label_col = "celltypist_majority_voting"
    else:
        label_col = "celltypist_prediction"

    # Summary
    type_counts = adata.obs[label_col].value_counts()
    n_types = type_counts.nunique()

    result = {
        "metrics": {
            "n_cell_types": n_types,
            "label_column": label_col,
            "cell_type_counts": {str(k): int(v) for k, v in type_counts.items()},
            "model_used": model,
        },
        "plots": [],
        "provenance": {
            "tool_id": "celltypist_annotation",
            "parameters": {
                "model": model,
                "majority_voting": majority_voting,
                "over_clustering": over_clustering,
            },
            "n_cell_types_found": n_types,
        },
        "warnings": [],
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_annotation")

    return result


def annotate_manual(
    adata: ad.AnnData,
    *,
    marker_dict: dict[str, list[str]],
    groupby: str = "leiden",
    layer: str | None = None,
    checkpoint_dir: str | None = None,
) -> dict:
    """Score clusters against user-provided marker genes.

    For each cluster, computes a score for each cell type based on
    mean expression of its marker genes. Returns ranked annotation
    suggestions per cluster for the user to confirm or override.

    Parameters
    ----------
    marker_dict
        Maps cell type names to lists of marker genes, e.g.::

            {
                "CD4 T": ["CD3D", "CD4", "IL7R"],
                "B cell": ["CD79A", "MS4A1", "CD19"],
                "NK": ["NKG7", "GNLY", "KLRB1"],
            }
    groupby
        obs column with cluster labels.
    layer
        Expression layer to use for scoring. *None* uses ``adata.X``
        (log-normalized). For raw counts use ``'counts'``.
    checkpoint_dir
        Save checkpoint after annotation.

    Returns
    -------
    Result dict with suggested annotations per cluster.
    """
    if groupby not in adata.obs.columns:
        raise ValueError(f"Column '{groupby}' not in adata.obs.")

    # Validate markers against available genes
    # Use adata.raw if available (has all genes), else adata
    if adata.raw is not None:
        available_genes = set(adata.raw.var_names)
    else:
        available_genes = set(adata.var_names)

    valid_markers = {}
    missing_markers = {}
    for cell_type, genes in marker_dict.items():
        found = [g for g in genes if g in available_genes]
        not_found = [g for g in genes if g not in available_genes]
        valid_markers[cell_type] = found
        if not_found:
            missing_markers[cell_type] = not_found

    # Score each cell type using sc.tl.score_genes
    score_keys = {}
    for cell_type, genes in valid_markers.items():
        if len(genes) == 0:
            continue
        score_key = f"score_{cell_type.replace(' ', '_')}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=score_key)
        score_keys[cell_type] = score_key

    # For each cluster, find the best-matching cell type
    groups = adata.obs[groupby].cat.categories.tolist()
    suggestions = {}
    for group in groups:
        mask = adata.obs[groupby] == group
        scores = {}
        for cell_type, key in score_keys.items():
            scores[cell_type] = float(adata.obs.loc[mask, key].mean())
        # Sort by score descending
        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        suggestions[str(group)] = [
            {"cell_type": ct, "score": round(s, 4)} for ct, s in ranked
        ]

    warnings = []
    if missing_markers:
        missing_str = "; ".join(
            f"{ct}: {genes}" for ct, genes in missing_markers.items()
        )
        warnings.append(f"Marker genes not found in data: {missing_str}")

    # Check for clusters where top score is very low (ambiguous)
    ambiguous = []
    for group, ranked in suggestions.items():
        if ranked and ranked[0]["score"] < 0.1:
            ambiguous.append(group)
    if ambiguous:
        warnings.append(
            f"Clusters with weak marker scores (top < 0.1): {ambiguous}. "
            "These may not match any provided cell type."
        )

    result = {
        "metrics": {
            "n_clusters": len(groups),
            "groupby": groupby,
            "suggestions": suggestions,
            "n_cell_types_provided": len(marker_dict),
            "missing_markers": missing_markers,
        },
        "plots": [],
        "provenance": {
            "tool_id": "manual_annotation",
            "parameters": {
                "groupby": groupby,
                "marker_dict": marker_dict,
                "layer": layer,
            },
        },
        "warnings": warnings,
    }

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_annotation")

    return result


def apply_annotation(
    adata: ad.AnnData,
    *,
    mapping: dict[str, str],
    groupby: str = "leiden",
    key_added: str = "cell_type",
) -> dict:
    """Apply a cluster-to-cell-type mapping.

    Parameters
    ----------
    mapping
        Maps cluster labels to cell type names, e.g.::

            {"0": "CD4 T", "1": "B cell", "2": "CD14 Mono", ...}

        Clusters not in the mapping are labeled ``'Unknown'``.
    groupby
        obs column with cluster labels.
    key_added
        obs column for the cell type annotations.

    Returns
    -------
    Result dict with cell type counts.
    """
    adata.obs[key_added] = (
        adata.obs[groupby]
        .map(mapping)
        .fillna("Unknown")
        .astype("category")
    )

    type_counts = adata.obs[key_added].value_counts()

    warnings = []
    n_unknown = int((adata.obs[key_added] == "Unknown").sum())
    if n_unknown > 0:
        warnings.append(
            f"{n_unknown} cells labeled 'Unknown' — some clusters are not "
            "in the mapping."
        )

    return {
        "metrics": {
            "n_cell_types": type_counts.nunique(),
            "cell_type_counts": {str(k): int(v) for k, v in type_counts.items()},
        },
        "plots": [],
        "provenance": {
            "tool_id": "apply_annotation",
            "parameters": {"mapping": mapping, "groupby": groupby, "key_added": key_added},
        },
        "warnings": warnings,
    }


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
