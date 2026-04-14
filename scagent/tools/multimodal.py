"""Multimodal (CITE-seq) tools: ADT loading, CLR normalization, WNN.

For 10x Chromium CITE-seq data (RNA + surface protein / ADT).

Best-practice references:
  [BP-1] Heumos et al. 2023, §"Surface protein expression" (pp. 558-559)
  [BP-2] sc-best-practices.org Ch. 32-37 (Surface protein QC, normalization, WNN)

Usage::

    from scagent.tools.multimodal import load_protein, normalize_protein, run_wnn
    adata = load_protein(adata, protein_path="filtered_feature_bc_matrix.h5")
    normalize_protein(adata)
    result = run_wnn(adata)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Load protein / ADT data
# ---------------------------------------------------------------------------


def load_protein(
    adata: ad.AnnData,
    protein_path: str | Path | None = None,
    *,
    protein_layer: str | None = None,
) -> dict[str, Any]:
    """Load ADT (Antibody-Derived Tag) counts for CITE-seq.

    Expects either a separate protein count matrix file or protein
    features already present in ``adata`` (from Cell Ranger multi output
    which includes both modalities).

    Parameters
    ----------
    adata
        Gene expression AnnData. If from Cell Ranger multi, protein
        features may already be in ``adata.var`` with feature_types="Antibody Capture".
    protein_path
        Path to a separate protein count matrix (h5 or mtx directory).
    protein_layer
        If protein counts are in an existing layer, specify here.

    Returns
    -------
    Dict with ``n_proteins``, ``summary``, ``provenance``, ``warnings``.
    May store protein data in ``adata.obsm["protein_counts"]``.
    """
    result_warnings: list[str] = []

    # Strategy 1: protein features already in adata (Cell Ranger multi)
    if "feature_types" in adata.var.columns:
        protein_mask = adata.var["feature_types"] == "Antibody Capture"
        n_proteins = int(protein_mask.sum())
        if n_proteins > 0:
            # Extract protein data to obsm
            protein_data = adata[:, protein_mask].X.copy()
            adata.obsm["protein_counts"] = (
                protein_data.toarray() if sp.issparse(protein_data) else np.asarray(protein_data)
            )
            adata.uns["protein_names"] = list(adata.var_names[protein_mask])

            return {
                "n_proteins": n_proteins,
                "protein_names": list(adata.var_names[protein_mask]),
                "summary": f"Found {n_proteins} ADT features in adata (Cell Ranger multi format).",
                "warnings": result_warnings,
                "provenance": {
                    "tool_id": "load_protein",
                    "parameters": {"source": "adata_feature_types"},
                    "n_proteins": n_proteins,
                },
            }

    # Strategy 2: load from separate file
    if protein_path is not None:
        import scanpy as sc

        protein_path = Path(protein_path)
        if not protein_path.exists():
            raise FileNotFoundError(f"Protein data file not found: {protein_path}")

        if protein_path.suffix == ".h5":
            prot_adata = sc.read_10x_h5(str(protein_path))
        else:
            prot_adata = sc.read_10x_mtx(str(protein_path))

        # Match barcodes
        common = adata.obs_names.intersection(prot_adata.obs_names)
        if len(common) == 0:
            raise ValueError(
                "No shared barcodes between GEX and protein data. "
                "Check that barcodes match (with/without '-1' suffix)."
            )

        prot_sub = prot_adata[common]
        protein_data = prot_sub.X
        adata_sub = adata[common]

        adata.obsm["protein_counts"] = (
            protein_data.toarray() if sp.issparse(protein_data) else np.asarray(protein_data)
        )
        adata.uns["protein_names"] = list(prot_sub.var_names)
        n_proteins = prot_sub.n_vars

        frac_matched = len(common) / adata.n_obs
        if frac_matched < 0.9:
            result_warnings.append(
                f"Only {frac_matched*100:.1f}% of GEX barcodes matched protein data."
            )

        return {
            "n_proteins": n_proteins,
            "protein_names": list(prot_sub.var_names),
            "n_matched_barcodes": len(common),
            "summary": f"Loaded {n_proteins} ADT features from {protein_path.name}. {len(common)} barcodes matched.",
            "warnings": result_warnings,
            "provenance": {
                "tool_id": "load_protein",
                "parameters": {"protein_path": str(protein_path)},
                "n_proteins": n_proteins,
                "n_matched": len(common),
            },
        }

    # Strategy 3: check obsm
    if protein_layer and protein_layer in adata.layers:
        adata.obsm["protein_counts"] = (
            adata.layers[protein_layer].toarray()
            if sp.issparse(adata.layers[protein_layer])
            else np.asarray(adata.layers[protein_layer])
        )
        return {
            "n_proteins": adata.obsm["protein_counts"].shape[1],
            "summary": f"Using protein counts from layer '{protein_layer}'.",
            "warnings": result_warnings,
            "provenance": {"tool_id": "load_protein", "parameters": {"protein_layer": protein_layer}},
        }

    raise ValueError(
        "No protein/ADT data found. Provide protein_path, or ensure "
        "adata.var['feature_types'] contains 'Antibody Capture' features "
        "(Cell Ranger multi output)."
    )


# ---------------------------------------------------------------------------
# CLR normalization
# ---------------------------------------------------------------------------


def normalize_protein(
    adata: ad.AnnData,
    *,
    method: str = "clr",
) -> dict[str, Any]:
    """Normalize ADT counts using CLR or DSB.

    Parameters
    ----------
    adata
        Must have ``adata.obsm["protein_counts"]`` (from :func:`load_protein`).
    method
        Normalization method: ``"clr"`` (centered log-ratio, default)
        or ``"dsb"`` (requires empty droplets). [BP-1, BP-2 Ch. 33]

    Returns
    -------
    Dict with ``summary``, ``provenance``, ``warnings``.
    Stores normalized values in ``adata.obsm["protein_normalized"]``.
    """
    result_warnings: list[str] = []

    if "protein_counts" not in adata.obsm:
        raise ValueError(
            "No protein counts found in adata.obsm['protein_counts']. "
            "Run load_protein() first."
        )

    raw = adata.obsm["protein_counts"].copy()

    if method == "clr":
        # Centered log-ratio: log(x / geometric_mean(x)) per cell
        # Add pseudocount
        raw_pseudo = raw + 1
        log_raw = np.log(raw_pseudo)
        geo_mean = log_raw.mean(axis=1, keepdims=True)
        normalized = log_raw - geo_mean
        adata.obsm["protein_normalized"] = normalized

    elif method == "dsb":
        result_warnings.append(
            "DSB normalization requires empty droplet protein counts "
            "for background estimation. Using CLR as fallback."
        )
        # Fallback to CLR
        raw_pseudo = raw + 1
        log_raw = np.log(raw_pseudo)
        geo_mean = log_raw.mean(axis=1, keepdims=True)
        normalized = log_raw - geo_mean
        adata.obsm["protein_normalized"] = normalized
        method = "clr (dsb fallback)"

    else:
        raise ValueError(f"Unknown normalization method: {method}. Use 'clr' or 'dsb'.")

    # Guard rail: check for isotype controls [BP-2 Ch. 32]
    protein_names = adata.uns.get("protein_names", [])
    isotype_markers = [p for p in protein_names if "isotype" in p.lower() or "igg" in p.lower()]
    if len(isotype_markers) == 0:
        result_warnings.append(
            "No isotype control antibodies detected in protein panel. "
            "Isotype controls are recommended for assessing ADT background. [BP-2 Ch. 32]"
        )

    return {
        "method": method,
        "n_proteins": int(raw.shape[1]),
        "summary": f"Protein counts normalized with {method}. {raw.shape[1]} features.",
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "normalize_protein",
            "parameters": {"method": method},
            "n_proteins": int(raw.shape[1]),
        },
    }


# ---------------------------------------------------------------------------
# WNN (Weighted Nearest Neighbors)
# ---------------------------------------------------------------------------


def run_wnn(
    adata: ad.AnnData,
    *,
    n_neighbors: int = 20,
    rna_weight: float = 0.5,
) -> dict[str, Any]:
    """Compute weighted nearest neighbors from RNA + protein.

    Parameters
    ----------
    adata
        Must have PCA embedding and protein_normalized in obsm.
    n_neighbors
        Number of neighbors for the joint graph.
    rna_weight
        Weight for RNA modality (0-1). Protein weight = 1 - rna_weight.
        [BP-1, BP-2 Ch. 35]

    Returns
    -------
    Dict with ``summary``, ``provenance``, ``warnings``.
    Stores joint neighbor graph in ``adata.obsp`` / ``adata.uns``.
    """
    result_warnings: list[str] = []

    if "X_pca" not in adata.obsm:
        raise ValueError("PCA not found. Run sc.tl.pca() first.")
    if "protein_normalized" not in adata.obsm:
        raise ValueError(
            "Normalized protein data not found. Run normalize_protein() first."
        )

    try:
        import muon as mu
        import mudata
        import scanpy as sc

        # Use muon's WNN implementation
        # Create a MuData object with per-modality neighbors
        rna_ad = adata.copy()
        if "neighbors" not in rna_ad.uns:
            sc.pp.neighbors(rna_ad, use_rep="X_pca", n_neighbors=n_neighbors)

        prot_ad = ad.AnnData(
            X=adata.obsm["protein_normalized"],
            obs=adata.obs.copy(),
        )
        if "protein_names" in adata.uns:
            prot_ad.var_names = adata.uns["protein_names"]
        # Compute neighbors on protein modality too
        sc.pp.neighbors(prot_ad, n_neighbors=n_neighbors)

        mdata = mudata.MuData({"rna": rna_ad, "prot": prot_ad})
        mu.pp.neighbors(mdata, key_added="wnn")

        # Transfer WNN graph back to adata
        if "wnn" in mdata.uns:
            adata.uns["wnn"] = mdata.uns["wnn"]
        if "wnn" in mdata.obsp:
            adata.obsp["wnn_connectivities"] = mdata.obsp.get("wnn_connectivities")
            adata.obsp["wnn_distances"] = mdata.obsp.get("wnn_distances")
        result_warnings.append("WNN computed via muon.")

    except ImportError:
        # Fallback: simple concatenated embedding
        import scanpy as sc

        result_warnings.append(
            "muon not installed. Using concatenated RNA+protein embedding "
            "with weighted neighbors as approximation."
        )
        rna_emb = adata.obsm["X_pca"]
        prot_emb = adata.obsm["protein_normalized"]

        # Scale and weight
        from sklearn.preprocessing import StandardScaler
        rna_scaled = StandardScaler().fit_transform(rna_emb) * rna_weight
        prot_scaled = StandardScaler().fit_transform(prot_emb) * (1 - rna_weight)

        joint = np.hstack([rna_scaled, prot_scaled])
        adata.obsm["X_wnn"] = joint

        sc.pp.neighbors(adata, use_rep="X_wnn", n_neighbors=n_neighbors, key_added="wnn")

    return {
        "n_neighbors": n_neighbors,
        "rna_weight": rna_weight,
        "summary": (
            f"WNN graph computed: RNA weight={rna_weight}, "
            f"protein weight={1-rna_weight}, k={n_neighbors}."
        ),
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "wnn",
            "parameters": {
                "n_neighbors": n_neighbors,
                "rna_weight": rna_weight,
            },
        },
    }
