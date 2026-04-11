"""Load 10x Genomics Chromium data."""

from __future__ import annotations

import hashlib
from pathlib import Path

import anndata as ad
import scanpy as sc


def load_10x_h5(
    filename: str,
    *,
    gex_only: bool = True,
    genome: str | None = None,
    checkpoint_dir: str | None = None,
) -> tuple[ad.AnnData, dict]:
    """Load a 10x Genomics filtered_feature_bc_matrix.h5 file.

    Parameters
    ----------
    filename
        Path to the .h5 file.
    gex_only
        If True, load only gene expression features (exclude antibody capture, etc.).
    genome
        Genome annotation. Required only for multi-genome references.
    checkpoint_dir
        If provided, save checkpoint to ``{checkpoint_dir}/after_load.h5ad``.

    Returns
    -------
    (adata, result) where result contains metrics, provenance, and warnings.
    """
    path = Path(filename)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filename}")

    # Compute file hash for provenance
    file_hash = _sha256_file(path)

    # Load
    adata = sc.read_10x_h5(str(path), gex_only=gex_only, genome=genome)

    # Fix duplicate gene names (common in 10x references)
    n_dupes = adata.var_names.duplicated().sum()
    if n_dupes > 0:
        adata.var_names_make_unique()

    # Detect species from gene name case
    species = _detect_species(adata)
    adata.uns["species"] = species

    # Build result
    result = {
        "metrics": {
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "species": species,
            "duplicate_var_names_fixed": int(n_dupes),
        },
        "plots": [],
        "provenance": {
            "tool_id": "load_10x_h5",
            "parameters": {
                "filename": str(path.name),
                "gex_only": gex_only,
                "genome": genome,
            },
            "file_hash": file_hash,
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
        },
        "warnings": [],
    }

    if n_dupes > 0:
        result["warnings"].append(
            f"Fixed {n_dupes} duplicate gene names with var_names_make_unique()."
        )

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_load")

    return adata, result


def _detect_species(adata: ad.AnnData) -> str:
    """Detect species from gene name case convention.

    Human genes: uppercase (e.g., CD3E, MT-CO1)
    Mouse genes: title-case (e.g., Cd3e, mt-Co1)
    """
    sample = adata.var_names[:200].tolist()
    n_upper = sum(1 for g in sample if g == g.upper())
    # If >70% are all-uppercase → human
    if n_upper / len(sample) > 0.7:
        return "human"
    return "mouse"


def _sha256_file(path: Path, chunk_size: int = 65536) -> str:
    """Compute SHA-256 hash of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    """Save an AnnData checkpoint and return the file path."""
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
