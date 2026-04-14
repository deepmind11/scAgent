"""Immune repertoire tools: VDJ loading, clonotype analysis, diversity, overlap.

Uses Scirpy (scverse ecosystem) for TCR/BCR analysis from 10x Chromium VDJ.

Best-practice references:
  [BP-1] Heumos et al. 2023, §"Adaptive immune receptor repertoires" (pp. 559-560)
  [BP-2] sc-best-practices.org Ch. 38-39 (IR profiling, clonotype analysis)

Usage::

    from scagent.tools.repertoire import load_vdj, run_clonotype_analysis
    adata = load_vdj(adata, vdj_path="filtered_contig_annotations.csv")
    result = run_clonotype_analysis(adata)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _require_scirpy():
    try:
        import scirpy as ir
        return ir
    except ImportError:
        raise ImportError(
            "scirpy is required for immune repertoire analysis. "
            "Install with: pip install scirpy"
        )


# ---------------------------------------------------------------------------
# Load V(D)J data
# ---------------------------------------------------------------------------


def load_vdj(
    adata: ad.AnnData,
    vdj_path: str | Path,
    *,
    receptor_type: str = "auto",
) -> dict[str, Any]:
    """Load 10x VDJ contig data and merge with GEX AnnData.

    Parameters
    ----------
    adata
        Gene expression AnnData (must share barcodes with VDJ data).
    vdj_path
        Path to ``filtered_contig_annotations.csv`` from Cell Ranger VDJ.
    receptor_type
        ``"TCR"``, ``"BCR"``, or ``"auto"`` (detect from data).

    Returns
    -------
    Dict with ``n_cells_with_ir``, ``summary``, ``provenance``, ``warnings``.
    Modifies ``adata`` in-place: adds IR fields to ``.obs`` and ``.obsm``.
    """
    ir = _require_scirpy()

    result_warnings: list[str] = []
    vdj_path = Path(vdj_path)

    if not vdj_path.exists():
        raise FileNotFoundError(f"VDJ file not found: {vdj_path}")

    # Load VDJ data
    adata_vdj = ir.io.read_10x_vdj(str(vdj_path))

    # Merge with GEX
    ir.pp.merge_with_ir(adata, adata_vdj)

    # Count cells with IR data
    has_ir = adata.obs.get("has_ir", pd.Series(False, index=adata.obs.index))
    if "ir_status" in adata.obs.columns:
        has_ir = adata.obs["ir_status"] != "orphan"
    n_with_ir = int(has_ir.sum())
    frac_ir = n_with_ir / adata.n_obs if adata.n_obs > 0 else 0

    if frac_ir < 0.1:
        result_warnings.append(
            f"Only {frac_ir*100:.1f}% of cells have IR data. "
            "Check barcode matching between GEX and VDJ libraries."
        )

    return {
        "n_cells_with_ir": n_with_ir,
        "frac_with_ir": round(frac_ir, 3),
        "n_cells_total": int(adata.n_obs),
        "summary": (
            f"VDJ data loaded: {n_with_ir}/{adata.n_obs} cells "
            f"({frac_ir*100:.1f}%) have immune receptor data."
        ),
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "load_vdj",
            "parameters": {"vdj_path": str(vdj_path), "receptor_type": receptor_type},
            "n_cells_with_ir": n_with_ir,
        },
    }


# ---------------------------------------------------------------------------
# Clonotype analysis
# ---------------------------------------------------------------------------


def run_clonotype_analysis(
    adata: ad.AnnData,
    *,
    sequence: str = "aa",
    metric: str = "identity",
    receptor_arms: str = "all",
    dual_ir: str = "primary_only",
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run clonotype definition, expansion, and diversity analysis.

    Parameters
    ----------
    adata
        Must have IR data merged (from :func:`load_vdj`).
    sequence
        Clonotype definition: ``"aa"`` (amino acid) or ``"nt"`` (nucleotide).
    metric
        Distance metric for clonotype clustering.
    receptor_arms
        Which receptor arms to consider.
    dual_ir
        How to handle dual IR (cells with two productive chains).
    plot_dir
        Directory for clonotype plots.

    Returns
    -------
    Dict with ``diversity``, ``expansion``, ``summary``, ``provenance``,
    ``plots``, ``warnings``.
    """
    ir = _require_scirpy()

    result_warnings: list[str] = []
    plots: list[str] = []

    # Define clonotypes
    ir.pp.ir_dist(adata, metric=metric, sequence=sequence)
    ir.tl.define_clonotypes(adata, receptor_arms=receptor_arms, dual_ir=dual_ir)

    # Clonal expansion
    ir.tl.clonal_expansion(adata)

    # Diversity metrics
    try:
        diversity_df = ir.tl.alpha_diversity(adata, groupby="clone_id" if "clone_id" in adata.obs.columns else None)
    except Exception:
        diversity_df = None

    # Count clonotypes
    if "clone_id" in adata.obs.columns:
        clone_counts = adata.obs["clone_id"].value_counts()
        n_clonotypes = int(clone_counts.nunique())
        n_expanded = int((clone_counts > 1).sum())
        top_clonotypes = clone_counts.head(10).to_dict()
    else:
        n_clonotypes = 0
        n_expanded = 0
        top_clonotypes = {}
        result_warnings.append("No clone_id found after clonotype definition.")

    # Expansion categories
    expansion_stats = {}
    if "clonal_expansion" in adata.obs.columns:
        expansion_stats = adata.obs["clonal_expansion"].value_counts().to_dict()

    # Plots
    if plot_dir:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(8, 5))
            ir.pl.clonal_expansion(adata, ax=ax, show=False)
            path = str(pdir / "clonal_expansion.png")
            fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            plots.append(path)
        except Exception as e:
            logger.warning("Clonal expansion plot failed: %s", e)

    return {
        "n_clonotypes": n_clonotypes,
        "n_expanded": n_expanded,
        "expansion_stats": expansion_stats,
        "top_clonotypes": top_clonotypes,
        "summary": (
            f"Clonotype analysis: {n_clonotypes} clonotypes defined, "
            f"{n_expanded} expanded (>1 cell). "
            f"Expansion: {expansion_stats}."
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "clonotype_analysis",
            "parameters": {
                "sequence": sequence,
                "metric": metric,
                "receptor_arms": receptor_arms,
            },
            "n_clonotypes": n_clonotypes,
            "n_expanded": n_expanded,
        },
    }
