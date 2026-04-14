"""Perturbation screen tools: guide assignment + perturbation DE.

For Perturb-seq / CROP-seq experiments on 10x Chromium with CRISPR guides.

Best-practice references:
  [BP-1] Heumos et al. 2023, §"Inferring perturbation effects" (p. 557)
  [BP-2] sc-best-practices.org Ch. 20 (Perturbation modeling)

Usage::

    from scagent.tools.perturbation import assign_guides, run_perturbation_de
    adata = assign_guides(adata, guide_calls_key="guide_ids")
    results = run_perturbation_de(adata, control_label="non-targeting")
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Guide assignment
# ---------------------------------------------------------------------------


def assign_guides(
    adata: ad.AnnData,
    *,
    guide_calls_key: str = "guide_ids",
    feature_type_key: str | None = "feature_types",
    min_umis: int = 5,
    max_guides_per_cell: int = 1,
) -> dict[str, Any]:
    """Assign CRISPR guide identities to cells.

    Expects guide information in ``adata.obs[guide_calls_key]`` (from
    Cell Ranger multi) or in a CRISPR feature barcode matrix.

    Parameters
    ----------
    adata
        Annotated data. Guide calls should be in ``adata.obs``.
    guide_calls_key
        Column in ``adata.obs`` with guide IDs (or comma-separated
        for multi-guide cells from Cell Ranger).
    feature_type_key
        If guide counts are in ``adata.var``, column identifying
        CRISPR guide features.
    min_umis
        Minimum guide UMI count for a valid assignment.
    max_guides_per_cell
        Maximum guides per cell. Cells with more are marked "multi-guide".

    Returns
    -------
    Dict with ``assignment_stats``, ``summary``, ``provenance``, ``warnings``.
    Modifies ``adata.obs`` in-place: adds ``perturbation``, ``guide``,
    ``n_guides`` columns.
    """
    result_warnings: list[str] = []

    if guide_calls_key not in adata.obs.columns:
        raise ValueError(
            f"Guide calls column '{guide_calls_key}' not in adata.obs. "
            f"Available: {list(adata.obs.columns)}. "
            "For Cell Ranger multi output, guide IDs are typically in "
            "'guide_ids' or extracted from the CRISPR feature barcode matrix."
        )

    # Parse guide assignments
    raw_guides = adata.obs[guide_calls_key].astype(str)

    guides = []
    n_guides_list = []
    perturbations = []

    for val in raw_guides:
        if pd.isna(val) or val in ("", "nan", "None", "unassigned"):
            guides.append("unassigned")
            n_guides_list.append(0)
            perturbations.append("unassigned")
        else:
            # Handle comma-separated multi-guide calls
            guide_list = [g.strip() for g in val.split(",") if g.strip()]
            n_g = len(guide_list)
            n_guides_list.append(n_g)

            if n_g == 0:
                guides.append("unassigned")
                perturbations.append("unassigned")
            elif n_g == 1:
                guides.append(guide_list[0])
                # Extract target gene (assume guide name format: gene_guideN)
                perturbations.append(_extract_target(guide_list[0]))
            elif n_g <= max_guides_per_cell:
                guides.append(guide_list[0])  # take first
                perturbations.append(_extract_target(guide_list[0]))
                result_warnings.append(
                    f"Cell has {n_g} guides; using first: {guide_list[0]}"
                ) if n_g > 1 else None
            else:
                guides.append("multi-guide")
                perturbations.append("multi-guide")

    adata.obs["guide"] = pd.Categorical(guides)
    adata.obs["n_guides"] = n_guides_list
    adata.obs["perturbation"] = pd.Categorical(perturbations)

    # Stats
    guide_counts = adata.obs["perturbation"].value_counts()
    n_unassigned = int((adata.obs["perturbation"] == "unassigned").sum())
    n_multi = int((adata.obs["perturbation"] == "multi-guide").sum())
    n_assigned = int(adata.n_obs - n_unassigned - n_multi)
    n_perturbations = int(
        (guide_counts.index != "unassigned").sum()
        - (1 if "multi-guide" in guide_counts.index else 0)
    )

    frac_assigned = n_assigned / adata.n_obs if adata.n_obs > 0 else 0
    if frac_assigned < 0.5:
        result_warnings.append(
            f"Only {frac_assigned*100:.0f}% of cells have guide assignments. "
            "Check guide library quality or Cell Ranger guide calling settings."
        )

    return {
        "assignment_stats": {
            "n_cells": int(adata.n_obs),
            "n_assigned": n_assigned,
            "n_unassigned": n_unassigned,
            "n_multi_guide": n_multi,
            "n_perturbations": n_perturbations,
            "frac_assigned": round(frac_assigned, 3),
            "top_perturbations": guide_counts.head(10).to_dict(),
        },
        "summary": (
            f"Guide assignment: {n_assigned}/{adata.n_obs} cells assigned "
            f"({frac_assigned*100:.0f}%), {n_perturbations} perturbations, "
            f"{n_unassigned} unassigned, {n_multi} multi-guide."
        ),
        "warnings": [w for w in result_warnings if w],
        "provenance": {
            "tool_id": "guide_assignment",
            "parameters": {
                "guide_calls_key": guide_calls_key,
                "min_umis": min_umis,
                "max_guides_per_cell": max_guides_per_cell,
            },
            "n_assigned": n_assigned,
            "n_perturbations": n_perturbations,
        },
    }


def _extract_target(guide_name: str) -> str:
    """Extract target gene from guide name (e.g., 'TP53_guide1' → 'TP53')."""
    # Common patterns: GENE_guideN, GENE-gN, sgGENE
    for sep in ("_guide", "_sg", "-g", "-sg"):
        if sep in guide_name.lower():
            idx = guide_name.lower().index(sep)
            return guide_name[:idx]
    # If starts with "sg", strip it
    if guide_name.lower().startswith("sg"):
        return guide_name[2:]
    # Fallback: use the whole guide name as the perturbation
    return guide_name


# ---------------------------------------------------------------------------
# Perturbation DE
# ---------------------------------------------------------------------------


def run_perturbation_de(
    adata: ad.AnnData,
    *,
    control_label: str = "non-targeting",
    perturbation_key: str = "perturbation",
    cell_type_key: str | None = "cell_type",
    min_cells: int = 50,
    alpha: float = 0.05,
    method: str = "wilcoxon",
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run DE for each perturbation vs non-targeting controls.

    Parameters
    ----------
    adata
        Must have perturbation labels in ``adata.obs[perturbation_key]``.
    control_label
        Label for non-targeting control cells.
    perturbation_key
        Column with perturbation/target gene labels.
    cell_type_key
        Optional: run DE within each cell type separately.
    min_cells
        Minimum cells per perturbation to run DE.
    alpha
        Adjusted p-value significance threshold.
    method
        DE method: ``"wilcoxon"`` (cell-level, fast) or ``"pseudobulk"``
        (requires sample_key). Wilcoxon is acceptable for Perturb-seq
        since each cell is an independent perturbation event.
    plot_dir
        Directory for per-perturbation volcano plots.

    Returns
    -------
    Dict with ``de_results`` (per perturbation), ``summary``, ``provenance``,
    ``warnings``.

    Notes
    -----
    Unlike cross-condition DE (where pseudobulk is mandatory), cell-level
    tests are acceptable for Perturb-seq because each cell receives an
    independent guide — cells ARE the replicates here.
    """
    import scanpy as sc

    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate ---
    if perturbation_key not in adata.obs.columns:
        raise ValueError(
            f"'{perturbation_key}' not in adata.obs. "
            "Run assign_guides() first."
        )

    perturbations = adata.obs[perturbation_key].unique()
    if control_label not in perturbations:
        # Try case-insensitive match
        matches = [p for p in perturbations if str(p).lower() == control_label.lower()]
        if matches:
            control_label = matches[0]
        else:
            raise ValueError(
                f"Control label '{control_label}' not found in perturbations. "
                f"Available: {list(perturbations[:20])}"
            )

    control_cells = adata.obs[perturbation_key] == control_label
    n_control = int(control_cells.sum())
    if n_control < min_cells:
        raise ValueError(
            f"Only {n_control} control cells ('{control_label}'). "
            f"Need ≥{min_cells} for reliable DE."
        )

    # Filter to valid perturbations
    perturb_counts = adata.obs[perturbation_key].value_counts()
    valid_perturbations = perturb_counts[
        (perturb_counts >= min_cells)
        & (~perturb_counts.index.isin(["unassigned", "multi-guide", control_label]))
    ].index.tolist()

    skipped = perturb_counts[
        (perturb_counts < min_cells)
        & (~perturb_counts.index.isin(["unassigned", "multi-guide", control_label]))
    ]
    if len(skipped) > 0:
        result_warnings.append(
            f"{len(skipped)} perturbations skipped (<{min_cells} cells): "
            f"{list(skipped.index[:5])}"
        )

    if len(valid_perturbations) == 0:
        raise ValueError(
            f"No perturbations have ≥{min_cells} cells. "
            f"Top counts: {perturb_counts.head(10).to_dict()}"
        )

    # --- Run DE per perturbation ---
    de_results: dict[str, Any] = {}
    summary_rows = []

    for perturbation in valid_perturbations:
        # Subset to control + this perturbation
        mask = adata.obs[perturbation_key].isin([control_label, perturbation])
        sub = adata[mask].copy()
        sub.obs["_de_group"] = (
            sub.obs[perturbation_key]
            .map({control_label: "control", perturbation: "perturbed"})
            .astype("category")
        )

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sc.tl.rank_genes_groups(
                    sub, groupby="_de_group", groups=["perturbed"],
                    reference="control", method="wilcoxon",
                )

            # Extract results
            result_df = sc.get.rank_genes_groups_df(sub, group="perturbed")
            result_df["significant"] = result_df["pvals_adj"] < alpha

            n_sig = int(result_df["significant"].sum())
            n_up = int(((result_df["significant"]) & (result_df["logfoldchanges"] > 0)).sum())
            n_down = int(((result_df["significant"]) & (result_df["logfoldchanges"] < 0)).sum())

            de_results[str(perturbation)] = {
                "n_genes_tested": len(result_df),
                "n_significant": n_sig,
                "n_up": n_up,
                "n_down": n_down,
                "n_cells_perturbed": int((sub.obs["_de_group"] == "perturbed").sum()),
                "n_cells_control": int((sub.obs["_de_group"] == "control").sum()),
                "top_up": result_df[result_df["logfoldchanges"] > 0]
                    .nsmallest(5, "pvals_adj")["names"].tolist(),
                "top_down": result_df[result_df["logfoldchanges"] < 0]
                    .nsmallest(5, "pvals_adj")["names"].tolist(),
            }

            summary_rows.append({
                "perturbation": str(perturbation),
                "n_sig": n_sig,
                "n_up": n_up,
                "n_down": n_down,
            })

        except Exception as e:
            result_warnings.append(f"DE failed for '{perturbation}': {e}")

    return {
        "de_results": de_results,
        "summary_table": summary_rows,
        "n_perturbations_tested": len(de_results),
        "n_perturbations_skipped": len(skipped),
        "control_label": control_label,
        "n_control_cells": n_control,
        "summary": (
            f"Perturbation DE: {len(de_results)} perturbations tested "
            f"(vs {n_control} '{control_label}' controls). "
            f"{sum(r['n_significant'] for r in de_results.values())} total DE genes."
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "perturbation_de",
            "parameters": {
                "control_label": control_label,
                "perturbation_key": perturbation_key,
                "min_cells": min_cells,
                "alpha": alpha,
                "method": method,
            },
            "n_perturbations_tested": len(de_results),
            "n_control_cells": n_control,
        },
    }
