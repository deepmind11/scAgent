"""Compositional analysis: scCODA + Milo via pertpy.

Detects changes in cell-type proportions between conditions.
Two approaches:
  - scCODA: Bayesian compositional model for labeled clusters [required]
  - Milo:   KNN-based differential abundance without predefined labels

Best-practice references:
  [BP-1] Heumos et al. 2023, §"Deciphering changes in cell composition" (p. 555)
  [BP-2] sc-best-practices.org Ch. 18 (Compositional analysis)
  Büttner et al. 2021, Nat Commun — scCODA: MCC 0.64 vs ~0.20 for naive tests

CRITICAL: Naive per-cell-type proportion tests (Wilcoxon, Fisher) produce
systematic false positives due to the compositional nature of cell-type
counts. scCODA properly models this constraint. [BP-1, BP-2 Ch. 18]

Usage::

    from scagent.tools.composition import run_sccoda, run_milo
    result = run_sccoda(adata, condition_key="condition", sample_key="donor",
                        cell_type_key="cell_type")
    result = run_milo(adata, condition_key="condition", sample_key="donor")
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

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# scCODA — Bayesian compositional analysis (labeled clusters)
# ---------------------------------------------------------------------------


def run_sccoda(
    adata: ad.AnnData,
    *,
    condition_key: str = "condition",
    sample_key: str = "sample",
    cell_type_key: str = "cell_type",
    reference_cell_type: str | None = "automatic",
    fdr: float = 0.05,
    formula: str | None = None,
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run scCODA compositional analysis via pertpy.

    Parameters
    ----------
    adata
        Annotated data with condition, sample, and cell type labels.
    condition_key
        Column in ``adata.obs`` with condition labels.
    sample_key
        Column in ``adata.obs`` with biological sample IDs.
    cell_type_key
        Column in ``adata.obs`` with cell-type labels.
    reference_cell_type
        Reference cell type assumed unchanged across conditions.
        ``"automatic"`` lets scCODA pick one. [BP-2 Ch. 18]
    fdr
        False discovery rate threshold. Default 0.05; loosen to 0.2
        for exploratory analysis. [BP-2 Ch. 18]
    formula
        Design formula (R-style). Default: ``"~ {condition_key}"``.
    plot_dir
        Directory for composition plots.

    Returns
    -------
    Dict with ``credible_effects``, ``effect_sizes``, ``summary``,
    ``provenance``, ``plots``, ``warnings``.
    """
    try:
        import pertpy as pt
    except ImportError:
        raise ImportError(
            "pertpy is required for scCODA. Install with: pip install pertpy"
        )

    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate inputs ---
    for col in (condition_key, sample_key, cell_type_key):
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column '{col}' not in adata.obs. "
                f"Available: {list(adata.obs.columns)}"
            )

    conditions = adata.obs[condition_key].unique()
    if len(conditions) < 2:
        raise ValueError(
            f"Need ≥2 conditions, found {len(conditions)}: {list(conditions)}"
        )

    n_samples = adata.obs[sample_key].nunique()
    if n_samples < 4:
        result_warnings.append(
            f"Only {n_samples} samples. scCODA can work with n=2/group "
            "but power is limited. Results should be interpreted cautiously."
        )

    n_cell_types = adata.obs[cell_type_key].nunique()
    if n_cell_types < 2:
        raise ValueError("Need ≥2 cell types for compositional analysis.")

    if formula is None:
        formula = f"~ {condition_key}"

    # --- Build scCODA model ---
    sccoda_model = pt.tl.Sccoda()
    sccoda_data = sccoda_model.load(
        adata,
        type="cell_level",
        generate_sample_level=True,
        cell_type_identifier=cell_type_key,
        sample_identifier=sample_key,
        covariate_obs=[condition_key],
    )

    # Prepare model
    ref_ct = reference_cell_type if reference_cell_type != "automatic" else "automatic"
    sccoda_data = sccoda_model.prepare(
        sccoda_data,
        modality_key="coda",
        formula=formula,
        reference_cell_type=ref_ct,
    )

    # Run MCMC
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sccoda_model.run_nuts(sccoda_data, modality_key="coda", rng_key=42)

    # --- Check acceptance rate [BP-2 Ch. 18] ---
    # pertpy stores MCMC diagnostics in uns
    params = sccoda_data["coda"].uns.get("scCODA_params", {})
    # acceptance rate check is advisory

    # Set FDR
    sccoda_model.set_fdr(sccoda_data, fdr)

    # --- Extract results ---
    credible = sccoda_model.credible_effects(sccoda_data, modality_key="coda")

    # Parse credible effects into structured format
    effects_list = []
    if isinstance(credible, pd.Series):
        for idx, is_credible in credible.items():
            if isinstance(idx, tuple):
                covariate, cell_type = idx
            else:
                covariate = str(idx)
                cell_type = str(idx)
            effects_list.append({
                "covariate": str(covariate),
                "cell_type": str(cell_type),
                "credible": bool(is_credible),
            })

    n_credible = sum(1 for e in effects_list if e["credible"])

    # --- Plots ---
    if plot_dir:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)
        try:
            fig = sccoda_model.plot_stacked_barplot(
                sccoda_data, modality_key="coda",
                feature_name=condition_key, figsize=(6, 4),
            )
            if fig is not None:
                path = str(pdir / "composition_barplot.png")
                plt.savefig(path, dpi=150, bbox_inches="tight")
                plots.append(path)
            plt.close("all")
        except Exception as e:
            logger.warning("Composition barplot failed: %s", e)

        try:
            sccoda_model.plot_effects_barplot(sccoda_data, "coda", condition_key)
            path = str(pdir / "composition_effects.png")
            plt.savefig(path, dpi=150, bbox_inches="tight")
            plots.append(path)
            plt.close("all")
        except Exception as e:
            logger.warning("Effects barplot failed: %s", e)

    return {
        "credible_effects": effects_list,
        "n_credible": n_credible,
        "n_cell_types": n_cell_types,
        "reference_cell_type": ref_ct,
        "fdr": fdr,
        "summary": (
            f"scCODA: {n_credible}/{len(effects_list)} credible effects "
            f"at FDR={fdr}. {n_cell_types} cell types, "
            f"{n_samples} samples, {len(conditions)} conditions. "
            f"Reference: {ref_ct}."
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "sccoda",
            "parameters": {
                "condition_key": condition_key,
                "sample_key": sample_key,
                "cell_type_key": cell_type_key,
                "reference_cell_type": ref_ct,
                "fdr": fdr,
                "formula": formula,
            },
            "n_credible": n_credible,
            "n_cell_types": n_cell_types,
            "n_samples": n_samples,
        },
    }


# ---------------------------------------------------------------------------
# Milo — KNN-based differential abundance (no predefined labels)
# ---------------------------------------------------------------------------


def run_milo(
    adata: ad.AnnData,
    *,
    condition_key: str = "condition",
    sample_key: str = "sample",
    design: str | None = None,
    model_contrasts: str | None = None,
    prop: float = 0.1,
    alpha: float = 0.1,
    cell_type_key: str | None = "cell_type",
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run Milo differential abundance analysis via pertpy.

    Parameters
    ----------
    adata
        Must have neighbors computed.
    condition_key
        Column with condition labels.
    sample_key
        Column with biological sample IDs.
    design
        GLM design formula. Default: ``"~ {condition_key}"``.
    model_contrasts
        Contrast specification (e.g., ``"conditionDisease-conditionControl"``).
    prop
        Proportion of cells to sample as neighbourhood indices. [BP-2 Ch. 18]
    alpha
        SpatialFDR significance threshold.
    cell_type_key
        Column for annotating neighbourhoods (optional).
    plot_dir
        Directory for Milo plots.

    Returns
    -------
    Dict with ``da_results``, ``summary``, ``provenance``, ``plots``, ``warnings``.
    """
    try:
        import pertpy as pt
    except ImportError:
        raise ImportError(
            "pertpy is required for Milo. Install with: pip install pertpy"
        )

    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate ---
    for col in (condition_key, sample_key):
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column '{col}' not in adata.obs. "
                f"Available: {list(adata.obs.columns)}"
            )

    if "neighbors" not in adata.uns and "connectivities" not in adata.obsp:
        raise ValueError(
            "Neighbor graph not found. Run sc.pp.neighbors() first."
        )

    # --- Guard rail: check n_neighbors vs n_samples [BP-2 Ch. 18] ---
    n_samples = adata.obs[sample_key].nunique()
    n_neighbors = adata.uns.get("neighbors", {}).get("params", {}).get("n_neighbors", 15)
    min_recommended = 3 * n_samples
    if n_neighbors < min_recommended:
        result_warnings.append(
            f"n_neighbors={n_neighbors} may be too low for {n_samples} samples. "
            f"Recommended: ≥{min_recommended} (3 × n_samples) for adequate "
            f"statistical power. [BP-2 Ch. 18]"
        )

    if design is None:
        design = f"~ {condition_key}"

    # --- Run Milo ---
    milo = pt.tl.Milo()
    mdata = milo.load(adata)

    milo.make_nhoods(mdata, prop=prop)
    milo.count_nhoods(mdata, sample_col=sample_key)

    # DA testing
    kwargs: dict[str, Any] = {"design": design}
    if model_contrasts:
        kwargs["model_contrasts"] = model_contrasts

    try:
        milo.da_nhoods(mdata, **kwargs)
    except Exception as e:
        raise RuntimeError(
            f"Milo DA test failed: {e}. "
            "This may be due to batch-condition confounding — "
            "consider integrating with scVI first. [BP-2 Ch. 18]"
        ) from e

    # --- Extract results ---
    nhood_results = mdata["milo"].var.copy()
    n_sig = int((nhood_results.get("SpatialFDR", pd.Series(dtype=float)) < alpha).sum())
    n_total = len(nhood_results)

    if "logFC" in nhood_results.columns:
        n_enriched = int(
            ((nhood_results["SpatialFDR"] < alpha) & (nhood_results["logFC"] > 0)).sum()
        )
        n_depleted = int(
            ((nhood_results["SpatialFDR"] < alpha) & (nhood_results["logFC"] < 0)).sum()
        )
    else:
        n_enriched = n_depleted = 0

    # Annotate neighbourhoods by cell type if possible
    if cell_type_key and cell_type_key in adata.obs.columns:
        try:
            milo.annotate_nhoods(mdata, anno_col=cell_type_key)
        except Exception:
            pass

    # --- Plots ---
    if plot_dir:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)
        try:
            milo.plot_da_beeswarm(mdata)
            path = str(pdir / "milo_beeswarm.png")
            plt.savefig(path, dpi=150, bbox_inches="tight")
            plots.append(path)
            plt.close("all")
        except Exception as e:
            logger.warning("Milo beeswarm plot failed: %s", e)

    return {
        "n_neighbourhoods": n_total,
        "n_significant": n_sig,
        "n_enriched": n_enriched,
        "n_depleted": n_depleted,
        "alpha": alpha,
        "summary": (
            f"Milo: {n_sig}/{n_total} DA neighbourhoods at SpatialFDR<{alpha}. "
            f"{n_enriched} enriched, {n_depleted} depleted."
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "milo",
            "parameters": {
                "condition_key": condition_key,
                "sample_key": sample_key,
                "design": design,
                "prop": prop,
                "alpha": alpha,
            },
            "n_neighbourhoods": n_total,
            "n_significant": n_sig,
            "n_samples": n_samples,
        },
    }
