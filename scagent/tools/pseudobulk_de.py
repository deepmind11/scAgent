"""Pseudobulk differential expression via PyDESeq2.

The #1 statistical error in scRNA-seq cross-condition DE is treating cells
as independent replicates.  This module **aggregates cells into pseudobulk
samples** (sum of raw counts per cell type × biological sample) and then
runs DESeq2, which models biological replicate variance correctly.

Reference:
    Squair et al. (2021) *Nature Communications* 12:5692
    "Confronting false discoveries in single-cell differential expression"

Usage::

    from scagent.tools.pseudobulk_de import run_pseudobulk_de
    results = run_pseudobulk_de(
        adata,
        cell_type_key="cell_type",
        sample_key="donor",
        condition_key="condition",
    )
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
# Pseudobulk aggregation
# ---------------------------------------------------------------------------


def aggregate_pseudobulk(
    adata: ad.AnnData,
    cell_type_key: str,
    sample_key: str,
    condition_key: str,
    *,
    min_cells: int = 10,
    layer: str | None = None,
) -> dict[str, dict[str, Any]]:
    """Aggregate single-cell counts into pseudobulk samples.

    Parameters
    ----------
    adata
        Annotated data with raw counts in ``.X`` (or ``.raw.X``, or
        a specified *layer*).
    cell_type_key
        Column in ``adata.obs`` with cell-type labels.
    sample_key
        Column in ``adata.obs`` with biological-sample IDs.
    condition_key
        Column in ``adata.obs`` with condition labels.
    min_cells
        Minimum cells required to form a valid pseudobulk sample.
        Cell type × sample combinations below this are dropped.
    layer
        AnnData layer to use for raw counts.  *None* means ``.X``
        (falls back to ``.raw.X`` if ``.X`` looks log-normalized).

    Returns
    -------
    Dict keyed by cell type, each containing:
        - ``counts``:  pseudobulk count DataFrame (genes × samples)
        - ``metadata``: sample metadata DataFrame (samples × cols)
        - ``n_cells``:  dict of cell counts per sample
        - ``dropped``:  list of sample IDs dropped for too few cells
    """
    for col in (cell_type_key, sample_key, condition_key):
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column '{col}' not in adata.obs. Available: {list(adata.obs.columns)}"
            )

    # Decide which matrix has raw counts
    X = _get_raw_counts(adata, layer)
    var_names = _get_var_names(adata, layer)

    cell_types = adata.obs[cell_type_key].unique()
    result: dict[str, dict[str, Any]] = {}

    for ct in sorted(cell_types):
        ct_mask = adata.obs[cell_type_key] == ct
        ct_obs = adata.obs.loc[ct_mask]
        ct_X = X[ct_mask.values]

        samples = ct_obs[sample_key].unique()
        counts_list: list[pd.Series] = []
        meta_rows: list[dict] = []
        n_cells: dict[str, int] = {}
        dropped: list[str] = []

        for sample in sorted(samples):
            s_mask = ct_obs[sample_key] == sample
            n = int(s_mask.sum())

            if n < min_cells:
                dropped.append(str(sample))
                continue

            # Sum raw counts across cells → one pseudobulk vector
            s_X = ct_X[s_mask.values]
            if sp.issparse(s_X):
                summed = np.asarray(s_X.sum(axis=0)).flatten()
            else:
                summed = np.asarray(s_X.sum(axis=0)).flatten()

            col_name = f"{ct}__{sample}"
            counts_list.append(pd.Series(summed, index=var_names, name=col_name))
            n_cells[str(sample)] = n

            # One condition per sample — take the first
            cond = ct_obs.loc[s_mask, condition_key].iloc[0]
            meta_rows.append({
                "sample": str(sample),
                "condition": str(cond),
                "n_cells": n,
            })

        if not counts_list:
            logger.warning("Cell type '%s': no samples passed min_cells=%d filter", ct, min_cells)
            continue

        counts_df = pd.DataFrame(counts_list).T.astype(int)
        counts_df.index = var_names
        meta_df = pd.DataFrame(meta_rows).set_index("sample")

        result[str(ct)] = {
            "counts": counts_df,
            "metadata": meta_df,
            "n_cells": n_cells,
            "dropped": dropped,
        }

    return result


# ---------------------------------------------------------------------------
# DESeq2 wrapper
# ---------------------------------------------------------------------------


def run_pseudobulk_de(
    adata: ad.AnnData,
    *,
    cell_type_key: str = "cell_type",
    sample_key: str = "sample",
    condition_key: str = "condition",
    min_cells_per_pseudobulk: int = 10,
    design_formula: str = "~ condition",
    alpha: float = 0.05,
    min_replicates: int = 2,
    layer: str | None = None,
    plot_dir: str | None = None,
    lfc_threshold: float = 0.0,
) -> dict:
    """Run pseudobulk DE on all cell types.

    Parameters
    ----------
    adata
        AnnData with raw counts, cell type labels, sample IDs, and
        condition labels in ``.obs``.
    cell_type_key
        Column with cell-type annotations.
    sample_key
        Column with biological-sample IDs (the *real* replicates).
    condition_key
        Column with condition labels (e.g., ``"disease"``/``"healthy"``).
    min_cells_per_pseudobulk
        Minimum cells to form a pseudobulk sample.
    design_formula
        DESeq2 design formula (default ``"~ condition"``).
    alpha
        Adjusted p-value significance threshold.
    min_replicates
        Hard minimum biological replicates per condition per cell type.
        Cell types below this are skipped with a warning.
    layer
        AnnData layer with raw counts (*None* = ``.X``).
    plot_dir
        Directory for volcano plots (one per cell type).
    lfc_threshold
        Minimum |log2FC| for a gene to be reported as significant.

    Returns
    -------
    Result dict with ``de_results`` (per cell type), ``summary``,
    ``provenance``, ``warnings``, ``plots``.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # --- Validate inputs ---
    conditions = adata.obs[condition_key].unique()
    if len(conditions) < 2:
        raise ValueError(
            f"Need at least 2 conditions for DE, found {len(conditions)}: {list(conditions)}. "
            "For marker gene detection (one-vs-rest within clusters), use wilcoxon_markers instead."
        )

    n_samples = adata.obs[sample_key].nunique()
    if n_samples < 2 * min_replicates:
        raise ValueError(
            f"Need at least {min_replicates} biological replicates per condition "
            f"({2 * min_replicates} total samples), found {n_samples}. "
            "Pseudobulk DE is not appropriate with fewer replicates."
        )

    # --- Aggregate ---
    pb = aggregate_pseudobulk(
        adata, cell_type_key, sample_key, condition_key,
        min_cells=min_cells_per_pseudobulk, layer=layer,
    )

    de_results: dict[str, pd.DataFrame] = {}
    skipped: list[dict] = []
    all_warnings: list[str] = []
    plots: list[str] = []

    for ct, data in pb.items():
        meta = data["metadata"]
        counts = data["counts"]

        # Check replicates per condition
        cond_counts = meta["condition"].value_counts()
        low = cond_counts[cond_counts < min_replicates]
        if len(low) > 0:
            msg = (
                f"Cell type '{ct}': conditions {list(low.index)} have <{min_replicates} "
                f"replicates ({dict(low)}). Skipping."
            )
            all_warnings.append(msg)
            skipped.append({"cell_type": ct, "reason": msg})
            continue

        if data["dropped"]:
            all_warnings.append(
                f"Cell type '{ct}': samples {data['dropped']} dropped "
                f"(<{min_cells_per_pseudobulk} cells)"
            )

        # --- Run PyDESeq2 ---
        try:
            # PyDESeq2 expects genes × samples → transpose so it's samples × genes
            counts_t = counts.T
            # Ensure metadata index matches count columns
            meta_aligned = meta.loc[
                [c.split("__")[1] if "__" in c else c for c in counts_t.index]
            ]
            # Rebuild index to match
            counts_t.index = meta_aligned.index

            # Build contrast: [factor, test_level, reference_level]
            cond_levels = sorted(meta_aligned["condition"].unique())
            contrast = ["condition", cond_levels[0], cond_levels[-1]]

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                dds = DeseqDataSet(
                    counts=counts_t,
                    metadata=meta_aligned,
                    design=design_formula,
                )
                dds.deseq2()
                stat_res = DeseqStats(dds, contrast=contrast, alpha=alpha)
                stat_res.summary()

            res_df = stat_res.results_df.copy()
            res_df["cell_type"] = ct
            res_df["significant"] = (
                (res_df["padj"] < alpha) & (res_df["log2FoldChange"].abs() > lfc_threshold)
            )

            de_results[ct] = res_df

            if plot_dir:
                p = _volcano_plot(res_df, ct, alpha, lfc_threshold, plot_dir)
                if p:
                    plots.append(p)

        except Exception as e:
            msg = f"Cell type '{ct}': DESeq2 failed — {e}"
            all_warnings.append(msg)
            skipped.append({"cell_type": ct, "reason": msg})
            logger.warning(msg)

    # --- Summary ---
    summary_rows = []
    for ct, df in de_results.items():
        n_sig = int(df["significant"].sum())
        n_up = int(((df["significant"]) & (df["log2FoldChange"] > 0)).sum())
        n_down = int(((df["significant"]) & (df["log2FoldChange"] < 0)).sum())
        summary_rows.append({
            "cell_type": ct,
            "n_genes_tested": len(df),
            "n_significant": n_sig,
            "n_up": n_up,
            "n_down": n_down,
        })

    result = {
        "de_results": {ct: df.to_dict("records") for ct, df in de_results.items()},
        "de_dataframes": de_results,  # keep DataFrames for enrichment
        "summary": summary_rows,
        "skipped_cell_types": skipped,
        "plots": plots,
        "warnings": all_warnings,
        "metrics": {
            "n_cell_types_tested": len(de_results),
            "n_cell_types_skipped": len(skipped),
            "conditions": list(conditions),
            "n_samples": int(n_samples),
        },
        "provenance": {
            "tool_id": "deseq2_pseudobulk",
            "parameters": {
                "cell_type_key": cell_type_key,
                "sample_key": sample_key,
                "condition_key": condition_key,
                "min_cells_per_pseudobulk": min_cells_per_pseudobulk,
                "design_formula": design_formula,
                "alpha": alpha,
            },
            "n_cell_types_tested": len(de_results),
            "n_samples_per_condition": {
                str(c): int(adata.obs.loc[
                    adata.obs[condition_key] == c, sample_key
                ].nunique())
                for c in conditions
            },
        },
    }

    return result


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _volcano_plot(
    df: pd.DataFrame,
    cell_type: str,
    alpha: float,
    lfc_threshold: float,
    plot_dir: str,
) -> str | None:
    """Volcano plot for one cell type's DE results."""
    try:
        d = Path(plot_dir)
        d.mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(figsize=(8, 6))

        # Non-significant
        ns = df[~df["significant"]]
        ax.scatter(ns["log2FoldChange"], -np.log10(ns["pvalue"].clip(lower=1e-300)),
                   c="grey", alpha=0.4, s=4, label="NS")

        # Significant
        sig = df[df["significant"]]
        colors = ["#e74c3c" if x > 0 else "#3498db" for x in sig["log2FoldChange"]]
        ax.scatter(sig["log2FoldChange"], -np.log10(sig["pvalue"].clip(lower=1e-300)),
                   c=colors, alpha=0.7, s=8, label=f"Sig (n={len(sig)})")

        ax.axhline(-np.log10(alpha), ls="--", c="grey", lw=0.8)
        if lfc_threshold > 0:
            ax.axvline(lfc_threshold, ls="--", c="grey", lw=0.8)
            ax.axvline(-lfc_threshold, ls="--", c="grey", lw=0.8)

        ax.set_xlabel("log₂ fold change")
        ax.set_ylabel("-log₁₀(p-value)")
        ax.set_title(f"Pseudobulk DE — {cell_type} ({len(sig)} sig. genes)")
        ax.legend(loc="upper right", fontsize=8)

        safe_name = cell_type.replace("/", "_").replace(" ", "_")
        path = str(d / f"volcano_{safe_name}.png")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return path
    except Exception as e:
        logger.warning("Volcano plot failed for %s: %s", cell_type, e)
        return None


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _get_raw_counts(adata: ad.AnnData, layer: str | None) -> Any:
    """Return the raw-count matrix (sparse or dense)."""
    if layer is not None:
        return adata.layers[layer]
    # Heuristic: if .X has non-integer values, try .raw
    if sp.issparse(adata.X):
        sample = adata.X[:100].toarray()
    else:
        sample = np.asarray(adata.X[:100])
    is_integer = np.allclose(sample, sample.astype(int))
    if is_integer:
        return adata.X
    if adata.raw is not None:
        return adata.raw.X
    logger.warning(
        "adata.X appears log-normalized and no .raw found. "
        "Using .X anyway — DE results may be unreliable."
    )
    return adata.X


def _get_var_names(adata: ad.AnnData, layer: str | None) -> pd.Index:
    """Return gene names matching the count matrix."""
    if layer is not None:
        return adata.var_names
    # If we fell back to .raw, use its var_names
    if sp.issparse(adata.X):
        sample = adata.X[:100].toarray()
    else:
        sample = np.asarray(adata.X[:100])
    is_integer = np.allclose(sample, sample.astype(int))
    if is_integer:
        return adata.var_names
    if adata.raw is not None:
        return adata.raw.var_names
    return adata.var_names
