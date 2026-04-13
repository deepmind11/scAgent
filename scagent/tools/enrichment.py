"""Gene set enrichment analysis (GSEA) via GSEApy.

Takes DE results (from pseudobulk_de) and runs pre-ranked GSEA to
find enriched pathways and biological processes.

Usage::

    from scagent.tools.enrichment import run_gsea
    results = run_gsea(de_results_dict)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def run_gsea(
    de_results: dict[str, pd.DataFrame] | pd.DataFrame,
    *,
    gene_sets: str = "MSigDB_Hallmark_2020",
    ranking_metric: str = "log2fc",
    permutation_num: int = 1000,
    min_size: int = 15,
    max_size: int = 500,
    fdr_threshold: float = 0.25,
    organism: str = "human",
    plot_dir: str | None = None,
) -> dict:
    """Run GSEA on DE results from pseudobulk_de.

    Parameters
    ----------
    de_results
        Either a dict keyed by cell type (each value a DataFrame with
        ``log2FoldChange`` and ``pvalue`` columns), or a single DataFrame.
        Typically the ``de_dataframes`` field from ``run_pseudobulk_de()``.
    gene_sets
        Gene set library name (GSEApy builtin) or path to ``.gmt`` file.
        Common: ``'MSigDB_Hallmark_2020'``, ``'GO_Biological_Process_2023'``,
        ``'KEGG_2021_Human'``, ``'Reactome_2022'``.
    ranking_metric
        ``'log2fc'`` — rank by log2 fold change (default).
        ``'signal_to_noise'`` — rank by −log10(pval) × sign(log2fc).
    permutation_num
        Number of permutations for p-value estimation.
    min_size
        Minimum genes in a gene set to include it.
    max_size
        Maximum genes in a gene set.
    fdr_threshold
        FDR cutoff for reporting significant terms.
    organism
        ``'human'`` or ``'mouse'`` — affects gene set selection for
        some libraries.
    plot_dir
        Directory for enrichment plots.

    Returns
    -------
    Result dict with ``enrichment_results`` (per cell type),
    ``summary``, ``provenance``, ``warnings``, ``plots``.
    """
    import gseapy as gp

    # Normalize input
    if isinstance(de_results, pd.DataFrame):
        de_results = {"all": de_results}

    enrichment_results: dict[str, pd.DataFrame] = {}
    all_warnings: list[str] = []
    plots: list[str] = []
    summary_rows: list[dict] = []

    for ct, df in de_results.items():
        # Build ranked gene list
        rnk = _build_ranking(df, ranking_metric)
        if rnk is None or len(rnk) < 100:
            msg = f"Cell type '{ct}': too few ranked genes ({0 if rnk is None else len(rnk)}). Skipping GSEA."
            all_warnings.append(msg)
            continue

        try:
            pre_res = gp.prerank(
                rnk=rnk,
                gene_sets=gene_sets,
                min_size=min_size,
                max_size=max_size,
                permutation_num=permutation_num,
                seed=42,
                verbose=False,
                no_plot=True,
            )

            res_df = pre_res.res2d.copy()
            res_df["cell_type"] = ct

            # Ensure numeric types
            for col in ("NES", "FDR q-val", "NOM p-val", "FWER p-val"):
                if col in res_df.columns:
                    res_df[col] = pd.to_numeric(res_df[col], errors="coerce")

            fdr_col = "FDR q-val" if "FDR q-val" in res_df.columns else "fdr"
            nes_col = "NES" if "NES" in res_df.columns else "nes"

            n_sig = int((res_df[fdr_col] < fdr_threshold).sum()) if fdr_col in res_df.columns else 0
            n_up = 0
            n_down = 0
            if fdr_col in res_df.columns and nes_col in res_df.columns:
                sig_mask = res_df[fdr_col] < fdr_threshold
                n_up = int((sig_mask & (res_df[nes_col] > 0)).sum())
                n_down = int((sig_mask & (res_df[nes_col] < 0)).sum())

            enrichment_results[ct] = res_df

            summary_rows.append({
                "cell_type": ct,
                "n_terms_tested": len(res_df),
                "n_significant": n_sig,
                "n_up": n_up,
                "n_down": n_down,
            })

            if n_sig == 0:
                all_warnings.append(
                    f"Cell type '{ct}': no significant terms at FDR<{fdr_threshold}. "
                    "Biological signal may be weak or try different gene sets."
                )

        except Exception as e:
            msg = f"Cell type '{ct}': GSEA failed — {e}"
            all_warnings.append(msg)
            logger.warning(msg)

    result = {
        "enrichment_results": {
            ct: df.to_dict("records") for ct, df in enrichment_results.items()
        },
        "enrichment_dataframes": enrichment_results,
        "summary": summary_rows,
        "plots": plots,
        "warnings": all_warnings,
        "metrics": {
            "n_cell_types_tested": len(enrichment_results),
            "gene_sets": gene_sets,
        },
        "provenance": {
            "tool_id": "gsea",
            "parameters": {
                "gene_sets": gene_sets,
                "ranking_metric": ranking_metric,
                "permutation_num": permutation_num,
                "min_size": min_size,
                "max_size": max_size,
                "fdr_threshold": fdr_threshold,
            },
            "n_cell_types_tested": len(enrichment_results),
            "n_significant_terms": sum(r["n_significant"] for r in summary_rows),
        },
    }

    return result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _build_ranking(df: pd.DataFrame, metric: str) -> pd.Series | None:
    """Build a ranked gene list from DE results."""
    if "log2FoldChange" not in df.columns:
        logger.warning("DE results missing 'log2FoldChange' column")
        return None

    if metric == "log2fc":
        rnk = df["log2FoldChange"].dropna()
    elif metric == "signal_to_noise":
        if "pvalue" not in df.columns:
            logger.warning("signal_to_noise metric requires 'pvalue' column")
            return None
        pval = df["pvalue"].clip(lower=1e-300)
        rnk = -np.log10(pval) * np.sign(df["log2FoldChange"])
        rnk = rnk.dropna()
    else:
        logger.warning("Unknown ranking metric '%s', using log2fc", metric)
        rnk = df["log2FoldChange"].dropna()

    # GSEApy expects a Series with gene names as index
    rnk.index = df.loc[rnk.index].index
    rnk = rnk.sort_values(ascending=False)
    # Drop duplicates (keep first / highest)
    rnk = rnk[~rnk.index.duplicated(keep="first")]

    return rnk
