"""Cell-cell communication inference via LIANA+.

LIANA+ is a meta-framework that wraps multiple L-R methods (CellPhoneDB,
NATMI, Connectome, SingleCellSignalR, etc.) and provides consensus rankings.

Best-practice references:
  [BP-1] Heumos et al. 2023, §"Communication events across cells" (p. 557)
  [BP-2] sc-best-practices.org Ch. 22 (Cell-cell communication)
  Dimitrov et al. 2022/2024 — LIANA: L-R database bias, consensus approach

IMPORTANT: L-R interaction databases are biased toward specific pathways,
functional categories, and tissue-enriched proteins [BP-1]. Choice of method
AND database strongly affects predictions. LIANA's consensus mitigates this.

Usage::

    from scagent.tools.communication import run_liana
    result = run_liana(adata, cell_type_key="cell_type")
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
# LIANA+ — consensus cell-cell communication
# ---------------------------------------------------------------------------


def run_liana(
    adata: ad.AnnData,
    *,
    cell_type_key: str = "cell_type",
    methods: list[str] | None = None,
    resource_name: str = "consensus",
    n_perms: int = 1000,
    use_raw: bool = False,
    return_all_lrs: bool = False,
    min_cells: int = 10,
    top_n: int = 50,
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run LIANA+ cell-cell communication analysis.

    Parameters
    ----------
    adata
        Must have cell type annotations.
    cell_type_key
        Column in ``adata.obs`` with cell-type labels.
    methods
        LIANA methods to run. Default: all available (CellPhoneDB,
        NATMI, Connectome, etc.). The consensus is computed automatically.
    resource_name
        L-R resource/database. Default ``"consensus"`` (LIANA's curated set).
    n_perms
        Number of permutations for p-value estimation.
    use_raw
        Whether to use ``adata.raw`` for expression values.
    return_all_lrs
        Return all L-R pairs or only significant ones.
    min_cells
        Minimum cells per cell type to include in analysis.
    top_n
        Number of top interactions to report in summary.
    plot_dir
        Directory for communication plots.

    Returns
    -------
    Dict with ``interactions`` (top ranked), ``summary``, ``provenance``,
    ``plots``, ``warnings``.
    """
    try:
        import liana as li
    except ImportError:
        raise ImportError(
            "liana is required for cell-cell communication. "
            "Install with: pip install liana"
        )

    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate ---
    if cell_type_key not in adata.obs.columns:
        raise ValueError(
            f"Cell type key '{cell_type_key}' not in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    cell_types = adata.obs[cell_type_key].unique()
    n_cell_types = len(cell_types)
    if n_cell_types < 2:
        raise ValueError(
            "Need ≥2 cell types for communication analysis. "
            f"Found {n_cell_types} in '{cell_type_key}'."
        )

    # Check min_cells per cell type
    ct_counts = adata.obs[cell_type_key].value_counts()
    small_cts = ct_counts[ct_counts < min_cells]
    if len(small_cts) > 0:
        result_warnings.append(
            f"{len(small_cts)} cell types have <{min_cells} cells and "
            f"will have limited power: {list(small_cts.index[:5])}"
        )

    # --- L-R database bias warning [BP-1] ---
    result_warnings.append(
        "Note: L-R interaction databases are biased toward specific pathways "
        "and tissue-enriched proteins [BP-1]. LIANA's consensus resource "
        "mitigates but does not eliminate this bias."
    )

    # --- Run LIANA ---
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if methods is not None:
            # Run specific methods
            for method_name in methods:
                method_fn = getattr(li.method, method_name, None)
                if method_fn is None:
                    result_warnings.append(
                        f"Method '{method_name}' not found in liana.method"
                    )
                    continue
                method_fn(
                    adata,
                    groupby=cell_type_key,
                    resource_name=resource_name,
                    use_raw=use_raw,
                    n_perms=n_perms,
                    return_all_lrs=return_all_lrs,
                    verbose=False,
                )
        else:
            # Run all methods via rank_aggregate (consensus)
            li.mt.rank_aggregate(
                adata,
                groupby=cell_type_key,
                resource_name=resource_name,
                use_raw=use_raw,
                n_perms=n_perms,
                return_all_lrs=return_all_lrs,
                verbose=False,
            )

    # --- Extract results ---
    if hasattr(adata, "uns") and "liana_res" in adata.uns:
        liana_res = adata.uns["liana_res"]
    else:
        # Try obsm or other storage locations depending on liana version
        liana_res = getattr(adata, "uns", {}).get("liana_res", pd.DataFrame())

    if isinstance(liana_res, pd.DataFrame) and len(liana_res) > 0:
        # Sort by aggregate rank if available
        sort_col = None
        for candidate in ["magnitude_rank", "specificity_rank", "rank_aggregate"]:
            if candidate in liana_res.columns:
                sort_col = candidate
                break

        if sort_col:
            top_interactions = liana_res.nsmallest(top_n, sort_col)
        else:
            top_interactions = liana_res.head(top_n)

        # Convert to list of dicts for JSON serialization
        interaction_records = []
        for _, row in top_interactions.iterrows():
            record = {}
            for col in ["source", "target", "ligand_complex", "receptor_complex",
                        "magnitude_rank", "specificity_rank", "rank_aggregate"]:
                if col in row.index:
                    val = row[col]
                    record[col] = float(val) if isinstance(val, (np.floating, float)) else str(val)
            interaction_records.append(record)

        n_total = len(liana_res)
        n_pairs = len(
            liana_res[["source", "target"]].drop_duplicates()
        ) if "source" in liana_res.columns else 0
    else:
        interaction_records = []
        n_total = 0
        n_pairs = 0
        result_warnings.append(
            "LIANA returned no results. Check cell type annotations "
            "and ensure sufficient cells per type."
        )

    # --- Plots ---
    if plot_dir and len(interaction_records) > 0:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)
        try:
            li.pl.dotplot(
                adata,
                colour="magnitude_rank" if "magnitude_rank" in liana_res.columns else None,
                size="specificity_rank" if "specificity_rank" in liana_res.columns else None,
                source_labels=list(cell_types[:5]),
                target_labels=list(cell_types[:5]),
                top_n=min(top_n, 20),
            )
            path = str(pdir / "liana_dotplot.png")
            plt.savefig(path, dpi=150, bbox_inches="tight")
            plots.append(path)
            plt.close("all")
        except Exception as e:
            logger.warning("LIANA dotplot failed: %s", e)

    return {
        "interactions": interaction_records,
        "n_total_interactions": n_total,
        "n_cell_type_pairs": n_pairs,
        "n_cell_types": n_cell_types,
        "top_n": top_n,
        "summary": (
            f"LIANA: {n_total} L-R interactions identified across "
            f"{n_pairs} cell-type pairs ({n_cell_types} cell types). "
            f"Top {min(top_n, len(interaction_records))} shown."
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "liana",
            "parameters": {
                "cell_type_key": cell_type_key,
                "resource_name": resource_name,
                "n_perms": n_perms,
                "use_raw": use_raw,
            },
            "n_total_interactions": n_total,
            "n_cell_type_pairs": n_pairs,
            "n_cell_types": n_cell_types,
        },
    }
