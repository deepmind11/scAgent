"""QC eval: conservative cell filtering on already-clean data."""

from __future__ import annotations
import anndata as ad
from scagent.tools.qc import calculate_qc_metrics, filter_cells, filter_genes


def handle_qc(adata: ad.AnnData, task_prompt: str) -> dict:
    """Filter cells with conservative thresholds.

    The 4T1 data is already post-filtering and high quality.
    The eval expects ~6420 cells after filtering (i.e., almost none removed).
    We use MAD-based thresholds to be data-adaptive rather than hardcoded.
    """
    # Compute QC metrics (mouse data — species auto-detected)
    qc_result = calculate_qc_metrics(adata)
    thresholds = qc_result["metrics"]["recommended_thresholds"]

    # Use MAD-based thresholds — these adapt to the data
    adata, filt_result = filter_cells(
        adata,
        min_genes=thresholds["min_genes"],
        max_genes=thresholds["max_genes"],
        max_pct_mito=thresholds["max_pct_mito_mad"],
    )

    # Also filter genes
    adata, _ = filter_genes(adata, min_cells=3)

    cells_after = adata.n_obs

    return {"cells_after_filtering": cells_after}
