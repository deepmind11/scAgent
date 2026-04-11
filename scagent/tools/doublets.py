"""Doublet detection using Scanpy's built-in Scrublet wrapper."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


def detect_doublets(
    adata: ad.AnnData,
    *,
    expected_doublet_rate: float | None = None,
    n_prin_comps: int = 30,
    sim_doublet_ratio: float = 2.0,
    random_state: int = 0,
    remove: bool = True,
    plot_dir: str | None = None,
    checkpoint_dir: str | None = None,
) -> tuple[ad.AnnData, dict]:
    """Detect doublets using Scrublet (via ``sc.pp.scrublet``).

    Operates on raw counts — must be called BEFORE normalization.

    Parameters
    ----------
    adata
        AnnData with raw (unnormalized) counts in ``adata.X``.
    expected_doublet_rate
        Expected doublet fraction. If *None*, auto-calculated from cell count
        using the 10x Chromium rate (~0.8% per 1,000 cells captured).
    n_prin_comps
        Number of principal components for Scrublet's internal embedding.
    sim_doublet_ratio
        Ratio of simulated doublets to observed cells.
    random_state
        Random seed for reproducibility.
    remove
        If *True*, remove predicted doublets. If *False*, annotate only.
    plot_dir
        Directory to save the doublet score histogram.
    checkpoint_dir
        If provided, save checkpoint after this step.

    Returns
    -------
    (adata, result) — adata with doublets handled; result dict with metrics.
    """
    # Auto-calculate expected rate if not provided
    if expected_doublet_rate is None:
        expected_doublet_rate = min(adata.n_obs / 1000 * 0.008, 0.15)

    cells_before = adata.n_obs

    # Run Scrublet via Scanpy's wrapper (avoids direct annoy dependency issues)
    sc.pp.scrublet(
        adata,
        expected_doublet_rate=expected_doublet_rate,
        sim_doublet_ratio=sim_doublet_ratio,
        n_prin_comps=n_prin_comps,
        random_state=random_state,
        verbose=False,
    )

    # Scanpy stores results in adata.obs['predicted_doublet'] and adata.obs['doublet_score']
    doublet_scores = adata.obs["doublet_score"].values
    predicted_doublets = adata.obs["predicted_doublet"].values

    n_doublets = int(predicted_doublets.sum())
    doublet_rate = n_doublets / cells_before

    # Get threshold from uns if available
    threshold = float(adata.uns.get("scrublet", {}).get("threshold", np.nan))

    # Generate plot
    plots = []
    if plot_dir is not None:
        plots = _plot_doublet_histogram(
            doublet_scores, threshold, n_doublets, cells_before, plot_dir
        )

    # Optionally remove doublets
    if remove and n_doublets > 0:
        adata = adata[~adata.obs["predicted_doublet"], :].copy()

    cells_after = adata.n_obs

    result = {
        "metrics": {
            "cells_before": cells_before,
            "cells_after": cells_after,
            "n_doublets": n_doublets,
            "doublet_rate": round(doublet_rate, 4),
            "expected_doublet_rate": round(expected_doublet_rate, 4),
            "threshold": round(threshold, 4) if not np.isnan(threshold) else None,
            "doublets_removed": remove,
        },
        "plots": plots,
        "provenance": {
            "tool_id": "scrublet_doublets",
            "parameters": {
                "expected_doublet_rate": expected_doublet_rate,
                "n_prin_comps": n_prin_comps,
                "sim_doublet_ratio": sim_doublet_ratio,
                "random_state": random_state,
                "remove": remove,
            },
            "threshold_used": round(threshold, 4) if not np.isnan(threshold) else None,
            "n_doublets": n_doublets,
        },
        "warnings": [],
    }

    # Flag suspicious doublet rates
    if doublet_rate > 2 * expected_doublet_rate:
        result["warnings"].append(
            f"Detected doublet rate ({doublet_rate:.1%}) is >2× the expected rate "
            f"({expected_doublet_rate:.1%}). Inspect the score histogram — the threshold "
            "may need manual adjustment."
        )
    if doublet_rate < expected_doublet_rate / 3:
        result["warnings"].append(
            f"Detected doublet rate ({doublet_rate:.1%}) is much lower than expected "
            f"({expected_doublet_rate:.1%}). Scrublet may have failed to detect doublets. "
            "Check the bimodality of the score histogram."
        )

    if checkpoint_dir is not None:
        _save_checkpoint(adata, checkpoint_dir, "after_doublets")

    return adata, result


def _plot_doublet_histogram(
    scores: np.ndarray,
    threshold: float,
    n_doublets: int,
    n_cells: int,
    plot_dir: str,
) -> list[str]:
    d = Path(plot_dir)
    d.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(scores, bins=100, edgecolor="none", alpha=0.8)
    if not np.isnan(threshold):
        ax.axvline(
            threshold, color="red", linestyle="--", linewidth=1.5,
            label=f"threshold = {threshold:.3f}",
        )
    ax.set_xlabel("Doublet score")
    ax.set_ylabel("Number of cells")
    rate_pct = n_doublets / n_cells * 100
    ax.set_title(f"Scrublet Doublet Scores — {n_doublets} doublets ({rate_pct:.1f}%)")
    ax.legend()
    fig.tight_layout()

    path = str(d / "doublet_histogram.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return [path]


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
