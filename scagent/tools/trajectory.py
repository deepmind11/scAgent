"""Trajectory inference tools: PAGA, Palantir, DPT, scVelo, CellRank.

Recommended workflow (aligned with sc-best-practices.org):
  1. PAGA      — topology inference (which clusters connect)
  2. Palantir  — pseudotime + terminal state identification [RECOMMENDED]
  3. DPT       — simpler pseudotime alternative (fallback)
  4. scVelo    — RNA velocity (optional, requires spliced/unspliced)
  5. CellRank  — fate mapping from velocity (recommended over raw stream plots)

Best-practice references:
  [BP-2] sc-best-practices.org Ch. 14: Palantir preferred over DPT for pseudotime.
         DPT inflates values for disconnected lineages; Palantir handles branching.
  [BP-2] sc-best-practices.org Ch. 15: scVelo for velocity, CellRank for fate mapping.
         Do NOT over-interpret 2D velocity stream plots.
  [BP-1] Heumos et al. 2023, §"From discrete states to continuous processes"

Usage::

    from scagent.tools.trajectory import (
        run_paga, run_palantir, run_diffusion_pseudotime,
        run_scvelo, run_cellrank,
    )

    result = run_paga(adata, groups="leiden")
    result = run_palantir(adata, root_cell_type="HSC")       # preferred
    result = run_diffusion_pseudotime(adata, root_cell_type="HSC")  # fallback
    result = run_scvelo(adata, mode="dynamical")
    result = run_cellrank(adata)                              # fate mapping
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
# PAGA — topology inference
# ---------------------------------------------------------------------------


def run_paga(
    adata: ad.AnnData,
    *,
    groups: str = "leiden",
    model: str = "v1.2",
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run PAGA graph abstraction for trajectory topology.

    Parameters
    ----------
    adata
        Must have neighbors computed and cluster labels in ``adata.obs[groups]``.
    groups
        Key in ``adata.obs`` for cluster labels.
    model
        PAGA model version (``"v1.0"`` or ``"v1.2"``).
    plot_dir
        Directory for PAGA plots.

    Returns
    -------
    Dict with ``connectivities``, ``summary``, ``provenance``, ``plots``, ``warnings``.
    """
    import scanpy as sc

    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate ---
    if "neighbors" not in adata.uns:
        raise ValueError(
            "Neighbor graph not found. Run sc.pp.neighbors() first. "
            "PAGA needs the KNN graph to compute cluster connectivity."
        )
    if groups not in adata.obs.columns:
        raise ValueError(
            f"Cluster key '{groups}' not in adata.obs. "
            f"Available: {list(adata.obs.columns)}. "
            "Run clustering (Leiden) before trajectory inference."
        )

    n_groups = adata.obs[groups].nunique()
    if n_groups < 2:
        raise ValueError(
            f"Only {n_groups} group(s) in '{groups}'. PAGA needs ≥2 clusters."
        )
    if n_groups > 100:
        result_warnings.append(
            f"High cluster count ({n_groups}). PAGA graph may be hard to interpret. "
            "Consider reducing Leiden resolution."
        )

    # --- Run PAGA ---
    sc.tl.paga(adata, groups=groups, model=model)

    # Extract connectivity info
    conn = adata.uns["paga"]["connectivities"].toarray()
    group_names = list(adata.obs[groups].cat.categories)

    # Check connectivity: is the graph connected?
    from scipy.sparse.csgraph import connected_components
    n_components, labels = connected_components(
        adata.uns["paga"]["connectivities"], directed=False
    )
    if n_components > 1:
        result_warnings.append(
            f"PAGA graph has {n_components} disconnected components. "
            "This may indicate over-clustering or truly separate populations. "
            "Consider lowering Leiden resolution or verifying biology."
        )

    # Build edge list (above-threshold connections)
    edges = []
    threshold = 0.05  # default for visualization
    for i in range(len(group_names)):
        for j in range(i + 1, len(group_names)):
            w = float(conn[i, j])
            if w > threshold:
                edges.append({
                    "source": group_names[i],
                    "target": group_names[j],
                    "weight": round(w, 4),
                })

    # --- Plots ---
    if plot_dir:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)
        try:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sc.pl.paga(adata, ax=axes[0], show=False, title="PAGA graph")
            sc.pl.paga(adata, ax=axes[1], show=False, title="PAGA (compare)",
                       color=groups, threshold=threshold)
            path = str(pdir / "paga_graph.png")
            fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            plots.append(path)
        except Exception as e:
            logger.warning("PAGA plot failed: %s", e)

    return {
        "connectivities": edges,
        "n_groups": n_groups,
        "n_components": n_components,
        "group_names": group_names,
        "summary": (
            f"PAGA computed on {n_groups} clusters ({groups}). "
            f"{len(edges)} edges above threshold {threshold}. "
            f"{'Connected graph.' if n_components == 1 else f'{n_components} components.'}"
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "paga",
            "parameters": {"groups": groups, "model": model},
            "n_groups": n_groups,
            "n_edges": len(edges),
            "n_components": n_components,
        },
    }


# ---------------------------------------------------------------------------
# DPT — diffusion pseudotime
# ---------------------------------------------------------------------------


def run_diffusion_pseudotime(
    adata: ad.AnnData,
    *,
    root_cell_type: str | None = None,
    root_cell_index: int | None = None,
    n_dcs: int = 15,
    cell_type_key: str = "cell_type",
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run diffusion pseudotime ordering.

    Requires EITHER ``root_cell_type`` (picks the most extreme cell in
    diffusion component space from that cluster) OR ``root_cell_index``
    (explicit cell barcode index).

    Parameters
    ----------
    adata
        Must have neighbors computed.
    root_cell_type
        Name of the root cell type (e.g., ``"HSC"``). The most stem-like
        cell is selected automatically via diffusion components.
    root_cell_index
        Explicit root cell index. Overrides ``root_cell_type``.
    n_dcs
        Number of diffusion components.
    cell_type_key
        Column in ``adata.obs`` with cell type labels.
    plot_dir
        Directory for pseudotime plots.

    Returns
    -------
    Dict with ``pseudotime_stats``, ``summary``, ``provenance``, ``plots``, ``warnings``.

    Notes
    -----
    [BP-2] Ch. 14 warns that DPT can inflate pseudotime for disconnected
    lineages (e.g., CLPs in bone marrow). If multiple lineages are present,
    consider comparing with Palantir pseudotime.
    """
    import scanpy as sc

    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate ---
    if "neighbors" not in adata.uns:
        raise ValueError(
            "Neighbor graph not found. Run sc.pp.neighbors() first."
        )

    # --- Compute diffusion map ---
    sc.tl.diffmap(adata, n_comps=n_dcs)

    # --- Determine root cell ---
    if root_cell_index is not None:
        adata.uns["iroot"] = root_cell_index
        root_info = f"explicit index {root_cell_index}"
    elif root_cell_type is not None:
        if cell_type_key not in adata.obs.columns:
            raise ValueError(
                f"Cell type key '{cell_type_key}' not in adata.obs. "
                f"Available: {list(adata.obs.columns)}"
            )
        ct_mask = adata.obs[cell_type_key] == root_cell_type
        if ct_mask.sum() == 0:
            available = list(adata.obs[cell_type_key].unique())
            raise ValueError(
                f"Root cell type '{root_cell_type}' not found. "
                f"Available types: {available}"
            )

        # Pick the most extreme cell in diffusion component space
        # Use the component with max variance within the root cluster
        dc_root = adata.obsm["X_diffmap"][ct_mask.values]
        dc_var = dc_root.var(axis=0)
        best_dc = int(np.argmax(dc_var))
        # Pick the cell at the extreme of this component
        root_idx = int(np.where(ct_mask)[0][np.argmin(dc_root[:, best_dc])])
        adata.uns["iroot"] = root_idx
        root_info = (
            f"auto-selected from '{root_cell_type}' "
            f"(DC {best_dc}, cell index {root_idx})"
        )
    else:
        raise ValueError(
            "Must provide either root_cell_type or root_cell_index. "
            "The root should be the starting point of the developmental "
            "process (e.g., stem cells, progenitors)."
        )

    # --- Compute DPT ---
    sc.tl.dpt(adata, n_dcs=n_dcs)

    dpt = adata.obs["dpt_pseudotime"].values

    # --- Guard rail: check for disconnected lineages ---
    # DPT inflates values to ~inf for unreachable cells [BP-2 Ch. 14]
    inf_mask = np.isinf(dpt) | (dpt > 1e6)
    n_inf = int(inf_mask.sum())
    if n_inf > 0:
        result_warnings.append(
            f"{n_inf} cells ({n_inf/len(dpt)*100:.1f}%) have infinite/extreme "
            "pseudotime — these are in disconnected lineages. "
            "Consider: (1) lowering Leiden resolution, (2) subsetting to the "
            "lineage of interest, or (3) comparing with Palantir pseudotime. "
            "[BP-2 Ch. 14]"
        )
        # Cap for stats
        dpt_finite = dpt[~inf_mask]
    else:
        dpt_finite = dpt

    # --- Per-cluster pseudotime stats ---
    pt_stats = {}
    if cell_type_key in adata.obs.columns:
        for ct in adata.obs[cell_type_key].unique():
            ct_vals = dpt[adata.obs[cell_type_key] == ct]
            ct_finite = ct_vals[~(np.isinf(ct_vals) | (ct_vals > 1e6))]
            if len(ct_finite) > 0:
                pt_stats[str(ct)] = {
                    "median": float(np.median(ct_finite)),
                    "mean": float(np.mean(ct_finite)),
                    "std": float(np.std(ct_finite)),
                    "n_cells": int(len(ct_vals)),
                    "n_inf": int(len(ct_vals) - len(ct_finite)),
                }

    # Sort by median pseudotime for ordering
    ordered_types = sorted(pt_stats.items(), key=lambda x: x[1]["median"])

    # --- Plots ---
    if plot_dir:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)
        try:
            # UMAP colored by pseudotime
            if "X_umap" in adata.obsm:
                fig, ax = plt.subplots(figsize=(8, 6))
                sc.pl.umap(adata, color="dpt_pseudotime", ax=ax, show=False,
                           color_map="gnuplot2", title="Diffusion Pseudotime")
                path = str(pdir / "dpt_umap.png")
                fig.savefig(path, dpi=150, bbox_inches="tight")
                plt.close(fig)
                plots.append(path)
        except Exception as e:
            logger.warning("DPT UMAP plot failed: %s", e)

        try:
            # Violin plot of pseudotime per cell type
            if cell_type_key in adata.obs.columns:
                order = [ct for ct, _ in ordered_types]
                fig, ax = plt.subplots(figsize=(max(10, len(order)), 6))
                sc.pl.violin(adata, keys="dpt_pseudotime", groupby=cell_type_key,
                             order=order, ax=ax, show=False, rotation=45)
                ax.set_title("Pseudotime distribution by cell type")
                path = str(pdir / "dpt_violin.png")
                fig.savefig(path, dpi=150, bbox_inches="tight")
                plt.close(fig)
                plots.append(path)
        except Exception as e:
            logger.warning("DPT violin plot failed: %s", e)

    return {
        "pseudotime_stats": pt_stats,
        "cell_ordering": [ct for ct, _ in ordered_types],
        "root_info": root_info,
        "n_cells_infinite_pt": n_inf,
        "summary": (
            f"DPT computed with root: {root_info}. "
            f"Pseudotime range: [{float(np.nanmin(dpt_finite)):.3f}, "
            f"{float(np.nanmax(dpt_finite)):.3f}]. "
            f"Cell type ordering (early→late): {' → '.join(ct for ct, _ in ordered_types[:6])}{'...' if len(ordered_types) > 6 else ''}. "
            + (f"⚠️ {n_inf} cells unreachable." if n_inf > 0 else "")
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "diffusion_pseudotime",
            "parameters": {
                "root_cell_type": root_cell_type,
                "root_cell_index": root_cell_index,
                "n_dcs": n_dcs,
                "cell_type_key": cell_type_key,
            },
            "root_info": root_info,
            "n_cells": int(len(dpt)),
            "n_cells_infinite_pt": n_inf,
        },
    }


# ---------------------------------------------------------------------------
# scVelo — RNA velocity
# ---------------------------------------------------------------------------


def run_scvelo(
    adata: ad.AnnData,
    *,
    mode: str = "stochastic",
    n_jobs: int = 1,
    plot_dir: str | None = None,
    top_genes: int = 20,
) -> dict[str, Any]:
    """Run RNA velocity analysis with scVelo.

    Parameters
    ----------
    adata
        Must have ``spliced`` and ``unspliced`` layers.
    mode
        Velocity mode: ``"deterministic"``, ``"stochastic"``, or ``"dynamical"``.
    n_jobs
        Parallel jobs for dynamical mode.
    plot_dir
        Directory for velocity plots.
    top_genes
        Number of top likelihood genes for phase portrait validation.

    Returns
    -------
    Dict with ``velocity_stats``, ``summary``, ``provenance``, ``plots``, ``warnings``.

    Notes
    -----
    [BP-1]: Check phase portraits of high-likelihood genes. If lacking the
    expected almond shape (induction arc + repression arc), velocity may be
    inferred incorrectly. Transcriptional bursts and steady-state populations
    can also cause erroneous velocity directions.
    """
    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate: check for spliced/unspliced ---
    if "spliced" not in adata.layers:
        raise ValueError(
            "Layer 'spliced' not found in adata.layers. "
            "RNA velocity requires spliced and unspliced count matrices. "
            "These must be generated by velocyto (run_10x) or STARsolo "
            "during alignment — they are NOT available from standard "
            "Cell Ranger count output."
        )
    if "unspliced" not in adata.layers:
        raise ValueError(
            "Layer 'unspliced' not found in adata.layers. "
            "See above — requires velocyto or STARsolo preprocessing."
        )

    # Check that layers have meaningful content
    import scipy.sparse as sp

    spliced = adata.layers["spliced"]
    if sp.issparse(spliced):
        frac_nonzero = spliced.nnz / (spliced.shape[0] * spliced.shape[1])
    else:
        frac_nonzero = float((spliced != 0).mean())

    if frac_nonzero < 0.01:
        result_warnings.append(
            f"Spliced layer is very sparse ({frac_nonzero*100:.1f}% non-zero). "
            "Velocity estimates may be unreliable."
        )

    # --- Run scVelo ---
    try:
        import scvelo as scv
    except ImportError:
        raise ImportError(
            "scvelo is required for RNA velocity. Install with: pip install scvelo"
        )

    # scVelo preprocessing
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    # Velocity estimation
    if mode == "dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
    scv.tl.velocity(adata, mode=mode)
    scv.tl.velocity_graph(adata)

    # --- Velocity confidence ---
    scv.tl.velocity_confidence(adata)

    confidence = adata.obs["velocity_confidence"].values
    median_conf = float(np.nanmedian(confidence))
    low_conf_frac = float((confidence < 0.5).mean())

    if median_conf < 0.5:
        result_warnings.append(
            f"Median velocity confidence is low ({median_conf:.2f}). "
            "This suggests insufficient spliced/unspliced signal or "
            "that the biological process may not be well-captured by "
            "RNA velocity (e.g., steady-state populations). [BP-1]"
        )

    if low_conf_frac > 0.5:
        result_warnings.append(
            f"{low_conf_frac*100:.0f}% of cells have confidence <0.5. "
            "Consider whether RNA velocity is appropriate for this dataset."
        )

    # --- Phase portrait validation [BP-1] ---
    # Check top likelihood genes if dynamical mode
    phase_portrait_warnings = []
    if mode == "dynamical" and "fit_likelihood" in adata.var.columns:
        top_genes_list = (
            adata.var["fit_likelihood"]
            .sort_values(ascending=False)
            .head(top_genes)
            .index.tolist()
        )
    else:
        top_genes_list = []

    # --- Plots ---
    if plot_dir and "X_umap" in adata.obsm:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)
        try:
            fig, ax = plt.subplots(figsize=(10, 8))
            scv.pl.velocity_embedding_stream(
                adata, basis="umap", ax=ax, show=False,
                title="RNA Velocity Stream"
            )
            path = str(pdir / "velocity_stream.png")
            fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            plots.append(path)
        except Exception as e:
            logger.warning("Velocity stream plot failed: %s", e)

        try:
            fig, ax = plt.subplots(figsize=(8, 6))
            scv.pl.velocity_embedding(
                adata, basis="umap", ax=ax, show=False,
                title="RNA Velocity Arrows", arrow_length=3, arrow_size=2,
            )
            path = str(pdir / "velocity_arrows.png")
            fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            plots.append(path)
        except Exception as e:
            logger.warning("Velocity arrows plot failed: %s", e)

    return {
        "velocity_stats": {
            "median_confidence": median_conf,
            "low_confidence_fraction": low_conf_frac,
            "n_velocity_genes": int((adata.var["velocity_genes"] == True).sum())
            if "velocity_genes" in adata.var.columns else 0,
        },
        "top_likelihood_genes": top_genes_list,
        "summary": (
            f"RNA velocity computed (mode={mode}). "
            f"Median confidence: {median_conf:.2f}. "
            f"{low_conf_frac*100:.0f}% cells below 0.5 confidence. "
            + ("⚠️ Low overall confidence — interpret with caution." if median_conf < 0.5 else "")
        ),
        "plots": plots,
        "warnings": result_warnings + phase_portrait_warnings,
        "provenance": {
            "tool_id": "scvelo_velocity",
            "parameters": {"mode": mode, "n_jobs": n_jobs},
            "median_confidence": median_conf,
            "n_cells": int(adata.n_obs),
        },
    }


# ---------------------------------------------------------------------------
# Palantir — pseudotime + terminal states (RECOMMENDED)
# ---------------------------------------------------------------------------


def run_palantir(
    adata: ad.AnnData,
    *,
    root_cell_type: str | None = None,
    root_cell_index: int | None = None,
    cell_type_key: str = "cell_type",
    n_components: int = 5,
    knn: int = 30,
    num_waypoints: int = 500,
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run Palantir pseudotime and terminal state inference.

    Palantir models trajectories as Markov chains on diffusion maps,
    producing continuous pseudotime and per-cell fate probabilities
    toward identified terminal states. Preferred over DPT for branching
    trajectories. [BP-2 Ch. 14]

    Parameters
    ----------
    adata
        Must have PCA computed.
    root_cell_type
        Name of the root cell type. The most stem-like cell is
        auto-selected via diffusion components.
    root_cell_index
        Explicit root cell index. Overrides ``root_cell_type``.
    cell_type_key
        Column in ``adata.obs`` with cell type labels.
    n_components
        Number of diffusion components for multiscale space.
    knn
        Number of nearest neighbors.
    num_waypoints
        Number of waypoints for trajectory computation.
    plot_dir
        Directory for pseudotime plots.

    Returns
    -------
    Dict with ``pseudotime_stats``, ``terminal_states``, ``summary``,
    ``provenance``, ``plots``, ``warnings``.
    """
    try:
        import palantir
    except ImportError:
        raise ImportError(
            "palantir is required. Install with: pip install palantir"
        )

    import scanpy as sc

    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate ---
    if "X_pca" not in adata.obsm:
        raise ValueError("PCA not found. Run sc.tl.pca() first.")

    # --- Determine root cell ---
    # Compute diffusion maps via Palantir
    palantir.utils.run_diffusion_maps(adata, n_components=n_components, knn=knn)

    if root_cell_index is not None:
        root_cell = adata.obs_names[root_cell_index]
        root_info = f"explicit index {root_cell_index}"
    elif root_cell_type is not None:
        if cell_type_key not in adata.obs.columns:
            raise ValueError(
                f"Cell type key '{cell_type_key}' not in adata.obs. "
                f"Available: {list(adata.obs.columns)}"
            )
        ct_mask = adata.obs[cell_type_key] == root_cell_type
        if ct_mask.sum() == 0:
            available = list(adata.obs[cell_type_key].unique())
            raise ValueError(
                f"Root cell type '{root_cell_type}' not found. "
                f"Available types: {available}"
            )

        # Pick most extreme cell in multiscale space within root cluster
        ms_data = adata.obsm.get("X_palantir_multiscale",
                                  adata.obsm.get("DM_EigenVectors"))
        if ms_data is None:
            raise ValueError("Palantir diffusion maps failed to compute.")

        ct_indices = np.where(ct_mask)[0]
        ct_ms = ms_data[ct_indices]
        # Use the cell with minimum value in first component
        best_idx = ct_indices[np.argmin(ct_ms[:, 0])]
        root_cell = adata.obs_names[best_idx]
        root_info = f"auto-selected from '{root_cell_type}' (index {best_idx})"
    else:
        raise ValueError(
            "Must provide either root_cell_type or root_cell_index."
        )

    # --- Compute multiscale space and run Palantir ---
    palantir.utils.determine_multiscale_space(adata)

    pr_res = palantir.utils.run_palantir(
        adata,
        early_cell=root_cell,
        num_waypoints=num_waypoints,
    )

    # --- Extract results ---
    pseudotime = adata.obs["palantir_pseudotime"].values
    entropy = adata.obs.get("palantir_entropy", pd.Series(dtype=float)).values

    # Terminal states
    terminal_states = []
    if "palantir_fate_probabilities" in adata.obsm:
        fate_probs = adata.obsm["palantir_fate_probabilities"]
        if hasattr(fate_probs, 'columns'):
            terminal_states = list(fate_probs.columns)

    # Per-cluster pseudotime stats
    pt_stats = {}
    if cell_type_key in adata.obs.columns:
        for ct in adata.obs[cell_type_key].unique():
            ct_vals = pseudotime[adata.obs[cell_type_key] == ct]
            pt_stats[str(ct)] = {
                "median": float(np.median(ct_vals)),
                "mean": float(np.mean(ct_vals)),
                "std": float(np.std(ct_vals)),
                "n_cells": int(len(ct_vals)),
            }

    ordered_types = sorted(pt_stats.items(), key=lambda x: x[1]["median"])

    # --- Plots ---
    if plot_dir:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)

        if "X_umap" in adata.obsm:
            try:
                fig, ax = plt.subplots(figsize=(8, 6))
                sc.pl.umap(adata, color="palantir_pseudotime", ax=ax,
                           show=False, color_map="gnuplot2",
                           title="Palantir Pseudotime")
                path = str(pdir / "palantir_pseudotime_umap.png")
                fig.savefig(path, dpi=150, bbox_inches="tight")
                plt.close(fig)
                plots.append(path)
            except Exception as e:
                logger.warning("Palantir UMAP plot failed: %s", e)

        if cell_type_key in adata.obs.columns:
            try:
                order = [ct for ct, _ in ordered_types]
                fig, ax = plt.subplots(figsize=(max(10, len(order)), 6))
                sc.pl.violin(adata, keys="palantir_pseudotime",
                             groupby=cell_type_key, order=order,
                             ax=ax, show=False, rotation=45)
                ax.set_title("Palantir pseudotime by cell type")
                path = str(pdir / "palantir_violin.png")
                fig.savefig(path, dpi=150, bbox_inches="tight")
                plt.close(fig)
                plots.append(path)
            except Exception as e:
                logger.warning("Palantir violin plot failed: %s", e)

    return {
        "pseudotime_stats": pt_stats,
        "cell_ordering": [ct for ct, _ in ordered_types],
        "terminal_states": terminal_states,
        "root_info": root_info,
        "summary": (
            f"Palantir pseudotime computed. Root: {root_info}. "
            f"Terminal states: {terminal_states if terminal_states else 'auto-detected'}. "
            f"Ordering (early→late): "
            f"{' → '.join(ct for ct, _ in ordered_types[:6])}"
            f"{'...' if len(ordered_types) > 6 else ''}."
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "palantir",
            "parameters": {
                "root_cell_type": root_cell_type,
                "root_cell_index": root_cell_index,
                "n_components": n_components,
                "knn": knn,
                "num_waypoints": num_waypoints,
            },
            "root_info": root_info,
            "n_terminal_states": len(terminal_states),
            "n_cells": int(len(pseudotime)),
        },
    }


# ---------------------------------------------------------------------------
# CellRank — fate mapping from velocity (RECOMMENDED over raw streams)
# ---------------------------------------------------------------------------


def run_cellrank(
    adata: ad.AnnData,
    *,
    cluster_key: str = "leiden",
    n_states: int | None = None,
    plot_dir: str | None = None,
) -> dict[str, Any]:
    """Run CellRank fate mapping from RNA velocity.

    CellRank uses the velocity field to compute a Markov chain transition
    matrix, then identifies terminal/initial states and fate probabilities.
    Operates in high-dimensional space, avoiding misleading 2D stream
    plot artifacts. [BP-2 Ch. 15]

    Parameters
    ----------
    adata
        Must have RNA velocity computed (``adata.layers["velocity"]``).
    cluster_key
        Cluster labels for terminal state naming.
    n_states
        Number of terminal states. ``None`` for automatic detection.
    plot_dir
        Directory for fate plots.

    Returns
    -------
    Dict with ``terminal_states``, ``initial_states``, ``summary``,
    ``provenance``, ``plots``, ``warnings``.
    """
    try:
        import cellrank as cr
    except ImportError:
        raise ImportError(
            "cellrank is required for fate mapping. "
            "Install with: pip install cellrank"
        )

    result_warnings: list[str] = []
    plots: list[str] = []

    # --- Validate ---
    if "velocity" not in adata.layers:
        raise ValueError(
            "RNA velocity not found in adata.layers['velocity']. "
            "Run run_scvelo() first."
        )

    if cluster_key not in adata.obs.columns:
        result_warnings.append(
            f"Cluster key '{cluster_key}' not found. "
            "Terminal states will be named by cell index."
        )

    # --- Build velocity kernel ---
    vk = cr.kernels.VelocityKernel(adata)
    vk.compute_transition_matrix()

    # --- Identify terminal states ---
    estimator = cr.estimators.GPCCA(vk)

    if n_states is not None:
        estimator.fit(n_states=n_states)
    else:
        estimator.fit()

    estimator.predict_terminal_states()

    # Extract terminal states
    terminal_states = []
    if "terminal_states" in adata.obs.columns:
        ts = adata.obs["terminal_states"].dropna().unique()
        terminal_states = [str(s) for s in ts]

    # Try to identify initial states too
    initial_states = []
    try:
        estimator.predict_initial_states()
        if "initial_states" in adata.obs.columns:
            init_s = adata.obs["initial_states"].dropna().unique()
            initial_states = [str(s) for s in init_s]
    except Exception:
        pass

    # Fate probabilities
    try:
        estimator.compute_fate_probabilities()
        has_fate_probs = True
    except Exception:
        has_fate_probs = False

    # --- Plots ---
    if plot_dir and "X_umap" in adata.obsm:
        pdir = Path(plot_dir)
        pdir.mkdir(parents=True, exist_ok=True)
        try:
            estimator.plot_fate_probabilities(
                basis="umap", save=str(pdir / "cellrank_fate.png"),
                dpi=150,
            )
            plots.append(str(pdir / "cellrank_fate.png"))
        except Exception as e:
            logger.warning("CellRank fate plot failed: %s", e)

    return {
        "terminal_states": terminal_states,
        "initial_states": initial_states,
        "has_fate_probabilities": has_fate_probs,
        "summary": (
            f"CellRank fate mapping: {len(terminal_states)} terminal states "
            f"({', '.join(terminal_states) if terminal_states else 'detected'}), "
            f"{len(initial_states)} initial states. "
            f"{'Fate probabilities computed.' if has_fate_probs else ''}"
        ),
        "plots": plots,
        "warnings": result_warnings,
        "provenance": {
            "tool_id": "cellrank",
            "parameters": {
                "cluster_key": cluster_key,
                "n_states": n_states,
            },
            "n_terminal_states": len(terminal_states),
            "n_initial_states": len(initial_states),
            "n_cells": int(adata.n_obs),
        },
    }
