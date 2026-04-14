"""Trajectory eval: infer trajectory, find terminal markers.

Uses the scAgent trajectory stack:
  1. Preprocess (inspector + normalize)
  2. PAGA for topology
  3. Palantir for pseudotime (recommended by sc-best-practices.org Ch. 14)
  4. DE between terminal groups

Tries both extremes of DC0 as root cell and picks the trajectory with
stronger terminal group separation (higher mean |log fold change|).
"""

from __future__ import annotations
import warnings
import numpy as np
import anndata as ad
import scanpy as sc

warnings.filterwarnings("ignore")


def handle_trajectory(adata: ad.AnnData, task_prompt: str) -> dict:
    """Infer trajectory and find markers distinguishing the two terminal groups."""
    import palantir

    # --- Step 1: Preprocess ---
    from scagent.inspector import inspect_adata
    from scagent.dependencies import ensure_ready_for

    state = inspect_adata(adata)
    ensure_ready_for(adata, state, needs="log_normalized")
    adata.raw = adata.copy()

    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.tl.pca(adata, n_comps=30, random_state=0)
    sc.pp.neighbors(adata, n_neighbors=15, random_state=0)
    sc.tl.leiden(adata, resolution=0.5)

    # --- Step 2: PAGA topology ---
    from scagent.tools.trajectory import run_paga
    run_paga(adata, groups="leiden")

    # --- Step 3: Palantir pseudotime ---
    palantir.utils.run_diffusion_maps(adata, n_components=10)
    palantir.utils.determine_multiscale_space(adata)

    dc0 = adata.obsm["DM_EigenVectors"][:, 0]

    # Try both extremes as root — pick the one with stronger DE separation
    best_genes = []
    best_score = -1

    for root_idx in [int(np.argmin(dc0)), int(np.argmax(dc0))]:
        root_cell = adata.obs_names[root_idx]

        palantir.utils.run_palantir(
            adata, early_cell=root_cell, num_waypoints=500,
        )

        pt = adata.obs["palantir_pseudotime"].values
        low_thresh = np.percentile(pt, 10)
        high_thresh = np.percentile(pt, 90)

        adata.obs["terminal_group"] = "middle"
        adata.obs.loc[pt <= low_thresh, "terminal_group"] = "terminal_A"
        adata.obs.loc[pt >= high_thresh, "terminal_group"] = "terminal_B"
        adata.obs["terminal_group"] = adata.obs["terminal_group"].astype("category")

        # --- Step 4: DE between terminal groups ---
        terminal_mask = adata.obs["terminal_group"].isin(["terminal_A", "terminal_B"])
        adata_terminal = adata[terminal_mask].copy()

        sc.tl.rank_genes_groups(
            adata_terminal, groupby="terminal_group",
            method="wilcoxon", use_raw=True, n_genes=50,
        )

        # Score: mean absolute log fold change of top DE genes
        # Higher = better separation between terminal groups
        top_genes = []
        total_lfc = 0
        for group in ["terminal_A", "terminal_B"]:
            try:
                names = list(adata_terminal.uns["rank_genes_groups"]["names"][group][:25])
                lfcs = list(adata_terminal.uns["rank_genes_groups"]["logfoldchanges"][group][:25])
                top_genes.extend(names)
                total_lfc += sum(abs(float(x)) for x in lfcs)
            except (KeyError, IndexError):
                pass

        seen = set()
        unique_genes = [g for g in top_genes if g not in seen and not seen.add(g)]
        score = total_lfc

        if score > best_score:
            best_score = score
            best_genes = unique_genes

    return {"top_marker_genes": best_genes}
