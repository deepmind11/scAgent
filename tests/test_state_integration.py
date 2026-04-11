#!/usr/bin/env python3
"""Integration test: branched analysis on real PBMC 10k dataset.

Scenario:
  1. Load → QC → normalize → HVG → PCA → neighbors (on main)
  2. Fork "high_res" from main
  3. On high_res: leiden(2.0)
  4. Switch back to main
  5. On main: leiden(0.8)
  6. Verify: both branches have different clustering, correct state
  7. Verify: provenance has both branches with real hashes
  8. Verify: disk usage is ~2 snapshots, not 12

Requires: data/pbmc10k/filtered_feature_bc_matrix.h5
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from scagent.provenance import ProvenanceGraph, record_step
from scagent.state import StateManager

DATA_FILE = Path("data/pbmc10k/filtered_feature_bc_matrix.h5")

pytestmark = pytest.mark.skipif(
    not DATA_FILE.exists(),
    reason=f"Dataset not found: {DATA_FILE}",
)


def test_branched_analysis():
    from scagent.tools.loading import load_10x_h5
    from scagent.tools.qc import calculate_qc_metrics, filter_cells, filter_genes
    from scagent.tools.normalize import log_normalize
    from scagent.tools.feature_selection import select_hvg
    from scagent.tools.pca import run_pca
    from scagent.tools.neighbors import compute_neighbors
    from scagent.tools.clustering import run_leiden

    with tempfile.TemporaryDirectory() as tmp:
        proj = Path(tmp) / ".scagent"
        sm = StateManager(proj)
        graph = ProvenanceGraph(proj)

        # === Phase 1: Shared pipeline on main ===
        adata, r = load_10x_h5(str(DATA_FILE))
        sm.mark_dirty("after_load")
        h = sm.save_snapshot(adata, "after_load", step_index=0)
        record_step(graph, r, output_hash=h, user_prompt="load")

        r_qc = calculate_qc_metrics(adata)
        sm.mark_dirty("after_qc_metrics")
        record_step(graph, r_qc, user_prompt="QC metrics")

        adata, r = filter_cells(adata, min_genes=200, max_genes=5000, max_pct_mito=20.0)
        sm.mark_dirty("after_filter_cells")
        record_step(graph, r, user_prompt="filter cells")

        adata, r = filter_genes(adata, min_cells=3)
        sm.mark_dirty("after_filter_genes")
        record_step(graph, r, user_prompt="filter genes")

        r = log_normalize(adata)
        sm.mark_dirty("after_normalize")
        record_step(graph, r, user_prompt="normalize")

        r = select_hvg(adata, n_top_genes=2000)
        sm.mark_dirty("after_hvg")
        record_step(graph, r, user_prompt="select HVGs")

        r = run_pca(adata, n_comps=50)
        sm.mark_dirty("after_pca")
        record_step(graph, r, user_prompt="PCA")

        r = compute_neighbors(adata, n_neighbors=15, n_pcs=30)
        sm.mark_dirty("after_neighbors")
        h_neighbors = sm.save_snapshot(adata, "after_neighbors", step_index=6)
        record_step(graph, r, output_hash=h_neighbors, user_prompt="neighbors")

        # === Phase 2: Fork high_res ===
        fork_hash = sm.create_branch("high_res", adata=adata)
        graph.fork_branch("high_res", from_branch="main")
        adata_hr = sm.switch_branch("high_res", adata=adata)

        # leiden at resolution 2.0 on high_res
        r_hr = run_leiden(adata_hr, resolution=2.0)
        sm.mark_dirty("after_leiden_2.0")
        h_hr = sm.save_snapshot(adata_hr, "after_leiden_2.0", step_index=7)
        record_step(graph, r_hr, output_hash=h_hr, branch="high_res",
                     user_prompt="cluster at resolution 2.0")

        n_clusters_hr = adata_hr.obs["leiden"].nunique()

        # === Phase 3: Switch back to main, leiden at 0.8 ===
        adata_main = sm.switch_branch("main", adata=adata_hr)

        r_main = run_leiden(adata_main, resolution=0.8)
        sm.mark_dirty("after_leiden_0.8")
        h_main = sm.save_snapshot(adata_main, "after_leiden_0.8", step_index=7)
        record_step(graph, r_main, output_hash=h_main, branch="main",
                     user_prompt="cluster at resolution 0.8")

        n_clusters_main = adata_main.obs["leiden"].nunique()

        # === Verification ===

        # 1. Different clustering results
        assert n_clusters_hr > n_clusters_main, (
            f"high_res ({n_clusters_hr} clusters) should have more than "
            f"main ({n_clusters_main} clusters)"
        )
        print(f"\n  main: {n_clusters_main} clusters (resolution 0.8)")
        print(f"  high_res: {n_clusters_hr} clusters (resolution 2.0)")

        # 2. Both branches exist with correct state
        branches = {b.name: b for b in sm.list_branches()}
        assert "main" in branches
        assert "high_res" in branches
        assert branches["main"].head_step == "after_leiden_0.8"
        assert branches["high_res"].head_step == "after_leiden_2.0"
        assert branches["high_res"].parent_branch == "main"

        # 3. Provenance has both branches
        main_chain = graph.get_chain("main")
        hr_chain = graph.get_chain("high_res")
        assert len(main_chain) >= 7  # shared steps + leiden
        assert len(hr_chain) >= 1    # at least leiden

        # 4. Full chain for high_res includes parent prefix
        full_hr = graph.get_full_chain("high_res")
        assert len(full_hr) > len(hr_chain)

        # 5. Disk usage: should be ~2-3 snapshots, not 12
        snap_count = 0
        for b in sm.list_branches():
            snap_count += b.snapshot_count
        assert snap_count <= 4, f"Expected ≤4 snapshots, got {snap_count}"

        total_mb = sm.total_disk_usage() / (1024 * 1024)
        print(f"  Total disk: {total_mb:.0f} MB across {snap_count} snapshots")

        # 6. Switching back loads the right data
        adata_check = sm.switch_branch("high_res")
        assert adata_check.obs["leiden"].nunique() == n_clusters_hr

        adata_check2 = sm.switch_branch("main")
        assert adata_check2.obs["leiden"].nunique() == n_clusters_main

        print("  ✓ Branch switching loads correct state")
