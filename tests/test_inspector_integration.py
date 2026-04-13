"""Integration tests for scagent.inspector using real scBench eval data."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from scagent.inspector import (
    inspect_adata,
    find_raw_counts,
    RAW_COUNTS,
    LOG_NORMALIZED,
    SCALED,
)

# Eval data cache — these are the actual scBench .h5ad files
CACHE_DIR = Path(__file__).parent.parent / ".eval_cache" / "cache"

# Map: description → filename (from manifest)
_DATA_FILES = {
    "qc_raw": "332388f86e1d8227__157798549.node",
    "norm_raw": "bca9e0cc52adc8c7__157798550.node",
    "hvg_scaled": "e18a87d8c3f0e75d__157798551.node",
    "trajectory_caf": "243ae62a228e9a81__157798553.node",
    "clustering_annotated": "a3d3d7ddb0e39418__158379062.node",
}


def _load(key: str) -> ad.AnnData:
    path = CACHE_DIR / _DATA_FILES[key]
    if not path.exists():
        pytest.skip(f"Eval data not available: {path}")
    return ad.read_h5ad(path)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestInspectRealData:
    def test_inspect_raw_data(self):
        """QC eval data — raw integer counts, no preprocessing."""
        adata = _load("qc_raw")
        state = inspect_adata(adata)

        assert state.x_state == RAW_COUNTS
        assert state.raw_counts_location == "X"
        assert state.n_cells == 6420
        assert state.n_genes == 19973
        assert state.species == "mouse"
        assert not state.has_qc_metrics
        assert not state.has_normalized
        assert not state.has_pca
        assert not state.has_clusters

    def test_inspect_scaled_data(self):
        """HVG eval data — z-score scaled with raw_counts in layers."""
        adata = _load("hvg_scaled")
        state = inspect_adata(adata)

        assert state.x_state == SCALED
        assert state.raw_counts_location == "layers['raw_counts']"
        assert state.n_cells == 6420
        assert state.has_normalized  # uns['log1p'] is present
        assert not state.has_pca
        assert not state.has_clusters

    def test_inspect_annotated_data(self):
        """Clustering eval data — scaled, has cell_type labels."""
        adata = _load("clustering_annotated")
        state = inspect_adata(adata)

        assert state.x_state == SCALED
        assert state.raw_counts_location == "layers['raw_counts']"
        assert state.has_cell_types
        assert state.celltype_key == "cell_type"
        assert state.has_normalized

    def test_inspect_trajectory_data(self):
        """Trajectory eval data — raw counts, CAF subset."""
        adata = _load("trajectory_caf")
        state = inspect_adata(adata)

        assert state.x_state == RAW_COUNTS
        assert state.raw_counts_location == "X"
        assert state.n_cells == 1689
        assert not state.has_clusters
        assert not state.has_cell_types

    def test_ensure_ready_for_hvg(self):
        """Scaled data → recover raw counts → normalize → HVG works."""
        adata = _load("hvg_scaled")
        state = inspect_adata(adata)

        # Recover raw counts
        raw = find_raw_counts(adata, state)
        assert (raw >= 0).all()
        sample = raw[:100].flatten()[:1000]
        assert np.allclose(sample, np.round(sample))  # integers

        # Create fresh AnnData from raw counts and run HVG
        import scanpy as sc

        adata_fresh = ad.AnnData(
            X=raw,
            obs=adata.obs.copy(),
            var=adata.var.copy(),
        )
        sc.pp.normalize_total(adata_fresh, target_sum=1e4)
        sc.pp.log1p(adata_fresh)
        sc.pp.highly_variable_genes(adata_fresh, n_top_genes=2000)

        n_hvg = adata_fresh.var["highly_variable"].sum()
        assert n_hvg == 2000  # should work without crashing
