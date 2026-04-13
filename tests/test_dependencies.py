"""Unit tests for scagent.dependencies."""

from __future__ import annotations

import numpy as np
import pytest
import anndata as ad

from scagent.inspector import AnnDataState, RAW_COUNTS, LOG_NORMALIZED, SCALED, inspect_adata
from scagent.dependencies import check_prerequisites, plan_steps, ensure_ready_for


# ---------------------------------------------------------------------------
# Helpers — build states directly
# ---------------------------------------------------------------------------


def _state(**overrides) -> AnnDataState:
    """Build an AnnDataState with defaults (everything empty/False)."""
    defaults = dict(
        x_state=RAW_COUNTS,
        raw_counts_location="X",
        has_qc_metrics=False,
        has_normalized=False,
        has_raw_frozen=False,
        has_hvg=False,
        has_pca=False,
        has_neighbors=False,
        has_umap=False,
        has_tsne=False,
        has_clusters=False,
        has_cell_types=False,
        has_de_results=False,
        has_diffmap=False,
        has_velocity_layers=False,
        cluster_key=None,
        celltype_key=None,
        condition_key=None,
        batch_key=None,
        n_cells=1000,
        n_genes=5000,
        species="human",
        warnings=[],
    )
    defaults.update(overrides)
    return AnnDataState(**defaults)


# ---------------------------------------------------------------------------
# check_prerequisites
# ---------------------------------------------------------------------------


class TestCheckPrerequisites:
    def test_prerequisites_met(self):
        """Clustering with neighbors present → can run."""
        state = _state(has_neighbors=True)
        can_run, missing = check_prerequisites("clustering", state)
        assert can_run
        assert missing == []

    def test_prerequisites_missing(self):
        """Clustering without neighbors → cannot run."""
        state = _state(has_neighbors=False)
        can_run, missing = check_prerequisites("clustering", state)
        assert not can_run
        assert len(missing) > 0
        assert any("neighbors" in m for m in missing)

    def test_umap_needs_neighbors(self):
        state = _state(has_neighbors=False)
        can_run, missing = check_prerequisites("umap", state)
        assert not can_run

    def test_umap_with_neighbors(self):
        state = _state(has_neighbors=True)
        can_run, _ = check_prerequisites("umap", state)
        assert can_run

    def test_pseudobulk_needs_raw_counts(self):
        """Pseudobulk DE without raw counts → cannot run."""
        state = _state(
            x_state=SCALED,
            raw_counts_location=None,
            has_cell_types=True,
            condition_key="condition",
        )
        can_run, missing = check_prerequisites("pseudobulk_de", state)
        assert not can_run
        assert any("raw counts" in m for m in missing)

    def test_pseudobulk_needs_condition(self):
        """Pseudobulk DE without condition column → cannot run."""
        state = _state(
            has_cell_types=True,
            raw_counts_location="layers['counts']",
            condition_key=None,
        )
        can_run, missing = check_prerequisites("pseudobulk_de", state)
        assert not can_run
        assert any("condition" in m for m in missing)

    def test_pseudobulk_all_met(self):
        state = _state(
            has_cell_types=True,
            raw_counts_location="layers['counts']",
            condition_key="condition",
        )
        can_run, _ = check_prerequisites("pseudobulk_de", state)
        assert can_run

    def test_hvg_needs_log_normalized(self):
        state = _state(x_state=RAW_COUNTS, has_normalized=True)
        can_run, missing = check_prerequisites("hvg", state)
        assert not can_run
        assert any("log_normalized" in m for m in missing)

    def test_unknown_step_passes(self):
        """Unknown step ID → assume OK (don't block)."""
        state = _state()
        can_run, missing = check_prerequisites("some_future_step", state)
        assert can_run
        assert missing == []


# ---------------------------------------------------------------------------
# plan_steps
# ---------------------------------------------------------------------------


class TestPlanSteps:
    def test_nothing_needed(self):
        """Goal already completed → empty plan."""
        state = _state(has_umap=True)
        steps = plan_steps("umap", state)
        assert steps == []

    def test_fill_gaps_simple(self):
        """Need clustering, have neighbors → just clustering."""
        state = _state(has_neighbors=True)
        steps = plan_steps("clustering", state)
        assert steps == ["clustering"]

    def test_fill_gaps_chain(self):
        """Need clustering, have PCA → neighbors + clustering."""
        state = _state(has_pca=True)
        steps = plan_steps("clustering", state)
        assert steps == ["neighbors", "clustering"]

    def test_fill_gaps_long_chain(self):
        """Need annotation, have HVG → pca + neighbors + clustering + annotation."""
        state = _state(x_state=LOG_NORMALIZED, has_hvg=True)
        steps = plan_steps("annotation", state)
        assert steps == ["pca", "neighbors", "clustering", "annotation"]

    def test_skips_completed(self):
        """PCA done, need UMAP → neighbors + umap only."""
        state = _state(has_pca=True)
        steps = plan_steps("umap", state)
        assert steps == ["neighbors", "umap"]

    def test_everything_done(self):
        """All prerequisites met → just the goal."""
        state = _state(has_neighbors=True, has_clusters=True)
        steps = plan_steps("marker_genes", state)
        assert steps == ["marker_genes"]

    def test_normalize_from_raw(self):
        """Need HVG, have raw counts → normalize + hvg."""
        state = _state(x_state=RAW_COUNTS)
        steps = plan_steps("hvg", state)
        assert "normalize" in steps
        assert "hvg" in steps


# ---------------------------------------------------------------------------
# ensure_ready_for
# ---------------------------------------------------------------------------


class TestEnsureReadyFor:
    def test_already_log_normalized(self):
        """X already log-normalized → no-op."""
        rng = np.random.default_rng(42)
        X = rng.uniform(0, 10, size=(100, 50)).astype(np.float32)
        adata = ad.AnnData(X=X)
        adata.uns["log1p"] = {}
        state = inspect_adata(adata)
        result = ensure_ready_for(adata, state, needs="log_normalized")
        assert result is adata  # same object, not a copy
        # X should be unchanged
        assert np.allclose(adata.X, X)

    def test_raw_to_log_normalized(self):
        """Raw counts → normalize + log."""
        rng = np.random.default_rng(42)
        X = rng.poisson(2, size=(100, 50)).astype(np.int64)
        adata = ad.AnnData(X=X.copy())
        state = inspect_adata(adata)
        ensure_ready_for(adata, state, needs="log_normalized")
        # X should now be floats, non-negative, reasonable range
        assert adata.X.dtype in (np.float32, np.float64)
        assert (adata.X >= 0).all()
        assert adata.X.max() < 20
        # Raw counts preserved
        assert "counts" in adata.layers

    def test_scaled_to_log_normalized(self):
        """Scaled X with raw in layers → recover + normalize."""
        rng = np.random.default_rng(42)
        counts = rng.poisson(2, size=(100, 50)).astype(np.int64)
        X = rng.normal(0, 1, size=(100, 50)).astype(np.float32)
        adata = ad.AnnData(X=X)
        adata.layers["raw_counts"] = counts
        state = inspect_adata(adata)
        ensure_ready_for(adata, state, needs="log_normalized")
        assert (adata.X >= 0).all()
        assert adata.X.max() < 20

    def test_scaled_no_raw_raises(self):
        """Scaled X, no raw counts → raises."""
        rng = np.random.default_rng(42)
        X = rng.normal(0, 1, size=(100, 50)).astype(np.float32)
        adata = ad.AnnData(X=X)
        state = inspect_adata(adata)
        with pytest.raises(ValueError, match="No raw counts found"):
            ensure_ready_for(adata, state, needs="log_normalized")

    def test_ensure_raw_counts(self):
        """Recover raw counts from layers."""
        rng = np.random.default_rng(42)
        counts = rng.poisson(2, size=(100, 50)).astype(np.int64)
        X = rng.normal(0, 1, size=(100, 50)).astype(np.float32)
        adata = ad.AnnData(X=X)
        adata.layers["counts"] = counts
        state = inspect_adata(adata)
        ensure_ready_for(adata, state, needs="raw_counts")
        assert np.allclose(adata.X, counts)
