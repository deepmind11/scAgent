"""Unit tests for scagent.inspector."""

from __future__ import annotations

import numpy as np
import pytest
import anndata as ad
from scipy.sparse import csr_matrix

from scagent.inspector import (
    AnnDataState,
    inspect_adata,
    find_raw_counts,
    summarize_state,
    RAW_COUNTS,
    LOG_NORMALIZED,
    SCALED,
    TRANSFORMED_UNKNOWN,
)


# ---------------------------------------------------------------------------
# Helpers — build AnnData objects in known states
# ---------------------------------------------------------------------------


def _make_raw(n_cells=200, n_genes=100, sparse=False):
    """Raw integer UMI counts."""
    rng = np.random.default_rng(42)
    X = rng.poisson(2, size=(n_cells, n_genes)).astype(np.int64)
    if sparse:
        X = csr_matrix(X)
    genes = [f"GENE{i}" for i in range(n_genes)]
    return ad.AnnData(X=X, var={"_": genes})


def _make_log_normalized(n_cells=200, n_genes=100):
    """Log-normalized float data with uns['log1p'] marker."""
    rng = np.random.default_rng(42)
    counts = rng.poisson(2, size=(n_cells, n_genes)).astype(np.float32)
    # Simulate normalize_total + log1p
    row_sums = counts.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    X = np.log1p(counts / row_sums * 1e4)
    adata = ad.AnnData(X=X)
    adata.uns["log1p"] = {}
    return adata


def _make_scaled(n_cells=200, n_genes=100):
    """Z-score scaled data with raw counts in layers."""
    rng = np.random.default_rng(42)
    counts = rng.poisson(2, size=(n_cells, n_genes)).astype(np.int64)
    # Normalize, log, scale
    row_sums = counts.sum(axis=1, keepdims=True).astype(np.float64)
    row_sums[row_sums == 0] = 1
    X = np.log1p(counts / row_sums * 1e4)
    X = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-8)  # z-score
    adata = ad.AnnData(X=X.astype(np.float32))
    adata.layers["raw_counts"] = counts
    adata.uns["log1p"] = {}
    return adata


def _make_transformed_unknown(n_cells=200, n_genes=100):
    """Floats, non-negative, no log1p marker, max > 20."""
    rng = np.random.default_rng(42)
    X = rng.exponential(10, size=(n_cells, n_genes)).astype(np.float32)
    X = np.abs(X)  # ensure non-negative
    X[0, 0] = 50.0  # ensure max > 20
    return ad.AnnData(X=X)


# ---------------------------------------------------------------------------
# X state detection
# ---------------------------------------------------------------------------


class TestDetectXState:
    def test_detect_raw_counts(self):
        adata = _make_raw()
        state = inspect_adata(adata)
        assert state.x_state == RAW_COUNTS

    def test_detect_raw_counts_sparse(self):
        adata = _make_raw(sparse=True)
        state = inspect_adata(adata)
        assert state.x_state == RAW_COUNTS

    def test_detect_log_normalized(self):
        adata = _make_log_normalized()
        state = inspect_adata(adata)
        assert state.x_state == LOG_NORMALIZED

    def test_detect_scaled(self):
        adata = _make_scaled()
        state = inspect_adata(adata)
        assert state.x_state == SCALED

    def test_detect_transformed_unknown(self):
        adata = _make_transformed_unknown()
        state = inspect_adata(adata)
        assert state.x_state == TRANSFORMED_UNKNOWN

    def test_log_normalized_heuristic_no_uns(self):
        """Non-integer, non-negative, max < 20, but no uns['log1p']."""
        rng = np.random.default_rng(42)
        X = rng.uniform(0, 10, size=(200, 100)).astype(np.float32)
        adata = ad.AnnData(X=X)
        state = inspect_adata(adata)
        assert state.x_state == LOG_NORMALIZED
        assert any("Inferring log-normalized" in w for w in state.warnings)


# ---------------------------------------------------------------------------
# Raw counts location
# ---------------------------------------------------------------------------


class TestFindRawCounts:
    def test_find_raw_in_layers_counts(self):
        adata = _make_log_normalized()
        rng = np.random.default_rng(42)
        adata.layers["counts"] = rng.poisson(2, size=adata.shape).astype(np.int64)
        state = inspect_adata(adata)
        assert state.raw_counts_location == "layers['counts']"

    def test_find_raw_in_layers_raw_counts(self):
        adata = _make_scaled()
        state = inspect_adata(adata)
        assert state.raw_counts_location == "layers['raw_counts']"

    def test_find_raw_in_x(self):
        adata = _make_raw()
        state = inspect_adata(adata)
        assert state.raw_counts_location == "X"

    def test_find_raw_none(self):
        """Scaled X, no layers → None + warning."""
        rng = np.random.default_rng(42)
        X = rng.normal(0, 1, size=(200, 100)).astype(np.float32)
        adata = ad.AnnData(X=X)
        state = inspect_adata(adata)
        assert state.raw_counts_location is None
        assert any("no raw counts found" in w.lower() for w in state.warnings)

    def test_find_raw_counts_function(self):
        adata = _make_scaled()
        state = inspect_adata(adata)
        raw = find_raw_counts(adata, state)
        assert (raw >= 0).all()
        assert np.allclose(raw[:100].flatten()[:1000],
                           np.round(raw[:100].flatten()[:1000]))

    def test_find_raw_counts_raises_if_none(self):
        rng = np.random.default_rng(42)
        X = rng.normal(0, 1, size=(200, 100)).astype(np.float32)
        adata = ad.AnnData(X=X)
        state = inspect_adata(adata)
        with pytest.raises(ValueError, match="No raw counts found"):
            find_raw_counts(adata, state)

    def test_layers_counts_preferred_over_raw_counts(self):
        """layers['counts'] takes priority over layers['raw_counts']."""
        adata = _make_log_normalized()
        rng = np.random.default_rng(42)
        adata.layers["counts"] = rng.poisson(2, size=adata.shape).astype(np.int64)
        adata.layers["raw_counts"] = rng.poisson(3, size=adata.shape).astype(np.int64)
        state = inspect_adata(adata)
        assert state.raw_counts_location == "layers['counts']"


# ---------------------------------------------------------------------------
# Breadcrumb detection
# ---------------------------------------------------------------------------


class TestBreadcrumbs:
    def test_breadcrumbs_empty(self):
        adata = _make_raw()
        state = inspect_adata(adata)
        assert not state.has_qc_metrics
        assert not state.has_normalized
        assert not state.has_hvg
        assert not state.has_pca
        assert not state.has_neighbors
        assert not state.has_umap
        assert not state.has_clusters
        assert not state.has_cell_types
        assert not state.has_de_results

    def test_breadcrumbs_qc(self):
        adata = _make_raw()
        adata.obs["n_genes_by_counts"] = np.random.default_rng(0).integers(100, 5000, adata.n_obs)
        state = inspect_adata(adata)
        assert state.has_qc_metrics

    def test_breadcrumbs_normalized(self):
        adata = _make_log_normalized()
        state = inspect_adata(adata)
        assert state.has_normalized

    def test_breadcrumbs_pca(self):
        adata = _make_raw()
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 30))
        state = inspect_adata(adata)
        assert state.has_pca

    def test_breadcrumbs_neighbors(self):
        adata = _make_raw()
        adata.uns["neighbors"] = {"params": {}}
        adata.obsp["connectivities"] = csr_matrix((adata.n_obs, adata.n_obs))
        adata.obsp["distances"] = csr_matrix((adata.n_obs, adata.n_obs))
        state = inspect_adata(adata)
        assert state.has_neighbors

    def test_breadcrumbs_neighbors_partial(self):
        """uns['neighbors'] without obsp → not detected."""
        adata = _make_raw()
        adata.uns["neighbors"] = {"params": {}}
        state = inspect_adata(adata)
        assert not state.has_neighbors

    def test_breadcrumbs_umap(self):
        adata = _make_raw()
        adata.obsm["X_umap"] = np.random.default_rng(0).normal(size=(adata.n_obs, 2))
        state = inspect_adata(adata)
        assert state.has_umap

    def test_breadcrumbs_de(self):
        adata = _make_raw()
        adata.uns["rank_genes_groups"] = {"names": {}}
        state = inspect_adata(adata)
        assert state.has_de_results

    def test_breadcrumbs_hvg(self):
        adata = _make_raw()
        adata.var["highly_variable"] = np.random.default_rng(0).choice(
            [True, False], size=adata.n_vars
        )
        state = inspect_adata(adata)
        assert state.has_hvg

    def test_breadcrumbs_velocity(self):
        adata = _make_raw()
        adata.layers["spliced"] = adata.X.copy()
        adata.layers["unspliced"] = adata.X.copy()
        state = inspect_adata(adata)
        assert state.has_velocity_layers

    def test_breadcrumbs_diffmap(self):
        adata = _make_raw()
        adata.obsm["X_diffmap"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))
        state = inspect_adata(adata)
        assert state.has_diffmap


# ---------------------------------------------------------------------------
# Fuzzy column matching
# ---------------------------------------------------------------------------


class TestFuzzyMatching:
    def test_cluster_key_leiden(self):
        adata = _make_raw()
        adata.obs["leiden"] = "0"
        state = inspect_adata(adata)
        assert state.has_clusters
        assert state.cluster_key == "leiden"

    def test_cluster_key_leiden_resolution(self):
        adata = _make_raw()
        adata.obs["leiden_0.8"] = "0"
        state = inspect_adata(adata)
        assert state.has_clusters
        assert state.cluster_key == "leiden_0.8"

    def test_cluster_key_louvain(self):
        adata = _make_raw()
        adata.obs["louvain"] = "0"
        state = inspect_adata(adata)
        assert state.cluster_key == "louvain"

    def test_celltype_key_standard(self):
        adata = _make_raw()
        adata.obs["cell_type"] = "T cell"
        state = inspect_adata(adata)
        assert state.has_cell_types
        assert state.celltype_key == "cell_type"

    def test_celltype_key_celltypist(self):
        adata = _make_raw()
        adata.obs["CellTypist_annotation"] = "T cell"
        state = inspect_adata(adata)
        assert state.celltype_key == "CellTypist_annotation"

    def test_condition_key(self):
        adata = _make_raw()
        rng = np.random.default_rng(42)
        adata.obs["condition"] = rng.choice(["disease", "healthy"], size=adata.n_obs)
        state = inspect_adata(adata)
        assert state.condition_key == "condition"

    def test_condition_key_disease_status(self):
        adata = _make_raw()
        rng = np.random.default_rng(42)
        adata.obs["disease_status"] = rng.choice(["AML", "healthy"], size=adata.n_obs)
        state = inspect_adata(adata)
        assert state.condition_key == "disease_status"

    def test_condition_key_single_value_ignored(self):
        """Condition column with only 1 unique value → not detected."""
        adata = _make_raw()
        adata.obs["condition"] = "healthy"
        state = inspect_adata(adata)
        assert state.condition_key is None

    def test_batch_key(self):
        adata = _make_raw()
        rng = np.random.default_rng(42)
        adata.obs["batch"] = rng.choice(["batch1", "batch2"], size=adata.n_obs)
        state = inspect_adata(adata)
        assert state.batch_key == "batch"

    def test_batch_key_donor(self):
        adata = _make_raw()
        rng = np.random.default_rng(42)
        adata.obs["donor_id"] = rng.choice(["D1", "D2", "D3"], size=adata.n_obs)
        state = inspect_adata(adata)
        assert state.batch_key == "donor_id"


# ---------------------------------------------------------------------------
# Species detection
# ---------------------------------------------------------------------------


class TestSpecies:
    def test_species_human(self):
        """Uppercase gene names → human."""
        genes = [f"GENE{i}" for i in range(100)]
        adata = ad.AnnData(
            X=np.zeros((10, 100)),
            var={"_": genes},
        )
        adata.var_names = genes
        state = inspect_adata(adata)
        assert state.species == "human"

    def test_species_mouse(self):
        """Title-case gene names → mouse."""
        genes = [f"Gene{i}" for i in range(100)]
        adata = ad.AnnData(
            X=np.zeros((10, 100)),
            var={"_": genes},
        )
        adata.var_names = genes
        state = inspect_adata(adata)
        assert state.species == "mouse"


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


class TestSummary:
    def test_summarize_state_basic(self):
        adata = _make_raw()
        state = inspect_adata(adata)
        summary = summarize_state(state)
        assert "Data Inspection" in summary
        assert "200 cells" in summary
        assert "100 genes" in summary
        assert "raw integer counts" in summary

    def test_summarize_state_scaled_with_warnings(self):
        rng = np.random.default_rng(42)
        X = rng.normal(0, 1, size=(200, 100)).astype(np.float32)
        adata = ad.AnnData(X=X)
        state = inspect_adata(adata)
        summary = summarize_state(state)
        assert "z-score scaled" in summary
        assert "⚠️ Warnings" in summary

    def test_summarize_state_with_metadata(self):
        adata = _make_log_normalized()
        rng = np.random.default_rng(42)
        adata.obs["leiden"] = rng.choice(["0", "1", "2"], size=adata.n_obs)
        adata.obs["cell_type"] = rng.choice(["T", "B"], size=adata.n_obs)
        state = inspect_adata(adata)
        summary = summarize_state(state)
        assert "Metadata" in summary
        assert "leiden" in summary
        assert "cell_type" in summary
