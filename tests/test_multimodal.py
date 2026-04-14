"""Tests for multimodal (CITE-seq) tools."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import anndata as ad

from scagent.tools.multimodal import load_protein, normalize_protein, run_wnn


@pytest.fixture
def citeseq_adata():
    """Synthetic CITE-seq dataset with ADT features in feature_types."""
    np.random.seed(42)
    n_cells = 200
    n_rna_genes = 100
    n_proteins = 20

    # RNA + protein as single matrix (Cell Ranger multi format)
    rna_X = np.random.poisson(5, (n_cells, n_rna_genes)).astype(np.float32)
    prot_X = np.random.poisson(50, (n_cells, n_proteins)).astype(np.float32)
    X = np.hstack([rna_X, prot_X])

    adata = ad.AnnData(X)
    feature_types = (
        ["Gene Expression"] * n_rna_genes
        + ["Antibody Capture"] * n_proteins
    )
    adata.var["feature_types"] = feature_types
    adata.var_names = (
        [f"gene_{i}" for i in range(n_rna_genes)]
        + [f"CD{i}" for i in range(n_proteins)]
    )
    adata.obs["cell_type"] = pd.Categorical(
        np.random.choice(["Tcell", "Bcell"], n_cells)
    )

    return adata


@pytest.fixture
def adata_with_protein(citeseq_adata):
    """CITE-seq adata with protein already loaded + PCA computed."""
    result = load_protein(citeseq_adata)
    # Run PCA on RNA genes only
    rna_mask = citeseq_adata.var["feature_types"] == "Gene Expression"
    rna_adata = citeseq_adata[:, rna_mask].copy()
    sc.pp.normalize_total(rna_adata)
    sc.pp.log1p(rna_adata)
    sc.tl.pca(rna_adata, n_comps=20)
    citeseq_adata.obsm["X_pca"] = rna_adata.obsm["X_pca"]
    return citeseq_adata


class TestLoadProtein:
    def test_load_from_feature_types(self, citeseq_adata):
        result = load_protein(citeseq_adata)

        assert "protein_counts" in citeseq_adata.obsm
        assert result["n_proteins"] == 20
        assert citeseq_adata.obsm["protein_counts"].shape == (200, 20)
        assert result["provenance"]["tool_id"] == "load_protein"

    def test_protein_names_stored(self, citeseq_adata):
        load_protein(citeseq_adata)
        assert "protein_names" in citeseq_adata.uns
        assert len(citeseq_adata.uns["protein_names"]) == 20

    def test_no_protein_data_raises(self):
        adata = ad.AnnData(np.random.rand(10, 5).astype(np.float32))
        adata.var_names = [f"gene_{i}" for i in range(5)]
        with pytest.raises(ValueError, match="No protein"):
            load_protein(adata)


class TestNormalizeProtein:
    def test_clr_normalization(self, citeseq_adata):
        load_protein(citeseq_adata)
        result = normalize_protein(citeseq_adata, method="clr")

        assert "protein_normalized" in citeseq_adata.obsm
        assert citeseq_adata.obsm["protein_normalized"].shape == (200, 20)
        assert result["method"] == "clr"

    def test_clr_values_centered(self, citeseq_adata):
        load_protein(citeseq_adata)
        normalize_protein(citeseq_adata)

        # CLR should be approximately centered per cell (mean ~0)
        means = citeseq_adata.obsm["protein_normalized"].mean(axis=1)
        assert np.allclose(means, 0, atol=0.01)

    def test_requires_protein_counts(self):
        adata = ad.AnnData(np.random.rand(10, 5).astype(np.float32))
        with pytest.raises(ValueError, match="protein counts"):
            normalize_protein(adata)

    def test_isotype_control_warning(self, citeseq_adata):
        load_protein(citeseq_adata)
        result = normalize_protein(citeseq_adata)

        # No isotype controls in our synthetic data
        assert any("isotype" in w.lower() for w in result["warnings"])


class TestWNN:
    def test_wnn_basic(self, adata_with_protein):
        normalize_protein(adata_with_protein)
        result = run_wnn(adata_with_protein)

        assert result["provenance"]["tool_id"] == "wnn"
        assert result["rna_weight"] == 0.5

    def test_wnn_requires_pca(self, citeseq_adata):
        load_protein(citeseq_adata)
        normalize_protein(citeseq_adata)

        with pytest.raises(ValueError, match="PCA"):
            run_wnn(citeseq_adata)

    def test_wnn_requires_protein(self):
        adata = ad.AnnData(np.random.rand(10, 5).astype(np.float32))
        adata.obsm["X_pca"] = np.random.rand(10, 5)

        with pytest.raises(ValueError, match="protein"):
            run_wnn(adata)


class TestMultimodalDAG:
    def test_multimodal_dag_structure(self):
        from unittest.mock import MagicMock
        from scagent.dag import _multimodal_steps

        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        steps = _multimodal_steps(ctx)
        step_ids = [s.id for s in steps]

        assert "load_protein" in step_ids
        assert "normalize_protein" in step_ids
        assert "wnn" in step_ids
        assert "protein_markers" in step_ids

    def test_wnn_depends_on_pca_and_protein(self):
        from unittest.mock import MagicMock
        from scagent.dag import _multimodal_steps

        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        steps = _multimodal_steps(ctx)
        wnn = next(s for s in steps if s.id == "wnn")

        assert "pca" in wnn.depends_on
        assert "normalize_protein" in wnn.depends_on

    def test_clustering_on_wnn(self):
        from unittest.mock import MagicMock
        from scagent.dag import _multimodal_steps

        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        steps = _multimodal_steps(ctx)
        cluster = next(s for s in steps if s.id == "clustering")

        assert "wnn" in cluster.depends_on
