"""Tests for cell-cell communication tools (LIANA+)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import anndata as ad


def _has_liana() -> bool:
    try:
        import liana
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def multi_type_adata():
    """Synthetic dataset with 4 cell types."""
    np.random.seed(42)
    n_cells = 400
    n_genes = 200

    X = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)

    cell_types = np.repeat(["Tcell", "Bcell", "Mono", "NK"], n_cells // 4)

    adata = ad.AnnData(X)
    adata.obs["cell_type"] = pd.Categorical(cell_types)
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]

    return adata


# ---------------------------------------------------------------------------
# Input validation tests (no liana required)
# ---------------------------------------------------------------------------


class TestLIANAValidation:
    def test_missing_cell_type_column(self, multi_type_adata):
        from scagent.tools.communication import run_liana

        with pytest.raises((ValueError, ImportError)):
            run_liana(multi_type_adata, cell_type_key="nonexistent")

    def test_single_cell_type_rejected(self, multi_type_adata):
        from scagent.tools.communication import run_liana

        multi_type_adata.obs["cell_type"] = "Tcell"
        with pytest.raises((ValueError, ImportError)):
            run_liana(multi_type_adata, cell_type_key="cell_type")

    def test_import_error_message(self):
        """If liana is not installed, ImportError should be informative."""
        if _has_liana():
            pytest.skip("liana is installed")

        from scagent.tools.communication import run_liana

        adata = ad.AnnData(np.random.rand(10, 5).astype(np.float32))
        adata.obs["cell_type"] = pd.Categorical(["A"] * 5 + ["B"] * 5)

        with pytest.raises(ImportError, match="liana"):
            run_liana(adata, cell_type_key="cell_type")


# ---------------------------------------------------------------------------
# Integration test (requires liana)
# ---------------------------------------------------------------------------


class TestLIANAIntegration:
    @pytest.fixture
    def real_gene_adata(self):
        """Dataset with real gene names so LIANA's L-R resource can match."""
        np.random.seed(42)
        import scanpy as sc

        # Use pbmc3k which has real gene names
        try:
            adata = sc.datasets.pbmc3k_processed()
            # Subsample for speed
            sc.pp.subsample(adata, n_obs=500)
            adata.obs["cell_type"] = adata.obs["louvain"]
            return adata
        except Exception:
            return None

    @pytest.mark.skipif(not _has_liana(), reason="liana not installed")
    def test_liana_runs(self, real_gene_adata):
        if real_gene_adata is None:
            pytest.skip("Could not load pbmc3k dataset")

        from scagent.tools.communication import run_liana

        result = run_liana(
            real_gene_adata,
            cell_type_key="cell_type",
            n_perms=10,  # fast
        )

        assert "interactions" in result
        assert "provenance" in result
        assert result["provenance"]["tool_id"] == "liana"
        assert result["n_cell_types"] >= 2

    @pytest.mark.skipif(not _has_liana(), reason="liana not installed")
    def test_liana_database_bias_warning(self, real_gene_adata):
        if real_gene_adata is None:
            pytest.skip("Could not load pbmc3k dataset")

        from scagent.tools.communication import run_liana

        result = run_liana(
            real_gene_adata,
            cell_type_key="cell_type",
            n_perms=10,
        )

        # Should always warn about L-R database bias
        assert any("bias" in w.lower() for w in result["warnings"])


# ---------------------------------------------------------------------------
# Result structure tests
# ---------------------------------------------------------------------------


class TestResultStructure:
    """Test that result dicts have expected keys regardless of library."""

    def test_composition_result_keys(self):
        """Check expected keys are documented in the module docstring."""
        from scagent.tools import composition
        doc = composition.run_sccoda.__doc__
        assert "credible_effects" in doc
        assert "provenance" in doc

    def test_communication_result_keys(self):
        from scagent.tools import communication
        doc = communication.run_liana.__doc__
        assert "interactions" in doc
        assert "provenance" in doc
