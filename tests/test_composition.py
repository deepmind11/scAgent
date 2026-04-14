"""Tests for composition analysis tools (scCODA, Milo)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import anndata as ad


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _has_pertpy() -> bool:
    try:
        import pertpy
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def multi_condition_adata():
    """Synthetic multi-condition dataset with 6 donors (3 per condition)."""
    np.random.seed(42)
    n_cells = 600
    n_genes = 100

    X = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)

    # 3 cell types, 6 donors, 2 conditions
    cell_types = np.repeat(["TypeA", "TypeB", "TypeC"], n_cells // 3)
    donors = np.tile(np.repeat(["D1", "D2", "D3", "D4", "D5", "D6"], n_cells // 6), 1)
    # Truncate to match
    cell_types = cell_types[:n_cells]
    donors = donors[:n_cells]
    conditions = np.where(np.isin(donors, ["D1", "D2", "D3"]), "ctrl", "disease")

    adata = ad.AnnData(X)
    adata.obs["cell_type"] = pd.Categorical(cell_types)
    adata.obs["donor"] = pd.Categorical(donors)
    adata.obs["condition"] = pd.Categorical(conditions)
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]

    return adata


@pytest.fixture
def preprocessed_adata(multi_condition_adata):
    """Multi-condition data with neighbors computed."""
    import scanpy as sc

    adata = multi_condition_adata
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=15)
    return adata


# ---------------------------------------------------------------------------
# scCODA input validation tests (no pertpy required)
# ---------------------------------------------------------------------------


class TestScCODAValidation:
    def test_missing_condition_column(self, multi_condition_adata):
        from scagent.tools.composition import run_sccoda

        with pytest.raises((ValueError, ImportError)):
            run_sccoda(
                multi_condition_adata,
                condition_key="nonexistent",
                sample_key="donor",
                cell_type_key="cell_type",
            )

    def test_missing_sample_column(self, multi_condition_adata):
        from scagent.tools.composition import run_sccoda

        with pytest.raises((ValueError, ImportError)):
            run_sccoda(
                multi_condition_adata,
                condition_key="condition",
                sample_key="nonexistent",
                cell_type_key="cell_type",
            )

    def test_single_condition_rejected(self, multi_condition_adata):
        from scagent.tools.composition import run_sccoda

        multi_condition_adata.obs["condition"] = "ctrl"
        with pytest.raises((ValueError, ImportError)):
            run_sccoda(
                multi_condition_adata,
                condition_key="condition",
                sample_key="donor",
                cell_type_key="cell_type",
            )

    @pytest.mark.skipif(
        not _has_pertpy(), reason="pertpy not installed"
    )
    def test_sccoda_runs(self):
        """scCODA needs a proportion shift so auto-reference can find a stable type."""
        from scagent.tools.composition import run_sccoda

        np.random.seed(42)
        n_cells = 600
        n_genes = 50
        X = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)

        # Create imbalanced proportions: ctrl has equal, disease has more TypeA
        cell_types = (
            ["TypeA"] * 100 + ["TypeB"] * 100 + ["TypeC"] * 100  # ctrl: equal
            + ["TypeA"] * 200 + ["TypeB"] * 50 + ["TypeC"] * 50  # disease: TypeA enriched
        )
        donors = (
            ["D1"] * 50 + ["D2"] * 50 + ["D1"] * 50 + ["D2"] * 50 + ["D1"] * 50 + ["D2"] * 50
            + ["D3"] * 100 + ["D4"] * 100 + ["D3"] * 25 + ["D4"] * 25 + ["D3"] * 25 + ["D4"] * 25
        )
        conditions = ["ctrl"] * 300 + ["disease"] * 300

        adata = ad.AnnData(X)
        adata.obs["cell_type"] = pd.Categorical(cell_types)
        adata.obs["donor"] = pd.Categorical(donors)
        adata.obs["condition"] = pd.Categorical(conditions)
        adata.var_names = [f"gene_{i}" for i in range(n_genes)]

        result = run_sccoda(
            adata,
            condition_key="condition",
            sample_key="donor",
            cell_type_key="cell_type",
            reference_cell_type="TypeC",  # explicitly set a stable reference
        )
        assert "credible_effects" in result
        assert "provenance" in result
        assert result["provenance"]["tool_id"] == "sccoda"


# ---------------------------------------------------------------------------
# Milo input validation tests
# ---------------------------------------------------------------------------


class TestMiloValidation:
    def test_milo_requires_neighbors(self, multi_condition_adata):
        from scagent.tools.composition import run_milo

        with pytest.raises((ValueError, ImportError)):
            run_milo(
                multi_condition_adata,
                condition_key="condition",
                sample_key="donor",
            )

    def test_milo_missing_condition_column(self, preprocessed_adata):
        from scagent.tools.composition import run_milo

        with pytest.raises((ValueError, ImportError)):
            run_milo(
                preprocessed_adata,
                condition_key="nonexistent",
                sample_key="donor",
            )

    def test_milo_neighbor_power_warning(self, preprocessed_adata):
        """n_neighbors=15 < 3*6=18 should trigger warning."""
        from scagent.tools.composition import run_milo

        # This will either warn (if pertpy installed) or raise ImportError
        try:
            result = run_milo(
                preprocessed_adata,
                condition_key="condition",
                sample_key="donor",
            )
            # If pertpy is available, check warning
            assert any("n_neighbors" in w for w in result["warnings"])
        except ImportError:
            pytest.skip("pertpy not installed")



