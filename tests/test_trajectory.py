"""Tests for trajectory inference tools (PAGA, DPT, scVelo)."""

from __future__ import annotations

import numpy as np
import pytest
import scanpy as sc

from scagent.tools.trajectory import (
    run_diffusion_pseudotime,
    run_paga,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def branching_adata():
    """Synthetic dataset with a branching trajectory (3 branches from a root)."""
    np.random.seed(42)
    n_per_branch = 100
    n_genes = 200

    # Root cells (progenitors)
    root = np.random.poisson(5, (n_per_branch, n_genes)).astype(np.float32)
    # Branch A — upregulate genes 0-49
    branch_a = root.copy()
    branch_a[:, :50] += np.random.poisson(8, (n_per_branch, 50))
    # Branch B — upregulate genes 50-99
    branch_b = root.copy()
    branch_b[:, 50:100] += np.random.poisson(8, (n_per_branch, 50))
    # Branch C — upregulate genes 100-149
    branch_c = root.copy()
    branch_c[:, 100:150] += np.random.poisson(8, (n_per_branch, 50))

    X = np.vstack([root, branch_a, branch_b, branch_c])
    labels = (
        ["root"] * n_per_branch
        + ["branch_A"] * n_per_branch
        + ["branch_B"] * n_per_branch
        + ["branch_C"] * n_per_branch
    )

    import anndata as ad
    import pandas as pd

    adata = ad.AnnData(X)
    adata.obs["cell_type"] = pd.Categorical(labels)
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]

    # Standard preprocessing
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=100)
    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_pcs=20)
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata)

    return adata


@pytest.fixture
def simple_adata():
    """Minimal preprocessed dataset for quick tests."""
    np.random.seed(0)
    import anndata as ad
    import pandas as pd

    adata = ad.AnnData(np.random.poisson(3, (200, 100)).astype(np.float32))
    adata.obs["cell_type"] = pd.Categorical(
        ["type_A"] * 100 + ["type_B"] * 100
    )
    adata.var_names = [f"g{i}" for i in range(100)]

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=20)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    return adata


# ---------------------------------------------------------------------------
# PAGA tests
# ---------------------------------------------------------------------------


class TestRunPAGA:
    def test_basic_paga(self, branching_adata):
        result = run_paga(branching_adata, groups="leiden")

        assert "connectivities" in result
        assert "provenance" in result
        assert result["provenance"]["tool_id"] == "paga"
        assert result["n_groups"] >= 2
        assert "paga" in branching_adata.uns

    def test_paga_with_cell_type_groups(self, branching_adata):
        result = run_paga(branching_adata, groups="cell_type")

        assert result["n_groups"] == 4  # root + 3 branches
        assert set(result["group_names"]) == {"root", "branch_A", "branch_B", "branch_C"}

    def test_paga_connectivity_check(self, branching_adata):
        result = run_paga(branching_adata, groups="cell_type")
        # Should be connected for this synthetic data
        assert result["n_components"] >= 1
        # Should have edges
        assert len(result["connectivities"]) > 0

    def test_paga_requires_neighbors(self):
        import anndata as ad
        adata = ad.AnnData(np.random.rand(50, 10).astype(np.float32))
        adata.obs["leiden"] = "0"

        with pytest.raises(ValueError, match="Neighbor graph not found"):
            run_paga(adata)

    def test_paga_requires_valid_groups_key(self, simple_adata):
        with pytest.raises(ValueError, match="not in adata.obs"):
            run_paga(simple_adata, groups="nonexistent_key")

    def test_paga_summary_format(self, branching_adata):
        result = run_paga(branching_adata, groups="leiden")
        assert isinstance(result["summary"], str)
        assert len(result["summary"]) > 0

    def test_paga_provenance(self, branching_adata):
        result = run_paga(branching_adata, groups="cell_type", model="v1.2")
        prov = result["provenance"]
        assert prov["parameters"]["groups"] == "cell_type"
        assert prov["parameters"]["model"] == "v1.2"
        assert prov["n_groups"] == 4

    def test_paga_plot_dir(self, branching_adata, tmp_path):
        result = run_paga(
            branching_adata, groups="leiden",
            plot_dir=str(tmp_path / "plots"),
        )
        assert len(result["plots"]) >= 1
        import os
        assert os.path.exists(result["plots"][0])


# ---------------------------------------------------------------------------
# DPT tests
# ---------------------------------------------------------------------------


class TestRunDPT:
    def test_basic_dpt(self, branching_adata):
        result = run_diffusion_pseudotime(
            branching_adata, root_cell_type="root",
            cell_type_key="cell_type",
        )

        assert "dpt_pseudotime" in branching_adata.obs.columns
        assert "pseudotime_stats" in result
        assert "root" in result["pseudotime_stats"]
        assert result["provenance"]["tool_id"] == "diffusion_pseudotime"

    def test_dpt_root_has_lowest_pseudotime(self, branching_adata):
        result = run_diffusion_pseudotime(
            branching_adata, root_cell_type="root",
            cell_type_key="cell_type",
        )

        # Root should have lowest median pseudotime
        stats = result["pseudotime_stats"]
        root_median = stats["root"]["median"]
        for ct, ct_stats in stats.items():
            if ct != "root":
                assert ct_stats["median"] >= root_median - 0.05, (
                    f"Branch '{ct}' has lower pseudotime than root "
                    f"({ct_stats['median']:.3f} vs {root_median:.3f})"
                )

    def test_dpt_cell_ordering(self, branching_adata):
        result = run_diffusion_pseudotime(
            branching_adata, root_cell_type="root",
            cell_type_key="cell_type",
        )

        # Root should be first in the ordering
        ordering = result["cell_ordering"]
        assert ordering[0] == "root", f"Expected root first, got {ordering}"

    def test_dpt_explicit_root_index(self, branching_adata):
        result = run_diffusion_pseudotime(
            branching_adata, root_cell_index=0,
            cell_type_key="cell_type",
        )
        assert "explicit index 0" in result["root_info"]

    def test_dpt_requires_root(self, branching_adata):
        with pytest.raises(ValueError, match="Must provide either"):
            run_diffusion_pseudotime(branching_adata)

    def test_dpt_invalid_root_type(self, branching_adata):
        with pytest.raises(ValueError, match="not found"):
            run_diffusion_pseudotime(
                branching_adata, root_cell_type="nonexistent",
                cell_type_key="cell_type",
            )

    def test_dpt_requires_neighbors(self):
        import anndata as ad
        adata = ad.AnnData(np.random.rand(50, 10).astype(np.float32))

        with pytest.raises(ValueError, match="Neighbor graph not found"):
            run_diffusion_pseudotime(adata, root_cell_index=0)

    def test_dpt_summary(self, branching_adata):
        result = run_diffusion_pseudotime(
            branching_adata, root_cell_type="root",
            cell_type_key="cell_type",
        )
        assert "DPT computed" in result["summary"]
        assert "→" in result["summary"]

    def test_dpt_plot_dir(self, branching_adata, tmp_path):
        # Need UMAP for plotting
        result = run_diffusion_pseudotime(
            branching_adata, root_cell_type="root",
            cell_type_key="cell_type",
            plot_dir=str(tmp_path / "plots"),
        )
        assert len(result["plots"]) >= 1


# ---------------------------------------------------------------------------
# scVelo tests (layer validation only — no actual scvelo needed)
# ---------------------------------------------------------------------------


class TestRunPalantir:
    def test_palantir_basic(self, branching_adata):
        from scagent.tools.trajectory import run_palantir

        result = run_palantir(
            branching_adata, root_cell_type="root",
            cell_type_key="cell_type", num_waypoints=50,
        )

        assert "palantir_pseudotime" in branching_adata.obs.columns
        assert result["provenance"]["tool_id"] == "palantir"
        assert "pseudotime_stats" in result
        assert "root" in result["pseudotime_stats"]

    def test_palantir_root_has_lowest_pseudotime(self, branching_adata):
        from scagent.tools.trajectory import run_palantir

        result = run_palantir(
            branching_adata, root_cell_type="root",
            cell_type_key="cell_type", num_waypoints=50,
        )

        stats = result["pseudotime_stats"]
        root_median = stats["root"]["median"]
        for ct, ct_stats in stats.items():
            if ct != "root":
                assert ct_stats["median"] >= root_median - 0.1, (
                    f"'{ct}' has lower pseudotime than root"
                )

    def test_palantir_requires_root(self, branching_adata):
        from scagent.tools.trajectory import run_palantir

        with pytest.raises(ValueError, match="Must provide either"):
            run_palantir(branching_adata)

    def test_palantir_requires_pca(self):
        import anndata as ad
        adata = ad.AnnData(np.random.rand(50, 10).astype(np.float32))
        from scagent.tools.trajectory import run_palantir

        with pytest.raises(ValueError, match="PCA"):
            run_palantir(adata, root_cell_index=0)


class TestRunCellRankValidation:
    def test_cellrank_requires_velocity(self, simple_adata):
        from scagent.tools.trajectory import run_cellrank

        with pytest.raises(ValueError, match="velocity"):
            run_cellrank(simple_adata)


class TestRunScVeloValidation:
    def test_scvelo_refuses_without_spliced(self, simple_adata):
        from scagent.tools.trajectory import run_scvelo

        with pytest.raises(ValueError, match="spliced"):
            run_scvelo(simple_adata)

    def test_scvelo_refuses_without_unspliced(self, simple_adata):
        from scagent.tools.trajectory import run_scvelo

        simple_adata.layers["spliced"] = simple_adata.X.copy()
        with pytest.raises(ValueError, match="unspliced"):
            run_scvelo(simple_adata)


# ---------------------------------------------------------------------------
# DAG integration
# ---------------------------------------------------------------------------


class TestTrajectoryDAG:
    def test_trajectory_dag_has_steps(self):
        from scagent.dag import AnalysisDAG, _developmental_trajectory_steps
        from unittest.mock import MagicMock

        ctx = MagicMock()
        ctx.paradigm = "developmental_trajectory"
        ctx.needs_batch_correction.return_value = False

        steps = _developmental_trajectory_steps(ctx)
        step_ids = [s.id for s in steps]

        assert "paga" in step_ids
        assert "pseudotime" in step_ids
        assert "rna_velocity" in step_ids
        assert "fate_mapping" in step_ids

    def test_trajectory_dag_ordering(self):
        from scagent.dag import _developmental_trajectory_steps
        from unittest.mock import MagicMock

        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        steps = _developmental_trajectory_steps(ctx)
        step_ids = [s.id for s in steps]

        # PAGA must come after clustering
        assert step_ids.index("paga") > step_ids.index("clustering")
        # Pseudotime must come after PAGA and annotation
        assert step_ids.index("pseudotime") > step_ids.index("paga")
        assert step_ids.index("pseudotime") > step_ids.index("annotation")
        # Fate mapping must come after velocity
        assert step_ids.index("fate_mapping") > step_ids.index("rna_velocity")

    def test_trajectory_dag_uses_palantir(self):
        """DAG should use Palantir for pseudotime, not DPT. [BP-2 Ch. 14]"""
        from scagent.dag import _developmental_trajectory_steps
        from unittest.mock import MagicMock

        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        steps = _developmental_trajectory_steps(ctx)
        tool_map = {s.id: s.tool_id for s in steps}

        assert tool_map["paga"] == "paga"
        assert tool_map["pseudotime"] == "palantir"
        assert tool_map["rna_velocity"] == "scvelo_velocity"
        assert tool_map["fate_mapping"] == "cellrank"

    def test_velocity_and_fate_mapping_optional(self):
        from scagent.dag import _developmental_trajectory_steps
        from unittest.mock import MagicMock

        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        steps = _developmental_trajectory_steps(ctx)
        vel_step = next(s for s in steps if s.id == "rna_velocity")
        fate_step = next(s for s in steps if s.id == "fate_mapping")

        assert vel_step.required is False
        assert vel_step.conditional is True
        assert fate_step.required is False
        assert fate_step.conditional is True
