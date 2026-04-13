#!/usr/bin/env python3
"""Unit tests for scagent.tools.pseudobulk_de."""

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy.sparse import csr_matrix

from scagent.tools.pseudobulk_de import aggregate_pseudobulk, run_pseudobulk_de


# ---------------------------------------------------------------------------
# Fixtures: synthetic multi-sample, multi-condition data
# ---------------------------------------------------------------------------

def _make_multi_sample_adata(
    n_cells_per_sample: int = 200,
    n_genes: int = 500,
    n_samples_per_condition: int = 3,
    n_cell_types: int = 2,
    de_effect_size: float = 2.0,
    n_de_genes: int = 50,
    seed: int = 42,
) -> ad.AnnData:
    """Create a synthetic multi-sample AnnData with known DE signal.

    - Two conditions: 'disease' and 'healthy'
    - n_samples_per_condition samples per condition
    - n_cell_types cell types
    - First n_de_genes genes are upregulated in 'disease' for cell type 0
    """
    rng = np.random.default_rng(seed)
    total_samples = n_samples_per_condition * 2
    total_cells = n_cells_per_sample * total_samples * n_cell_types

    conditions = ["disease"] * n_samples_per_condition + ["healthy"] * n_samples_per_condition
    sample_ids = [f"sample_{i}" for i in range(total_samples)]

    obs_rows = []
    X_blocks = []

    for ct_idx in range(n_cell_types):
        ct_name = f"CellType_{ct_idx}"
        for s_idx, (sid, cond) in enumerate(zip(sample_ids, conditions)):
            # Base expression: Poisson
            base_mean = rng.uniform(5, 20, size=n_genes)
            counts = rng.poisson(base_mean, size=(n_cells_per_sample, n_genes))

            # Inject DE signal: cell type 0 + disease condition
            if ct_idx == 0 and cond == "disease":
                counts[:, :n_de_genes] = rng.poisson(
                    base_mean[:n_de_genes] * de_effect_size,
                    size=(n_cells_per_sample, n_de_genes),
                )

            X_blocks.append(counts)
            for _ in range(n_cells_per_sample):
                obs_rows.append({
                    "cell_type": ct_name,
                    "sample": sid,
                    "condition": cond,
                })

    X = np.vstack(X_blocks).astype(np.float32)
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes)])

    adata = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    # Ensure categoricals
    for col in ("cell_type", "sample", "condition"):
        adata.obs[col] = adata.obs[col].astype("category")
    return adata


@pytest.fixture
def multi_sample_adata():
    return _make_multi_sample_adata()


@pytest.fixture
def small_adata():
    """Tiny dataset for fast tests."""
    return _make_multi_sample_adata(
        n_cells_per_sample=30, n_genes=100, n_samples_per_condition=2,
        n_cell_types=1, n_de_genes=20,
    )


# ---------------------------------------------------------------------------
# Tests: aggregation
# ---------------------------------------------------------------------------

class TestAggregatePseudobulk:
    def test_returns_all_cell_types(self, multi_sample_adata):
        pb = aggregate_pseudobulk(
            multi_sample_adata, "cell_type", "sample", "condition"
        )
        assert "CellType_0" in pb
        assert "CellType_1" in pb

    def test_counts_shape(self, multi_sample_adata):
        pb = aggregate_pseudobulk(
            multi_sample_adata, "cell_type", "sample", "condition"
        )
        ct0 = pb["CellType_0"]
        # 6 samples, 500 genes
        assert ct0["counts"].shape == (500, 6)

    def test_counts_are_integers(self, multi_sample_adata):
        pb = aggregate_pseudobulk(
            multi_sample_adata, "cell_type", "sample", "condition"
        )
        counts = pb["CellType_0"]["counts"]
        assert (counts == counts.astype(int)).all().all()

    def test_metadata_has_conditions(self, multi_sample_adata):
        pb = aggregate_pseudobulk(
            multi_sample_adata, "cell_type", "sample", "condition"
        )
        meta = pb["CellType_0"]["metadata"]
        assert "condition" in meta.columns
        assert set(meta["condition"].unique()) == {"disease", "healthy"}

    def test_min_cells_filter(self, multi_sample_adata):
        # Set min_cells very high — should drop everything
        pb = aggregate_pseudobulk(
            multi_sample_adata, "cell_type", "sample", "condition",
            min_cells=10000,
        )
        assert len(pb) == 0

    def test_missing_column_raises(self, multi_sample_adata):
        with pytest.raises(ValueError, match="not in adata.obs"):
            aggregate_pseudobulk(
                multi_sample_adata, "nonexistent", "sample", "condition"
            )


# ---------------------------------------------------------------------------
# Tests: DE
# ---------------------------------------------------------------------------

class TestRunPseudobulkDE:
    def test_detects_de_genes(self, multi_sample_adata):
        results = run_pseudobulk_de(
            multi_sample_adata,
            cell_type_key="cell_type",
            sample_key="sample",
            condition_key="condition",
        )
        # CellType_0 should have DE genes (we injected signal)
        assert "CellType_0" in results["de_dataframes"]
        df = results["de_dataframes"]["CellType_0"]
        n_sig = df["significant"].sum()
        assert n_sig > 0, "Expected some DE genes in CellType_0"

    def test_cell_type_without_signal_has_fewer_de(self, multi_sample_adata):
        results = run_pseudobulk_de(
            multi_sample_adata,
            cell_type_key="cell_type",
            sample_key="sample",
            condition_key="condition",
        )
        # CellType_1 has no injected signal
        if "CellType_1" in results["de_dataframes"]:
            df1 = results["de_dataframes"]["CellType_1"]
            n_sig_1 = df1["significant"].sum()
            df0 = results["de_dataframes"]["CellType_0"]
            n_sig_0 = df0["significant"].sum()
            assert n_sig_0 > n_sig_1, (
                f"CellType_0 ({n_sig_0} DE) should have more DE than CellType_1 ({n_sig_1})"
            )

    def test_summary_format(self, small_adata):
        results = run_pseudobulk_de(
            small_adata,
            cell_type_key="cell_type",
            sample_key="sample",
            condition_key="condition",
        )
        assert len(results["summary"]) > 0
        row = results["summary"][0]
        assert "cell_type" in row
        assert "n_significant" in row
        assert "n_up" in row
        assert "n_down" in row

    def test_provenance_recorded(self, small_adata):
        results = run_pseudobulk_de(
            small_adata,
            cell_type_key="cell_type",
            sample_key="sample",
            condition_key="condition",
        )
        prov = results["provenance"]
        assert prov["tool_id"] == "deseq2_pseudobulk"
        assert prov["parameters"]["cell_type_key"] == "cell_type"
        assert prov["n_cell_types_tested"] > 0

    def test_too_few_conditions_raises(self):
        # Only one condition
        adata = _make_multi_sample_adata(n_samples_per_condition=2)
        adata.obs["condition"] = "healthy"  # overwrite to single condition
        with pytest.raises(ValueError, match="at least 2 conditions"):
            run_pseudobulk_de(
                adata,
                cell_type_key="cell_type",
                sample_key="sample",
                condition_key="condition",
            )

    def test_too_few_replicates_raises(self):
        # Only 1 sample per condition
        adata = _make_multi_sample_adata(n_samples_per_condition=1)
        with pytest.raises(ValueError, match="biological replicates"):
            run_pseudobulk_de(
                adata,
                cell_type_key="cell_type",
                sample_key="sample",
                condition_key="condition",
                min_replicates=2,
            )

    def test_volcano_plots(self, small_adata, tmp_path):
        results = run_pseudobulk_de(
            small_adata,
            cell_type_key="cell_type",
            sample_key="sample",
            condition_key="condition",
            plot_dir=str(tmp_path),
        )
        assert len(results["plots"]) > 0
        from pathlib import Path
        for p in results["plots"]:
            assert Path(p).exists()
