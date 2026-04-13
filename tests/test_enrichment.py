#!/usr/bin/env python3
"""Unit tests for scagent.tools.enrichment."""

import numpy as np
import pandas as pd
import pytest

from scagent.tools.enrichment import run_gsea, _build_ranking


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def de_dataframe():
    """Synthetic DE results with known structure."""
    rng = np.random.default_rng(42)
    n_genes = 300
    genes = [f"GENE{i}" for i in range(n_genes)]
    df = pd.DataFrame({
        "log2FoldChange": rng.normal(0, 1.5, n_genes),
        "pvalue": rng.uniform(0, 1, n_genes),
        "padj": rng.uniform(0, 1, n_genes),
    }, index=genes)
    # Inject a few strong signals
    df.loc["GENE0":"GENE9", "log2FoldChange"] = rng.uniform(2, 4, 10)
    df.loc["GENE0":"GENE9", "pvalue"] = rng.uniform(1e-10, 1e-5, 10)
    return df


@pytest.fixture
def de_dict(de_dataframe):
    """DE results as a dict keyed by cell type."""
    return {"Monocytes": de_dataframe}


# ---------------------------------------------------------------------------
# Tests: ranking
# ---------------------------------------------------------------------------

class TestBuildRanking:
    def test_log2fc_ranking(self, de_dataframe):
        rnk = _build_ranking(de_dataframe, "log2fc")
        assert rnk is not None
        assert len(rnk) > 0
        # Should be sorted descending
        assert rnk.iloc[0] >= rnk.iloc[-1]

    def test_signal_to_noise_ranking(self, de_dataframe):
        rnk = _build_ranking(de_dataframe, "signal_to_noise")
        assert rnk is not None
        assert len(rnk) > 0

    def test_missing_column_returns_none(self):
        df = pd.DataFrame({"some_col": [1, 2, 3]})
        rnk = _build_ranking(df, "log2fc")
        assert rnk is None


# ---------------------------------------------------------------------------
# Tests: GSEA
# ---------------------------------------------------------------------------

class TestRunGSEA:
    def test_runs_on_single_df(self, de_dataframe):
        """GSEA should accept a single DataFrame (not just a dict)."""
        result = run_gsea(
            de_dataframe,
            gene_sets="KEGG_2021_Human",
            permutation_num=100,  # fast for tests
        )
        assert "enrichment_results" in result
        assert "summary" in result
        assert "provenance" in result
        assert result["provenance"]["tool_id"] == "gsea"

    def test_runs_on_dict(self, de_dict):
        result = run_gsea(
            de_dict,
            gene_sets="KEGG_2021_Human",
            permutation_num=100,
        )
        assert "enrichment_results" in result
        if "Monocytes" in result["enrichment_results"]:
            assert len(result["enrichment_results"]["Monocytes"]) > 0

    def test_provenance_recorded(self, de_dict):
        result = run_gsea(
            de_dict,
            gene_sets="KEGG_2021_Human",
            permutation_num=100,
        )
        prov = result["provenance"]
        assert prov["parameters"]["gene_sets"] == "KEGG_2021_Human"
        assert prov["parameters"]["permutation_num"] == 100

    def test_too_few_genes_warns(self):
        """DataFrame with very few genes should produce a warning."""
        df = pd.DataFrame({
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.01, 0.5],
        }, index=["A", "B"])
        result = run_gsea(df, gene_sets="KEGG_2021_Human", permutation_num=100)
        assert len(result["warnings"]) > 0
