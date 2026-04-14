"""Tests for perturbation screen tools (guide assignment, perturbation DE)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import scanpy as sc

import anndata as ad

from scagent.tools.perturbation import assign_guides, run_perturbation_de


@pytest.fixture
def perturb_adata():
    """Synthetic Perturb-seq dataset with guide assignments."""
    np.random.seed(42)
    n_cells = 500
    n_genes = 200

    X = np.random.poisson(5, (n_cells, n_genes)).astype(np.float32)
    adata = ad.AnnData(X)
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]

    # Assign guides: 100 non-targeting, 3 perturbations of 100 each, 100 unassigned
    guides = (
        ["non-targeting"] * 100
        + ["TP53_guide1"] * 100
        + ["BRCA1_guide1"] * 100
        + ["MYC_guide1"] * 100
        + [""] * 100  # unassigned
    )
    adata.obs["guide_ids"] = guides
    adata.obs["cell_type"] = pd.Categorical(
        np.random.choice(["Tcell", "Mono"], n_cells)
    )

    # Make TP53 perturbation have visible DE (upregulate genes 0-19)
    tp53_mask = np.array(guides) == "TP53_guide1"
    X[tp53_mask, :20] += np.random.poisson(15, (100, 20))

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    return adata


class TestAssignGuides:
    def test_basic_assignment(self, perturb_adata):
        result = assign_guides(perturb_adata, guide_calls_key="guide_ids")

        assert "perturbation" in perturb_adata.obs.columns
        assert "guide" in perturb_adata.obs.columns
        assert result["assignment_stats"]["n_assigned"] == 400
        assert result["assignment_stats"]["n_unassigned"] == 100
        assert result["assignment_stats"]["n_perturbations"] > 0

    def test_target_extraction(self, perturb_adata):
        assign_guides(perturb_adata, guide_calls_key="guide_ids")

        perturbations = perturb_adata.obs["perturbation"].unique()
        assert "TP53" in perturbations
        assert "BRCA1" in perturbations
        assert "MYC" in perturbations

    def test_missing_column_raises(self, perturb_adata):
        with pytest.raises(ValueError, match="not in adata.obs"):
            assign_guides(perturb_adata, guide_calls_key="nonexistent")

    def test_provenance(self, perturb_adata):
        result = assign_guides(perturb_adata, guide_calls_key="guide_ids")
        assert result["provenance"]["tool_id"] == "guide_assignment"

    def test_multi_guide_cells(self):
        adata = ad.AnnData(np.random.rand(10, 5).astype(np.float32))
        adata.obs["guide_ids"] = [
            "TP53_g1,BRCA1_g1",  # multi-guide
            "TP53_g1",
            "",
        ] * 3 + ["MYC_g1"]

        result = assign_guides(adata, guide_calls_key="guide_ids")
        assert result["assignment_stats"]["n_multi_guide"] == 3


class TestPerturbationDE:
    def test_basic_de(self, perturb_adata):
        assign_guides(perturb_adata, guide_calls_key="guide_ids")
        result = run_perturbation_de(
            perturb_adata,
            control_label="non-targeting",
            min_cells=50,
        )

        assert "de_results" in result
        assert "TP53" in result["de_results"]
        assert result["n_perturbations_tested"] == 3
        assert result["provenance"]["tool_id"] == "perturbation_de"

    def test_tp53_has_de_genes(self, perturb_adata):
        assign_guides(perturb_adata, guide_calls_key="guide_ids")
        result = run_perturbation_de(
            perturb_adata,
            control_label="non-targeting",
            min_cells=50,
        )

        # TP53 perturbation should have significant DE genes
        tp53_result = result["de_results"]["TP53"]
        assert tp53_result["n_significant"] > 0

    def test_missing_control_raises(self, perturb_adata):
        assign_guides(perturb_adata, guide_calls_key="guide_ids")
        with pytest.raises(ValueError, match="not found"):
            run_perturbation_de(
                perturb_adata,
                control_label="nonexistent_control",
            )

    def test_min_cells_filter(self, perturb_adata):
        """High min_cells should raise because control itself is too small."""
        assign_guides(perturb_adata, guide_calls_key="guide_ids")
        with pytest.raises(ValueError, match="control cells"):
            run_perturbation_de(
                perturb_adata,
                control_label="non-targeting",
                min_cells=200,  # above control's 100 cells
            )

    def test_missing_perturbation_key_raises(self, perturb_adata):
        with pytest.raises(ValueError, match="not in adata.obs"):
            run_perturbation_de(perturb_adata)
