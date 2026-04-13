#!/usr/bin/env python3
"""Unit tests for scagent.context."""

import json

import pytest

from scagent.context import ExperimentContext, VALID_PARADIGMS


@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path / ".scagent"


@pytest.fixture
def ctx(tmp_dir):
    return ExperimentContext(tmp_dir)


class TestCreateContext:
    def test_minimal_valid(self, ctx):
        ctx.paradigm = "cell_atlas"
        ctx.organism = {"species": "Homo sapiens", "ncbi_taxon": 9606}
        ctx.tissue = {"name": "PBMCs"}

        errors = ctx.validate()
        assert errors == []
        assert ctx.is_complete()

    def test_empty_context_invalid(self, ctx):
        errors = ctx.validate()
        assert len(errors) >= 2  # species and tissue missing (paradigm is optional)

    def test_invalid_paradigm_raises(self, ctx):
        with pytest.raises(ValueError, match="Invalid paradigm"):
            ctx.paradigm = "nonexistent"


class TestValidation:
    def test_missing_paradigm_is_valid(self, ctx):
        """Paradigm is optional — missing paradigm is not an error."""
        ctx.organism = {"species": "Homo sapiens"}
        ctx.tissue = {"name": "brain"}
        errors = ctx.validate()
        assert not any("paradigm" in e for e in errors)

    def test_missing_tissue(self, ctx):
        ctx.paradigm = "cell_atlas"
        ctx.organism = {"species": "Homo sapiens"}
        errors = ctx.validate()
        assert any("tissue" in e for e in errors)

    def test_missing_species(self, ctx):
        ctx.paradigm = "cell_atlas"
        ctx.tissue = {"name": "brain"}
        errors = ctx.validate()
        assert any("species" in e for e in errors)

    def test_conditions_required_for_disease_vs_healthy(self, ctx):
        ctx.paradigm = "disease_vs_healthy"
        ctx.organism = {"species": "Homo sapiens"}
        ctx.tissue = {"name": "PBMCs"}
        # No conditions set
        errors = ctx.validate()
        assert any("conditions" in e for e in errors)

    def test_conditions_not_required_for_cell_atlas(self, ctx):
        ctx.paradigm = "cell_atlas"
        ctx.organism = {"species": "Homo sapiens"}
        ctx.tissue = {"name": "PBMCs"}
        # No conditions — should be fine for atlas
        errors = ctx.validate()
        assert not any("conditions" in e for e in errors)

    def test_conditions_valid_when_set(self, ctx):
        ctx.paradigm = "disease_vs_healthy"
        ctx.organism = {"species": "Homo sapiens"}
        ctx.tissue = {"name": "PBMCs"}
        ctx.design = {"conditions": ["disease", "healthy"]}
        errors = ctx.validate()
        assert errors == []


class TestInferFromData:
    def test_infer_human(self):
        import anndata as ad
        import numpy as np
        from scipy.sparse import csr_matrix

        # Human genes: uppercase
        var_names = [f"GENE{i}" for i in range(50)]
        var_names[0] = "MT-CO1"  # mitochondrial
        X = csr_matrix(np.ones((10, 50), dtype=np.float32))
        adata = ad.AnnData(X=X)
        adata.var_names = var_names

        inferred = ExperimentContext.infer_from_data(adata)
        assert inferred["organism"]["species"] == "Homo sapiens"
        assert inferred["organism"]["inferred"] is True
        assert inferred["_data_shape"]["n_cells"] == 10

    def test_infer_mouse(self):
        import anndata as ad
        import numpy as np
        from scipy.sparse import csr_matrix

        # Mouse genes: lowercase-start
        var_names = [f"Gene{i}" for i in range(50)]
        X = csr_matrix(np.ones((10, 50), dtype=np.float32))
        adata = ad.AnnData(X=X)
        adata.var_names = var_names

        inferred = ExperimentContext.infer_from_data(adata)
        assert inferred["organism"]["species"] == "Mus musculus"

    def test_infer_umi_platform(self):
        import anndata as ad
        import numpy as np
        from scipy.sparse import csr_matrix

        # Integer counts → UMI-based
        X = csr_matrix(np.array([[1, 0, 3], [0, 2, 0]], dtype=np.float32))
        adata = ad.AnnData(X=X, var={"gene_id": ["A", "B", "C"]})
        adata.var_names = ["GENE1", "GENE2", "GENE3"]

        inferred = ExperimentContext.infer_from_data(adata)
        assert inferred["library"]["umi"] is True
        assert inferred["platform"]["vendor"] == "10x Genomics"


class TestDecisionHelpers:
    def test_needs_batch_correction_multi_sample(self, ctx):
        ctx.samples = [{"id": "a"}, {"id": "b"}]
        assert ctx.needs_batch_correction()

    def test_needs_batch_correction_single_sample(self, ctx):
        ctx.samples = [{"id": "a"}]
        ctx.design = {}
        assert not ctx.needs_batch_correction()

    def test_needs_pseudobulk_de(self, ctx):
        ctx.paradigm = "disease_vs_healthy"
        assert ctx.needs_pseudobulk_de()

        ctx.paradigm = "cell_atlas"
        assert not ctx.needs_pseudobulk_de()

    def test_needs_trajectory(self, ctx):
        ctx.paradigm = "developmental_trajectory"
        assert ctx.needs_trajectory()

        ctx.paradigm = "cell_atlas"
        assert not ctx.needs_trajectory()


class TestSaveLoadRoundtrip:
    def test_roundtrip(self, tmp_dir):
        ctx1 = ExperimentContext(tmp_dir)
        ctx1.paradigm = "cell_atlas"
        ctx1.organism = {"species": "Homo sapiens", "ncbi_taxon": 9606}
        ctx1.tissue = {"name": "PBMCs"}
        ctx1.hypotheses = ["test hypothesis"]
        ctx1.save()

        ctx2 = ExperimentContext(tmp_dir)
        assert ctx2.paradigm == "cell_atlas"
        assert ctx2.organism["species"] == "Homo sapiens"
        assert ctx2.tissue["name"] == "PBMCs"
        assert ctx2.hypotheses == ["test hypothesis"]

    def test_file_is_valid_json(self, tmp_dir):
        ctx = ExperimentContext(tmp_dir)
        ctx.paradigm = "cell_atlas"
        ctx.organism = {"species": "Homo sapiens"}
        ctx.tissue = {"name": "brain"}
        path = ctx.save()
        data = json.loads(path.read_text())
        assert data["paradigm"] == "cell_atlas"


class TestIsComplete:
    def test_complete(self, ctx):
        ctx.paradigm = "cell_atlas"
        ctx.organism = {"species": "Homo sapiens"}
        ctx.tissue = {"name": "PBMCs"}
        assert ctx.is_complete()

    def test_incomplete(self, ctx):
        ctx.paradigm = "cell_atlas"
        # Missing tissue and organism
        assert not ctx.is_complete()
