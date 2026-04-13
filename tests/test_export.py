#!/usr/bin/env python3
"""Unit tests for scagent.export (methods + reproducibility package)."""

import json
import pytest
from pathlib import Path
from scagent.provenance import ProvenanceGraph
from scagent.context import ExperimentContext
from scagent.export import generate_methods, generate_repro_package


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def graph(tmp_path):
    """Provenance graph with a realistic PBMC pipeline."""
    g = ProvenanceGraph(tmp_path / "prov")

    steps = [
        {
            "tool_id": "load_10x_h5",
            "parameters": {"file_path": "pbmc_10k_v3_filtered.h5"},
            "extras": {"effects": {"n_cells": 11769, "n_genes": 33538},
                       "description": "Load PBMC 10k dataset"},
        },
        {
            "tool_id": "qc_metrics",
            "parameters": {},
            "extras": {},
        },
        {
            "tool_id": "filter_cells",
            "parameters": {"min_genes": 200, "max_genes": 5000, "max_pct_mito": 20},
            "extras": {"effects": {"n_cells_before": 11769, "n_cells_after": 10834}},
        },
        {
            "tool_id": "filter_genes",
            "parameters": {"min_cells": 3},
            "extras": {"effects": {"n_genes_before": 33538, "n_genes_after": 18952}},
        },
        {
            "tool_id": "detect_doublets",
            "parameters": {"expected_doublet_rate": 0.06},
            "extras": {"effects": {"n_doublets": 347, "n_cells": 10834}},
        },
        {
            "tool_id": "normalize",
            "parameters": {"target_sum": 10000},
            "extras": {},
        },
        {
            "tool_id": "log_transform",
            "parameters": {},
            "extras": {},
        },
        {
            "tool_id": "highly_variable_genes",
            "parameters": {"flavor": "seurat_v3", "n_top_genes": 2000},
            "extras": {"effects": {"n_hvg": 2000}},
        },
        {
            "tool_id": "pca",
            "parameters": {"n_comps": 50},
            "extras": {},
        },
        {
            "tool_id": "neighbors",
            "parameters": {"n_neighbors": 15, "n_pcs": 50},
            "extras": {},
        },
        {
            "tool_id": "leiden",
            "parameters": {"resolution": 0.8},
            "extras": {"effects": {"n_clusters": 21}},
        },
        {
            "tool_id": "umap",
            "parameters": {},
            "extras": {},
        },
    ]
    for s in steps:
        g.record(
            tool_id=s["tool_id"],
            parameters=s["parameters"],
            extras=s.get("extras"),
            input_hash="",
            output_hash="",
            user_prompt="",
        )
    return g


@pytest.fixture
def context(tmp_path):
    ctx = ExperimentContext(tmp_path / "ctx")
    ctx.paradigm = "disease_vs_healthy"
    ctx.organism = {"species": "Homo sapiens", "ncbi_taxon": 9606}
    ctx.tissue = {"name": "PBMCs"}
    ctx.platform = {"name": "10x Chromium"}
    ctx.library = {"type": "3prime"}
    return ctx


# ---------------------------------------------------------------------------
# Methods section
# ---------------------------------------------------------------------------

class TestMethods:
    def test_generates_text(self, graph):
        text = generate_methods(graph)
        assert len(text) > 100
        assert isinstance(text, str)

    def test_includes_tool_steps(self, graph):
        text = generate_methods(graph)
        assert "filtered" in text.lower() or "filter" in text.lower()
        assert "200" in text  # min_genes
        assert "Leiden" in text or "leiden" in text
        assert "0.8" in text  # resolution

    def test_includes_context_preamble(self, graph, context):
        text = generate_methods(graph, context)
        assert "Homo sapiens" in text
        assert "PBMCs" in text

    def test_includes_normalization(self, graph):
        text = generate_methods(graph)
        assert "10000" in text or "10,000" in text  # target_sum

    def test_includes_pca(self, graph):
        text = generate_methods(graph)
        assert "50" in text  # n_comps

    def test_includes_hvg(self, graph):
        text = generate_methods(graph)
        assert "2000" in text or "2,000" in text

    def test_includes_doublets(self, graph):
        text = generate_methods(graph)
        assert "Scrublet" in text or "doublet" in text.lower()

    def test_without_context(self, graph):
        text = generate_methods(graph, context=None)
        assert len(text) > 50
        assert "Homo sapiens" not in text

    def test_empty_graph(self, tmp_path):
        g = ProvenanceGraph(tmp_path / "empty")
        text = generate_methods(g)
        assert isinstance(text, str)


# ---------------------------------------------------------------------------
# Reproducibility package
# ---------------------------------------------------------------------------

class TestReproPackage:
    def test_creates_directory(self, graph, context, tmp_path):
        out = generate_repro_package(graph, context, out_dir=tmp_path / "repro")
        assert out.exists()
        assert out.is_dir()

    def test_all_files_present(self, graph, context, tmp_path):
        out = generate_repro_package(graph, context, out_dir=tmp_path / "repro")
        expected = ["replay.py", "params.json", "environment.json",
                    "provenance.json", "methods.md", "context.json", "README.md"]
        for f in expected:
            assert (out / f).exists(), f"Missing {f}"

    def test_params_json_valid(self, graph, context, tmp_path):
        out = generate_repro_package(graph, context, out_dir=tmp_path / "repro")
        params = json.loads((out / "params.json").read_text())
        assert isinstance(params, list)
        assert len(params) == 12  # 12 steps
        assert params[0]["tool_id"] == "load_10x_h5"
        assert params[-1]["tool_id"] == "umap"

    def test_replay_script_is_python(self, graph, tmp_path):
        out = generate_repro_package(graph, out_dir=tmp_path / "repro")
        script = (out / "replay.py").read_text()
        assert script.startswith("#!/usr/bin/env python3")
        assert "import scanpy" in script
        assert "leiden" in script.lower()
        # Should be valid Python syntax
        compile(script, "replay.py", "exec")

    def test_methods_md_matches(self, graph, context, tmp_path):
        out = generate_repro_package(graph, context, out_dir=tmp_path / "repro")
        methods_file = (out / "methods.md").read_text()
        methods_direct = generate_methods(graph, context)
        assert methods_file == methods_direct

    def test_provenance_json_valid(self, graph, tmp_path):
        out = generate_repro_package(graph, out_dir=tmp_path / "repro")
        prov = json.loads((out / "provenance.json").read_text())
        assert "@context" in prov or "@graph" in prov

    def test_environment_json_has_timestamp(self, graph, tmp_path):
        out = generate_repro_package(graph, out_dir=tmp_path / "repro")
        env = json.loads((out / "environment.json").read_text())
        assert "exported_at" in env
        assert "branch" in env

    def test_without_context_no_context_json(self, graph, tmp_path):
        out = generate_repro_package(graph, context=None, out_dir=tmp_path / "repro")
        assert not (out / "context.json").exists()

    def test_readme_has_step_count(self, graph, tmp_path):
        out = generate_repro_package(graph, out_dir=tmp_path / "repro")
        readme = (out / "README.md").read_text()
        assert "12" in readme  # 12 steps
        assert "main" in readme  # branch

    def test_custom_branch(self, graph, tmp_path):
        out = generate_repro_package(graph, out_dir=tmp_path / "repro",
                                     branch="main")
        env = json.loads((out / "environment.json").read_text())
        assert env["branch"] == "main"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_custom_tool_in_methods(self, tmp_path):
        g = ProvenanceGraph(tmp_path / "prov")
        g.record(
            tool_id="custom",
            parameters={},
            extras={"description": "Manually adjusted cluster 5 labels"},
            input_hash="", output_hash="", user_prompt="",
        )
        text = generate_methods(g)
        assert "Manually adjusted" in text

    def test_unknown_tool_in_methods(self, tmp_path):
        g = ProvenanceGraph(tmp_path / "prov")
        g.record(
            tool_id="some_future_tool",
            parameters={"x": 1},
            extras={"description": "Novel analysis step"},
            input_hash="", output_hash="", user_prompt="",
        )
        text = generate_methods(g)
        assert "Novel analysis" in text

    def test_replay_script_handles_unknown_tool(self, tmp_path):
        g = ProvenanceGraph(tmp_path / "prov")
        g.record(
            tool_id="alien_tool",
            parameters={"x": 1},
            extras={},
            input_hash="", output_hash="", user_prompt="",
        )
        out = generate_repro_package(g, out_dir=tmp_path / "repro")
        script = (out / "replay.py").read_text()
        assert "TODO" in script
        assert "alien_tool" in script
