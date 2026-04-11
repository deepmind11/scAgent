#!/usr/bin/env python3
"""Unit tests for scagent.provenance."""

import json
import tempfile
from pathlib import Path

import pytest

from scagent.provenance import ProvenanceGraph, record_step, record_custom


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_dir(tmp_path):
    """Provide a temp .scagent directory."""
    d = tmp_path / ".scagent"
    d.mkdir()
    return d


@pytest.fixture
def graph(tmp_dir):
    return ProvenanceGraph(tmp_dir)


def _fake_tool_result(tool_id: str = "filter_cells", **extra_prov) -> dict:
    """Return a minimal tool result dict matching scagent/tools/ conventions."""
    prov = {
        "tool_id": tool_id,
        "parameters": {"min_genes": 200, "max_genes": 5000, "max_pct_mito": 10.0},
        **extra_prov,
    }
    return {"metrics": {}, "plots": [], "provenance": prov, "warnings": []}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestRecordSingleStep:
    def test_basic_record(self, graph):
        aid = graph.record("filter_cells", {"min_genes": 200})
        assert aid.startswith("sca:activity_filter_cells_")
        assert graph.n_activities == 1
        # 2 entities: synthetic input + output
        assert graph.n_entities == 2

    def test_activity_retrievable(self, graph):
        aid = graph.record("filter_cells", {"min_genes": 200}, user_prompt="filter")
        a = graph.get_activity(aid)
        assert a is not None
        assert a["tool_id"] == "filter_cells"
        assert a["parameters"]["min_genes"] == 200
        assert a["user_prompt"] == "filter"

    def test_extras_stored(self, graph):
        aid = graph.record(
            "filter_cells",
            {"min_genes": 200},
            extras={"cells_before": 11769, "cells_after": 10834},
        )
        a = graph.get_activity(aid)
        assert a["extras"]["cells_before"] == 11769
        assert a["extras"]["cells_after"] == 10834


class TestRecordChain:
    def test_five_step_chain(self, graph):
        tools = [
            ("load_10x_h5", {"filename": "test.h5"}),
            ("filter_cells", {"min_genes": 200}),
            ("filter_genes", {"min_cells": 3}),
            ("log_normalize", {"target_sum": 10000}),
            ("pca", {"n_comps": 50}),
        ]
        for tool_id, params in tools:
            graph.record(tool_id, params)

        assert graph.n_activities == 5
        # 5 outputs + 1 synthetic input = 6
        assert graph.n_entities == 6

    def test_entity_chain_linked(self, graph):
        graph.record("load_10x_h5", {"filename": "test.h5"})
        graph.record("filter_cells", {"min_genes": 200})
        graph.record("pca", {"n_comps": 50})

        # Each output entity should derive from the previous
        output_entities = [e for e in graph._entities if e.step_index >= 0]
        for i in range(1, len(output_entities)):
            assert output_entities[i].derived_from == output_entities[i - 1].id


class TestSerializeJsonld:
    def test_valid_structure(self, graph):
        graph.record("filter_cells", {"min_genes": 200})
        doc = graph.serialize()

        assert "@context" in doc
        assert doc["@context"]["prov"] == "http://www.w3.org/ns/prov#"
        assert "@graph" in doc
        assert isinstance(doc["@graph"], list)

    def test_node_types_present(self, graph):
        graph.record("filter_cells", {"min_genes": 200})
        doc = graph.serialize()
        types = set()
        for node in doc["@graph"]:
            t = node["@type"]
            if isinstance(t, list):
                types.update(t)
            else:
                types.add(t)
        assert "prov:Entity" in types
        assert "prov:Activity" in types
        assert "prov:Agent" in types

    def test_activity_has_required_fields(self, graph):
        graph.record("filter_cells", {"min_genes": 200}, user_prompt="do filter")
        doc = graph.serialize()
        activities = [n for n in doc["@graph"] if n.get("@type") == "prov:Activity"]
        assert len(activities) == 1
        a = activities[0]
        assert "prov:used" in a
        assert "prov:wasAssociatedWith" in a
        assert "prov:startedAtTime" in a
        assert "sca:parameters" in a
        assert "sca:user_prompt" in a


class TestSaveLoadRoundtrip:
    def test_roundtrip(self, tmp_dir):
        g1 = ProvenanceGraph(tmp_dir)
        g1.record("load_10x_h5", {"filename": "test.h5"}, user_prompt="load")
        g1.record("filter_cells", {"min_genes": 200}, extras={"cells_before": 100, "cells_after": 90})
        g1.record("pca", {"n_comps": 50})
        path = g1.save()

        g2 = ProvenanceGraph(tmp_dir)
        assert g2.n_activities == 3
        assert g2.n_entities == g1.n_entities

        # Check activity content survived
        chain = g2.get_chain()
        assert chain[0]["tool_id"] == "load_10x_h5"
        assert chain[1]["tool_id"] == "filter_cells"
        assert chain[1]["extras"]["cells_before"] == 100
        assert chain[2]["tool_id"] == "pca"

    def test_load_classmethod(self, tmp_dir):
        g1 = ProvenanceGraph(tmp_dir)
        g1.record("pca", {"n_comps": 50})

        g2 = ProvenanceGraph.load(tmp_dir)
        assert g2.n_activities == 1

    def test_file_is_valid_json(self, tmp_dir):
        g = ProvenanceGraph(tmp_dir)
        g.record("pca", {"n_comps": 50})
        path = g.save()
        data = json.loads(path.read_text())
        assert "@context" in data


class TestListActivities:
    def test_filter_by_tool_id(self, graph):
        graph.record("filter_cells", {"min_genes": 200})
        graph.record("pca", {"n_comps": 50})
        graph.record("filter_cells", {"min_genes": 300})

        results = graph.list_activities(tool_id="filter_cells")
        assert len(results) == 2
        assert all(r["tool_id"] == "filter_cells" for r in results)

    def test_filter_by_branch(self, graph):
        graph.record("pca", {"n_comps": 50}, branch="main")
        graph.record("pca", {"n_comps": 30}, branch="alt")

        assert len(graph.list_activities(branch="main")) == 1
        assert len(graph.list_activities(branch="alt")) == 1

    def test_no_filter_returns_all(self, graph):
        graph.record("load_10x_h5", {"filename": "a.h5"})
        graph.record("pca", {"n_comps": 50})
        assert len(graph.list_activities()) == 2


class TestSummary:
    def test_produces_markdown(self, graph):
        graph.record("load_10x_h5", {"filename": "test.h5"})
        graph.record("filter_cells", {"min_genes": 200, "max_genes": 5000})

        md = graph.summary()
        assert "## Analysis Provenance" in md
        assert "load_10x_h5" in md
        assert "filter_cells" in md
        assert "| # |" in md  # table header

    def test_empty_graph(self, graph):
        md = graph.summary()
        assert "No provenance" in md


class TestReplayPlan:
    def test_returns_tool_param_tuples(self, graph):
        graph.record("load_10x_h5", {"filename": "test.h5"})
        graph.record("filter_cells", {"min_genes": 200})
        graph.record("pca", {"n_comps": 50})

        plan = graph.replay_plan()
        assert len(plan) == 3
        assert plan[0] == ("load_10x_h5", {"filename": "test.h5"})
        assert plan[1] == ("filter_cells", {"min_genes": 200})
        assert plan[2] == ("pca", {"n_comps": 50})


class TestDiffBranches:
    def test_detects_divergence(self, graph):
        # Shared prefix
        graph.record("load_10x_h5", {"filename": "test.h5"}, branch="main")
        graph.record("filter_cells", {"min_genes": 200}, branch="main")

        graph.record("load_10x_h5", {"filename": "test.h5"}, branch="alt")
        graph.record("filter_cells", {"min_genes": 200}, branch="alt")

        # Diverge
        graph.record("pca", {"n_comps": 50}, branch="main")
        graph.record("pca", {"n_comps": 30}, branch="alt")

        d = graph.diff("main", "alt")
        assert d["shared_steps"] == 2
        assert len(d["branch_a_only"]) == 1
        assert len(d["branch_b_only"]) == 1
        assert d["parameter_diffs"][0]["param"] == "n_comps"

    def test_identical_branches(self, graph):
        for b in ("main", "alt"):
            graph.record("pca", {"n_comps": 50}, branch=b)

        d = graph.diff("main", "alt")
        assert d["shared_steps"] == 1
        assert len(d["branch_a_only"]) == 0
        assert len(d["branch_b_only"]) == 0


class TestRecordStepAdapter:
    def test_extracts_provenance(self, graph):
        result = _fake_tool_result(
            "filter_cells",
            cells_before=11769,
            cells_after=10834,
        )
        aid = record_step(graph, result, user_prompt="filter cells")

        a = graph.get_activity(aid)
        assert a["tool_id"] == "filter_cells"
        assert a["parameters"]["min_genes"] == 200
        assert a["extras"]["cells_before"] == 11769
        assert a["extras"]["cells_after"] == 10834

    def test_raises_without_provenance_key(self, graph):
        with pytest.raises(ValueError, match="no 'provenance' key"):
            record_step(graph, {"metrics": {}})


class TestRecordCustom:
    def test_basic_custom_step(self, graph):
        code = "adata.obs['mito_ribo'] = adata.obs['pct_counts_mt'] / adata.obs['pct_counts_ribo']"
        aid = record_custom(
            graph,
            description="mito/ribo ratio column",
            code=code,
            user_prompt="compute the ratio of mito to ribo genes",
        )
        a = graph.get_activity(aid)
        assert a is not None
        assert a["tool_id"] == "custom"
        assert a["parameters"]["description"] == "mito/ribo ratio column"
        assert a["extras"]["code"] == code

    def test_custom_with_effects(self, graph):
        aid = record_custom(
            graph,
            description="remove low-quality cells",
            code="adata = adata[adata.obs['score'] > 0.5].copy()",
            effects={"cells_before": 1000, "cells_after": 850, "added_obs_columns": []},
            user_prompt="filter cells with score below 0.5",
        )
        a = graph.get_activity(aid)
        assert a["extras"]["effects"]["cells_before"] == 1000
        assert a["extras"]["effects"]["cells_after"] == 850

    def test_custom_in_replay_plan(self, graph):
        graph.record("pca", {"n_comps": 50})
        record_custom(graph, description="custom filter", code="adata = adata[mask]")
        graph.record("leiden_clustering", {"resolution": 1.0})

        plan = graph.replay_plan()
        assert len(plan) == 3
        assert plan[1] == ("custom", {"description": "custom filter"})

    def test_custom_in_summary(self, graph):
        record_custom(graph, description="add score column", code="x = 1")
        md = graph.summary()
        assert "custom" in md
        assert "add score column" in md


class TestSizeBudget:
    def test_15_steps_under_50kb(self, tmp_dir):
        graph = ProvenanceGraph(tmp_dir)
        tools = [
            ("load_10x_h5", {"filename": "filtered_feature_bc_matrix.h5", "gex_only": True}),
            ("calculate_qc_metrics", {"species": "human"}),
            ("filter_cells", {"min_genes": 200, "max_genes": 5000, "max_pct_mito": 10.0, "min_counts": 500}),
            ("filter_genes", {"min_cells": 3}),
            ("scrublet_doublets", {"expected_doublet_rate": 0.06, "random_state": 0}),
            ("log_normalize", {"target_sum": 10000, "exclude_highly_expressed": False}),
            ("highly_variable_genes", {"n_top_genes": 2000, "flavor": "seurat_v3"}),
            ("pca", {"n_comps": 50, "random_state": 0, "svd_solver": "arpack"}),
            ("neighbor_graph", {"n_neighbors": 15, "n_pcs": 30, "metric": "euclidean"}),
            ("umap", {"min_dist": 0.5, "spread": 1.0, "random_state": 0}),
            ("harmony", {"key": "batch", "basis": "X_pca"}),
            ("leiden_clustering", {"resolution": 1.0, "random_state": 42, "key_added": "leiden_1.0"}),
            ("wilcoxon_markers", {"groupby": "leiden_1.0", "n_genes": 100, "method": "wilcoxon"}),
            ("celltypist_annotation", {"model": "Immune_All_Low.pkl", "majority_voting": True}),
            ("manual_annotation", {"groupby": "leiden_1.0", "marker_dict": {"T cell": ["CD3D"], "B cell": ["CD79A"]}}),
        ]
        for tool_id, params in tools:
            graph.record(tool_id, params, user_prompt=f"run {tool_id}")

        path = graph.save()
        size = path.stat().st_size
        assert size < 50_000, f"provenance.jsonld is {size} bytes, expected < 50KB"
        print(f"  15-step provenance file: {size:,} bytes")


class TestIdempotentSave:
    def test_double_save_same_content(self, tmp_dir):
        graph = ProvenanceGraph(tmp_dir)
        graph.record("pca", {"n_comps": 50})
        p1 = graph.save()
        content1 = p1.read_text()
        p2 = graph.save()
        content2 = p2.read_text()
        assert content1 == content2


class TestEmptyGraph:
    def test_no_crashes(self, graph):
        assert graph.n_activities == 0
        assert graph.n_entities == 0
        assert graph.list_activities() == []
        assert graph.get_chain() == []
        assert graph.replay_plan() == []
        assert "No provenance" in graph.summary()
        assert graph.get_activity("nonexistent") is None
        assert graph.get_entity("nonexistent") is None


class TestSessions:
    def test_session_created_on_init(self, graph):
        assert len(graph.sessions) == 1
        s = graph.sessions[0]
        assert s["session_number"] == 1
        assert s["ended_at"] is None  # still active
        assert "python" in s["software_versions"]

    def test_activities_tagged_to_session(self, graph):
        graph.record("pca", {"n_comps": 50})
        graph.record("leiden_clustering", {"resolution": 1.0})
        s = graph.sessions[0]
        assert len(s["activities"]) == 2

    def test_end_session(self, graph):
        graph.record("pca", {"n_comps": 50})
        graph.end_session()
        s = graph.sessions[0]
        assert s["ended_at"] is not None

    def test_multi_session_roundtrip(self, tmp_dir):
        # Session 1
        g1 = ProvenanceGraph(tmp_dir)
        g1.record("load_10x_h5", {"filename": "test.h5"})
        g1.record("filter_cells", {"min_genes": 200})
        g1.end_session()

        # Session 2 — new graph instance loads the file
        g2 = ProvenanceGraph(tmp_dir)
        g2.record("pca", {"n_comps": 50})

        assert len(g2.sessions) == 2
        assert g2.sessions[0]["session_number"] == 1
        assert g2.sessions[0]["ended_at"] is not None
        assert len(g2.sessions[0]["activities"]) == 2
        assert g2.sessions[1]["session_number"] == 2
        assert g2.sessions[1]["ended_at"] is None  # still active
        assert len(g2.sessions[1]["activities"]) == 1

    def test_session_boundary_in_summary(self, tmp_dir):
        g1 = ProvenanceGraph(tmp_dir)
        g1.record("load_10x_h5", {"filename": "test.h5"})
        g1.end_session()

        g2 = ProvenanceGraph(tmp_dir)
        g2.record("pca", {"n_comps": 50})

        md = g2.summary()
        assert "Session 1" in md
        assert "Session 2" in md
        assert "2 sessions" in md

    def test_session_serialized_in_jsonld(self, tmp_dir):
        g = ProvenanceGraph(tmp_dir)
        g.record("pca", {"n_comps": 50})
        doc = g.serialize()
        session_nodes = [n for n in doc["@graph"] if n.get("@type") == "sca:Session"]
        assert len(session_nodes) == 1
        assert session_nodes[0]["sca:session_number"] == 1
        assert "prov:startedAtTime" in session_nodes[0]


class TestForkBranch:
    def test_fork_creates_branch(self, graph):
        graph.record("pca", {"n_comps": 50}, branch="main")
        graph.fork_branch("experiment", from_branch="main")

        assert "experiment" in graph.branches
        # Recording on new branch works
        graph.record("leiden_clustering", {"resolution": 2.0}, branch="experiment")
        assert len(graph.get_chain("experiment")) == 1

    def test_fork_nonexistent_branch_raises(self, graph):
        with pytest.raises(ValueError, match="does not exist"):
            graph.fork_branch("new", from_branch="nonexistent")

    def test_fork_duplicate_branch_raises(self, graph):
        graph.record("pca", {"n_comps": 50}, branch="main")
        with pytest.raises(ValueError, match="already exists"):
            graph.fork_branch("main", from_branch="main")


class TestFullChain:
    def test_unfrosted_branch_same_as_chain(self, graph):
        """A branch that was never forked: get_full_chain == get_chain."""
        graph.record("load_10x_h5", {"filename": "test.h5"})
        graph.record("pca", {"n_comps": 50})

        assert graph.get_full_chain("main") == graph.get_chain("main")

    def test_forked_branch_includes_parent_prefix(self, graph):
        """A forked branch's full chain includes the parent's steps."""
        # Build main: load → filter → pca → neighbors
        graph.record("load_10x_h5", {"filename": "test.h5"}, branch="main")
        graph.record("filter_cells", {"min_genes": 200}, branch="main")
        graph.record("pca", {"n_comps": 50}, branch="main")
        graph.record("neighbor_graph", {"n_neighbors": 15}, branch="main")

        # Fork at this point
        graph.fork_branch("high_res", from_branch="main")

        # Continue on the fork
        graph.record("leiden_clustering", {"resolution": 2.0}, branch="high_res")
        graph.record("wilcoxon_markers", {"groupby": "leiden"}, branch="high_res")

        # Full chain of fork should include main's 4 steps + fork's 2
        full = graph.get_full_chain("high_res")
        tools = [a["tool_id"] for a in full]
        assert tools == [
            "load_10x_h5", "filter_cells", "pca", "neighbor_graph",
            "leiden_clustering", "wilcoxon_markers",
        ]

        # Plain get_chain should only have the fork's 2 steps
        own = graph.get_chain("high_res")
        assert len(own) == 2

    def test_replay_plan_full_vs_own(self, graph):
        graph.record("load_10x_h5", {"filename": "test.h5"}, branch="main")
        graph.record("pca", {"n_comps": 50}, branch="main")
        graph.fork_branch("alt", from_branch="main")
        graph.record("leiden_clustering", {"resolution": 0.5}, branch="alt")

        plan_full = graph.replay_plan("alt", full=True)
        plan_own = graph.replay_plan("alt", full=False)

        assert len(plan_full) == 3  # load + pca + leiden
        assert len(plan_own) == 1   # just leiden
        assert plan_full[0][0] == "load_10x_h5"
        assert plan_full[-1][0] == "leiden_clustering"


class TestPromoteBranch:
    def test_promote_and_export(self, graph):
        graph.record("load_10x_h5", {"filename": "test.h5"}, branch="main")
        graph.record("pca", {"n_comps": 50}, branch="main")
        graph.fork_branch("final", from_branch="main")
        graph.record("leiden_clustering", {"resolution": 1.0}, branch="final")

        assert graph.promoted_branch is None

        graph.promote_branch("final")
        assert graph.promoted_branch == "final"

        # export_plan should give the full chain for the promoted branch
        plan = graph.export_plan()
        assert len(plan) == 3
        assert plan[-1][0] == "leiden_clustering"

    def test_promote_survives_roundtrip(self, tmp_dir):
        g1 = ProvenanceGraph(tmp_dir)
        g1.record("pca", {"n_comps": 50}, branch="main")
        g1.fork_branch("final", from_branch="main")
        g1.record("leiden_clustering", {"resolution": 1.0}, branch="final")
        g1.promote_branch("final")
        g1.end_session()

        g2 = ProvenanceGraph(tmp_dir)
        assert g2.promoted_branch == "final"

    def test_export_falls_back_to_main(self, graph):
        graph.record("pca", {"n_comps": 50})
        plan = graph.export_plan()
        assert len(plan) == 1  # main's single step

    def test_promote_nonexistent_raises(self, graph):
        with pytest.raises(ValueError, match="does not exist"):
            graph.promote_branch("nonexistent")
