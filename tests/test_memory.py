#!/usr/bin/env python3
"""Unit tests for scagent.memory (MemPalace integration)."""

import pytest
from scagent.memory import (
    ProjectMemory, normalize_room, room_from_tool, _make_id,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mem(tmp_path):
    return ProjectMemory(tmp_path / "palace", project_name="test_proj")


@pytest.fixture
def populated(mem):
    """Palace with content on two branches."""
    mem.store_exchange("Why resolution 0.8?",
                       "NK markers GNLY/NKG7 clean at 0.8, fragmented at 1.2",
                       branch="main", room="clustering")
    mem.store_exchange("Should we filter doublets?",
                       "Yes, Scrublet found 3.2%. Remove before normalization.",
                       branch="main", room="qc")
    mem.store_step({
        "tool_id": "leiden",
        "parameters": {"resolution": 0.8},
        "effects": {"n_clusters": 21},
    }, branch="main")
    mem.store_decision("Use pseudobulk DE",
                       "Cell-level Wilcoxon inflates FP. 6 donors = replicates.",
                       branch="main")
    # Second branch
    mem.store_decision("Try resolution 1.2",
                       "Want finer CD4 T-cell subclusters",
                       branch="high_res")
    return mem


# ---------------------------------------------------------------------------
# Store + recall basics
# ---------------------------------------------------------------------------

class TestStore:
    def test_store_returns_drawer_id(self, mem):
        did = mem.store("hello", room="qc", branch="main")
        assert did.startswith("drawer_")

    def test_store_exchange(self, mem):
        did = mem.store_exchange("q?", "a.", branch="main")
        assert did.startswith("drawer_")

    def test_store_step(self, mem):
        did = mem.store_step({"tool_id": "pca", "parameters": {"n_comps": 50}},
                             branch="main")
        assert did.startswith("drawer_")

    def test_store_decision(self, mem):
        did = mem.store_decision("choice", "reason", branch="main")
        assert did.startswith("drawer_")


class TestRecall:
    def test_finds_relevant(self, populated):
        hits = populated.recall("resolution clustering", branch="main")
        assert len(hits) > 0
        assert any("0.8" in h["text"] for h in hits)

    def test_room_filter(self, populated):
        hits = populated.recall("filter", room="qc", branch="main")
        for h in hits:
            assert h["room"] == "qc"

    def test_empty_palace(self, mem):
        assert mem.recall("anything") == []

    def test_hit_has_branch(self, populated):
        hits = populated.recall("resolution", branch="main")
        assert all(h["branch"] == "main" for h in hits)

    def test_hit_has_similarity(self, populated):
        hits = populated.recall("doublet", branch="main")
        assert len(hits) > 0
        assert 0 <= hits[0]["similarity"] <= 1


# ---------------------------------------------------------------------------
# Branch as tag — the core contract
# ---------------------------------------------------------------------------

class TestBranching:
    def test_branch_isolation(self, populated):
        """main recall does NOT return high_res content."""
        hits = populated.recall("resolution", branch="main")
        assert all(h["branch"] == "main" for h in hits)

    def test_cross_branch_with_none(self, populated):
        """branch=None searches everything."""
        hits = populated.recall("resolution", branch=None)
        branches = {h["branch"] for h in hits}
        assert "main" in branches
        assert "high_res" in branches

    def test_specific_other_branch(self, populated):
        """Can explicitly query another branch."""
        hits = populated.recall("resolution", branch="high_res")
        assert len(hits) > 0
        assert all(h["branch"] == "high_res" for h in hits)

    def test_same_content_different_branches(self, mem):
        """Same text on two branches → separate drawer IDs."""
        id1 = mem.store("same text", branch="main")
        id2 = mem.store("same text", branch="dev")
        assert id1 != id2


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------

class TestStatus:
    def test_status_after_stores(self, populated):
        s = populated.status()
        assert s["project"] == "test_proj"
        assert s["project_drawers"] == 5
        assert "main" in s["branches"]
        assert "high_res" in s["branches"]
        assert s["branches"]["main"] == 4
        assert s["branches"]["high_res"] == 1
        assert s["sessions"] >= 1


# ---------------------------------------------------------------------------
# Room normalization
# ---------------------------------------------------------------------------

class TestNormalizeRoom:
    @pytest.mark.parametrize("alias,expected", [
        ("qc", "qc"), ("quality_control", "qc"), ("filtering", "qc"),
        ("normalization", "preprocessing"), ("hvg", "preprocessing"),
        ("leiden", "clustering"), ("umap", "clustering"),
        ("markers", "annotation"), ("cell_type", "annotation"),
        ("pseudobulk", "de"), ("deseq2", "de"),
        ("gsea", "enrichment"), ("pathway", "enrichment"),
        ("something_unknown", "exploration"),
    ])
    def test_aliases(self, alias, expected):
        assert normalize_room(alias) == expected


class TestRoomFromTool:
    @pytest.mark.parametrize("tool_id,expected", [
        ("load_10x_h5", "qc"), ("filter_cells", "qc"),
        ("normalize", "preprocessing"), ("pca", "preprocessing"),
        ("leiden", "clustering"), ("neighbors", "clustering"),
        ("deseq2_pseudobulk", "de"), ("gsea", "enrichment"),
        ("custom", "exploration"), ("unknown_tool", "exploration"),
    ])
    def test_mapping(self, tool_id, expected):
        assert room_from_tool(tool_id) == expected


# ---------------------------------------------------------------------------
# Deterministic IDs
# ---------------------------------------------------------------------------

class TestMakeId:
    def test_deterministic(self):
        assert _make_id("w", "r", "c") == _make_id("w", "r", "c")

    def test_different_content(self):
        assert _make_id("w", "r", "a") != _make_id("w", "r", "b")

    def test_different_branch(self):
        assert _make_id("w", "r", "c", "main") != _make_id("w", "r", "c", "dev")

    def test_format(self):
        assert _make_id("proj", "qc", "x").startswith("drawer_proj_qc_")
