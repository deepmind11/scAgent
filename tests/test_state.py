#!/usr/bin/env python3
"""Unit tests for scagent.state."""

import json

import anndata as ad
import numpy as np
import pytest
from scipy.sparse import csr_matrix

from scagent.state import StateManager, _hash_file


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tmp_dir(tmp_path):
    """Provide a temp .scagent directory."""
    return tmp_path / ".scagent"


@pytest.fixture
def sm(tmp_dir):
    return StateManager(tmp_dir)


def _make_adata(n_obs: int = 100, n_vars: int = 50, seed: int = 0) -> ad.AnnData:
    """Create a small test AnnData with sparse X."""
    rng = np.random.default_rng(seed)
    X = csr_matrix(rng.poisson(3, size=(n_obs, n_vars)).astype(np.float32))
    obs_names = [f"cell_{i}" for i in range(n_obs)]
    var_names = [f"gene_{i}" for i in range(n_vars)]
    return ad.AnnData(X=X, obs={"cell_id": obs_names}, var={"gene_id": var_names})


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestSaveLoadRoundtrip:
    def test_basic_roundtrip(self, sm):
        adata = _make_adata()
        h = sm.save_snapshot(adata, step_name="after_load")
        assert len(h) == 8

        loaded = sm.load_snapshot()
        assert loaded.shape == adata.shape
        np.testing.assert_array_equal(loaded.X.toarray(), adata.X.toarray())

    def test_step_info_in_head(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="after_pca", step_index=5)
        state = sm.current_state()
        assert state.step_name == "after_pca"
        assert state.step_index == 5
        assert not state.has_unsaved_changes


class TestHashing:
    def test_deterministic(self, sm):
        adata = _make_adata()
        h1 = sm.save_snapshot(adata, step_name="s1")
        h2 = sm.save_snapshot(adata, step_name="s2")
        assert h1 == h2  # same data → same hash

    def test_changes_on_mutation(self, sm):
        adata1 = _make_adata(seed=0)
        adata2 = _make_adata(seed=42)
        h1 = sm.save_snapshot(adata1, step_name="s1")
        h2 = sm.save_snapshot(adata2, step_name="s2")
        assert h1 != h2


class TestCreateBranch:
    def test_creates_directory_structure(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="after_pca")

        fork_hash = sm.create_branch("experiment")

        branch_dir = sm._branch_dir("experiment")
        assert branch_dir.is_dir()
        assert (branch_dir / "parent.json").exists()
        assert (branch_dir / "head.json").exists()
        assert (branch_dir / "snapshots" / f"{fork_hash}.h5ad").exists()

    def test_parent_metadata(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="after_neighbors", step_index=6)

        sm.create_branch("high_res")
        parent = json.loads(
            (sm._branch_dir("high_res") / "parent.json").read_text()
        )
        assert parent["parent_branch"] == "main"
        assert parent["fork_step"] == "after_neighbors"
        assert parent["fork_step_index"] == 6

    def test_duplicate_name_raises(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="s")
        sm.create_branch("test")
        with pytest.raises(ValueError, match="already exists"):
            sm.create_branch("test")

    def test_nonexistent_source_raises(self, sm):
        with pytest.raises(ValueError, match="does not exist"):
            sm.create_branch("new", from_branch="nonexistent")

    def test_does_not_switch(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="s")
        sm.create_branch("experiment")
        assert sm.active_branch == "main"  # still on main

    def test_create_from_unsaved_state(self, sm):
        adata = _make_adata()
        # No snapshot on main yet — must pass adata
        fork_hash = sm.create_branch("experiment", adata=adata)
        assert fork_hash is not None
        assert len(fork_hash) == 8


class TestSwitchBranch:
    def test_switch_loads_correct_state(self, sm):
        adata_main = _make_adata(seed=0)
        sm.save_snapshot(adata_main, step_name="main_state")

        sm.create_branch("alt")
        adata_alt = _make_adata(seed=99)
        # Save different data on alt
        sm._active_branch = "alt"  # direct set for test
        sm.save_snapshot(adata_alt, step_name="alt_state", branch="alt")
        sm._active_branch = "main"

        loaded = sm.switch_branch("alt")
        assert sm.active_branch == "alt"
        np.testing.assert_array_equal(loaded.X.toarray(), adata_alt.X.toarray())

    def test_switch_saves_dirty(self, sm):
        adata = _make_adata(seed=0)
        sm.save_snapshot(adata, step_name="initial")
        sm.create_branch("alt")

        # Modify in-memory state
        adata_modified = _make_adata(seed=42)
        sm.mark_dirty("modified")

        # Switch away — should save the dirty state
        sm.switch_branch("alt", adata=adata_modified)

        # Switch back and verify
        loaded = sm.switch_branch("main")
        np.testing.assert_array_equal(
            loaded.X.toarray(), adata_modified.X.toarray()
        )

    def test_switch_nonexistent_raises(self, sm):
        with pytest.raises(ValueError, match="does not exist"):
            sm.switch_branch("nonexistent")


class TestDeleteBranch:
    def test_delete_removes_directory(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="s")
        sm.create_branch("temp")
        assert sm._branch_exists("temp")

        sm.delete_branch("temp")
        assert not sm._branch_exists("temp")

    def test_delete_active_raises(self, sm):
        with pytest.raises(ValueError, match="active branch"):
            sm.delete_branch("main")

    def test_delete_main_raises(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="s")
        sm.create_branch("alt")
        sm.switch_branch("alt")
        with pytest.raises(ValueError, match="Cannot delete.*main"):
            sm.delete_branch("main")

    def test_delete_with_children_raises(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="s")
        sm.create_branch("parent_branch")

        # Create a child of parent_branch
        sm.switch_branch("parent_branch")
        adata2 = _make_adata(seed=1)
        sm.save_snapshot(adata2, step_name="s2", branch="parent_branch")
        sm.create_branch("child", from_branch="parent_branch")
        sm.switch_branch("main")

        with pytest.raises(ValueError, match="child branches"):
            sm.delete_branch("parent_branch")

        # Force works
        sm.delete_branch("parent_branch", force=True)
        assert not sm._branch_exists("parent_branch")


class TestListBranches:
    def test_lists_all_branches(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="s")
        sm.create_branch("alpha")
        sm.create_branch("beta")

        branches = sm.list_branches()
        names = [b.name for b in branches]
        assert "main" in names
        assert "alpha" in names
        assert "beta" in names
        assert len(names) == 3

    def test_active_branch_marked(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="s")
        branches = sm.list_branches()
        active = [b for b in branches if b.is_active]
        assert len(active) == 1
        assert active[0].name == "main"

    def test_branch_info_fields(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="after_pca", step_index=5)
        sm.create_branch("test")

        branches = {b.name: b for b in sm.list_branches()}
        test_branch = branches["test"]
        assert test_branch.parent_branch == "main"
        assert test_branch.fork_step == "after_pca"
        assert test_branch.snapshot_count >= 1
        assert test_branch.disk_size_bytes > 0


class TestCurrentState:
    def test_after_save(self, sm):
        adata = _make_adata()
        h = sm.save_snapshot(adata, step_name="after_leiden", step_index=8)

        state = sm.current_state()
        assert state.branch == "main"
        assert state.hash == h
        assert state.step_name == "after_leiden"
        assert state.step_index == 8
        assert not state.has_unsaved_changes

    def test_dirty_tracking(self, sm):
        adata = _make_adata()
        sm.save_snapshot(adata, step_name="s")
        sm.mark_dirty("after_clustering", step_index=9)

        state = sm.current_state()
        assert state.has_unsaved_changes


class TestSessionSaveLoad:
    def test_roundtrip(self, tmp_dir):
        sm1 = StateManager(tmp_dir)
        adata = _make_adata()
        sm1.save_on_exit(adata, step_name="after_annotation")

        sm2 = StateManager(tmp_dir)
        loaded = sm2.load_on_start()
        assert loaded is not None
        assert loaded.shape == adata.shape

    def test_load_on_start_empty(self, tmp_dir):
        sm = StateManager(tmp_dir)
        result = sm.load_on_start()
        assert result is None


class TestGcSnapshots:
    def test_removes_unreferenced(self, sm):
        adata1 = _make_adata(seed=0)
        adata2 = _make_adata(seed=1)
        adata3 = _make_adata(seed=2)

        sm.save_snapshot(adata1, step_name="s1")
        sm.save_snapshot(adata2, step_name="s2")
        h3 = sm.save_snapshot(adata3, step_name="s3")

        # Only h3 should be referenced (latest HEAD)
        snap_dir = sm._branch_dir("main") / "snapshots"
        # After save_snapshot, gc_branch already ran — only HEAD should remain
        remaining = list(snap_dir.glob("*.h5ad"))
        assert len(remaining) == 1
        assert remaining[0].stem == h3


class TestMultipleSavesOverwrite:
    def test_old_snapshot_replaced(self, sm):
        adata1 = _make_adata(seed=0)
        adata2 = _make_adata(seed=42)
        h1 = sm.save_snapshot(adata1, step_name="s1")
        h2 = sm.save_snapshot(adata2, step_name="s2")
        assert h1 != h2

        # Only latest snapshot should remain
        snap_dir = sm._branch_dir("main") / "snapshots"
        remaining = [f.stem for f in snap_dir.glob("*.h5ad")]
        assert h2 in remaining
        assert h1 not in remaining


class TestProvenanceIntegration:
    def test_hash_usable_in_record_step(self, sm):
        """save_snapshot returns a hash that can be passed to record_step."""
        from scagent.provenance import ProvenanceGraph, record_step

        adata = _make_adata()
        file_hash = sm.save_snapshot(adata, step_name="after_pca")
        assert isinstance(file_hash, str)
        assert len(file_hash) == 8

        # Use it in provenance
        graph = ProvenanceGraph(sm._dir)
        result = {
            "metrics": {},
            "plots": [],
            "provenance": {"tool_id": "pca", "parameters": {"n_comps": 50}},
            "warnings": [],
        }
        aid = record_step(graph, result, output_hash=file_hash)
        activity = graph.get_activity(aid)
        assert activity is not None
