#!/usr/bin/env python3
"""Unit tests for scagent.dag."""

import json

import pytest

from scagent.context import ExperimentContext
from scagent.dag import AnalysisDAG


@pytest.fixture
def tmp_dir(tmp_path):
    return tmp_path / ".scagent"


def _make_ctx(
    tmp_dir,
    paradigm="cell_atlas",
    conditions=None,
    samples=None,
    batch_key=None,
):
    """Helper: create a valid ExperimentContext."""
    ctx = ExperimentContext(tmp_dir)
    ctx.paradigm = paradigm
    ctx.organism = {"species": "Homo sapiens", "ncbi_taxon": 9606}
    ctx.tissue = {"name": "PBMCs"}
    if conditions:
        ctx.design = {"conditions": conditions, "batch_key": batch_key}
    if samples:
        ctx.samples = samples
    else:
        ctx.samples = [{"id": "sample1"}]
    return ctx


class TestCellAtlasDAG:
    def test_has_expected_steps(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)

        step_ids = [s.id for s in dag.steps]
        assert "load" in step_ids
        assert "qc_metrics" in step_ids
        assert "normalize" in step_ids
        assert "clustering" in step_ids
        assert "annotation" in step_ids

    def test_no_pseudobulk_de(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        step_ids = [s.id for s in dag.steps]
        assert "pseudobulk_de" not in step_ids

    def test_no_trajectory(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        step_ids = [s.id for s in dag.steps]
        assert "trajectory" not in step_ids


class TestDiseaseVsHealthyDAG:
    def test_has_pseudobulk_de(self, tmp_dir):
        ctx = _make_ctx(tmp_dir, paradigm="disease_vs_healthy",
                        conditions=["disease", "healthy"])
        dag = AnalysisDAG.from_context(ctx)
        step_ids = [s.id for s in dag.steps]
        assert "pseudobulk_de" in step_ids
        assert "pathway_enrichment" in step_ids

    def test_has_composition(self, tmp_dir):
        ctx = _make_ctx(tmp_dir, paradigm="disease_vs_healthy",
                        conditions=["disease", "healthy"])
        dag = AnalysisDAG.from_context(ctx)
        step_ids = [s.id for s in dag.steps]
        assert "composition" in step_ids


class TestTrajectoryDAG:
    def test_has_trajectory_steps(self, tmp_dir):
        ctx = _make_ctx(tmp_dir, paradigm="developmental_trajectory")
        dag = AnalysisDAG.from_context(ctx)
        step_ids = [s.id for s in dag.steps]
        assert "paga" in step_ids
        assert "pseudotime" in step_ids

    def test_no_pseudobulk_de(self, tmp_dir):
        ctx = _make_ctx(tmp_dir, paradigm="developmental_trajectory")
        dag = AnalysisDAG.from_context(ctx)
        step_ids = [s.id for s in dag.steps]
        assert "pseudobulk_de" not in step_ids


class TestConditionalBatchCorrection:
    def test_included_when_multi_sample(self, tmp_dir):
        ctx = _make_ctx(
            tmp_dir,
            samples=[{"id": "a"}, {"id": "b"}],
            batch_key="sample",
        )
        dag = AnalysisDAG.from_context(ctx)
        step_ids = [s.id for s in dag.steps]
        assert "batch_correction" in step_ids

    def test_excluded_when_single_sample(self, tmp_dir):
        ctx = _make_ctx(tmp_dir, samples=[{"id": "a"}])
        dag = AnalysisDAG.from_context(ctx)
        step_ids = [s.id for s in dag.steps]
        assert "batch_correction" not in step_ids


class TestNextStep:
    def test_first_step(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        nxt = dag.next_step()
        assert nxt is not None
        assert nxt.id == "load"

    def test_advances_after_complete(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        dag.complete_step("load")
        nxt = dag.next_step()
        assert nxt.id == "qc_metrics"

    def test_none_when_all_done(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        for s in dag.steps:
            dag.complete_step(s.id)
        assert dag.next_step() is None


class TestSkipStep:
    def test_skip_with_reason(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        dag.skip_step("load", reason="data already loaded")
        step = dag.get_step("load")
        assert step.status == "skipped"
        assert step.skip_reason == "data already loaded"

    def test_skip_unblocks_dependents(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        dag.skip_step("load")
        # qc_metrics depends on load — should be valid now
        assert dag.is_valid_step("qc_metrics")


class TestIsValidStep:
    def test_deps_not_met(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        # clustering depends on neighbors → not valid yet
        assert not dag.is_valid_step("clustering")

    def test_deps_met(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        # load has no deps → always valid
        assert dag.is_valid_step("load")

    def test_nonexistent_step_raises(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        with pytest.raises(ValueError, match="No step"):
            dag.is_valid_step("nonexistent")


class TestProgress:
    def test_progress_count(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        total_steps = len(dag.steps)

        assert dag.progress == (0, total_steps)
        dag.complete_step("load")
        assert dag.progress == (1, total_steps)

    def test_skipped_excluded_from_total(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        total = len(dag.steps)
        dag.skip_step("load")
        done, new_total = dag.progress
        assert new_total == total - 1


class TestSummary:
    def test_produces_markdown(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        dag.complete_step("load")

        md = dag.summary()
        assert "## Analysis Plan" in md
        assert "cell_atlas" in md
        assert "✅" in md
        assert "⬜" in md
        assert "Next step" in md


class TestSaveLoadRoundtrip:
    def test_roundtrip(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag1 = AnalysisDAG.from_context(ctx)
        dag1.complete_step("load")
        dag1.complete_step("qc_metrics")
        dag1.skip_step("doublet_detection", reason="low cell count")
        tmp_dir.mkdir(parents=True, exist_ok=True)
        dag1.save(tmp_dir)

        dag2 = AnalysisDAG.load(tmp_dir)
        assert dag2.paradigm == "cell_atlas"
        assert dag2.get_step("load").status == "done"
        assert dag2.get_step("qc_metrics").status == "done"
        assert dag2.get_step("doublet_detection").status == "skipped"
        assert dag2.get_step("doublet_detection").skip_reason == "low cell count"
        assert dag2.get_step("normalize").status == "pending"

    def test_file_is_valid_json(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        path = dag.save(tmp_dir)
        data = json.loads(path.read_text())
        assert data["paradigm"] == "cell_atlas"
        assert isinstance(data["steps"], list)


class TestAddStep:
    def test_add_step_at_end(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        from scagent.dag import DAGStep
        new_step = DAGStep(
            id="trajectory", name="Trajectory", category="trajectory",
            tool_id="paga", depends_on=["clustering", "annotation"],
        )
        dag.add_step(new_step)
        assert dag.get_step("trajectory") is not None
        assert dag.steps[-1].id == "trajectory"

    def test_add_step_after(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        from scagent.dag import DAGStep
        new_step = DAGStep(
            id="custom_filter", name="Custom filter", category="qc",
            depends_on=["filter_cells"],
        )
        dag.add_step(new_step, after="filter_cells")
        ids = [s.id for s in dag.steps]
        assert ids.index("custom_filter") == ids.index("filter_cells") + 1

    def test_add_duplicate_raises(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        from scagent.dag import DAGStep
        dup = DAGStep(id="load", name="Load again", category="loading")
        with pytest.raises(ValueError, match="already exists"):
            dag.add_step(dup)

    def test_add_after_unknown_raises(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        from scagent.dag import DAGStep
        step = DAGStep(id="new", name="New", category="misc")
        with pytest.raises(ValueError, match="not found"):
            dag.add_step(step, after="nonexistent")


class TestMarkPrecomputed:
    def test_mark_precomputed(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        dag.mark_precomputed("load")
        assert dag.get_step("load").status == "done"

    def test_mark_precomputed_unknown_ignored(self, tmp_dir):
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)
        dag.mark_precomputed("nonexistent")  # should not raise

    def test_mark_precomputed_from_state(self, tmp_dir):
        from scagent.inspector import AnnDataState
        ctx = _make_ctx(tmp_dir)
        dag = AnalysisDAG.from_context(ctx)

        state = AnnDataState(
            x_state="log_normalized",
            raw_counts_location="layers['counts']",
            has_qc_metrics=True,
            has_normalized=True,
            has_hvg=True,
            has_pca=True,
            has_neighbors=False,
            n_cells=1000, n_genes=5000,
        )
        count = dag.mark_precomputed_from_state(state)
        assert count > 0
        assert dag.get_step("load").status == "done"
        assert dag.get_step("normalize").status == "done"
        assert dag.get_step("pca").status == "done"
        # Neighbors not done → clustering still pending
        assert dag.get_step("clustering").status == "pending"


class TestFromUnknownParadigm:
    def test_raises(self, tmp_dir):
        ctx = ExperimentContext(tmp_dir)
        ctx._data["paradigm"] = "nonexistent"
        with pytest.raises(ValueError, match="No DAG definition"):
            AnalysisDAG.from_context(ctx)

    def test_none_paradigm_raises(self, tmp_dir):
        ctx = ExperimentContext(tmp_dir)
        with pytest.raises(ValueError, match="paradigm is not set"):
            AnalysisDAG.from_context(ctx)
