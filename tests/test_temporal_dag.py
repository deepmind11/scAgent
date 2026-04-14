"""Tests for temporal/longitudinal DAG and perturbation screen DAG."""

from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from scagent.dag import (
    AnalysisDAG,
    _perturbation_screen_steps,
    _temporal_longitudinal_steps,
    _PARADIGM_BUILDERS,
)


# ---------------------------------------------------------------------------
# Temporal DAG tests
# ---------------------------------------------------------------------------


class TestTemporalDAG:
    @pytest.fixture
    def ctx(self):
        ctx = MagicMock()
        ctx.paradigm = "temporal_longitudinal"
        ctx.needs_batch_correction.return_value = True
        return ctx

    def test_temporal_dag_has_required_steps(self, ctx):
        steps = _temporal_longitudinal_steps(ctx)
        step_ids = [s.id for s in steps]

        # Must have standard prefix
        for required in ["load", "qc_metrics", "filter_cells", "normalize", "pca"]:
            assert required in step_ids, f"Missing required step: {required}"

        # Must have batch correction (always required for temporal)
        assert "batch_correction" in step_ids

        # Must have temporal-specific steps
        assert "pseudobulk_de" in step_ids
        assert "composition" in step_ids

    def test_batch_correction_always_required(self, ctx):
        steps = _temporal_longitudinal_steps(ctx)
        bc_step = next(s for s in steps if s.id == "batch_correction")
        assert bc_step.required is True

    def test_batch_correction_before_neighbors(self, ctx):
        steps = _temporal_longitudinal_steps(ctx)
        step_ids = [s.id for s in steps]
        assert step_ids.index("batch_correction") < step_ids.index("neighbors")

    def test_pseudobulk_de_after_annotation(self, ctx):
        steps = _temporal_longitudinal_steps(ctx)
        de_step = next(s for s in steps if s.id == "pseudobulk_de")
        assert "annotation" in de_step.depends_on

    def test_composition_after_annotation(self, ctx):
        steps = _temporal_longitudinal_steps(ctx)
        comp_step = next(s for s in steps if s.id == "composition")
        assert "annotation" in comp_step.depends_on
        # Composition is optional
        assert comp_step.required is False

    def test_temporal_dag_from_context(self, ctx):
        dag = AnalysisDAG.from_context(ctx)
        assert dag.paradigm == "temporal_longitudinal"
        assert len(dag.steps) > 10

    def test_temporal_dag_summary(self, ctx):
        dag = AnalysisDAG.from_context(ctx)
        summary = dag.summary()
        assert "temporal_longitudinal" in summary
        assert "Pseudobulk DE" in summary


# ---------------------------------------------------------------------------
# Perturbation DAG tests
# ---------------------------------------------------------------------------


class TestPerturbationDAG:
    @pytest.fixture
    def ctx(self):
        ctx = MagicMock()
        ctx.paradigm = "perturbation_screen"
        ctx.needs_batch_correction.return_value = False
        return ctx

    def test_perturbation_dag_has_required_steps(self, ctx):
        steps = _perturbation_screen_steps(ctx)
        step_ids = [s.id for s in steps]

        assert "guide_assignment" in step_ids
        assert "perturbation_de" in step_ids

    def test_guide_assignment_after_annotation(self, ctx):
        steps = _perturbation_screen_steps(ctx)
        ga_step = next(s for s in steps if s.id == "guide_assignment")
        assert "annotation" in ga_step.depends_on

    def test_perturbation_de_after_guide_assignment(self, ctx):
        steps = _perturbation_screen_steps(ctx)
        de_step = next(s for s in steps if s.id == "perturbation_de")
        assert "guide_assignment" in de_step.depends_on

    def test_perturbation_dag_from_context(self, ctx):
        dag = AnalysisDAG.from_context(ctx)
        assert dag.paradigm == "perturbation_screen"


# ---------------------------------------------------------------------------
# Full DAG coverage test
# ---------------------------------------------------------------------------


class TestDAGCoverage:
    """Every paradigm in VALID_PARADIGMS must have a DAG builder."""

    def test_all_paradigms_have_builders(self):
        from scagent.context import VALID_PARADIGMS

        for paradigm in VALID_PARADIGMS:
            assert paradigm in _PARADIGM_BUILDERS, (
                f"Paradigm '{paradigm}' is in VALID_PARADIGMS but has no DAG builder. "
                f"Available builders: {sorted(_PARADIGM_BUILDERS.keys())}"
            )

    def test_all_builders_produce_valid_dags(self):
        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        for paradigm, builder in _PARADIGM_BUILDERS.items():
            ctx.paradigm = paradigm
            steps = builder(ctx)
            assert len(steps) > 0, f"Builder for '{paradigm}' returned empty steps"

            # Every step should have an id and name
            for step in steps:
                assert step.id, f"Step in '{paradigm}' has empty id"
                assert step.name, f"Step '{step.id}' in '{paradigm}' has empty name"

            # Check no duplicate step IDs
            ids = [s.id for s in steps]
            assert len(ids) == len(set(ids)), (
                f"Duplicate step IDs in '{paradigm}': "
                f"{[x for x in ids if ids.count(x) > 1]}"
            )

    def test_all_dags_have_load_step(self):
        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        for paradigm, builder in _PARADIGM_BUILDERS.items():
            ctx.paradigm = paradigm
            steps = builder(ctx)
            step_ids = [s.id for s in steps]
            assert "load" in step_ids, (
                f"Paradigm '{paradigm}' DAG is missing 'load' step"
            )

    def test_all_dags_dependencies_exist(self):
        """Every step's depends_on references must be valid step IDs."""
        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        for paradigm, builder in _PARADIGM_BUILDERS.items():
            ctx.paradigm = paradigm
            steps = builder(ctx)
            step_ids = {s.id for s in steps}

            for step in steps:
                for dep in step.depends_on:
                    assert dep in step_ids, (
                        f"Step '{step.id}' in '{paradigm}' depends on "
                        f"'{dep}' which is not in the DAG"
                    )
