"""Paradigm-aware analysis DAG for scAgent.

Generates an ordered, dependency-checked plan of analysis steps based
on the experiment context.  The agent uses this to suggest next steps,
validate ordering, and track progress.

Usage::

    from scagent.context import ExperimentContext
    from scagent.dag import AnalysisDAG

    ctx = ExperimentContext(Path(".scagent"))
    dag = AnalysisDAG.from_context(ctx)
    print(dag.summary())
    step = dag.next_step()
"""

from __future__ import annotations

import json
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from scagent.context import ExperimentContext


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class DAGStep:
    """A single step in the analysis DAG."""

    id: str
    name: str  # human-readable
    category: str  # qc, normalization, clustering, de, trajectory, etc.
    tool_id: str | None = None  # maps to tool registry; None for meta-steps
    depends_on: list[str] = field(default_factory=list)
    status: str = "pending"  # pending, done, skipped
    skip_reason: str | None = None
    required: bool = True
    conditional: bool = False  # True if step is included only when condition is met


# ---------------------------------------------------------------------------
# DAG definitions per paradigm
# ---------------------------------------------------------------------------

def _cell_atlas_steps(ctx: ExperimentContext) -> list[DAGStep]:
    """Cell atlas: single-tissue deep profiling."""
    steps = [
        DAGStep("load", "Load data", "loading", "load_10x_h5"),
        DAGStep("qc_metrics", "QC metrics", "qc", "calculate_qc_metrics", ["load"]),
        DAGStep("filter_cells", "Filter cells", "qc", "filter_cells", ["qc_metrics"]),
        DAGStep("filter_genes", "Filter genes", "qc", "filter_genes", ["filter_cells"]),
        DAGStep("doublet_detection", "Doublet detection", "qc", "scrublet_doublets", ["filter_genes"]),
        DAGStep("normalize", "Normalize", "normalization", "log_normalize", ["doublet_detection"]),
        DAGStep("hvg", "Highly variable genes", "feature_selection", "highly_variable_genes", ["normalize"]),
        DAGStep("pca", "PCA", "dimensionality_reduction", "pca", ["hvg"]),
    ]

    if ctx.needs_batch_correction():
        steps.append(DAGStep(
            "batch_correction", "Batch correction", "integration", "harmony",
            ["pca"], conditional=True,
        ))
        neighbor_dep = "batch_correction"
    else:
        neighbor_dep = "pca"

    steps.extend([
        DAGStep("neighbors", "Neighbor graph", "neighbors", "neighbor_graph", [neighbor_dep]),
        DAGStep("clustering", "Clustering", "clustering", "leiden_clustering", ["neighbors"]),
        DAGStep("umap", "UMAP embedding", "embedding", "umap", ["neighbors"]),
        DAGStep("markers", "Marker genes", "differential_expression", "wilcoxon_markers", ["clustering"]),
        DAGStep("annotation", "Cell type annotation", "annotation", "celltypist_annotation", ["clustering", "markers"]),
    ])

    return steps


def _disease_vs_healthy_steps(ctx: ExperimentContext) -> list[DAGStep]:
    """Disease vs. healthy: cross-condition comparison."""
    steps = _cell_atlas_steps(ctx)  # shared prefix

    # Add DE + enrichment + composition
    steps.extend([
        DAGStep(
            "pseudobulk_de", "Pseudobulk DE", "differential_expression",
            "deseq2_pseudobulk", ["annotation"], required=True,
        ),
        DAGStep(
            "pathway_enrichment", "Pathway enrichment", "enrichment",
            "gsea", ["pseudobulk_de"], required=False,
        ),
        DAGStep(
            "composition", "Composition analysis", "composition",
            "sccoda", ["annotation"], required=False, conditional=True,
        ),
    ])

    return steps


def _developmental_trajectory_steps(ctx: ExperimentContext) -> list[DAGStep]:
    """Developmental trajectory: pseudotime / lineage analysis.

    Workflow: shared prefix → PAGA topology → DPT pseudotime → (optional) scVelo.
    Best-practice refs: [BP-1] pp. 553-554, [BP-2] Ch. 14.
    """
    steps = _cell_atlas_steps(ctx)  # shared prefix

    # Add trajectory-specific steps
    steps.extend([
        DAGStep(
            "paga", "PAGA trajectory topology", "trajectory",
            "paga", ["clustering"],
        ),
        DAGStep(
            "pseudotime", "Diffusion pseudotime", "trajectory",
            "diffusion_pseudotime", ["paga", "annotation"],
        ),
        DAGStep(
            "rna_velocity", "RNA velocity", "trajectory",
            "scvelo_velocity", ["neighbors"],
            required=False, conditional=True,
        ),
    ])

    return steps


def _perturbation_screen_steps(ctx: ExperimentContext) -> list[DAGStep]:
    """Perturbation screen: CRISPR / Perturb-seq analysis.

    Workflow: shared prefix → guide assignment → perturbation DE → enrichment.
    Best-practice refs: [BP-1] p. 557, [BP-2] Ch. 20.
    """
    steps = _cell_atlas_steps(ctx)  # shared prefix

    steps.extend([
        DAGStep(
            "guide_assignment", "Guide assignment", "perturbation",
            "guide_assignment", ["annotation"],
        ),
        DAGStep(
            "perturbation_de", "Perturbation DE", "differential_expression",
            "perturbation_de", ["guide_assignment"],
        ),
        DAGStep(
            "perturbation_enrichment", "Perturbation enrichment", "enrichment",
            "gsea", ["perturbation_de"], required=False,
        ),
    ])

    return steps


def _temporal_longitudinal_steps(ctx: ExperimentContext) -> list[DAGStep]:
    """Temporal / longitudinal: time-series cross-condition analysis.

    Reuses existing tools: pseudobulk DE for pairwise timepoint contrasts,
    composition for proportion trends, trajectory for temporal ordering.
    Batch correction is critical — different timepoints = different batches.
    Best-practice refs: [BP-1] pp. 552-553, 555, [BP-2] Ch. 18.
    """
    steps = [
        DAGStep("load", "Load data", "loading", "load_10x_h5"),
        DAGStep("qc_metrics", "QC metrics", "qc", "calculate_qc_metrics", ["load"]),
        DAGStep("filter_cells", "Filter cells", "qc", "filter_cells", ["qc_metrics"]),
        DAGStep("filter_genes", "Filter genes", "qc", "filter_genes", ["filter_cells"]),
        DAGStep("doublet_detection", "Doublet detection", "qc", "scrublet_doublets", ["filter_genes"]),
        DAGStep("normalize", "Normalize", "normalization", "log_normalize", ["doublet_detection"]),
        DAGStep("hvg", "Highly variable genes", "feature_selection", "highly_variable_genes", ["normalize"]),
        DAGStep("pca", "PCA", "dimensionality_reduction", "pca", ["hvg"]),
    ]

    # Batch correction is ALWAYS required for temporal data
    # (different timepoints = different batches) [BP-1]
    steps.append(DAGStep(
        "batch_correction", "Batch correction", "integration", "harmony",
        ["pca"], required=True,
    ))

    steps.extend([
        DAGStep("neighbors", "Neighbor graph", "neighbors", "neighbor_graph", ["batch_correction"]),
        DAGStep("clustering", "Clustering", "clustering", "leiden_clustering", ["neighbors"]),
        DAGStep("umap", "UMAP embedding", "embedding", "umap", ["neighbors"]),
        DAGStep("markers", "Marker genes", "differential_expression", "wilcoxon_markers", ["clustering"]),
        DAGStep("annotation", "Cell type annotation", "annotation", "celltypist_annotation", ["clustering", "markers"]),
        # Temporal-specific analysis steps
        DAGStep(
            "pseudobulk_de", "Pseudobulk DE (timepoints)", "differential_expression",
            "deseq2_pseudobulk", ["annotation"], required=True,
        ),
        DAGStep(
            "composition", "Composition over time", "composition",
            "sccoda", ["annotation"], required=False,
        ),
        DAGStep(
            "pathway_enrichment", "Pathway enrichment", "enrichment",
            "gsea", ["pseudobulk_de"], required=False,
        ),
    ])

    return steps


def _immune_repertoire_steps(ctx: ExperimentContext) -> list[DAGStep]:
    """Immune repertoire: V(D)J / TCR/BCR clonotype analysis.

    Workflow: shared GEX prefix → load VDJ → clonotype analysis → overlap.
    Best-practice refs: [BP-1] pp. 559-560, [BP-2] Ch. 38-39.
    """
    steps = _cell_atlas_steps(ctx)  # shared GEX prefix

    steps.extend([
        DAGStep(
            "load_vdj", "Load V(D)J data", "loading",
            "load_vdj", ["load"],
        ),
        DAGStep(
            "clonotype_analysis", "Clonotype analysis", "repertoire",
            "clonotype_analysis", ["load_vdj", "annotation"],
        ),
        DAGStep(
            "repertoire_overlap", "Repertoire overlap", "repertoire",
            "repertoire_overlap", ["clonotype_analysis"],
            required=False,
        ),
    ])

    return steps


def _multimodal_steps(ctx: ExperimentContext) -> list[DAGStep]:
    """Multimodal (CITE-seq): RNA + surface protein analysis.

    Workflow: GEX prefix → load protein → normalize protein → WNN → cluster → annotate.
    Best-practice refs: [BP-1] pp. 558-559, [BP-2] Ch. 32-37.
    """
    steps = [
        DAGStep("load", "Load data", "loading", "load_10x_h5"),
        DAGStep("qc_metrics", "QC metrics", "qc", "calculate_qc_metrics", ["load"]),
        DAGStep("filter_cells", "Filter cells", "qc", "filter_cells", ["qc_metrics"]),
        DAGStep("filter_genes", "Filter genes", "qc", "filter_genes", ["filter_cells"]),
        DAGStep("doublet_detection", "Doublet detection", "qc", "scrublet_doublets", ["filter_genes"]),
        DAGStep("normalize", "Normalize", "normalization", "log_normalize", ["doublet_detection"]),
        DAGStep("hvg", "Highly variable genes", "feature_selection", "highly_variable_genes", ["normalize"]),
        DAGStep("pca", "PCA", "dimensionality_reduction", "pca", ["hvg"]),
    ]

    # Protein modality
    steps.extend([
        DAGStep(
            "load_protein", "Load protein data", "loading",
            "load_protein", ["load"],
        ),
        DAGStep(
            "normalize_protein", "Normalize protein (CLR)", "normalization",
            "normalize_protein", ["load_protein"],
        ),
        DAGStep(
            "wnn", "Weighted nearest neighbors", "neighbors",
            "wnn", ["pca", "normalize_protein"],
        ),
        DAGStep(
            "clustering", "Clustering (WNN)", "clustering",
            "leiden_clustering", ["wnn"],
        ),
        DAGStep(
            "umap", "UMAP embedding", "embedding",
            "umap", ["wnn"],
        ),
        DAGStep(
            "markers", "Marker genes", "differential_expression",
            "wilcoxon_markers", ["clustering"],
        ),
        DAGStep(
            "protein_markers", "Protein markers", "differential_expression",
            "protein_markers", ["clustering", "normalize_protein"],
        ),
        DAGStep(
            "annotation", "Cell type annotation", "annotation",
            "celltypist_annotation", ["clustering", "markers", "protein_markers"],
        ),
    ])

    return steps


_PARADIGM_BUILDERS = {
    "cell_atlas": _cell_atlas_steps,
    "disease_vs_healthy": _disease_vs_healthy_steps,
    "developmental_trajectory": _developmental_trajectory_steps,
    "perturbation_screen": _perturbation_screen_steps,
    "temporal_longitudinal": _temporal_longitudinal_steps,
    "immune_repertoire": _immune_repertoire_steps,
    "multimodal": _multimodal_steps,
}


# ---------------------------------------------------------------------------
# AnalysisDAG
# ---------------------------------------------------------------------------

class AnalysisDAG:
    """Paradigm-aware analysis DAG with progress tracking.

    Parameters
    ----------
    paradigm
        Experiment paradigm (e.g., ``"cell_atlas"``).
    steps
        Ordered list of :class:`DAGStep` objects.
    """

    FILENAME = "dag.json"

    def __init__(self, paradigm: str, steps: list[DAGStep]) -> None:
        self.paradigm = paradigm
        self.steps = steps
        self._step_index: dict[str, DAGStep] = {s.id: s for s in steps}

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def from_context(cls, ctx: ExperimentContext) -> "AnalysisDAG":
        """Generate a DAG appropriate for the experiment context."""
        paradigm = ctx.paradigm
        if paradigm is None:
            raise ValueError("Cannot generate DAG: paradigm is not set in experiment context")
        builder = _PARADIGM_BUILDERS.get(paradigm)
        if builder is None:
            raise ValueError(
                f"No DAG definition for paradigm '{paradigm}'. "
                f"Supported: {sorted(_PARADIGM_BUILDERS.keys())}"
            )
        steps = builder(ctx)
        return cls(paradigm, steps)

    # ------------------------------------------------------------------
    # Navigation
    # ------------------------------------------------------------------

    def next_step(self) -> DAGStep | None:
        """Return the first pending step whose dependencies are all met."""
        for step in self.steps:
            if step.status != "pending":
                continue
            if self._deps_met(step):
                return step
        return None

    def complete_step(self, step_id: str) -> None:
        """Mark a step as done."""
        step = self._get(step_id)
        step.status = "done"

    def skip_step(self, step_id: str, reason: str = "") -> None:
        """Mark a step as skipped."""
        step = self._get(step_id)
        step.status = "skipped"
        step.skip_reason = reason

    def is_valid_step(self, step_id: str) -> bool:
        """Check if a step can be executed (all dependencies met)."""
        step = self._get(step_id)
        return self._deps_met(step)

    def get_step(self, step_id: str) -> DAGStep | None:
        """Return a step by ID, or *None*."""
        return self._step_index.get(step_id)

    # ------------------------------------------------------------------
    # Dynamic modification
    # ------------------------------------------------------------------

    def add_step(
        self,
        step: DAGStep,
        after: str | None = None,
    ) -> None:
        """Add a step to the DAG dynamically.

        Parameters
        ----------
        step
            The new step to add.
        after
            Insert after this step ID. If *None*, appends at the end.

        Raises
        ------
        ValueError
            If *step.id* already exists or *after* is not found.
        """
        if step.id in self._step_index:
            raise ValueError(f"Step '{step.id}' already exists in DAG")

        if after is not None:
            if after not in self._step_index:
                raise ValueError(f"Step '{after}' not found in DAG")
            idx = next(i for i, s in enumerate(self.steps) if s.id == after)
            self.steps.insert(idx + 1, step)
        else:
            self.steps.append(step)

        self._step_index[step.id] = step

    def mark_precomputed(self, step_id: str) -> None:
        """Mark a step as already done (e.g., detected by inspector).

        Silently ignores unknown step IDs so callers can mark steps
        that may or may not exist in this particular DAG.
        """
        step = self._step_index.get(step_id)
        if step is not None:
            step.status = "done"

    def mark_precomputed_from_state(self, state) -> int:
        """Mark DAG steps as done based on an :class:`AnnDataState`.

        Parameters
        ----------
        state
            Inspector output (:class:`~scagent.inspector.AnnDataState`).

        Returns
        -------
        Number of steps marked as precomputed.
        """
        # Map state flags → DAG step IDs
        _FLAG_TO_STEPS = {
            "has_qc_metrics": ["load", "qc_metrics", "filter_cells", "filter_genes"],
            "has_normalized": ["normalize"],
            "has_hvg": ["hvg"],
            "has_pca": ["pca"],
            "has_neighbors": ["neighbors"],
            "has_clusters": ["clustering"],
            "has_umap": ["umap"],
            "has_cell_types": ["annotation"],
            "has_de_results": ["markers"],
        }
        count = 0
        for flag, step_ids in _FLAG_TO_STEPS.items():
            if getattr(state, flag, False):
                for sid in step_ids:
                    step = self._step_index.get(sid)
                    if step is not None and step.status == "pending":
                        step.status = "done"
                        count += 1
        return count

    @property
    def done_steps(self) -> list[DAGStep]:
        return [s for s in self.steps if s.status == "done"]

    @property
    def pending_steps(self) -> list[DAGStep]:
        return [s for s in self.steps if s.status == "pending"]

    @property
    def progress(self) -> tuple[int, int]:
        """Return (done_count, total_count)."""
        total = sum(1 for s in self.steps if s.status != "skipped")
        done = sum(1 for s in self.steps if s.status == "done")
        return done, total

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Human-readable Markdown progress table."""
        done, total = self.progress
        lines = [
            f"## Analysis Plan — {self.paradigm} ({done}/{total} steps done)\n",
            "| # | Step | Category | Status |",
            "|---|------|----------|--------|",
        ]
        for i, s in enumerate(self.steps):
            status_icon = {"done": "✅", "skipped": "⏭️", "pending": "⬜"}[s.status]
            extra = ""
            if s.conditional:
                extra = " _(conditional)_"
            if s.skip_reason:
                extra += f" — {s.skip_reason}"
            lines.append(f"| {i} | {s.name} | {s.category} | {status_icon} {s.status}{extra} |")

        nxt = self.next_step()
        if nxt:
            lines.append(f"\n**Next step:** {nxt.name} (`{nxt.id}`)")

        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, project_dir: Path | str) -> Path:
        """Write the DAG to ``.scagent/dag.json``."""
        p = Path(project_dir)
        p.mkdir(parents=True, exist_ok=True)
        path = p / self.FILENAME
        data = {
            "paradigm": self.paradigm,
            "steps": [asdict(s) for s in self.steps],
        }
        path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        return path

    @classmethod
    def load(cls, project_dir: Path | str) -> "AnalysisDAG":
        """Load a DAG from ``.scagent/dag.json``."""
        p = Path(project_dir) / cls.FILENAME
        data = json.loads(p.read_text(encoding="utf-8"))
        steps = [DAGStep(**s) for s in data["steps"]]
        return cls(data["paradigm"], steps)

    # ------------------------------------------------------------------
    # Private
    # ------------------------------------------------------------------

    def _get(self, step_id: str) -> DAGStep:
        step = self._step_index.get(step_id)
        if step is None:
            raise ValueError(f"No step with id '{step_id}' in DAG")
        return step

    def _deps_met(self, step: DAGStep) -> bool:
        """Check if all dependencies are done or skipped."""
        for dep_id in step.depends_on:
            dep = self._step_index.get(dep_id)
            if dep is None:
                return False
            if dep.status not in ("done", "skipped"):
                return False
        return True
