#!/usr/bin/env python3
"""Integration test: init flow on real PBMC 10k dataset.

1. Load data, infer context
2. Set paradigm + tissue, validate
3. Generate DAG, verify steps
4. Walk through 3 steps, track progress
5. Verify project.json + dag.json on disk

Requires: data/pbmc10k/filtered_feature_bc_matrix.h5
"""

import json
import tempfile
from pathlib import Path

import pytest

from scagent.context import ExperimentContext
from scagent.dag import AnalysisDAG

DATA_FILE = Path("data/pbmc10k/filtered_feature_bc_matrix.h5")

pytestmark = pytest.mark.skipif(
    not DATA_FILE.exists(),
    reason=f"Dataset not found: {DATA_FILE}",
)


def test_init_flow():
    from scagent.tools.loading import load_10x_h5

    with tempfile.TemporaryDirectory() as tmp:
        proj = Path(tmp) / ".scagent"

        # === Step 1: Load data and infer ===
        adata, _ = load_10x_h5(str(DATA_FILE))
        inferred = ExperimentContext.infer_from_data(adata)

        assert inferred["organism"]["species"] == "Homo sapiens"
        assert inferred["organism"]["inferred"] is True
        assert inferred["_data_shape"]["n_cells"] == 11769
        assert inferred["_data_shape"]["n_genes"] == 33538
        assert inferred["library"]["umi"] is True
        print(f"\n  Inferred: {inferred['organism']['species']}, "
              f"{inferred['_data_shape']['n_cells']} cells, "
              f"UMI={inferred['library']['umi']}")

        # === Step 2: Create context ===
        ctx = ExperimentContext(proj)
        ctx.paradigm = "cell_atlas"
        ctx.organism = inferred["organism"]
        ctx.tissue = {"name": "peripheral blood mononuclear cells"}
        ctx.platform = inferred.get("platform", {})
        ctx.library = inferred.get("library", {"type": "3prime", "umi": True})
        ctx.samples = [{"id": "pbmc10k", "file": str(DATA_FILE)}]

        errors = ctx.validate()
        assert errors == [], f"Validation errors: {errors}"
        assert ctx.is_complete()
        ctx.save()

        # Verify project.json on disk
        pj = proj / "project.json"
        assert pj.exists()
        data = json.loads(pj.read_text())
        assert data["paradigm"] == "cell_atlas"
        assert data["organism"]["species"] == "Homo sapiens"

        # === Step 3: Generate DAG ===
        dag = AnalysisDAG.from_context(ctx)
        dag.save(proj)

        step_ids = [s.id for s in dag.steps]
        assert "load" in step_ids
        assert "normalize" in step_ids
        assert "clustering" in step_ids
        assert "annotation" in step_ids
        # Cell atlas should NOT have pseudobulk DE
        assert "pseudobulk_de" not in step_ids
        # Single sample — no batch correction
        assert "batch_correction" not in step_ids

        print(f"  DAG: {len(dag.steps)} steps for {dag.paradigm}")
        print(f"  Steps: {' → '.join(step_ids)}")

        # Verify dag.json on disk
        dj = proj / "dag.json"
        assert dj.exists()

        # === Step 4: Walk through 3 steps ===
        assert dag.next_step().id == "load"
        dag.complete_step("load")
        assert dag.next_step().id == "qc_metrics"
        dag.complete_step("qc_metrics")
        assert dag.next_step().id == "filter_cells"
        dag.complete_step("filter_cells")

        done, total = dag.progress
        assert done == 3
        print(f"  Progress: {done}/{total} steps done")

        # === Step 5: Summary ===
        md = dag.summary()
        assert "✅" in md
        assert "3/" in md
        print(f"  Summary:\n{md}")

        # === Step 6: Save updated state, reload from disk ===
        dag.save(proj)
        dag3 = AnalysisDAG.load(proj)
        assert dag3.get_step("load").status == "done"
        assert dag3.get_step("filter_cells").status == "done"
        assert dag3.get_step("normalize").status == "pending"
        print("  ✓ Roundtrip: DAG state persists across save/load")
