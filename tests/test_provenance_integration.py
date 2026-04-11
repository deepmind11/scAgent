#!/usr/bin/env python3
"""Integration test: run real PBMC pipeline and verify provenance recording.

Requires: data/pbmc10k/filtered_feature_bc_matrix.h5
"""

import json
import tempfile
from pathlib import Path

import pytest

from scagent.provenance import ProvenanceGraph, record_step

DATA_FILE = Path("data/pbmc10k/filtered_feature_bc_matrix.h5")

pytestmark = pytest.mark.skipif(
    not DATA_FILE.exists(),
    reason=f"Dataset not found: {DATA_FILE}",
)


def test_full_pipeline_provenance():
    """Run load → qc → norm → hvg → pca → neighbors → leiden, record each step."""
    from scagent.tools.loading import load_10x_h5
    from scagent.tools.qc import calculate_qc_metrics, filter_cells, filter_genes
    from scagent.tools.normalize import log_normalize
    from scagent.tools.feature_selection import select_hvg
    from scagent.tools.pca import run_pca
    from scagent.tools.neighbors import compute_neighbors
    from scagent.tools.clustering import run_leiden

    with tempfile.TemporaryDirectory() as tmp:
        proj = Path(tmp) / ".scagent"
        graph = ProvenanceGraph(proj)

        # Step 0: Load
        adata, r0 = load_10x_h5(str(DATA_FILE))
        record_step(graph, r0, user_prompt="load the PBMC 10k dataset")

        # Step 1: QC metrics
        r1 = calculate_qc_metrics(adata)
        record_step(graph, r1, user_prompt="compute QC metrics")

        # Step 2: Filter cells
        r2 = filter_cells(adata, min_genes=200, max_genes=5000, max_pct_mito=20.0)
        record_step(graph, r2, user_prompt="filter cells")

        # Step 3: Filter genes
        r3 = filter_genes(adata, min_cells=3)
        record_step(graph, r3, user_prompt="filter genes")

        # Step 4: Normalize
        r4 = log_normalize(adata)
        record_step(graph, r4, user_prompt="log-normalize")

        # Step 5: HVG
        r5 = select_hvg(adata, n_top_genes=2000)
        record_step(graph, r5, user_prompt="select highly variable genes")

        # Step 6: PCA
        r6 = run_pca(adata, n_comps=50)
        record_step(graph, r6, user_prompt="run PCA")

        # Step 7: Neighbors
        r7 = compute_neighbors(adata, n_neighbors=15, n_pcs=30)
        record_step(graph, r7, user_prompt="compute neighbor graph")

        # Step 8: Leiden
        r8 = run_leiden(adata, resolution=1.0)
        record_step(graph, r8, user_prompt="cluster at resolution 1.0")

        # -- Verify --

        # 1. File exists and is valid JSON-LD
        prov_path = proj / "provenance.jsonld"
        assert prov_path.exists()
        doc = json.loads(prov_path.read_text())
        assert "@context" in doc
        assert "@graph" in doc

        # 2. Correct number of activities
        assert graph.n_activities == 9

        # 3. Chain is ordered
        chain = graph.get_chain()
        expected_tools = [
            "load_10x_h5",
            "calculate_qc_metrics",
            "filter_cells",
            "filter_genes",
            "log_normalize",
            "highly_variable_genes",
            "pca",
            "neighbor_graph",
            "leiden_clustering",
        ]
        actual_tools = [a["tool_id"] for a in chain]
        assert actual_tools == expected_tools

        # 4. Entity chain is linked
        output_entities = [e for e in graph._entities if e.step_index >= 0]
        for i in range(1, len(output_entities)):
            assert output_entities[i].derived_from is not None

        # 5. Parameters match what we actually passed
        filter_act = [a for a in chain if a["tool_id"] == "filter_cells"][0]
        assert filter_act["parameters"]["min_genes"] == 200
        assert filter_act["parameters"]["max_pct_mito"] == 20.0

        pca_act = [a for a in chain if a["tool_id"] == "pca"][0]
        assert pca_act["parameters"]["n_comps"] == 50

        leiden_act = [a for a in chain if a["tool_id"] == "leiden_clustering"][0]
        assert leiden_act["parameters"]["resolution"] == 1.0

        # 6. Software versions captured
        last_activity = graph._activities[-1]
        assert "scanpy" in last_activity.software_versions
        assert "python" in last_activity.software_versions
        assert "scagent" in last_activity.software_versions

        # 7. Replay plan matches call sequence
        plan = graph.replay_plan()
        assert len(plan) == 9
        assert [t for t, _ in plan] == expected_tools

        # 8. Summary is readable
        md = graph.summary()
        assert "## Analysis Provenance" in md
        assert "9 steps" in md

        # 9. File size within budget
        size = prov_path.stat().st_size
        assert size < 50_000, f"provenance.jsonld is {size:,} bytes"
        print(f"\n  Integration test: 9-step provenance = {size:,} bytes")
        print(f"  Summary:\n{md}")
