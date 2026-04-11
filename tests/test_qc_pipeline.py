#!/usr/bin/env python3
"""End-to-end smoke test of the preprocessing pipeline on PBMC 10k data.

Exercises: load → QC metrics → filter cells → filter genes → doublets → normalize.
Verifies expected dimensions, reproducibility, and checkpoint/resume.
"""

import shutil
import tempfile
from pathlib import Path

import numpy as np

from scagent.tools.loading import load_10x_h5
from scagent.tools.qc import calculate_qc_metrics, filter_cells, filter_genes
from scagent.tools.doublets import detect_doublets
from scagent.tools.normalize import log_normalize
from scagent.tools.pipeline import run_preprocessing

DATA_FILE = "data/pbmc10k/filtered_feature_bc_matrix.h5"


def test_load():
    print("=== test_load ===")
    adata, info = load_10x_h5(DATA_FILE)
    assert adata.n_obs == 11769, f"Expected 11769 cells, got {adata.n_obs}"
    assert adata.n_vars == 33538, f"Expected 33538 genes, got {adata.n_vars}"
    assert info["metrics"]["species"] == "human"
    assert info["metrics"]["duplicate_var_names_fixed"] > 0
    print(f"  ✓ {adata.n_obs} cells × {adata.n_vars} genes, species={info['metrics']['species']}")
    return adata


def test_qc_metrics(adata):
    print("=== test_qc_metrics ===")
    result = calculate_qc_metrics(adata)
    m = result["metrics"]

    assert "pct_counts_mt" in adata.obs.columns
    assert "n_genes_by_counts" in adata.obs.columns
    assert m["n_mito_genes"] > 0, "No mitochondrial genes found"
    assert 0 < m["median_pct_mito"] < 50, f"Unexpected median mito: {m['median_pct_mito']}"
    assert "recommended_thresholds" in m

    print(f"  ✓ Median genes={m['median_genes_per_cell']:.0f}, "
          f"counts={m['median_counts_per_cell']:.0f}, mito={m['median_pct_mito']:.1f}%")
    print(f"  ✓ MAD thresholds: {m['recommended_thresholds']}")


def test_filter_cells(adata):
    print("=== test_filter_cells ===")
    n_before = adata.n_obs
    adata_f, result = filter_cells(adata, min_genes=200, max_genes=5000, max_pct_mito=20.0)
    m = result["metrics"]

    assert m["cells_before"] == n_before
    assert m["cells_after"] < n_before, "No cells were filtered"
    assert m["cells_after"] > n_before * 0.5, "More than 50% cells removed"
    assert adata_f.n_obs == m["cells_after"]

    print(f"  ✓ {m['cells_before']} → {m['cells_after']} cells ({m['pct_removed']}% removed)")
    return adata_f


def test_filter_genes(adata):
    print("=== test_filter_genes ===")
    n_before = adata.n_vars
    adata_f, result = filter_genes(adata, min_cells=3)
    m = result["metrics"]

    assert m["genes_before"] == n_before
    assert m["genes_after"] < n_before
    assert m["genes_after"] > n_before * 0.3

    print(f"  ✓ {m['genes_before']} → {m['genes_after']} genes ({m['genes_removed']} removed)")
    return adata_f


def test_doublets(adata):
    print("=== test_doublets ===")
    n_before = adata.n_obs
    adata_d, result = detect_doublets(adata, random_state=0)
    m = result["metrics"]

    assert 0.01 < m["doublet_rate"] < 0.15, f"Unusual doublet rate: {m['doublet_rate']}"
    assert m["cells_after"] < n_before
    assert m["cells_after"] > n_before * 0.8  # shouldn't lose >20% to doublets

    print(f"  ✓ {m['n_doublets']} doublets ({m['doublet_rate']:.1%}), "
          f"expected ~{m['expected_doublet_rate']:.1%}")
    return adata_d


def test_normalize(adata):
    print("=== test_normalize ===")
    result = log_normalize(adata)
    m = result["metrics"]

    assert m["max_normalized_value"] < 15, f"Max value too high: {m['max_normalized_value']}"
    assert not m["has_nan"], "NaN values after normalization"
    assert m["raw_counts_preserved"], "Raw counts not preserved"
    assert m["raw_frozen"], "adata.raw not set"
    assert m["counts_are_integers"], "Raw counts are not integers"
    assert adata.raw is not None
    assert "counts" in adata.layers

    print(f"  ✓ Max value={m['max_normalized_value']:.1f}, "
          f"counts_layer={'counts' in adata.layers}, raw={adata.raw is not None}")


def test_reproducibility():
    print("=== test_reproducibility ===")
    # Run doublet detection twice with same seed
    adata1, _ = load_10x_h5(DATA_FILE)
    calculate_qc_metrics(adata1)
    adata1, _ = filter_cells(adata1, min_genes=200, max_genes=5000, max_pct_mito=20.0)
    adata1, _ = filter_genes(adata1, min_cells=3)
    adata1, r1 = detect_doublets(adata1, random_state=42)

    adata2, _ = load_10x_h5(DATA_FILE)
    calculate_qc_metrics(adata2)
    adata2, _ = filter_cells(adata2, min_genes=200, max_genes=5000, max_pct_mito=20.0)
    adata2, _ = filter_genes(adata2, min_cells=3)
    adata2, r2 = detect_doublets(adata2, random_state=42)

    assert r1["metrics"]["n_doublets"] == r2["metrics"]["n_doublets"], "Same seed produced different results!"
    print(f"  ✓ Same seed=42 → same {r1['metrics']['n_doublets']} doublets both times")


def test_pipeline_full():
    print("=== test_pipeline_full ===")
    tmpdir = tempfile.mkdtemp(prefix="scagent_test_")
    try:
        adata, results = run_preprocessing(
            DATA_FILE,
            max_pct_mito=20.0,
            random_state=0,
            checkpoint_dir=f"{tmpdir}/checkpoints",
            plot_dir=f"{tmpdir}/plots",
        )

        # Check all steps ran
        assert "load" in results
        assert "qc_metrics" in results
        assert "filter_cells" in results
        assert "filter_genes" in results
        assert "doublets" in results
        assert "normalize" in results

        # Check final state
        assert adata.raw is not None
        assert "counts" in adata.layers
        assert adata.n_obs > 5000  # shouldn't lose too many cells overall

        # Check checkpoints exist
        cp_dir = Path(tmpdir) / "checkpoints"
        checkpoints = sorted(p.stem for p in cp_dir.glob("*.h5ad"))
        assert "after_load" in checkpoints
        assert "after_normalize" in checkpoints

        print(f"  ✓ Full pipeline: {results['load']['metrics']['n_cells']} → {adata.n_obs} cells")
        print(f"  ✓ Checkpoints: {checkpoints}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def test_pipeline_resume():
    print("=== test_pipeline_resume ===")
    tmpdir = tempfile.mkdtemp(prefix="scagent_test_")
    try:
        # Run full pipeline first to create checkpoints
        run_preprocessing(
            DATA_FILE,
            max_pct_mito=20.0,
            random_state=0,
            checkpoint_dir=f"{tmpdir}/checkpoints",
        )

        # Resume from after_filter_cells with different mito threshold
        # This simulates "go back to after cell filtering and try 15% mito"
        adata2, results2 = run_preprocessing(
            DATA_FILE,  # ignored since start_from is set
            max_pct_mito=15.0,  # different threshold (but won't re-filter since starting after)
            random_state=0,
            checkpoint_dir=f"{tmpdir}/checkpoints",
            start_from="after_filter_cells",
        )

        # Should only have results for steps after filter_cells
        assert "load" not in results2
        assert "filter_cells" not in results2
        assert "filter_genes" in results2
        assert "doublets" in results2
        assert "normalize" in results2
        assert adata2.raw is not None

        print(f"  ✓ Resumed from after_filter_cells: {adata2.n_obs} cells")
        print(f"  ✓ Steps run: {list(results2.keys())}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    adata = test_load()
    test_qc_metrics(adata)
    adata = test_filter_cells(adata)
    adata = test_filter_genes(adata)
    adata = test_doublets(adata)
    test_normalize(adata)
    test_reproducibility()
    test_pipeline_full()
    test_pipeline_resume()
    print("\n✅ All tests passed.")
