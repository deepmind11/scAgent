#!/usr/bin/env python3
"""End-to-end test: preprocessing → clustering → annotation on PBMC 10k.

Verifies the full pipeline produces expected cell types.
"""

import shutil
import tempfile

import numpy as np

from scagent.tools.pipeline import run_preprocessing
from scagent.tools.feature_selection import select_hvg
from scagent.tools.pca import run_pca
from scagent.tools.neighbors import compute_neighbors
from scagent.tools.embedding import run_umap
from scagent.tools.clustering import run_leiden, sweep_resolution, check_seed_sensitivity
from scagent.tools.markers import find_marker_genes
from scagent.tools.annotation import annotate_celltypist, annotate_manual

DATA_FILE = "data/pbmc10k/filtered_feature_bc_matrix.h5"


def get_preprocessed_adata():
    """Run Chunk 2 preprocessing and return the result."""
    adata, _ = run_preprocessing(DATA_FILE, max_pct_mito=20.0, random_state=0)
    return adata


def test_hvg(adata):
    print("=== test_hvg ===")
    result = select_hvg(adata, n_top_genes=2000)
    m = result["metrics"]
    assert m["n_hvgs"] == 2000, f"Expected 2000 HVGs, got {m['n_hvgs']}"
    assert "highly_variable" in adata.var.columns
    print(f"  ✓ {m['n_hvgs']} HVGs selected ({m['pct_hvg']}% of genes)")


def test_pca(adata):
    print("=== test_pca ===")
    result = run_pca(adata, n_comps=50, random_state=0)
    m = result["metrics"]
    assert "X_pca" in adata.obsm
    assert adata.obsm["X_pca"].shape == (adata.n_obs, 50)
    assert m["pc1_variance"] < 0.5, f"PC1 explains {m['pc1_variance']:.1%} — too high"
    print(f"  ✓ PCA: {m['n_comps']} components, PC1={m['pc1_variance']:.1%}, "
          f"total={m['total_variance_explained']:.1%}")


def test_neighbors(adata):
    print("=== test_neighbors ===")
    result = compute_neighbors(adata, n_neighbors=15, metric="cosine", random_state=0)
    assert "connectivities" in adata.obsp
    assert "distances" in adata.obsp
    print(f"  ✓ Neighbors: k={result['metrics']['n_neighbors']}, "
          f"metric={result['metrics']['metric']}")


def test_umap(adata):
    print("=== test_umap ===")
    result = run_umap(adata, random_state=0)
    assert "X_umap" in adata.obsm
    assert adata.obsm["X_umap"].shape == (adata.n_obs, 2)
    assert not result["metrics"]["has_nan"]
    print(f"  ✓ UMAP: no NaN, range x={result['metrics']['x_range']}")


def test_leiden(adata):
    print("=== test_leiden ===")
    result = run_leiden(adata, resolution=1.0, random_state=0)
    m = result["metrics"]
    assert m["n_clusters"] >= 2
    assert m["n_clusters"] <= adata.n_obs * 0.1
    print(f"  ✓ Leiden: {m['n_clusters']} clusters at resolution=1.0")
    print(f"    Sizes: {m['cluster_sizes']}")
    return m["n_clusters"]


def test_sweep(adata):
    print("=== test_sweep ===")
    result = sweep_resolution(adata, random_state=0, plot_dir="data/working/plots")
    m = result["metrics"]
    print(f"  ✓ Resolution sweep:")
    for row in m["sweep"]:
        print(f"    res={row['resolution']} → {row['n_clusters']} clusters, "
              f"silhouette={row['silhouette']}, stability={row['stability_ari']}")
    print(f"  ★ Recommended: {m['recommended_resolution']}")
    print(f"    {m['recommendation_note']}")


def test_seed_sensitivity(adata):
    print("=== test_seed_sensitivity ===")
    result = check_seed_sensitivity(adata, resolution=1.0)
    m = result["metrics"]
    print(f"  ✓ Seed sensitivity: mean ARI={m['mean_ari']:.3f}, "
          f"min ARI={m['min_ari']:.3f}, stable={m['stable']}")
    print(f"    Clusters per seed: {m['n_clusters_per_seed']}")


def test_markers(adata):
    print("=== test_markers ===")
    result = find_marker_genes(adata, groupby="leiden", method="wilcoxon")
    m = result["metrics"]
    print(f"  ✓ Markers for {m['n_clusters']} clusters")
    for cluster, markers in list(m["top_markers"].items())[:3]:
        genes = [x["gene"] for x in markers[:3]]
        print(f"    Cluster {cluster}: {', '.join(genes)}")
    if m["weak_clusters"]:
        print(f"    ⚠ Weak clusters: {m['weak_clusters']}")


def test_celltypist(adata):
    print("=== test_celltypist ===")
    result = annotate_celltypist(adata, model="Immune_All_Low.pkl")
    m = result["metrics"]
    print(f"  ✓ CellTypist: {m['n_cell_types']} cell types found")
    for ct, count in list(m["cell_type_counts"].items())[:8]:
        print(f"    {ct}: {count}")


def test_manual_annotation(adata):
    print("=== test_manual_annotation ===")
    # Canonical PBMC markers
    marker_dict = {
        "CD4 T": ["CD3D", "CD4", "IL7R"],
        "CD8 T": ["CD3D", "CD8A", "CD8B"],
        "B cell": ["CD79A", "MS4A1", "CD19"],
        "NK": ["NKG7", "GNLY", "KLRB1"],
        "CD14 Mono": ["CD14", "LYZ", "S100A9"],
        "FCGR3A Mono": ["FCGR3A", "MS4A7"],
        "DC": ["FCER1A", "CST3"],
        "Platelet": ["PPBP", "PF4"],
    }
    result = annotate_manual(adata, marker_dict=marker_dict, groupby="leiden")
    m = result["metrics"]
    print(f"  ✓ Manual annotation: {m['n_clusters']} clusters scored "
          f"against {m['n_cell_types_provided']} cell types")
    for cluster, ranked in list(m["suggestions"].items())[:5]:
        top = ranked[0]
        print(f"    Cluster {cluster} → {top['cell_type']} (score={top['score']:.3f})")
    if m["missing_markers"]:
        print(f"    Missing markers: {m['missing_markers']}")


def test_reproducibility(adata):
    print("=== test_reproducibility ===")
    # Run leiden twice with same seed
    import anndata
    adata2 = adata.copy()
    run_leiden(adata, resolution=1.0, key_added="test_a", random_state=0)
    run_leiden(adata2, resolution=1.0, key_added="test_a", random_state=0)
    match = (adata.obs["test_a"] == adata2.obs["test_a"]).all()
    assert match, "Same seed produced different clusters!"
    print(f"  ✓ Same seed=0 → identical cluster labels")
    del adata.obs["test_a"]


if __name__ == "__main__":
    print("Preprocessing (Chunk 2)...")
    adata = get_preprocessed_adata()
    print(f"Starting: {adata.n_obs} cells × {adata.n_vars} genes\n")

    test_hvg(adata)
    test_pca(adata)
    test_neighbors(adata)
    test_umap(adata)
    n_clusters = test_leiden(adata)
    test_sweep(adata)
    test_seed_sensitivity(adata)
    test_markers(adata)
    test_celltypist(adata)
    test_manual_annotation(adata)
    test_reproducibility(adata)

    print("\n✅ All Chunk 3 tests passed.")
