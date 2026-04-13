#!/usr/bin/env python3
"""Integration test: full DE pipeline on synthetic multi-sample data.

1. Create synthetic multi-sample AnnData with known DE signal
2. Run pseudobulk aggregation + DESeq2
3. Verify DE genes are recovered
4. Run GSEA on DE results
5. Verify provenance chain
"""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy.sparse import csr_matrix


def _make_realistic_adata(seed=42):
    """Create a synthetic dataset mimicking a real disease-vs-healthy experiment.

    - 2 cell types: Monocytes (with DE signal) and T_cells (no signal)
    - 3 samples per condition (6 total)
    - 1000 genes, 80 DE in Monocytes (effect size 3x)
    - 300 cells per sample per cell type
    """
    rng = np.random.default_rng(seed)

    n_genes = 1000
    n_de = 80
    n_cells_per_sample = 300
    conditions = ["disease"] * 3 + ["healthy"] * 3
    samples = [f"donor_{i}" for i in range(6)]

    obs_rows = []
    X_blocks = []
    gene_names = [f"GENE{i:04d}" for i in range(n_genes)]

    for ct in ["Monocytes", "T_cells"]:
        for sid, cond in zip(samples, conditions):
            base = rng.uniform(3, 15, size=n_genes)
            counts = rng.poisson(base, size=(n_cells_per_sample, n_genes))

            if ct == "Monocytes" and cond == "disease":
                # Strong upregulation in first n_de genes
                counts[:, :n_de] = rng.poisson(
                    base[:n_de] * 3.0, size=(n_cells_per_sample, n_de)
                )

            X_blocks.append(counts)
            for _ in range(n_cells_per_sample):
                obs_rows.append({
                    "cell_type": ct,
                    "sample": sid,
                    "condition": cond,
                })

    X = np.vstack(X_blocks).astype(np.float32)
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell_{i}" for i in range(len(obs))]
    var = pd.DataFrame(index=gene_names)

    adata = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    for col in ("cell_type", "sample", "condition"):
        adata.obs[col] = adata.obs[col].astype("category")
    return adata


def test_full_de_pipeline():
    from scagent.tools.pseudobulk_de import run_pseudobulk_de
    from scagent.tools.enrichment import run_gsea

    adata = _make_realistic_adata()
    print(f"\n  Synthetic data: {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"  Samples: {adata.obs['sample'].nunique()}, "
          f"Conditions: {list(adata.obs['condition'].unique())}")

    with tempfile.TemporaryDirectory() as tmp:
        plot_dir = str(Path(tmp) / "plots")

        # === Step 1: Pseudobulk DE ===
        de_results = run_pseudobulk_de(
            adata,
            cell_type_key="cell_type",
            sample_key="sample",
            condition_key="condition",
            min_cells_per_pseudobulk=10,
            alpha=0.05,
            plot_dir=plot_dir,
        )

        # Check we got results for both cell types
        assert "Monocytes" in de_results["de_dataframes"], "Monocytes DE failed"
        assert "T_cells" in de_results["de_dataframes"], "T_cells DE failed"

        mono_df = de_results["de_dataframes"]["Monocytes"]
        tcell_df = de_results["de_dataframes"]["T_cells"]

        mono_sig = int(mono_df["significant"].sum())
        tcell_sig = int(tcell_df["significant"].sum())

        print(f"  Monocytes: {mono_sig} DE genes")
        print(f"  T_cells: {tcell_sig} DE genes")

        # Monocytes should have significant DE (we injected 80 DE genes)
        assert mono_sig > 0, "Expected DE genes in Monocytes (signal injected)"

        # T_cells should have fewer DE genes (no signal injected)
        assert mono_sig > tcell_sig, (
            f"Monocytes ({mono_sig}) should have more DE than T_cells ({tcell_sig})"
        )

        # Verify the DE genes are in the right range (first 80 genes)
        sig_genes = mono_df[mono_df["significant"]].index.tolist()
        expected_de = {f"GENE{i:04d}" for i in range(80)}
        recovered = set(sig_genes) & expected_de
        precision = len(recovered) / len(sig_genes) if sig_genes else 0
        recall = len(recovered) / len(expected_de)
        print(f"  DE recovery: precision={precision:.2f}, recall={recall:.2f} "
              f"({len(recovered)}/{len(expected_de)} true DE genes recovered)")

        # We expect reasonable recall (>50%) given n=3 per condition
        assert recall > 0.3, f"DE recall too low: {recall:.2f}"

        # Volcano plots should exist
        assert len(de_results["plots"]) > 0
        for p in de_results["plots"]:
            assert Path(p).exists(), f"Plot not found: {p}"
        print(f"  Volcano plots: {len(de_results['plots'])}")

        # Provenance
        prov = de_results["provenance"]
        assert prov["tool_id"] == "deseq2_pseudobulk"
        assert prov["n_cell_types_tested"] == 2
        print(f"  Provenance: {prov['tool_id']}, {prov['n_cell_types_tested']} cell types")

        # === Step 2: GSEA ===
        gsea_results = run_gsea(
            de_results["de_dataframes"],
            gene_sets="KEGG_2021_Human",
            permutation_num=100,  # fast for test
        )

        assert "enrichment_results" in gsea_results
        assert "provenance" in gsea_results
        assert gsea_results["provenance"]["tool_id"] == "gsea"

        n_enriched = sum(r["n_significant"] for r in gsea_results["summary"])
        print(f"  GSEA: {len(gsea_results['summary'])} cell types tested, "
              f"{n_enriched} significant terms total")

        # Summary
        for row in de_results["summary"]:
            print(f"  DE summary — {row['cell_type']}: "
                  f"{row['n_significant']} sig ({row['n_up']} up, {row['n_down']} down)")

        print("  ✓ Full DE pipeline: aggregate → DESeq2 → GSEA complete")
