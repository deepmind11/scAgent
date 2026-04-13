# Best Practice Reference: Normalization

> Distilled from Heumos et al. 2023, 10x Genomics Analysis Guide (2025), and sc-best-practices.org.

## Why Normalize

Cells have different total counts due to differences in cell size, capture efficiency, and sequencing depth. Normalization makes expression profiles comparable across cells.

## Method Comparison

A benchmark of 22 transformations (Ahlmann-Eltze & Huber 2023, cited by Heumos et al.) showed:

| Method | Best for | Notes |
|--------|----------|-------|
| **Shifted logarithm** (log1p with size factors) | Variance stabilization for dimensionality reduction | Simple, fast, widely used. Do NOT use CPM (1e6) as input — use CP10K (1e4) or scran size factors. |
| **Scran** (pooling-based size factors) | Batch correction tasks, compositional bias | Pools cells with similar depth, estimates size factors via linear regression. Better handles heterogeneous cell sizes. |
| **Analytical Pearson residuals** | HVG selection, rare cell type identification | Fits a GLM with sequencing depth as covariate. Not sensitive to normalization choice. |
| **SCTransform** (Seurat) | Variance stabilization + HVG selection in one step | Models counts as a function of sequencing depth, outputs Pearson residuals. Used by 10x Cloud pipeline. |

## Key Guidance

- "The normalization method should be chosen carefully and based on the subsequent analysis task" (Heumos et al.).
- **For standard Scanpy workflows**: CP10K + log1p is the pragmatic default. It works well for most downstream tasks.
- **For complex integration**: Scran size factors perform well (Luecken et al. 2022 benchmark).
- **For HVG selection on raw counts**: Analytical Pearson residuals or deviance-based selection bypass the normalization choice entirely.
- Do NOT use counts-per-million (target_sum=1e6) — "it reflects an unrealistically large overdispersion" (Heumos et al.).
- Preserve raw counts in `adata.layers['counts']` — needed later for pseudobulk DE (DESeq2/edgeR operate on raw counts).

## 10x-Specific Guidance

- 10x Cloud pipeline uses **SCTransform** from Seurat, which performs normalization, variance stabilization, and HVG selection in one step.
- If following a Scanpy workflow, CP10K + log1p is equivalent in practice for most analyses.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- 10x Genomics 2025: https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data
- Ahlmann-Eltze & Huber, "Comparison of transformations for single-cell RNA-seq data." Nat Methods 20, 665–672 (2023).
