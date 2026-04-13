# Best Practice Reference: Batch Integration

> Distilled from Heumos et al. 2023, 10x Genomics Analysis Guide (2025), and sc-best-practices.org.

## When to Integrate

- Only for **multi-sample or multi-batch** experiments.
- First check: after merging samples and computing UMAP, do cells cluster by sample/batch instead of cell type? If yes, integration is needed.
- **Caution**: "If your biological condition of interest perfectly correlates with batch... distinguishing true biological effects from technical batch effects becomes inherently more challenging" (10x Genomics).

## Method Recommendations (from Luecken et al. 2022 benchmark, 16 methods × 14 metrics)

| Method | Best for | Notes |
|--------|----------|-------|
| **Harmony** | Simple-to-moderate batch structures | Fast, operates on PCA space, does not alter counts. "Identified as a high-performance batch correction technique" (10x Genomics). Recommended default. |
| **scVI / scANVI** | Complex integration, atlas-scale | Deep learning. scANVI can use cell-type labels to better conserve biology. Best for complex tasks. |
| **Scanorama** | Atlas integration | Mutual nearest neighbors approach, performed well for complex tasks. |
| **CCA (Seurat)** | Simple batch structures | Linear embedding, fast. |
| **BBKNN** | Replaces neighbor graph directly | Does not produce a corrected embedding — modifies the graph. |

## Key Principles

- "Linear-embedding models such as CCA and Harmony were shown to perform well for simpler integration tasks with distinct batch structures. For complex tasks such as atlas integration, deep-learning approaches such as scANVI, scVI and Scanorama performed best" (Heumos et al.).
- Use **scIB** metrics to evaluate integration quality (batch mixing vs. biological conservation).
- Integration typically operates on embeddings (PCA), not on raw counts. **Harmony does not alter the underlying UMI counts** (10x Genomics) — important for downstream DE.

## Cell Cycle Effects

- Cell cycle can be a confounding biological variable rather than a technical batch.
- Removing cell cycle effects "can be favourable for downstream analysis; however, knowing whether cells are cycling may provide valuable insights" (Heumos et al.).
- Recommendation: use Scanpy/Seurat built-in cell cycle scoring first; if needed, apply Tricycle for more complex datasets.

## 10x-Specific Guidance

- 10x Cloud provides an SCTransform + Harmony workflow.
- Toggle "skip batch correction" for first-time merging to assess if batch effects exist before correcting.
- Harmony corrects the PCA embedding; downstream UMAP and clustering use the corrected embedding.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- 10x Genomics 2025: https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data
- Luecken et al. "Benchmarking atlas-level data integration in single-cell genomics." Nat Methods 19, 41–50 (2022).
