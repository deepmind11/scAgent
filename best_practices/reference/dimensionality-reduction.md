# Best Practice Reference: Dimensionality Reduction

> Distilled from Heumos et al. 2023, 10x Genomics Analysis Guide (2025), and sc-best-practices.org.

## Two Distinct Purposes

1. **Data summarization** → PCA (linear, preserves global structure)
2. **Visualization** → UMAP / t-SNE / PHATE (nonlinear, 2D)

These serve different purposes and should not be conflated.

## PCA

- Standard first step after HVG selection. Reduces from ~2,000 HVGs to 30–50 principal components.
- Use the **elbow plot** to identify how many PCs capture meaningful variance. Typically 15–30 PCs contain the biology; beyond that is noise.
- Default of 50 PCs is safe — downstream methods (neighbors, clustering) implicitly select the useful ones.

## UMAP

- Nonlinear projection for **visualization only**.
- "Generally faster than t-SNE and can scale better with larger datasets" (10x Genomics).
- "Often produces embeddings that appear to better represent the global arrangement of clusters" vs. t-SNE (10x Genomics).
- **Critical**: "Relying only on 2D embeddings can lead to misinterpretation of the relationships between cells, and results should not be formulated only on the basis of visual inspection" (Heumos et al.).
- **Never cluster on UMAP coordinates.** Clustering must happen on the KNN graph built from PCA space.
- **Distances between clusters are not meaningful** on UMAP — "must be interpreted with caution" (sc-best-practices.org).

## t-SNE

- "Excellent at preserving local (neighborhood) structure. However, it may distort global relationships, making distances between clusters difficult to interpret" (10x Genomics).
- Less commonly used now that UMAP is standard, but still useful for publication figures.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- 10x Genomics 2025: https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data
- sc-best-practices.org, Clustering chapter.
- Chari et al. "The specious art of single-cell genomics." PLoS Comput Biol (2023) — on UMAP misinterpretation.
