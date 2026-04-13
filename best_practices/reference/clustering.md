# Best Practice Reference: Clustering

> Distilled from Heumos et al. 2023, 10x Genomics Analysis Guide (2025), and sc-best-practices.org.

## Algorithm

- **Use Leiden** — it yields guaranteed connected communities and is computationally more efficient than Louvain (Heumos et al.).
- "Since the Louvain algorithm is no longer maintained, using Leiden instead is preferred" (sc-best-practices.org).
- Both operate on the **KNN graph** computed on the PCA-reduced (or integration-corrected) representation.
- **Never cluster on UMAP coordinates.**

## Resolution Parameter

- Controls the granularity of clustering. Higher resolution → more clusters.
- "We recommend using the Leiden algorithm at different resolutions to obtain an ideal clustering for annotating cells" (Heumos et al.).
- Typical ranges: 0.3–0.5 for major lineages, 0.8–1.2 for standard analysis, 1.5+ for fine subtypes.
- **Sub-clustering** on specific clusters allows focusing on substructure without affecting other clusters.

## Graph-Based vs. K-Means

- Graph-based (Leiden/Louvain) is preferred for scRNA-seq: "can better identify clusters that are complex and not necessarily round" and "determines cluster numbers from the data's structure, not from a user-predefined number" (10x Genomics).
- K-means is useful for quick initial assessment with small k (e.g., 3–5) to identify major populations, or when k is confidently known.

## Practical Tips

- K for the KNN graph: typically 5–100 depending on dataset size. Default of 15 works for most cases.
- Run `n_iterations=2` in Leiden for good quality at reasonable speed; use `-1` for full convergence on smaller datasets (sc-best-practices.org).
- Save clustering results under resolution-specific keys (e.g., `leiden_res0_5`, `leiden_res1`) for comparison.
- Check seed sensitivity — if clustering changes substantially with different random seeds (ARI < 0.9), boundaries are not well-defined.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- 10x Genomics 2025: https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data
- sc-best-practices.org, Clustering chapter: https://www.sc-best-practices.org/cellular_structure/clustering.html
- Traag et al. "From Louvain to Leiden: guaranteeing well-connected communities." Sci Rep 9, 5233 (2019).
