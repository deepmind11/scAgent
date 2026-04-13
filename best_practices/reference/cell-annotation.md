# Best Practice Reference: Cell Type Annotation

> Distilled from Heumos et al. 2023, 10x Genomics Analysis Guide (2025), and sc-best-practices.org.

## Recommended Approach: Three-Step Process (Heumos et al.)

1. **Automated annotation** (first pass) — CellTypist, Clustifyr, or reference mapping
2. **Expert manual annotation** — cross-check against known marker genes
3. **Verification** — validate final annotations, especially for rare or novel populations

## Automated Methods

### Classifier-based (pre-trained models)
- **CellTypist**: trained on cross-tissue immune cell atlas (Conde et al. 2022). 61 models available covering immune cells, brain, lung, gut, etc. Uses thousands of genes for classification.
- **Clustifyr**: correlates cluster expression profiles to reference datasets.
- 10x Cell Ranger (v9.0+): built-in reference-based annotation for human/mouse.
- "Annotation results obtained with pre-trained classifiers are strongly affected by the classifier type and the quality of the training data" (Heumos et al.).

### Reference mapping (label transfer)
- **scArches**: maps query data to a reference atlas using transfer learning. Extends an existing variational autoencoder model.
- **Symphony**: efficient reference mapping via compressed sketching.
- **Azimuth** (Seurat): supervised PCA-based query-to-reference mapping.
- "The quality of the transferred annotations depends on the quality of the reference data, the model and the suitability to the data set" (Heumos et al.).
- Use **uncertainty scores** to flag low-confidence annotations — cells with high uncertainty may represent novel states.

## Manual Annotation

- Based on known marker genes and clustering. Two complementary approaches:
  1. **Markers → clusters**: check where known markers are expressed across clusters
  2. **Clusters → markers**: find DE genes per cluster, look them up in databases
- "Marker-based annotation can be sensitive to the cluster resolution you choose, the robustness and uniqueness of the marker sets, and your knowledge of expected cell types" (sc-best-practices.org).
- Markers from protein-level studies (FACS) may not transfer perfectly to transcriptomic data.

## Cross-Validation (Critical)

- **Always validate automated annotations against marker genes** (10x Genomics, sc-best-practices.org).
- Use dotplots and UMAP expression overlays to check canonical markers match predicted cell types.
- "Automated annotation algorithms should be used with caution and should be regarded as a starting point... expression of known marker genes is still the most accepted support" (sc-best-practices.org).
- Check dendrograms — if a B cell subtype doesn't cluster with other B cells, the label is likely wrong.

## 10x-Specific Guidance

- Cell Ranger v9+ and 10x Cloud perform automated annotation during pipeline runs.
- Validate in Loupe Browser: (1) visualize known marker expression on UMAP, (2) check DE gene lists for canonical markers per annotated type.
- XIST enrichment in female samples is a useful positive control for cross-condition DE validation.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- 10x Genomics 2025: https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data
- sc-best-practices.org, Annotation chapter: https://www.sc-best-practices.org/cellular_structure/annotation.html
- Conde et al. "Cross-tissue immune cell analysis reveals tissue-specific features." Science 376 (2022).
