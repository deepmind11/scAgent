# Best Practice Reference: Feature Selection (HVG)

> Distilled from Heumos et al. 2023 and sc-best-practices.org.

## Purpose

Select genes that explain biological variation (between subpopulations) rather than technical noise (within subpopulations). The count matrix typically has 20,000–30,000 genes; most are uninformative.

## Method Recommendation

- **Deviance** (from the `scry` package) identifies highly informative genes by fitting a gene-wise model that assumes constant expression across all cells, then quantifying which genes violate this assumption. "It performed favourably for identifying genes with high variance across subpopulations" in an independent comparison (Heumos et al.).
- **Key advantage**: Deviance works on raw counts and is therefore **not sensitive to the normalization choice**.
- **Seurat v3 method** (`flavor='seurat_v3'`) selects HVGs by fitting a mean-variance relationship — the standard default in most Scanpy tutorials.
- **Analytical Pearson residuals** simultaneously normalize and select features — an alternative when you want to skip explicit normalization.

## Practical Guidance

- **2,000 genes** is the standard default for most analyses.
- **3,000–5,000** for heterogeneous tissues with many cell types (whole-organ atlases).
- **1,000** for conservative selection or datasets with few expected cell types.
- For multi-batch datasets, use `batch_key` to select HVGs per batch and take the union — this avoids batch-specific genes dominating the selection.
- After feature selection, dimensionality is further reduced by PCA (typically 30–50 components).

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- Townes et al. "Feature selection and dimension reduction for single-cell RNA-Seq based on a multinomial model." Genome Biology 20, 295 (2019).
