# Best Practice Reference: Pseudobulk Differential Expression

> Distilled from Heumos et al. 2023, 10x Genomics Analysis Guide (2025), and sc-best-practices.org.

## The Pseudoreplication Problem (THE Critical Issue)

- Individual cells from the same sample are NOT independent replicates.
- "Pseudoreplication leads to an inflated false discovery rate (FDR) as DGE methods do not account for the inherent correlation of replicates (cells from the same individual)" (Heumos et al.).
- This is confirmed by Squair et al. 2021, Zimmerman et al. 2021, and Murphy & Skene 2022.
- **Never use cell-level Wilcoxon/t-test for cross-condition comparisons.** Even if the user asks.

## Recommended Methods

| Method | Type | Notes |
|--------|------|-------|
| **edgeR** (quasi-likelihood) | Pseudobulk | Accounts for dispersion uncertainty. Flexible design matrices. |
| **DESeq2** | Pseudobulk | Robust, widely used. Recommended for complex experimental designs. |
| **Limma-voom** | Pseudobulk | Fast, flexible. Good for large sample sizes. |
| **MAST** (with random effects) | Mixed model | Models cells individually but accounts for donor as random effect. Slower but preserves cell-level info. |

- "Pseudobulk methods with sum aggregation and mixed models such as MAST with random effect setting were found to be superior to naive methods" (Heumos et al.).
- Murphy & Skene 2022 confirmed pseudobulk methods perform best; whether sum or mean aggregation is better needs further study.

## Practical Requirements

- **Minimum 2 biological replicates per condition** — refuse to run with fewer.
- **3+ replicates recommended** for adequate statistical power.
- "The best way to increase statistical power is to increase the number of independent experimental samples" — NOT more cells per sample (Zimmerman et al. 2021).
- **Use raw counts** as input to edgeR/DESeq2, NOT normalized values.
- Aggregate per cell type × sample before testing.
- Filter: require minimum cells per pseudobulk (typically ≥10–30 cells).

## Design Matrices

- "The validity of DGE results strongly depends on the capture of the major axis of variation in the statistical model" (Heumos et al.).
- Run PCA on pseudobulk samples first to identify confounders (sc-best-practices.org).
- Include relevant covariates (donor, sex, age) in the design matrix.
- Use flexible methods (edgeR, DESeq2, limma) that support complex designs.
- Always correct p-values for multiple testing (Benjamini-Hochberg).

## 10x-Specific Guidance

- Loupe Browser supports pseudobulk multi-sample comparison (treats each sample as a replicate, not each cell).
- Harmony corrects embeddings but does NOT alter counts — for DE on counts, include batch as a covariate in the statistical model if it's a technical variable.
- "If the 'batch' variable used by Harmony is a biological variable (e.g., patient ID)... there is no need to adjust during DE as this is biological signal rather than artifact" (10x Genomics).

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- 10x Genomics 2025: https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data
- sc-best-practices.org, DE chapter: https://www.sc-best-practices.org/conditions/differential_gene_expression.html
- Squair et al. "Confronting false discoveries in single-cell DE." Nat Commun 12, 5692 (2021).
- Zimmerman et al. "A practical solution to pseudoreplication bias." Nat Commun 12, 738 (2021).
