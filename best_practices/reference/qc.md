# Best Practice Reference: Quality Control

> Distilled from Heumos et al. 2023 (Nat Rev Genetics), 10x Genomics Analysis Guide (2025),
> and sc-best-practices.org. Not dogma — use as informed defaults.

## Core QC Metrics

Three metrics are evaluated **jointly** per sample — never use one in isolation:

1. **Count depth** (total UMIs per cell) — low values indicate empty droplets or debris
2. **Genes detected** (n_genes_by_counts) — low values indicate poor-quality cells; very high values may indicate doublets
3. **Mitochondrial fraction** (pct_counts_mt) — high values indicate dying/lysed cells where cytoplasmic RNA leaked out

## Thresholds

- **Use MAD-based filtering** as a starting point (Heumos et al.: "sample-wise automatic filtering based on the number of median absolute deviations"). Then adjust based on tissue context.
- **QC must be per-sample**, not across the merged dataset — thresholds can vary substantially between samples.
- **If the mito distribution is bimodal**, the high-mito peak may represent stressed or apoptotic cells with biological meaning — discuss with the user before removing.
- **Set permissive thresholds initially** and potentially remove more cells during re-analysis (Heumos et al.).

## Ambient RNA

- Cell-free RNA contaminates genuine cells, causing marker genes to "bleed" across populations.
- **SoupX** estimates contamination from empty droplet profiles and cell clusters.
- **CellBender** uses an unsupervised Bayesian model requiring no prior knowledge of cell types.
- Consider ambient RNA removal as an initial step, especially for solid tissues (Heumos et al.).
- 10x Genomics: ambient RNA correction is important when looking for subtle expression patterns or rare cell types.

## Doublets (see also doublet-detection reference)

- Doublets violate the one-cell-per-droplet assumption. Heterotypic doublets (different cell types) are particularly problematic.
- **scDblFinder** outperforms other methods in accuracy and speed across multiple benchmarks (Heumos et al.).
- Apply multiple doublet detection methods and compare results to increase accuracy.
- QC strategy "often needs to be reassessed during downstream analysis when low-quality cells and doublets cluster together" (Heumos et al.).

## 10x-Specific Guidance

- Always check the Cell Ranger `web_summary.html` first — look for alerts, the barcode rank plot ("cliff-and-knee" shape), and cell recovery vs. target.
- Metrics to verify: "Confidently mapped reads in cells" should be high (>90%); median genes/cell should match tissue expectations.
- Use Loupe Browser for interactive QC exploration and filtering.

## Sources

- Heumos et al. "Best practices for single-cell analysis across modalities." Nat Rev Genetics 24, 550–572 (2023). https://doi.org/10.1038/s41576-023-00586-w
- 10x Genomics. "Best Practices for Analysis of 10x Genomics Single Cell RNA-seq Data." (2025). https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data
- sc-best-practices.org — Preprocessing chapters. https://www.sc-best-practices.org
