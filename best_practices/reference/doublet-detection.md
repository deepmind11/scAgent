# Best Practice Reference: Doublet Detection

> Distilled from Heumos et al. 2023, 10x Genomics Analysis Guide (2025), and sc-best-practices.org.

## Method Recommendation

- **scDblFinder** is the top-performing method across multiple independent benchmarks for both accuracy and computational efficiency (Heumos et al.).
- Scrublet (used in scAgent) is a reasonable alternative that generates artificial doublets and compares them to measured cells — the same core approach as scDblFinder.
- "It can be beneficial to apply multiple doublet detection methods and compare the results to increase the accuracy of doublet detection" (Heumos et al.).

## When to Run

- After initial QC filtering, before normalization.
- Must run on **raw counts** (not normalized).
- For multi-sample experiments, run per sample before merging — doublet rates and profiles differ across samples.

## Expected Rates

- 10x Chromium: ~0.8% per 1,000 cells captured (so 10K cells → ~8%).
- Rates much higher than expected may indicate over-calling; rates much lower may indicate a failed threshold.

## Interpretation Pitfalls

- Doublets from the same cell type (homotypic) are very hard to detect and often remain in the data.
- Heterotypic doublets (different cell types) are the main concern — they create artificial intermediate populations.
- The QC strategy "often needs to be reassessed during downstream analysis when low-quality cells and doublets cluster together" (Heumos et al.).

## 10x-Specific Guidance

- 10x Genomics considers algorithmic doublet detection an "optional advanced QC refinement" step — manual review of diagnostic plots is often sufficient as an initial method.
- Commonly recommended tools: scDblFinder and DoubletFinder.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- 10x Genomics 2025: https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data
