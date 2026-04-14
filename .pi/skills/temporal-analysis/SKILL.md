---
name: temporal-analysis
description: Analyze time-series / longitudinal scRNA-seq experiments. Use when the paradigm is temporal_longitudinal, or when data includes multiple timepoints from the same biological process.
---

# Skill: Temporal / Longitudinal Analysis

Analyze scRNA-seq data with multiple timepoints (e.g., disease progression, drug response over time, developmental time courses).

## When to Use

- **Paradigm is `temporal_longitudinal`** — the DAG requires this workflow.
- **User has samples from multiple timepoints** (Day 0, Day 7, Day 14, etc.).
- **Interest in how cell composition or gene expression changes over time.**

## Key Principle: Timepoints = Batches

Different timepoints are almost always processed in different batches. **Batch correction is mandatory** before any cross-timepoint analysis. The DAG enforces this — batch correction is required, not optional.

## Workflow

### 1. Standard Preprocessing (with mandatory batch correction)

The DAG runs the standard prefix (QC → normalize → HVG → PCA) then **always** applies batch correction (Harmony or scVI). This is not optional because timepoint-batch confounding is expected. [BP-1 pp. 552-553]

### 2. Pseudobulk DE (Timepoint Contrasts)

Pairwise comparison between timepoints using the same pseudobulk DE pipeline as `disease_vs_healthy`.

```python
from scagent.tools.pseudobulk_de import run_pseudobulk_de

# Treat timepoint as the condition
results = run_pseudobulk_de(
    adata,
    cell_type_key="cell_type",
    sample_key="donor",
    condition_key="timepoint",  # e.g., "Day0" vs "Day7"
)
```

**Critical:** Pseudobulk aggregation by sample is mandatory. Cell-level tests on temporal data inflate FDR even more than in standard cross-condition DE because of correlated time effects within samples. [BP-1 p. 555]

### 3. Composition Analysis (Proportion Trends)

How do cell type proportions change over time?

- **scCODA** for discrete timepoint comparisons (labeled clusters)
- **Milo** for continuous time covariates — Milo supports GLM testing with continuous covariates, directly applicable to time-course DA. [BP-2 Ch. 18]

```python
# Milo with continuous timepoint (via pertpy)
# design="~ timepoint" where timepoint is numeric
```

### 4. Pathway Enrichment (Optional)

Run GSEA on DE results to identify pathways that change over time.

## Guard Rails

1. **Batch correction is always required.** The DAG enforces this. Do not skip it.
2. **Pseudobulk aggregation is mandatory.** Never do cell-level DE for cross-timepoint comparisons.
3. **Minimum replicates per timepoint.** At least 2 biological replicates per timepoint for valid DE. Warn if <3.
4. **Cell cycle effects should be assessed.** Time-course data often has varying cell cycle composition. Consider regressing out cell cycle effects before DE. [BP-1]
5. **If only 2 timepoints, treat as `disease_vs_healthy`.** Temporal analysis adds value primarily with ≥3 timepoints where trends can be assessed.

## Interpretation Guidance

- **Composition changes ≠ gene expression changes.** A cell type can change proportion without changing its expression profile, and vice versa. Both should be checked.
- **Pseudobulk DE between consecutive timepoints** may miss gradual trends. Consider regression on time if ≥3 timepoints.
- **Batch effects can mimic temporal effects.** If each timepoint was processed in a separate batch, some "temporal" signals may be batch artifacts. Check with PCA on pseudobulk samples. [BP-1]

## Best-Practice References

- [BP-1] §"Removing confounding sources of variation" (pp. 552-553)
- [BP-1] §"Differential gene expression analysis" (p. 555)
- [BP-2] Ch. 18 §"Without labeled clusters" — Milo time-course example
- [BP-3] QC/normalization per-sample before integration
