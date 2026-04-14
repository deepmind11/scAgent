---
name: composition-analysis
description: Detect changes in cell-type proportions between conditions using scCODA or Milo. Use when the paradigm involves multi-condition data and the user asks about proportion changes, differential abundance, or cell-type composition shifts.
---

# Skill: Compositional Analysis

Detect which cell types change in relative abundance between experimental conditions.

## When to Use

- **Paradigm is `disease_vs_healthy` or `temporal_longitudinal`** — composition is a DAG step.
- **User asks "which cell types change in proportion?"**
- **Multi-condition experiment with biological replicates.**

## CRITICAL: Why Not Simple Proportion Tests?

Cell-type counts from scRNA-seq are **compositional data** — proportions sum to 1 per sample. If one cell type increases, others MUST decrease proportionally, even if unchanged biologically. Naive per-type tests (Wilcoxon, Fisher, Poisson regression) produce **systematic false positives** due to this induced negative correlation. [BP-1, BP-2 Ch. 18]

Example: if disease doubles cell type A but leaves B and C unchanged, sequencing a fixed sample of 600 cells makes B and C APPEAR to decrease. scCODA correctly identifies only A as changed; naive tests flag B and C as false positives.

## Two Approaches

### 1. scCODA — Labeled clusters (DEFAULT)

Use when cell types are well-defined and you want to know which types change.

```python
from scagent.tools.composition import run_sccoda

result = run_sccoda(
    adata,
    condition_key="condition",
    sample_key="donor",
    cell_type_key="cell_type",
    reference_cell_type="automatic",  # or specify a stable type
    fdr=0.05,
    plot_dir="plots/composition",
)
```

**Key decisions:**
- **Reference cell type:** Must be one believed unchanged. Use `"automatic"` if unsure. [BP-2 Ch. 18]
- **FDR threshold:** Default 0.05. Loosen to 0.2 for exploratory analysis. [BP-2 Ch. 18]
- **NUTS acceptance rate:** Should be 0.4–0.9. Outside range → suspect sampling issues.

### 2. Milo — KNN-based (no predefined labels)

Use when clusters are unclear, transitional states exist, or you want sub-cell-type resolution.

```python
from scagent.tools.composition import run_milo

result = run_milo(
    adata,
    condition_key="condition",
    sample_key="donor",
    cell_type_key="cell_type",  # optional, for annotation
    plot_dir="plots/composition",
)
```

**Key decisions:**
- **n_neighbors:** Must be ≥ 3 × n_samples for adequate power. [BP-2 Ch. 18]
- **Batch correction first:** If conditions and batches are confounded, integrate with scVI before Milo.
- **Continuous covariates:** Milo supports `design="~ timepoint"` for time-course DA testing.

## Guard Rails

1. **Never use naive proportion tests for cross-condition comparisons.** Always use scCODA or Milo.
2. **scCODA requires ≥2 samples per condition.** Warn if <3.
3. **Milo requires sufficient neighbours.** Check n_neighbors ≥ 3 × n_samples.
4. **Composition ≠ DE.** Composition finds proportion changes; DE finds gene expression changes. Both should be checked — they answer different questions.
5. **If conditions and batch are confounded, warn the user** and recommend scVI integration before Milo. [BP-2 Ch. 18]

## Best-Practice References

- [BP-1] §"Deciphering changes in cell composition" (p. 555)
- [BP-2] Ch. 18 — Compositional analysis (scCODA + Milo tutorial)
- Büttner et al. 2021, Nat Commun — scCODA: MCC 0.64 vs ~0.20 naive
