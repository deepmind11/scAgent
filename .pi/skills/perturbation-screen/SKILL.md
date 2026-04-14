---
name: perturbation-screen
description: Analyze Perturb-seq / CROP-seq CRISPR screen data. Use when the paradigm is perturbation_screen or when data includes CRISPR guide assignments.
---

# Skill: Perturbation Screen Analysis

Analyze single-cell CRISPR perturbation screens (Perturb-seq, CROP-seq) on 10x Chromium.

## When to Use

- **Paradigm is `perturbation_screen`** — the DAG includes guide assignment and perturbation DE.
- **User has CRISPR guide/feature barcode data** alongside gene expression.
- **User asks about perturbation effects, gene knockdowns, or screen results.**

## Workflow

### Step 1: Guide Assignment

Assign CRISPR guide identities to cells from Cell Ranger multi output.

```python
from scagent.tools.perturbation import assign_guides

result = assign_guides(adata, guide_calls_key="guide_ids")
# Adds: adata.obs["guide"], adata.obs["perturbation"], adata.obs["n_guides"]
```

### Step 2: Perturbation DE

Compare each perturbation to non-targeting controls.

```python
from scagent.tools.perturbation import run_perturbation_de

results = run_perturbation_de(
    adata,
    control_label="non-targeting",
    min_cells=50,
    alpha=0.05,
)
```

### Step 3: Pathway Enrichment (reuse existing tool)

Run GSEA on per-perturbation DE results to identify affected pathways.

## Key Difference from Cross-Condition DE

In Perturb-seq, **cells ARE the replicates** — each cell receives an independent guide. Cell-level Wilcoxon tests are acceptable here, unlike cross-condition DE where pseudobulk is mandatory. [BP-1, BP-2 Ch. 20]

## Guard Rails

1. **Check guide assignment rate.** Warn if <50% of cells have guides.
2. **Minimum cells per perturbation: 50.** Skip perturbations below this threshold.
3. **Non-targeting controls must be present** and have sufficient cells.
4. **Perturbation ≠ complete knockout.** Knockdown efficiency varies — effect sizes may be attenuated.
5. **Multi-guide cells** should be flagged and excluded from per-perturbation analysis.

## Best-Practice References

- [BP-1] §"Inferring perturbation effects" (p. 557)
- [BP-2] Ch. 20 — Perturbation modeling
