---
name: pseudobulk-de
description: Run cross-condition differential expression using pseudobulk aggregation + DESeq2. Use when comparing conditions (disease vs. healthy, treated vs. control) in scRNA-seq data, or when the user asks for DE between conditions.
---

# Skill: Pseudobulk Differential Expression

Run cross-condition DE using pseudobulk aggregation + DESeq2. This is the **only** correct method for comparing conditions (disease vs. healthy, treated vs. control) in scRNA-seq data.

## When to Use

- **Paradigm is `disease_vs_healthy`, `perturbation_screen`, or `temporal_longitudinal`** — the DAG requires this step.
- **User asks for DE between conditions** — always use pseudobulk, never cell-level Wilcoxon.
- **After cell type annotation** — DE is run per cell type.

## CRITICAL: Why Not Cell-Level Wilcoxon?

Individual cells from the same sample are NOT independent replicates. Treating them as such massively inflates false positives (p-values in the 10⁻⁵⁰ range for noise). The biological replicates are **samples/donors**, not cells.

- **Wilcoxon/t-test**: for finding **cluster markers** (one cluster vs. rest)
- **Pseudobulk DESeq2**: for comparing **conditions** (disease vs. healthy)

This distinction is the #1 statistical error in the field (Squair et al. 2021).

## How to Use

```python
from scagent.tools.pseudobulk_de import run_pseudobulk_de

results = run_pseudobulk_de(
    adata,
    cell_type_key="cell_type",   # from annotation step
    sample_key="donor",          # biological replicates
    condition_key="condition",   # disease vs healthy
    min_cells_per_pseudobulk=10,
    alpha=0.05,
    plot_dir="plots/de",
)

# Results include per-cell-type DE tables + volcano plots
for row in results["summary"]:
    print(f"{row['cell_type']}: {row['n_significant']} DE genes "
          f"({row['n_up']} up, {row['n_down']} down)")
```

## Required Metadata

The agent needs these from the experiment context:
- **cell_type_key**: which `adata.obs` column has cell type labels (from annotation)
- **sample_key**: which column identifies biological samples (from `context.design.batch_key` or `context.samples`)
- **condition_key**: which column has condition labels (from `context.design.conditions`)

## Guard Rails

- **Refuse if <2 replicates per condition.** Explain why.
- **Warn if <3 replicates.** DESeq2 can run with n=2 but power is low.
- **Warn if cell types are dropped** for too few cells — this is expected for rare types.
- **Never fall back to cell-level DE for cross-condition comparisons.** Even if the user asks.

## What to Tell the User

**Before running:**
> "Running pseudobulk DE: aggregating cells by cell type × sample, then using DESeq2. This correctly accounts for biological replicate variance. Found 8 cell types × 4 samples per condition."

**After running:**
> "DE results for 8 cell types. CD14 Monocytes: 342 DE genes (187 up, 155 down). CD4 T cells: 45 DE genes. See volcano plots in plots/de/."

**If they ask for Wilcoxon between conditions:**
> "Wilcoxon rank-sum treats each cell as an independent replicate, which inflates false positives dramatically. For cross-condition comparisons, pseudobulk DESeq2 is the correct method. Shall I run that instead?"

## Best Practice Reference

Load `best_practices/reference/pseudobulk-de.md` for detailed literature-backed guidance on the pseudoreplication problem, method comparison (edgeR vs. DESeq2 vs. MAST), and design matrix construction.
