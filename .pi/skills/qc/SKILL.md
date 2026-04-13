---
name: qc
description: Run quality control on single-cell data. Compute QC metrics, visualise distributions, recommend filtering thresholds, and apply cell and gene filtering. Use when the user asks for QC, quality check, filtering, or after loading data.
---

# Quality Control & Filtering

Compute quality metrics, show distributions, recommend thresholds, and filter low-quality cells and rarely-expressed genes.

## When to Use

- After loading data (always the first analysis step)
- User asks for QC, quality check, or cell filtering
- User wants to adjust filtering thresholds and re-filter

## Step 1: Compute QC Metrics

```python
from scagent.tools.qc import calculate_qc_metrics

result = calculate_qc_metrics(adata, plot_dir="data/working/plots")
```

This adds `pct_counts_mt`, `n_genes_by_counts`, and `total_counts` to `adata.obs`. It also generates three diagnostic plots:

1. **Violin plots** — distribution of genes/cell, UMIs/cell, % mito
2. **Scatter plot** — total counts vs genes detected, colored by % mito
3. **Mito histogram** — distribution of mitochondrial content

**Show all three plots to the user.** Always show distributions before recommending thresholds.

## Step 2: Discuss Thresholds

Present the MAD-based recommendations from `result["metrics"]["recommended_thresholds"]` alongside tissue-specific guidance:

| Tissue | Typical max % mito | Notes |
|--------|-------------------|-------|
| PBMCs (clean) | 5–10% | Low mito expected |
| PBMCs (variable) | 10–20% | Some datasets have higher baseline |
| Solid tissue (lung, gut) | 10–15% | More ambient RNA |
| Heart, kidney | 15–25% | Naturally high mito content |
| Tumour samples | 10–20% | Mixed viability |

**Key rules:**
- If the mito distribution is bimodal, the high-mito population may be real biology (stressed or apoptotic cells), not just debris. Ask the user.
- Always show the actual distribution. Never apply a threshold without discussing it.
- The MAD-based thresholds are a starting point, not gospel.

## Step 3: Filter Cells

After the user confirms thresholds:

```python
from scagent.tools.qc import filter_cells

adata, result = filter_cells(
    adata,
    min_genes=200,
    max_genes=5000,
    max_pct_mito=20.0,    # adjust based on discussion
    checkpoint_dir="data/working/checkpoints",
)
```

**Show the result:** "Filtered 11,769 → 10,991 cells (6.6% removed)."

**If >50% of cells are removed**, the function returns a warning. Stop and discuss with the user — the thresholds are likely too aggressive.

## Step 4: Filter Genes

```python
from scagent.tools.qc import filter_genes

adata, result = filter_genes(
    adata,
    min_cells=3,
    checkpoint_dir="data/working/checkpoints",
)
```

`min_cells=3` is the standard default (gene must be expressed in at least 3 cells).

## What to Do Next

Suggest doublet detection: "QC filtering is done. Next I'll check for cell doublets using Scrublet — this should be done before normalization."

Or if the user wants to skip doublets, suggest normalization.

## Going Back

If the user wants to try different thresholds, load the checkpoint:

```python
import anndata
adata = anndata.read_h5ad("data/working/checkpoints/after_load.h5ad")
# Re-run calculate_qc_metrics and filter_cells with new thresholds
```

Or use the pipeline runner:

```python
from scagent.tools.pipeline import run_preprocessing
adata, results = run_preprocessing(
    filename,
    max_pct_mito=15.0,  # new threshold
    checkpoint_dir="data/working/checkpoints",
    start_from="after_load",
    stop_after="filter_genes",
)
```

## Best Practice Reference

Before recommending thresholds or making QC decisions, load and consult `best_practices/reference/qc.md` for literature-backed guidance on filtering strategies, ambient RNA correction, and per-sample QC.

## Parameters Reference

See `tools/filter_cells.json` and `tools/filter_genes.json` for full parameter schemas.
