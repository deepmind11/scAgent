---
name: feature-selection
description: Select highly variable genes (HVGs) for downstream dimensionality reduction and clustering. Use after normalization and before PCA.
---

# Highly Variable Gene Selection

Select the most informative genes — those with the highest biological variability across cells.

## When to Use

- After normalization (log-normalized data in `adata.X`)
- Before PCA

## How to Run

```python
from scagent.tools.feature_selection import select_hvg

result = select_hvg(adata, n_top_genes=2000, plot_dir="data/working/plots")
```

The default is `flavor='seurat_v3'`, which uses raw counts. The function automatically finds raw counts in `adata.layers['counts']`, `adata.layers['raw_counts']`, or `adata.X` itself. If no raw counts are found, it falls back to `flavor='seurat'` on log-normalized data.

## Critical: Use Raw Counts for HVG Selection

**Always use `seurat_v3` on raw counts when available.** This is the standard best practice (Hafemeister & Satija 2019).

- `seurat_v3` uses a variance-stabilizing transformation on count data — it recovers biologically meaningful markers (cell-type markers, signaling genes) that the older `seurat` method misses.
- The older `seurat` flavor on log-normalized data is biased toward lowly-expressed genes and misses important markers.
- Before running HVG, **inspect the dataset** for raw counts: check `adata.layers` (look for `'counts'` or `'raw_counts'`), check whether `adata.X` contains integers, and check `adata.raw`.

## What to Show the User

- Number of HVGs selected and % of total genes
- Which flavor and raw counts source was used
- The mean-vs-dispersion plot (`hvg_selection.png`)

## Parameter Guidance

- **2000** is the standard default
- **3000–5000** for heterogeneous tissues with many cell types (e.g., whole-organ atlases)
- **1000** if you want a very conservative set

## Important Notes

- HVGs are *marked* in `adata.var['highly_variable']` — the data is NOT subset. PCA uses `use_highly_variable=True` to restrict automatically.
- `adata.raw` (set during normalization) retains all genes for marker detection and plotting.
- If the dataset has multiple batches, use `batch_key` to select HVGs per batch and take the union.

## What to Do Next

Suggest PCA: "HVGs selected. Next I'll run PCA to reduce the dimensionality — this captures the main axes of variation in 50 components."

## Best Practice Reference

Load `best_practices/reference/feature-selection.md` for benchmark-backed guidance on deviance vs. Seurat v3 HVG methods and gene count recommendations.

## Parameters Reference

See `tools/highly_variable_genes.json`.
