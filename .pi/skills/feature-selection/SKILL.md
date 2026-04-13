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

## What to Show the User

- Number of HVGs selected and % of total genes
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
