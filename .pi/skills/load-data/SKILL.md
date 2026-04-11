---
name: load-data
description: Load 10x Genomics Chromium data from .h5 files into AnnData. Use when the user provides a data file, asks to load data, or starts a new analysis.
---

# Load 10x Data

Load a 10x Genomics `filtered_feature_bc_matrix.h5` file and prepare it for analysis.

## When to Use

- User provides a `.h5` file path
- User says "load my data", "start analysis", "open dataset"
- Beginning of any new analysis pipeline

## How to Run

```python
from scagent.tools.loading import load_10x_h5

adata, result = load_10x_h5(
    "path/to/filtered_feature_bc_matrix.h5",
    checkpoint_dir="data/working/checkpoints",
)
```

The function handles:
- Loading the sparse matrix from HDF5
- Fixing duplicate gene names (`var_names_make_unique()`)
- Detecting species (human/mouse) from gene name case
- Computing a SHA-256 file hash for provenance

## What to Show the User

From `result["metrics"]`:
- Number of cells and genes
- Detected species
- Whether duplicate gene names were fixed

Example output:
```
Loaded 11,769 cells × 33,538 genes
Species: human (auto-detected from gene names)
Fixed 24 duplicate gene names.
```

## What to Do Next

Suggest running QC: "The data is loaded. Let's look at quality metrics — I'll compute QC statistics and show you the distributions so we can set filtering thresholds."

## Parameters Reference

See `tools/load_10x_h5.json` for the full parameter schema.

## Notes

- The file must be the **filtered** matrix from Cell Ranger (not raw). The filtered matrix contains only cell-containing barcodes.
- If the user has a `.h5ad` file instead, use `anndata.read_h5ad()` directly — it's already processed.
- For multi-sample experiments, load each sample separately and concatenate later (during the integration step).
