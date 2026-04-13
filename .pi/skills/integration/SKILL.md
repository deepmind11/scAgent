---
name: integration
description: Batch correction and data integration for multi-sample experiments. Use after PCA when data contains multiple batches, samples, or donors. Skip for single-sample experiments.
---

# Batch Integration

Correct batch effects so that cells from different samples mix by cell type, not by batch.

## When to Use

**Only for multi-sample or multi-batch experiments.** Check if `adata.obs` contains a column with batch/sample/donor labels.

- If the data has a batch column → run integration
- If single-sample (like PBMC 10k) → skip this step entirely

## How to Decide

Ask the user: "Does your dataset have multiple samples or batches?" If yes, ask which obs column contains the batch labels.

## How to Run (Harmony)

```python
from scagent.tools.integration import run_harmony

result = run_harmony(adata, key="batch", random_state=0)
```

After Harmony, the corrected embedding is in `adata.obsm['X_pca_harmony']`. Use it for the neighbor graph:

```python
from scagent.tools.neighbors import compute_neighbors

compute_neighbors(adata, use_rep="X_pca_harmony")
```

## When NOT to Use

- Single-sample experiments (no batch to correct)
- If samples were processed identically and UMAP shows good mixing already
- If you want to preserve batch differences (e.g., comparing conditions)

## Evaluating Integration

After integration, generate a UMAP colored by batch. Batches should overlap (mix by cell type). If batches are still separated, integration may have been insufficient.

## Alternative Methods

Harmony is the default (fast, good for most cases). Other options:
- **scVI** — deep learning, better for complex batch structures, slower
- **Scanorama** — mutual nearest neighbors, fast
- **BBKNN** — replaces the neighbor graph step entirely

These are registered in the tool registry but not yet implemented as wrapper functions.

## Best Practice Reference

Load `best_practices/reference/integration.md` for benchmark-backed method selection (Harmony vs. scVI vs. Scanorama) based on integration complexity.

## Parameters Reference

See `tools/harmony.json`.
