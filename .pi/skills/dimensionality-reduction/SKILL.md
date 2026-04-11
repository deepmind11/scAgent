---
name: dimensionality-reduction
description: Run PCA for dimensionality reduction and UMAP for visualization. PCA runs after HVG selection; UMAP runs after the neighbor graph is built.
---

# Dimensionality Reduction: PCA and UMAP

Two distinct steps that appear at different points in the pipeline.

## PCA — After HVG Selection

`run_pca()` handles two steps: z-score scaling (each gene to mean=0, var=1) then PCA. Scaling is standard and required so that highly-expressed genes don't dominate the principal components.

```python
from scagent.tools.pca import run_pca

result = run_pca(adata, n_comps=50, random_state=0, plot_dir="data/working/plots")
```

After this call, `adata.X` contains z-scores (not log-counts). `adata.raw` still holds log-normalized values for all genes — used for plotting and marker detection.

**Show the user:**
- The elbow plot (`pca_elbow.png`) — shows variance explained per PC
- PC1 variance (should be <50%; if higher, suspect technical artifact)
- Total variance explained by all 50 PCs

**Elbow plot interpretation:** Look for a "knee" where additional PCs stop adding meaningful variance. Typically 15–30 PCs capture the biology; beyond that is noise. The default 50 is safe — downstream steps (neighbors, clustering) can use fewer.

## UMAP — After Neighbor Graph

UMAP runs AFTER the neighbor graph is built (and optionally after integration). Do not run it before neighbors.

```python
from scagent.tools.embedding import run_umap

result = run_umap(adata, random_state=0)
```

**Critical warning: UMAP is for VISUALIZATION ONLY.**
- Never cluster on UMAP coordinates
- Never interpret distances between clusters as meaningful
- The topology (which clusters are near each other) is informative; metric distances are not

**Show the user:** A UMAP scatter plot colored by cluster labels after Leiden clustering is run.

## Parameters Reference

See `tools/pca.json` and `tools/umap.json`.
