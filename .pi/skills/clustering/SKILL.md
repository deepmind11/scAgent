---
name: clustering
description: Cluster cells using the Leiden algorithm. Includes resolution sweep to help choose cluster granularity, and seed sensitivity analysis to verify clustering stability.
---

# Leiden Clustering

Group cells into clusters based on the neighbor graph.

## When to Use

- After the neighbor graph is built
- This is a required step in every analysis

## Step 1: Build the Neighbor Graph

If not already done:

```python
from scagent.tools.neighbors import compute_neighbors

compute_neighbors(adata, n_neighbors=15, metric="cosine", random_state=0)
```

Use `use_rep="X_pca_harmony"` after integration, or `"X_pca"` for single-sample data.

## Step 2: Resolution Sweep

Run Leiden at multiple resolutions to see how cluster count varies:

```python
from scagent.tools.clustering import sweep_resolution

result = sweep_resolution(adata, random_state=0, plot_dir="data/working/plots")
```

**Show the user** the resolution-vs-clusters table and plot. Discuss:
- 0.3–0.5 → major lineages (few broad clusters)
- 0.8–1.2 → standard analysis
- 1.5–3.0 → fine-grained subtypes

Help the user pick based on their biological question. More clusters isn't always better.

## Step 3: Run Leiden

```python
from scagent.tools.clustering import run_leiden

result = run_leiden(adata, resolution=1.0, random_state=0)
```

**Show:** Number of clusters, cluster sizes, any warnings about tiny clusters.

## Step 4: Seed Sensitivity Check

After choosing a resolution, check if the clustering is stable:

```python
from scagent.tools.clustering import check_seed_sensitivity

result = check_seed_sensitivity(adata, resolution=1.0)
```

This runs Leiden with 5 different random seeds and computes the Adjusted Rand Index (ARI) between each pair.

- **ARI ≥ 0.9:** Clustering is stable. Proceed with confidence.
- **ARI < 0.9:** Clustering is fragile at this resolution. The boundaries between some clusters are not well-defined. Options:
  - Try a different resolution
  - Increase/decrease n_neighbors
  - Inspect which clusters are splitting differently

**Always mention the ARI result to the user.** It's important for them to know how reproducible their clusters are.

## Step 5: UMAP Visualization

After clustering, generate UMAP and show it colored by clusters:

```python
from scagent.tools.embedding import run_umap
import scanpy as sc

run_umap(adata, random_state=0)
sc.pl.umap(adata, color="leiden", save="_clusters.png")
```

## What to Do Next

Suggest marker gene detection: "Now let's find the genes that define each cluster — this will help us identify what cell types they are."

## Best Practice Reference

Load `best_practices/reference/clustering.md` for Leiden vs. Louvain rationale, resolution selection guidance, and graph-based vs. K-means comparison.

## Parameters Reference

See `tools/leiden.json` and `tools/neighbors.json`.
