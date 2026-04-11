---
name: marker-genes
description: Find marker genes that define each cluster using differential expression (Wilcoxon rank-sum test). Use after clustering to characterize clusters before cell type annotation.
---

# Marker Gene Detection

Find genes that are significantly upregulated in each cluster compared to all other cells. These markers characterize each cluster and are the basis for cell type annotation.

## When to Use

- After Leiden clustering
- Before cell type annotation (provides evidence for annotation)
- When the user asks "what genes define cluster X?"

## How to Run

```python
from scagent.tools.markers import find_marker_genes

result = find_marker_genes(
    adata,
    groupby="leiden",
    method="wilcoxon",
    plot_dir="data/working/plots",
)
```

## What to Show the User

From `result["metrics"]["top_markers"]` — a table of top 5 markers per cluster:

| Cluster | Top markers | Scores |
|---------|------------|--------|
| 0 | IL7R, CD3D, IL32 | 22.1, 19.8, 18.5 |
| 1 | LYZ, S100A9, S100A8 | 35.2, 32.1, 30.8 |
| ... | ... | ... |

Also show the dotplot (`marker_dotplot.png`) if generated.

## Interpreting Markers

The agent should recognize common markers for well-known cell types:
- **T cells:** CD3D, CD3E, TRAC
- **CD4 T:** CD4, IL7R, CCR7 (naive), IL2RA (Treg)
- **CD8 T:** CD8A, CD8B, GZMK (memory), GZMB (effector)
- **B cells:** CD79A, MS4A1, CD19
- **NK cells:** NKG7, GNLY, KLRB1
- **CD14 Monocytes:** CD14, LYZ, S100A9, S100A8
- **FCGR3A Monocytes:** FCGR3A, MS4A7
- **Dendritic cells:** FCER1A, CST3
- **Platelets:** PPBP, PF4

Use these to suggest cluster identities, but always present as suggestions — the user makes the final call.

## Warning Signs

- **Clusters with no clear markers** (<3 significant genes): May be over-clustered. Suggest merging with a neighboring cluster or lowering resolution.
- **Multiple clusters with the same top markers**: May be split versions of the same cell type. Suggest lowering resolution.
- **Ribosomal or mitochondrial genes as top markers**: Likely technical artifact, not biological signal.

## What to Do Next

Suggest cell type annotation: "Based on these markers, I can suggest cell type identities. Would you like me to run automated annotation with CellTypist, or shall we annotate manually using marker genes?"

## Important Note

This finds **cluster markers** — genes distinguishing one cluster from all others. For **cross-condition differential expression** (disease vs healthy), use pseudobulk methods (DESeq2/edgeR) in a later step. Do not confuse the two.

## Parameters Reference

See `tools/wilcoxon_markers.json`.
