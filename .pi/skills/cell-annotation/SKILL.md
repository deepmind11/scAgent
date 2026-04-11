---
name: cell-annotation
description: Annotate cell types using CellTypist (automated) or marker-gene scoring (manual). Use after clustering and marker gene detection. Both methods are always available.
---

# Cell Type Annotation

Assign cell type identities to clusters. Two methods are always available:

1. **CellTypist** — automated, pre-trained classifier (when a suitable model exists)
2. **Marker-based manual** — score clusters against known marker genes (always works)

Both can be used together: run CellTypist for a first pass, then cross-validate against marker genes.

## Method 1: CellTypist (Automated)

```python
from scagent.tools.annotation import annotate_celltypist

result = annotate_celltypist(adata, model="Immune_All_Low.pkl")
```

**Model selection:**
- `Immune_All_Low.pkl` — broad immune subtypes (best for PBMCs)
- `Immune_All_High.pkl` — fine-grained immune subtypes
- Tissue-specific models available for brain, lung, gut, etc.
- Run `celltypist.models.models_description()` to list all 61 models

**When CellTypist works well:** Human/mouse immune cells, well-studied tissues with a matching model.

**When it does NOT work:** Exotic organisms, novel cell types, tissues without a pre-trained model. In these cases, use marker-based annotation instead.

## Method 2: Marker-Based Manual Annotation

Always available regardless of species, tissue, or organism.

```python
from scagent.tools.annotation import annotate_manual

result = annotate_manual(
    adata,
    marker_dict={
        "CD4 T": ["CD3D", "CD4", "IL7R"],
        "CD8 T": ["CD3D", "CD8A", "CD8B"],
        "B cell": ["CD79A", "MS4A1", "CD19"],
        "NK": ["NKG7", "GNLY", "KLRB1"],
        "CD14 Mono": ["CD14", "LYZ", "S100A9"],
        "FCGR3A Mono": ["FCGR3A", "MS4A7"],
        "DC": ["FCER1A", "CST3"],
        "Platelet": ["PPBP", "PF4"],
    },
    groupby="leiden",
)
```

The function scores each cluster against each cell type and returns ranked suggestions per cluster. The user confirms or overrides.

**Where do markers come from?**
- The agent's knowledge of canonical markers for the tissue
- Published marker databases (CellMarker 2.0, PanglaoDB)
- The user's domain expertise
- The marker gene detection step (cluster markers from Wilcoxon test)

## Applying the Final Annotation

Once the user approves a cluster-to-cell-type mapping:

```python
from scagent.tools.annotation import apply_annotation

result = apply_annotation(
    adata,
    mapping={"0": "CD4 T", "1": "CD14 Mono", "2": "CD4 T", ...},
    groupby="leiden",
    key_added="cell_type",
)
```

## What to Show the User

- Cell type counts table
- UMAP colored by cell type annotations
- For CellTypist: compare predicted types against top marker genes per cluster
- For manual: show the scoring results and ask for confirmation

## Cross-Validation

When using CellTypist, always cross-check against markers:
- Does the CellTypist label match the top marker genes?
- Are there clusters where CellTypist disagrees with the markers?
- Flag any discrepancies and discuss with the user

## Parameters Reference

See `tools/celltypist.json` and `tools/wilcoxon_markers.json`.
