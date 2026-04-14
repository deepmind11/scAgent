---
name: immune-repertoire
description: Analyze TCR/BCR immune receptor data from 10x Chromium V(D)J. Use when the paradigm is immune_repertoire or when V(D)J data is available alongside gene expression.
---

# Skill: Immune Repertoire Analysis

Analyze adaptive immune receptor (TCR/BCR) repertoires from 10x Chromium V(D)J data using Scirpy.

## When to Use

- **Paradigm is `immune_repertoire`** — the DAG includes VDJ loading and clonotype analysis.
- **User has 10x Chromium V(D)J data** (filtered_contig_annotations.csv).
- **Interest in clonotype diversity, clonal expansion, or repertoire overlap.**

## Workflow

### Step 1: Load V(D)J Data

```python
from scagent.tools.repertoire import load_vdj
result = load_vdj(adata, vdj_path="filtered_contig_annotations.csv")
```

### Step 2: Clonotype Analysis

```python
from scagent.tools.repertoire import run_clonotype_analysis
result = run_clonotype_analysis(adata, sequence="aa")
```

## Guard Rails

1. **Check barcode matching** between GEX and VDJ — warn if <10% cells have IR data.
2. **Clonotype definition** depends on sequence type (aa vs nt) and receptor arms.
3. **Diversity metrics** (Shannon, Simpson) require sufficient clonotype counts per group.

## Best-Practice References

- [BP-1] §"Adaptive immune receptor repertoires" (pp. 559-560)
- [BP-2] Ch. 38-39 (Scirpy tutorials)
