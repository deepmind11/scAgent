---
name: multimodal
description: Analyze CITE-seq (RNA + surface protein) data from 10x Chromium. Use when the paradigm is multimodal or when ADT / antibody-derived tag data is available.
---

# Skill: Multimodal (CITE-seq) Analysis

Analyze joint RNA + surface protein data from 10x Chromium CITE-seq experiments.

## When to Use

- **Paradigm is `multimodal`** — the DAG includes protein loading, normalization, and WNN.
- **User has CITE-seq data** with Antibody-Derived Tag (ADT) counts.
- **Interest in surface protein markers, joint RNA+protein clustering, or protein-informed annotation.**

## Workflow

### Step 1: Load Protein Data

```python
from scagent.tools.multimodal import load_protein
result = load_protein(adata)  # auto-detects from Cell Ranger multi output
```

### Step 2: Normalize Protein (CLR)

```python
from scagent.tools.multimodal import normalize_protein
result = normalize_protein(adata, method="clr")
```

### Step 3: WNN (Joint Graph)

```python
from scagent.tools.multimodal import run_wnn
result = run_wnn(adata, n_neighbors=20, rna_weight=0.5)
```

### Step 4: Cluster + Annotate on WNN graph

Use standard Leiden clustering on the WNN graph, then annotate using both RNA markers and protein markers.

## Guard Rails

1. **Check for isotype controls** in the ADT panel — flag if missing. [BP-2 Ch. 32]
2. **CLR normalization** is the default and most robust. DSB requires empty droplet data.
3. **WNN weight** between RNA and protein should be tuned (default 0.5). More informative modality gets higher weight.
4. **Protein panel size** affects interpretation — small panels (<50 antibodies) limit protein-only analysis.

## Best-Practice References

- [BP-1] §"Surface protein expression" (pp. 558-559)
- [BP-2] Ch. 32-37 (CITE-seq QC, CLR, WNN, annotation)
