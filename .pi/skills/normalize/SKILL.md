---
name: normalize
description: Normalize and log-transform single-cell count data using CP10K + log1p. Use after QC and doublet detection. Preserves raw counts for downstream differential expression.
---

# Normalization

Normalize raw UMI counts to make cells comparable, then log-transform to reduce skewness.

## When to Use

- After QC filtering and doublet detection
- User asks to normalize or prepare data for clustering

## Prerequisites

- Data must contain raw counts (not already normalized)
- QC filtering and doublet detection should be complete

## How to Run

```python
from scagent.tools.normalize import log_normalize

result = log_normalize(
    adata,
    target_sum=1e4,
    checkpoint_dir="data/working/checkpoints",
)
```

## What Happens Internally

The function performs four operations in a specific order:

1. **Save raw counts** → `adata.layers['counts']` (integer UMIs, needed for pseudobulk DE later)
2. **Normalize** → `sc.pp.normalize_total(adata, target_sum=1e4)` (each cell sums to 10,000 = CP10K)
3. **Log-transform** → `sc.pp.log1p(adata)` (natural log of counts + 1)
4. **Freeze as raw** → `adata.raw = adata.copy()` (saves the full gene set with normalized values)

**This order is critical.** It ensures:
- `adata.layers['counts']` has raw integers for pseudobulk DE (DESeq2, edgeR)
- `adata.raw` has log-normalized values for ALL genes (used by `sc.tl.rank_genes_groups` with `use_raw=True`)
- `adata.X` has log-normalized values (will be subset to HVGs in the next step, but raw keeps all genes)

## What to Show the User

From `result["metrics"]`:
- Max normalized value (should be < 15, typically ~8–12)
- Confirmation that raw counts are preserved
- Confirmation that `adata.raw` is set

Example: "Normalized to CP10K + log1p. Max value: 8.5. Raw counts preserved in `adata.layers['counts']` for downstream DE analysis."

## What NOT to Do

- Do NOT run normalization before Scrublet — Scrublet needs raw counts.
- Do NOT call `adata.raw = adata` again after HVG subsetting — it's already set here with all genes.
- Do NOT use CPM (target_sum=1e6) for scRNA-seq — CP10K is the standard.

## Alternatives

The default CP10K + log1p works well for most analyses. Alternatives exist but are not the default:
- **SCTransform** (Seurat-style variance stabilization) — available via `scvi-tools`
- **scran pooling normalization** — better handles compositional bias
- **Pearson residuals** — analytic normalization, no need for HVG selection

These can be added as tool options later. For now, CP10K + log1p is the recommended default.

## What to Do Next

Suggest feature selection: "Data is normalized. Next step is selecting highly variable genes (HVGs) — these are the genes with the most biological signal, which we'll use for PCA and clustering."

## Parameters Reference

See `tools/log_normalize.json` for the full parameter schema.
