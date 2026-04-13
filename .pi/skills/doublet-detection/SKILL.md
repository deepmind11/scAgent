---
name: doublet-detection
description: Detect and optionally remove cell doublets using Scanpy's built-in scrublet wrapper (sc.pp.scrublet). Use after QC filtering and before normalization. Must run on raw (unnormalized) counts.
---

# Doublet Detection

Identify cell doublets (two cells captured in one droplet) using `sc.pp.scrublet()` — Scanpy's built-in reimplementation of the Scrublet algorithm.

> **⚠️ WARNING: Do NOT `import scrublet` directly.**
> The `scrublet` PyPI package contains native C extensions that cause **Bus Error (SIGBUS) crashes on Python 3.13** due to incompatible compiled code. Always use `sc.pp.scrublet()` from Scanpy instead — it is a pure-Python reimplementation of the same algorithm with an identical interface and does not depend on the broken native library.

## When to Use

- After QC filtering, before normalization
- User asks about doublets or data quality
- Part of the standard preprocessing pipeline

## Prerequisites

- Data must contain raw counts (NOT normalized). If `adata.X.max() > 100`, counts are likely raw. If `max < 15`, data may already be log-normalized — Scrublet will not work correctly.
- QC filtering should be done first (remove low-quality cells before doublet detection).

## How to Run

### Single-sample experiment

```python
import scanpy as sc

# Calculate expected doublet rate from cell count
# 10x Chromium rate: ~0.8% per 1,000 cells captured
n_cells = adata.n_obs
expected_doublet_rate = 0.008 * (n_cells / 1000)

sc.pp.scrublet(adata, expected_doublet_rate=expected_doublet_rate, random_state=0)
```

This adds two columns to `adata.obs` in-place:
- `predicted_doublet` — `bool`, True for predicted doublets
- `doublet_score` — `float`, continuous doublet score (0–1)

### Multi-sample experiment (required)

For datasets with multiple samples, **run scrublet per sample** then concatenate. Running on the merged dataset confuses the algorithm because inter-sample variation looks like doublets.

```python
import scanpy as sc

sample_adatas = []
for sample_id in adata.obs["sample"].unique():
    adata_sub = adata[adata.obs["sample"] == sample_id].copy()
    n_cells = adata_sub.n_obs
    expected_doublet_rate = 0.008 * (n_cells / 1000)
    sc.pp.scrublet(adata_sub, expected_doublet_rate=expected_doublet_rate, random_state=0)
    sample_adatas.append(adata_sub)

adata = sc.concat(sample_adatas)
```

### Removing doublets

After running `sc.pp.scrublet`, filter doublets explicitly:

```python
n_before = adata.n_obs
adata = adata[~adata.obs["predicted_doublet"]].copy()
n_removed = n_before - adata.n_obs
```

To annotate only (without removing), simply skip the filtering step.

## What to Show the User

1. **Doublet score histogram** — plot with `sc.pl.scrublet_score_distribution(adata)` or build a custom histogram from `adata.obs['doublet_score']`. Save to `data/working/plots/doublet_histogram.png` and show it.
2. **Key numbers** — compute directly from `adata.obs`:
   - `n_doublets = adata.obs['predicted_doublet'].sum()`
   - `doublet_rate = n_doublets / adata.n_obs`
   - Expected rate (calculated from cell count as above)
   - Whether doublets were removed or just annotated

Example: "Detected 545 doublets (5.0% of cells). Expected rate for 10,991 cells: ~8.8%. Removed doublets, leaving 10,446 cells."

## Interpreting Results

- **Rate close to expected:** Normal. Proceed.
- **Rate much higher than expected (>2×):** Scrublet may be over-calling. Inspect the histogram — is there a clear bimodal separation? If not, the threshold may need manual adjustment.
- **Rate much lower than expected (<1/3):** Scrublet may have failed. The score distribution may be unimodal (no clear doublet peak). Consider adjusting `n_prin_comps`.

## Reproducibility

The `random_state` parameter (default 0) ensures identical results on re-run. Same data + same seed = same doublets, always.

## What to Do Next

Suggest normalization: "Doublet detection is complete. Next I'll normalize the data — this converts raw counts to log-transformed values for downstream analysis."

## Parameters Reference

See `tools/scrublet_doublets.json` for the full parameter schema.

## Best Practice Reference

Load `best_practices/reference/doublet-detection.md` for benchmark-backed method recommendations and interpretation guidance.

## Notes

- **Never `import scrublet`** — see warning at top. Use `sc.pp.scrublet()` exclusively.
- For multi-sample experiments, always run scrublet per sample before merging (see multi-sample code above).
- `sc.pp.scrublet` performs its own internal normalization and PCA — this is separate from the analysis pipeline's normalization step.
