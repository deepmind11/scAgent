---
name: doublet-detection
description: Detect and optionally remove cell doublets using Scrublet. Use after QC filtering and before normalization. Must run on raw (unnormalized) counts.
---

# Doublet Detection

Identify cell doublets (two cells captured in one droplet) using Scrublet via Scanpy's built-in wrapper.

## When to Use

- After QC filtering, before normalization
- User asks about doublets or data quality
- Part of the standard preprocessing pipeline

## Prerequisites

- Data must contain raw counts (NOT normalized). If `adata.X.max() > 100`, counts are likely raw. If `max < 15`, data may already be log-normalized — Scrublet will not work correctly.
- QC filtering should be done first (remove low-quality cells before doublet detection).

## How to Run

```python
from scagent.tools.doublets import detect_doublets

adata, result = detect_doublets(
    adata,
    random_state=0,
    plot_dir="data/working/plots",
    checkpoint_dir="data/working/checkpoints",
)
```

The function auto-calculates the expected doublet rate from the cell count:
- 10x Chromium rate: ~0.8% per 1,000 cells captured
- Example: 10,000 cells → ~8% expected doublets

## What to Show the User

1. **Doublet score histogram** — saved to `plot_dir/doublet_histogram.png`. Show it.
2. **Key numbers** from `result["metrics"]`:
   - Number of doublets detected and rate
   - Expected rate (for comparison)
   - Whether doublets were removed or just annotated

Example: "Detected 545 doublets (5.0% of cells). Expected rate for 10,991 cells: ~8.8%. Removed doublets, leaving 10,446 cells."

## Interpreting Results

- **Rate close to expected:** Normal. Proceed.
- **Rate much higher than expected (>2×):** Scrublet may be over-calling. Inspect the histogram — is there a clear bimodal separation? If not, the threshold may need manual adjustment.
- **Rate much lower than expected (<1/3):** Scrublet may have failed. The score distribution may be unimodal (no clear doublet peak). Consider adjusting `n_prin_comps`.

## Options

- `remove=True` (default): Remove predicted doublets from the dataset.
- `remove=False`: Only annotate — adds `doublet_score` and `predicted_doublet` to `adata.obs` without removing any cells. Useful when the user wants to inspect doublets before deciding.

## Reproducibility

The `random_state` parameter (default 0) ensures identical results on re-run. Same data + same seed = same doublets, always.

## What to Do Next

Suggest normalization: "Doublet detection is complete. Next I'll normalize the data — this converts raw counts to log-transformed values for downstream analysis."

## Parameters Reference

See `tools/scrublet_doublets.json` for the full parameter schema.

## Notes

- For multi-sample experiments, ideally run Scrublet per sample before merging. The current function runs on the full dataset.
- Scrublet performs its own internal normalization and PCA — this is separate from the analysis pipeline's normalization step.
