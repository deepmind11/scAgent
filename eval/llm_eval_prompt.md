# scAgent LLM Benchmark Prompt

Paste this into a fresh Claude Code session launched from `/Users/hgz/Projects/scAgent/`.

---

## Prompt

You are scAgent, a single-cell RNA-seq analysis assistant. You have Python tools available in the `scagent.tools` package.

**Activate the venv first:** `source /Users/hgz/Projects/scAgent/.venv/bin/activate`

### Available tools (import from `scagent.tools.*`):

```python
# Data loading
from scagent.tools.loading import load_10x_h5  
# load_10x_h5(filename) → (adata, result_dict)

# QC
from scagent.tools.qc import calculate_qc_metrics, filter_cells, filter_genes
# calculate_qc_metrics(adata) → result with metrics["recommended_thresholds"]
# filter_cells(adata, min_genes=200, max_genes=5000, max_pct_mito=10.0) → (adata, result)
# filter_genes(adata, min_cells=3) → (adata, result)

# Normalization
from scagent.tools.normalize import log_normalize
# log_normalize(adata, target_sum=1e4) → result (modifies in place, saves raw to layers['counts'])

# Feature selection
from scagent.tools.feature_selection import select_hvg
# select_hvg(adata, n_top_genes=2000) → result (marks adata.var['highly_variable'])

# PCA (includes z-score scaling)
from scagent.tools.pca import run_pca
# run_pca(adata, n_comps=50, random_state=0) → result

# Neighbors
from scagent.tools.neighbors import compute_neighbors
# compute_neighbors(adata, n_neighbors=15, metric="cosine", random_state=0) → result

# UMAP
from scagent.tools.embedding import run_umap
# run_umap(adata, random_state=0) → result

# Clustering
from scagent.tools.clustering import run_leiden, sweep_resolution
# run_leiden(adata, resolution=1.0, random_state=0) → result
# sweep_resolution(adata) → result with recommended_resolution

# Markers
from scagent.tools.markers import find_marker_genes
# find_marker_genes(adata, groupby="leiden", method="wilcoxon", use_raw=True) → result

# Annotation
from scagent.tools.annotation import annotate_celltypist, annotate_manual, apply_annotation

# You can also use scanpy directly: import scanpy as sc
```

### Important notes:
- Some datasets have a `raw_counts` layer — check `adata.layers` after loading. If present, you may need to set `adata.X = adata.layers['raw_counts'].copy()` before normalizing
- Check `adata.X.max()` to determine if data is raw counts (>100) or already normalized (<15)
- Species is auto-detected from gene names (uppercase=human, title-case=mouse)
- For QC: `calculate_qc_metrics()` returns MAD-based recommended thresholds — use them for conservative filtering
- z-score scaling is included in `run_pca()` — don't call `sc.pp.scale()` separately unless you're NOT using `run_pca()`
- For trajectory analysis, use `sc.tl.diffmap()` for diffusion pseudotime

---

### Tasks

Run ALL 7 scBench evals below. For each one:
1. Load the data, inspect it (shape, layers, obs columns, X.max)
2. Write and run a Python script to solve the task
3. Save the answer as JSON to `eval/results/llm_{eval_id}.json`
4. Then run the grader to check if it passes

**Grading command** (run after each eval, replacing EVAL_ID):
```bash
cd /Users/hgz/Projects/scAgent && source .venv/bin/activate
python3 -c "
import json, sys
sys.path.insert(0, '/Users/hgz/Projects/scbench-eval')
from latch_eval_tools.graders import GRADER_REGISTRY
answer = json.load(open('eval/results/llm_EVAL_ID.json'))
grader_config = json.load(open('eval/results/grader_EVAL_ID.json'))
grader = GRADER_REGISTRY[grader_config['type']]()
result = grader.evaluate_answer(answer, grader_config['config'])
print(f'PASSED: {result.passed}')
print(result.reasoning)
"
```

---

### Eval 1: QC (`chromium_qc_4T1_filter_cells`)
**Data:** `/Users/hgz/Projects/scAgent/.eval_cache/cache/332388f86e1d8227__157798549.node`
**Task:** Use the given .h5ad containing all viable cells from 4T1 tumors. Choose conservative quality-control thresholds to remove only clearly low-quality cells and apply the filtering.
**Answer format:** `{"cells_after_filtering": <int>}`
**Save to:** `eval/results/llm_chromium_qc_4T1_filter_cells.json`
**Grader config:** `eval/results/grader_chromium_qc_4T1_filter_cells.json`

### Eval 2: Normalization (`chromium_4t1_normalization`)
**Data:** `/Users/hgz/Projects/scAgent/.eval_cache/cache/bca9e0cc52adc8c7__157798550.node`
**Task:** Apply global scale log normalization and z-score scale to the given .h5ad. Identify the cell with the highest RAW (pre-normalization) count for the gene 'Mrc1'. For that cell, report the final normalized value of 'Mrc1'.
**Answer format:** `{"gene_value": <float>}`
**Save to:** `eval/results/llm_chromium_4t1_normalization.json`
**Grader config:** `eval/results/grader_chromium_4t1_normalization.json`

### Eval 3: HVG (`chromium_4t1_hvg_gene_sets`)
**Data:** `/Users/hgz/Projects/scAgent/.eval_cache/cache/e18a87d8c3f0e75d__157798551.node`
**Task:** Identify 2000 highly variable genes in given .h5ad, report all highly variable genes.
**Answer format:** `{"top_marker_genes": ["Gene1", "Gene2", ...]}`
**Save to:** `eval/results/llm_chromium_4t1_hvg_gene_sets.json`
**Grader config:** `eval/results/grader_chromium_4t1_hvg_gene_sets.json`

### Eval 4: Cell typing (`chromium_celltyping_01_4t1_compartment_fractions`)
**Data:** `/Users/hgz/Projects/scAgent/.eval_cache/cache/e18a87d8c3f0e75d__157798551.node`
**Task:** Annotate the major cell populations using canonical marker genes. Classify all cells into three compartments (exact strings): Immune, Epithelial/Cancer, CAF. Treat gene names case-insensitively. Report cell_type_distribution as percentage of total cells per compartment (should sum to ~100%).
**Answer format:** `{"cell_type_distribution": {"Immune": <pct>, "Epithelial/Cancer": <pct>, "CAF": <pct>}}`
**Save to:** `eval/results/llm_chromium_celltyping_01_4t1_compartment_fractions.json`
**Grader config:** `eval/results/grader_chromium_celltyping_01_4t1_compartment_fractions.json`

### Eval 5: Clustering (`chromium_clustering_01_4t1_pericyte_adjacent_to_caf`)
**Data:** `/Users/hgz/Projects/scAgent/.eval_cache/cache/a3d3d7ddb0e39418__158379062.node`
**Task:** Using the provided dataset with cell_type labels in .obs['cell_type'], preprocess and compute PCA. Build KNN in PCA space. For pericyte cells, compute fraction of KNN edges to each non-pericyte type. Compute: CAF adjacency fraction / max(other non-pericyte fractions). Report this ratio.
**Answer format:** `{"pericyte_caf_relative_pca_knn_adjacency": <float>}`
**Save to:** `eval/results/llm_chromium_clustering_01_4t1_pericyte_adjacent_to_caf.json`
**Grader config:** `eval/results/grader_chromium_clustering_01_4t1_pericyte_adjacent_to_caf.json`

### Eval 6: DE (`chromium_differential_expression_01_contractile_caf_marker_recovery`)
**Data:** `/Users/hgz/Projects/scAgent/.eval_cache/cache/a3d3d7ddb0e39418__158379062.node`
**Task:** Using the dataset with cell_type labels in .obs['cell_type']. Cluster CAFs into 6 subclusters. Identify the CAF subcluster with the strongest contractile/cytoskeletal program. Perform DE comparing this subcluster against all other CAFs subtypes. Report top 50 marker genes.
**Answer format:** `{"top_marker_genes": ["Gene1", ...]}`
**Save to:** `eval/results/llm_chromium_differential_expression_01_contractile_caf_marker_recovery.json`
**Grader config:** `eval/results/grader_chromium_differential_expression_01_contractile_caf_marker_recovery.json`

### Eval 7: Trajectory (`chromium_trajectory_01_caf_terminal_marker_recovery`)
**Data:** `/Users/hgz/Projects/scAgent/.eval_cache/cache/243ae62a228e9a81__157798553.node`
**Task:** Perform expression-based preprocessing and infer the dominant trajectory capturing the major axis of CAF heterogeneity. Identify two terminal cell groups (the two ends of the dominant trajectory). Compute genes that most strongly distinguish these two terminal groups. Report top marker genes.
**Answer format:** `{"top_marker_genes": ["Gene1", ...]}`
**Save to:** `eval/results/llm_chromium_trajectory_01_caf_terminal_marker_recovery.json`
**Grader config:** `eval/results/grader_chromium_trajectory_01_caf_terminal_marker_recovery.json`

---

After running all 7, summarize results in a table:
```
| Eval | Passed? | Details |
|------|---------|---------|
| QC | ? | ? |
| Normalization | ? | ? |
| HVG | ? | ? |
| Cell typing | ? | ? |
| Clustering | ? | ? |
| DE | ? | ? |
| Trajectory | ? | ? |

Score: X/7
```
