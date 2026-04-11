# Chunk 7 Design Plan: Eval Adapter & First Benchmark

**Status:** Draft for review
**Depends on:** Chunks 2 & 3 (QC, clustering, annotation — done)
**Goal:** Get a baseline scBench score on the 7 Chromium canonical evals

---

## 1. How scBench Works (From Source Code Analysis)

### Flow
```
EvalRunner(eval_json)
  → download_data(): latch cp {data_node} → .eval_cache/, symlink to work_dir/
  → agent_function(task_prompt, work_dir) → returns {"answer": {...}}
  → grader.evaluate_answer(answer, config) → passed/failed
```

### Agent Interface
The agent_function receives:
- `task_prompt`: The task text PLUS contextual data about files, enhanced with local file paths
- `work_dir`: Path where data files are symlinked (e.g., `work_dir/157798549.node`)

The agent must:
1. Read the .h5ad file(s) from work_dir
2. Perform the analysis described in the task
3. Return `{"answer": {the_answer_dict}}` — OR write `eval_answer.json` to work_dir

### The 7 Chromium Evals

| Eval | Grader Type | What It Tests | Answer Format |
|------|------------|---------------|---------------|
| **QC** (filter_cells) | numeric_tolerance | Conservative QC → report cells remaining | `{"cells_after_filtering": int}` ±100 of 6420 |
| **Normalization** | numeric_tolerance | LogNorm + z-score → report specific gene value | `{"gene_value": float}` ±0.1 of 2.3 |
| **HVG** | marker_gene_precision_recall | Find 2000 HVGs, must recover canonical markers | `{"top_marker_genes": [...]}` recall ≥50% of 9 markers |
| **Cell typing** | distribution_comparison | Annotate Immune/Epithelial/CAF fractions | `{"cell_type_distribution": {...}}` ±5% per type |
| **Clustering** | numeric_tolerance | Compute pericyte-CAF adjacency in PCA space | `{"pericyte_caf_relative_pca_knn_adjacency": float}` ≥1.1 |
| **DE** | marker_gene_precision_recall | Sub-cluster CAFs, find contractile markers | `{"top_marker_genes": [...]}` recover ≥50% of 9 markers |
| **Trajectory** | marker_gene_precision_recall | Infer trajectory, find terminal markers | `{"top_marker_genes": [...]}` recover Ly6c1 AND Acta2 |

### Data: 5 Unique Datasets (All 4T1 Mouse Tumor)
- `157798549.node` — QC eval (6,420 cells, already filtered)
- `157798550.node` — Normalization eval
- `157798551.node` — HVG + cell typing evals
- `157798553.node` — Trajectory eval (CAF subset)
- `158379062.node` — Clustering + DE evals (with cell_type labels)

All are **mouse** (gene names like Cd74, Epcam, Ly6c1). Our tools handle this (species detection from gene case).

---

## 2. Architecture Decision: Direct Python vs. CLI Agent

Two options for the adapter:

### Option A: Launch `pi` CLI (like Claude Code)
```python
def agent_fn(task_prompt, work_dir):
    subprocess.run(["pi", "--print", ...], stdin=task_prompt, cwd=work_dir)
    return json.load(open(work_dir / "eval_answer.json"))
```
- **Pro:** Tests the real agent end-to-end (LLM reads skills, writes code, etc.)
- **Con:** Expensive (each eval costs ~$1-5 in API calls), slow (~2-5 min per eval), non-deterministic (LLM may make mistakes)
- **Con:** We'd need scAgent running as a pi agent, which requires API credits

### Option B: Direct Python adapter (no LLM)
```python
def agent_fn(task_prompt, work_dir):
    # Parse what the task asks for
    # Call our scagent.tools functions directly
    # Return the answer
    ...
```
- **Pro:** Free, fast, deterministic, testable
- **Pro:** Measures whether our TOOLS are correct, independent of LLM reliability
- **Con:** Doesn't test the agent's ability to interpret the task

### Decision: **Both, in this chunk.**

1. **Phase 1: Direct adapter (Option B)** — Tests our tool functions against the benchmark. This tells us: "If our tools are called correctly, do they produce the right answers?" Fast, free, immediate feedback. Also reveals any tool bugs before we spend API credits.

2. **Phase 2: LLM agent adapter (Option A)** — Tests the full agent end-to-end via `pi --print`. This tells us: "Does the agent read the task, choose the right tools, and produce the right answer?" Costs API credits (~$1-5 per eval) but is the real benchmark.

We run both and compare: the gap between Phase 1 and Phase 2 tells us how much the LLM is helping vs. hurting.

---

## 3. Direct Adapter Design

### `scbench_adapter.py`

Each eval type needs a handler function that:
1. Reads the .h5ad from work_dir
2. Calls our scagent.tools functions
3. Returns the answer in the expected format

```python
def handle_qc(adata, task_prompt):
    """QC: filter cells, return count."""
    calculate_qc_metrics(adata)
    adata, result = filter_cells(adata, ...)
    return {"cells_after_filtering": result["metrics"]["cells_after"]}

def handle_normalization(adata, task_prompt):
    """Norm: log-normalize + z-score, return specific gene value."""
    # This one is tricky: task asks for z-score scaling too (sc.pp.scale)
    # AND asks for the value of a specific gene in a specific cell
    ...

def handle_hvg(adata, task_prompt):
    """HVG: select 2000 HVGs, return gene list."""
    log_normalize(adata)
    select_hvg(adata, n_top_genes=2000)
    hvg_names = adata.var_names[adata.var["highly_variable"]].tolist()
    return {"top_marker_genes": hvg_names}

def handle_celltyping(adata, task_prompt):
    """Cell typing: annotate compartments, return fractions."""
    # Needs manual annotation with specific markers for
    # Immune, Epithelial/Cancer, CAF compartments
    ...

def handle_clustering(adata, task_prompt):
    """Clustering: compute pericyte-CAF adjacency ratio."""
    # Data already has cell_type labels
    # Need to compute KNN in PCA, measure adjacency
    ...

def handle_de(adata, task_prompt):
    """DE: sub-cluster CAFs, find contractile markers."""
    # Filter to CAFs, cluster into 6, find contractile sub-cluster
    # Run DE on that vs. other CAFs
    ...

def handle_trajectory(adata, task_prompt):
    """Trajectory: infer trajectory, find terminal markers."""
    # Need trajectory inference (PAGA or diffusion pseudotime)
    # Find genes distinguishing the two ends
    ...
```

### Task Routing

The adapter needs to figure out which handler to call. Since we know the 7 eval IDs, we can map directly:

```python
EVAL_HANDLERS = {
    "chromium_qc_4T1_filter_cells": handle_qc,
    "chromium_4t1_normalization": handle_normalization,
    "chromium_4t1_hvg_gene_sets": handle_hvg,
    "chromium_celltyping_01_4t1_compartment_fractions": handle_celltyping,
    "chromium_clustering_01_4t1_pericyte_adjacent_to_caf": handle_clustering,
    "chromium_differential_expression_01_contractile_caf_marker_recovery": handle_de,
    "chromium_trajectory_01_caf_terminal_marker_recovery": handle_trajectory,
}
```

But that's overfitting to the eval set. Better: parse the task prompt to determine what kind of analysis to run. Still, for the 7 canonical evals, we know exactly what's asked.

**Pragmatic decision:** Use eval ID routing for the direct adapter (Option B). The LLM adapter (Option A, later) will parse the prompt naturally.

---

## 4. Per-Eval Analysis

### 4.1 QC (`chromium_qc_4T1_filter_cells`)

**Task:** "Conservative quality-control thresholds to remove only clearly low-quality cells."
**Ground truth:** 6420 cells (tolerance ±100, so 6320–6520 passes)
**Note from eval:** "The data are already post cell-filtering and expected to be high quality; a reasonable QC procedure should remove few cells."

**Strategy:** This data is ALREADY clean. We need very conservative thresholds.
- Run `calculate_qc_metrics()` to see the distributions
- Use very permissive thresholds (data is already filtered)
- The key insight: "conservative" means DON'T remove many cells

**Risk:** Our default thresholds might remove too many cells. Need to adapt to the data.

### 4.2 Normalization (`chromium_4t1_normalization`)

**Task:** "Log normalization and z-score scale. Find cell with highest raw count for 'Mrc1'. Report final normalized value."
**Ground truth:** 2.3 (±0.1)

**Strategy:**
1. Save raw counts
2. Find cell with max raw Mrc1 expression
3. `sc.pp.normalize_total(adata, target_sum=1e4)`
4. `sc.pp.log1p(adata)`
5. `sc.pp.scale(adata)` — **NOTE: our normalize skill doesn't do z-score scaling!**
6. Report the Mrc1 value for that specific cell

**Gap:** Our `log_normalize()` doesn't include `sc.pp.scale()`. The eval explicitly asks for z-score scaling. Need to handle this.

### 4.3 HVG (`chromium_4t1_hvg_gene_sets`)

**Task:** "Identify 2000 highly variable genes."
**Ground truth:** Must recover ≥50% of 9 canonical markers (Cd74, Crabp1, Epcam, Thy1, Pdpn, Pdgfra, Mcam, Rgs5 — note mouse gene casing)
**Grader:** precision=0 (any extra genes fine), recall must be ≥0.5

**Strategy:**
1. `log_normalize(adata)`
2. `select_hvg(adata, n_top_genes=2000)`
3. Return all HVG names

**Risk:** Low — these canonical markers should naturally be in the top 2000 HVGs for a tumor dataset. But we need to handle mouse gene casing.

### 4.4 Cell Typing (`chromium_celltyping_01_4t1_compartment_fractions`)

**Task:** Classify into Immune, Epithelial/Cancer, CAF. Report percentages.
**Ground truth:** Immune=66.4%, Epithelial/Cancer=24.5%, CAF=8.0% (±5% each)

**Strategy:** Use `annotate_manual()` with marker genes for these compartments:
- Immune: Ptprc (Cd45), Cd3d, Cd79a, etc.
- Epithelial/Cancer: Epcam, Krt8, Krt18
- CAF: Pdgfra, Thy1, Col1a1, Dcn

Then compute percentages. This is a MOUSE dataset so gene names are title-case.

### 4.5 Clustering (`chromium_clustering_01_4t1_pericyte_adjacent_to_caf`)

**Task:** Compute pericyte-CAF adjacency in PCA KNN space.
**Ground truth:** ≥1.1 (the ratio of CAF adjacency to max other adjacency for pericytes)
**Data already has:** cell_type labels in obs

**Strategy:**
1. Standard preprocessing
2. PCA
3. KNN in PCA space
4. For pericyte cells, count edges to each cell type
5. Compute CAF fraction / max other fraction

This is a custom computation — no existing wrapper. Need to write the logic.

### 4.6 DE (`chromium_differential_expression_01_contractile_caf_marker_recovery`)

**Task:** Cluster CAFs into 6 subclusters, find contractile subcluster markers.
**Ground truth:** Must recover ≥50% of 9 contractile markers (Tpm1, Tpm2, Myl9, Tagln, Cnn2, Cnn3, Igfbp3, Tnc, Tmem119)

**Strategy:**
1. Filter to CAFs only (using cell_type labels)
2. Standard preprocessing on CAF subset
3. Leiden clustering at a resolution giving ~6 clusters
4. Find which cluster has highest expression of contractile genes
5. Run `find_marker_genes()` on that cluster vs. rest
6. Return top 50 genes

### 4.7 Trajectory (`chromium_trajectory_01_caf_terminal_marker_recovery`)

**Task:** Infer trajectory, find terminal markers (Ly6c1 and Acta2).
**Ground truth:** Must recover BOTH Ly6c1 AND Acta2 (recall=1.0 at k=2)

**Strategy:**
1. Standard preprocessing
2. Diffusion pseudotime or PAGA
3. Identify two terminal groups
4. DE between them
5. Return top markers

**Gap:** We don't have trajectory wrappers implemented yet (Chunk 11). Options:
- Write a lightweight diffusion pseudotime handler using `sc.tl.diffmap` + `sc.tl.dpt`
- Use PAGA (already in tool registry but no wrapper)
- Accept this eval may fail for now

---

## 5. Gaps to Fill

| Gap | Needed for | Solution |
|-----|-----------|----------|
| `sc.pp.scale()` (z-score) | Normalization eval + PCA correctness | Fix in our pipeline: add `scale()` step between HVG and PCA. This IS standard Scanpy — it fell through the gap between `log_normalize()` and `run_pca()`. |
| Pericyte adjacency computation | Clustering eval | Custom code in handler — KNN neighbor type composition analysis |
| Trajectory inference | Trajectory eval | Lightweight `sc.tl.diffmap` + `sc.tl.dpt` handler |
| Mouse gene casing | All evals (4T1 is mouse) | Already handled — species auto-detection works |
| LLM agent adapter | Phase 2 | Wire `pi --print` with scAgent system prompt as agent_function |

---

## 6. File List

| File | Purpose |
|------|---------|
| `eval/adapter.py` | Main adapter: EvalRunner → handler dispatch → answer |
| `eval/handlers/qc.py` | QC eval handler |
| `eval/handlers/normalization.py` | Normalization eval handler |
| `eval/handlers/hvg.py` | HVG eval handler |
| `eval/handlers/celltyping.py` | Cell typing eval handler |
| `eval/handlers/clustering.py` | Clustering/adjacency eval handler |
| `eval/handlers/de.py` | Differential expression eval handler |
| `eval/handlers/trajectory.py` | Trajectory eval handler |
| `eval/run_benchmark.py` | Script to run all 7 evals and report results |

---

## 7. Expected Results (Honest Prediction)

| Eval | Can We Pass? | Confidence | Notes |
|------|-------------|-----------|-------|
| QC | **Likely yes** | 80% | Need conservative thresholds for already-clean data |
| Normalization | **Likely yes** | 85% | Straightforward computation, just need to add scale() |
| HVG | **Likely yes** | 90% | Standard HVG selection, canonical markers should be in top 2000 |
| Cell typing | **Maybe** | 60% | Manual marker annotation on mouse tumor — need right markers |
| Clustering | **Maybe** | 50% | Custom adjacency computation, less standard |
| DE | **Maybe** | 55% | Need to find the right CAF subcluster |
| Trajectory | **Unlikely** | 25% | No trajectory wrapper yet |

**Predicted baseline: 3-4 out of 7 (43-57%)**

This is a realistic first score. The best existing agent (per scBench paper) scores 52.8% overall. If we hit 43-57% on our first attempt with a direct adapter (no LLM), that's a strong foundation.

---

## 8. Build Order

1. Set up `eval/` directory structure
2. Build the adapter framework (loader, dispatcher)
3. Implement handlers one at a time, testing each against its eval
4. Run all 7, record results
5. Analyze failures
6. Commit results + analysis

---

## 9. Important: Not Overfitting

The direct adapter (Option B) knows which eval it's running via the eval ID. This is arguably "overfitting" — but it's actually the right approach for a BASELINE:

- We're testing whether our **tools produce correct results**, not whether an LLM can interpret the task
- The handlers use our standard tool functions with reasonable parameters
- We're NOT tuning parameters to match the ground truth — we're using principled defaults
- The LLM adapter (Option A, later) will be the true test of the agent

The distinction: the adapter routes to the right handler, but the handler uses our tools with the same defaults a biologist would choose.
