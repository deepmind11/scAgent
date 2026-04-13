# Plan: AnnData Inspector + Architectural Shift

**Status:** Draft for review  
**Depends on:** All existing chunks (tools, provenance, state, context, DAG, knowledge)  
**Changes:** New module, DAG refactor, architecture doc update, eval handler fixes

---

## 0. What This Plan Does

Three connected changes:

1. **New module: `scagent/inspector.py`** — inspects foreign AnnData to determine what preprocessing has been done and what's available to work with.
2. **DAG refactor** — shift from "mandatory paradigm-first pipeline" to "optional plan + always-on dependency graph." The DAG becomes advisory, not blocking.
3. **Architecture doc update** — reflect the new design philosophy: question-driven analysis, not paradigm-driven.

---

## 1. The Core Problem

scAgent currently assumes it controls the data from the start. When it encounters foreign AnnData (scBench evals, a researcher's pre-processed data), it has no way to know:

- What state `adata.X` is in (raw counts? log-normalized? z-score scaled?)
- Where raw counts are stored (if anywhere)
- Which analysis steps have already been completed
- What metadata is available (conditions, batches, cell types)

This caused the HVG eval crash: the handler checked `X.max() > 50`, concluded "raw counts," and ran `log_normalize()` on z-score scaled data.

---

## 2. Design: `scagent/inspector.py`

### 2.1 `inspect_adata(adata) → dict`

Single entry point. Returns a comprehensive state report. Runs once when data is loaded.

**X matrix state detection** (priority-ordered):

```
1. (X < 0).any()                          → 'scaled'
2. np.allclose(sample, round(sample))     → 'raw_counts'  
3. 'log1p' in adata.uns                   → 'log_normalized'
4. max_val < 20                           → 'log_normalized' (heuristic)
5. else                                   → 'transformed_unknown'
```

**Completed steps detection** (breadcrumb checks):

| Step | Breadcrumb | Key to check |
|------|-----------|--------------|
| QC metrics | `obs['n_genes_by_counts']` | Exact column name from Scanpy |
| Normalization | `uns['log1p']` | Written by `sc.pp.log1p()` |
| Raw frozen | `adata.raw is not None` | Standard Scanpy convention |
| HVG selection | `var['highly_variable']` | Exact column name |
| PCA | `obsm['X_pca']` | Exact key |
| Neighbors | `uns['neighbors']` + `obsp['connectivities']` | Both must exist |
| UMAP | `obsm['X_umap']` | Exact key |
| tSNE | `obsm['X_tsne']` | Exact key |
| Clustering | `obs` columns starting with `leiden`/`louvain`/`cluster` | Fuzzy match |
| Cell types | `obs` columns matching `cell_type`/`celltype`/`annotation`/etc. | Fuzzy match |
| DE results | `uns['rank_genes_groups']` | Exact key |
| Diffusion map | `obsm['X_diffmap']` | Exact key |
| RNA velocity | `layers['spliced']` + `layers['unspliced']` | Both must exist |

**Raw counts location** (search order):

```
1. layers['counts']
2. layers['raw_counts']  
3. X itself (if x_state == 'raw_counts')
4. None → flag as missing
```

**Metadata detection** (fuzzy obs column matching):

| Looking for | Column name candidates | Extra validation |
|-------------|----------------------|------------------|
| Cluster key | `leiden*`, `louvain*`, `cluster*` | Categorical, 2–100 unique values |
| Cell type key | `cell_type`, `celltype`, `cell_label`, `annotation`, `celltypist*`, `azimuth*` | Categorical |
| Condition key | `condition`, `disease`, `treatment`, `group`, `status` | ≥2 unique values |
| Batch key | `batch`, `sample`, `donor`, `patient`, `library_id`, `channel` | ≥2 unique values |
| Sample key | `sample`, `sample_id`, `library_id` | ≥1 unique values |

**Species detection** (existing logic from `context.py`, moved here):

Gene name casing: >70% uppercase → human, <30% → mouse. Also check for `mt-` vs `MT-` prefix.

### 2.2 Return format

```python
@dataclass
class AnnDataState:
    # Matrix state
    x_state: str               # 'raw_counts' | 'log_normalized' | 'scaled' | 'transformed_unknown'
    raw_counts_location: str | None  # "layers['counts']", "X", or None
    
    # Completed steps (booleans)
    has_qc_metrics: bool
    has_normalized: bool
    has_raw_frozen: bool
    has_hvg: bool
    has_pca: bool
    has_neighbors: bool
    has_umap: bool
    has_tsne: bool
    has_clusters: bool
    has_cell_types: bool
    has_de_results: bool
    has_diffmap: bool
    has_velocity_layers: bool
    
    # Detected keys (actual column/key names found)
    cluster_key: str | None     # e.g., 'leiden', 'leiden_0.8'
    celltype_key: str | None    # e.g., 'cell_type'
    condition_key: str | None   # e.g., 'condition'
    batch_key: str | None       # e.g., 'sample'
    
    # Data shape
    n_cells: int
    n_genes: int
    species: str | None         # 'human', 'mouse', or None
    
    # Warnings
    warnings: list[str]         # e.g., "X is z-score scaled but no raw counts found"
```

### 2.3 `find_raw_counts(adata, state) → AnnData or matrix`

Given the inspector output, retrieve the best available raw counts for re-processing:

```
1. If raw_counts_location is not None → return that layer/X
2. If adata.raw exists and looks like log-normalized → return adata.raw (usable for some steps)
3. Else → raise with clear error message
```

### 2.4 `summarize_state(state) → str`

Human-readable Markdown summary for the agent to present:

```
## Data Inspection

**Shape:** 6,420 cells × 19,918 genes (human)
**X state:** z-score scaled (negative values detected, per-gene mean ≈ 0, std ≈ 1)
**Raw counts:** available in `layers['raw_counts']` ✅

### Completed Steps
✅ QC metrics  ✅ Normalized  ✅ HVG  ✅ PCA  ✅ Neighbors  ✅ Clustered  ✅ Cell types
❌ UMAP  ❌ DE results

### Metadata
- Clusters: `obs['leiden']` (12 clusters)
- Cell types: `obs['cell_type']` (8 types)  
- Condition: `obs['condition']` (disease, healthy)
- Batch: `obs['sample']` (4 samples)

### ⚠️ Warnings
- X is z-score scaled. Will use `layers['raw_counts']` as starting point for any re-processing.
```

---

## 3. Dependency Graph (replaces rigid DAG)

### 3.1 New: `scagent/dependencies.py`

Static knowledge about what each analysis step requires. Not a plan — just facts.

```python
DEPENDENCIES = {
    'qc_metrics':       {'requires_x_state': ['raw_counts']},
    'filter_cells':     {'requires': ['qc_metrics']},
    'filter_genes':     {'requires': ['filter_cells']},
    'doublet_detection':{'requires': ['filter_genes']},
    'normalize':        {'requires_x_state': ['raw_counts']},
    'hvg':              {'requires_x_state': ['log_normalized']},
    'pca':              {'requires': ['hvg']},
    'neighbors':        {'requires': ['pca']},
    'clustering':       {'requires': ['neighbors']},
    'umap':             {'requires': ['neighbors']},
    'marker_genes':     {'requires': ['clustering']},
    'annotation':       {'requires': ['clustering']},
    'pseudobulk_de':    {'requires': ['annotation'], 'requires_data': ['raw_counts', 'condition_key']},
    'trajectory':       {'requires': ['neighbors']},
    'gsea':             {'requires': ['pseudobulk_de']},
    'composition':      {'requires': ['annotation'], 'requires_data': ['condition_key']},
}
```

### 3.2 `check_prerequisites(goal, state) → (can_run, missing)`

Given a goal (what the user wants) and the inspector state, returns what's missing:

```python
check_prerequisites('pseudobulk_de', state)
# → (False, ['annotation not found', 'condition_key not found'])

check_prerequisites('clustering', state)  
# → (True, [])  # neighbors already exist

check_prerequisites('marker_genes', state)
# → (True, [])  # clusters already exist
```

### 3.3 `plan_steps(goal, state) → list[str]`

Work backwards from the goal, skip steps already completed, return only what's needed:

```python
plan_steps('pseudobulk_de', state)
# If state has clusters but no annotation:
# → ['annotation', 'pseudobulk_de']

plan_steps('clustering', state)
# If state has PCA but no neighbors:
# → ['neighbors', 'clustering']

plan_steps('umap', state)
# If state already has UMAP:
# → []  (nothing to do)
```

---

## 4. Refactor: `scagent/dag.py`

The existing `AnalysisDAG` is not deleted. It changes role:

### Current behavior (remove)
- Generated at init, required before analysis
- Gates all steps — "you can't do this, it's not in the DAG"
- Fixed per paradigm

### New behavior
- **Optional.** Generated when the user states a clear goal ("I want to build a cell atlas")
- **Advisory.** Shows a suggested plan, tracks progress, but doesn't block other actions
- **Appendable.** User can ask for trajectory mid-atlas and it gets added
- Still useful for: showing progress, suggesting next steps, reproducibility export

### Changes to `AnalysisDAG`

```python
class AnalysisDAG:
    # Keep: from_context(), next_step(), complete_step(), skip_step(), summary(), save/load
    # Add:
    def add_step(self, step: DAGStep, after: str | None = None) -> None:
        """Add a step to the DAG dynamically."""
    
    def mark_precomputed(self, step_id: str) -> None:
        """Mark a step as already done (from inspector)."""
    
    # Change: is_valid_step() consults dependencies.py, not just the DAG
```

### Changes to DAG skill (`.pi/skills/dag/SKILL.md`)

- Remove: "Before each tool call, check what the DAG says to do next"
- Add: "If no DAG exists, use dependencies to validate the step is possible"
- Add: "If the user asks for something not in the DAG, check dependencies, execute if valid, optionally add to DAG"

---

## 5. Refactor: Init Flow

### Changes to `ExperimentContext`

- **Paradigm becomes optional.** Context stores it if provided, but doesn't require it.
- Add: `context.data_state: dict | None` — stores the inspector output at init time
- `validate()` no longer fails on missing paradigm — only warns

### Changes to init skill (`.pi/skills/init/SKILL.md`)

Two modes:

**Mode A: Fresh start (no pre-processed data)**
Same as today: load raw data, infer species/platform, ask paradigm/tissue/conditions, generate DAG.

**Mode B: Onboarding (pre-processed data)**  
New flow:
1. Run `inspect_adata()` on the provided data
2. Show the researcher what was detected (summarize_state)
3. Ask: "What would you like to do with this data?" (no paradigm menu — free-form question)
4. Use `plan_steps()` to figure out what's needed
5. Optionally generate a DAG if the researcher wants a structured plan
6. Record the inspection result in `project.json` for future reference

---

## 6. Eval Handler Fixes

Update all eval handlers to use the inspector instead of ad-hoc checks:

### `eval/handlers/hvg.py`
```python
# Before (broken):
if adata.X.max() > 50:
    log_normalize(adata)

# After:
state = inspect_adata(adata)
if state.x_state != 'log_normalized':
    adata.X = find_raw_counts(adata, state)
    log_normalize(adata)
```

### `eval/handlers/celltyping.py`
Same pattern — use inspector to determine starting point.

### `eval/handlers/trajectory.py`
Same pattern — check if data needs preprocessing before trajectory inference.

### General pattern for all handlers
```python
def handle_X(adata, task_prompt):
    state = inspect_adata(adata)
    adata = ensure_ready_for(adata, state, needs='log_normalized')
    # ... actual analysis ...
```

`ensure_ready_for()` is a convenience function that:
- If X is in the right state → do nothing
- If raw counts are available → re-process from those
- If nothing is available → raise with clear error

---

## 7. Architecture Doc Updates

### Section 1: Design Philosophy
Update the three pillars:
```
1. Paradigm-aware routing → Question-driven analysis with dependency awareness
2. Provenance-first → (unchanged)
3. Branched exploration → (unchanged)
```

Add a new principle:
> **State-aware onboarding** — the system inspects foreign data to determine what has been done, 
> what's available, and what's needed. It never assumes it controls the data from the start.

### Section 3: Experiment Context
- Note that paradigm is now optional
- Add subsection on data inspection at onboarding

### Section 5: Analysis DAG
Major rewrite:
- Rename to "Analysis Planning"
- Explain the two-layer system: static dependency graph (always active) + optional DAG (when user wants a plan)
- Show how the agent handles "user asks for something not in the DAG"
- Show the onboarding flow for foreign data

### New Section (after Section 5): Data Inspection
- The `inspect_adata()` function and what it checks
- The `AnnDataState` dataclass
- How the agent uses it: "inspect once, reason per question"
- Examples of each X state and how to handle them

### Section 9: REPL Interface
Add example of onboarding flow:
```
You: [drops a pre-processed .h5ad]

scAgent: I've inspected your data:
  6,420 cells × 19,918 genes (mouse)
  X is z-score scaled, raw counts in layers['raw_counts']
  Clusters: obs['leiden'] (12 clusters)
  Cell types: obs['cell_type'] (8 types)
  
  What would you like to do?

You: Find DE genes between the CAF subtypes

scAgent: I see 'CAF' in your cell type labels. I'll:
  1. Subset to CAF cells
  2. Re-cluster to find subtypes (using raw counts → normalize → PCA → cluster)
  3. Run DE between subtypes
  
  Proceed?
```

### Section 11: How This Beats scBench
Update the table:

| scBench Failure Mode | scAgent Design Response |
|---------------------|--------------------------|
| Agent doesn't know data is pre-processed | **Inspector** detects X state, finds raw counts in layers |
| Agent applies wrong normalization | **Inspector** checks X state before any transformation |
| (keep existing rows) | |

---

## 8. Testing

### 8.1 Unit Tests: `tests/test_inspector.py`

| Test | What it checks |
|------|---------------|
| `test_detect_raw_counts` | Integer non-negative X → 'raw_counts' |
| `test_detect_log_normalized` | Float non-negative X + uns['log1p'] → 'log_normalized' |
| `test_detect_scaled` | Negative values → 'scaled' |
| `test_detect_transformed_unknown` | Float non-negative, no uns['log1p'], ambiguous range → 'transformed_unknown' |
| `test_find_raw_in_layers_counts` | layers['counts'] found → correct location |
| `test_find_raw_in_layers_raw_counts` | layers['raw_counts'] found → correct location |
| `test_find_raw_in_x` | X is raw counts, no layers → 'X' |
| `test_find_raw_none` | Scaled X, no layers → None + warning |
| `test_breadcrumbs_qc` | obs has QC columns → has_qc_metrics=True |
| `test_breadcrumbs_pca` | obsm has X_pca → has_pca=True |
| `test_breadcrumbs_neighbors` | uns has neighbors + obsp → has_neighbors=True |
| `test_breadcrumbs_clusters` | obs has leiden column → has_clusters=True, cluster_key='leiden' |
| `test_breadcrumbs_empty` | Fresh AnnData → all False |
| `test_fuzzy_cluster_key` | obs has 'leiden_0.8' → finds it |
| `test_fuzzy_celltype_key` | obs has 'CellTypist_annotation' → finds it |
| `test_fuzzy_condition_key` | obs has 'disease_status' → finds it |
| `test_fuzzy_batch_key` | obs has 'donor_id' → finds it |
| `test_species_human` | Uppercase genes → 'human' |
| `test_species_mouse` | Lowercase-initial genes → 'mouse' |
| `test_summarize_state` | Returns readable Markdown with correct status icons |

### 8.2 Unit Tests: `tests/test_dependencies.py`

| Test | What it checks |
|------|---------------|
| `test_check_prerequisites_met` | Clustering with neighbors present → (True, []) |
| `test_check_prerequisites_missing` | Clustering without neighbors → (False, ['neighbors']) |
| `test_plan_steps_nothing_needed` | Goal already completed → [] |
| `test_plan_steps_fill_gaps` | DE needs annotation needs clustering → correct chain |
| `test_plan_steps_skips_completed` | PCA done, need clustering → [neighbors, clustering] only |
| `test_pseudobulk_needs_raw_counts` | Pseudobulk without raw counts → error |
| `test_pseudobulk_needs_condition` | Pseudobulk without condition column → error |

### 8.3 Integration Tests: `tests/test_inspector_integration.py`

Use the actual scBench eval data:

| Test | What it checks |
|------|---------------|
| `test_inspect_raw_data` | QC eval data (raw counts) → correct detection |
| `test_inspect_scaled_data` | HVG eval data (z-scored) → detects scaled + finds layers['raw_counts'] |
| `test_inspect_annotated_data` | Clustering eval data → detects cell types + clusters |
| `test_ensure_ready_for_hvg` | Scaled data → recovers raw counts → normalizes → HVG works |

### 8.4 Updated Tests: `tests/test_dag.py`

| Test | What it checks |
|------|---------------|
| `test_dag_optional` | System works without a DAG (dependency graph only) |
| `test_add_step_dynamic` | Can add trajectory step to atlas DAG mid-analysis |
| `test_mark_precomputed` | Inspector state → DAG steps marked as already done |

---

## 9. File Changes Summary

| File | Action | Description |
|------|--------|-------------|
| `scagent/inspector.py` | **New** | `inspect_adata()`, `find_raw_counts()`, `summarize_state()`, `AnnDataState` |
| `scagent/dependencies.py` | **New** | `DEPENDENCIES`, `check_prerequisites()`, `plan_steps()`, `ensure_ready_for()` |
| `scagent/dag.py` | **Modify** | Add `add_step()`, `mark_precomputed()`. Relax `is_valid_step()` to use dependencies |
| `scagent/context.py` | **Modify** | Make paradigm optional in `validate()`. Add `data_state` field |
| `scagent/__init__.py` | **Modify** | Export `inspect_adata`, `AnnDataState`, `check_prerequisites`, `plan_steps` |
| `outputs/architecture.md` | **Modify** | Update pillars, add inspector section, rewrite DAG section, update REPL examples |
| `.pi/skills/init/SKILL.md` | **Modify** | Add Mode B (onboarding foreign data) |
| `.pi/skills/dag/SKILL.md` | **Modify** | DAG is advisory, dependency graph validates, handle out-of-DAG requests |
| `eval/handlers/hvg.py` | **Modify** | Use inspector instead of `X.max() > 50` |
| `eval/handlers/celltyping.py` | **Modify** | Use inspector |
| `eval/handlers/trajectory.py` | **Modify** | Use inspector |
| `tests/test_inspector.py` | **New** | 19 unit tests |
| `tests/test_dependencies.py` | **New** | 7 unit tests |
| `tests/test_inspector_integration.py` | **New** | 4 integration tests |
| `tests/test_dag.py` | **Modify** | Add 3 tests for new DAG behavior |

---

## 10. Implementation Order

```
Phase 1: Inspector (core, no dependencies on other changes)
  1. scagent/inspector.py — AnnDataState, inspect_adata(), find_raw_counts(), summarize_state()
  2. tests/test_inspector.py — unit tests
  3. tests/test_inspector_integration.py — integration tests with eval data

Phase 2: Dependencies (builds on inspector)
  4. scagent/dependencies.py — DEPENDENCIES, check_prerequisites(), plan_steps(), ensure_ready_for()
  5. tests/test_dependencies.py — unit tests

Phase 3: DAG refactor (builds on dependencies)
  6. scagent/dag.py — add_step(), mark_precomputed(), relax is_valid_step()
  7. scagent/context.py — make paradigm optional
  8. tests/test_dag.py — update + add new tests
  9. .pi/skills/init/SKILL.md — add Mode B
  10. .pi/skills/dag/SKILL.md — advisory behavior

Phase 4: Eval fixes (uses inspector)
  11. eval/handlers/*.py — use inspector in all handlers
  12. Run scBench — verify HVG eval no longer crashes

Phase 5: Documentation
  13. outputs/architecture.md — full update
  14. scagent/__init__.py — exports
```

---

## 11. What This Does NOT Change

| Component | Status | Why |
|-----------|--------|-----|
| Provenance | Unchanged | Still records every step. Now also records inspector output at onboarding. |
| State manager (branching) | Unchanged | Branches still work the same way |
| Memory (MemPalace) | Unchanged | Still stores conversations and decisions |
| Knowledge (MarkerDB) | Unchanged | Still provides marker lookups |
| Export (repro package) | Unchanged | Still exports from provenance |
| Tool wrappers | Unchanged | Individual tools don't change — the inspector sits above them |

---

## 12. Open Questions

| # | Question | Options | Recommendation |
|---|----------|---------|----------------|
| 1 | Should `inspect_adata()` be expensive (full matrix scan) or cheap (sample-based)? | Full scan, sample first 200 rows, configurable | **Sample first 200 rows.** The negative check and integer check only need a sample. Per-gene mean/std for scaled detection can use the sample too. Full scan is unnecessary. |
| 2 | Should the inspector store its result in `project.json`? | Yes (persistent), No (re-run each session) | **Yes.** Store as `data_state` in project.json so the agent remembers across sessions. Re-inspect only if the data file changes. |
| 3 | How to handle edge case: data is log-normalized but `uns['log1p']` is missing (R pipeline)? | Trust heuristic (max < 20), warn, ask user | **Heuristic + warn.** Label as 'log_normalized' if max < 20 and values are non-negative non-integer, but add a warning: "No Scanpy log-transform marker found — inferring from value range." |
| 4 | Should `ensure_ready_for()` modify adata in place or return a copy? | In place, copy, configurable | **In place** (consistent with Scanpy convention). Document that it modifies adata. |
| 5 | Should the existing `ExperimentContext.infer_from_data()` be merged into the inspector? | Merge, keep separate, inspector calls it | **Inspector calls it.** `infer_from_data()` handles species/platform inference. Inspector handles everything else. No duplication. |
