# Chunk 2 Design Plan: Skills — QC & Preprocessing (Revised)

**Status:** Draft for review (v2 — revised based on feedback)
**Depends on:** Chunk 1 (tool registry — done)
**Feeds into:** Chunk 3 (clustering & annotation skills)

---

## Design Change from v1

**v1 approach (rejected):** Skills contain Scanpy code examples. Agent reads skill, writes Python scripts from scratch each time. Problems: hallucination risk, untestable, unreproducible.

**v2 approach:** We write **deterministic Python wrapper functions** in `scagent/tools/`. The agent just calls them. Skills teach the agent WHEN to call each function and WHAT parameters to choose — not how to implement the underlying Scanpy logic.

This means Chunk 2 now builds two things:
1. **Python module** (`scagent/tools/`) — tested, deterministic functions
2. **Skills** (`.pi/skills/`) — instruction documents for the agent

---

## 1. Python Module Design

### 1.1 Module Structure

```
scagent/
  __init__.py
  tools/
    __init__.py
    loading.py       # load_10x_h5()
    qc.py            # calculate_qc_metrics(), filter_cells(), filter_genes()
    doublets.py      # detect_doublets()
    normalize.py     # log_normalize()
```

### 1.2 API Conventions

Every function follows these rules:

1. **Explicit parameters** — no `**kwargs`, no hidden state. Every meaningful parameter is a named argument with a default.

2. **Per-step random seeds** — any function using randomness accepts `random_state: int = 0`. This seed is passed to the underlying library call, NOT set globally via `np.random.seed()`. Re-running with the same seed always gives the same result, regardless of prior execution history.

3. **Return a result dict** — every function returns a dict with:
   - `metrics`: key numbers (cells_before, cells_after, etc.)
   - `plots`: list of file paths to generated plots (empty if no plot_dir)
   - `provenance`: dict of parameters used + results (for future PROV-JSONLD recording)
   - `warnings`: list of warning strings (e.g., ">50% cells removed")

4. **AnnData handling**:
   - Functions that subset (filter_cells, filter_genes, detect_doublets) return a NEW adata + result dict
   - Functions that modify in place (calculate_qc_metrics, log_normalize) return only result dict
   - This makes it explicit when adata identity changes (important for checkpointing)

5. **Checkpoint saving** — each function accepts an optional `checkpoint_dir: str = None`. If provided, saves adata after the step to `{checkpoint_dir}/{step_name}.h5ad`. This enables going back to any intermediate state.

6. **No side effects** beyond the adata and files — no global state, no singletons, no class instances to manage.

### 1.3 Function Signatures

```python
# === loading.py ===

def load_10x_h5(
    filename: str,
    *,
    gex_only: bool = True,
    genome: str | None = None,
    checkpoint_dir: str | None = None,
) -> tuple[AnnData, dict]:
    """
    Load 10x Genomics filtered matrix.
    
    Returns:
        (adata, result) where result contains metrics, plots, provenance, warnings
    """

# === qc.py ===

def calculate_qc_metrics(
    adata: AnnData,
    *,
    species: str | None = None,   # 'human' or 'mouse'; auto-detected if None
    plot_dir: str | None = None,
) -> dict:
    """
    Compute QC metrics and generate diagnostic plots. Modifies adata in place.
    
    Adds to adata.obs: n_genes_by_counts, total_counts, pct_counts_mt
    Adds to adata.var: mt (boolean), n_cells_by_counts
    Stores detected species in adata.uns['species']
    
    Returns result dict with metrics (medians, MAD-based threshold recommendations)
    and plot paths.
    """

def filter_cells(
    adata: AnnData,
    *,
    min_genes: int = 200,
    max_genes: int = 5000,
    max_pct_mito: float = 10.0,
    min_counts: int | None = None,
    checkpoint_dir: str | None = None,
) -> tuple[AnnData, dict]:
    """
    Filter low-quality cells.
    
    Returns:
        (filtered_adata, result) — new AnnData object (subset)
    """

def filter_genes(
    adata: AnnData,
    *,
    min_cells: int = 3,
    checkpoint_dir: str | None = None,
) -> tuple[AnnData, dict]:
    """
    Filter rarely-expressed genes.
    
    Returns:
        (filtered_adata, result) — new AnnData object (subset)
    """

# === doublets.py ===

def detect_doublets(
    adata: AnnData,
    *,
    expected_doublet_rate: float | None = None,  # auto-calculate if None
    n_prin_comps: int = 30,
    sim_doublet_ratio: float = 2.0,
    random_state: int = 0,
    remove: bool = True,          # if False, annotate only (don't filter)
    plot_dir: str | None = None,
    checkpoint_dir: str | None = None,
) -> tuple[AnnData, dict]:
    """
    Detect doublets using Scrublet. Operates on raw counts.
    
    If expected_doublet_rate is None, auto-calculates as n_cells / 1000 * 0.008
    (10x Chromium rate: ~0.8% per 1,000 cells captured).
    
    Returns:
        (adata_with_doublets_handled, result)
    """

# === normalize.py ===

def log_normalize(
    adata: AnnData,
    *,
    target_sum: float = 1e4,
    exclude_highly_expressed: bool = False,
    checkpoint_dir: str | None = None,
) -> dict:
    """
    Log-normalize counts (CP10K + log1p). Modifies adata in place.
    
    Order of operations:
    1. Copy raw counts to adata.layers['counts']
    2. sc.pp.normalize_total(adata, target_sum=target_sum)
    3. sc.pp.log1p(adata)
    4. adata.raw = adata.copy()  (freezes full gene set, normalized)
    
    This order ensures:
    - adata.layers['counts'] = raw integer UMIs (for pseudobulk DE)
    - adata.raw.X = log-normalized, all genes (for Wilcoxon markers, plotting)
    - adata.X = log-normalized (will be subset to HVGs later for PCA)
    
    Returns result dict.
    """
```

### 1.4 Random Seed Strategy

**Principle:** Each function owns its own seed. No global RNG state dependency.

| Function | Uses randomness? | How seed is applied |
|----------|-----------------|---------------------|
| `load_10x_h5` | No | N/A |
| `calculate_qc_metrics` | No | N/A |
| `filter_cells` | No | N/A |
| `filter_genes` | No | N/A |
| `detect_doublets` | **Yes** (Scrublet simulates doublets) | `np.random.RandomState(random_state)` passed to Scrublet constructor |
| `log_normalize` | No | N/A |

For Chunk 3's functions (PCA, neighbors, UMAP, Leiden), randomness becomes more important. Same pattern applies: each function takes `random_state`, passes it to the Scanpy call, records it in provenance.

**Re-running a step:** If the user goes back to a checkpoint and re-runs a step with the same `random_state`, they get identical results. If they want to explore sensitivity, they change the seed — and the result dict records which seed was used.

**NOT doing:**
- Global `np.random.seed()` — breaks reproducibility on partial re-runs
- Derived seeds (`hash(master + step_name)`) — clever but unnecessarily complex. Just use explicit per-step seeds. Default 0 is fine.

### 1.5 Checkpoint Strategy

When `checkpoint_dir` is provided:
```python
# After filter_cells completes:
adata.write(f"{checkpoint_dir}/after_filter_cells.h5ad")
```

Naming convention: `after_{step_name}.h5ad`

Checkpoints created by the pipeline:
```
after_load.h5ad           # raw loaded data
after_qc_metrics.h5ad     # with QC annotations (no filtering yet)
after_filter_cells.h5ad   # after cell filtering
after_filter_genes.h5ad   # after gene filtering
after_doublets.h5ad       # after doublet removal
after_normalize.h5ad      # after normalization
```

The agent can tell the user "Going back to the state after QC filtering" and load the appropriate checkpoint. This is a lightweight precursor to the full state manager (Chunk 5) but immediately usable.

---

## 2. Skill Design

Skills become shorter and cleaner — they teach the agent WHEN and WHY, not HOW.

### 2.1 `load-data`

**Frontmatter:**
```yaml
name: load-data
description: Load 10x Genomics Chromium data from .h5 files. Use when the user provides a data file, asks to load data, or starts a new analysis.
```

**Body teaches the agent:**
- Call `scagent.tools.loading.load_10x_h5(filename)` 
- Report the returned metrics (n_cells, n_genes, species)
- Warn if duplicate var_names were fixed
- Suggest running QC next
- Save checkpoint

**~40 lines**

### 2.2 `qc`

**Frontmatter:**
```yaml
name: qc
description: Run quality control on single-cell data. Compute QC metrics, show distributions, recommend and apply cell/gene filtering thresholds. Use when the user asks for QC, quality check, filtering, or after loading data.
```

**Body teaches the agent:**
- Step 1: Call `calculate_qc_metrics(adata)` — show the returned plots to the user, present the MAD-based threshold recommendations from the result dict
- Step 2: Discuss thresholds with the user. Key guidance:
  - Mito: 5% for PBMCs, 10-15% for solid tissue, up to 20% for heart/kidney
  - Look at the distribution plots — are there clear outlier populations?
  - If bimodal mito distribution, the high-mito cells might be real biology
- Step 3: Call `filter_cells(adata, ...)` with confirmed thresholds
- Step 4: Call `filter_genes(adata, min_cells=3)`
- Step 5: Show before/after summary from result dicts
- Step 6: Suggest doublet detection next

**Key instruction to agent:** Always show the distributions before recommending thresholds. Never filter without asking.

**~80 lines** (shorter than v1 because code is in the module)

### 2.3 `doublet-detection`

**Frontmatter:**
```yaml
name: doublet-detection
description: Detect cell doublets using Scrublet. Use after QC filtering and before normalization. Must run on raw (unnormalized) counts.
```

**Body teaches the agent:**
- Call `detect_doublets(adata)` — auto-calculates expected rate from cell count
- Show the doublet score histogram to the user
- Compare detected rate to expected rate — flag if >2× different
- If remove=True (default), report how many cells were removed
- If user wants to explore, suggest remove=False to just annotate
- Suggest normalization next

**~50 lines**

### 2.4 `normalize`

**Frontmatter:**
```yaml
name: normalize
description: Normalize and log-transform single-cell count data. Use after QC and doublet detection. Preserves raw counts for downstream differential expression.
```

**Body teaches the agent:**
- Call `log_normalize(adata)` — handles the full order of operations internally
- Explain to user what happened: raw counts saved, CP10K normalization, log-transform
- Show validation results from the result dict
- Mention that alternatives exist (SCTransform, scran) but CP10K + log1p is the standard default
- Suggest feature selection (HVG) as the next step

**~40 lines**

---

## 3. Verification Plan

### 3.1 Unit-level tests for Python functions

`tests/test_qc_pipeline.py` — exercises each function on the PBMC 10k data:

```python
def test_load():
    adata, info = load_10x_h5("data/pbmc10k/filtered_feature_bc_matrix.h5")
    assert adata.n_obs == 11769
    assert adata.n_vars == 33538
    assert info['metrics']['species'] == 'human'

def test_qc_metrics():
    result = calculate_qc_metrics(adata)
    assert 'pct_counts_mt' in adata.obs.columns
    assert result['metrics']['median_pct_mito'] < 20  # PBMCs should be low

def test_filter_cells():
    adata_f, info = filter_cells(adata, min_genes=200, max_genes=5000, max_pct_mito=5.0)
    assert adata_f.n_obs < adata.n_obs
    assert adata_f.n_obs > adata.n_obs * 0.5  # didn't lose >50%
    assert info['metrics']['cells_after'] == adata_f.n_obs

def test_doublets():
    adata_d, info = detect_doublets(adata, random_state=0)
    rate = info['metrics']['doublet_rate']
    assert 0.01 < rate < 0.15

def test_normalize():
    info = log_normalize(adata)
    assert adata.X.max() < 15
    assert 'counts' in adata.layers
    assert adata.raw is not None

def test_reproducibility():
    """Same seed = same result"""
    adata1, info1 = detect_doublets(adata_copy1, random_state=42)
    adata2, info2 = detect_doublets(adata_copy2, random_state=42)
    assert info1['metrics']['n_doublets'] == info2['metrics']['n_doublets']
    
    adata3, info3 = detect_doublets(adata_copy3, random_state=99)
    # Different seed may give different count (but both should be reasonable)
```

### 3.2 Expected results for PBMC 10k

| Step | Expected |
|------|----------|
| Load | 11,769 cells × 33,538 genes |
| QC metrics | Median mito ~3-5%, median genes ~1500-2500 |
| Filter (5% mito) | ~10,500-11,500 cells retained |
| Filter genes | ~18,000-20,000 genes retained (min_cells=3) |
| Doublets | ~3-6% detected (~400-700 cells) |
| Normalize | Max value ~10-12 |

---

## 4. File List

| File | Purpose |
|------|---------|
| `scagent/__init__.py` | Package init |
| `scagent/tools/__init__.py` | Tools subpackage init |
| `scagent/tools/loading.py` | `load_10x_h5()` |
| `scagent/tools/qc.py` | `calculate_qc_metrics()`, `filter_cells()`, `filter_genes()` |
| `scagent/tools/doublets.py` | `detect_doublets()` |
| `scagent/tools/normalize.py` | `log_normalize()` |
| `.pi/skills/load-data/SKILL.md` | Agent instructions for data loading |
| `.pi/skills/qc/SKILL.md` | Agent instructions for QC |
| `.pi/skills/doublet-detection/SKILL.md` | Agent instructions for doublet detection |
| `.pi/skills/normalize/SKILL.md` | Agent instructions for normalization |
| `tests/test_qc_pipeline.py` | Smoke test on PBMC 10k |

---

## 5. Open Questions (Resolved)

| Question | Decision |
|----------|----------|
| Auto-chain vs. wait between steps? | Auto-chain by default. Checkpoints at each step allow going back to any intermediate state. |
| Agent writes Scanpy code or calls our functions? | Calls our functions. Python module handles all boilerplate. |
| Global seed or per-step seed? | Per-step. Each function takes `random_state`, passes to underlying call. Default 0. No global `np.random.seed()`. |
| Where to set `adata.raw`? | In `log_normalize()`, after log-transform. Documented so Chunk 3 respects it. |
| How to go back to intermediate state? | Checkpoints: `after_{step}.h5ad` saved when `checkpoint_dir` is provided. Agent loads the appropriate checkpoint. |

---

## 6. Build Order

1. `scagent/__init__.py` + `scagent/tools/__init__.py` (trivial)
2. `scagent/tools/loading.py` → test on PBMC data
3. `scagent/tools/qc.py` → test on PBMC data
4. `scagent/tools/doublets.py` → test on PBMC data
5. `scagent/tools/normalize.py` → test on PBMC data
6. `tests/test_qc_pipeline.py` — full pipeline smoke test
7. `.pi/skills/` — all 4 SKILL.md files (these are fast once the module works)
8. Git commit

Each step is testable independently. We don't move to step N+1 until step N passes its tests.
