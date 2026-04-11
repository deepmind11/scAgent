# Chunk 3 Design Plan: Skills — Clustering & Annotation

**Status:** Draft for review
**Depends on:** Chunk 2 (QC & preprocessing — done)
**Feeds into:** Chunk 7 (first scBench benchmark)

---

## 1. What This Chunk Builds

The second half of the core analysis pipeline: from normalized data to annotated cell types.

**DAG for this chunk:**
```
[normalized adata from Chunk 2]
    │
    ▼
  HVG selection
    │
    ▼
   PCA
    │
    ├──────────────────────┐
    ▼                      ▼
  [Integration]       (skip if single-sample)
    │
    ▼
  Neighbor graph
    │
    ├──────────────┐
    ▼              ▼
  Leiden         UMAP
    │              │
    ▼              │
  Marker genes ◄───┘
    │
    ▼
  Cell annotation
```

**Milestone:** Run the full pipeline on PBMC 10k and verify we get the expected cell types:
CD4 T, CD8 T, B, NK, CD14 Monocytes, FCGR3A Monocytes, Dendritic cells, Platelets.

---

## 2. Python Module: 7 New Files

Same pattern as Chunk 2: deterministic wrapper functions, per-step random seeds, checkpoint support.

### 2.1 `scagent/tools/feature_selection.py`

```python
def select_hvg(
    adata,
    *,
    n_top_genes: int = 2000,
    flavor: str = "seurat",
    batch_key: str | None = None,
    min_mean: float = 0.0125,
    max_mean: float = 3.0,
    min_disp: float = 0.5,
    checkpoint_dir: str | None = None,
) -> dict:
    """Select highly variable genes. Modifies adata in place.
    
    Subsets adata.var to mark highly_variable genes.
    Does NOT subset adata.X yet — PCA handles that via use_highly_variable=True.
    """
```

- No randomness involved.
- `flavor="seurat"` works on log-normalized data (our default).
- `flavor="seurat_v3"` works on raw counts — would need different placement in the DAG.
- Returns: n_hvgs selected, plot of mean vs dispersion.

### 2.2 `scagent/tools/pca.py`

```python
def run_pca(
    adata,
    *,
    n_comps: int = 50,
    use_highly_variable: bool = True,
    svd_solver: str = "arpack",
    random_state: int = 0,
    plot_dir: str | None = None,
    checkpoint_dir: str | None = None,
) -> dict:
    """Run PCA. Modifies adata in place.
    
    Adds adata.obsm['X_pca'] and adata.uns['pca'].
    Generates elbow plot (variance ratio) to help choose n_comps.
    """
```

- `random_state` matters for `svd_solver="randomized"`. For `"arpack"` it's deterministic regardless, but we pass it anyway for consistency.
- Elbow plot is important — helps the user decide how many PCs to keep.
- Validation: PC1 should explain <50% variance (if it does, likely a technical artifact).

### 2.3 `scagent/tools/neighbors.py`

```python
def compute_neighbors(
    adata,
    *,
    n_neighbors: int = 15,
    n_pcs: int = 50,
    metric: str = "cosine",
    use_rep: str = "X_pca",
    method: str = "umap",
    random_state: int = 0,
    checkpoint_dir: str | None = None,
) -> dict:
    """Compute k-nearest neighbor graph. Modifies adata in place.
    
    Uses Sciaraffa et al. 2025 optimal defaults: cosine metric, UMAP method.
    After Harmony, set use_rep='X_pca_harmony'.
    After scVI, set use_rep='X_scVI'.
    """
```

- `random_state` affects the approximate nearest neighbor search.
- `use_rep` is critical — must point to the correct embedding (raw PCA vs. integrated).

### 2.4 `scagent/tools/embedding.py`

```python
def run_umap(
    adata,
    *,
    min_dist: float = 0.5,
    spread: float = 1.0,
    random_state: int = 0,
    checkpoint_dir: str | None = None,
) -> dict:
    """Compute UMAP embedding for visualization. Modifies adata in place.
    
    UMAP is for VISUALIZATION ONLY. Never cluster on UMAP coordinates.
    """
```

- `random_state` matters — UMAP is stochastic. Different seeds give different layouts.
- The layout can look very different with different seeds, but the topology should be stable.

### 2.5 `scagent/tools/clustering.py`

```python
def run_leiden(
    adata,
    *,
    resolution: float = 1.0,
    key_added: str = "leiden",
    random_state: int = 0,
    n_iterations: int = -1,
    checkpoint_dir: str | None = None,
) -> dict:
    """Run Leiden clustering. Modifies adata in place.
    
    Adds cluster labels to adata.obs[key_added].
    """

def sweep_resolution(
    adata,
    *,
    resolutions: list[float] = [0.3, 0.5, 0.8, 1.0, 1.5, 2.0],
    random_state: int = 0,
    plot_dir: str | None = None,
) -> dict:
    """Run Leiden at multiple resolutions and report cluster counts.
    
    Helps the user choose the right resolution. Each resolution gets its own
    key in adata.obs (e.g., leiden_0.3, leiden_0.5, ...).
    
    Returns: table of resolution → n_clusters, plus a plot.
    Does NOT pick a resolution — that's the user's decision.
    """
```

**Sensitivity testing:** `sweep_resolution()` is where the user's earlier concern about "does my clustering change with a different seed?" gets addressed. The sweep shows how cluster count varies with resolution. We could also add a seed sensitivity check:

```python
def check_seed_sensitivity(
    adata,
    *,
    resolution: float = 1.0,
    seeds: list[int] = [0, 1, 2, 3, 4],
) -> dict:
    """Run Leiden at the same resolution with different seeds.
    
    Reports ARI (Adjusted Rand Index) between pairs to measure stability.
    If ARI < 0.9 across seeds, the clustering is fragile at this resolution.
    """
```

This is a natural place for it — after the user picks a resolution, we can optionally check whether it's stable.

### 2.6 `scagent/tools/markers.py`

```python
def find_marker_genes(
    adata,
    *,
    groupby: str = "leiden",
    method: str = "wilcoxon",
    n_genes: int = 100,
    use_raw: bool = True,
    corr_method: str = "benjamini-hochberg",
    plot_dir: str | None = None,
    checkpoint_dir: str | None = None,
) -> dict:
    """Find marker genes per cluster using Wilcoxon rank-sum test.
    
    Returns top markers per cluster with scores, logfc, pvals.
    Generates dotplot and/or heatmap of top markers.
    """
```

- No randomness.
- `use_raw=True` uses `adata.raw` (set in normalize skill) — tests on all genes, not just HVGs.
- Returns structured results: dict mapping cluster → list of (gene, score, logfc, pval_adj).
- Plot: dotplot of top 3-5 markers per cluster.

### 2.7 `scagent/tools/annotation.py`

```python
def annotate_celltypist(
    adata,
    *,
    model: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
    over_clustering: str | None = None,
    checkpoint_dir: str | None = None,
) -> dict:
    """Annotate cell types using CellTypist.
    
    Downloads model on first use.
    Adds predictions to adata.obs['celltypist_prediction'] and
    adata.obs['celltypist_majority_voting'].
    
    Returns: cell type counts, confidence summary.
    """
```

- No randomness (CellTypist is a logistic regression classifier).
- Model download happens automatically on first use — need to handle gracefully.
- `majority_voting=True` assigns cell type per cluster, more robust than per-cell.

### 2.8 Integration: `scagent/tools/integration.py`

```python
def run_harmony(
    adata,
    *,
    key: str,
    basis: str = "X_pca",
    adjusted_basis: str = "X_pca_harmony",
    max_iter_harmony: int = 10,
    random_state: int = 0,
    checkpoint_dir: str | None = None,
) -> dict:
    """Run Harmony batch correction. Modifies adata in place.
    
    Adds adata.obsm['X_pca_harmony'].
    Only needed for multi-sample/multi-batch experiments.
    For single-sample data, skip and use X_pca directly.
    """
```

- **Conditional step**: The PBMC 10k is single-sample, so we won't run integration on it. But the code and skill need to exist for multi-sample experiments.
- `random_state` for Harmony's iterative optimization.
- The skill needs to teach the agent HOW TO DECIDE whether integration is needed (check if `batch_key` exists in obs).

---

## 3. Pipeline Runner Extension

Extend `pipeline.py` with a clustering pipeline:

```python
CLUSTERING_STEPS = [
    "hvg",
    "pca",
    "integration",  # conditional — skipped if no batch_key
    "neighbors",
    "leiden",
    "umap",
    "markers",
    "annotation",
]

def run_clustering(
    adata,
    *,
    # HVG params
    n_top_genes: int = 2000,
    # PCA params
    n_comps: int = 50,
    # Neighbor params
    n_neighbors: int = 15,
    metric: str = "cosine",
    # Integration params
    batch_key: str | None = None,    # if None, skip integration
    # Clustering params
    resolution: float = 1.0,
    # Control
    random_state: int = 0,
    checkpoint_dir: str | None = None,
    plot_dir: str | None = None,
    stop_after: str | None = None,
    start_from: str | None = None,
) -> tuple[AnnData, dict]:
```

---

## 4. Random Seeds in This Chunk

| Function | Uses randomness? | Seed handling |
|----------|-----------------|---------------|
| `select_hvg` | No | N/A |
| `run_pca` | Depends on solver | Pass to `sc.pp.pca(random_state=...)` |
| `compute_neighbors` | Yes (ANN search) | Pass to `sc.pp.neighbors(random_state=...)` |
| `run_umap` | Yes | Pass to `sc.tl.umap(random_state=...)` |
| `run_leiden` | Yes | Pass to `sc.tl.leiden(random_state=...)` |
| `run_harmony` | Yes | Pass to `harmonypy` |
| `find_marker_genes` | No | N/A |
| `annotate_celltypist` | No | N/A |

**All stochastic functions default to `random_state=0`.** Same data + same params + same seed = same result.

**Sensitivity check for clustering:** The `sweep_resolution()` and `check_seed_sensitivity()` functions let the user explicitly test whether their clustering is robust. The skill should suggest this after the initial Leiden run.

---

## 5. Skills: 6 New Files

### 5.1 `feature-selection` (~40 lines)
- Call `select_hvg()`, show mean-vs-dispersion plot
- Default 2000 genes. Suggest higher (3000-5000) for heterogeneous tissues.
- Mention `flavor="seurat_v3"` as alternative (but requires raw counts, different DAG position).

### 5.2 `dimensionality-reduction` (~60 lines)
- **PCA section:** Call `run_pca()`, show elbow plot, discuss how many PCs
- **UMAP section:** Call `run_umap()` AFTER neighbors step
- Emphasize: UMAP is visualization only, never cluster on it, never interpret distances

### 5.3 `integration` (~60 lines)
- **Decision logic:** Does `adata.obs` contain a batch/sample column? If single-sample, skip.
- Call `run_harmony()` if needed
- After integration, set `use_rep='X_pca_harmony'` for neighbors
- How to evaluate: UMAP colored by batch — batches should overlap

### 5.4 `clustering` (~80 lines)
- Call `sweep_resolution()` first — show resolution vs n_clusters table
- Help user pick resolution based on expected biology
- Call `run_leiden()` with chosen resolution
- Optionally run `check_seed_sensitivity()` — if ARI < 0.9, warn
- Sciaraffa et al. 2025 defaults: resolution=2.0 + k=10 + cosine — mention but don't force

### 5.5 `marker-genes` (~50 lines)
- Call `find_marker_genes()`, show dotplot of top markers
- Present as a table: cluster, top 5 genes, scores
- Use this to validate clusters (clusters without clear markers may need merging)
- This feeds directly into cell annotation

### 5.6 `cell-annotation` (~60 lines)
- Call `annotate_celltypist()` with appropriate model
- Model selection guidance: `Immune_All_Low.pkl` for PBMCs, tissue-specific models for others
- Cross-validate: compare CellTypist predictions against top marker genes
- Show: UMAP colored by cell type, cell type counts table
- For PBMCs, check that expected types are present: T cells, B cells, NK, monocytes, DCs

---

## 6. Verification Plan

### 6.1 Unit tests per function
Test each function individually on PBMC 10k data (already preprocessed in Chunk 2).

### 6.2 Full pipeline test
Run load → QC → normalize → HVG → PCA → neighbors → Leiden → UMAP → markers → annotation on PBMC 10k.

**Expected results:**
| Step | Expected |
|------|----------|
| HVG | ~2000 genes selected |
| PCA | 50 components, PC1 < 50% variance |
| Neighbors | k=15, cosine metric |
| Leiden (res=1.0) | ~10-15 clusters |
| UMAP | 2D embedding, no NaN |
| Markers | Each cluster has ≥5 significant genes |
| CellTypist | Should find: T cells (CD4/CD8), B cells, NK, Monocytes (CD14/FCGR3A), DCs, Platelets |

### 6.3 Reproducibility test
Run clustering twice with same seed → identical cluster labels.
Run with different seed → ARI > 0.9 (stable).

### 6.4 Seed sensitivity test
Run Leiden at resolution=1.0 with seeds 0-4. Compute pairwise ARI. All should be >0.9.

---

## 7. File List

| File | Purpose |
|------|---------|
| `scagent/tools/feature_selection.py` | `select_hvg()` |
| `scagent/tools/pca.py` | `run_pca()` |
| `scagent/tools/neighbors.py` | `compute_neighbors()` |
| `scagent/tools/embedding.py` | `run_umap()` |
| `scagent/tools/clustering.py` | `run_leiden()`, `sweep_resolution()`, `check_seed_sensitivity()` |
| `scagent/tools/markers.py` | `find_marker_genes()` |
| `scagent/tools/annotation.py` | `annotate_celltypist()` |
| `scagent/tools/integration.py` | `run_harmony()` |
| `scagent/tools/pipeline.py` | Extend with `run_clustering()` |
| `.pi/skills/feature-selection/SKILL.md` | HVG skill |
| `.pi/skills/dimensionality-reduction/SKILL.md` | PCA + UMAP skill |
| `.pi/skills/integration/SKILL.md` | Batch correction skill |
| `.pi/skills/clustering/SKILL.md` | Leiden + resolution sweep skill |
| `.pi/skills/marker-genes/SKILL.md` | Marker gene detection skill |
| `.pi/skills/cell-annotation/SKILL.md` | CellTypist annotation skill |
| `tests/test_clustering_pipeline.py` | Full pipeline test |

### Dependencies to add to pyproject.toml
- `celltypist` (already installed, need to add to deps)
- `scikit-learn` (for ARI in seed sensitivity check — likely already a transitive dep)

---

## 8. Open Questions

1. **Should `check_seed_sensitivity()` be automatic or opt-in?**
   Decision: Opt-in. The clustering skill mentions it and the agent can suggest it, but don't force it on every run. It adds ~5× the clustering time.

2. **UMAP seed sensitivity?**
   UMAP with different seeds gives different-looking plots but the topology should be stable. Not worth testing by default — it's just visualization.

3. **CellTypist model download — how to handle first-run?**
   CellTypist downloads models on first use (~50MB). The function should handle this gracefully — print a message, download, cache.

4. **Should we subset adata to HVGs or just mark them?**
   Decision: Mark them in `adata.var['highly_variable']`. PCA uses `use_highly_variable=True` to automatically restrict to HVGs without subsetting. This keeps all genes available in `adata.raw` for marker detection and plotting.

---

## 9. Build Order

1. `feature_selection.py` → test
2. `pca.py` → test  
3. `neighbors.py` → test
4. `embedding.py` → test
5. `clustering.py` (Leiden + sweep + seed sensitivity) → test
6. `markers.py` → test
7. `annotation.py` → test
8. `integration.py` → test (on synthetic multi-batch or just API test)
9. Extend `pipeline.py` with `run_clustering()`
10. Full pipeline test on PBMC 10k
11. Write 6 skills
12. Git commit
