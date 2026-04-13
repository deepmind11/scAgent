# scAgent: Architecture for an Agentic scRNA-seq Analysis System

**System Design Document — April 2026**

---

## 1. Design Philosophy

The goal is **not** to build an agent that passes scBench by memorizing answers. It's to build a system so well-structured that correctness follows from design. The hypothesis:

> An agent that deeply understands experimental paradigms, enforces proper tool sequencing, tracks provenance, and presents biologists with verifiable intermediate results will outperform general-purpose coding agents on scRNA-seq tasks — without ever seeing the benchmark.

Four pillars:
1. **Question-driven analysis with dependency awareness** — the system responds to what the researcher asks, validates steps against a dependency graph, and optionally generates a structured plan when the researcher has a clear goal. Paradigm selection is helpful but not required.
2. **State-aware onboarding** — the system inspects foreign data to determine what has been done, what's available, and what's needed. It never assumes it controls the data from the start.
3. **Provenance-first** — every result is traceable to inputs, parameters, and tool versions
4. **Branched exploration** — biologists can fork, compare, and merge analysis paths

---

## 2. System Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                       scAgent CLI (REPL)                        │
│  Researcher types natural language ←→ Agent responds + acts     │
└──────────────────────────┬──────────────────────────────────────┘
                           │
              ┌────────────┴────────────┐
              │     ORCHESTRATOR        │
              │  ┌──────────────────┐   │
              │  │ Paradigm Router  │   │  ← Identifies experiment type
              │  │ Plan Generator   │   │  ← Generates analysis DAG
              │  │ Step Executor    │   │  ← Runs tools with params
              │  │ Verifier         │   │  ← Checks outputs for sanity
              │  └──────────────────┘   │
              └────────────┬────────────┘
                           │
         ┌─────────────────┼─────────────────┐
         │                 │                 │
    ┌────┴────┐     ┌──────┴──────┐   ┌─────┴─────┐
    │  TOOL   │     │  PROVENANCE │   │  STATE    │
    │ REGISTRY│     │  GRAPH      │   │  MANAGER  │
    │         │     │ (PROV-JSONLD│   │           │
    │ Scanpy  │     │  lineage)   │   │ AnnData   │
    │ Seurat  │     │             │   │ snapshots │
    │ CellR.  │     │             │   │ branch    │
    │ DESeq2  │     │             │   │ tree      │
    │ ...     │     │             │   │           │
    └─────────┘     └─────────────┘   └───────────┘
         │                 │                 │
         └─────────────────┼─────────────────┘
                           │
         ┌─────────────────┼─────────────────┐
         │                 │                 │
    ┌────┴────┐     ┌──────┴──────┐   ┌─────┴─────┐
    │ MEMORY  │     │  KNOWLEDGE  │   │ EXPERIMENT│
    │ PALACE  │     │  MarkerDB + │   │ CONTEXT   │
    │(mempal.)│     │  best prac- │   │ (minSCe)  │
    │         │     │  tice refs  │   │           │
    │ Chat    │     │             │   │ Species   │
    │ history │     │ Canonical + │   │ Platform  │
    │ across  │     │ CellTypist +│   │ Tissue    │
    │ months  │     │ User DBs   │   │ Chemistry │
    └─────────┘     └─────────────┘   └───────────┘
```

---

## 3. Experiment Context (minSCe)

Based on Füllgrabe et al. (2020) "Guidelines for reporting single-cell RNA-seq experiments" (*Nature Biotechnology* 38:1384–1386), every scAgent project is initialized with a structured context object. This is the system's "ground truth" about the experiment and gates all downstream decisions.

```jsonc
// scagent-project.json
{
  "@context": "https://scagent.dev/context/v1",
  "@type": "SingleCellExperiment",

  // ── Identity ──
  "project_id": "proj_2026_04_pbmc_aging",
  "paradigm": "disease_vs_healthy",         // ← from the 9 paradigms
  "created": "2026-04-11T10:00:00Z",

  // ── minSCe metadata ──
  "organism": {
    "species": "Homo sapiens",
    "ncbi_taxon": 9606
  },
  "tissue": {
    "name": "peripheral blood mononuclear cells",
    "uberon_id": "UBERON:0000178"
  },
  "platform": {
    "vendor": "10x Genomics",
    "instrument": "Chromium",
    "chemistry": "GEM-X Single Cell 3' v4",
    "cell_ranger_version": "9.0.0"
  },
  "library": {
    "type": "3prime",
    "modalities": ["gene_expression"],      // could include: vdj, protein, crispr, atac
    "umi": true
  },
  "samples": [
    {
      "id": "donor_1", "condition": "young", "sex": "male", "age_range": "18-35",
      "target_cells": 5000, "gem_well": "GW1"
    },
    {
      "id": "donor_2", "condition": "old", "sex": "male", "age_range": "65-80",
      "target_cells": 5000, "gem_well": "GW2"
    }
  ],

  // ── Design ──
  "design": {
    "type": "case_control",
    "conditions": ["young", "old"],
    "biological_replicates_per_condition": 2,
    "batch_structure": "one_gem_well_per_sample"
  },

  // ── Hypotheses (optional, for the agent to reason about) ──
  "hypotheses": [
    "Age-associated changes in monocyte inflammatory gene programs",
    "Decline in naive T cell proportion with aging"
  ]
}
```

**Why this matters for the agent:** When the paradigm is `disease_vs_healthy`, the orchestrator knows that:
- Cross-condition DE (pseudobulk) is a required downstream step
- Batch correction is likely needed
- Differential abundance testing is relevant
- Trajectory inference is probably NOT the primary goal

When the paradigm is `developmental_trajectory`, the orchestrator knows:
- Pseudotime/trajectory inference IS the primary goal
- RNA velocity may be relevant
- Cross-condition DE may not apply

---

## 4. Tool Registry

Each tool is registered with its capabilities, valid paradigms, step placement, parameter schemas, and validation rules.

```jsonc
// tools/leiden_clustering.json
{
  "tool_id": "leiden_clustering",
  "name": "Leiden Clustering",
  "category": "clustering",
  "framework": "scanpy",
  "function": "sc.tl.leiden",

  "valid_after": ["neighbor_graph"],          // must come after this step
  "valid_before": ["differential_expression", "cell_annotation"],
  "paradigms": ["all"],                       // valid in all paradigms

  "parameters": {
    "resolution": {
      "type": "float",
      "default": 1.0,
      "range": [0.1, 5.0],
      "guidance": "Higher values → more clusters. Use 0.3-0.5 for major lineages, 0.8-1.2 for standard, 1.5-3.0 for subtypes. Sciaraffa et al. 2025 found resolution=2.0 optimal across 3 datasets.",
      "auto_strategy": "sweep",               // agent should try multiple values
      "sweep_values": [0.3, 0.5, 0.8, 1.0, 1.5, 2.0]
    },
    "key_added": {
      "type": "string",
      "default": "leiden_{resolution}"
    }
  },

  "outputs": {
    "cluster_labels": "adata.obs['leiden_{resolution}']",
    "n_clusters": "int"
  },

  "validation": {
    "min_clusters": 2,
    "max_clusters_ratio": 0.1,                // no more than 10% of cells as clusters
    "check": "all clusters have ≥ 10 cells"
  },

  "provenance_captures": ["resolution", "n_neighbors_used", "n_pcs_used", "random_state"]
}
```

### Tool Registry Structure

All 31 tool schemas live flat in `tools/` — no subdirectories. Categories are encoded in each schema's `category` field, not the filesystem. This keeps the registry simple to enumerate and extend.

```
tools/
├── load_10x_h5.json          # loading
├── filter_cells.json          # qc
├── filter_genes.json          # qc
├── scrublet_doublets.json     # qc
├── log_normalize.json         # normalization
├── highly_variable_genes.json # feature_selection
├── pca.json                   # dimensionality_reduction
├── umap.json                  # dimensionality_reduction
├── neighbors.json             # graph
├── harmony.json               # integration
├── scvi.json                  # integration
├── scanorama.json             # integration
├── bbknn.json                 # integration
├── leiden.json                # clustering
├── louvain.json               # clustering
├── celltypist.json            # annotation
├── validate_annotation.json   # annotation
├── query_markers.json         # annotation
├── wilcoxon_markers.json      # differential_expression
├── deseq2_pseudobulk.json     # differential_expression
├── edger_pseudobulk.json      # differential_expression
├── monocle3.json              # trajectory
├── paga.json                  # trajectory
├── slingshot.json             # trajectory
├── scvelo.json                # trajectory
├── cellchat.json              # communication
├── cellphonedb.json           # communication
├── gsea.json                  # enrichment
├── cluster_profiler.json      # enrichment
├── sccoda.json                # composition
└── milor.json                 # composition
```

---

## 5. Analysis Planning (DAG + Dependency Graph)

Analysis planning has two layers:

### 5.1 Static Dependency Graph (always active)

The dependency graph encodes what each step requires — prerequisite steps and data conditions. It is **not a plan** — it's static knowledge used to validate any step the user requests.

```python
# scagent/dependencies.py
DEPENDENCIES = {
    'qc_metrics':       {'requires': [],              'requires_x': ['raw_counts']},
    'normalize':        {'requires': [],              'requires_x': ['raw_counts']},
    'hvg':              {'requires': ['normalize'],   'requires_x': ['log_normalized']},
    'pca':              {'requires': ['hvg']},
    'neighbors':        {'requires': ['pca']},
    'clustering':       {'requires': ['neighbors']},
    'umap':             {'requires': ['neighbors']},
    'annotation':       {'requires': ['clustering']},
    'pseudobulk_de':    {'requires': ['annotation'],  'requires_data': ['raw_counts', 'condition_key']},
    'trajectory':       {'requires': ['neighbors']},
    ...
}
```

When the researcher asks for something, the agent:
1. Checks `check_prerequisites(goal, state)` — can this run now?
2. If not, calls `plan_steps(goal, state)` — what's the minimal set of steps to get there?
3. Fills the gaps, then runs the requested step.

This works with or without a DAG. A researcher can ask "run trajectory" mid-atlas-analysis and the dependency graph validates it.

### 5.2 Optional Analysis DAG (when the user wants a plan)

When the researcher states a clear goal ("I want to build a cell atlas"), the system generates a paradigm-specific DAG. This DAG is **advisory** — it tracks progress and suggests next steps, but doesn't block out-of-plan requests.

The DAG can be modified dynamically:
- `dag.add_step(step, after="clustering")` — add trajectory to an atlas DAG
- `dag.mark_precomputed_from_state(state)` — mark already-completed steps when onboarding foreign data

### Default DAG for `disease_vs_healthy` paradigm:

For each paradigm, the orchestrator generates an analysis DAG. Steps are nodes; edges are data dependencies. The researcher can modify the DAG interactively.

### Default DAG for `disease_vs_healthy` paradigm:

```
cellranger_multi
       │
       ▼
  qc_filtering ──────────────────┐
       │                         │
       ▼                         ▼
  normalization            doublet_detection
       │                         │
       ▼                         │
  feature_selection ◄────────────┘
       │
       ▼
      pca
       │
       ▼
  batch_correction (if multi-batch)
       │
       ▼
  neighbor_graph
       │
       ├──────────────────┐
       ▼                  ▼
   clustering           umap
       │                  │
       ▼                  │
  marker_de ◄─────────────┘
       │
       ▼
  cell_annotation
       │
       ├──────────────────┐
       ▼                  ▼
  cross_condition_de   composition_analysis
       │                  │
       ▼                  │
  pathway_enrichment     │
       │                  │
       └──────┬───────────┘
              ▼
         report_generation
```

### Default DAG for `developmental_trajectory` paradigm:

```
cellranger_multi → qc → normalization → hvg → pca → integration
       │
       ├→ neighbor_graph → clustering → annotation
       │         │
       │         ├→ paga_trajectory
       │         ├→ monocle3_pseudotime
       │         └→ rna_velocity (scvelo)
       │
       └→ marker_de → tf_analysis (scenic)
```

---

## 5b. Data Inspection (`scagent/inspector.py`)

When scAgent encounters data it didn't create — a researcher's pre-processed `.h5ad`, scBench eval data, or a collaborator's shared file — it needs to figure out what state the data is in before doing anything.

### The Problem

`adata.X` could be raw integer counts, log-normalized floats, or z-score scaled values with negatives. Attempting to normalize already-scaled data produces NaN values and crashes. The inspector runs **once** at data load and produces a comprehensive state report.

### X State Detection (priority-ordered)

| Check | Signal | State |
|-------|--------|-------|
| `(X < 0).any()` | Negative values | `scaled` |
| `allclose(sample, round(sample))` | Integer values | `raw_counts` |
| `'log1p' in adata.uns` | Scanpy breadcrumb | `log_normalized` |
| `X.max() < 20` | Value range heuristic | `log_normalized` (with warning) |
| else | Unknown transform | `transformed_unknown` |

### Breadcrumb Detection

Each Scanpy step leaves specific keys. The inspector checks for all of them:

| Step | Breadcrumb | Reliable? |
|------|-----------|-----------|
| QC metrics | `obs['n_genes_by_counts']` | ✅ |
| Normalization | `uns['log1p']` | ✅ |
| HVG selection | `var['highly_variable']` | ✅ |
| PCA | `obsm['X_pca']` | ✅ |
| Neighbors | `uns['neighbors']` + `obsp['connectivities']` | ✅ |
| UMAP | `obsm['X_umap']` | ✅ |
| Clustering | `obs['leiden*']` / `obs['louvain*']` | ⚠️ fuzzy match |
| Cell types | `obs['cell_type']` / `obs['annotation']` / etc. | ⚠️ fuzzy match |
| DE results | `uns['rank_genes_groups']` | ✅ |

### Raw Counts Recovery

The inspector also locates raw counts for re-processing:
1. `layers['counts']` (Scanpy convention)
2. `layers['raw_counts']` (some pipelines)
3. `X` itself (if state is `raw_counts`)

### Usage

```python
from scagent.inspector import inspect_adata, summarize_state
from scagent.dependencies import ensure_ready_for

state = inspect_adata(adata)
print(summarize_state(state))
# → "6,420 cells × 19,918 genes (mouse), X: z-score scaled,
#    raw counts in layers['raw_counts'], clusters in obs['leiden']..."

# Bring X to the right state for HVG selection
ensure_ready_for(adata, state, needs="log_normalized")
```

### Onboarding Flow

```
User: [drops a pre-processed .h5ad]

scAgent: I've inspected your data:
  6,420 cells × 19,918 genes (mouse)
  X is z-score scaled, raw counts in layers['raw_counts']
  Clusters: obs['leiden'] (12 clusters)
  Cell types: obs['cell_type'] (8 types)
  
  What would you like to do?

User: Find DE genes between the CAF subtypes

scAgent: I see 'CAF' in your cell type labels. I'll:
  1. Subset to CAF cells
  2. Re-cluster to find subtypes (raw counts → normalize → PCA → cluster)
  3. Run DE between subtypes
  Proceed?
```

---

## 6. Provenance Graph (PROV-JSONLD)

Every step executed produces a provenance record. The format follows W3C PROV-O, serialized as JSON-LD. This is the "portable knowledge format" — a colleague can reconstruct exactly what happened.

```jsonc
{
  "@context": {
    "prov": "http://www.w3.org/ns/prov#",
    "sca": "https://scagent.dev/ontology/v1#",
    "xsd": "http://www.w3.org/2001/XMLSchema#"
  },
  "@graph": [
    {
      // The output data entity
      "@id": "sca:adata_clustered_v3",
      "@type": "prov:Entity",
      "prov:wasGeneratedBy": "sca:activity_leiden_001",
      "prov:wasDerivedFrom": "sca:adata_neighbors_v2",
      "sca:snapshot_hash": "sha256:a3f8c2...",
      "sca:branch": "main",
      "sca:step_index": 7
    },
    {
      // The activity that produced it
      "@id": "sca:activity_leiden_001",
      "@type": "prov:Activity",
      "prov:used": "sca:adata_neighbors_v2",
      "prov:wasAssociatedWith": "sca:tool_leiden_clustering",
      "prov:startedAtTime": "2026-04-11T14:23:00Z",
      "prov:endedAtTime": "2026-04-11T14:23:12Z",

      "sca:parameters": {
        "resolution": 1.0,
        "random_state": 42,
        "n_iterations": -1,
        "key_added": "leiden_1.0"
      },
      "sca:tool_version": "scanpy==1.10.4",
      "sca:python_version": "3.11.8",
      "sca:triggered_by": "user_request",
      "sca:user_prompt": "cluster the cells at resolution 1.0"
    },
    {
      // The tool (agent)
      "@id": "sca:tool_leiden_clustering",
      "@type": ["prov:Agent", "prov:SoftwareAgent"],
      "sca:tool_id": "leiden_clustering",
      "sca:registry_version": "0.1.0"
    }
  ]
}
```

**Portability:** This file is ~2KB per step. A full analysis of 15 steps is ~30KB. Send it to a colleague → they can see exactly what was done, with what parameters, in what order, and reproduce it.

The provenance graph is append-only. Branches create new provenance chains. The graph is queryable:
- "What parameters were used for clustering?" → look up all activities with tool_id=leiden_clustering
- "What changed between branch A and branch B?" → diff the provenance chains
- "Reproduce this analysis" → replay all activities in order with recorded parameters

---

## 7. Branched Analysis (Version Control for Data)

The researcher can fork any point in the analysis to explore alternatives. This is NOT git for code — it's version control for analysis state.

```
main:     qc → norm → hvg → pca → neighbors → leiden(0.8) → annotation → DE
                                                    │
                                                    ├── branch: "high_res"
                                                    │   leiden(2.0) → annotation → DE
                                                    │
                                                    └── branch: "harmony_corrected"
                                                        harmony → neighbors → leiden(1.0) → annotation → DE
```

### Implementation

Each "state" is an AnnData snapshot identified by a content hash. We don't copy the entire AnnData each time — we use a **copy-on-write** strategy:

```
.scagent/
├── project.json                    # experiment context (minSCe)
├── provenance.jsonld               # full provenance graph
├── dag.json                        # current analysis DAG
├── branches/
│   ├── main/
│   │   ├── HEAD                    # pointer to current state hash
│   │   ├── states/
│   │   │   ├── a3f8c2.h5ad        # snapshot after clustering
│   │   │   ├── b7d1e4.h5ad        # snapshot after annotation
│   │   │   └── ...
│   │   └── log.jsonl               # step-by-step execution log
│   ├── high_res/
│   │   ├── HEAD
│   │   ├── parent_branch: "main"
│   │   ├── fork_point: "a3f8c2"   # forked from this state
│   │   └── states/
│   │       └── c9e3f1.h5ad
│   └── harmony_corrected/
│       └── ...
├── tools/                          # tool registry (local overrides)
├── knowledge/                      # MarkerDB user databases
│   ├── cellmarker2.xlsx            # optional CellMarker 2.0
│   └── panglaodb.tsv               # optional PanglaoDB
└── memory/                         # mempalace data
    └── palace.db
```

### Branch Commands (REPL)

```
> /branch list
  main (active) — 12 steps, last: differential_expression
  high_res — 3 steps, forked from main@leiden_clustering
  harmony_corrected — 5 steps, forked from main@pca

> /branch create "no_mito_filter" from main@qc_filtering
  Created branch 'no_mito_filter' from state a1b2c3

> /branch compare high_res main
  Diverges at: leiden_clustering
  high_res: resolution=2.0 → 24 clusters → ...
  main: resolution=0.8 → 12 clusters → ...
  Shared cell types: 9/12 overlap
  high_res splits: CD4_naive → CD4_naive_1, CD4_naive_2

> /branch merge high_res into main
  ⚠ Cannot auto-merge analysis branches. Use /branch switch to activate.
```

---

## 8. Memory & Context Management

### 8.1 MemPalace Integration (Long-term Chat History)

MemPalace stores all researcher–agent conversations verbatim in ChromaDB, organized by the palace metaphor:

- **Wing:** The project (e.g., "PBMC aging study")
- **Hall:** Type of interaction (analysis decisions, parameter discussions, interpretation)
- **Room:** Specific topic (e.g., "clustering resolution debate", "B cell marker validation")

On each REPL turn, the agent can search MemPalace for relevant past context:
```python
# Before generating a response, the agent searches
results = mempalace.search("why did we choose resolution 1.0", wing="pbmc_aging")
# → Returns the conversation from 3 weeks ago where the researcher
#   discussed marker genes and decided 0.8 was too coarse
```

### 8.2 MarkerDB + Best Practice References (Domain Knowledge)

Knowledge is structured, not a generic knowledge graph. Two components:

**MarkerDB** (`knowledge.py`) — 3-tier marker gene lookup:
1. Built-in canonical markers (always available, no download)
2. CellTypist model extraction (61 trained models, ~50 MB each)
3. User-supplied databases (CellMarker 2.0 xlsx, PanglaoDB TSV in `.scagent/knowledge/`)

```
> /knowledge query "NK cells"
  Canonical: NKG7, GNLY, PRF1, GZMB, KLRD1, NCAM1 (CD56)
  CellTypist (Immune_All_Low): NKG7 (rank=1), GNLY (rank=2), GZMB (rank=3)
  Confidence: HIGH (sources agree)

> /knowledge validate --markers GNLY,NKG7,GZMB,PRF1 --label "NK cells"
  Match: 4/6 canonical markers, overlap=67%, confidence=HIGH
```

**Best practice references** (`best_practices/reference/`, 12 Markdown files) — distilled literature-backed guidance per analysis step, sourced from Heumos et al. 2023, 10x Genomics Analysis Guide 2025, and sc-best-practices.org. Loaded into agent context via skills on demand. Backed by two full PDFs.

### 8.3 Context Window Management

Modern LLMs provide 200K+ token context windows (Claude: 200K, Gemini: 1M+). scAgent should
use this capacity aggressively — the whole point is that biologists accumulate context over months.

The strategy is **tiered**: a small always-pinned layer, a large working memory layer that grows
with the session, and a compaction mechanism that fires only when genuinely necessary.

```
┌──────────────────────────────────────────────────────────────┐
│              CONTEXT WINDOW (200K tokens available)          │
│                                                              │
│  TIER 1 — ALWAYS PINNED (~3K tokens):                        │
│  ├─ System prompt + paradigm-specific rules        ~1,200    │
│  ├─ Experiment context (minSCe JSON)                 ~500    │
│  ├─ Current branch + state summary                   ~500    │
│  ├─ Current DAG position + next steps                ~400    │
│  └─ Active tool schema (for current step)            ~400    │
│                                                              │
│  TIER 2 — WORKING MEMORY (grows up to ~150K tokens):         │
│  ├─ Full conversation history (all turns this session)       │
│  │   Typical: ~1,500 tokens/turn × 40 turns = ~60K          │
│  ├─ Provenance chain for current branch              ~2-5K   │
│  ├─ Last tool outputs (QC plots, DE tables, UMAPs)   ~5-10K  │
│  ├─ Knowledge graph retrievals (papers, markers)     ~3-5K   │
│  └─ MemPalace retrievals (prior session context)    ~3-5K    │
│                                                              │
│  TIER 3 — OVERFLOW / RETRIEVAL (~47K reserved):              │
│  ├─ Full provenance history (all branches)                   │
│  ├─ Extended marker gene tables                              │
│  ├─ Paper excerpts from knowledge graph                      │
│  └─ Earlier conversation turns (retrieved on demand)         │
│                                                              │
│  COMPACTION TRIGGER:                                         │
│  When Tier 2 exceeds 150K tokens:                            │
│  → Summarize conversation turns older than 20 turns          │
│  → Mine full verbatim history to MemPalace                   │
│  → Compress tool outputs to key findings only                │
│  → Keep Tier 1 + last 20 turns + summaries of older turns    │
│                                                              │
│  CROSS-SESSION CONTINUITY:                                   │
│  On session start:                                           │
│  → Load Tier 1 (always)                                      │
│  → Load last session summary from MemPalace (~2K)            │
│  → Load recent provenance chain (~2K)                        │
│  → Retrieve relevant past discussions via MemPalace search   │
│  Total cold-start context: ~10K tokens                       │
│  Researcher can /recall "topic" to pull more past context    │
└──────────────────────────────────────────────────────────────┘
```

**Key principle:** Don't compress prematurely. With 200K tokens, a researcher can have ~80-100 back-and-forths
before compaction is needed. For a typical analysis session (20-30 turns), the full conversation fits
comfortably. Compaction is a safety valve, not a routine operation.

For 1M-token models (Gemini), compaction may never be needed within a single session.

---

## 9. REPL Interface Design

```
$ scagent init --paradigm disease_vs_healthy

╔══════════════════════════════════════════════════════════════╗
║  scAgent v0.1 — Single-Cell Analysis Agent                  ║
║  Project: pbmc_aging_study                                  ║
║  Paradigm: Disease vs. Healthy                              ║
║  Platform: 10x Chromium GEM-X 3' v4                         ║
║  Branch: main (step 0/12)                                   ║
╚══════════════════════════════════════════════════════════════╝

Commands:
  /status          Current analysis state
  /plan            Show analysis DAG
  /branch          Branch management
  /provenance      View lineage
  /knowledge       Query knowledge graph
  /tools           List available tools
  /export          Export provenance / report
  /help            Full command list

You: I've uploaded my Cell Ranger outputs for 4 donors. Can you run QC?

scAgent: I see 4 filtered_feature_bc_matrix directories. Let me load them
and run QC. Here's what I'll do:

  1. Load all 4 samples into a merged AnnData
  2. Compute QC metrics (n_genes, n_counts, pct_mito, pct_ribo)
  3. Show you the distributions so we can set thresholds together

[Executing... ████████████ 100%]

Here are the QC distributions:

  ┌─ Genes per cell ──────────────────┐
  │ Donor 1: median=3,274 (IQR: 2,100-4,500)
  │ Donor 2: median=3,102 (IQR: 1,900-4,200)
  │ ...
  └────────────────────────────────────┘

  ┌─ % Mitochondrial ────────────────┐
  │ All donors: median=4.2% (95th pctl: 8.7%)
  │ Suggested threshold: 10%
  └────────────────────────────────────┘

Based on your PBMC samples, I suggest:
  • Min genes: 200  |  Max genes: 5000
  • Max mito: 10%   |  Min cells per gene: 3

These are conservative defaults. Want to adjust?

You: Looks good, apply those. Also use scrublet for doublet detection.

scAgent: Applying QC filters + Scrublet doublet detection...

  ┌─ Filtering Results ──────────────┐
  │ Donor 1: 5,710 → 5,423 cells (95.0% retained)
  │ Donor 2: 5,200 → 4,891 cells (94.1% retained)
  │ ...
  │ Doublets removed: 312 total (2.8%)
  └────────────────────────────────────┘

  State saved: branch=main, step=qc_filtering, hash=a1b2c3
  Provenance recorded: filter_cells(min_genes=200, max_genes=5000, max_mito=0.10)
                       + scrublet(expected_doublet_rate=0.06)

  Next suggested step: normalization (SCTransform recommended for Chromium)
```

---

## 10. scGPT Integration (Future)

Leave a clean interface for foundation model embeddings:

```jsonc
// tools/foundation_models/scgpt_embeddings.json
{
  "tool_id": "scgpt_embeddings",
  "name": "scGPT Cell Embeddings",
  "category": "dimensionality_reduction",
  "framework": "scgpt",
  "status": "planned",                       // not yet active

  "requires": {
    "gpu": true,
    "vram_gb": 16,
    "model_checkpoint": "scGPT_human"
  },

  "valid_after": ["normalization"],
  "valid_before": ["clustering", "cell_annotation"],

  "parameters": {
    "model": {
      "type": "string",
      "options": ["scGPT_human", "scGPT_mouse", "scGPT_brain"],
      "default": "scGPT_human"
    },
    "batch_size": { "type": "int", "default": 64 }
  },

  "outputs": {
    "cell_embeddings": "adata.obsm['X_scgpt']"
  },

  "notes": "Alternative to PCA for downstream clustering/annotation. scCluBench (2025) found foundation model embeddings underperform task-specific clustering methods but excel at classification. Consider using for annotation transfer rather than clustering."
}
```

When hardware becomes available, the tool just gets activated in the registry. No architecture changes needed.

---

## 11. How This Beats scBench (By Design, Not By Overfitting)

| scBench Failure Mode | scAgent Design Response |
|---------------------|--------------------------|
| Agent doesn't know what platform the data is from | **Experiment context** explicitly encodes platform, chemistry, kit |
| Agent doesn't know data is pre-processed | **Inspector** detects X state (raw/normalized/scaled), finds raw counts in layers |
| Agent uses wrong normalization for the data type | **Inspector + dependency graph** check X state before any transformation; `ensure_ready_for()` handles recovery |
| Agent treats cells as independent replicates in cross-condition DE | **Dependency graph** for pseudobulk requires condition column + raw counts |
| Agent picks arbitrary clustering resolution | **Parameter guidance** from literature baked into tool schema; resolution sweep is default |
| Agent can't identify cell types | **Knowledge graph** loaded with marker databases; annotation validated against DE markers |
| Agent forgets what it did 5 steps ago | **Provenance graph** always available; current state summary always in context |
| Agent hallucinates gene names | **Validation rules** check gene names against the loaded AnnData var_names |
| Agent crashes on unfamiliar data structures | **All data is AnnData** — the system enforces a single data format |

### The Key Insight

scBench's hardest tasks are **cell typing** (48% best) and **DE** (41% best). These fail because general agents:
1. Don't know what tissue they're working with
2. Don't know what markers to expect
3. Don't enforce pseudobulk for cross-sample comparisons
4. Pick arbitrary parameters without justification

scAgent addresses all four by construction — the experiment context tells the agent the tissue, the knowledge graph provides markers, the DAG enforces pseudobulk, and the tool registry constrains parameters with literature-backed defaults.

---

## 12. Technology Stack

| Component | Technology | Role |
|-----------|-----------|------|
| CLI/REPL | pi (platform) + scAgent system prompt | User interface |
| LLM backbone | Claude via API (through pi) | Reasoning, code generation, interpretation |
| Data format | AnnData (.h5ad) | Universal single-cell data container |
| Analysis core | Scanpy (Python) | Primary analysis framework |
| R bridge | rpy2 or subprocess | For Seurat, DESeq2, edgeR when needed |
| Chat memory | MemPalace (ChromaDB) | Long-term conversation storage |
| Knowledge | MarkerDB + best_practices/reference/ | Marker genes + literature-backed guidance |
| Provenance | Custom PROV-JSONLD writer | Lineage tracking |
| State snapshots | AnnData .h5ad + SHA256 | Branched analysis states |
| Tool registry | JSON schemas | Tool capabilities and constraints |
| Visualization | Matplotlib/Scanpy plots + Rich terminal | Inline plot display |

---

## 13. Implementation Roadmap

### Phase 1: Foundation (Weeks 1-3)
- [ ] Project initialization (`scagent init`)
- [ ] Experiment context schema (minSCe JSON)
- [ ] Tool registry with 10 core tools (QC through clustering)
- [ ] Basic REPL with single-branch execution
- [ ] Provenance recording (PROV-JSONLD)
- [ ] AnnData state management (save/load snapshots)

### Phase 2: Intelligence (Weeks 4-6)
- [ ] Paradigm router (maps experiment type → analysis DAG)
- [ ] LLM integration for natural language → tool invocation
- [ ] Parameter guidance from tool schemas
- [ ] Validation checks on tool outputs
- [ ] MemPalace integration for chat persistence

### Phase 3: Exploration (Weeks 7-9)
- [ ] Branch/fork/compare workflow
- [ ] Multiple resolution sweep with comparison view
- [ ] Expand MarkerDB with tissue-specific databases
- [ ] Cell annotation pipeline with validation
- [ ] DE pipeline (markers + pseudobulk)

### Phase 4: Polish (Weeks 10-12)
- [ ] Full DAG for all 9 paradigms
- [ ] Context compaction strategy
- [ ] Export: provenance file, methods section, reproducibility package
- [ ] Tool registry expanded to 30+ tools
- [ ] Documentation and onboarding flow

### Phase 5: Validation (Post-launch)
- [ ] Run against scBench canonical 30 (Chromium subset: 7 evals)
- [ ] Request full 394-eval access from LatchBio
- [ ] Identify paradigm/task gaps and iterate
- [ ] scGPT integration when GPU available

---

## 14. Open Design Questions

| Question | Options | Recommendation |
|----------|---------|----------------|
| Should branching copy full AnnData or use deltas? | Full copy (simple), delta/diff (space-efficient), lazy checkpointing (only save on branch) | Start with lazy checkpointing — only snapshot when branching. Use full copy for simplicity. |
| How to handle R tools (DESeq2, Seurat)? | rpy2 in-process, subprocess with file exchange, separate R server | Subprocess with .h5ad → .rds conversion. Cleanest isolation. |
| Should the agent auto-advance the DAG or wait for user? | Auto-advance, always wait, suggest + confirm | Suggest + confirm. Show next step, ask "proceed?" |
| How granular should provenance be? | Per-tool, per-function-call, per-line-of-code | Per-tool invocation. Each tool is one provenance Activity. |
| Context compaction: summarize or discard? | Summarize to MemPalace, discard old turns, hierarchical summary | Summarize + mine to MemPalace. Keep summaries queryable. |
