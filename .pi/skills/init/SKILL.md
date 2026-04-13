---
name: init
description: Set up the experiment context for a new scAgent project. Use when no project.json exists, the user says init or new project, or the user jumps straight to analysis without context.
---

# Skill: Project Initialization

Set up the experiment context for a new scAgent project. This is the first thing that happens — before any analysis.

## When to Trigger

- **No `.scagent/project.json` exists** and the user wants to start analysis → run init.
- **User says** "init", "new project", "set up", or drops data files → run init.
- **User jumps straight to analysis** without context → pause, run init first. Never silently default.

## The Flow: Infer First, Ask Only What You Must

### Step 1: Load data and auto-detect

```python
from pathlib import Path
from scagent.context import ExperimentContext

ctx = ExperimentContext(Path(".scagent"))

# If data is loaded, infer what we can
inferred = ExperimentContext.infer_from_data(adata)
# inferred may contain: organism, library, platform, _data_shape, _batch_hint
```

### Step 2: Confirm inferred fields

Tell the user what you detected:
> "I've loaded 11,769 cells × 33,538 genes. Gene names suggest **human**. The data appears to be **10x Chromium 3' UMI-based**. Is that correct?"

### Step 3: Ask analysis-critical fields (only what's missing)

**Always ask — these cannot be inferred:**

1. **Tissue**: "What tissue are these cells from?" (e.g., PBMCs, brain cortex, tumor)
2. **Paradigm**: "What's your experimental question?"
   - Characterize cell types in a tissue → `cell_atlas`
   - Compare conditions (disease vs. healthy, treated vs. control) → `disease_vs_healthy`
   - Track a developmental or differentiation trajectory → `developmental_trajectory`
3. **Conditions** (only if paradigm requires them): "Do your samples belong to different conditions? Which samples belong to which group?"
4. **Batch structure** (only if multi-sample): "Were samples processed in separate batches?"

**Don't ask** about donor demographics, dissociation method, sequencing instrument, etc. Record them if the researcher volunteers the info, but don't block on them.

### Step 4: Build and save context

```python
ctx.paradigm = "disease_vs_healthy"
ctx.organism = inferred.get("organism", {"species": "Homo sapiens", "ncbi_taxon": 9606})
ctx.tissue = {"name": "PBMCs"}
ctx.design = {"conditions": ["young", "old"], "batch_key": "sample"}
ctx.samples = [{"id": "donor_1", "condition": "young"}, ...]
ctx.save()
```

### Step 5: Generate and show the DAG

```python
from scagent.dag import AnalysisDAG

dag = AnalysisDAG.from_context(ctx)
dag.save(Path(".scagent"))
print(dag.summary())
```

Show the researcher the plan:
> "Here's your analysis plan for **disease vs. healthy**. It includes pseudobulk DE and batch correction. Want to modify anything before we start?"

## minSCe Checklist Reference

Based on Füllgrabe et al. 2020 — analysis-critical fields:

| Field | How to get it |
|-------|--------------|
| Species | Inferred from gene names, confirm |
| Tissue | Ask |
| Platform / chemistry | Inferred from data format, confirm if unclear |
| Library type (3'/5'/full-length) | Inferred from platform |
| Modalities (GEX only? multiome?) | Ask if unclear |
| Paradigm | Always ask |
| Conditions | Ask if paradigm requires them |
| Batch structure | Ask if multi-sample |
| Cell enrichment / sorting | Ask only if relevant (e.g., immune panel) |

## Mode B: Onboarding Pre-Processed Data

When the user provides a `.h5ad` file that has already been (partially) processed:

### Step 1: Inspect the data

```python
from scagent.inspector import inspect_adata, summarize_state

state = inspect_adata(adata)
print(summarize_state(state))
```

### Step 2: Show what was detected

Tell the user what you found:
> "I've inspected your data: 6,420 cells × 19,918 genes (mouse). X is z-score scaled, raw counts are available in `layers['raw_counts']`. I see Leiden clusters and cell type annotations. What would you like to do?"

### Step 3: Respond to the question

Use the dependency system to figure out what's needed:

```python
from scagent.dependencies import plan_steps, check_prerequisites

# User asks for DE → check what's needed
can_run, missing = check_prerequisites("pseudobulk_de", state)
steps = plan_steps("pseudobulk_de", state)
```

### Step 4: Optionally generate a DAG

If the user wants a structured plan, generate one and mark pre-computed steps:

```python
dag = AnalysisDAG.from_context(ctx)
dag.mark_precomputed_from_state(state)
```

## Memory (automatic)

When ``ExperimentContext`` loads or saves, it automatically creates a ``ProjectMemory`` backed by a ChromaDB palace at ``.scagent/palace/``.  **No manual init is needed.**

Access memory through the context object:

```python
ctx = ExperimentContext(Path(".scagent"))
# ctx.memory is now a ready-to-use ProjectMemory instance

ctx.memory.store_decision("Use resolution 0.8", "NK markers clean", branch="main")
hits = ctx.memory.recall("clustering resolution")
```

On session start, if a palace exists, search for prior context:

```python
hits = ctx.memory.recall("last session summary")
if hits:
    # Tell the user what was recovered
```

## Rules

- **Never silently default to cell_atlas.** If paradigm is unknown, ask — or let the user's question drive the analysis.
- **Paradigm is optional.** The user can work question-by-question without declaring a paradigm. The dependency graph validates steps regardless.
- **Don't re-ask every session.** Once `project.json` is saved, load it. Only ask again if critical fields are missing.
- **Accept partial info.** If the researcher gives a rich prompt, extract what you can and confirm. Don't re-ask what they already told you.
- **Inspect before assuming.** When receiving a foreign `.h5ad`, always run `inspect_adata()` before any preprocessing. Never assume X contains raw counts.
