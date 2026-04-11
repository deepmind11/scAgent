# Chunk 4 Design Plan: Provenance Writer

**Status:** Draft for review
**Depends on:** Chunks 1–3 (tool registry + tool wrappers with `provenance` dicts)
**Feeds into:** Chunk 5 (branching — provenance chains per branch), Chunk 6 (DAG — provenance validates step ordering)

---

## 0. Current State (What We Already Have)

The good news: **every tool wrapper already returns provenance data.** Every function in `scagent/tools/` returns a result dict with a `"provenance"` key:

```python
# Example from scagent/tools/qc.py — filter_cells()
return {
    "metrics": { ... },
    "plots": [],
    "provenance": {
        "tool_id": "filter_cells",
        "parameters": {
            "min_genes": min_genes,
            "max_genes": max_genes,
            "max_pct_mito": max_pct_mito,
            "min_counts": min_counts,
        },
        "cells_before": cells_before,
        "cells_after": cells_after,
    },
    "warnings": [],
}
```

This exists across **all 15 tool wrapper functions** (loading, qc ×3, doublets, normalize, hvg, pca, neighbors, umap, harmony, leiden ×3, markers, annotation ×2).

The bad news: **nothing collects them.** The provenance dicts are returned, printed by the agent, and lost. No persistence, no linking, no serialization.

Chunk 4 bridges this gap.

---

## 1. What Chunk 4 Builds

One Python module (`scagent/provenance.py`) that:

1. **Accepts** the raw provenance dicts that tool wrappers already return
2. **Enriches** them with timestamps, software versions, entity linking (input → output), and step ordering
3. **Serializes** to PROV-JSONLD (W3C PROV-O compliant)
4. **Persists** to `.scagent/provenance.jsonld` (append-only)
5. **Queries** — "what parameters for clustering?", "list all steps", "show full lineage"
6. **Diffs** — compare two provenance chains (prep for Chunk 5 branching)
7. **Replays** — emit an ordered list of (tool_id, parameters) for reproduction

Plus one skill (`.pi/skills/provenance/SKILL.md`) that teaches the agent to record provenance after each step and respond to `/provenance` queries.

---

## 2. Module Design

### 2.1 File: `scagent/provenance.py`

```
ProvenanceGraph
├── __init__(project_dir: Path)
├── record(tool_result: dict, *, input_hash: str, output_hash: str, user_prompt: str, branch: str) -> str  # returns activity_id
├── get_activity(activity_id: str) -> dict
├── get_entity(entity_id: str) -> dict
├── list_activities(branch: str = None, tool_id: str = None) -> list[dict]
├── get_chain(branch: str = "main") -> list[dict]  # ordered steps
├── diff(branch_a: str, branch_b: str) -> dict  # what diverges
├── replay_plan(branch: str = "main") -> list[tuple[str, dict]]  # [(tool_id, params), ...]
├── serialize() -> dict  # full PROV-JSONLD
├── save() -> Path  # writes to .scagent/provenance.jsonld
├── load(path: Path) -> "ProvenanceGraph"  # class method
└── summary(branch: str = "main") -> str  # human-readable summary
```

### 2.2 Internal Data Model

Internally, the graph stores three types of nodes (matching W3C PROV-O):

```python
@dataclass
class Entity:
    """A data state — an AnnData snapshot at a point in the analysis."""
    id: str                    # e.g., "sca:adata_filtered_v1"
    hash: str                  # SHA-256 of the .h5ad file (from Chunk 5, placeholder for now)
    branch: str                # "main", "high_res", etc.
    step_index: int            # position in the branch's chain
    generated_by: str          # activity_id that created this
    derived_from: str | None   # entity_id of the input (None for first entity)

@dataclass
class Activity:
    """A tool invocation — one step in the analysis."""
    id: str                    # e.g., "sca:activity_filter_cells_001"
    tool_id: str               # from the provenance dict
    parameters: dict           # from the provenance dict
    extra: dict                # tool-specific extras (cells_before, n_clusters, etc.)
    used: str                  # entity_id of input
    generated: str             # entity_id of output
    started_at: str            # ISO timestamp
    ended_at: str              # ISO timestamp
    user_prompt: str           # what the user asked
    software_versions: dict    # scanpy, python, scagent versions

@dataclass
class Agent:
    """A software tool in the registry."""
    id: str                    # e.g., "sca:tool_filter_cells"
    tool_id: str
    registry_version: str
```

### 2.3 ID Generation

IDs follow the architecture doc's convention:

```
Entity:   sca:adata_{step_name}_v{step_index}       e.g., sca:adata_filtered_v2
Activity: sca:activity_{tool_id}_{counter:03d}       e.g., sca:activity_leiden_001
Agent:    sca:tool_{tool_id}                         e.g., sca:tool_leiden_clustering
```

The counter is global across the graph (monotonically increasing). This ensures IDs are unique even if the same tool is called multiple times (e.g., leiden at different resolutions).

---

## 3. PROV-JSONLD Serialization

Output format matches architecture §6 exactly:

```jsonc
{
  "@context": {
    "prov": "http://www.w3.org/ns/prov#",
    "sca": "https://scagent.dev/ontology/v1#",
    "xsd": "http://www.w3.org/2001/XMLSchema#"
  },
  "@graph": [
    // Entity nodes
    {
      "@id": "sca:adata_loaded_v0",
      "@type": "prov:Entity",
      "prov:wasGeneratedBy": "sca:activity_load_10x_h5_001",
      "sca:snapshot_hash": "sha256:...",
      "sca:branch": "main",
      "sca:step_index": 0
    },
    // Activity nodes
    {
      "@id": "sca:activity_load_10x_h5_001",
      "@type": "prov:Activity",
      "prov:used": "sca:input_file",
      "prov:wasAssociatedWith": "sca:tool_load_10x_h5",
      "prov:startedAtTime": "2026-04-11T14:00:00Z",
      "prov:endedAtTime": "2026-04-11T14:00:03Z",
      "sca:parameters": { "filename": "filtered_feature_bc_matrix.h5", "gex_only": true },
      "sca:tool_version": "scanpy==1.12.1",
      "sca:python_version": "3.12.x",
      "sca:scagent_version": "0.1.0",
      "sca:triggered_by": "user_request",
      "sca:user_prompt": "load the PBMC dataset"
    },
    // Agent (tool) nodes
    {
      "@id": "sca:tool_load_10x_h5",
      "@type": ["prov:Agent", "prov:SoftwareAgent"],
      "sca:tool_id": "load_10x_h5",
      "sca:registry_version": "0.1.0"
    }
    // ... more nodes ...
  ]
}
```

### 3.1 Size Budget

Per the architecture: ~2KB per step, ~30KB for a full 15-step analysis. We should validate this during testing.

---

## 4. Integration with Existing Tool Wrappers

**No changes to tool wrappers needed.** The wrappers already return the `provenance` dicts. We just need a function that extracts the provenance from a tool result and feeds it to the graph.

```python
# Convenience function — the agent (or a future orchestrator) calls this after each tool
def record_step(
    graph: ProvenanceGraph,
    tool_result: dict,         # the full return value from any scagent/tools/ function
    *,
    input_hash: str = "",      # from StateManager (Chunk 5), empty for now
    output_hash: str = "",     # from StateManager (Chunk 5), empty for now
    user_prompt: str = "",     # what the user asked
    branch: str = "main",
) -> str:
    """Extract provenance from a tool result dict and record it in the graph.
    
    Returns the activity_id.
    """
```

This is a thin adapter that:
1. Pulls `tool_result["provenance"]["tool_id"]` and `tool_result["provenance"]["parameters"]`
2. Separates parameters from tool-specific extras (cells_before, n_clusters, etc.)
3. Stamps with current time + software versions
4. Calls `graph.record()`

### 4.1 Software Version Capture

Captured once at graph initialization and attached to every activity:

```python
def _capture_versions() -> dict:
    return {
        "scanpy": sc.__version__,
        "anndata": ad.__version__,
        "python": platform.python_version(),
        "scagent": scagent.__version__,
    }
```

### 4.2 Separating Parameters from Extras

The raw provenance dicts from tool wrappers mix parameters with metadata. For example, `filter_cells` returns:

```python
"provenance": {
    "tool_id": "filter_cells",
    "parameters": { "min_genes": 200, "max_genes": 5000, ... },
    "cells_before": 11769,    # ← this is NOT a parameter, it's a metric
    "cells_after": 10834,
}
```

The `record_step()` function should separate these:
- **`parameters`** → goes into `sca:parameters` (needed for replay)
- **Everything else** (`cells_before`, `cells_after`, `n_clusters`, `file_hash`, etc.) → goes into `sca:extras` (informational, not needed for replay)

Rule: `tool_id` and `parameters` are reserved keys. Everything else in the provenance dict is an extra.

---

## 5. Storage

### 5.1 File Location

```
.scagent/
└── provenance.jsonld          # the full PROV-JSONLD graph
```

Note: `.scagent/` doesn't exist yet. Chunk 4 will create it when first recording. Chunk 5 will add `branches/`, Chunk 6 will add `project.json` and `dag.json`.

### 5.2 Persistence Strategy

- **Save on every record.** Each `record()` call appends to the internal graph, then `save()` writes the full serialized JSONLD. This is simple and safe — a 30KB file rewrite is cheap.
- **Append-only semantics.** Activities and entities are never modified or deleted. A branch fork creates new entities that reference the fork point.
- **Load on init.** `ProvenanceGraph(project_dir)` loads existing `.scagent/provenance.jsonld` if it exists, or starts empty.

### 5.3 No Database

The architecture mentions this is ~30KB for 15 steps. A flat JSON-LD file is fine. No SQLite, no ChromaDB, no over-engineering. If we ever hit 1000+ steps, we can revisit.

---

## 6. Query API

These methods are for the agent (via skill) and for future `/provenance` REPL commands.

### 6.1 `list_activities(branch=None, tool_id=None)`

```python
>>> graph.list_activities()
[
    {"step": 0, "tool": "load_10x_h5", "time": "14:00:03", "output": "adata_loaded_v0"},
    {"step": 1, "tool": "calculate_qc_metrics", "time": "14:00:05", "output": "adata_qc_v1"},
    {"step": 2, "tool": "filter_cells", "time": "14:00:06", "output": "adata_filtered_v2"},
    ...
]

>>> graph.list_activities(tool_id="leiden_clustering")
[
    {"step": 7, "tool": "leiden_clustering", "params": {"resolution": 1.0}, ...},
]
```

### 6.2 `get_chain(branch)`

Returns the ordered sequence of activities for a branch — the full lineage from data load to current state.

### 6.3 `summary(branch)`

Returns a human-readable Markdown summary the agent can include in its response:

```
## Analysis Provenance (branch: main, 8 steps)

| # | Step | Tool | Key Params | Duration |
|---|------|------|------------|----------|
| 0 | Load data | load_10x_h5 | file=pbmc10k.h5 | 2.8s |
| 1 | QC metrics | calculate_qc_metrics | species=human | 1.2s |
| 2 | Filter cells | filter_cells | min_genes=200, max_mito=10% | 0.4s |
| 3 | Filter genes | filter_genes | min_cells=3 | 0.2s |
| 4 | Doublets | scrublet | rate=0.06 | 3.1s |
| 5 | Normalize | log_normalize | target_sum=10000 | 0.5s |
| 6 | HVG | highly_variable_genes | n=2000 | 0.8s |
| 7 | PCA | pca | n_comps=50 | 1.4s |

Total: 10.4s | Software: scanpy 1.12.1, scagent 0.1.0
```

### 6.4 `diff(branch_a, branch_b)`

Compares two branches. Returns:

```python
{
    "diverge_point": "sca:adata_neighbors_v6",   # last shared entity
    "diverge_step": 6,
    "branch_a_steps": [                           # steps unique to branch_a
        {"tool": "leiden_clustering", "params": {"resolution": 0.8}},
        ...
    ],
    "branch_b_steps": [                           # steps unique to branch_b
        {"tool": "leiden_clustering", "params": {"resolution": 2.0}},
        ...
    ],
    "parameter_diffs": [                          # same tool, different params
        {"tool": "leiden_clustering", "param": "resolution", "a": 0.8, "b": 2.0},
    ],
}
```

This is a Chunk 5 feature in spirit, but the diff logic lives in the provenance module since it compares provenance chains.

### 6.5 `replay_plan(branch)`

Returns an ordered list of `(tool_id, parameters)` tuples — everything needed to reproduce the analysis from scratch:

```python
>>> graph.replay_plan()
[
    ("load_10x_h5", {"filename": "filtered_feature_bc_matrix.h5", "gex_only": True}),
    ("calculate_qc_metrics", {"species": "human"}),
    ("filter_cells", {"min_genes": 200, "max_genes": 5000, "max_pct_mito": 10.0}),
    ("filter_genes", {"min_cells": 3}),
    ("scrublet", {"expected_doublet_rate": 0.06, "random_state": 0}),
    ("log_normalize", {"target_sum": 10000}),
    ("highly_variable_genes", {"n_top_genes": 2000, "flavor": "seurat_v3"}),
    ("pca", {"n_comps": 50, "random_state": 0}),
]
```

---

## 7. Skill: `.pi/skills/provenance/SKILL.md`

Teaches the agent two things:

### 7.1 Recording (After Each Tool Call)

After calling any `scagent.tools.*` function, the agent should record provenance:

```python
from scagent.provenance import ProvenanceGraph, record_step

graph = ProvenanceGraph(Path(".scagent"))
# ... call tool ...
result = filter_cells(adata, min_genes=200, max_genes=5000, max_pct_mito=10.0)
record_step(graph, result, user_prompt="filter cells with standard thresholds")
```

The skill should make clear:
- **Always record.** Even if the user didn't ask. Provenance is automatic.
- **Pass the user prompt.** This links the activity to the researcher's intent.
- **No need to construct the PROV-JSONLD yourself.** `record_step()` handles everything.

### 7.2 Querying (When the User Asks)

When the user asks about provenance (`/provenance`, "what did we do?", "what parameters did we use?", "how can I reproduce this?"):

```python
# Show analysis history
print(graph.summary())

# Show specific tool usage
activities = graph.list_activities(tool_id="leiden_clustering")

# Generate reproduction plan
plan = graph.replay_plan()
```

---

## 8. Testing

### 8.1 Unit Tests: `tests/test_provenance.py`

| Test | What It Checks |
|------|----------------|
| `test_record_single_step` | Record one tool result, verify activity/entity/agent are created |
| `test_record_chain` | Record 5 steps in sequence, verify entity linking (each derived_from previous) |
| `test_serialize_jsonld` | Serialize to PROV-JSONLD, validate structure against the architecture spec |
| `test_save_load_roundtrip` | Save to file, load from file, verify graph is identical |
| `test_list_activities` | Record multiple steps, query by tool_id, verify filtering |
| `test_summary` | Record a chain, check that summary() produces a markdown table |
| `test_replay_plan` | Record a chain, verify replay_plan() returns correct (tool_id, params) tuples |
| `test_diff_branches` | Record two divergent branches, verify diff identifies the divergence point |
| `test_record_step_adapter` | Use `record_step()` with a real tool result dict (from qc.py), verify extraction |
| `test_size_budget` | Record 15 steps, check serialized size is ≤ 50KB |
| `test_idempotent_save` | Save twice without recording, verify file doesn't change |
| `test_empty_graph` | Query an empty graph, verify no crashes |

### 8.2 Integration Test: `tests/test_provenance_integration.py`

Run the real pipeline (load → qc → normalize → hvg → pca → neighbors → leiden) on the PBMC 10k dataset, recording provenance after each step. Verify:

1. The provenance.jsonld file is valid JSON-LD
2. It contains 7+ activities in correct order
3. Entity chain is linked (each wasDerivedFrom the previous)
4. Parameters match what was actually passed to the tools
5. Software versions are captured
6. `replay_plan()` output matches the actual call sequence
7. `summary()` produces a readable table
8. File size is within budget

---

## 9. What This Chunk Does NOT Build

Keeping scope clear — these are deferred:

| Feature | Deferred To | Notes |
|---------|-------------|-------|
| `.scagent/branches/` directory structure | Chunk 5 | Provenance supports `branch` field, but branch management is Chunk 5 |
| AnnData snapshot hashing | Chunk 5 | `input_hash` and `output_hash` are placeholders (empty strings) until Chunk 5 |
| DAG validation ("is this step valid here?") | Chunk 6 | Provenance records what happened; DAG validates what should happen |
| `/provenance` REPL command | Chunk 6 or later | The skill covers agent-driven queries; REPL commands are future |
| Methods section export ("generate a Methods paragraph") | Chunk 12 | Uses `replay_plan()` + templates. Polish phase. |

---

## 10. Implementation Checklist

```
Chunk 4: Provenance Writer
├── [ ] 4.1  scagent/provenance.py — ProvenanceGraph class
│     ├── [ ] Entity, Activity, Agent dataclasses
│     ├── [ ] __init__(project_dir) — load or create
│     ├── [ ] record() — add an activity + entity pair
│     ├── [ ] record_step() — adapter for tool result dicts
│     ├── [ ] _capture_versions() — scanpy, anndata, python, scagent
│     ├── [ ] get_activity(), get_entity()
│     ├── [ ] list_activities(branch, tool_id)
│     ├── [ ] get_chain(branch)
│     ├── [ ] summary(branch) — markdown table
│     ├── [ ] diff(branch_a, branch_b)
│     ├── [ ] replay_plan(branch)
│     ├── [ ] serialize() — PROV-JSONLD dict
│     ├── [ ] save() / load()
│     └── [ ] _generate_id() — entity/activity/agent ID generation
├── [ ] 4.2  .pi/skills/provenance/SKILL.md
│     ├── [ ] Recording instructions (after every tool call)
│     ├── [ ] Query instructions (summary, list, replay)
│     └── [ ] Example agent interaction
├── [ ] 4.3  tests/test_provenance.py — 12 unit tests
├── [ ] 4.4  tests/test_provenance_integration.py — end-to-end with real pipeline
├── [ ] 4.5  Update scagent/__init__.py — export ProvenanceGraph
├── [ ] 4.6  Create .scagent/ directory on first use (mkdir -p)
└── [ ] 4.7  Verify: run full PBMC pipeline → valid provenance.jsonld
```

**Estimated output:**
- `scagent/provenance.py` — ~350-450 lines
- `.pi/skills/provenance/SKILL.md` — ~100 lines
- `tests/test_provenance.py` — ~250 lines
- `tests/test_provenance_integration.py` — ~80 lines

---

## 11. Open Questions (For Review)

| # | Question | Options | Recommendation |
|---|----------|---------|----------------|
| 1 | **Where does `record_step()` live?** | Method on `ProvenanceGraph`, or standalone function in the module | **Standalone function** — it's an adapter between tool results and the graph. Keeps the graph class clean. |
| 2 | **Should we auto-save on every `record()`?** | Yes (safe, simple), No (explicit `save()` call required) | **Yes.** The file is small, and we never want to lose provenance. The agent might crash mid-session. |
| 3 | **How to handle `input_hash` / `output_hash` before Chunk 5?** | Skip hashes, use empty string, use a placeholder like "not_tracked" | **Empty string.** Chunk 5 will populate. The schema allows it. Don't over-engineer. |
| 4 | **Should the provenance graph be a singleton per session?** | Singleton, passed around, created per use | **Passed around** (or created once per session in the skill). No global state. |
| 5 | **Step naming: use `tool_id` or a human-readable step name?** | tool_id only, step name from DAG (Chunk 6), both | **`tool_id` for now.** Chunk 6 can add DAG-derived step names later. |
