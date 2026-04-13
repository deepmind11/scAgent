---
name: provenance
description: Record and query the full lineage of every analysis step. Use after every tool call to record provenance, or when the user asks what was done, what parameters were used, or how to reproduce the analysis.
---

# Skill: Provenance Tracking

Record and query the full lineage of every analysis step.

## When to Use

- **After every tool call** — record provenance automatically.
- **When the user asks** "what did we do?", "what parameters?", "how to reproduce?", `/provenance`, `/history`.

## Recording — After Each Tool Call

After calling **any** `scagent.tools.*` function, immediately record provenance:

```python
from pathlib import Path
from scagent.provenance import ProvenanceGraph, record_step

graph = ProvenanceGraph(Path(".scagent"))

# Call the tool
result = filter_cells(adata, min_genes=200, max_genes=5000, max_pct_mito=10.0)

# Record — always pass the user's original prompt
record_step(graph, result, user_prompt="filter cells with standard thresholds")
```

Rules:
1. **Always record.** Even if the user didn't ask. Provenance is automatic.
2. **Pass the user prompt.** This links the activity to the researcher's intent.
3. **One graph per session.** Create `ProvenanceGraph(Path(".scagent"))` once at session start. Reuse it for every `record_step()` call.
4. **Don't construct PROV-JSONLD yourself.** `record_step()` handles extraction, timestamping, versioning, and serialization.

## Querying — When the User Asks

```python
# Show full analysis history as a Markdown table
print(graph.summary())

# List all activities (optionally filter)
activities = graph.list_activities(tool_id="leiden_clustering")

# Generate a reproduction plan
plan = graph.replay_plan()
# → [("load_10x_h5", {"filename": "..."}), ("filter_cells", {"min_genes": 200, ...}), ...]

# Compare branches (when branching is available)
diff = graph.diff("main", "high_res")
```

### What to Show the User

- For `/provenance` or "what did we do?" → use `graph.summary()` which returns a Markdown table.
- For "what parameters did we use for clustering?" → use `graph.list_activities(tool_id="leiden_clustering")`.
- For "how can I reproduce this?" → use `graph.replay_plan()` and present it as an ordered list.
- For "what's different between these branches?" → use `graph.diff(branch_a, branch_b)`.

## Custom / Ad-Hoc Analysis Steps

When executing code that does NOT come from a `scagent.tools.*` function — custom computations, external packages, researcher-provided scripts — use `record_custom()`:

```python
from scagent.provenance import record_custom

record_custom(
    graph,
    description="mito/ribo ratio column",
    code="adata.obs['mito_ribo'] = adata.obs['pct_counts_mt'] / adata.obs['pct_counts_ribo']",
    user_prompt="compute the ratio of mito to ribo genes",
    effects={"added_obs_columns": ["mito_ribo"]},
)
```

Rules:
- **Always record custom steps.** If you wrote and executed code that modifies the AnnData, record it.
- **Include the exact code.** Paste the code string verbatim — this is what makes it reproducible.
- **Describe the effects.** What columns were added? How many cells were removed? This goes in `effects`.
- Custom steps appear as `tool_id="custom"` in the provenance graph and replay plan.

## File Location

Provenance is stored at `.scagent/provenance.jsonld`. This file is:
- **Append-only** — steps are never modified or deleted.
- **Auto-saved** — written after every `record_step()` call.
- **Portable** — send it to a colleague for full analysis reproducibility.
- **~2KB per step** — a 15-step analysis is ~30KB.
