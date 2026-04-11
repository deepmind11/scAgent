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

## File Location

Provenance is stored at `.scagent/provenance.jsonld`. This file is:
- **Append-only** — steps are never modified or deleted.
- **Auto-saved** — written after every `record_step()` call.
- **Portable** — send it to a colleague for full analysis reproducibility.
- **~2KB per step** — a 15-step analysis is ~30KB.
