# Skill: Analysis DAG

Follow the paradigm-aware analysis plan, track progress, suggest next steps.

## When to Use

- **Before each tool call** — check what the DAG says to do next.
- **After each tool call** — mark the step as done.
- **When the user asks** "what's next?", "show progress", `/plan`, `/status`.

## How to Use

```python
from pathlib import Path
from scagent.dag import AnalysisDAG

dag = AnalysisDAG.load(Path(".scagent"))

# Suggest next step
step = dag.next_step()
# → DAGStep(id="normalize", name="Normalize", tool_id="log_normalize", ...)

# After executing the tool, mark done
dag.complete_step("normalize")
dag.save(Path(".scagent"))

# Show progress
print(dag.summary())

# Skip a step with reason
dag.skip_step("batch_correction", reason="single-sample experiment")
```

## Behavior Rules

1. **Suggest, don't auto-run.** Show the next step and ask: "Next step is normalization. Proceed?" The researcher confirms.
2. **Out-of-order requests**: if the researcher asks for something that has unmet dependencies, explain: "Clustering requires a neighbor graph first. Want me to run that?"
3. **Skipping**: if the researcher wants to skip a step, record the reason. The DAG continues from the next available step.
4. **Custom steps**: if the researcher asks for something not in the DAG, execute it and record with `record_custom()`. The DAG position doesn't change.
5. **Show progress on request**: use `dag.summary()` for a markdown table with ✅/⬜/⏭️ status.

## What to Tell the User

- **On suggest**: "Next step: **Normalize** (log-normalization with target_sum=10,000). Proceed?"
- **On complete**: "✅ Normalization done. 3/12 steps complete. Next: highly variable genes."
- **On skip**: "⏭️ Skipping batch correction (single-sample experiment). Next: neighbor graph."
- **On progress**: show the full summary table.
