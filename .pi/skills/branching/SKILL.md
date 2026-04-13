---
name: branching
description: Manage parallel analysis paths — fork, switch, compare, and delete branches. Use when the user wants to try different parameters, compare methods, backtrack, or explore alternative analysis paths.
---

# Skill: Branched Analysis

Manage parallel analysis paths — fork, switch, compare, delete.

## When to Suggest Branching

- User wants to **try different parameters**: "try resolution 2.0" → fork, try on new branch.
- User wants to **compare methods**: "compare Harmony vs. Scanorama" → fork, run each on a branch.
- User wants to **backtrack**: "go back to before clustering" → fork from that point, never rewind in-place.

**Rule: never rewind on the same branch. Always branch first.** The current state is preserved on the original branch.

## How to Use

```python
from pathlib import Path
from scagent.state import StateManager

sm = StateManager(Path(".scagent"))

# --- Session start: load previous state ---
adata = sm.load_on_start()  # returns None if no prior state

# --- After each tool call: mark dirty ---
result = run_leiden(adata, resolution=1.0)
sm.mark_dirty(step_name="after_leiden", step_index=8)
# No disk write happens here — just tracking.

# --- Create a branch (does NOT switch to it) ---
sm.create_branch("high_res", adata=adata)

# --- Switch to the new branch ---
adata = sm.switch_branch("high_res", adata=adata)
# This saves current branch to disk, loads high_res snapshot.

# --- Work on the new branch ---
result = run_leiden(adata, resolution=2.0)
sm.mark_dirty(step_name="after_leiden_high_res")

# --- Switch back ---
adata = sm.switch_branch("main", adata=adata)

# --- List branches ---
for b in sm.list_branches():
    marker = "* " if b.is_active else "  "
    parent = f", forked from {b.parent_branch}@{b.fork_step}" if b.parent_branch else ""
    mb = b.disk_size_bytes / 1024 / 1024
    print(f"{marker}{b.name:20s} — step {b.head_step_index} ({b.head_step}), {mb:.0f} MB{parent}")

# --- Session end: persist state ---
sm.save_on_exit(adata, step_name="after_annotation")
```

## What to Tell the User

- **After fork**: "Created branch 'high_res' from your current state at step 7 (neighbors)."
- **After switch**: "Switched to 'main'. You're at step 9 (annotation)."
- **On list**: show the branch table.
- **On backtrack request**: "I'll create a branch from step 5 to explore from there. Your current work stays on main."
- **On disk warning**: "Branches are using X GB. Consider deleting unused branches."

## Key Constraints

- AnnData is **NOT** saved to disk after every tool call. Only on fork, switch, session end, or explicit save.
- `create_branch` does **NOT** auto-switch. Call `switch_branch` explicitly after.
- Pass `adata=` to `create_branch` and `switch_branch` when the in-memory state has unsaved changes.
- Never delete the active branch or `main`.
