---
name: memory
description: Cross-session memory for scAgent using MemPalace. Store and recall analysis decisions, tool results, and discussions across sessions via natural language search.
---

# Skill: Memory (MemPalace)

Cross-session memory for scAgent.  Every analysis decision, tool result, and
discussion is stored in a local ChromaDB palace and searchable by natural
language.

## Prerequisites

MemPalace ships with scAgent.  Verify:

```bash
mempalace --version
```

## Palace layout

| Concept | scAgent mapping |
|---------|----------------|
| **Wing** | Project name (e.g. ``PBMC_aging``) |
| **Room** | Analysis phase: qc, preprocessing, clustering, annotation, de, enrichment, exploration, decisions |
| **Branch** | Metadata tag — the analysis branch (``main``, ``high_res``, …) |

## When to store

- **After every tool run** — call ``store_step(result["provenance"], branch=current_branch)``
- **After an analysis decision** — call ``store_decision(what, why, branch=current_branch)``
- **After a meaningful exchange** — call ``store_exchange(question, response, branch=current_branch)``
- **When the save hook fires** — save everything it asks for, tagged with the current branch

## When to recall

- **Session start** — search for recent context: ``recall("last session summary")``
- **Before answering** something that might have been discussed — ``recall("topic")``
- **When the user asks "didn't we already…"** — always search first
- **Cross-branch** — pass ``branch=None`` to search across all branches

## How to use

```python
from scagent.memory import ProjectMemory

mem = ProjectMemory(".scagent/palace", project_name="PBMC_aging")

# Store
mem.store_step(result["provenance"], branch="main")
mem.store_decision("Use resolution 0.8", "NK markers clean", branch="main")
mem.store_exchange("Why not Wilcoxon?", "Inflates false positives", branch="main")

# Recall — current branch (default)
hits = mem.recall("clustering resolution", branch="main")

# Recall — cross-branch (see what other branches found)
hits = mem.recall("clustering resolution", branch=None)

# Recall — specific other branch
hits = mem.recall("clustering resolution", branch="high_res")
```

## Branching

Branch is just a metadata tag.  The agent decides what to filter:

- **Default**: search the current branch only
- **Cross-branch**: pass ``branch=None`` to see everything
- **Specific branch**: pass ``branch="other"`` to look at one branch

This lets two branches communicate: branch B can check what branch A tried
by searching with ``branch="branch_a"`` or ``branch=None``.

## Hooks

Two hooks auto-trigger saves (install via ``.claude/settings.local.json``):

| Hook | When | What |
|------|------|------|
| ``scagent_save_hook.sh`` | Every 15 human messages | Blocks agent, tells it to save key decisions/steps |
| ``scagent_precompact_hook.sh`` | Before context compaction | Emergency save — save EVERYTHING |

When a hook fires, save to MemPalace using the methods above.  Always tag
with the current branch.

## Compaction protocol

When the **precompact hook** fires, save in this priority order:

1. **Experiment context** — species, tissue, paradigm, conditions, platform
2. **Current branch + DAG position** — what's done, what step is next
3. **Data state** — cell count, cluster count, key `.obs` columns
4. **Key decisions + rationale** — why resolution X, why pseudobulk, rejected alternatives
5. **Open questions** — what the user asked that isn't resolved yet

Use ``store_decision()`` for items 4–5.  Use ``store_exchange()`` for a
short summary of the session so far.

When the **save hook** fires (every 15 exchanges), save only items 4 and
any new tool results (``store_step()``).  Don't re-save context that hasn't
changed.

**After compaction / on session cold-start:**

1. ``recall("experiment context", branch=current)``
2. ``recall("current DAG position", branch=current)``
3. ``recall("recent decisions", branch=current, n_results=5)``
4. Tell the user what was recovered

## What to tell the user

**On session start (if memory exists):**
> "Found context from your previous session. You were clustering at resolution 0.8 with 21 clusters on branch main."

**When memory answers a question:**
> "From our earlier discussion: we tried resolution 1.2 but it over-split NK cells."

**Cross-branch insight:**
> "Branch high_res tried resolution 1.2 — it found 33 clusters but NK markers fragmented."
