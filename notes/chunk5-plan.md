# Chunk 5 Design Plan: State Manager (Branching)

**Status:** Draft for review
**Depends on:** Chunk 4 (provenance — done, provides entity hashes + fork_branch)
**Feeds into:** Chunk 6 (DAG — needs to know current branch + step position)

---

## 0. The Size Problem

AnnData file sizes for the PBMC 10k dataset (measured):

| Stage | Size | Why |
|-------|------|-----|
| After load (sparse) | **193 MB** | 11,769 × 33,538 sparse matrix |
| After QC/filter (sparse) | **193 MB** | Slightly fewer cells, still sparse |
| After HVG (sparse) | **194 MB** | Added `.var` columns |
| After PCA + neighbors (dense) | **1,804 MB** | `sc.pp.scale()` densifies the matrix |

If we snapshot after every step on a 12-step pipeline, that's **5–15 GB per branch.** With 3 branches exploring different parameters, we'd need **15–45 GB** of disk for one project.

**Conclusion:** We cannot snapshot every step. The architecture doc already anticipated this — it recommends "lazy checkpointing: only snapshot when branching."

---

## 1. Design: Lazy Checkpointing

### When to Save

Snapshots are written to disk only when necessary:

| Trigger | Why |
|---------|-----|
| **Branch fork** | Freeze the fork point so the new branch has a base to work from |
| **Branch switch** | Save current branch state before loading another |
| **Session end** | Persist state so the next session can resume |
| **Explicit save** | User says "save" or `/save` |

During normal single-branch work, AnnData lives **in memory only.** No disk writes between steps.

### When to Load

| Trigger | What's Loaded |
|---------|---------------|
| **Session start** | Latest snapshot of the active branch |
| **Branch switch** | Latest snapshot of the target branch |
| **Explicit load** | A specific snapshot by hash |

### What This Means for the Researcher

- Normal workflow: no disk overhead between steps. Fast.
- Creating a branch: one snapshot write (~3–5 seconds for 1.8 GB). Acceptable.
- Switching branches: one save + one load (~6–10 seconds). Acceptable.
- Session end/start: one save + one load. Acceptable.

---

## 2. Storage Layout

```
.scagent/
├── provenance.jsonld              # (Chunk 4, exists)
├── state.json                     # global state metadata
├── branches/
│   ├── main/
│   │   ├── head.json              # pointer to current snapshot
│   │   └── snapshots/
│   │       └── a3f8c2de.h5ad     # latest snapshot (content-hashed name)
│   ├── high_res/
│   │   ├── head.json
│   │   ├── parent.json            # {"branch": "main", "fork_hash": "a3f8c2de"}
│   │   └── snapshots/
│   │       └── b7d1e4f9.h5ad
│   └── harmony_corrected/
│       ├── head.json
│       ├── parent.json
│       └── snapshots/
│           └── c9e3f1a2.h5ad
```

### File Descriptions

**`state.json`** — global metadata:
```jsonc
{
  "active_branch": "main",
  "created": "2026-04-11T10:00:00Z"
}
```

**`branches/{name}/head.json`** — per-branch pointer:
```jsonc
{
  "hash": "a3f8c2de",          // first 8 chars of SHA-256
  "step_name": "leiden_clustering",
  "step_index": 7,
  "timestamp": "2026-04-11T14:30:00Z"
}
```

**`branches/{name}/parent.json`** — only for forked branches:
```jsonc
{
  "parent_branch": "main",
  "fork_hash": "a3f8c2de",     // which snapshot was forked from
  "fork_step": "neighbor_graph",
  "fork_step_index": 6,
  "forked_at": "2026-04-11T15:00:00Z"
}
```

### Content-Addressed Naming

Snapshot files are named by the first 8 characters of their SHA-256 hash. This gives:
- **Identity verification:** you can always re-hash the file to confirm it's the right state.
- **Dedup potential:** if two branches happen to reach the same state, they'd share a hash (unlikely in practice but correct by construction).
- **Provenance linking:** the hash in `head.json` matches the `sca:snapshot_hash` in `provenance.jsonld`.

Computing SHA-256 on a 1.8 GB file takes ~3–5 seconds. This happens only during saves, which are already infrequent.

### Snapshot Retention

Each branch keeps **only its latest snapshot** plus the fork point snapshot (if it's referenced by a child branch). Old snapshots are deleted on save.

Exception: if the user explicitly requests "save a checkpoint here" (future feature), that snapshot is pinned and not garbage-collected.

---

## 3. Module Design

### 3.1 File: `scagent/state.py`

```
StateManager
├── __init__(project_dir: Path)
│
│   # Branch lifecycle
├── create_branch(name: str, from_branch: str = None) -> None
├── switch_branch(name: str) -> AnnData
├── delete_branch(name: str) -> None
├── list_branches() -> list[BranchInfo]
├── active_branch -> str
│
│   # Snapshot operations
├── save_snapshot(adata: AnnData, step_name: str) -> str   # returns hash
├── load_snapshot(branch: str = None) -> AnnData           # loads HEAD
├── load_snapshot_by_hash(hash: str) -> AnnData            # loads specific
│
│   # State queries
├── current_state() -> StateInfo                           # hash, step, branch
├── branch_info(name: str) -> BranchInfo                   # parent, steps, size
│
│   # Session lifecycle
├── save_on_exit(adata: AnnData, step_name: str) -> None   # called on session end
├── load_on_start() -> AnnData | None                      # called on session start
│
│   # Cleanup
└── gc_snapshots() -> int                                  # remove unreferenced snapshots, return bytes freed
```

### 3.2 Data Types

```python
@dataclass
class BranchInfo:
    name: str
    is_active: bool
    head_hash: str | None          # None if no snapshot yet
    head_step: str | None
    head_step_index: int
    parent_branch: str | None      # None for root branch (main)
    fork_hash: str | None
    fork_step: str | None
    snapshot_count: int
    disk_size_bytes: int

@dataclass
class StateInfo:
    branch: str
    hash: str | None
    step_name: str | None
    step_index: int
    has_unsaved_changes: bool      # True if in-memory adata differs from disk
```

### 3.3 Hashing

```python
def _hash_adata(adata: AnnData) -> str:
    """SHA-256 hash of the AnnData, using a temp file write.
    
    Returns the first 8 hex characters (enough to avoid collisions
    within a project — 4 billion possible values).
    """
```

We hash the serialized `.h5ad` file, not the in-memory object. This ensures the hash is deterministic and matches what's on disk.

**Why 8 characters?** Within a single project, we'll have at most ~100 snapshots across all branches. 8 hex chars = 4.3 billion possible values. Collision probability is negligible.

---

## 4. Integration with Provenance (Chunk 4)

Currently, provenance entities have `hash: ""` (empty placeholder). The StateManager fills these in:

```python
# After saving a snapshot
hash = state_manager.save_snapshot(adata, "after_leiden")

# Record provenance with the real hash
record_step(graph, result, output_hash=hash, ...)
```

The `input_hash` comes from the previous step's `output_hash` (already tracked by the provenance entity chain). For the first step, it comes from the input file hash (already captured by `load_10x_h5`).

### Wiring (in the skill or future orchestrator)

```python
# After each tool call:
# 1. Record provenance (already doing this)
record_step(graph, result, output_hash=current_hash, ...)
# 2. Track in-memory state (no disk write unless triggered)
state_manager.mark_dirty(step_name="after_leiden")
```

The `mark_dirty()` call just notes that the in-memory AnnData has changed since the last save. It doesn't write to disk. The actual save happens only at the triggers listed in §1.

---

## 5. Integration with Existing Checkpoints

The existing tool wrappers all have `checkpoint_dir` parameters. These are a simpler, flat checkpointing system (`after_load.h5ad`, `after_filter_cells.h5ad`, etc.) used by `pipeline.py`.

**Decision: keep both, don't break existing.**

- `checkpoint_dir` stays for the pipeline's simple use case (sequential debugging, eval harness).
- `StateManager` is the branch-aware, content-addressed system for interactive use.
- They don't conflict because they write to different directories (`checkpoint_dir` is user-specified; StateManager writes to `.scagent/branches/`).

Long-term, the skill should prefer StateManager. But the tool wrappers don't need to change.

---

## 6. Branch Operations — Detailed Behavior

### 6.1 `create_branch(name, from_branch="main")`

1. Validate: name doesn't exist, from_branch does exist.
2. **Save current state** of `from_branch` to disk (if not already saved). This is the fork point snapshot.
3. Create `.scagent/branches/{name}/` directory.
4. Write `parent.json` with fork metadata.
5. **Copy** (or hardlink) the fork point snapshot to the new branch's `snapshots/` dir.
6. Write `head.json` pointing to that snapshot.
7. Call `provenance.fork_branch(name, from_branch)` to record in provenance graph.
8. **Do NOT switch.** The user is still on `from_branch`.

### 6.2 `switch_branch(name)`

1. Validate: target branch exists.
2. **Save** current branch's in-memory AnnData to disk (if dirty).
3. Update `state.json` → `active_branch = name`.
4. **Load** target branch's HEAD snapshot into memory.
5. Return the loaded AnnData.

### 6.3 `delete_branch(name)`

1. Validate: name is not active branch, name is not "main".
2. Check if any other branch was forked from this one. If yes, refuse (or force flag).
3. Remove `.scagent/branches/{name}/` directory and all snapshots.
4. Note: provenance records for this branch are NOT deleted (append-only graph).

### 6.4 `list_branches()`

Returns `BranchInfo` for each branch:

```
> /branch list
  * main          — step 7 (neighbors), 1 snapshot, 1.8 GB
    high_res      — step 9 (annotation), forked from main@neighbors, 1 snapshot, 1.8 GB
    harmony_test  — step 8 (clustering), forked from main@pca, 1 snapshot, 1.8 GB
```

---

## 7. Skill: `.pi/skills/branching/SKILL.md`

Teaches the agent:

1. **When to suggest branching:**
   - User wants to try different parameters ("try resolution 2.0")
   - User wants to compare normalization methods
   - User wants to backtrack ("let's go back to before clustering")

2. **How to use the StateManager:**
   ```python
   from scagent.state import StateManager
   sm = StateManager(Path(".scagent"))
   
   # Create a branch for exploring high resolution
   sm.create_branch("high_res")
   adata = sm.switch_branch("high_res")
   # ... run tools on adata ...
   
   # Switch back to main
   adata = sm.switch_branch("main")
   
   # Compare
   print(sm.list_branches())
   ```

3. **What to tell the user:**
   - After fork: "I've created branch 'high_res' from your current state. Switching to it now."
   - After switch: "Switched to branch 'main'. You're at step 7 (neighbors)."
   - On list: show the table from §6.4.

---

## 8. Testing

### 8.1 Unit Tests: `tests/test_state.py`

| Test | What It Checks |
|------|----------------|
| `test_save_load_roundtrip` | Save AnnData → load → verify shapes and values match |
| `test_hash_deterministic` | Same AnnData → same hash twice |
| `test_hash_changes_on_mutation` | Modify AnnData → hash changes |
| `test_create_branch` | Create a branch, verify directory structure + parent.json |
| `test_switch_branch` | Save on branch A, switch to B, verify different state |
| `test_switch_saves_dirty` | Modify in-memory adata, switch away → snapshot is updated |
| `test_delete_branch` | Delete branch, verify directory removed |
| `test_delete_active_branch_raises` | Can't delete the branch you're on |
| `test_delete_main_raises` | Can't delete main |
| `test_list_branches` | Create 3 branches, verify list output |
| `test_current_state` | After save, current_state() returns correct hash + step |
| `test_session_save_load` | save_on_exit → load_on_start roundtrip |
| `test_gc_snapshots` | Create old snapshots, gc_snapshots() removes them |
| `test_provenance_integration` | save_snapshot returns hash, usable in record_step |
| `test_branch_with_parent_info` | Forked branch has correct parent + fork metadata |
| `test_multiple_saves_overwrite` | Saving twice on same branch replaces old snapshot |

### 8.2 Integration Test: `tests/test_state_integration.py`

Real PBMC dataset branching scenario:

1. Load data → QC → normalize → PCA → neighbors (on main)
2. Fork "high_res" from main
3. On high_res: leiden(2.0) → markers
4. Switch back to main
5. On main: leiden(0.8) → markers
6. Verify: both branches have different clustering, correct step counts
7. Verify: provenance graph has both branches with correct hashes
8. Verify: disk usage is ~2 snapshots (not 12)

---

## 9. What This Chunk Does NOT Build

| Feature | Deferred To | Notes |
|---------|-------------|-------|
| `/branch compare` with cell-type overlap analysis | Chunk 11+ | Needs annotation + metrics |
| `/branch merge` | Never? | Architecture says "cannot auto-merge." Branch switch is the answer. |
| Incremental / delta snapshots | Future | Would save disk but adds complexity. Full snapshots are simple and correctfor now. |
| Automatic branching ("try 3 resolutions") | Chunk 6+ | DAG-driven parameter sweeps. StateManager provides the primitive. |
| Undo / rollback within a branch | Future | Would require keeping intermediate snapshots. Conflicts with lazy approach. |

---

## 10. Open Questions

| # | Question | Options | Recommendation |
|---|----------|---------|----------------|
| 1 | **Hash the .h5ad file or hash in-memory arrays?** | File hash (deterministic, simple), array hash (faster, no temp file) | **File hash.** Guarantees what's on disk matches the hash. Array hashing has float precision and ordering issues. |
| 2 | **Snapshot on every `record_step` or only on triggers?** | Every step (safe, expensive), triggers only (lazy, efficient) | **Triggers only.** 193MB–1.8GB per save rules out every-step snapshots. |
| 3 | **How to handle "I want to go back 3 steps on main"?** | Re-run from the last snapshot (replay), keep intermediate snapshots, refuse | **Branch, then replay.** Always snapshot the current state first (as a branch or auto-save), then replay from provenance to the target step on a new branch. Never rewind in-place — that would discard the current state. The skill should enforce this: "I'll create a branch from step 7 to explore from there. Your current work stays on main." |
| 4 | **Should `create_branch` auto-switch to the new branch?** | Yes (more intuitive), No (explicit is better) | **No.** Explicit `switch_branch()` after. Keeps operations atomic and predictable. |
| 5 | **Disk budget: should we warn/refuse when snapshots exceed a threshold?** | Hard limit, soft warning, no limit | **Soft warning** at 10 GB total. Print: "Branches are using X GB. Consider deleting unused branches with `/branch delete`." |

---

## 11. Implementation Checklist

```
Chunk 5: State Manager
├── [ ] 5.1  scagent/state.py — StateManager class
│     ├── [ ] BranchInfo, StateInfo dataclasses
│     ├── [ ] __init__(project_dir) — load or create state.json
│     ├── [ ] save_snapshot(adata, step_name) → hash
│     ├── [ ] load_snapshot(branch) → AnnData
│     ├── [ ] load_snapshot_by_hash(hash) → AnnData
│     ├── [ ] _hash_adata(adata) → str (8-char hex)
│     ├── [ ] create_branch(name, from_branch)
│     ├── [ ] switch_branch(name) → AnnData
│     ├── [ ] delete_branch(name)
│     ├── [ ] list_branches() → list[BranchInfo]
│     ├── [ ] active_branch property
│     ├── [ ] current_state() → StateInfo
│     ├── [ ] mark_dirty(step_name) — track unsaved changes
│     ├── [ ] save_on_exit(adata, step_name)
│     ├── [ ] load_on_start() → AnnData | None
│     └── [ ] gc_snapshots() → int (bytes freed)
├── [ ] 5.2  .pi/skills/branching/SKILL.md
│     ├── [ ] When to suggest branching
│     ├── [ ] How to use StateManager
│     └── [ ] What to tell the user
├── [ ] 5.3  tests/test_state.py — 16 unit tests
├── [ ] 5.4  tests/test_state_integration.py — real PBMC branching scenario
├── [ ] 5.5  Wire provenance: save_snapshot hash → record_step output_hash
├── [ ] 5.6  Update scagent/__init__.py — export StateManager
└── [ ] 5.7  Verify: PBMC pipeline with 2 branches → correct state on disk
```

**Estimated output:**
- `scagent/state.py` — ~300–400 lines
- `.pi/skills/branching/SKILL.md` — ~80 lines
- `tests/test_state.py` — ~300 lines
- `tests/test_state_integration.py` — ~100 lines
