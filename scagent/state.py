"""State management for scAgent: lazy-checkpointed, branch-aware AnnData snapshots.

AnnData files are large (193 MB sparse, 1.8 GB after PCA). This module
only writes to disk on specific triggers: branch fork, branch switch,
session end, or explicit save.  During normal single-branch work, AnnData
lives in memory only.

Usage::

    from pathlib import Path
    from scagent.state import StateManager

    sm = StateManager(Path(".scagent"))
    sm.save_snapshot(adata, step_name="after_pca")
    sm.create_branch("high_res")
    adata = sm.switch_branch("high_res")
"""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import anndata as ad

# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

DISK_WARNING_BYTES = 10 * 1024 * 1024 * 1024  # 10 GB


@dataclass
class BranchInfo:
    """Summary of a branch's state."""

    name: str
    is_active: bool
    head_hash: str | None  # None if no snapshot yet
    head_step: str | None
    head_step_index: int
    parent_branch: str | None  # None for root (main)
    fork_hash: str | None
    fork_step: str | None
    snapshot_count: int
    disk_size_bytes: int


@dataclass
class StateInfo:
    """Current in-memory state summary."""

    branch: str
    hash: str | None
    step_name: str | None
    step_index: int
    has_unsaved_changes: bool


# ---------------------------------------------------------------------------
# StateManager
# ---------------------------------------------------------------------------


class StateManager:
    """Branch-aware, lazy-checkpointed AnnData state manager.

    Parameters
    ----------
    project_dir
        Path to the ``.scagent`` directory.  Created on first save if
        it does not exist.
    """

    def __init__(self, project_dir: Path | str) -> None:
        self._dir = Path(project_dir)
        self._branches_dir = self._dir / "branches"
        self._state_path = self._dir / "state.json"

        # In-memory tracking
        self._active_branch: str = "main"
        self._dirty: bool = False
        self._dirty_step: str | None = None
        self._head_hash: str | None = None
        self._head_step: str | None = None
        self._head_step_index: int = -1

        # Load existing state or initialise
        if self._state_path.exists():
            self._load_state()
        else:
            # Ensure main branch directory exists
            self._branch_dir("main").mkdir(parents=True, exist_ok=True)
            (self._branch_dir("main") / "snapshots").mkdir(exist_ok=True)
            self._save_state()

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def active_branch(self) -> str:
        return self._active_branch

    def current_state(self) -> StateInfo:
        """Return a summary of the current in-memory state."""
        return StateInfo(
            branch=self._active_branch,
            hash=self._head_hash,
            step_name=self._head_step,
            step_index=self._head_step_index,
            has_unsaved_changes=self._dirty,
        )

    # ------------------------------------------------------------------
    # Snapshot operations
    # ------------------------------------------------------------------

    def save_snapshot(
        self,
        adata: ad.AnnData,
        step_name: str,
        *,
        step_index: int | None = None,
        branch: str | None = None,
    ) -> str:
        """Write AnnData to disk and return the content hash (8-char hex).

        Parameters
        ----------
        adata
            The AnnData object to snapshot.
        step_name
            Human-readable step name (e.g., ``"after_pca"``).
        step_index
            Numeric step position.  Auto-incremented from HEAD if *None*.
        branch
            Target branch.  Defaults to the active branch.
        """
        branch = branch or self._active_branch
        snap_dir = self._branch_dir(branch) / "snapshots"
        snap_dir.mkdir(parents=True, exist_ok=True)

        # Write to a temp file first, then hash, then rename
        with tempfile.NamedTemporaryFile(
            dir=snap_dir, suffix=".h5ad", delete=False
        ) as tmp:
            tmp_path = Path(tmp.name)

        adata.write(tmp_path)
        file_hash = _hash_file(tmp_path)

        final_path = snap_dir / f"{file_hash}.h5ad"
        if final_path.exists():
            tmp_path.unlink()  # identical snapshot already on disk
        else:
            tmp_path.rename(final_path)

        # Resolve step index
        if step_index is None:
            head = self._read_head(branch)
            step_index = (head.get("step_index", -1) + 1) if head else 0

        # Update head pointer
        self._write_head(
            branch,
            hash=file_hash,
            step_name=step_name,
            step_index=step_index,
        )

        # Clean up old snapshots on this branch (keep only current + fork refs)
        self._gc_branch(branch, keep={file_hash})

        # Update in-memory tracking if this is the active branch
        if branch == self._active_branch:
            self._head_hash = file_hash
            self._head_step = step_name
            self._head_step_index = step_index
            self._dirty = False
            self._dirty_step = None

        self._check_disk_warning()
        return file_hash

    def load_snapshot(self, branch: str | None = None) -> ad.AnnData:
        """Load the HEAD snapshot of a branch.

        Parameters
        ----------
        branch
            Branch to load from.  Defaults to the active branch.

        Raises
        ------
        FileNotFoundError
            If the branch has no snapshots.
        """
        branch = branch or self._active_branch
        head = self._read_head(branch)
        if head is None or head.get("hash") is None:
            raise FileNotFoundError(
                f"No snapshot on branch '{branch}'. Run some analysis first."
            )
        return self.load_snapshot_by_hash(head["hash"], branch=branch)

    def load_snapshot_by_hash(
        self, hash: str, *, branch: str | None = None
    ) -> ad.AnnData:
        """Load a specific snapshot by its content hash."""
        # Search across branches if no branch specified
        branches_to_search = (
            [branch] if branch else self._list_branch_names()
        )
        for b in branches_to_search:
            path = self._branch_dir(b) / "snapshots" / f"{hash}.h5ad"
            if path.exists():
                return ad.read_h5ad(path)
        raise FileNotFoundError(f"Snapshot '{hash}' not found.")

    # ------------------------------------------------------------------
    # Dirty tracking
    # ------------------------------------------------------------------

    def mark_dirty(self, step_name: str, step_index: int | None = None) -> None:
        """Mark that the in-memory AnnData has unsaved changes."""
        self._dirty = True
        self._dirty_step = step_name
        if step_index is not None:
            self._head_step_index = step_index

    # ------------------------------------------------------------------
    # Branch lifecycle
    # ------------------------------------------------------------------

    def create_branch(
        self,
        name: str,
        *,
        from_branch: str | None = None,
        adata: ad.AnnData | None = None,
    ) -> str:
        """Create a new branch forked from *from_branch*.

        If the fork source has unsaved changes, *adata* must be provided
        so the fork point can be snapshotted.

        Returns the fork-point snapshot hash.
        """
        from_branch = from_branch or self._active_branch

        if not self._branch_exists(name):
            pass  # good, doesn't exist yet
        else:
            raise ValueError(f"Branch '{name}' already exists")

        if not self._branch_exists(from_branch):
            raise ValueError(f"Source branch '{from_branch}' does not exist")

        # Ensure the fork point is on disk
        head = self._read_head(from_branch)
        if head is None or head.get("hash") is None:
            # No snapshot yet — we need the adata to create one
            if adata is None:
                raise ValueError(
                    f"Branch '{from_branch}' has no snapshot. "
                    "Pass adata= to snapshot the current state."
                )
            fork_hash = self.save_snapshot(
                adata, step_name="fork_point", branch=from_branch
            )
            head = self._read_head(from_branch)
        elif self._dirty and from_branch == self._active_branch and adata is not None:
            # Active branch has unsaved changes — save before forking
            fork_hash = self.save_snapshot(
                adata,
                step_name=self._dirty_step or "fork_point",
                branch=from_branch,
            )
            head = self._read_head(from_branch)
        else:
            fork_hash = head["hash"]

        # Create branch directory
        new_branch_dir = self._branch_dir(name)
        new_branch_dir.mkdir(parents=True, exist_ok=True)
        (new_branch_dir / "snapshots").mkdir(exist_ok=True)

        # Copy fork-point snapshot to new branch
        src = self._branch_dir(from_branch) / "snapshots" / f"{fork_hash}.h5ad"
        dst = new_branch_dir / "snapshots" / f"{fork_hash}.h5ad"
        if src.exists() and not dst.exists():
            shutil.copy2(src, dst)

        # Write parent.json
        parent_meta = {
            "parent_branch": from_branch,
            "fork_hash": fork_hash,
            "fork_step": head.get("step_name", ""),
            "fork_step_index": head.get("step_index", -1),
            "forked_at": datetime.now(timezone.utc).isoformat(),
        }
        (new_branch_dir / "parent.json").write_text(
            json.dumps(parent_meta, indent=2), encoding="utf-8"
        )

        # Write head.json pointing to fork snapshot
        self._write_head(
            name,
            hash=fork_hash,
            step_name=head.get("step_name", "fork_point"),
            step_index=head.get("step_index", 0),
        )

        return fork_hash

    def switch_branch(
        self, name: str, *, adata: ad.AnnData | None = None
    ) -> ad.AnnData:
        """Switch to a different branch.

        Saves the current branch's state if *adata* is provided and
        there are unsaved changes, then loads the target branch's HEAD.

        Returns the loaded AnnData.
        """
        if not self._branch_exists(name):
            raise ValueError(f"Branch '{name}' does not exist")

        # Save current branch if dirty
        if self._dirty and adata is not None:
            self.save_snapshot(
                adata,
                step_name=self._dirty_step or self._head_step or "unknown",
            )

        # Switch
        self._active_branch = name
        head = self._read_head(name)
        if head:
            self._head_hash = head.get("hash")
            self._head_step = head.get("step_name")
            self._head_step_index = head.get("step_index", -1)
        else:
            self._head_hash = None
            self._head_step = None
            self._head_step_index = -1

        self._dirty = False
        self._dirty_step = None
        self._save_state()

        return self.load_snapshot(name)

    def delete_branch(self, name: str, *, force: bool = False) -> None:
        """Delete a branch and its snapshots.

        Raises
        ------
        ValueError
            If *name* is the active branch, is ``"main"``, or has child
            branches (unless *force* is True).
        """
        if name == self._active_branch:
            raise ValueError(f"Cannot delete the active branch '{name}'")
        if name == "main":
            raise ValueError("Cannot delete the 'main' branch")
        if not self._branch_exists(name):
            raise ValueError(f"Branch '{name}' does not exist")

        # Check for child branches
        if not force:
            children = self._child_branches(name)
            if children:
                raise ValueError(
                    f"Branch '{name}' has child branches: {children}. "
                    "Use force=True to delete anyway."
                )

        shutil.rmtree(self._branch_dir(name))

    def list_branches(self) -> list[BranchInfo]:
        """List all branches with their metadata."""
        result: list[BranchInfo] = []
        for name in self._list_branch_names():
            head = self._read_head(name)
            parent = self._read_parent(name)
            snap_dir = self._branch_dir(name) / "snapshots"
            snapshots = list(snap_dir.glob("*.h5ad")) if snap_dir.exists() else []
            disk_size = sum(f.stat().st_size for f in snapshots)

            result.append(
                BranchInfo(
                    name=name,
                    is_active=(name == self._active_branch),
                    head_hash=head.get("hash") if head else None,
                    head_step=head.get("step_name") if head else None,
                    head_step_index=head.get("step_index", -1) if head else -1,
                    parent_branch=parent.get("parent_branch") if parent else None,
                    fork_hash=parent.get("fork_hash") if parent else None,
                    fork_step=parent.get("fork_step") if parent else None,
                    snapshot_count=len(snapshots),
                    disk_size_bytes=disk_size,
                )
            )
        return result

    # ------------------------------------------------------------------
    # Session lifecycle
    # ------------------------------------------------------------------

    def save_on_exit(self, adata: ad.AnnData, step_name: str, **kw: Any) -> str:
        """Save current state for cross-session persistence. Returns hash."""
        return self.save_snapshot(adata, step_name, **kw)

    def load_on_start(self) -> ad.AnnData | None:
        """Load the active branch's HEAD snapshot, or *None* if empty."""
        try:
            adata = self.load_snapshot(self._active_branch)
            head = self._read_head(self._active_branch)
            if head:
                self._head_hash = head.get("hash")
                self._head_step = head.get("step_name")
                self._head_step_index = head.get("step_index", -1)
            return adata
        except FileNotFoundError:
            return None

    # ------------------------------------------------------------------
    # Garbage collection
    # ------------------------------------------------------------------

    def gc_snapshots(self) -> int:
        """Remove unreferenced snapshots across all branches. Returns bytes freed."""
        # Collect all referenced hashes
        referenced: set[str] = set()
        for name in self._list_branch_names():
            head = self._read_head(name)
            if head and head.get("hash"):
                referenced.add(head["hash"])
            parent = self._read_parent(name)
            if parent and parent.get("fork_hash"):
                referenced.add(parent["fork_hash"])

        # Also keep fork hashes referenced by child branches
        for name in self._list_branch_names():
            parent = self._read_parent(name)
            if parent and parent.get("fork_hash"):
                referenced.add(parent["fork_hash"])

        freed = 0
        for name in self._list_branch_names():
            snap_dir = self._branch_dir(name) / "snapshots"
            if not snap_dir.exists():
                continue
            for f in snap_dir.glob("*.h5ad"):
                file_hash = f.stem
                if file_hash not in referenced:
                    freed += f.stat().st_size
                    f.unlink()
        return freed

    def total_disk_usage(self) -> int:
        """Total bytes used by all snapshots."""
        total = 0
        for name in self._list_branch_names():
            snap_dir = self._branch_dir(name) / "snapshots"
            if snap_dir.exists():
                total += sum(f.stat().st_size for f in snap_dir.glob("*.h5ad"))
        return total

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _branch_dir(self, name: str) -> Path:
        return self._branches_dir / name

    def _branch_exists(self, name: str) -> bool:
        return self._branch_dir(name).is_dir()

    def _list_branch_names(self) -> list[str]:
        if not self._branches_dir.exists():
            return ["main"]
        names = sorted(
            d.name for d in self._branches_dir.iterdir() if d.is_dir()
        )
        # Ensure main is first
        if "main" in names:
            names.remove("main")
            names.insert(0, "main")
        return names

    def _child_branches(self, name: str) -> list[str]:
        """Return branches that were forked from *name*."""
        children = []
        for b in self._list_branch_names():
            parent = self._read_parent(b)
            if parent and parent.get("parent_branch") == name:
                children.append(b)
        return children

    def _read_head(self, branch: str) -> dict | None:
        path = self._branch_dir(branch) / "head.json"
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
        return None

    def _write_head(
        self,
        branch: str,
        *,
        hash: str,
        step_name: str,
        step_index: int,
    ) -> None:
        branch_dir = self._branch_dir(branch)
        branch_dir.mkdir(parents=True, exist_ok=True)
        data = {
            "hash": hash,
            "step_name": step_name,
            "step_index": step_index,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        (branch_dir / "head.json").write_text(
            json.dumps(data, indent=2), encoding="utf-8"
        )

    def _read_parent(self, branch: str) -> dict | None:
        path = self._branch_dir(branch) / "parent.json"
        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))
        return None

    def _save_state(self) -> None:
        self._dir.mkdir(parents=True, exist_ok=True)
        data = {
            "active_branch": self._active_branch,
            "updated": datetime.now(timezone.utc).isoformat(),
        }
        self._state_path.write_text(
            json.dumps(data, indent=2), encoding="utf-8"
        )

    def _load_state(self) -> None:
        data = json.loads(self._state_path.read_text(encoding="utf-8"))
        self._active_branch = data.get("active_branch", "main")

        # Load head info for active branch
        head = self._read_head(self._active_branch)
        if head:
            self._head_hash = head.get("hash")
            self._head_step = head.get("step_name")
            self._head_step_index = head.get("step_index", -1)

    def _gc_branch(self, branch: str, keep: set[str]) -> None:
        """Remove snapshots on *branch* that aren't in *keep* or referenced."""
        snap_dir = self._branch_dir(branch) / "snapshots"
        if not snap_dir.exists():
            return

        # Also keep fork hashes referenced by child branches
        for child in self._child_branches(branch):
            parent = self._read_parent(child)
            if parent and parent.get("fork_hash"):
                keep.add(parent["fork_hash"])

        for f in snap_dir.glob("*.h5ad"):
            if f.stem not in keep:
                f.unlink()

    def _check_disk_warning(self) -> None:
        """Print a warning if total disk usage exceeds the threshold."""
        total = self.total_disk_usage()
        if total > DISK_WARNING_BYTES:
            gb = total / (1024 ** 3)
            print(
                f"⚠ Branches are using {gb:.1f} GB of disk. "
                "Consider deleting unused branches."
            )


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------


def _hash_file(path: Path) -> str:
    """SHA-256 hash of a file, returning the first 8 hex characters."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(8192)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()[:8]
