"""Cross-session memory for scAgent via MemPalace.

Thin wrapper around MemPalace.  Stores conversation exchanges, provenance
steps, and analysis decisions in a ChromaDB palace.  Branch is a metadata
tag — the caller decides what to filter.

Usage::

    from scagent.memory import ProjectMemory
    mem = ProjectMemory(".scagent/palace", project_name="PBMC_aging")
    mem.store("We chose resolution 0.8 because NK markers were cleaner",
              room="clustering", branch="main")
    hits = mem.recall("resolution choice")
"""

from __future__ import annotations

import hashlib
import logging
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)

COLLECTION_NAME = "mempalace_drawers"

# Canonical rooms for scRNA-seq analysis phases
ROOMS = (
    "qc", "preprocessing", "clustering", "annotation",
    "de", "enrichment", "exploration", "decisions",
)

# Aliases → canonical room
_ROOM_ALIASES: dict[str, str] = {
    "quality_control": "qc", "quality control": "qc", "filtering": "qc",
    "normalize": "preprocessing", "normalization": "preprocessing",
    "hvg": "preprocessing", "pca": "preprocessing", "batch_correction": "preprocessing",
    "leiden": "clustering", "louvain": "clustering", "umap": "clustering",
    "markers": "annotation", "cell_type": "annotation", "celltypist": "annotation",
    "differential_expression": "de", "pseudobulk": "de", "deseq2": "de",
    "gsea": "enrichment", "pathway": "enrichment",
    "decision": "decisions",
}

# tool_id → room
_TOOL_ROOMS: dict[str, str] = {
    "load_10x_h5": "qc", "qc_metrics": "qc",
    "filter_cells": "qc", "filter_genes": "qc", "detect_doublets": "qc",
    "normalize": "preprocessing", "log_transform": "preprocessing",
    "highly_variable_genes": "preprocessing", "scale": "preprocessing",
    "pca": "preprocessing", "batch_correction": "preprocessing",
    "neighbors": "clustering", "leiden": "clustering", "umap": "clustering",
    "rank_genes_groups": "annotation", "annotate_celltypist": "annotation",
    "wilcoxon_markers": "annotation",
    "deseq2_pseudobulk": "de",
    "gsea": "enrichment",
    "custom": "exploration",
}


def normalize_room(room: str) -> str:
    """Map a room name or alias to its canonical form."""
    r = room.lower().strip()
    if r in ROOMS:
        return r
    return _ROOM_ALIASES.get(r, "exploration")


def room_from_tool(tool_id: str) -> str:
    """Infer the room from a tool_id."""
    return _TOOL_ROOMS.get(tool_id, "exploration")


class ProjectMemory:
    """MemPalace-backed project memory.

    Parameters
    ----------
    palace_path
        Directory for the ChromaDB store (e.g. ``.scagent/palace``).
    project_name
        Wing name — typically the project/experiment name.
    """

    def __init__(self, palace_path: str | Path, project_name: str = "scagent"):
        self.palace_path = str(Path(palace_path).resolve())
        self.project_name = project_name
        self._collection = None
        self._session_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    # ------------------------------------------------------------------
    # Lazy collection
    # ------------------------------------------------------------------

    def _col(self):
        if self._collection is None:
            from mempalace.palace import get_collection
            self._collection = get_collection(
                self.palace_path, collection_name=COLLECTION_NAME,
            )
        return self._collection

    # ------------------------------------------------------------------
    # Store (one method — all variants)
    # ------------------------------------------------------------------

    def store(
        self,
        content: str,
        *,
        room: str = "exploration",
        branch: str = "main",
        source_type: str = "general",
        metadata: dict | None = None,
    ) -> str:
        """Store text in the palace.

        Parameters
        ----------
        content
            Verbatim text to store.
        room
            Analysis phase (auto-normalized).
        branch
            Current analysis branch — just a tag, no filtering logic.
        source_type
            ``exchange``, ``provenance``, ``decision``, or ``general``.
        metadata
            Extra key-value pairs to attach.

        Returns
        -------
        Drawer ID.
        """
        room = normalize_room(room)
        col = self._col()
        drawer_id = _make_id(self.project_name, room, content, branch)
        meta = {
            "wing": self.project_name,
            "room": room,
            "branch": branch,
            "source_type": source_type,
            "session_id": self._session_id,
            "filed_at": datetime.now().isoformat(),
            "source_file": f"scagent_{source_type}",
        }
        if metadata:
            meta.update({k: str(v) for k, v in metadata.items()})
        col.upsert(documents=[content], ids=[drawer_id], metadatas=[meta])
        return drawer_id

    # ------------------------------------------------------------------
    # Convenience store helpers
    # ------------------------------------------------------------------

    def store_exchange(
        self, question: str, response: str, *,
        room: str = "exploration", branch: str = "main",
    ) -> str:
        """Store a conversation exchange (Q+A pair)."""
        content = f"> {question}\n\n{response}"
        return self.store(content, room=room, branch=branch, source_type="exchange")

    def store_step(self, step: dict, *, branch: str = "main") -> str:
        """Store a provenance step as searchable text."""
        tool_id = step.get("tool_id", "unknown")
        room = room_from_tool(tool_id)
        lines = [f"Tool: {tool_id}"]
        params = step.get("parameters", {})
        if params:
            lines.append(f"Parameters: {_fmt(params)}")
        effects = step.get("effects", {})
        if effects:
            lines.append(f"Effects: {_fmt(effects)}")
        desc = step.get("description", "")
        if desc:
            lines.append(f"Description: {desc}")
        return self.store(
            "\n".join(lines), room=room, branch=branch,
            source_type="provenance", metadata={"tool_id": tool_id},
        )

    def store_decision(
        self, decision: str, rationale: str, *,
        room: str = "decisions", branch: str = "main",
    ) -> str:
        """Store an analysis decision with rationale."""
        content = f"Decision: {decision}\nRationale: {rationale}"
        return self.store(content, room=room, branch=branch, source_type="decision")

    # ------------------------------------------------------------------
    # Recall
    # ------------------------------------------------------------------

    def recall(
        self,
        query: str,
        *,
        room: str | None = None,
        branch: str | None = "main",
        n_results: int = 5,
    ) -> list[dict]:
        """Search memory.

        Parameters
        ----------
        query
            Natural language search.
        room
            Filter to one room (or ``None`` for all rooms).
        branch
            Filter to one branch.  Pass ``None`` to search across all
            branches (cross-branch communication).
        n_results
            Max results.

        Returns
        -------
        List of hits with ``text``, ``room``, ``branch``, ``similarity``,
        ``metadata``.
        """
        col = self._col()
        conditions: list[dict] = [{"wing": self.project_name}]
        if branch is not None:
            conditions.append({"branch": branch})
        if room is not None:
            conditions.append({"room": normalize_room(room)})
        where = {"$and": conditions} if len(conditions) > 1 else conditions[0]

        try:
            results = col.query(
                query_texts=[query], n_results=n_results,
                where=where,
                include=["documents", "metadatas", "distances"],
            )
        except Exception as e:
            logger.warning("Memory recall failed: %s", e)
            return []

        docs = results.get("documents", [[]])[0]
        metas = results.get("metadatas", [[]])[0]
        dists = results.get("distances", [[]])[0]
        return [
            {
                "text": doc,
                "room": m.get("room", "unknown"),
                "branch": m.get("branch", "main"),
                "source_type": m.get("source_type", "unknown"),
                "similarity": round(max(0.0, 1 - d), 3),
                "metadata": m,
            }
            for doc, m, d in zip(docs, metas, dists)
        ]

    # ------------------------------------------------------------------
    # Status
    # ------------------------------------------------------------------

    def status(self) -> dict:
        """Palace statistics."""
        col = self._col()
        total = col.count()
        try:
            r = col.get(where={"wing": self.project_name},
                        limit=10_000, include=["metadatas"])
            metas = r.get("metadatas", [])
        except Exception:
            metas = []

        rooms: dict[str, int] = {}
        branches: dict[str, int] = {}
        types: dict[str, int] = {}
        sessions: set[str] = set()
        for m in metas:
            rooms[m.get("room", "?")] = rooms.get(m.get("room", "?"), 0) + 1
            branches[m.get("branch", "?")] = branches.get(m.get("branch", "?"), 0) + 1
            types[m.get("source_type", "?")] = types.get(m.get("source_type", "?"), 0) + 1
            s = m.get("session_id")
            if s:
                sessions.add(s)
        return {
            "project": self.project_name,
            "total_drawers": total,
            "project_drawers": len(metas),
            "rooms": rooms,
            "branches": branches,
            "types": types,
            "sessions": len(sessions),
        }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_id(wing: str, room: str, content: str, branch: str = "main") -> str:
    raw = f"{wing}:{branch}:{room}:{content}"
    return f"drawer_{wing}_{room}_{hashlib.sha256(raw.encode()).hexdigest()[:24]}"


def _fmt(d: dict) -> str:
    return ", ".join(
        f"{k}={v}" if not isinstance(v, (dict, list)) else f"{k}=..."
        for k, v in d.items()
    )
