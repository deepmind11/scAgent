"""Provenance tracking for scAgent analysis steps.

Records every tool invocation as a W3C PROV-O compliant graph serialized
to JSON-LD.  The graph is append-only: activities and entities are never
modified after creation.

Usage::

    from pathlib import Path
    from scagent.provenance import ProvenanceGraph, record_step

    graph = ProvenanceGraph(Path(".scagent"))
    result = filter_cells(adata, min_genes=200, max_genes=5000)
    record_step(graph, result, user_prompt="filter with standard thresholds")
"""

from __future__ import annotations

import json
import platform
import uuid
from copy import deepcopy
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import scagent

# ---------------------------------------------------------------------------
# Internal data model (mirrors W3C PROV-O)
# ---------------------------------------------------------------------------

@dataclass
class Entity:
    """A data state — an AnnData snapshot at a point in the analysis."""

    id: str
    hash: str  # SHA-256 of h5ad; empty until Chunk 5
    branch: str
    step_index: int
    generated_by: str  # activity id
    derived_from: str | None  # previous entity id (None for first)


@dataclass
class Activity:
    """A single tool invocation."""

    id: str
    tool_id: str
    parameters: dict[str, Any]
    extras: dict[str, Any]  # tool-specific metadata (cells_before, n_clusters, …)
    used: str  # input entity id
    generated: str  # output entity id
    started_at: str  # ISO-8601
    ended_at: str  # ISO-8601
    user_prompt: str
    software_versions: dict[str, str]
    branch: str


@dataclass
class ToolAgent:
    """A software tool registered in the tool registry."""

    id: str
    tool_id: str
    registry_version: str = "0.1.0"


@dataclass
class Session:
    """A researcher interaction session — marks when scAgent was opened/closed."""

    id: str  # e.g., "sca:session_001"
    session_number: int
    started_at: str  # ISO-8601
    ended_at: str | None  # None while session is active
    software_versions: dict[str, str]
    activities: list[str] = field(default_factory=list)  # activity IDs recorded in this session


# ---------------------------------------------------------------------------
# Version capture helper
# ---------------------------------------------------------------------------

def _capture_versions() -> dict[str, str]:
    """Snapshot current software versions."""
    from importlib.metadata import version as pkg_version

    versions: dict[str, str] = {"python": platform.python_version(), "scagent": scagent.__version__}
    for pkg in ("scanpy", "anndata"):
        try:
            versions[pkg] = pkg_version(pkg)
        except Exception:
            pass
    return versions


# ---------------------------------------------------------------------------
# ProvenanceGraph
# ---------------------------------------------------------------------------

class ProvenanceGraph:
    """Append-only provenance graph persisted as PROV-JSONLD.

    Parameters
    ----------
    project_dir
        Path to the ``.scagent`` directory.  Created on first ``save()``
        if it does not exist.
    """

    FILENAME = "provenance.jsonld"

    def __init__(self, project_dir: Path | str) -> None:
        self._dir = Path(project_dir)
        self._path = self._dir / self.FILENAME

        self._entities: list[Entity] = []
        self._activities: list[Activity] = []
        self._agents: dict[str, ToolAgent] = {}  # keyed by tool_id
        self._sessions: list[Session] = []
        self._promoted_branch: str | None = None

        # Monotonic counter for unique activity IDs
        self._counter: int = 0

        # Per-branch step counters
        self._branch_steps: dict[str, int] = {}

        # Load existing graph if present
        if self._path.exists():
            self._load()

        # Start a new session
        self._start_session()

    # ------------------------------------------------------------------
    # Session management
    # ------------------------------------------------------------------

    def _start_session(self) -> None:
        """Start a new session.  Called automatically on __init__."""
        n = len(self._sessions) + 1
        self._current_session = Session(
            id=f"sca:session_{n:03d}",
            session_number=n,
            started_at=datetime.now(timezone.utc).isoformat(),
            ended_at=None,
            software_versions=_capture_versions(),
        )
        self._sessions.append(self._current_session)

    def end_session(self) -> None:
        """Mark the current session as ended.  Auto-saves."""
        if self._current_session and self._current_session.ended_at is None:
            self._current_session.ended_at = datetime.now(timezone.utc).isoformat()
            self.save()

    @property
    def current_session(self) -> Session:
        return self._current_session

    @property
    def sessions(self) -> list[dict[str, Any]]:
        """Return all sessions as dicts."""
        return [asdict(s) for s in self._sessions]

    # ------------------------------------------------------------------
    # Recording
    # ------------------------------------------------------------------

    def record(
        self,
        tool_id: str,
        parameters: dict[str, Any],
        extras: dict[str, Any] | None = None,
        *,
        input_hash: str = "",
        output_hash: str = "",
        user_prompt: str = "",
        branch: str = "main",
        started_at: str | None = None,
        ended_at: str | None = None,
    ) -> str:
        """Record a tool invocation.  Returns the activity ID.

        Automatically creates the input entity for the first step in a
        branch, links entities in chain order, registers the tool agent,
        and auto-saves.
        """
        now = datetime.now(timezone.utc).isoformat()
        if started_at is None:
            started_at = now
        if ended_at is None:
            ended_at = now

        self._counter += 1
        step = self._branch_steps.get(branch, 0)

        # Ensure a tool agent exists
        agent_id = f"sca:tool_{tool_id}"
        if tool_id not in self._agents:
            self._agents[tool_id] = ToolAgent(id=agent_id, tool_id=tool_id)

        # Resolve input entity
        if step == 0:
            # First step on this branch — create a synthetic input entity
            input_entity_id = f"sca:input_{branch}"
            if not self._find_entity(input_entity_id):
                self._entities.append(
                    Entity(
                        id=input_entity_id,
                        hash=input_hash,
                        branch=branch,
                        step_index=-1,
                        generated_by="",
                        derived_from=None,
                    )
                )
        else:
            # Link to previous entity on this branch
            input_entity_id = self._latest_entity_id(branch)

        # Create output entity
        # Derive a readable suffix from the tool_id
        step_name = tool_id.replace("_", "_")
        output_entity_id = f"sca:adata_{step_name}_v{step}"

        output_entity = Entity(
            id=output_entity_id,
            hash=output_hash,
            branch=branch,
            step_index=step,
            generated_by="",  # filled below
            derived_from=input_entity_id,
        )

        # Create activity
        activity_id = f"sca:activity_{tool_id}_{self._counter:03d}"
        output_entity.generated_by = activity_id

        activity = Activity(
            id=activity_id,
            tool_id=tool_id,
            parameters=deepcopy(parameters),
            extras=deepcopy(extras) if extras else {},
            used=input_entity_id,
            generated=output_entity_id,
            started_at=started_at,
            ended_at=ended_at,
            user_prompt=user_prompt,
            software_versions=_capture_versions(),
            branch=branch,
        )

        self._entities.append(output_entity)
        self._activities.append(activity)
        self._branch_steps[branch] = step + 1

        # Tag this activity to the current session
        if self._current_session:
            self._current_session.activities.append(activity_id)

        self.save()
        return activity_id

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def get_activity(self, activity_id: str) -> dict[str, Any] | None:
        """Return an activity by ID, or *None*."""
        for a in self._activities:
            if a.id == activity_id:
                return asdict(a)
        return None

    def get_entity(self, entity_id: str) -> dict[str, Any] | None:
        """Return an entity by ID, or *None*."""
        for e in self._entities:
            if e.id == entity_id:
                return asdict(e)
        return None

    def list_activities(
        self,
        branch: str | None = None,
        tool_id: str | None = None,
    ) -> list[dict[str, Any]]:
        """List activities, optionally filtered by branch and/or tool_id."""
        out: list[dict[str, Any]] = []
        for a in self._activities:
            if branch is not None and a.branch != branch:
                continue
            if tool_id is not None and a.tool_id != tool_id:
                continue
            out.append(
                {
                    "activity_id": a.id,
                    "tool_id": a.tool_id,
                    "parameters": a.parameters,
                    "extras": a.extras,
                    "started_at": a.started_at,
                    "ended_at": a.ended_at,
                    "branch": a.branch,
                    "input": a.used,
                    "output": a.generated,
                }
            )
        return out

    def get_chain(self, branch: str = "main") -> list[dict[str, Any]]:
        """Return the ordered activity chain for *branch* (this branch only)."""
        return [a for a in self.list_activities(branch=branch)]

    def get_full_chain(self, branch: str = "main") -> list[dict[str, Any]]:
        """Return the complete lineage for *branch*, including inherited parent steps.

        When a branch was forked from another, this traces back through the
        fork point and prepends the parent's steps.  The result is the full
        sequence of tool invocations needed to reproduce the branch's current
        state from raw data.
        """
        # Find the fork entity for this branch (if any)
        fork_entity = None
        for e in self._entities:
            if e.branch == branch and e.step_index == -1 and e.id.startswith("sca:fork_"):
                fork_entity = e
                break

        if fork_entity is None or fork_entity.derived_from is None:
            # Not forked — just return this branch's chain
            return self.get_chain(branch)

        # Find which branch + step the fork came from
        parent_entity = self._find_entity(fork_entity.derived_from)
        if parent_entity is None:
            return self.get_chain(branch)

        # Recursively get the parent's full chain up to the fork point
        parent_chain = self.get_full_chain(parent_entity.branch)
        # Keep parent steps up to and including the fork point
        prefix = [a for a in parent_chain
                  if a["output"] == parent_entity.id
                  or parent_chain.index(a) < len(parent_chain)]

        # The fork point is the parent entity's step_index
        # Keep all parent activities up through the entity that was forked from
        cutoff_ids: set[str] = set()
        for e in self._entities:
            if e.branch == parent_entity.branch and e.step_index <= parent_entity.step_index:
                cutoff_ids.add(e.generated_by)
        prefix = [a for a in parent_chain if a["activity_id"] in cutoff_ids]

        return prefix + self.get_chain(branch)

    def summary(self, branch: str = "main") -> str:
        """Human-readable Markdown table of the analysis provenance."""
        chain = self.get_chain(branch)
        if not chain:
            return f"_No provenance recorded for branch `{branch}`._"

        # Build a lookup: activity_id → session_number
        activity_to_session: dict[str, int] = {}
        for s in self._sessions:
            for aid in s.activities:
                activity_to_session[aid] = s.session_number

        lines = [
            f"## Analysis Provenance (branch: {branch}, {len(chain)} steps, {len(self._sessions)} sessions)\n",
            "| # | Tool | Key Params | Extras |",
            "|---|------|------------|--------|",
        ]
        current_session_num = None
        for i, a in enumerate(chain):
            # Insert session boundary marker
            sess_num = activity_to_session.get(a["activity_id"])
            if sess_num is not None and sess_num != current_session_num:
                sess = self._sessions[sess_num - 1]
                ts = sess.started_at[:19].replace("T", " ")  # trim to readable
                lines.append(f"| | **── Session {sess_num} ({ts}) ──** | | |")
                current_session_num = sess_num

            params_str = ", ".join(
                f"{k}={_fmt(v)}" for k, v in list(a["parameters"].items())[:4]
            )
            extras_str = ", ".join(
                f"{k}={_fmt(v)}" for k, v in a["extras"].items()
            )
            lines.append(f"| {i} | {a['tool_id']} | {params_str} | {extras_str} |")

        # Software versions from last activity
        if chain:
            last = self._activities[-1]
            vers = ", ".join(f"{k} {v}" for k, v in last.software_versions.items())
            lines.append(f"\n_Software: {vers}_")

        return "\n".join(lines)

    def diff(self, branch_a: str, branch_b: str) -> dict[str, Any]:
        """Compare two branches.  Returns divergence point and per-branch deltas."""
        chain_a = self.get_chain(branch_a)
        chain_b = self.get_chain(branch_b)

        # Find divergence point — last shared entity
        shared = 0
        for a, b in zip(chain_a, chain_b):
            if a["tool_id"] == b["tool_id"] and a["parameters"] == b["parameters"]:
                shared += 1
            else:
                break

        diverge_entity = chain_a[shared - 1]["input"] if shared > 0 else None

        return {
            "diverge_point": diverge_entity,
            "diverge_step": shared,
            "shared_steps": shared,
            "branch_a_only": chain_a[shared:],
            "branch_b_only": chain_b[shared:],
            "parameter_diffs": _find_param_diffs(chain_a[shared:], chain_b[shared:]),
        }

    def replay_plan(self, branch: str = "main", *, full: bool = True) -> list[tuple[str, dict[str, Any]]]:
        """Return ``[(tool_id, parameters), …]`` for reproducing the analysis.

        Parameters
        ----------
        branch
            Which branch to generate the plan for.
        full
            If *True* (default), includes inherited parent steps for forked
            branches — i.e. the complete pipeline from raw data.
            If *False*, only steps on this branch.
        """
        chain = self.get_full_chain(branch) if full else self.get_chain(branch)
        return [(a["tool_id"], a["parameters"]) for a in chain]

    def promote_branch(self, branch: str) -> None:
        """Mark *branch* as the canonical/final analysis pipeline.

        Creates a metadata marker so ``export_plan()`` and colleagues
        know which branch represents the settled result vs. exploration.
        """
        if branch not in self._branch_steps:
            raise ValueError(f"Branch '{branch}' does not exist")
        self._promoted_branch = branch
        self.save()

    @property
    def promoted_branch(self) -> str | None:
        """The branch marked as the canonical result, or *None*."""
        return self._promoted_branch

    def export_plan(self) -> list[tuple[str, dict[str, Any]]]:
        """Return the replay plan for the promoted (canonical) branch.

        If no branch has been promoted, falls back to ``main``.
        This is what you send to a colleague for reproduction.
        """
        branch = self._promoted_branch or "main"
        return self.replay_plan(branch, full=True)

    @property
    def n_activities(self) -> int:
        return len(self._activities)

    @property
    def n_entities(self) -> int:
        return len(self._entities)

    @property
    def branches(self) -> list[str]:
        return sorted(self._branch_steps.keys())

    # ------------------------------------------------------------------
    # Serialization (PROV-JSONLD)
    # ------------------------------------------------------------------

    def serialize(self) -> dict[str, Any]:
        """Return the full graph as a W3C PROV-O JSON-LD document."""
        graph_nodes: list[dict[str, Any]] = []

        # Entities
        for e in self._entities:
            node: dict[str, Any] = {
                "@id": e.id,
                "@type": "prov:Entity",
                "sca:snapshot_hash": e.hash,
                "sca:branch": e.branch,
                "sca:step_index": e.step_index,
            }
            if e.generated_by:
                node["prov:wasGeneratedBy"] = e.generated_by
            if e.derived_from:
                node["prov:wasDerivedFrom"] = e.derived_from
            graph_nodes.append(node)

        # Activities
        for a in self._activities:
            node = {
                "@id": a.id,
                "@type": "prov:Activity",
                "prov:used": a.used,
                "prov:wasAssociatedWith": f"sca:tool_{a.tool_id}",
                "prov:startedAtTime": a.started_at,
                "prov:endedAtTime": a.ended_at,
                "sca:parameters": a.parameters,
                "sca:triggered_by": "user_request",
                "sca:user_prompt": a.user_prompt,
                "sca:branch": a.branch,
            }
            if a.extras:
                node["sca:extras"] = a.extras
            for k, v in a.software_versions.items():
                node[f"sca:{k}_version"] = v
            graph_nodes.append(node)

        # Tool agents
        for ta in self._agents.values():
            graph_nodes.append(
                {
                    "@id": ta.id,
                    "@type": ["prov:Agent", "prov:SoftwareAgent"],
                    "sca:tool_id": ta.tool_id,
                    "sca:registry_version": ta.registry_version,
                }
            )

        # Sessions
        for s in self._sessions:
            snode: dict[str, Any] = {
                "@id": s.id,
                "@type": "sca:Session",
                "sca:session_number": s.session_number,
                "prov:startedAtTime": s.started_at,
                "sca:activities": s.activities,
            }
            if s.ended_at:
                snode["prov:endedAtTime"] = s.ended_at
            for k, v in s.software_versions.items():
                snode[f"sca:{k}_version"] = v
            graph_nodes.append(snode)

        doc: dict[str, Any] = {
            "@context": {
                "prov": "http://www.w3.org/ns/prov#",
                "sca": "https://scagent.dev/ontology/v1#",
                "xsd": "http://www.w3.org/2001/XMLSchema#",
            },
            "@graph": graph_nodes,
        }
        if self._promoted_branch:
            doc["sca:promoted_branch"] = self._promoted_branch
        return doc

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self) -> Path:
        """Write the graph to ``.scagent/provenance.jsonld``."""
        self._dir.mkdir(parents=True, exist_ok=True)
        data = self.serialize()
        self._path.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
        return self._path

    def _load(self) -> None:
        """Populate internal state from an existing JSONLD file."""
        data = json.loads(self._path.read_text(encoding="utf-8"))
        self._promoted_branch = data.get("sca:promoted_branch")
        graph_nodes = data.get("@graph", [])

        # Separate by type
        entity_nodes = []
        activity_nodes = []
        agent_nodes = []
        session_nodes = []

        for node in graph_nodes:
            ntype = node.get("@type", "")
            if ntype == "prov:Entity":
                entity_nodes.append(node)
            elif ntype == "prov:Activity":
                activity_nodes.append(node)
            elif ntype == "sca:Session":
                session_nodes.append(node)
            elif isinstance(ntype, list) and "prov:Agent" in ntype:
                agent_nodes.append(node)

        # Reconstruct agents
        for n in agent_nodes:
            ta = ToolAgent(
                id=n["@id"],
                tool_id=n.get("sca:tool_id", ""),
                registry_version=n.get("sca:registry_version", "0.1.0"),
            )
            self._agents[ta.tool_id] = ta

        # Reconstruct entities
        for n in entity_nodes:
            e = Entity(
                id=n["@id"],
                hash=n.get("sca:snapshot_hash", ""),
                branch=n.get("sca:branch", "main"),
                step_index=n.get("sca:step_index", -1),
                generated_by=n.get("prov:wasGeneratedBy", ""),
                derived_from=n.get("prov:wasDerivedFrom"),
            )
            self._entities.append(e)

        # Reconstruct activities
        for n in activity_nodes:
            # Extract software versions
            sw: dict[str, str] = {}
            for k, v in n.items():
                if k.startswith("sca:") and k.endswith("_version"):
                    sw_key = k[4:-8]  # strip "sca:" and "_version"
                    sw[sw_key] = v

            a = Activity(
                id=n["@id"],
                tool_id=n.get("prov:wasAssociatedWith", "").replace("sca:tool_", ""),
                parameters=n.get("sca:parameters", {}),
                extras=n.get("sca:extras", {}),
                used=n.get("prov:used", ""),
                generated="",  # resolved below
                started_at=n.get("prov:startedAtTime", ""),
                ended_at=n.get("prov:endedAtTime", ""),
                user_prompt=n.get("sca:user_prompt", ""),
                software_versions=sw,
                branch=n.get("sca:branch", "main"),
            )
            # Find generated entity
            for e in self._entities:
                if e.generated_by == a.id:
                    a.generated = e.id
                    break
            self._activities.append(a)

        # Reconstruct sessions
        for n in sorted(session_nodes, key=lambda x: x.get("sca:session_number", 0)):
            sw: dict[str, str] = {}
            for k, v in n.items():
                if k.startswith("sca:") and k.endswith("_version"):
                    sw[k[4:-8]] = v
            s = Session(
                id=n["@id"],
                session_number=n.get("sca:session_number", 0),
                started_at=n.get("prov:startedAtTime", ""),
                ended_at=n.get("prov:endedAtTime"),
                software_versions=sw,
                activities=n.get("sca:activities", []),
            )
            self._sessions.append(s)

        # Restore counters
        self._counter = len(self._activities)
        for a in self._activities:
            step_idx = -1
            for e in self._entities:
                if e.generated_by == a.id:
                    step_idx = e.step_index
                    break
            if step_idx >= 0:
                branch = a.branch
                self._branch_steps[branch] = max(
                    self._branch_steps.get(branch, 0), step_idx + 1
                )

    @classmethod
    def load(cls, path: Path | str) -> "ProvenanceGraph":
        """Load a provenance graph from a directory containing ``provenance.jsonld``."""
        p = Path(path)
        if p.name == cls.FILENAME:
            p = p.parent
        return cls(p)

    # ------------------------------------------------------------------
    # Branch support (minimal — full management is Chunk 5)
    # ------------------------------------------------------------------

    def fork_branch(self, new_branch: str, from_branch: str = "main") -> None:
        """Create a new branch forked from the current head of *from_branch*.

        The new branch will continue from the last entity of *from_branch*.
        """
        if from_branch not in self._branch_steps:
            raise ValueError(f"Branch '{from_branch}' does not exist")
        if new_branch in self._branch_steps:
            raise ValueError(f"Branch '{new_branch}' already exists")

        # The next record on new_branch will derive from from_branch's latest entity
        head_entity_id = self._latest_entity_id(from_branch)

        # Create a synthetic link entity on the new branch
        fork_entity = Entity(
            id=f"sca:fork_{new_branch}_from_{from_branch}",
            hash="",
            branch=new_branch,
            step_index=-1,
            generated_by="",
            derived_from=head_entity_id,
        )
        self._entities.append(fork_entity)
        self._branch_steps[new_branch] = 0
        self.save()

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _find_entity(self, entity_id: str) -> Entity | None:
        for e in self._entities:
            if e.id == entity_id:
                return e
        return None

    def _latest_entity_id(self, branch: str) -> str:
        """Return the ID of the most recent entity on *branch*."""
        candidates = [e for e in self._entities if e.branch == branch]
        if not candidates:
            raise ValueError(f"No entities on branch '{branch}'")
        return max(candidates, key=lambda e: e.step_index).id


def _fmt(v: Any) -> str:
    """Format a value for the summary table."""
    if isinstance(v, float):
        return f"{v:g}"
    if isinstance(v, list) and len(v) > 3:
        return f"[{v[0]}, …, {v[-1]}] ({len(v)} items)"
    return str(v)


def _find_param_diffs(
    chain_a: list[dict[str, Any]],
    chain_b: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Find same-tool invocations with different parameters across two chains."""
    diffs: list[dict[str, Any]] = []
    tools_a = {a["tool_id"]: a["parameters"] for a in chain_a}
    tools_b = {b["tool_id"]: b["parameters"] for b in chain_b}

    for tool in set(tools_a) & set(tools_b):
        pa, pb = tools_a[tool], tools_b[tool]
        for k in set(pa) | set(pb):
            va, vb = pa.get(k), pb.get(k)
            if va != vb:
                diffs.append({"tool_id": tool, "param": k, "a": va, "b": vb})
    return diffs


# ---------------------------------------------------------------------------
# Convenience adapter — bridges tool result dicts to the graph
# ---------------------------------------------------------------------------

def record_step(
    graph: ProvenanceGraph,
    tool_result: dict[str, Any],
    *,
    input_hash: str = "",
    output_hash: str = "",
    user_prompt: str = "",
    branch: str = "main",
    started_at: str | None = None,
    ended_at: str | None = None,
) -> str:
    """Extract provenance from a tool result and record it.

    Accepts either a plain result dict (most tools) or an
    ``(adata, result_dict)`` tuple (load_10x_h5, filter_cells,
    filter_genes, detect_doublets).  Pulls ``tool_id`` and
    ``parameters`` from the ``"provenance"`` key, separates
    tool-specific extras, and forwards to ``graph.record()``.

    Returns the activity ID.
    """
    # Handle (adata, result_dict) tuples
    if isinstance(tool_result, tuple):
        tool_result = tool_result[-1]  # last element is the result dict

    prov = tool_result.get("provenance")
    if prov is None:
        raise ValueError("tool_result has no 'provenance' key")

    tool_id = prov["tool_id"]
    parameters = prov.get("parameters", {})

    # Everything in the provenance dict that is NOT tool_id or parameters
    # is an extra (cells_before, n_clusters, file_hash, etc.)
    extras = {k: v for k, v in prov.items() if k not in ("tool_id", "parameters")}

    return graph.record(
        tool_id=tool_id,
        parameters=parameters,
        extras=extras,
        input_hash=input_hash,
        output_hash=output_hash,
        user_prompt=user_prompt,
        branch=branch,
        started_at=started_at,
        ended_at=ended_at,
    )


def record_custom(
    graph: ProvenanceGraph,
    *,
    description: str,
    code: str,
    user_prompt: str = "",
    effects: dict[str, Any] | None = None,
    input_hash: str = "",
    output_hash: str = "",
    branch: str = "main",
    started_at: str | None = None,
    ended_at: str | None = None,
) -> str:
    """Record a custom (non-tool) analysis step in the provenance graph.

    Use this when the agent executes code that doesn't come from a
    registered ``scagent.tools.*`` function — e.g., ad-hoc computations,
    custom filtering, external package calls, or researcher-provided
    scripts.

    Parameters
    ----------
    description
        Short human-readable label (e.g., ``"mito/ribo ratio column"``).
    code
        The Python code that was executed, verbatim.
    user_prompt
        The researcher's original request.
    effects
        Optional dict describing what changed — e.g.,
        ``{"added_obs_columns": ["mito_ribo_ratio"], "cells_removed": 0}``.
        Free-form; stored as extras in provenance.
    """
    extras = {"code": code}
    if effects:
        extras["effects"] = effects

    return graph.record(
        tool_id="custom",
        parameters={"description": description},
        extras=extras,
        input_hash=input_hash,
        output_hash=output_hash,
        user_prompt=user_prompt,
        branch=branch,
        started_at=started_at,
        ended_at=ended_at,
    )
