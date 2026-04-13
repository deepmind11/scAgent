"""scAgent — Agentic single-cell RNA-seq analysis system."""

__version__ = "0.1.0"

from scagent.provenance import ProvenanceGraph, record_step, record_custom
from scagent.state import StateManager
from scagent.context import ExperimentContext
from scagent.dag import AnalysisDAG
from scagent.memory import ProjectMemory
from scagent.knowledge import MarkerDB
from scagent.export import generate_methods, generate_repro_package
from scagent.inspector import AnnDataState, inspect_adata, find_raw_counts, summarize_state
from scagent.dependencies import check_prerequisites, plan_steps, ensure_ready_for

__all__ = [
    "ProvenanceGraph", "record_step", "record_custom",
    "StateManager",
    "ExperimentContext", "AnalysisDAG",
    "ProjectMemory",
    "MarkerDB",
    "generate_methods", "generate_repro_package",
    "AnnDataState", "inspect_adata", "find_raw_counts", "summarize_state",
    "check_prerequisites", "plan_steps", "ensure_ready_for",
]
