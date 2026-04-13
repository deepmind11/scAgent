"""Static dependency graph for scRNA-seq analysis steps.

Encodes what each analysis step requires — both prerequisite steps
and data conditions.  Used to validate whether a step can run given
the current AnnData state, and to plan the minimal set of steps
needed to reach a goal.

Usage::

    from scagent.inspector import inspect_adata
    from scagent.dependencies import check_prerequisites, plan_steps, ensure_ready_for

    state = inspect_adata(adata)
    can_run, missing = check_prerequisites("clustering", state)
    steps = plan_steps("pseudobulk_de", state)
    adata = ensure_ready_for(adata, state, needs="log_normalized")
"""

from __future__ import annotations

import logging
from copy import copy
from typing import Any

import numpy as np
from scipy.sparse import issparse

from scagent.inspector import (
    AnnDataState,
    RAW_COUNTS,
    LOG_NORMALIZED,
    SCALED,
    TRANSFORMED_UNKNOWN,
    inspect_adata,
    find_raw_counts,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Dependency definitions
# ---------------------------------------------------------------------------

# Maps step_id → requirements.
#
# "requires"       : list of step_ids that must be completed first
# "requires_x"     : list of acceptable x_state values
# "requires_data"  : list of data conditions that must be met:
#                     "raw_counts"    → raw counts must be available somewhere
#                     "condition_key" → a condition column must exist in obs
#                     "batch_key"     → a batch column must exist in obs
# "provides"       : what state flag this step sets to True
#                     (used by plan_steps to know when a dep is satisfied)

DEPENDENCIES: dict[str, dict[str, Any]] = {
    "qc_metrics": {
        "requires": [],
        "requires_x": [RAW_COUNTS],
        "provides": "has_qc_metrics",
    },
    "filter_cells": {
        "requires": ["qc_metrics"],
        "provides": "has_qc_metrics",  # still has QC after filtering
    },
    "filter_genes": {
        "requires": ["filter_cells"],
        "provides": "has_qc_metrics",
    },
    "doublet_detection": {
        "requires": ["filter_genes"],
        "provides": "has_qc_metrics",
    },
    "normalize": {
        "requires": [],
        "requires_x": [RAW_COUNTS],
        "provides": "has_normalized",
    },
    "hvg": {
        "requires": ["normalize"],
        "requires_x": [LOG_NORMALIZED],
        "provides": "has_hvg",
    },
    "pca": {
        "requires": ["hvg"],
        "provides": "has_pca",
    },
    "neighbors": {
        "requires": ["pca"],
        "provides": "has_neighbors",
    },
    "clustering": {
        "requires": ["neighbors"],
        "provides": "has_clusters",
    },
    "umap": {
        "requires": ["neighbors"],
        "provides": "has_umap",
    },
    "tsne": {
        "requires": ["neighbors"],
        "provides": "has_tsne",
    },
    "marker_genes": {
        "requires": ["clustering"],
        "provides": "has_de_results",
    },
    "annotation": {
        "requires": ["clustering"],
        "provides": "has_cell_types",
    },
    "pseudobulk_de": {
        "requires": ["annotation"],
        "requires_data": ["raw_counts", "condition_key"],
        "provides": "has_de_results",
    },
    "trajectory": {
        "requires": ["neighbors"],
        "provides": "has_diffmap",
    },
    "gsea": {
        "requires": ["pseudobulk_de"],
        "provides": None,
    },
    "composition": {
        "requires": ["annotation"],
        "requires_data": ["condition_key"],
        "provides": None,
    },
}

# Map provides → state field for checking if a step is already done
_PROVIDES_TO_STATE: dict[str, str] = {
    v["provides"]: k
    for k, v in DEPENDENCIES.items()
    if v.get("provides")
}


# ---------------------------------------------------------------------------
# Check prerequisites
# ---------------------------------------------------------------------------


def check_prerequisites(
    goal: str, state: AnnDataState
) -> tuple[bool, list[str]]:
    """Check if a goal step can be run given the current state.

    Parameters
    ----------
    goal
        Step ID (e.g., ``"clustering"``, ``"pseudobulk_de"``).
    state
        Inspector output.

    Returns
    -------
    ``(can_run, missing)`` where *missing* lists human-readable reasons
    why the step cannot run yet.
    """
    if goal not in DEPENDENCIES:
        return True, []  # unknown step — assume OK

    dep = DEPENDENCIES[goal]
    missing: list[str] = []

    # Check required predecessor steps
    for req_step in dep.get("requires", []):
        if not _step_is_done(req_step, state):
            req_dep = DEPENDENCIES.get(req_step, {})
            provides = req_dep.get("provides", req_step)
            missing.append(f"{req_step} not completed (need {provides})")

    # Check X state requirements
    req_x = dep.get("requires_x", [])
    if req_x and state.x_state not in req_x:
        missing.append(
            f"X must be {' or '.join(req_x)}, but is '{state.x_state}'"
        )

    # Check data conditions
    for data_req in dep.get("requires_data", []):
        if data_req == "raw_counts" and state.raw_counts_location is None:
            missing.append("raw counts not available (needed for pseudobulk)")
        elif data_req == "condition_key" and state.condition_key is None:
            missing.append("no condition column found in obs")
        elif data_req == "batch_key" and state.batch_key is None:
            missing.append("no batch column found in obs")

    return len(missing) == 0, missing


# ---------------------------------------------------------------------------
# Plan steps
# ---------------------------------------------------------------------------


def plan_steps(goal: str, state: AnnDataState) -> list[str]:
    """Work backwards from a goal, skip completed steps, return what's needed.

    Parameters
    ----------
    goal
        Target step ID.
    state
        Inspector output.

    Returns
    -------
    Ordered list of step IDs to execute (may be empty if goal is already
    done or has no missing prerequisites).
    """
    if goal not in DEPENDENCIES:
        return [goal]

    # If the goal is already done, nothing to do
    if _step_is_done(goal, state):
        return []

    # BFS backwards through dependencies
    needed: list[str] = []
    visited: set[str] = set()

    def _resolve(step_id: str) -> None:
        if step_id in visited:
            return
        visited.add(step_id)

        if _step_is_done(step_id, state):
            return  # already done, skip

        dep = DEPENDENCIES.get(step_id, {})
        for req in dep.get("requires", []):
            _resolve(req)

        needed.append(step_id)

    _resolve(goal)
    return needed


# ---------------------------------------------------------------------------
# ensure_ready_for — convenience function
# ---------------------------------------------------------------------------


def ensure_ready_for(
    adata,
    state: AnnDataState | None = None,
    *,
    needs: str = "log_normalized",
) -> "ad.AnnData":
    """Ensure adata.X is in the required state, modifying in place.

    Parameters
    ----------
    adata
        The AnnData object.
    state
        Pre-computed inspector state. Computed if *None*.
    needs
        Required X state: ``"raw_counts"`` or ``"log_normalized"``.

    Returns
    -------
    The same AnnData (modified in place) with X in the requested state.

    Raises
    ------
    ValueError
        If the data cannot be brought to the required state.
    """
    import scanpy as sc

    if state is None:
        state = inspect_adata(adata)

    if needs == "raw_counts":
        if state.x_state == RAW_COUNTS:
            return adata
        # Try to recover raw counts
        raw = find_raw_counts(adata, state)
        adata.X = raw
        return adata

    if needs == "log_normalized":
        if state.x_state == LOG_NORMALIZED:
            return adata

        if state.x_state == RAW_COUNTS:
            # Preserve raw counts before transforming
            if "counts" not in adata.layers:
                adata.layers["counts"] = adata.X.copy()
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            return adata

        if state.x_state in (SCALED, TRANSFORMED_UNKNOWN):
            # Try to recover raw counts from layers
            raw = find_raw_counts(adata, state)
            adata.X = raw
            if "counts" not in adata.layers:
                adata.layers["counts"] = raw.copy()
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            return adata

    raise ValueError(
        f"Cannot bring X to '{needs}' state. "
        f"Current state: '{state.x_state}', "
        f"raw counts: {state.raw_counts_location}"
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------


def _step_is_done(step_id: str, state: AnnDataState) -> bool:
    """Check if a step has already been completed based on inspector state."""
    dep = DEPENDENCIES.get(step_id)
    if dep is None:
        return False

    provides = dep.get("provides")
    if provides is None:
        return False

    return getattr(state, provides, False)
