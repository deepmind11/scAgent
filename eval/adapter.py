"""scAgent direct adapter for scBench.

Routes eval tasks to the appropriate handler function.
Each handler uses our scagent.tools wrappers with principled defaults.
"""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad

from .handlers.qc import handle_qc
from .handlers.normalization import handle_normalization
from .handlers.hvg import handle_hvg
from .handlers.celltyping import handle_celltyping
from .handlers.clustering import handle_clustering
from .handlers.de import handle_de
from .handlers.trajectory import handle_trajectory


# Map eval IDs to handler functions
EVAL_HANDLERS = {
    "chromium_qc_4T1_filter_cells": handle_qc,
    "chromium_4t1_normalization": handle_normalization,
    "chromium_4t1_hvg_gene_sets": handle_hvg,
    "chromium_celltyping_01_4t1_compartment_fractions": handle_celltyping,
    "chromium_clustering_01_4t1_pericyte_adjacent_to_caf": handle_clustering,
    "chromium_differential_expression_01_contractile_caf_marker_recovery": handle_de,
    "chromium_trajectory_01_caf_terminal_marker_recovery": handle_trajectory,
}


def scagent_direct_agent(task_prompt: str, work_dir: Path) -> dict:
    """Direct Python adapter: no LLM, calls tools directly.

    This is used as the ``agent_function`` for scBench's ``EvalRunner``.

    Parameters
    ----------
    task_prompt
        The eval task prompt (with contextual data appended).
    work_dir
        Working directory containing the symlinked .h5ad data files.

    Returns
    -------
    dict with ``"answer"`` key containing the eval answer.
    """
    # Find the .h5ad file(s) in work_dir (may be in data/ subdirectory)
    h5ad_files = list(Path(work_dir).glob("*.h5ad"))
    if not h5ad_files:
        h5ad_files = list(Path(work_dir).glob("**/*.h5ad"))
    if not h5ad_files:
        h5ad_files = list(Path(work_dir).glob("*.node"))
    if not h5ad_files:
        h5ad_files = list(Path(work_dir).glob("**/*.node"))
    if not h5ad_files:
        raise FileNotFoundError(f"No data files found in {work_dir}")

    # Load the first data file
    adata = ad.read_h5ad(h5ad_files[0])

    # Determine which eval this is by matching the task prompt
    eval_id = _identify_eval(task_prompt)
    if eval_id is None:
        raise ValueError("Could not identify eval type from task prompt")

    handler = EVAL_HANDLERS.get(eval_id)
    if handler is None:
        raise ValueError(f"No handler for eval: {eval_id}")

    print(f"  Handler: {eval_id}")
    print(f"  Data: {adata.n_obs} cells × {adata.n_vars} genes")

    # Run the handler
    answer = handler(adata, task_prompt)

    # Write eval_answer.json (scBench also checks for this file)
    answer_file = Path(work_dir) / "eval_answer.json"
    answer_file.write_text(json.dumps(answer, indent=2))

    return {"answer": answer}


def _identify_eval(task_prompt: str) -> str | None:
    """Identify the eval type from keywords in the task prompt."""
    prompt_lower = task_prompt.lower()

    if "cells_after_filtering" in prompt_lower:
        return "chromium_qc_4T1_filter_cells"
    if "gene_value" in prompt_lower and "mrc1" in prompt_lower:
        return "chromium_4t1_normalization"
    if "highly variable genes" in prompt_lower and "top_marker_genes" in prompt_lower:
        return "chromium_4t1_hvg_gene_sets"
    if "cell_type_distribution" in prompt_lower and "immune" in prompt_lower:
        return "chromium_celltyping_01_4t1_compartment_fractions"
    if "pericyte" in prompt_lower and "adjacency" in prompt_lower:
        return "chromium_clustering_01_4t1_pericyte_adjacent_to_caf"
    if "contractile" in prompt_lower and "caf" in prompt_lower:
        return "chromium_differential_expression_01_contractile_caf_marker_recovery"
    if "trajectory" in prompt_lower and "terminal" in prompt_lower:
        return "chromium_trajectory_01_caf_terminal_marker_recovery"

    return None
