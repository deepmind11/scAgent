"""Composable preprocessing pipeline."""

from __future__ import annotations

from pathlib import Path

import anndata as ad

from .loading import load_10x_h5
from .qc import calculate_qc_metrics, filter_cells, filter_genes
from .doublets import detect_doublets
from .normalize import log_normalize


# Ordered list of step names for checkpoint resolution
PREPROCESSING_STEPS = [
    "load",
    "qc_metrics",
    "filter_cells",
    "filter_genes",
    "doublets",
    "normalize",
]


def run_preprocessing(
    filename: str,
    *,
    # QC params
    min_genes: int = 200,
    max_genes: int = 5000,
    max_pct_mito: float = 20.0,
    min_counts: int | None = None,
    min_cells: int = 3,
    # Doublet params
    expected_doublet_rate: float | None = None,
    n_prin_comps: int = 30,
    random_state: int = 0,
    remove_doublets: bool = True,
    # Normalize params
    target_sum: float = 1e4,
    # Control
    checkpoint_dir: str | None = None,
    plot_dir: str | None = None,
    stop_after: str | None = None,
    start_from: str | None = None,
) -> tuple[ad.AnnData, dict]:
    """Run the full preprocessing pipeline: load → QC → doublets → normalize.

    Parameters
    ----------
    filename
        Path to the 10x .h5 file. Ignored if ``start_from`` is set.
    checkpoint_dir
        If provided, save ``after_{step}.h5ad`` at each step. Also used
        to load checkpoints when ``start_from`` is set.
    plot_dir
        Directory for diagnostic plots.
    stop_after
        Stop after this step. One of: ``'load'``, ``'qc_metrics'``,
        ``'filter_cells'``, ``'filter_genes'``, ``'doublets'``, ``'normalize'``.
    start_from
        Resume from this checkpoint (e.g., ``'after_filter_cells'``).
        Requires ``checkpoint_dir`` to be set and the checkpoint file to exist.

    Returns
    -------
    (adata, results) where results is a dict mapping step names to their
    individual result dicts.
    """
    results: dict[str, dict] = {}

    # Determine which steps to run
    if start_from is not None:
        if checkpoint_dir is None:
            raise ValueError("checkpoint_dir must be set when using start_from.")
        adata = _load_checkpoint(checkpoint_dir, start_from)
        start_idx = _step_index(start_from.replace("after_", "")) + 1
    else:
        start_idx = 0
        adata = None  # will be loaded in step 0

    stop_idx = _step_index(stop_after) if stop_after else len(PREPROCESSING_STEPS) - 1

    # Execute steps
    for i in range(start_idx, stop_idx + 1):
        step = PREPROCESSING_STEPS[i]

        if step == "load":
            adata, result = load_10x_h5(
                filename, checkpoint_dir=checkpoint_dir,
            )
            results["load"] = result

        elif step == "qc_metrics":
            result = calculate_qc_metrics(adata, plot_dir=plot_dir)
            if checkpoint_dir:
                _save_checkpoint(adata, checkpoint_dir, "after_qc_metrics")
            results["qc_metrics"] = result

        elif step == "filter_cells":
            adata, result = filter_cells(
                adata,
                min_genes=min_genes,
                max_genes=max_genes,
                max_pct_mito=max_pct_mito,
                min_counts=min_counts,
                checkpoint_dir=checkpoint_dir,
            )
            results["filter_cells"] = result

        elif step == "filter_genes":
            adata, result = filter_genes(
                adata, min_cells=min_cells, checkpoint_dir=checkpoint_dir,
            )
            results["filter_genes"] = result

        elif step == "doublets":
            adata, result = detect_doublets(
                adata,
                expected_doublet_rate=expected_doublet_rate,
                n_prin_comps=n_prin_comps,
                random_state=random_state,
                remove=remove_doublets,
                plot_dir=plot_dir,
                checkpoint_dir=checkpoint_dir,
            )
            results["doublets"] = result

        elif step == "normalize":
            result = log_normalize(
                adata,
                target_sum=target_sum,
                checkpoint_dir=checkpoint_dir,
            )
            results["normalize"] = result

    return adata, results


def _step_index(step_name: str) -> int:
    """Get the index of a step in the pipeline."""
    try:
        return PREPROCESSING_STEPS.index(step_name)
    except ValueError:
        raise ValueError(
            f"Unknown step '{step_name}'. Valid steps: {PREPROCESSING_STEPS}"
        )


def _load_checkpoint(checkpoint_dir: str, name: str) -> ad.AnnData:
    """Load an AnnData checkpoint."""
    path = Path(checkpoint_dir) / f"{name}.h5ad"
    if not path.exists():
        raise FileNotFoundError(
            f"Checkpoint not found: {path}. "
            f"Available checkpoints: {_list_checkpoints(checkpoint_dir)}"
        )
    return ad.read_h5ad(path)


def _list_checkpoints(checkpoint_dir: str) -> list[str]:
    """List available checkpoint names."""
    d = Path(checkpoint_dir)
    if not d.exists():
        return []
    return sorted(p.stem for p in d.glob("after_*.h5ad"))


def _save_checkpoint(adata: ad.AnnData, checkpoint_dir: str, name: str) -> str:
    d = Path(checkpoint_dir)
    d.mkdir(parents=True, exist_ok=True)
    path = d / f"{name}.h5ad"
    adata.write(path)
    return str(path)
