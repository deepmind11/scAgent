"""HVG eval: select 2000 highly variable genes."""

from __future__ import annotations
import anndata as ad
from scagent.inspector import inspect_adata
from scagent.dependencies import ensure_ready_for
from scagent.tools.feature_selection import select_hvg


def handle_hvg(adata: ad.AnnData, task_prompt: str) -> dict:
    """Select 2000 HVGs and return the full gene list.

    The eval checks whether canonical biological markers are recovered
    in the HVG set (precision=0, so extra genes are fine).
    """
    state = inspect_adata(adata)
    ensure_ready_for(adata, state, needs="log_normalized")

    # Select HVGs
    select_hvg(adata, n_top_genes=2000)

    hvg_names = adata.var_names[adata.var["highly_variable"]].tolist()

    return {"top_marker_genes": hvg_names}
