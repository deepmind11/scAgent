"""Cell typing eval: annotate Immune / Epithelial-Cancer / CAF compartments."""

from __future__ import annotations
import numpy as np
import anndata as ad
import scanpy as sc


def handle_celltyping(adata: ad.AnnData, task_prompt: str) -> dict:
    """Classify cells into three compartments using canonical markers.

    Ground truth: Immune ~66.4%, Epithelial/Cancer ~24.5%, CAF ~8.0%
    """
    # Normalize if needed
    if adata.X.max() > 50:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    elif "raw_counts" in adata.layers:
        # Already normalized but raw counts available
        pass

    # Define compartment markers (mouse gene casing)
    # Case-insensitive matching as specified in the task
    var_names_lower = [g.lower() for g in adata.var_names]
    var_name_map = {g.lower(): g for g in adata.var_names}

    compartment_markers = {
        "Immune": ["ptprc", "cd3d", "cd3e", "cd79a", "cd14", "lyz2",
                    "nkg7", "s100a8", "s100a9", "csf1r", "itgam"],
        "Epithelial/Cancer": ["epcam", "krt8", "krt18", "krt19", "cdh1"],
        "CAF": ["pdgfra", "thy1", "col1a1", "col1a2", "dcn", "pdpn",
                "acta2", "fap"],
    }

    # Score each cell for each compartment
    scores = {}
    for compartment, markers in compartment_markers.items():
        # Find which markers exist in the data
        valid_markers = [var_name_map[m] for m in markers if m in var_name_map]
        if not valid_markers:
            scores[compartment] = np.zeros(adata.n_obs)
            continue
        sc.tl.score_genes(adata, gene_list=valid_markers,
                          score_name=f"_score_{compartment}")
        scores[compartment] = adata.obs[f"_score_{compartment}"].values.copy()

    # Assign each cell to highest-scoring compartment
    score_matrix = np.column_stack([scores[c] for c in compartment_markers])
    compartment_names = list(compartment_markers.keys())
    assignments = [compartment_names[i] for i in np.argmax(score_matrix, axis=1)]

    # Compute percentages
    from collections import Counter
    counts = Counter(assignments)
    total = sum(counts.values())
    distribution = {
        comp: round(counts.get(comp, 0) / total * 100, 1)
        for comp in compartment_names
    }

    return {"cell_type_distribution": distribution}
