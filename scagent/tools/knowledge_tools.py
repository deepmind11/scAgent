"""Tool wrappers for marker gene knowledge queries.

These are the callable functions that the agent invokes via tool registry.
"""

from __future__ import annotations

from scagent.knowledge import MarkerDB, MarkerHit, ValidationResult


def query_markers(
    cell_type: str,
    *,
    species: str = "human",
    tissue: str = "",
    top_n: int = 15,
    sources: list[str] | None = None,
    knowledge_dir: str | None = None,
    celltypist_models: list[str] | None = None,
) -> dict:
    """Query marker genes for a cell type.

    Parameters
    ----------
    cell_type
        Cell type name (e.g. ``"NK cells"``, ``"CD8+ T cells"``).
    species
        ``"human"`` or ``"mouse"``.
    tissue
        Optional tissue filter.
    top_n
        Max markers per source.
    sources
        Restrict to ``["canonical"]``, ``["celltypist"]``, etc.
    knowledge_dir
        Path to directory with CellMarker 2.0 / PanglaoDB files.
    celltypist_models
        CellTypist model names to use.

    Returns
    -------
    dict with ``cell_type``, ``markers`` list, ``sources`` set.
    """
    db = MarkerDB(knowledge_dir=knowledge_dir, celltypist_models=celltypist_models)
    hits = db.query(cell_type, species=species, tissue=tissue,
                    top_n=top_n, sources=sources)
    return _format_query(cell_type, hits)


def validate_annotation(
    cluster_markers: list[str],
    label: str,
    *,
    species: str = "human",
    top_n: int = 20,
    knowledge_dir: str | None = None,
    celltypist_models: list[str] | None = None,
) -> dict:
    """Validate a cluster's cell type label against known markers.

    Parameters
    ----------
    cluster_markers
        Top DE genes for the cluster.
    label
        Proposed cell type label.
    species
        ``"human"`` or ``"mouse"``.
    top_n
        How many known markers to compare.
    knowledge_dir
        Path to external marker databases.
    celltypist_models
        CellTypist model names to use.

    Returns
    -------
    dict with ``label``, ``confidence``, ``overlap_ratio``,
    ``matched_markers``, ``evidence``, ``alternatives``.
    """
    db = MarkerDB(knowledge_dir=knowledge_dir, celltypist_models=celltypist_models)
    result = db.validate_annotation(
        cluster_markers, label, species=species, top_n=top_n,
    )
    return _format_validation(result)


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _format_query(cell_type: str, hits: list[MarkerHit]) -> dict:
    by_source: dict[str, list[dict]] = {}
    for h in hits:
        src = h.source
        by_source.setdefault(src, []).append({
            "gene": h.gene, "rank": h.rank,
        })
    return {
        "cell_type": cell_type,
        "markers": [{"gene": h.gene, "source": h.source, "rank": h.rank} for h in hits],
        "sources": sorted(by_source.keys()),
        "by_source": by_source,
        "n_total": len(hits),
    }


def _format_validation(result: ValidationResult) -> dict:
    return {
        "label": result.proposed_label,
        "confidence": result.confidence,
        "overlap_ratio": result.overlap_ratio,
        "matched_markers": result.matched_markers,
        "total_known": result.total_known,
        "evidence": result.evidence,
        "alternatives": result.alternative_labels,
    }
