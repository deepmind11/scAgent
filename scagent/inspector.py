"""AnnData state inspector for scAgent.

Inspects any AnnData object to determine what preprocessing has been done,
what's available to work with, and what state ``adata.X`` is in.  Runs once
when foreign data is loaded — the agent uses the result to reason about
what steps are needed for any given question.

Usage::

    from scagent.inspector import inspect_adata, find_raw_counts, summarize_state

    state = inspect_adata(adata)
    print(summarize_state(state))
    raw = find_raw_counts(adata, state)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import issparse

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# X-state constants
# ---------------------------------------------------------------------------

RAW_COUNTS = "raw_counts"
LOG_NORMALIZED = "log_normalized"
SCALED = "scaled"
TRANSFORMED_UNKNOWN = "transformed_unknown"

# ---------------------------------------------------------------------------
# AnnDataState
# ---------------------------------------------------------------------------


@dataclass
class AnnDataState:
    """Comprehensive state report for an AnnData object."""

    # Matrix state
    x_state: str  # RAW_COUNTS | LOG_NORMALIZED | SCALED | TRANSFORMED_UNKNOWN
    raw_counts_location: str | None  # "layers['counts']", "X", or None

    # Completed steps (booleans)
    has_qc_metrics: bool = False
    has_normalized: bool = False
    has_raw_frozen: bool = False
    has_hvg: bool = False
    has_pca: bool = False
    has_neighbors: bool = False
    has_umap: bool = False
    has_tsne: bool = False
    has_clusters: bool = False
    has_cell_types: bool = False
    has_de_results: bool = False
    has_diffmap: bool = False
    has_velocity_layers: bool = False

    # Detected keys (actual column/key names found)
    cluster_key: str | None = None
    celltype_key: str | None = None
    condition_key: str | None = None
    batch_key: str | None = None

    # Data shape
    n_cells: int = 0
    n_genes: int = 0
    species: str | None = None

    # Warnings
    warnings: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def inspect_adata(adata) -> AnnDataState:
    """Inspect an AnnData object and return a comprehensive state report.

    Parameters
    ----------
    adata
        An :class:`~anndata.AnnData` object, potentially from an external
        source with unknown preprocessing history.

    Returns
    -------
    :class:`AnnDataState` with matrix state, completed steps, metadata
    keys, and warnings.
    """
    warnings: list[str] = []

    # ── X matrix state ──────────────────────────────────────────────
    x_state = _detect_x_state(adata, warnings)

    # ── Raw counts location ─────────────────────────────────────────
    raw_location = _find_raw_counts_location(adata, x_state)
    if x_state == SCALED and raw_location is None:
        warnings.append(
            "X is z-score scaled but no raw counts found in layers or .raw. "
            "Re-processing from raw counts will not be possible."
        )

    # ── Completed steps (breadcrumbs) ───────────────────────────────
    has_qc = "n_genes_by_counts" in getattr(adata, "obs", {}).keys()
    has_norm = "log1p" in getattr(adata, "uns", {})
    has_raw_frozen = adata.raw is not None
    has_hvg = "highly_variable" in getattr(adata, "var", {}).keys()
    has_pca = "X_pca" in getattr(adata, "obsm", {})
    has_neighbors = (
        "neighbors" in getattr(adata, "uns", {})
        and "connectivities" in getattr(adata, "obsp", {})
    )
    has_umap = "X_umap" in getattr(adata, "obsm", {})
    has_tsne = "X_tsne" in getattr(adata, "obsm", {})
    has_de = "rank_genes_groups" in getattr(adata, "uns", {})
    has_diffmap = "X_diffmap" in getattr(adata, "obsm", {})
    has_velocity = (
        "spliced" in getattr(adata, "layers", {})
        and "unspliced" in getattr(adata, "layers", {})
    )

    # ── Metadata detection (fuzzy column matching) ──────────────────
    obs_cols = list(getattr(adata.obs, "columns", []))
    cluster_key = _find_cluster_key(adata, obs_cols)
    celltype_key = _find_celltype_key(obs_cols)
    condition_key = _find_condition_key(adata, obs_cols)
    batch_key = _find_batch_key(adata, obs_cols)

    # ── Species detection ───────────────────────────────────────────
    species = _detect_species(adata)

    return AnnDataState(
        x_state=x_state,
        raw_counts_location=raw_location,
        has_qc_metrics=has_qc,
        has_normalized=has_norm,
        has_raw_frozen=has_raw_frozen,
        has_hvg=has_hvg,
        has_pca=has_pca,
        has_neighbors=has_neighbors,
        has_umap=has_umap,
        has_tsne=has_tsne,
        has_clusters=cluster_key is not None,
        has_cell_types=celltype_key is not None,
        has_de_results=has_de,
        has_diffmap=has_diffmap,
        has_velocity_layers=has_velocity,
        cluster_key=cluster_key,
        celltype_key=celltype_key,
        condition_key=condition_key,
        batch_key=batch_key,
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        species=species,
        warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Raw counts retrieval
# ---------------------------------------------------------------------------


def find_raw_counts(adata, state: AnnDataState | None = None):
    """Retrieve the best available raw count matrix.

    Parameters
    ----------
    adata
        The AnnData object.
    state
        Pre-computed inspector state. If *None*, runs ``inspect_adata()``
        first.

    Returns
    -------
    A matrix (sparse or dense) of raw counts.

    Raises
    ------
    ValueError
        If no raw counts can be found anywhere.
    """
    if state is None:
        state = inspect_adata(adata)

    loc = state.raw_counts_location
    if loc is None:
        raise ValueError(
            "No raw counts found. Checked: layers['counts'], "
            "layers['raw_counts'], adata.X. "
            f"X is '{state.x_state}', which cannot be used as raw counts."
        )

    if loc == "X":
        return adata.X.copy()
    elif loc.startswith("layers["):
        # Extract layer name from "layers['counts']"
        layer_name = loc.split("'")[1]
        return adata.layers[layer_name].copy()

    raise ValueError(f"Unknown raw counts location: {loc}")


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def summarize_state(state: AnnDataState) -> str:
    """Human-readable Markdown summary of the AnnData state."""
    lines: list[str] = []

    # Header
    species_str = f" ({state.species})" if state.species else ""
    lines.append("## Data Inspection\n")
    lines.append(
        f"**Shape:** {state.n_cells:,} cells × {state.n_genes:,} genes{species_str}"
    )

    # X state
    x_desc = {
        RAW_COUNTS: "raw integer counts",
        LOG_NORMALIZED: "log-normalized",
        SCALED: "z-score scaled (negative values detected)",
        TRANSFORMED_UNKNOWN: "transformed (unknown method)",
    }
    lines.append(f"**X state:** {x_desc.get(state.x_state, state.x_state)}")

    # Raw counts
    if state.raw_counts_location:
        lines.append(f"**Raw counts:** available in `{state.raw_counts_location}` ✅")
    else:
        lines.append("**Raw counts:** not found ❌")

    # Completed steps
    lines.append("\n### Completed Steps")
    steps = [
        ("QC metrics", state.has_qc_metrics),
        ("Normalized", state.has_normalized),
        ("HVG", state.has_hvg),
        ("PCA", state.has_pca),
        ("Neighbors", state.has_neighbors),
        ("UMAP", state.has_umap),
        ("Clustered", state.has_clusters),
        ("Cell types", state.has_cell_types),
        ("DE results", state.has_de_results),
    ]
    step_strs = [f"{'✅' if done else '❌'} {name}" for name, done in steps]
    # Two rows for readability
    mid = (len(step_strs) + 1) // 2
    lines.append("  ".join(step_strs[:mid]))
    lines.append("  ".join(step_strs[mid:]))

    # Metadata
    metadata_items = []
    if state.cluster_key:
        n = _safe_nunique(state, "cluster_key")
        metadata_items.append(
            f"- Clusters: `obs['{state.cluster_key}']`"
            + (f" ({n} clusters)" if n else "")
        )
    if state.celltype_key:
        metadata_items.append(f"- Cell types: `obs['{state.celltype_key}']`")
    if state.condition_key:
        metadata_items.append(f"- Condition: `obs['{state.condition_key}']`")
    if state.batch_key:
        metadata_items.append(f"- Batch: `obs['{state.batch_key}']`")

    if metadata_items:
        lines.append("\n### Metadata")
        lines.extend(metadata_items)

    # Warnings
    if state.warnings:
        lines.append("\n### ⚠️ Warnings")
        for w in state.warnings:
            lines.append(f"- {w}")

    return "\n".join(lines)


def _safe_nunique(state: AnnDataState, key_attr: str) -> int | None:
    """Try to get nunique for a key — but we don't store adata, so return None."""
    # The summary doesn't have access to adata, so we can't compute nunique.
    # This is a presentation-only helper; the caller can enrich the summary
    # if they have adata available.
    return None


# ---------------------------------------------------------------------------
# X state detection
# ---------------------------------------------------------------------------


def _detect_x_state(adata, warnings: list[str]) -> str:
    """Classify the state of adata.X.

    Priority order:
    1. Negative values → 'scaled'
    2. Integer values → 'raw_counts'
    3. uns['log1p'] present → 'log_normalized'
    4. max < 20 heuristic → 'log_normalized' (with warning)
    5. else → 'transformed_unknown'
    """
    X = adata.X
    n_rows = min(200, adata.n_obs)

    if issparse(X):
        sample = X[:n_rows].toarray()
    else:
        sample = X[:n_rows]

    sample_flat = sample.flatten()

    # Check 1: negatives → scaled
    if (sample_flat < 0).any():
        return SCALED

    # Check 2: integers → raw counts
    # Use a reasonable sample to avoid scanning the whole matrix
    check_vals = sample_flat[sample_flat != 0]
    if len(check_vals) == 0:
        # All zeros — treat as raw counts (empty/trivial data)
        return RAW_COUNTS
    if np.allclose(check_vals[:5000], np.round(check_vals[:5000])):
        return RAW_COUNTS

    # Check 3: Scanpy breadcrumb
    if "log1p" in getattr(adata, "uns", {}):
        return LOG_NORMALIZED

    # Check 4: value range heuristic
    if issparse(X):
        max_val = float(X.max())
    else:
        max_val = float(np.max(X))

    if max_val < 20:
        warnings.append(
            "No Scanpy log-transform marker (uns['log1p']) found. "
            "Inferring log-normalized from value range (max < 20). "
            "Data may have been processed outside Scanpy."
        )
        return LOG_NORMALIZED

    # Check 5: unknown
    return TRANSFORMED_UNKNOWN


# ---------------------------------------------------------------------------
# Raw counts location
# ---------------------------------------------------------------------------


def _find_raw_counts_location(adata, x_state: str) -> str | None:
    """Find where raw counts are stored.

    Search order:
    1. layers['counts']
    2. layers['raw_counts']
    3. X itself (if x_state is 'raw_counts')
    """
    layers = getattr(adata, "layers", {})

    for layer_name in ("counts", "raw_counts"):
        if layer_name in layers:
            # Quick validation: should be non-negative integers
            layer = layers[layer_name]
            if issparse(layer):
                sample = layer[:min(100, adata.n_obs)].toarray().flatten()
            else:
                sample = layer[:min(100, adata.n_obs)].flatten()

            non_neg = (sample >= 0).all()
            is_int = np.allclose(
                sample[sample != 0][:1000],
                np.round(sample[sample != 0][:1000]),
            ) if (sample != 0).any() else True

            if non_neg and is_int:
                return f"layers['{layer_name}']"

    if x_state == RAW_COUNTS:
        return "X"

    return None


# ---------------------------------------------------------------------------
# Fuzzy column matching
# ---------------------------------------------------------------------------

_CLUSTER_PREFIXES = ("leiden", "louvain", "cluster")
_CELLTYPE_PATTERNS = (
    "cell_type", "celltype", "cell_label", "cell_annotation",
    "annotation", "celltypist", "azimuth", "singler",
)
_CONDITION_PATTERNS = (
    "condition", "disease", "treatment", "group", "status",
    "phenotype", "diagnosis",
)
_BATCH_PATTERNS = (
    "batch", "sample", "donor", "patient", "library_id",
    "channel", "library", "subject",
)


def _find_cluster_key(adata, obs_cols: list[str]) -> str | None:
    """Find clustering column — could be leiden, louvain, custom name."""
    for col in obs_cols:
        col_lower = col.lower()
        if any(col_lower.startswith(p) for p in _CLUSTER_PREFIXES):
            return col
    return None


def _find_celltype_key(obs_cols: list[str]) -> str | None:
    """Find cell type annotation column."""
    for col in obs_cols:
        col_lower = col.lower().replace(" ", "_")
        if any(pat in col_lower for pat in _CELLTYPE_PATTERNS):
            return col
    return None


def _find_condition_key(adata, obs_cols: list[str]) -> str | None:
    """Find experimental condition column."""
    for col in obs_cols:
        col_lower = col.lower()
        if any(pat in col_lower for pat in _CONDITION_PATTERNS):
            try:
                if adata.obs[col].nunique() >= 2:
                    return col
            except Exception:
                pass
    return None


def _find_batch_key(adata, obs_cols: list[str]) -> str | None:
    """Find batch/sample column."""
    for col in obs_cols:
        col_lower = col.lower()
        if any(pat in col_lower for pat in _BATCH_PATTERNS):
            try:
                if adata.obs[col].nunique() >= 2:
                    return col
            except Exception:
                pass
    return None


# ---------------------------------------------------------------------------
# Species detection
# ---------------------------------------------------------------------------


def _detect_species(adata) -> str | None:
    """Detect species from gene name casing.

    Human genes: uppercase (BRCA1, TP53, MT-CO1)
    Mouse genes: title case (Brca1, Tp53, mt-Co1)

    Some gene names are contig/locus identifiers (e.g., RP23-317L18.4,
    AC139064.2) which are uppercase regardless of species.  We filter
    these out by only considering short alphabetic names for the casing
    check, and also look for the mitochondrial prefix convention.
    """
    var_names = list(adata.var_names[:1000])
    if not var_names:
        return None

    # Check mitochondrial prefix — most reliable single signal
    has_mt_human = any(g.startswith("MT-") for g in var_names)
    has_mt_mouse = any(g.startswith("mt-") for g in var_names)

    if has_mt_human and not has_mt_mouse:
        return "human"
    if has_mt_mouse and not has_mt_human:
        return "mouse"

    # Filter to "normal" gene names: alphabetic, 2-10 chars, no digits
    # This excludes contig IDs like RP23-317L18.4 that skew the count
    normal = [g for g in var_names if g.isalpha() and 2 <= len(g) <= 10]
    if len(normal) < 50:
        # Not enough clean gene names to judge; fall back to full list
        normal = var_names

    n_upper = sum(1 for g in normal if g == g.upper())
    frac_upper = n_upper / len(normal)

    if frac_upper > 0.7:
        return "human"
    elif frac_upper < 0.3:
        return "mouse"
    return None
