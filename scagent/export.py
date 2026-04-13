"""Export analysis provenance as methods text or reproducibility package.

Two export modes:

1. **Methods section** — human-readable paragraph suitable for a paper's
   Methods section, generated from the provenance chain.
2. **Reproducibility package** — a directory with a replay script,
   parameters JSON, software versions, and data manifest.

Usage::

    from scagent.export import generate_methods, generate_repro_package

    methods = generate_methods(provenance_graph, experiment_context)
    generate_repro_package(provenance_graph, experiment_context, out_dir="repro/")
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scagent.provenance import ProvenanceGraph
from scagent.context import ExperimentContext


# ---------------------------------------------------------------------------
# Tool → prose templates
# ---------------------------------------------------------------------------

_TOOL_PROSE: dict[str, str] = {
    "load_10x_h5": (
        "Raw count matrices were loaded from 10x Genomics Cell Ranger output"
        " ({file_path})."
    ),
    "qc_metrics": (
        "Quality control metrics were computed, including the number of genes"
        " per cell, total UMI counts, and the percentage of mitochondrial reads."
    ),
    "filter_cells": (
        "Cells were filtered based on: {filter_desc},"
        " retaining {n_cells_after} cells from {n_cells_before}."
    ),
    "filter_genes": (
        "Genes were filtered, requiring expression in at least"
        " {min_cells} cells, retaining {n_genes_after} genes"
        " from {n_genes_before}."
    ),
    "detect_doublets": (
        "Doublet detection was performed using Scrublet"
        " (expected_doublet_rate={expected_doublet_rate}),"
        " identifying {n_doublets} predicted doublets ({doublet_pct}%)."
    ),
    "normalize": (
        "Data was normalized using size-factor normalization"
        " (target_sum={target_sum})."
    ),
    "log_transform": "Counts were log1p-transformed.",
    "highly_variable_genes": (
        "Highly variable genes were selected using the {flavor} method"
        " (n_top_genes={n_top_genes})."
    ),
    "scale": "Gene expression was scaled to unit variance and zero mean.",
    "pca": (
        "Principal component analysis was performed retaining"
        " {n_comps} components."
    ),
    "batch_correction": (
        "Batch correction was applied using {method} on the"
        " '{batch_key}' variable."
    ),
    "neighbors": (
        "A neighborhood graph was computed using {n_neighbors} neighbors"
        " and {n_pcs} principal components."
    ),
    "leiden": (
        "Cells were clustered using the Leiden algorithm at resolution"
        " {resolution}, yielding {n_clusters} clusters."
    ),
    "umap": "UMAP embedding was computed for visualization.",
    "rank_genes_groups": (
        "Marker genes were identified using the {method} test"
        " ({groupby})."
    ),
    "annotate_celltypist": (
        "Automated cell type annotation was performed using CellTypist"
        " (model: {model})."
    ),
    "deseq2_pseudobulk": (
        "Differential expression analysis was performed using pseudobulk"
        " aggregation with DESeq2 (design: {design}, contrast:"
        " {contrast})."
    ),
    "gsea": (
        "Gene set enrichment analysis was performed using GSEApy"
        " (gene_sets: {gene_sets})."
    ),
    "custom": "Custom analysis: {description}",
}

# Chain dicts use 'extras' for effects/description. These helpers extract them.
def _get_effects(step: dict) -> dict:
    extras = step.get("extras", {})
    return extras.get("effects", extras)

def _get_description(step: dict) -> str:
    extras = step.get("extras", {})
    return extras.get("description", step.get("description", ""))


# ---------------------------------------------------------------------------
# Methods section generator
# ---------------------------------------------------------------------------

def generate_methods(
    graph: ProvenanceGraph,
    context: ExperimentContext | None = None,
    branch: str | None = None,
) -> str:
    """Generate a methods section from provenance.

    Parameters
    ----------
    graph
        The provenance graph.
    context
        Experiment context (adds organism/tissue/paradigm header).
    branch
        Branch to export. Default: promoted branch or ``main``.

    Returns
    -------
    Methods text as a string (Markdown-compatible).
    """
    branch = branch or graph._promoted_branch or "main"
    chain = graph.get_full_chain(branch)

    paragraphs: list[str] = []

    # Preamble from experiment context
    if context:
        preamble = _context_preamble(context)
        if preamble:
            paragraphs.append(preamble)

    # Software
    sw = _software_line(graph)
    if sw:
        paragraphs.append(sw)

    # Tool steps → prose
    step_lines: list[str] = []
    for step in chain:
        line = _step_to_prose(step)
        if line:
            step_lines.append(line)
    if step_lines:
        paragraphs.append(" ".join(step_lines))

    return "\n\n".join(paragraphs)


def _context_preamble(ctx: ExperimentContext) -> str:
    parts = []
    org = ctx.organism
    if org:
        species = org.get("species", "")
        if species:
            parts.append(f"Single-cell RNA sequencing was performed on {species}")
    tissue = ctx.tissue
    if tissue:
        name = tissue.get("name", "")
        if name:
            parts.append(f"{name} samples" if parts else f"Samples from {name}")
    platform = ctx.platform
    if platform:
        pname = platform.get("name", "") if isinstance(platform, dict) else str(platform)
        if pname:
            parts.append(f"using the {pname} platform")
    library = ctx.library
    if library:
        ltype = library.get("type", "") if isinstance(library, dict) else str(library)
        if ltype:
            parts.append(f"({ltype} library preparation)")

    if not parts:
        return ""
    return " ".join(parts) + "."


def _software_line(graph: ProvenanceGraph) -> str:
    # Get software versions from the provenance session
    sessions = graph.sessions
    if not sessions:
        return ""
    sw = sessions[0].get("software_versions", {})
    if not sw:
        return ""
    parts = [f"{k} v{v}" for k, v in sw.items() if v]
    if parts:
        return f"Analysis was performed using {', '.join(parts)}."
    return ""


def _step_to_prose(step: dict) -> str:
    tool_id = step.get("tool_id", "unknown")
    template = _TOOL_PROSE.get(tool_id)
    if template is None:
        desc = _get_description(step)
        if desc:
            return desc
        return ""

    params = step.get("parameters", {})
    effects = _get_effects(step)
    merged = {**params, **effects}

    # Build filter description for filter_cells
    if tool_id == "filter_cells":
        filter_parts = []
        for k in ("min_genes", "max_genes", "min_counts", "max_counts",
                   "max_pct_mito", "max_pct_ribo"):
            if k in params:
                filter_parts.append(f"{k}={params[k]}")
        merged["filter_desc"] = ", ".join(filter_parts) if filter_parts else "custom thresholds"

    # Doublet percentage
    if tool_id == "detect_doublets":
        n_doub = effects.get("n_doublets", 0)
        n_cells = effects.get("n_cells", 1)
        merged["doublet_pct"] = f"{100 * n_doub / n_cells:.1f}" if n_cells else "?"

    try:
        return template.format(**merged)
    except KeyError:
        # Template has placeholders we don't have data for — partial fill
        import re
        result = template
        for match in re.finditer(r"\{(\w+)\}", template):
            key = match.group(1)
            val = merged.get(key, "?")
            result = result.replace(f"{{{key}}}", str(val))
        return result


# ---------------------------------------------------------------------------
# Reproducibility package
# ---------------------------------------------------------------------------

def generate_repro_package(
    graph: ProvenanceGraph,
    context: ExperimentContext | None = None,
    out_dir: str | Path = "repro",
    branch: str | None = None,
) -> Path:
    """Generate a reproducibility package directory.

    Creates:
    - ``replay.py`` — Python script that replays the analysis
    - ``params.json`` — all parameters used
    - ``environment.json`` — software versions and platform
    - ``provenance.json`` — full PROV-O JSON-LD
    - ``methods.md`` — generated methods section
    - ``context.json`` — experiment context (if available)

    Parameters
    ----------
    graph
        The provenance graph.
    context
        Experiment context.
    out_dir
        Output directory.
    branch
        Branch to export.

    Returns
    -------
    Path to the output directory.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    branch = branch or graph._promoted_branch or "main"

    chain = graph.get_full_chain(branch)

    # 1. params.json — flat list of all tool calls + params
    params = []
    for step in chain:
        params.append({
            "step": len(params) + 1,
            "tool_id": step.get("tool_id", "unknown"),
            "parameters": step.get("parameters", {}),
            "extras": step.get("extras", {}),
        })
    _write_json(out / "params.json", params)

    # 2. environment.json
    sessions = graph.sessions
    env = dict(sessions[0]) if sessions else {}
    env["exported_at"] = datetime.now(timezone.utc).isoformat()
    env["branch"] = branch
    _write_json(out / "environment.json", env)

    # 3. provenance.json — full graph
    _write_json(out / "provenance.json", graph.serialize())

    # 4. methods.md
    methods = generate_methods(graph, context, branch)
    (out / "methods.md").write_text(methods, encoding="utf-8")

    # 5. context.json
    if context:
        _write_json(out / "context.json", context.raw)

    # 6. replay.py
    script = _generate_replay_script(chain)
    (out / "replay.py").write_text(script, encoding="utf-8")

    # 7. README
    readme = _generate_readme(branch, len(chain), env)
    (out / "README.md").write_text(readme, encoding="utf-8")

    return out


def _generate_replay_script(chain: list[dict]) -> str:
    """Generate a Python script that replays the analysis chain."""
    lines = [
        '#!/usr/bin/env python3',
        '"""Auto-generated replay script from scAgent provenance."""',
        '',
        'import scanpy as sc',
        'from pathlib import Path',
        '',
        '# ─── Replay ────────────────────────────────────────────',
        '',
    ]

    for i, step in enumerate(chain):
        tool_id = step.get("tool_id", "unknown")
        params = step.get("parameters", {})
        desc = step.get("description", "")
        lines.append(f"# Step {i + 1}: {tool_id}")
        if desc:
            lines.append(f"# {desc}")
        lines.append(f"# Parameters: {json.dumps(params, default=str)}")
        code = _tool_to_code(tool_id, params, i)
        lines.append(code)
        lines.append("")

    lines.append('print("Replay complete.")')
    return "\n".join(lines)


# Mapping tool_id → code template
_TOOL_CODE: dict[str, str] = {
    "load_10x_h5": 'adata = sc.read_10x_h5("{file_path}")\nadata.var_names_make_unique()',
    "qc_metrics": 'sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)',
    "filter_cells": (
        "adata = adata[\n"
        "    (adata.obs['n_genes_by_counts'] >= {min_genes})\n"
        "    & (adata.obs['n_genes_by_counts'] <= {max_genes})\n"
        "    & (adata.obs['pct_counts_mt'] <= {max_pct_mito})\n"
        "].copy()"
    ),
    "filter_genes": 'sc.pp.filter_genes(adata, min_cells={min_cells})',
    "detect_doublets": (
        "import scrublet as scr\n"
        "scrub = scr.Scrublet(adata.X, expected_doublet_rate={expected_doublet_rate})\n"
        "doublet_scores, predicted_doublets = scrub.scrub_doublets()\n"
        "adata = adata[~predicted_doublets].copy()"
    ),
    "normalize": "sc.pp.normalize_total(adata, target_sum={target_sum})",
    "log_transform": "sc.pp.log1p(adata)",
    "highly_variable_genes": 'sc.pp.highly_variable_genes(adata, flavor="{flavor}", n_top_genes={n_top_genes})',
    "scale": "sc.pp.scale(adata, max_value={max_value})",
    "pca": "sc.tl.pca(adata, n_comps={n_comps})",
    "batch_correction": 'sc.external.pp.harmony_integrate(adata, key="{batch_key}")',
    "neighbors": "sc.pp.neighbors(adata, n_neighbors={n_neighbors}, n_pcs={n_pcs})",
    "leiden": 'sc.tl.leiden(adata, resolution={resolution}, key_added="leiden")',
    "umap": "sc.tl.umap(adata)",
    "rank_genes_groups": 'sc.tl.rank_genes_groups(adata, groupby="{groupby}", method="{method}")',
    "annotate_celltypist": (
        "import celltypist\n"
        'model = celltypist.models.Model.load(model="{model}")\n'
        "predictions = celltypist.annotate(adata, model=model, majority_voting=True)"
    ),
}


def _tool_to_code(tool_id: str, params: dict, step_idx: int) -> str:
    template = _TOOL_CODE.get(tool_id)
    if template is None:
        return f"# TODO: {tool_id}({json.dumps(params, default=str)})"
    import re
    result = template
    for match in re.finditer(r"\{(\w+)\}", template):
        key = match.group(1)
        val = params.get(key, "?")
        result = result.replace(f"{{{key}}}", str(val))
    return result


def _generate_readme(branch: str, n_steps: int, env: dict) -> str:
    sw = env.get("software_versions", {})
    sw_lines = "\n".join(f"- {k}: {v}" for k, v in sw.items()) if sw else "- (see environment.json)"
    return f"""# Reproducibility Package

Auto-generated by scAgent from provenance records.

## Analysis

- **Branch:** {branch}
- **Steps:** {n_steps}
- **Exported:** {env.get("exported_at", "?")}

## Contents

| File | Description |
|------|-------------|
| `replay.py` | Python script that replays the analysis |
| `params.json` | All tool parameters |
| `environment.json` | Software versions and platform |
| `provenance.json` | Full W3C PROV-O provenance graph |
| `methods.md` | Generated methods section text |
| `context.json` | Experiment metadata (if available) |
| `README.md` | This file |

## Software

{sw_lines}

## Usage

```bash
pip install scagent
python replay.py
```
"""


def _write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
