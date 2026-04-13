---
name: knowledge
description: Query marker genes from curated databases and validate cell type annotations. Use when the user asks about marker genes, wants to verify annotations, or needs canonical markers for a cell type.
---

# Skill: Knowledge (Marker Gene Database)

Query marker genes from curated databases and validate cell type annotations.

## Sources (priority order)

| Source | Coverage | Availability |
|--------|----------|-------------|
| **Canonical** | ~30 PBMC/immune cell types, curated markers | Always available, no download |
| **CellTypist** | 98 immune subtypes (Low model), 61 models total | Needs model download (~5MB each) |
| **CellMarker 2.0** | 400+ cell types, multi-tissue | User drops xlsx in `.scagent/knowledge/` |
| **PanglaoDB** | 200+ cell types | User drops TSV in `.scagent/knowledge/` |

## When to use

### Query markers

Before or during cell type annotation:

```python
from scagent.tools.knowledge_tools import query_markers
result = query_markers("NK cells", species="human")
# → {markers: [{gene: "NKG7", source: "canonical", rank: 1}, ...]}
```

### Validate annotation

After `rank_genes_groups` or manual labelling — check if DE markers support
the label:

```python
from scagent.tools.knowledge_tools import validate_annotation
result = validate_annotation(
    cluster_markers=["NKG7", "GNLY", "GZMB", "PRF1", "FCGR3A"],
    label="NK cells",
)
# → {confidence: "high", overlap_ratio: 0.625, matched_markers: ["FCGR3A", "GNLY", ...]}
```

### During annotation workflow

1. Run `rank_genes_groups` to get top markers per cluster
2. For each cluster with a proposed label, run `validate_annotation`
3. If confidence is "low", check the `alternatives` field for better matches
4. If confidence is "medium", present both the label and alternatives to the user

## Cell type aliases

Common aliases are resolved automatically:
- "Tregs" → "Regulatory T cells"
- "CD14 monocytes" → "Classical monocytes"
- "NK" → "NK cells"
- "pDC" → "pDC"

## Adding external databases

Tell the user to download and place files in `.scagent/knowledge/`:
- **CellMarker 2.0**: `Human_cell_markers.xlsx` from http://bio-bigdata.hrbmu.edu.cn/CellMarker/
- **PanglaoDB**: `PanglaoDB_markers.tsv` from https://panglaodb.se/markers.html

Files are auto-detected by name pattern on first query.
