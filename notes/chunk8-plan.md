# Chunk 8 Design Plan: DE Pipeline (Pseudobulk + Enrichment)

**Status:** Building
**Depends on:** Chunk 6 (DAG enforces pseudobulk for disease_vs_healthy), tool registry
**Core insight:** scBench shows 41% DE accuracy for best LLM — because models default to cell-level Wilcoxon. We structurally prevent that.

---

## What We Build

1. **`scagent/tools/pseudobulk_de.py`** — Pseudobulk aggregation + PyDESeq2 differential expression
2. **`scagent/tools/enrichment.py`** — GSEA via GSEApy using DE results
3. **`.pi/skills/pseudobulk-de/SKILL.md`** — Skill teaching the agent when & how to use pseudobulk DE
4. **`.pi/skills/pathway-enrichment/SKILL.md`** — Skill for enrichment analysis

## Key Design Decisions

- **PyDESeq2** (pure Python) instead of R DESeq2 via subprocess — no R dependency, same algorithm
- **Aggregation is inside the tool** — the agent calls one function, it handles sum-per-celltype-per-sample internally
- **Per-cell-type DE** — runs DESeq2 separately for each cell type, returns results keyed by cell type
- **Guard rails:** refuse to run if <2 replicates per condition, warn if <3

## Implementation Checklist

```
├── [ ] scagent/tools/pseudobulk_de.py
│     ├── aggregate_pseudobulk() — sum raw counts per cell_type × sample
│     ├── run_deseq2() — PyDESeq2 per cell type
│     ├── Validation: min replicates, min cells per pseudobulk
│     └── Output: DE results + volcano plots + provenance dict
├── [ ] scagent/tools/enrichment.py
│     ├── run_gsea() — GSEApy prerank from DE results
│     └── Output: enrichment table + plots + provenance dict
├── [ ] .pi/skills/pseudobulk-de/SKILL.md
├── [ ] .pi/skills/pathway-enrichment/SKILL.md
├── [ ] tests/test_pseudobulk_de.py — unit tests with synthetic data
├── [ ] tests/test_enrichment.py — unit tests
└── [ ] tests/test_de_integration.py — end-to-end on synthetic multi-sample data
```
