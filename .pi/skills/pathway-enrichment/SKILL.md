# Skill: Pathway Enrichment (GSEA)

Run gene set enrichment analysis on DE results to find enriched biological pathways.

## When to Use

- **After pseudobulk DE** — the DAG places this after `pseudobulk_de`.
- **User asks about pathways, GO terms, biological processes, or functional enrichment.**

## How to Use

```python
from scagent.tools.enrichment import run_gsea

# de_results from run_pseudobulk_de()
enrichment = run_gsea(
    de_results["de_dataframes"],
    gene_sets="MSigDB_Hallmark_2020",
    fdr_threshold=0.25,
)

for row in enrichment["summary"]:
    print(f"{row['cell_type']}: {row['n_significant']} enriched terms")
```

## Gene Set Choices

Recommend based on the biological question:

| Gene set | When to use |
|----------|------------|
| `MSigDB_Hallmark_2020` | Default — 50 curated cancer/immune hallmarks |
| `GO_Biological_Process_2023` | Broad biological processes |
| `KEGG_2021_Human` | Metabolic and signaling pathways |
| `Reactome_2022` | Detailed molecular pathways |

## What to Tell the User

**Before running:**
> "Running GSEA on DE results using Hallmark gene sets. This identifies pathways with coordinated up/down-regulation."

**After running:**
> "3 enriched pathways in CD14 Monocytes: TNF-alpha signaling (NES=2.1), Inflammatory response (NES=1.8), Interferon gamma response (NES=1.6). GSEA uses a permissive FDR<0.25 threshold by design."

## Notes

- GSEA's 0.25 FDR threshold is intentionally more permissive than typical DE thresholds — it detects coordinated small effects across gene sets.
- If no significant terms are found, suggest trying different gene sets or checking the DE results.
- For mouse data, some gene set libraries need gene name conversion (upper-case to title-case).
