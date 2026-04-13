# Best Practice Reference: Pathway / Gene Set Enrichment

> Distilled from Heumos et al. 2023.

## Purpose

Summarize many DE genes into interpretable biological terms — pathways, processes, transcription factor activities.

## Key Finding

- "Gene set enrichment analysis was found to be **more sensitive to the choice of gene sets rather than statistical methods**" (Heumos et al.).
- Therefore: **choose your database carefully** to ensure relevant pathways are covered. The statistical method matters less.

## Gene Set Databases

| Database | Scope |
|----------|-------|
| **MSigDB** (Hallmarks) | 50 curated cancer/immune hallmark pathways — good default |
| **Gene Ontology** | Broad biological processes, molecular functions, cellular components |
| **KEGG** | Metabolic and signaling pathways |
| **Reactome** | Detailed molecular pathways |
| **PROGENy** | Signalling pathway activities (weighted gene sets) |
| **DoRothEA** | Transcription factor activities (weighted gene sets) |

## Methods

- Common methods: hypergeometric/Fisher's exact test, GSEA, GSVA.
- **decoupleR** provides access to multiple databases and methods in a single tool — recommended for flexibility (Heumos et al.).
- Enrichment methods designed for bulk transcriptomics can be applied to scRNA-seq (Heumos et al.), but some single-cell methods like Pagoda2 may outperform them.
- GSEA uses a permissive FDR < 0.25 threshold by design — it detects coordinated small effects across gene sets.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
