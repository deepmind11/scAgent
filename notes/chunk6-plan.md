# Chunk 6 Design Plan: Experiment Context & DAG

**Status:** Draft for review
**Depends on:** Chunks 4–5 (provenance + state), tool registry (Chunk 1)
**Feeds into:** Chunk 8 (DE pipeline uses paradigm to enforce pseudobulk), Chunk 11 (full paradigm support)

---

## 0. What Chunk 6 Builds

Three things:

1. **Experiment context** (`scagent-project.json`) — structured metadata about the experiment, guided by minSCe (Füllgrabe et al. 2020). This is the agent's "ground truth."
2. **Analysis DAG** — a paradigm-specific directed acyclic graph of analysis steps. The agent uses this to suggest next steps and validate ordering.
3. **Init skill** — a conversational flow that gathers context from the researcher, infers what it can from the data, and asks about what it can't.

---

## 1. The Init UX: Infer First, Ask Only What You Must

The init flow is **not a form.** It's a conversation where the agent is smart about what it already knows.

### 1.1 What Can Be Inferred from Data

| Field | How | Confidence |
|-------|-----|------------|
| Species | Gene name prefixes: `MT-` = human, `mt-` = mouse | High (already in QC) |
| Number of cells | `.h5` file shape | Exact |
| Number of genes | `.h5` file shape | Exact |
| Platform | File format (`.h5` = Cell Ranger, `.loom` = other) | Medium |
| Chemistry | Cell Ranger `metrics_summary.csv` if present | Medium |
| Multi-sample | Number of input files, or `batch` column in AnnData | Medium |
| UMI-based | Presence of Cell Ranger output structure | High |

### 1.2 What Cannot Be Inferred — Must Ask

| Field | Why It Matters | Gates |
|-------|---------------|-------|
| **Paradigm** | The biological question determines the entire DAG | DAG selection |
| **Tissue** | Expected cell types, marker databases to load | Annotation strategy |
| **Conditions** | Which samples are disease vs. healthy, treated vs. control | Whether DE / differential abundance applies |
| **Sample ↔ condition mapping** | Which files belong to which group | Pseudobulk aggregation, batch correction |
| **Batch structure** | Are batches confounded with conditions? | Whether batch correction is needed, which method |

### 1.3 Three UX Scenarios

**Scenario A: Rich prompt**
> *"I have PBMCs from 4 donors — 2 young, 2 old — on 10x Chromium v3. Looking for age-related immune changes."*

Agent extracts everything, confirms, generates context + DAG.

**Scenario B: Minimal prompt**
> *"Analyze these files"*

Agent loads data, infers species + platform, then asks:
```
I've loaded your data: 11,769 cells × 33,538 genes. Gene names suggest human.

To plan your analysis, I need a few things:

1. What tissue are these cells from?
2. What's your experimental question?
   • Characterize cell types in a tissue (cell atlas)
   • Compare conditions — e.g., disease vs. healthy, treated vs. control
   • Track a developmental / differentiation trajectory
   • Something else?
3. Do your samples have different conditions? If so, which?
```

**Scenario C: No context at all (data just dropped)**

Agent detects no `scagent-project.json`, proactively asks before proceeding with analysis. Does NOT silently default to cell atlas.

### 1.4 minSCe Checklist — What the Skill Asks

Based on Füllgrabe et al. 2020, Figure 2. We group into **analysis-critical** (must ask) and **metadata** (nice to have, ask only if the researcher offers):

**Analysis-critical (gates decisions):**

| minSCe Component | Field | Question |
|------------------|-------|----------|
| Biosource | Species | (inferred, confirm) |
| Biosource | Tissue / organ | "What tissue are these cells from?" |
| Single cell isolation | Cell enrichment (FACS sorted?) | "Were cells enriched or sorted before capture?" |
| Library construction | End bias (3'/5'/full-length) | (inferred from platform, confirm if unclear) |
| Library construction | Modalities (GEX only? + protein? + VDJ?) | "Is this gene expression only, or multiome?" |
| Design | Paradigm | "What's your biological question?" |
| Design | Conditions | "Do your samples belong to different conditions?" |
| Design | Batch structure | "Were samples processed in separate batches?" |

**Metadata (record if provided, don't block on it):**

| minSCe Component | Field | Notes |
|------------------|-------|-------|
| Biosource | Donor info (sex, age) | Helpful for DE covariates |
| Biosource | Cell line vs. primary | Changes QC expectations |
| Single cell isolation | Dissociation method | Affects stress gene signatures |
| Sequencing | Instrument | Usually in Cell Ranger output |
| Single cell isolation | Target cell number | Useful for doublet rate estimation |

---

## 2. Experiment Context Schema

File: `.scagent/project.json`

```jsonc
{
  "@context": "https://scagent.dev/context/v1",
  "@type": "SingleCellExperiment",

  // ── Identity ──
  "project_id": "proj_2026_04_pbmc_aging",     // auto-generated
  "created": "2026-04-11T10:00:00Z",

  // ── Paradigm (gates the DAG) ──
  "paradigm": "disease_vs_healthy",
  // Valid: cell_atlas, disease_vs_healthy, developmental_trajectory,
  //        perturbation_screen, temporal_longitudinal

  // ── Biosource (minSCe) ──
  "organism": {
    "species": "Homo sapiens",
    "ncbi_taxon": 9606,
    "inferred": true                            // was this auto-detected?
  },
  "tissue": {
    "name": "peripheral blood mononuclear cells",
    "uberon_id": "UBERON:0000178"               // optional, filled if known
  },
  "biosource_type": "specimen_from_organism",   // or: cell_line, organoid

  // ── Platform (minSCe: library construction + sequencing) ──
  "platform": {
    "vendor": "10x Genomics",
    "instrument": "Chromium",
    "chemistry": "3' v3",
    "cell_ranger_version": "9.0.0",
    "inferred": true
  },
  "library": {
    "type": "3prime",                           // 3prime, 5prime, full_length
    "modalities": ["gene_expression"],
    "umi": true
  },

  // ── Samples ──
  "samples": [
    {"id": "donor_1", "condition": "young", "file": "donor1.h5"},
    {"id": "donor_2", "condition": "old", "file": "donor2.h5"}
  ],

  // ── Experimental design ──
  "design": {
    "type": "case_control",                     // case_control, time_series, single_condition
    "conditions": ["young", "old"],
    "biological_replicates_per_condition": 2,
    "batch_key": "sample",                      // column in adata.obs for batch correction
    "batch_confounded_with_condition": false
  },

  // ── Single cell isolation (minSCe) ──
  "isolation": {
    "method": "microfluidics_droplet",          // microfluidics_droplet, facs, nanowell
    "enrichment": null,                         // e.g., "CD45+ FACS sort"
    "expected_doublet_rate": 0.04               // from 10x loading guidelines
  },

  // ── Hypotheses (optional) ──
  "hypotheses": [],

  // ── Status ──
  "status": "initialized"                       // initialized, in_progress, complete
}
```

### 2.1 Validation Rules

The `ExperimentContext` class validates:
- `paradigm` must be one of the known paradigm strings
- `organism.species` must be set (at minimum inferred)
- `tissue.name` must be set
- If paradigm is `disease_vs_healthy`: `design.conditions` must have ≥ 2 entries
- If paradigm is `developmental_trajectory`: conditions not required
- `samples` list must be non-empty
- `library.type` must be one of `3prime`, `5prime`, `full_length`

### 2.2 `inferred` Flags

Fields that were auto-detected (not explicitly stated by researcher) are marked `inferred: true`. The agent should confirm these: *"Gene names suggest human — is that correct?"*

---

## 3. Analysis DAG

### 3.1 Data Model

```python
@dataclass
class DAGStep:
    id: str                     # e.g., "qc_filtering"
    tool_id: str | None         # maps to tool registry; None for meta-steps
    category: str               # e.g., "qc", "clustering", "de"
    status: str                 # "pending", "running", "done", "skipped"
    depends_on: list[str]       # step IDs this depends on
    required: bool              # can this step be skipped?
    condition: str | None       # e.g., "design.batch_key is not null" → batch correction only if multi-batch

@dataclass
class AnalysisDAG:
    paradigm: str
    steps: list[DAGStep]
    current_step: str | None    # which step the researcher is at
```

### 3.2 DAG Definitions for 3 Paradigms

**Cell Atlas** (single-tissue deep profiling):
```
load → qc → normalize → hvg → pca → neighbors → clustering → umap
  → markers → annotation → report
```
No DE, no batch correction (unless multi-sample), no trajectory.

**Disease vs. Healthy** (cross-condition comparison):
```
load → qc → normalize → hvg → pca
  → [batch_correction if multi-batch]
  → neighbors → clustering → umap → markers → annotation
  → pseudobulk_de → pathway_enrichment
  → [composition_analysis]
  → report
```
Pseudobulk DE is required. Batch correction conditional on batch structure.

**Developmental Trajectory**:
```
load → qc → normalize → hvg → pca
  → [batch_correction if multi-batch]
  → neighbors → clustering → umap → markers → annotation
  → trajectory_inference → pseudotime
  → [rna_velocity]
  → report
```
Trajectory inference is the primary goal. DE is within-trajectory, not cross-condition.

### 3.3 Conditional Steps

Some steps only apply under certain conditions:

| Step | Condition | How to check |
|------|-----------|-------------|
| `batch_correction` | Multi-sample with batch structure | `len(context.samples) > 1 and context.design.batch_key` |
| `pseudobulk_de` | Paradigm has conditions to compare | `paradigm in (disease_vs_healthy, perturbation_screen)` |
| `composition_analysis` | Paradigm has conditions + cell types | Same as above, post-annotation |
| `trajectory_inference` | Paradigm is developmental | `paradigm == developmental_trajectory` |
| `rna_velocity` | Trajectory paradigm + spliced/unspliced available | Check for `.loom` or velocity-ready data |

### 3.4 DAG Operations

```python
class AnalysisDAG:
    def next_step() -> DAGStep | None        # suggest the next pending step
    def complete_step(step_id: str)           # mark as done
    def skip_step(step_id: str)               # mark as skipped (with reason)
    def is_valid_step(step_id: str) -> bool   # are all dependencies met?
    def summary() -> str                      # markdown status table
    def to_json() -> dict                     # serialize for .scagent/dag.json
    def from_paradigm(paradigm, context) -> "AnalysisDAG"   # factory
```

---

## 4. Module Design

### 4.1 File: `scagent/context.py`

```
ExperimentContext
├── __init__(project_dir: Path)
├── load() → ExperimentContext               # from .scagent/project.json
├── save()                                    # write to .scagent/project.json
├── from_prompt(text: str) → dict            # extract fields from free-text (for the skill)
├── infer_from_data(adata: AnnData) → dict   # auto-detect species, platform, etc.
├── validate() → list[str]                    # return list of missing/invalid fields
├── is_complete() → bool                      # all analysis-critical fields set?
├── needs_batch_correction() → bool
├── needs_pseudobulk_de() → bool
├── needs_trajectory() → bool
├── paradigm, organism, tissue, platform, ... # properties
```

### 4.2 File: `scagent/dag.py`

```
AnalysisDAG
├── from_paradigm(paradigm: str, context: ExperimentContext) → AnalysisDAG
├── next_step() → DAGStep | None
├── complete_step(step_id)
├── skip_step(step_id, reason)
├── is_valid_step(step_id) → bool
├── summary() → str
├── save(path) / load(path)                  # .scagent/dag.json
├── steps, current_step, paradigm            # properties
```

### 4.3 Skill: `.pi/skills/init/SKILL.md`

Teaches the agent the conversational init flow:

1. Check if `.scagent/project.json` exists. If yes, load it. If no, start init.
2. If data files are present, run `infer_from_data()` to auto-detect what we can.
3. Confirm inferred fields with the researcher.
4. Ask about the analysis-critical fields that couldn't be inferred (paradigm, tissue, conditions).
5. Don't ask about metadata fields unless the researcher volunteers info.
6. Generate `project.json` and `dag.json`.
7. Show the DAG to the researcher: "Here's your analysis plan. Want to modify anything?"
8. On subsequent sessions: if no `project.json`, ask before proceeding. Never default silently.

### 4.4 Skill: `.pi/skills/dag/SKILL.md`

Teaches the agent to follow the DAG:

1. Before each step, check `dag.next_step()` — suggest it to the researcher.
2. After each step, call `dag.complete_step(step_id)`.
3. If the researcher asks for something out of order, check `dag.is_valid_step()`.
4. If they want to skip a step, record it with `dag.skip_step(step_id, reason)`.
5. Show `dag.summary()` when asked about progress.

---

## 5. Storage Layout

```
.scagent/
├── project.json        # experiment context (this chunk)
├── dag.json            # current analysis DAG (this chunk)
├── provenance.jsonld   # (Chunk 4, exists)
├── state.json          # (Chunk 5, exists)
└── branches/           # (Chunk 5, exists)
```

---

## 6. Testing

### 6.1 Unit Tests: `tests/test_context.py`

| Test | What It Checks |
|------|----------------|
| `test_create_minimal_context` | Create with just paradigm + species + tissue, validate passes |
| `test_validate_missing_paradigm` | Validation fails if paradigm is missing |
| `test_validate_missing_tissue` | Validation fails if tissue is missing |
| `test_validate_conditions_required_for_dv_h` | disease_vs_healthy requires conditions |
| `test_validate_conditions_not_required_for_atlas` | cell_atlas doesn't require conditions |
| `test_infer_species_human` | Gene names with MT- → human |
| `test_infer_species_mouse` | Gene names with mt- → mouse |
| `test_needs_batch_correction` | Multi-sample + batch_key → True |
| `test_needs_pseudobulk_de` | disease_vs_healthy → True, cell_atlas → False |
| `test_save_load_roundtrip` | Write to JSON, load back, identical |
| `test_is_complete` | Returns True only when all critical fields set |

### 6.2 Unit Tests: `tests/test_dag.py`

| Test | What It Checks |
|------|----------------|
| `test_cell_atlas_dag_steps` | Cell atlas DAG has correct steps in correct order |
| `test_disease_vs_healthy_dag_has_pseudobulk` | disease_vs_healthy DAG includes pseudobulk_de |
| `test_trajectory_dag_has_trajectory` | developmental_trajectory DAG includes trajectory_inference |
| `test_next_step` | Returns first pending step |
| `test_complete_step` | Marks step done, advances next |
| `test_skip_step` | Marks step skipped with reason |
| `test_is_valid_step_dependencies` | Can't run clustering before neighbors |
| `test_conditional_batch_correction` | batch_correction included only when multi-batch |
| `test_summary_markdown` | summary() produces a readable table |
| `test_save_load_roundtrip` | Serialize to JSON, load back |
| `test_from_unknown_paradigm_raises` | Invalid paradigm raises error |

### 6.3 Integration Test: `tests/test_init_integration.py`

Real PBMC dataset:
1. Load data, run `infer_from_data()` → species=human, n_cells, n_genes
2. Create context with paradigm=cell_atlas, tissue=PBMC
3. Generate DAG, verify steps match expected cell_atlas pipeline
4. Walk through 3 steps, verify DAG tracking (complete, next_step)
5. Verify `project.json` and `dag.json` on disk

---

## 7. What This Chunk Does NOT Build

| Feature | Deferred To | Notes |
|---------|-------------|-------|
| All 9 paradigms | Chunk 11 | We build 3 now (atlas, disease_vs_healthy, trajectory) |
| Auto-execution of DAG steps | Future orchestrator | DAG suggests; the agent acts; researcher confirms |
| NLP extraction from free-text prompts | Future | `from_prompt()` is a stub — the skill handles extraction conversationally |
| Ontology lookups (UBERON, NCBI Taxonomy) | Future | Free-text tissue name is fine for now |

---

## 8. Open Questions

| # | Question | Options | Recommendation |
|---|----------|---------|----------------|
| 1 | **Should the DAG be editable by the researcher?** | Yes (add/remove/reorder steps), No (fixed per paradigm) | **Yes, but simple.** The researcher can skip steps or add custom steps. They can't reorder dependencies. |
| 2 | **Should init block until all critical fields are set?** | Hard block, soft warning, allow proceeding | **Soft block.** The agent asks, but if the researcher says "just analyze it", proceed with cell_atlas defaults and warn about assumptions. |
| 3 | **How to handle multi-modal data (GEX + protein + VDJ)?** | Full support now, detect + record only, ignore | **Detect + record only.** Note the modalities in `project.json` but only generate DAGs for GEX. Multi-modal pipelines are Chunk 11. |
| 4 | **Should the agent re-ask on every session if context is incomplete?** | Yes (every time), No (once is enough), Only if analysis-critical | **Only if analysis-critical fields are missing.** If paradigm is set but donor age isn't, don't nag. |
| 5 | **DAG file format?** | JSON, YAML, Python dataclass only | **JSON.** Consistent with project.json, readable, serializable. |

---

## 9. Implementation Checklist

```
Chunk 6: Experiment Context & DAG
├── [ ] 6.1  scagent/context.py — ExperimentContext class
│     ├── [ ] Schema as dataclass / typed dict
│     ├── [ ] validate() — check required fields per paradigm
│     ├── [ ] infer_from_data(adata) — species, platform, cell count
│     ├── [ ] is_complete() → bool
│     ├── [ ] needs_batch_correction() / needs_pseudobulk_de() / needs_trajectory()
│     ├── [ ] save() / load() — .scagent/project.json
│     └── [ ] inferred flags on auto-detected fields
├── [ ] 6.2  scagent/dag.py — AnalysisDAG class
│     ├── [ ] DAGStep dataclass
│     ├── [ ] from_paradigm(paradigm, context) — factory for 3 paradigms
│     ├── [ ] next_step() / complete_step() / skip_step()
│     ├── [ ] is_valid_step() — dependency checking
│     ├── [ ] summary() — markdown table
│     ├── [ ] save() / load() — .scagent/dag.json
│     └── [ ] Conditional step logic (batch correction, pseudobulk, trajectory)
├── [ ] 6.3  .pi/skills/init/SKILL.md — conversational init flow
│     ├── [ ] Infer-first, ask-only-what-you-must logic
│     ├── [ ] minSCe checklist (analysis-critical vs. metadata)
│     ├── [ ] Three UX scenarios (rich prompt, minimal, no context)
│     └── [ ] Proactive asking if no project.json exists
├── [ ] 6.4  .pi/skills/dag/SKILL.md — DAG-following behavior
│     ├── [ ] Suggest next step before each action
│     ├── [ ] Complete / skip step tracking
│     └── [ ] Show progress on request
├── [ ] 6.5  tests/test_context.py — 11 unit tests
├── [ ] 6.6  tests/test_dag.py — 11 unit tests
├── [ ] 6.7  tests/test_init_integration.py — real PBMC end-to-end
├── [ ] 6.8  Update scagent/__init__.py — export ExperimentContext, AnalysisDAG
└── [ ] 6.9  Verify: init on PBMC dataset → valid project.json + dag.json
```

**Estimated output:**
- `scagent/context.py` — ~250–350 lines
- `scagent/dag.py` — ~250–350 lines
- `.pi/skills/init/SKILL.md` — ~120 lines
- `.pi/skills/dag/SKILL.md` — ~80 lines
- `tests/test_context.py` — ~200 lines
- `tests/test_dag.py` — ~250 lines
- `tests/test_init_integration.py` — ~80 lines
