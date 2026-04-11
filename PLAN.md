# scAgent Build Plan

**Created:** 2026-04-11
**Status:** Active

---

## Prerequisites

### Prereq 1: Reference Dataset

**Chosen dataset: 10x Genomics 10k PBMCs, v3 chemistry**

| Field | Value |
|-------|-------|
| Source | https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0 |
| Cells | ~10,000 PBMCs from healthy donor |
| Chemistry | Chromium v3 (3' GEX) |
| Cell Ranger | 3.0.0 |
| Files needed | `filtered_feature_bc_matrix.h5` (~17 MB) |
| Paradigm | Cell atlas (single-tissue deep profiling) |

**Why this dataset:**
- The canonical Scanpy PBMC 3k is too small to stress-test anything. 10k is big enough to find real cell types but small enough to iterate fast.
- PBMC is the most well-characterized tissue in scRNA-seq — we have ground-truth expectations for what cell types should appear (CD4/CD8 T cells, B cells, NK, monocytes, DCs, platelets).
- This exact dataset is used across many tutorials (Scanpy, Seurat, Cell Ranger). We can cross-check our analysis against all of them.
- v3 chemistry is the current standard for Chromium.
- Filtered matrix means we skip Cell Ranger — exactly what you wanted.

**Additional datasets to add later (for different paradigms):**
- **Disease vs. healthy:** The scBench Chromium evals use a 4T1 tumor dataset. We should use that same one for eval compatibility.
- **Trajectory:** An organoid or hematopoiesis dataset with known lineage structure.

**Setup:**
```bash
mkdir -p ~/Projects/scAgent/data/pbmc10k
cd ~/Projects/scAgent/data/pbmc10k
# Download filtered_feature_bc_matrix.h5 from 10x website (requires email registration)
# OR: use scanpy to download programmatically later
```

### Prereq 2: Eval Harness

**Location: `~/Projects/scbench-eval/`** (separate from scAgent)

The eval system wraps [latchbio/scbench](https://github.com/latchbio/scbench). scBench provides:
- 7 canonical Chromium evals (QC, normalization, HVG, clustering, cell typing, DE, trajectory)
- `latch-eval-tools` Python package with graders and `EvalRunner`
- Agent interface: `agent_function(task_prompt, work_dir) -> answer_dict`
- Grader types: `numeric_tolerance`, `distribution_comparison`, `marker_gene_precision_recall`, `label_set_jaccard`, `multiple_choice`

We need an adapter that launches `scagent` (via `pi` in print mode or RPC mode) as the `agent_function`.

**Setup:**
```bash
mkdir -p ~/Projects/scbench-eval
cd ~/Projects/scbench-eval
git clone https://github.com/latchbio/scbench.git
cd scbench
pip install -e .

# Create our adapter
# scagent_adapter.py — bridges scbench's agent_function interface to pi --mode rpc
```

**Key design decision:** Use `pi --mode rpc` (stdin/stdout JSON-lines) to drive scAgent programmatically from the eval harness. This means the eval runner (Python) sends the task prompt to scAgent via RPC, scAgent runs the analysis, and writes `eval_answer.json` to the work_dir.

---

## Build Chunks

### Chunk 0: Environment & Scaffolding
**Goal:** Python environment, dependencies, data download, verify scBench runs with a dummy agent.

- [ ] Create `pyproject.toml` or `requirements.txt` for scAgent's Python deps (scanpy, anndata, scrublet, harmonypy, etc.)
- [ ] Set up a virtual environment or conda env
- [ ] Download PBMC 10k dataset to `data/pbmc10k/`
- [ ] Verify: `import scanpy as sc; adata = sc.read_10x_h5('data/pbmc10k/filtered_feature_bc_matrix.h5'); print(adata)` works
- [ ] Clone scBench, install it, run `scbench validate` on one canonical eval
- [ ] Run scBench with a trivial dummy agent to confirm the harness works end-to-end
- [ ] Decide: conda vs. uv vs. pip for dependency management

**Deliverable:** Working Python env, loadable dataset, scBench harness confirmed working.

---

### Chunk 1: Tool Registry (Core 10)
**Goal:** JSON schema files for the 10 core tools, covering QC through clustering.

Tools to register:
1. `load_10x_h5` — Load filtered_feature_bc_matrix.h5 into AnnData
2. `filter_cells` — QC filtering (min_genes, max_genes, max_mito)
3. `filter_genes` — Gene filtering (min_cells)
4. `scrublet_doublets` — Doublet detection
5. `normalize` — Log-normalize + scale (or SCTransform)
6. `highly_variable_genes` — Feature selection
7. `pca` — Dimensionality reduction
8. `harmony` — Batch correction
9. `neighbors` — Neighbor graph construction
10. `leiden` — Leiden clustering

Each JSON file follows the schema from the architecture doc (§4): tool_id, function, valid_after, valid_before, parameters with types/defaults/guidance, outputs, validation rules, provenance_captures.

**Deliverable:** `tools/` directory with 10 JSON files. Parseable and internally consistent.

---

### Chunk 2: Skills — QC & Preprocessing
**Goal:** scAgent skills that teach the LLM how to execute the first analysis steps.

Skills to create:
- `.pi/skills/load-data/SKILL.md` — How to load 10x data, check shape, report basic stats
- `.pi/skills/qc/SKILL.md` — QC workflow: compute metrics, show distributions, set thresholds, filter, report
- `.pi/skills/doublet-detection/SKILL.md` — Scrublet workflow with parameter guidance
- `.pi/skills/normalize/SKILL.md` — Normalization options, when to use which, execute

Each skill is a Markdown file that instructs the agent on:
1. What this step does biologically
2. What code to run (with Scanpy)
3. What parameters to use and why
4. What to check/validate after
5. What to show the user
6. What provenance to record

**Deliverable:** 4 skills in `.pi/skills/`. scAgent can load data, run QC, detect doublets, and normalize — with correct parameters and informative outputs — just by reading the skills.

---

### Chunk 3: Skills — Clustering & Annotation
**Goal:** Skills for HVG → PCA → neighbors → clustering → marker genes → cell type annotation.

Skills to create:
- `.pi/skills/feature-selection/SKILL.md` — HVG selection
- `.pi/skills/dimensionality-reduction/SKILL.md` — PCA + UMAP/tSNE
- `.pi/skills/clustering/SKILL.md` — Leiden clustering with resolution guidance, sweep strategy
- `.pi/skills/marker-genes/SKILL.md` — Differential expression for cluster markers (Wilcoxon)
- `.pi/skills/cell-annotation/SKILL.md` — Cell type annotation (manual + CellTypist)

**Deliverable:** 5 more skills. scAgent can now run a complete QC → annotation pipeline on the PBMC 10k dataset. This is the first end-to-end test.

**Milestone test:** Run scAgent on the PBMC 10k dataset interactively. Ask it to "analyze this dataset." Verify it produces sensible cell types (CD4 T, CD8 T, B, NK, Mono CD14, Mono FCGR3A, DC, Platelet).

---

### Chunk 4: Provenance Writer
**Goal:** Python module that emits PROV-JSONLD after each tool invocation.

- [ ] `scagent/provenance.py` — ProvenanceGraph class
  - `record_activity(tool_id, params, input_hash, output_hash, duration, user_prompt)`
  - `serialize() -> PROV-JSONLD dict`
  - `save(path)` / `load(path)`
  - `diff(branch_a, branch_b)` — compare two provenance chains
- [ ] Integrate with skills: each skill appends a provenance record after execution
- [ ] Provenance file stored at `.scagent/provenance.jsonld`

**Deliverable:** After running an analysis, `.scagent/provenance.jsonld` contains a complete, valid PROV-JSONLD graph. Each step is traceable.

---

### Chunk 5: State Manager (Branching)
**Goal:** AnnData snapshot save/load with content-addressed storage and branching.

- [ ] `scagent/state.py` — StateManager class
  - `save_snapshot(adata, step_name) -> hash`
  - `load_snapshot(hash) -> adata`
  - `create_branch(name, from_hash)`
  - `switch_branch(name)`
  - `list_branches()`
  - `current_state() -> hash, step_name, branch`
- [ ] Storage at `.scagent/branches/{name}/states/{hash}.h5ad`
- [ ] HEAD pointer file per branch
- [ ] Skill or REPL commands: `/branch create`, `/branch switch`, `/branch list`, `/branch compare`

**Deliverable:** Researcher can fork analysis at any point, try different parameters, switch back. State is persisted on disk.

---

### Chunk 6: Experiment Context & DAG
**Goal:** `scagent init` flow + paradigm-aware DAG generation.

- [ ] `scagent-project.json` schema (minSCe-based, from architecture §3)
- [ ] Skill: `.pi/skills/init/SKILL.md` — interactive project initialization
  - Ask: species, tissue, platform, chemistry, paradigm, samples
  - Generate `scagent-project.json`
  - Generate initial DAG based on paradigm
- [ ] DAG definitions for 3 paradigms:
  - `dags/cell_atlas.json`
  - `dags/disease_vs_healthy.json`
  - `dags/developmental_trajectory.json`
- [ ] DAG executor logic: track current position, suggest next step, validate ordering

**Deliverable:** `scagent init` creates a project with experiment context and a paradigm-appropriate analysis plan.

---

### Chunk 7: Eval Adapter & First Benchmark
**Goal:** Run scAgent against the 7 canonical scBench Chromium evals.

- [ ] `~/Projects/scbench-eval/scagent_adapter.py` — bridges scBench to scAgent via `pi --mode rpc` or `pi -p`
  - Receives `(task_prompt, work_dir)` from scBench
  - Copies the .h5ad into a temp scAgent project
  - Launches scAgent with the task prompt
  - Extracts `<EVAL_ANSWER>` from the response
  - Returns the answer dict
- [ ] Run all 7 Chromium canonical evals
- [ ] Record results: pass/fail per eval, agent output, time, cost
- [ ] Analyze failures: what went wrong, which skills/tools were missing

**Deliverable:** First scBench score. Baseline to improve from.

---

### Chunk 8: DE Pipeline (Pseudobulk)
**Goal:** Cross-condition differential expression using pseudobulk.

- [ ] Skills:
  - `.pi/skills/pseudobulk-de/SKILL.md` — aggregate to pseudobulk, run DESeq2/edgeR
  - `.pi/skills/pathway-enrichment/SKILL.md` — GSEA, clusterProfiler
- [ ] Tool registry entries for DESeq2 (via rpy2 or subprocess) and edgeR
- [ ] DAG enforcement: in `disease_vs_healthy` paradigm, DE step MUST use pseudobulk

This is where scBench's DE evals (41% accuracy for best model) should improve — by structurally preventing the cell-as-replicate error.

**Deliverable:** scAgent correctly runs pseudobulk DE on multi-condition data.

---

### Chunk 9: Memory Integration
**Goal:** Cross-session context persistence.

Two options:
- **MemPalace** — as designed in architecture. Requires checking if [github.com/milla-jovovich/mempalace](https://github.com/milla-jovovich/mempalace) actually exists and works.
- **Fallback** — Use pi's built-in memory package (`@samfp/pi-memory`) which is already installed in Feynman. Simpler, already works.

- [ ] Evaluate MemPalace: does it exist? Is it usable? Does it add value over pi-memory?
- [ ] If MemPalace works: integrate it for conversation history across sessions
- [ ] If not: use pi-memory for key facts + session summaries stored in `.scagent/memory/`
- [ ] Cross-session continuity: on startup, load last session summary + experiment context

**Deliverable:** Researcher can close scAgent and come back days later. The agent remembers what was done and why.

---

### Chunk 10: Knowledge Graph
**Goal:** Paper and marker gene database integration.

- [ ] Evaluate Graphify: does [github.com/safishamsi/graphify](https://github.com/safishamsi/graphify) exist and work?
- [ ] If Graphify works: integrate for paper processing
- [ ] If not: build a simpler marker database system
  - Load CellMarker 2.0, PanglaoDB, CellTypist models
  - Skill: `.pi/skills/knowledge/SKILL.md` — query marker databases, add papers
- [ ] `/knowledge add paper.pdf` command
- [ ] `/knowledge query "markers for CD8 effector T cells"` command

**Deliverable:** scAgent can look up marker genes from databases and papers. Cell annotation quality improves.

---

### Chunk 11: Full Paradigm Support
**Goal:** DAGs and specialized skills for remaining paradigms.

- [ ] Trajectory inference skills (Monocle 3, PAGA, scVelo)
- [ ] Perturbation screen skills (guide assignment, perturbation effects)
- [ ] V(D)J / immune repertoire skills
- [ ] CITE-seq / protein + RNA skills (WNN)
- [ ] Composition analysis skills (scCODA, miloR)
- [ ] Cell communication skills (CellChat, CellPhoneDB)

**Deliverable:** scAgent handles all 9 paradigms from the architecture.

---

### Chunk 12: Polish & Full Eval
**Goal:** Context management, export, and full scBench validation.

- [ ] Context compaction strategy (Tier 2 overflow at 150K → summarize)
- [ ] Export: provenance PDF, methods section generator, reproducibility package
- [ ] Request full 394-eval scBench access from kenny@latch.bio
- [ ] Run full eval, analyze per-task and per-platform breakdowns
- [ ] Iterate on skills/tools based on failure analysis
- [ ] Documentation and onboarding flow

**Deliverable:** Production-ready scAgent with validated scBench score.

---

## Dependency Graph

```
Prereq 1 (Dataset) ──┐
                      ├→ Chunk 0 (Environment) → Chunk 1 (Tool Registry) → Chunk 2 (Skills: QC)
Prereq 2 (Eval) ─────┘                                                          │
                                                                                 ▼
                                                                        Chunk 3 (Skills: Clustering)
                                                                                 │
                                                              ┌──────────────────┼──────────────────┐
                                                              ▼                  ▼                  ▼
                                                     Chunk 4 (Provenance)  Chunk 5 (Branching)  Chunk 6 (DAG/Init)
                                                              │                  │                  │
                                                              └──────────────────┴──────────────────┘
                                                                                 │
                                                                                 ▼
                                                                     Chunk 7 (First Benchmark)
                                                                                 │
                                                              ┌──────────────────┼──────────────────┐
                                                              ▼                  ▼                  ▼
                                                     Chunk 8 (DE Pipeline)  Chunk 9 (Memory)  Chunk 10 (Knowledge)
                                                              │                  │                  │
                                                              └──────────────────┴──────────────────┘
                                                                                 │
                                                                                 ▼
                                                                     Chunk 11 (Full Paradigms)
                                                                                 │
                                                                                 ▼
                                                                     Chunk 12 (Polish & Eval)
```

## Immediate Next Steps

1. **Set up Python environment** (Chunk 0)
2. **Download PBMC 10k dataset** (Prereq 1)
3. **Clone scBench and verify it runs** (Prereq 2)
4. **Write the first 10 tool registry JSONs** (Chunk 1)
5. **Write the QC skill** (Chunk 2) — this is where scAgent starts doing real work
