# scAgent

An agentic single-cell RNA-seq analysis system for 10x Genomics Chromium data.

scAgent is an AI assistant that helps wet-lab biologists perform, understand, and interpret scRNA-seq experiments through natural language conversation. It enforces best practices, tracks full provenance, and produces reproducible analyses — no programming knowledge required.

## Prerequisites

- **Python ≥ 3.11**
- **[Feynman](https://github.com/getcompanion-ai/feynman)** — the AI agent runtime
  ```bash
  curl -fsSL https://feynman.is/install | bash
  feynman setup   # authenticate with your Claude subscription
  ```
- *Optional:* **[Claude Code](https://docs.anthropic.com/en/docs/claude-code)** — also works as the agent runtime (scAgent's `.pi/` config is compatible with any pi-based agent)

## Installation

```bash
git clone https://github.com/deepmind11/scAgent.git
cd scAgent
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Quick Start

```bash
# Launch scAgent (works from anywhere)
scagent

# With specific model/thinking settings
scagent --model opus --thinking max

# Continue a previous session
scagent -c

# Resume/pick a session
scagent -r
```

On first launch, scAgent will ask you about your experiment:
1. **Load your data** — point it at a `filtered_feature_bc_matrix.h5` from Cell Ranger
2. **Describe your experiment** — tissue, organism, experimental design
3. **Start analyzing** — ask questions in natural language

Example prompts:
```
> Load data/pbmc10k/filtered_feature_bc_matrix.h5 and run QC
> What cell types are in my data?
> Compare gene expression between clusters 0 and 3
> Run pseudobulk DE between treatment and control
> Generate a methods section for my paper
```

## Evaluation: SC-Bench

scAgent is evaluated on [SC-Bench](https://github.com/latchbio/scbench) (Workman et al., 2026, [LatchBio](https://latch.bio)), a benchmark of 394 verifiable problems derived from practical scRNA-seq workflows. The current top baseline model on SC-Bench scores 52.8%. The 7 canonical Chromium evaluations are bundled in [`eval/evals_canonical_chromium/`](eval/evals_canonical_chromium/).

| Task | Result |
|------|--------|
| Normalization | ✅ Pass |
| HVG / Feature Selection | ✅ Pass |
| Clustering | ✅ Pass |
| Cell Type Annotation | ✅ Pass |
| Differential Expression | ✅ Pass |
| QC (cell filtering) | ❌ Fail |
| Trajectory Analysis | ❌ Fail |

**QC** fails because the agent applies standard textbook cutoffs (max 5,000 genes) to an already-cleaned dataset — it needs to inspect distributions before filtering and recognize when data is pre-filtered. **Trajectory** achieves 0.5 recall on terminal marker recovery — the agent identifies the correct CAF differentiation axis (universal fibroblast → myCAF) and recovers `Ly6c1` but misses `Acta2`, a canonical myCAF marker. Improving pseudotime terminal group detection and marker ranking is next.

The eval runs the full LLM agent end-to-end: the agent receives a task prompt, reasons about what analysis to perform, calls tools, and produces a structured answer that is graded automatically.

```bash
pip install -e ".[eval]"
python eval/run_llm_benchmark.py                        # default: claude-opus-4-6
python eval/run_llm_benchmark.py --model claude-sonnet-4-5  # or any model
```

Results are saved to `eval/results/`.

The evaluations use [SC-Bench](https://github.com/latchbio/scbench) by [LatchBio](https://latch.bio) ([eval-graders](https://github.com/latchbio/eval-graders)). The canonical eval JSONs are included under Apache 2.0.

```bibtex
@article{scbench2026,
  title={scBench: Evaluating AI Agents on Single-Cell RNA-seq Analysis},
  author={Workman, Kenny and Yang, Zhen and Muralidharan, Harihara and Abdulali, Aidan and Le, Hannah},
  year={2026},
  note={LatchBio}
}
```

## Key Features

### 34 Analysis Tools Across 7 Paradigms

Full pipeline from raw counts to publication: QC → normalization → HVG → PCA → batch integration (Harmony, scVI, BBKNN, Scanorama) → clustering (Leiden/Louvain) → cell type annotation (CellTypist) → differential expression (pseudobulk DESeq2/edgeR, Wilcoxon) → pathway enrichment (GSEA, ClusterProfiler) → trajectory inference (PAGA + DPT + scVelo) → compositional analysis (scCODA + Milo) → cell communication (LIANA+) → perturbation analysis (guide assignment + DE) → immune repertoire (Scirpy: clonotype, diversity) → multimodal CITE-seq (CLR + WNN).

Each tool is defined by a JSON schema in [`tools/`](tools/) with parameter types, constraints, and literature-backed defaults. 21 Python tool wrappers in [`scagent/tools/`](scagent/tools/) implement the analysis logic with input validation, guard rails, plotting, and structured provenance output. Default parameters and analysis guidelines are derived from [Best Practices for Single Cell Analysis across Modalities](https://www.nature.com/articles/s41576-023-00586-w) (Heumos et al., 2023), the [sc-best-practices.org online book](https://www.sc-best-practices.org/preamble.html) (Theis Lab), and the [10x Genomics Analysis Guide](https://www.10xgenomics.com/analysis-guides/best-practices-analysis-10x-single-cell-rnaseq-data). Experiment metadata collection follows the [minSCe guidelines](https://doi.org/10.1038/s41587-020-00744-z) (Füllgrabe et al., 2020). Per-step reference summaries are in [`best_practices/reference/`](best_practices/reference/).

### State-Aware Data Inspector

When you load data — whether raw from Cell Ranger or a half-processed `.h5ad` from a collaborator — the inspector ([`scagent/inspector.py`](scagent/inspector.py)) automatically determines what has already been done: is `adata.X` raw counts, log-normalized, or scaled? Are there PCA embeddings? Clustering labels? Batch columns? The agent uses this to reason about what steps are needed rather than assuming it controls the data from the start.

### Dependency Resolution

The dependency module ([`scagent/dependencies.py`](scagent/dependencies.py)) encodes what each analysis step requires — both prerequisite steps and data conditions. If you ask for differential expression but haven't clustered yet, scAgent can check what's missing (`check_prerequisites`), plan the minimal steps to get there (`plan_steps`), or auto-run the prerequisites (`ensure_ready_for`).

### Paradigm-Aware Analysis DAG

Every experiment has a paradigm — one of 7 supported analysis types. The analysis DAG in [`scagent/dag.py`](scagent/dag.py) generates a paradigm-specific step ordering with validated dependencies:

| Paradigm | Key Steps |
|---|---|
| `cell_atlas` | Standard QC → clustering → annotation |
| `disease_vs_healthy` | + pseudobulk DE + composition + enrichment |
| `developmental_trajectory` | + PAGA topology → DPT pseudotime → scVelo |
| `perturbation_screen` | + guide assignment → perturbation DE |
| `temporal_longitudinal` | + mandatory batch correction + time-course DE |
| `immune_repertoire` | + VDJ loading → clonotype → diversity |
| `multimodal` | + protein loading → CLR → WNN joint graph |

The DAG prevents invalid operations (pseudobulk DE on single-condition data, clustering on UMAP coordinates) and tracks progress through each step.

### State Management & Branching

You can run an analysis, then go back to any checkpoint and branch off in a different direction. [`scagent/state.py`](scagent/state.py) implements lazy-checkpointed branching — AnnData objects (200MB–2GB) live in memory during normal work and only get written to disk when you fork, switch branches, or end a session. The branch tree lives in `.scagent/branches/`.

### Long-Term Memory

Real-world scRNA-seq analysis happens over weeks — you run QC on Monday, come back to clustering on Thursday, and a reviewer asks about your normalization choice a month later. scAgent maintains cross-session memory via [MemPalace](https://github.com/AidanCooper/mempalace) (ChromaDB-backed) so nothing is lost between sessions. Every analysis decision, parameter choice, and conversation is stored with metadata (branch, timestamp, analysis phase). When the agent needs to recall a past decision — "why did we choose resolution 0.8?" — it does a semantic search and retrieves the relevant context, even from weeks ago.

### Provenance & Reproducibility

Every tool invocation is recorded as a W3C PROV-O graph in JSON-LD — inputs, parameters, outputs, software versions, timestamps. Full traceability.

The export module ([`scagent/export.py`](scagent/export.py)) generates from the provenance chain:
- **Methods section** — camera-ready prose for a paper, auto-generated from provenance records
- **Reproducibility package** — a self-contained directory with `methods.md`, `params.json`, `README.md`, and `replay.py` (a script that re-runs the entire analysis from raw data using the recorded parameters)

Share the reproducibility package with a collaborator or reviewer and they can reproduce your exact analysis independently — same parameters, same tool versions, same results.

### Biologist-Friendly

scAgent explains what it's doing and why at every step. It presents QC distributions and plots, shows parameter choices with rationale, displays top marker genes per cluster with supporting evidence, and asks for confirmation before advancing. No programming knowledge is assumed.

## Project Structure

```
scAgent/
├── scagent/              # Core Python package
│   ├── cli.py            #   CLI entry point (scagent command)
│   ├── state.py          #   Branch & snapshot management
│   ├── memory.py         #   Long-term memory (MemPalace/ChromaDB)
│   ├── provenance.py     #   W3C PROV-O provenance tracking (JSON-LD)
│   ├── dag.py            #   Paradigm-aware analysis DAG (7 paradigms)
│   ├── context.py        #   Experiment context & metadata
│   ├── knowledge.py      #   Marker gene database
│   ├── inspector.py      #   AnnData state inspection
│   ├── dependencies.py   #   Prerequisite checking
│   ├── export.py         #   Methods section & repro package generation
│   └── tools/            #   21 analysis tool implementations
│       ├── trajectory.py     # PAGA + DPT + scVelo
│       ├── composition.py    # scCODA + Milo (via pertpy)
│       ├── communication.py  # LIANA+ consensus L-R
│       ├── perturbation.py   # Guide assignment + perturbation DE
│       ├── repertoire.py     # Scirpy: VDJ + clonotype analysis
│       ├── multimodal.py     # CITE-seq: CLR + WNN
│       └── ...               # QC, clustering, DE, annotation, etc.
├── tools/                # Tool registry (34 JSON schemas)
├── .pi/                  # Agent configuration
│   ├── SYSTEM.md         #   System prompt (identity + rules)
│   ├── settings.json     #   Model defaults
│   └── skills/           #   26 step-specific instruction sets
├── best_practices/       # Literature-backed reference guides
│   ├── sc_best_practices.pdf        # Heumos et al. 2023
│   ├── best_practice_10xGenomics_scRNAseq.pdf  # 10x official guide
│   └── reference/        #   Per-step best practice summaries
├── eval/                 # SC-Bench evaluation framework & results
├── tests/                # 25 test files, 128+ tests
└── pyproject.toml        # Package metadata & dependencies
```

## How It Works

scAgent is built on [Feynman](https://github.com/getcompanion-ai/feynman), an open-source AI research agent. The `.pi/` directory configures the agent:

- **System prompt** ([`SYSTEM.md`](.pi/SYSTEM.md)) — defines scAgent's identity, rules (never cluster on UMAP, always pseudobulk for cross-condition DE, etc.), and interaction style
- **Skills** — 19 step-specific instruction sets loaded contextually when the agent reaches each analysis phase
- **Tool schemas** — JSON definitions for every tool with parameter constraints and defaults

When you chat with scAgent, it:
1. Identifies what you're asking for
2. Checks prerequisites against the analysis DAG
3. Loads the relevant skill for detailed instructions
4. Executes the analysis step with validated parameters
5. Records provenance and presents results with interpretation

## Running Tests

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

Integration tests require a dataset at `data/pbmc10k/filtered_feature_bc_matrix.h5`. Download from [10x Genomics](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0).

## Architecture

See [outputs/architecture.md](outputs/architecture.md) for the full system design document.

## Acknowledgments

- [SC-Bench](https://github.com/latchbio/scbench) by [LatchBio](https://latch.bio) — evaluation framework
- [Feynman](https://github.com/getcompanion-ai/feynman) — agent runtime
- [Scanpy](https://scanpy.readthedocs.io/) — core analysis engine
- [CellTypist](https://www.celltypist.org/) — cell type annotation
- [minSCe guidelines](https://doi.org/10.1038/s41587-020-00744-z) (Füllgrabe et al., 2020) — experiment metadata standards for scRNA-seq reporting

## Coming Soon

- **Full SC-Bench evaluation** — run against all 394 tasks (currently limited to 7 canonical Chromium evals)
- **Programmatic agent** — replace the subprocess-based LLM runner with a direct Anthropic API agent loop (`scagent/agent.py`)
- **scGPT / foundation model embeddings** — alternative to PCA for annotation transfer (tool schema defined, awaiting GPU support)

## License

[MIT](LICENSE)
