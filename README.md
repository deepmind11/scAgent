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

scAgent passes **6 out of 7** canonical tasks (85.7%) on [SC-Bench](https://github.com/latchbio/scbench) (Workman et al., 2026, [LatchBio](https://latch.bio)), a benchmark of 394 verifiable problems derived from practical scRNA-seq workflows. The current top baseline model on SC-Bench scores 52.8%.

| Task | Result |
|------|--------|
| QC (cell filtering) | ✅ Pass |
| Normalization | ✅ Pass |
| HVG / Feature Selection | ✅ Pass |
| Clustering | ✅ Pass |
| Cell Type Annotation | ✅ Pass |
| Differential Expression | ✅ Pass |
| Trajectory Analysis | ❌ Fail (tools not yet implemented) |

<details>
<summary>Reproduce the evaluation</summary>

The 7 canonical Chromium evaluations from SC-Bench are bundled in [`eval/evals_canonical_chromium/`](eval/evals_canonical_chromium/).

**Tool-level benchmark** (no LLM, tests analysis logic directly):
```bash
pip install -e ".[eval]"
python eval/run_benchmark.py
```

**Full agent benchmark** (LLM interprets task → chooses tools → produces answer):
```bash
pip install -e ".[eval]"
python eval/run_llm_benchmark.py                        # default: claude-opus-4-6
python eval/run_llm_benchmark.py --model claude-sonnet-4-5  # or any model
```

Both produce pass/fail summaries and save results to `eval/results/`.

The evaluations use [SC-Bench](https://github.com/latchbio/scbench) by [LatchBio](https://latch.bio) ([eval-graders](https://github.com/latchbio/eval-graders)). The canonical eval JSONs are included under Apache 2.0.

```bibtex
@article{scbench2026,
  title={scBench: Evaluating AI Agents on Single-Cell RNA-seq Analysis},
  author={Workman, Kenny and Yang, Zhen and Muralidharan, Harihara and Abdulali, Aidan and Le, Hannah},
  year={2026},
  note={LatchBio}
}
```

</details>

## Key Features

### 30+ Analysis Tools

Full pipeline from raw counts to publication: QC → normalization → HVG → PCA → batch integration (Harmony, scVI, BBKNN, Scanorama) → clustering (Leiden/Louvain) → cell type annotation (CellTypist) → differential expression (pseudobulk DESeq2/edgeR, Wilcoxon) → pathway enrichment (GSEA, ClusterProfiler) → cell communication (CellChat, CellPhoneDB). Trajectory analysis tools (Monocle3, Slingshot, scVelo, PAGA) have JSON schemas defined but are not yet implemented.

Each tool is defined by a JSON schema in [`tools/`](tools/) with parameter types, constraints, and literature-backed defaults.

### Paradigm-Aware Analysis DAG

Every experiment has a paradigm (cell atlas, disease vs. healthy, developmental trajectory, perturbation). The analysis DAG in [`scagent/dag.py`](scagent/dag.py) adapts valid step ordering based on the paradigm — preventing invalid operations like running pseudobulk DE on a single-condition atlas, or clustering on UMAP coordinates.

### State Management & Branching

You can run an analysis, then go back to any checkpoint and branch off in a different direction. [`scagent/state.py`](scagent/state.py) implements lazy-checkpointed branching — AnnData objects (200MB–2GB) live in memory during normal work and only get written to disk when you fork, switch branches, or end a session. The branch tree lives in `.scagent/branches/`.

### Long-Term Memory

Cross-session memory via [MemPalace](https://github.com/AidanCooper/mempalace) (ChromaDB-backed). Every analysis decision, parameter choice, and conversation is stored with metadata (branch, timestamp, analysis phase). When the agent needs to recall a past decision — "why did we choose resolution 0.8?" — it does a semantic search and retrieves the relevant context, even across sessions.

### Provenance & Reproducibility

Every tool invocation is recorded as a W3C PROV-O graph in JSON-LD — inputs, parameters, outputs, software versions, timestamps. Full traceability.

The export module ([`scagent/export.py`](scagent/export.py)) generates from the provenance chain:
- **Methods section** — camera-ready prose for a paper, auto-generated from provenance records
- **Reproducibility package** — `methods.md`, `params.json`, `README.md`, and `replay.py` (a script that re-runs the entire analysis from raw data using the recorded parameters)

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
│   ├── dag.py            #   Paradigm-aware analysis DAG
│   ├── context.py        #   Experiment context & metadata
│   ├── knowledge.py      #   Marker gene database
│   ├── inspector.py      #   AnnData state inspection
│   ├── dependencies.py   #   Prerequisite checking
│   ├── export.py         #   Methods section & repro package generation
│   └── tools/            #   Scanpy/analysis tool implementations
├── tools/                # Tool registry (~30 JSON schemas)
├── .pi/                  # Agent configuration
│   ├── SYSTEM.md         #   System prompt (identity + rules)
│   ├── settings.json     #   Model defaults
│   └── skills/           #   19 step-specific instruction sets
├── best_practices/       # Literature-backed reference guides
│   └── reference/        #   Per-step best practice summaries
├── eval/                 # SC-Bench evaluation framework & results
├── tests/                # Unit and integration tests
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

## License

[MIT](LICENSE)
