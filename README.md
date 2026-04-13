# scAgent

An agentic single-cell RNA-seq analysis system for 10x Genomics Chromium data.

scAgent is an AI assistant that helps wet-lab biologists perform, understand, and interpret scRNA-seq experiments through natural language conversation. It enforces best practices, tracks full provenance, and produces reproducible analyses.

## What it does

- **Full pipeline support:** QC → normalization → HVG selection → PCA → batch integration → clustering → cell type annotation → differential expression → pathway enrichment
- **Paradigm-aware:** Adapts the analysis DAG to your experimental design (cell atlas, disease vs. healthy, trajectory, perturbation)
- **Provenance tracking:** Every step is recorded with parameters, tool versions, inputs, and outputs (PROV-JSONLD)
- **Branched exploration:** Fork at any step, compare different parameter choices, merge back
- **Knowledge-backed:** Built-in marker gene database + CellTypist integration for cell type annotation
- **Biologist-friendly:** Explains what it's doing and why — no programming knowledge required

## Prerequisites

- **Python ≥ 3.11**
- **Node.js ≥ 20.19** — [install](https://nodejs.org) (the agent runtime needs it)
- A **Claude subscription** (Pro/Max) — no API key needed, authenticate via OAuth

## Installation

```bash
git clone https://github.com/deepmind11/scAgent.git
cd scAgent
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

That's it. The `scagent` command auto-installs the agent runtime ([pi-coding-agent](https://www.npmjs.com/package/@mariozechner/pi-coding-agent)) on first run if Node.js is available.

## Quick start

```bash
# Launch scAgent (works from anywhere)
scagent
```

On first launch, authenticate with your Claude subscription:
```
> /login
# Opens your browser → sign in with your Claude account → done
```

This is a one-time step. Your credentials are stored locally in `~/.pi/agent/auth.json` and auto-refresh.

```bash
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

## Project structure

```
scAgent/
├── scagent/              # Python package — core logic
│   ├── context.py        #   Experiment context (minSCe metadata)
│   ├── dag.py            #   Analysis DAG (step ordering + dependencies)
│   ├── dependencies.py   #   Prerequisite checking
│   ├── inspector.py      #   AnnData state inspection
│   ├── knowledge.py      #   Marker gene database
│   ├── memory.py         #   Cross-session project memory
│   ├── provenance.py     #   PROV-JSONLD provenance tracking
│   ├── state.py          #   Branch & snapshot management
│   ├── export.py         #   Methods section & repro package generation
│   └── tools/            #   Scanpy/analysis tool implementations
├── tools/                # Tool registry (JSON schemas)
├── .pi/                  # Agent configuration
│   ├── SYSTEM.md         #   System prompt (scAgent identity + rules)
│   ├── settings.json     #   Model defaults
│   └── skills/           #   Step-specific agent skills (19 skills)
├── best_practices/       # Literature-backed reference guides
│   └── reference/        #   Per-step best practice summaries
├── eval/                 # Evaluation & benchmarking framework
├── tests/                # Unit and integration tests
├── scagent/cli.py        # CLI entry point (`scagent` command)
└── pyproject.toml        # Package metadata & dependencies
```

## How it works

scAgent is built on [pi-coding-agent](https://www.npmjs.com/package/@mariozechner/pi-coding-agent), an open-source AI coding agent CLI. The `.pi/` directory contains:

- **System prompt** (`SYSTEM.md`) — defines scAgent's identity, rules, and constraints
- **Skills** — 19 step-specific instruction sets (QC, clustering, annotation, etc.) that the agent loads contextually
- **Tool schemas** — JSON definitions for every analysis tool with parameter constraints

The Python package (`scagent/`) provides the backend logic: provenance tracking, state management, dependency resolution, and the actual Scanpy-based analysis tools.

When you chat with scAgent, it:
1. Identifies what you're asking for
2. Checks prerequisites against the analysis DAG
3. Loads the relevant skill for detailed instructions
4. Executes the analysis step with validated parameters
5. Records provenance and presents results with interpretation

## Running tests

```bash
pip install -e .
pytest tests/ -v
```

Note: Integration tests require a dataset at `data/pbmc10k/filtered_feature_bc_matrix.h5`. Download from [10x Genomics](https://www.10xgenomics.com/datasets/10-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0).

## Architecture

See [outputs/architecture.md](outputs/architecture.md) for the full system design document.

## License

[MIT](LICENSE)
