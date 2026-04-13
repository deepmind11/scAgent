---
name: export
description: Generate paper-ready methods text and reproducibility packages from provenance. Use when the user asks to write up methods, share analysis, or create a reproducibility package or replay script.
---

# Skill: Export

Generate paper-ready methods text and reproducibility packages from provenance.

## When to use

- User says "write up my methods" or "generate methods section"
- User wants to share analysis with a colleague
- User asks for a reproducibility package or replay script

## Methods section

```python
from scagent.export import generate_methods

text = generate_methods(provenance_graph, experiment_context)
```

Reads the provenance chain and writes natural English prose suitable for a
paper's Methods section.  Includes software versions, filtering thresholds,
clustering parameters, etc.

## Reproducibility package

```python
from scagent.export import generate_repro_package

generate_repro_package(provenance_graph, experiment_context, out_dir="repro/")
```

Creates a directory with:

| File | What |
|------|------|
| `replay.py` | Python script replaying every tool step |
| `params.json` | All parameters per step |
| `environment.json` | Software versions + platform |
| `provenance.json` | Full W3C PROV-O graph |
| `methods.md` | Generated methods text |
| `context.json` | Experiment metadata |
| `README.md` | Package description |

## What branch to export

By default exports the **promoted branch** (the canonical pipeline).  If no
branch is promoted, falls back to ``main``.  The user can specify a branch:

```python
generate_methods(graph, ctx, branch="high_res")
```

## After export

Tell the user:
- Where the files are
- That `replay.py` can be run standalone with `pip install scagent && python replay.py`
- That `methods.md` is ready to paste into a paper
