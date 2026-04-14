---
name: trajectory
description: Run trajectory inference using PAGA + Palantir + (optional) scVelo/CellRank. Use when the paradigm is developmental_trajectory or when cells represent a continuous biological process.
---

# Skill: Trajectory Inference

Infer developmental trajectories following the sc-best-practices.org recommended workflow:
**PAGA → Palantir → (optional) scVelo → CellRank**

## When to Use

- **Paradigm is `developmental_trajectory`** — the DAG includes these steps.
- **User asks about differentiation, development, lineages, pseudotime, or cell fate.**
- **Data contains cells from a continuous process** (e.g., hematopoiesis, neurogenesis).

## When NOT to Use

- Static atlas without developmental context → `cell_atlas` paradigm.
- Disease vs. healthy comparison → `disease_vs_healthy` paradigm.
- Cells are all from the same terminal state → trajectory has no biological meaning.

## Workflow

### Step 1: PAGA (Topology)

PAGA finds the coarse-grained connectivity between clusters.

```python
from scagent.tools.trajectory import run_paga
result = run_paga(adata, groups="leiden", plot_dir="plots/trajectory")
```

### Step 2: Palantir (Pseudotime + Terminal States) — RECOMMENDED

Palantir models trajectories as Markov chains, producing smooth pseudotime and fate probabilities toward terminal states. **Preferred over DPT.** [BP-2 Ch. 14]

```python
from scagent.tools.trajectory import run_palantir
result = run_palantir(adata, root_cell_type="HSC", plot_dir="plots/trajectory")
# Outputs: palantir_pseudotime, terminal states, fate probabilities
```

**Why Palantir over DPT:** The sc-best-practices.org tutorial directly compares them on bone marrow data. DPT inflates pseudotime for disconnected lineages (CLPs get extreme values). Palantir increases continuously with developmental maturity. Their conclusion: "we would conclude to continue working with the Palantir pseudotime."

### Step 2 (alt): DPT (Fallback)

Use DPT only if Palantir is unavailable or for simple linear trajectories.

```python
from scagent.tools.trajectory import run_diffusion_pseudotime
result = run_diffusion_pseudotime(adata, root_cell_type="HSC")
```

### Step 3: scVelo (RNA Velocity — optional)

Adds directionality using spliced/unspliced RNA ratios. Only if data has the layers.

```python
from scagent.tools.trajectory import run_scvelo
result = run_scvelo(adata, mode="dynamical", plot_dir="plots/trajectory")
```

**Prerequisite:** `adata.layers["spliced"]` and `adata.layers["unspliced"]` must exist.
Standard Cell Ranger output does NOT include these — needs velocyto or STARsolo.

### Step 4: CellRank (Fate Mapping — recommended after velocity)

CellRank uses velocity to build a transition matrix and identify terminal/initial states with fate probabilities. **Recommended over interpreting raw velocity stream plots.** [BP-2 Ch. 15]

```python
from scagent.tools.trajectory import run_cellrank
result = run_cellrank(adata, plot_dir="plots/trajectory")
```

## Guard Rails

1. **PAGA first** for topology validation.
2. **Palantir is the default pseudotime method.** Use DPT only as fallback.
3. **Root cell type is required.** Prompt the user if not specified.
4. **Validate scVelo phase portraits.** Almond shape expected. If absent, velocity is unreliable. [BP-1]
5. **Do NOT over-interpret 2D velocity streams.** Use CellRank for rigorous fate analysis. [BP-2 Ch. 15]
6. **scVelo on steady-state populations is misleading.** Warn the user. [BP-1]

## Interpretation Guidance

- **Pseudotime is relative, not absolute.** Orders cells but doesn't measure real time.
- **PAGA edges ≠ lineage commitment.** They show transcriptomic similarity.
- **RNA velocity ≠ speed.** Direction of gene expression change, not rate.
- **Fate probabilities from CellRank are more reliable than stream plots.** [BP-2 Ch. 15]
- **Multiple methods should agree** on major branches to consider them robust.

## Best-Practice References

- [BP-2] sc-best-practices.org Ch. 14 — Palantir preferred over DPT
- [BP-2] sc-best-practices.org Ch. 15 — scVelo + CellRank for velocity/fate
- [BP-1] Heumos et al. 2023, §"From discrete states to continuous processes"
