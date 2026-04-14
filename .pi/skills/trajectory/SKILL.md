---
name: trajectory
description: Run trajectory inference and pseudotime analysis using PAGA + DPT (+ optional scVelo). Use when the paradigm is developmental_trajectory or when cells represent a continuous biological process (differentiation, development, activation).
---

# Skill: Trajectory Inference

Infer developmental trajectories using the Scanpy-native workflow: **PAGA → DPT → (optional) scVelo**.

## When to Use

- **Paradigm is `developmental_trajectory`** — the DAG includes this step.
- **User asks about differentiation, development, lineages, or pseudotime.**
- **Data contains cells from a continuous process** (e.g., hematopoiesis, neurogenesis, T cell activation).

## When NOT to Use

- Static atlas without developmental context → use `cell_atlas` paradigm.
- Disease vs. healthy comparison → use `disease_vs_healthy` paradigm (composition + DE, not trajectory).
- Cells are all from the same terminal state → trajectory has no biological meaning.

## Workflow

### Step 1: PAGA (Topology)

PAGA finds the coarse-grained connectivity between clusters — which cell types connect to which.

```python
from scagent.tools.trajectory import run_paga

result = run_paga(adata, groups="leiden", plot_dir="plots/trajectory")
# Check: is the graph connected? Are edges biologically plausible?
```

**Dynverse benchmark evidence:** PAGA/PAGA_Tree beats Slingshot on complex topologies (tree: 0.770 vs 0.630, multifurcation: 0.775 vs 0.719). Slingshot only wins on simpler linear/cycle topologies.

### Step 2: DPT (Pseudotime ordering)

DPT orders cells along the topology inferred by PAGA. Requires specifying a **root cell type** — the starting point of the developmental process.

```python
from scagent.tools.trajectory import run_diffusion_pseudotime

result = run_diffusion_pseudotime(
    adata,
    root_cell_type="HSC",          # biological starting point
    cell_type_key="cell_type",
    n_dcs=15,
    plot_dir="plots/trajectory",
)
# Outputs: adata.obs["dpt_pseudotime"], cell ordering, violin plots
```

**Choosing the root cell type:**
- Ask the user. Common roots: HSCs (blood), neural progenitors (brain), basal cells (epithelium).
- If unsure, inspect diffusion components: the most extreme cell in the component with highest variance within stem/progenitor clusters is a good candidate.
- Pseudotime is RELATIVE — it has no absolute timescale.

### Step 3: scVelo (RNA Velocity — optional)

RNA velocity adds DIRECTIONALITY to the trajectory using spliced/unspliced RNA ratios. Only available if the data was processed with velocyto or STARsolo.

```python
from scagent.tools.trajectory import run_scvelo

result = run_scvelo(adata, mode="stochastic", plot_dir="plots/trajectory")
# Outputs: velocity stream/arrow plots, confidence scores
```

**Prerequisite check:** The agent MUST verify that `adata.layers["spliced"]` and `adata.layers["unspliced"]` exist before attempting scVelo. If absent, explain to the user that standard Cell Ranger output doesn't include these — they need velocyto or STARsolo reprocessing.

## Guard Rails

1. **PAGA before DPT.** Always run PAGA first to validate the topology, then DPT for ordering.
2. **Root cell type is required for DPT.** Prompt the user if not specified.
3. **DPT can inflate pseudotime for disconnected lineages.** If some cells have infinite pseudotime, warn the user and suggest subsetting to the lineage of interest. [BP-2 Ch. 14]
4. **Compare pseudotime methods if multiple lineages exist.** DPT may not capture branching well — suggest Palantir as a complement. [BP-2 Ch. 14]
5. **Validate scVelo phase portraits.** If the expected almond shape is absent in top-likelihood genes, velocity inference may be incorrect. [BP-1]
6. **scVelo on steady-state populations is misleading.** Erroneous velocity directions between independent terminal populations. Warn the user. [BP-1]

## Interpretation Guidance

- **Pseudotime is relative, not absolute.** It orders cells but doesn't measure real time.
- **PAGA edges ≠ lineage commitment.** They show transcriptomic similarity, not fate decisions.
- **RNA velocity ≠ speed.** It shows direction of gene expression change, not the rate of differentiation.
- **Trajectory inference on static data is an estimation.** It assumes cells captured at different stages represent a continuous process. This may not be true (e.g., cycling cells can create circular artifacts).
- **Multiple trajectory methods should agree on major branches.** If PAGA and DPT disagree with velocity, the trajectory may not be robust.

## Best-Practice References

- [BP-1] Heumos et al. 2023, §"From discrete states to continuous processes" (pp. 553-554)
- [BP-2] sc-best-practices.org Ch. 14 — Pseudotemporal ordering
- Saelens et al. 2019, Nat Biotechnol — dynverse benchmark (45 methods × 339 datasets)
