# Best Practice Reference: Trajectory Inference

> Distilled from Heumos et al. 2023.

## When to Use

- For non-stationary biological processes: differentiation, development, reprogramming.
- When cells traverse a continuous space of states rather than discrete clusters.

## Method Selection

- Trajectories can be linear, cyclic, tree-like, or graph-shaped.
- **Slingshot** performed better for simple topologies (Saelens et al. 2019 benchmark).
- **PAGA** and **RaceID/StemID** scored better for complex trajectories.
- Use **dynguidelines** (dynverse) to select an applicable method based on expected topology.
- "When the expected topology is unknown, trajectories and downstream hypotheses should be confirmed by multiple trajectory inference methods using different underlying assumptions" (Heumos et al.).
- "Inferred trajectories might not necessarily have biological meaning" (Heumos et al.) — always seek independent validation.

## RNA Velocity

- **scVelo** models splicing kinetics (unspliced → spliced) to infer direction of state changes.
- **Check phase portraits** of high-likelihood genes: they should show an almond-shaped induction/repression pattern. If not, velocity may be incorrectly inferred.
- Limitations: assumes gene independence and constant rates. Fails with transcriptional bursts, steady-state populations, or processes outside the splicing kinetics timescale.
- **CellRank** uses velocity fields to estimate cellular fates — the recommended downstream tool.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- Saelens et al. "A comparison of single-cell trajectory inference methods." Nat Biotechnol 37, 547–554 (2019).
