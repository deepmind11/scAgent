# Best Practice Reference: Marker Gene Detection

> Distilled from Heumos et al. 2023 and sc-best-practices.org.

## Method

- **Wilcoxon rank-sum test** performed best for identifying cluster-specific markers (Heumos et al., confirmed by Pullin & McCarthy 2022).
- Apply the test to compare each cluster against all other cells.

## Important Caveats

- **P-values from cluster markers are inflated** because the same data was used to define the clusters being tested (Zhang et al. 2019, cited in Heumos et al. and sc-best-practices.org). This is expected — the goal is ranking, not formal hypothesis testing.
- These markers are for **cluster characterization** (what defines cluster X?), NOT for cross-condition DE (disease vs. healthy). The two are fundamentally different analyses.
- Filter markers for specificity: require minimum fraction expressed in-group (e.g., 20%) and maximum fraction in out-group (e.g., 20%) to get cluster-specific genes (sc-best-practices.org).

## Interpreting Markers

- Markers are often sparsely expressed — "only a subset of cells of a cell type in which a marker was detected" due to dropout (sc-best-practices.org). This is why we cluster first and then annotate clusters, not individual cells.
- If a cluster has no strong markers (<3 significant genes), it may be over-clustered — suggest merging or lowering resolution.
- If multiple clusters share the same top markers, they may be fragments of the same population.
- Ribosomal or mitochondrial genes as top markers usually indicate technical artifacts.

## Sources

- Heumos et al. 2023: https://doi.org/10.1038/s41576-023-00586-w
- sc-best-practices.org, Annotation chapter: https://www.sc-best-practices.org/cellular_structure/annotation.html
- Pullin & McCarthy, "A comparison of marker gene selection methods." bioRxiv (2022).
