---
name: cell-communication
description: Infer cell-cell communication via ligand-receptor analysis using LIANA+. Use when analyzing tumor microenvironment, immune interactions, or any multi-cell-type system where signaling is of interest.
---

# Skill: Cell-Cell Communication

Infer ligand-receptor interactions between cell types using LIANA+ (consensus meta-method).

## When to Use

- **Multi-cell-type systems** — tumor microenvironment, immune niches, developing tissues.
- **User asks about cell signaling, ligand-receptor interactions, or intercellular communication.**
- **After cell type annotation** — requires labeled cell types.

## How to Use

```python
from scagent.tools.communication import run_liana

result = run_liana(
    adata,
    cell_type_key="cell_type",
    resource_name="consensus",  # LIANA's curated L-R database
    n_perms=1000,
    top_n=50,
    plot_dir="plots/communication",
)

# Top interactions
for interaction in result["interactions"][:10]:
    print(f"{interaction['source']} → {interaction['target']}: "
          f"{interaction['ligand_complex']} — {interaction['receptor_complex']}")
```

## Why LIANA+ Over CellChat?

- **LIANA+ wraps 8 methods** (CellPhoneDB, NATMI, Connectome, etc.) + consensus ranking.
- **Python-native**, scverse-maintained. CellChat is R-only with a single scoring function.
- **Multiple databases.** CellChat uses only its own; LIANA provides consensus, CellPhoneDB, CellChatDB, etc.
- **[BP-1]:** "Owing to the lack of consensus between tools, we recommend using LIANA."

## Key Caveats

1. **L-R databases are biased** toward specific pathways, functional categories, and tissue-enriched proteins. [BP-1] Choose the resource carefully — `"consensus"` is the safest default.
2. **Statistical enrichment ≠ biological activity.** A significant L-R pair means co-expression, not proven signaling.
3. **Method + database choice strongly affects results.** [BP-1] Use LIANA's consensus to reduce method-specific artifacts.
4. **For higher confidence,** consider NicheNet for complementary intracellular activity estimates. [BP-1]

## Guard Rails

1. **Requires ≥2 cell types.** Refuse if only one cell type is annotated.
2. **Warn if cell types have <10 cells.** Limited statistical power.
3. **Always mention database bias** in results interpretation.
4. **Don't over-interpret individual interactions.** Focus on patterns across cell-type pairs.

## Best-Practice References

- [BP-1] §"Communication events across cells" (p. 557)
- [BP-2] Ch. 22 — Cell-cell communication (LIANA tutorial)
- Dimitrov et al. 2022, Nat Commun & 2024, Nat Cell Biol — LIANA/LIANA+
