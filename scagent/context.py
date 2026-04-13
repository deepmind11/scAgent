"""Experiment context for scAgent projects.

Stores structured metadata about a single-cell experiment following
the minSCe guidelines (Füllgrabe et al. 2020, *Nature Biotechnology*
38:1384–1386).  Gates all downstream analysis decisions — which DAG
to generate, whether batch correction is needed, whether DE should
be pseudobulk, etc.

Usage::

    from pathlib import Path
    from scagent.context import ExperimentContext

    ctx = ExperimentContext(Path(".scagent"))
    ctx.paradigm = "disease_vs_healthy"
    ctx.organism = {"species": "Homo sapiens", "ncbi_taxon": 9606}
    ctx.tissue = {"name": "PBMCs"}
    ctx.save()
"""

from __future__ import annotations

import json
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Valid values
# ---------------------------------------------------------------------------

VALID_PARADIGMS = frozenset({
    "cell_atlas",
    "disease_vs_healthy",
    "developmental_trajectory",
    "perturbation_screen",
    "temporal_longitudinal",
})

VALID_LIBRARY_TYPES = frozenset({"3prime", "5prime", "full_length"})

VALID_BIOSOURCE_TYPES = frozenset({
    "specimen_from_organism",
    "cell_line",
    "organoid",
    "cell_suspension",
})

# Paradigms that require cross-condition design
CONDITION_REQUIRED_PARADIGMS = frozenset({
    "disease_vs_healthy",
    "perturbation_screen",
    "temporal_longitudinal",
})


# ---------------------------------------------------------------------------
# ExperimentContext
# ---------------------------------------------------------------------------

class ExperimentContext:
    """Structured experiment metadata persisted as ``.scagent/project.json``.

    Parameters
    ----------
    project_dir
        Path to the ``.scagent`` directory.
    """

    FILENAME = "project.json"

    def __init__(self, project_dir: Path | str) -> None:
        self._dir = Path(project_dir)
        self._path = self._dir / self.FILENAME
        self._data: dict[str, Any] = self._empty()

        if self._path.exists():
            self._load()

    # ------------------------------------------------------------------
    # Properties (typed access to the context fields)
    # ------------------------------------------------------------------

    @property
    def paradigm(self) -> str | None:
        return self._data.get("paradigm")

    @paradigm.setter
    def paradigm(self, value: str) -> None:
        if value not in VALID_PARADIGMS:
            raise ValueError(
                f"Invalid paradigm '{value}'. Must be one of: {sorted(VALID_PARADIGMS)}"
            )
        self._data["paradigm"] = value

    @property
    def organism(self) -> dict[str, Any]:
        return self._data.get("organism", {})

    @organism.setter
    def organism(self, value: dict[str, Any]) -> None:
        self._data["organism"] = value

    @property
    def tissue(self) -> dict[str, Any]:
        return self._data.get("tissue", {})

    @tissue.setter
    def tissue(self, value: dict[str, Any]) -> None:
        self._data["tissue"] = value

    @property
    def platform(self) -> dict[str, Any]:
        return self._data.get("platform", {})

    @platform.setter
    def platform(self, value: dict[str, Any]) -> None:
        self._data["platform"] = value

    @property
    def library(self) -> dict[str, Any]:
        return self._data.get("library", {})

    @library.setter
    def library(self, value: dict[str, Any]) -> None:
        self._data["library"] = value

    @property
    def samples(self) -> list[dict[str, Any]]:
        return self._data.get("samples", [])

    @samples.setter
    def samples(self, value: list[dict[str, Any]]) -> None:
        self._data["samples"] = value

    @property
    def design(self) -> dict[str, Any]:
        return self._data.get("design", {})

    @design.setter
    def design(self, value: dict[str, Any]) -> None:
        self._data["design"] = value

    @property
    def isolation(self) -> dict[str, Any]:
        return self._data.get("isolation", {})

    @isolation.setter
    def isolation(self, value: dict[str, Any]) -> None:
        self._data["isolation"] = value

    @property
    def hypotheses(self) -> list[str]:
        return self._data.get("hypotheses", [])

    @hypotheses.setter
    def hypotheses(self, value: list[str]) -> None:
        self._data["hypotheses"] = value

    @property
    def project_id(self) -> str:
        return self._data.get("project_id", "")

    @property
    def raw(self) -> dict[str, Any]:
        """Return the full raw context dict."""
        return deepcopy(self._data)

    # ------------------------------------------------------------------
    # Inference from data
    # ------------------------------------------------------------------

    @staticmethod
    def infer_from_data(adata: Any) -> dict[str, Any]:
        """Auto-detect experiment metadata from an AnnData object.

        Returns a dict of inferred fields (each with ``inferred: True``).
        Does NOT modify the AnnData.
        """
        inferred: dict[str, Any] = {}

        # Species — from gene name prefixes
        var_names = list(adata.var_names[:1000])
        n_upper = sum(1 for g in var_names if g == g.upper())
        frac_upper = n_upper / len(var_names) if var_names else 0

        if frac_upper > 0.7:
            inferred["organism"] = {
                "species": "Homo sapiens",
                "ncbi_taxon": 9606,
                "inferred": True,
            }
        elif frac_upper < 0.3:
            inferred["organism"] = {
                "species": "Mus musculus",
                "ncbi_taxon": 10090,
                "inferred": True,
            }

        # Cell / gene counts
        inferred["_data_shape"] = {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        }

        # Platform hints — UMI-based if sparse integer counts
        try:
            import scipy.sparse as sp

            if sp.issparse(adata.X):
                sample = adata.X[:100].toarray()
            else:
                sample = adata.X[:100]
            is_integer = (sample == sample.astype(int)).all()
            if is_integer:
                inferred["library"] = {
                    "umi": True,
                    "type": "3prime",
                    "inferred": True,
                }
                inferred["platform"] = {
                    "vendor": "10x Genomics",
                    "instrument": "Chromium",
                    "inferred": True,
                }
        except Exception:
            pass

        # Batch key hints — look for common column names
        if hasattr(adata, "obs") and hasattr(adata.obs, "columns"):
            obs_cols = set(adata.obs.columns)
            for candidate in ("batch", "sample", "donor", "library_id", "channel"):
                if candidate in obs_cols:
                    n_batches = adata.obs[candidate].nunique()
                    if n_batches > 1:
                        inferred["_batch_hint"] = {
                            "column": candidate,
                            "n_values": int(n_batches),
                        }
                        break

        return inferred

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate(self) -> list[str]:
        """Return a list of validation errors.  Empty list = valid.

        Paradigm is optional — if set, it must be a valid value, but
        analysis can proceed without one (question-driven mode).
        """
        errors: list[str] = []

        if self.paradigm is not None and self.paradigm not in VALID_PARADIGMS:
            errors.append(f"paradigm '{self.paradigm}' is not valid")

        if not self.organism.get("species"):
            errors.append("organism.species is required")

        if not self.tissue.get("name"):
            errors.append("tissue.name is required")

        lib_type = self.library.get("type")
        if lib_type and lib_type not in VALID_LIBRARY_TYPES:
            errors.append(f"library.type '{lib_type}' is not valid")

        # Paradigm-specific validation
        if self.paradigm in CONDITION_REQUIRED_PARADIGMS:
            conditions = self.design.get("conditions", [])
            if len(conditions) < 2:
                errors.append(
                    f"paradigm '{self.paradigm}' requires at least 2 conditions "
                    f"in design.conditions (got {len(conditions)})"
                )

        return errors

    def is_complete(self) -> bool:
        """Return *True* if all analysis-critical fields are set."""
        return len(self.validate()) == 0

    # ------------------------------------------------------------------
    # Analysis decision helpers
    # ------------------------------------------------------------------

    def needs_batch_correction(self) -> bool:
        """True if the experiment has multiple batches that need correction."""
        batch_key = self.design.get("batch_key")
        if not batch_key:
            return len(self.samples) > 1
        return True

    def needs_pseudobulk_de(self) -> bool:
        """True if the paradigm involves cross-condition comparison."""
        return self.paradigm in CONDITION_REQUIRED_PARADIGMS

    def needs_trajectory(self) -> bool:
        """True if the paradigm involves trajectory / pseudotime analysis."""
        return self.paradigm == "developmental_trajectory"

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self) -> Path:
        """Write to ``.scagent/project.json``."""
        self._dir.mkdir(parents=True, exist_ok=True)
        self._path.write_text(
            json.dumps(self._data, indent=2, default=str), encoding="utf-8"
        )
        return self._path

    def _load(self) -> None:
        self._data = json.loads(self._path.read_text(encoding="utf-8"))

    @classmethod
    def load(cls, path: Path | str) -> "ExperimentContext":
        p = Path(path)
        if p.name == cls.FILENAME:
            p = p.parent
        return cls(p)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _empty() -> dict[str, Any]:
        return {
            "@context": "https://scagent.dev/context/v1",
            "@type": "SingleCellExperiment",
            "project_id": f"proj_{datetime.now(timezone.utc).strftime('%Y_%m_%d_%H%M%S')}",
            "created": datetime.now(timezone.utc).isoformat(),
            "paradigm": None,
            "organism": {},
            "tissue": {},
            "biosource_type": "specimen_from_organism",
            "platform": {},
            "library": {"type": "3prime", "modalities": ["gene_expression"], "umi": True},
            "samples": [],
            "design": {},
            "isolation": {"method": "microfluidics_droplet"},
            "hypotheses": [],
            "status": "initialized",
        }
