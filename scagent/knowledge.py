"""Marker gene knowledge base for scAgent.

Three sources, queried in priority order:

1. **Built-in canonical markers** — curated PBMC/immune markers, always
   available, no download needed.
2. **CellTypist models** — extract top discriminative genes from any of
   61 trained models via ``extract_top_markers()``.
3. **User-supplied databases** — load CellMarker 2.0 xlsx or PanglaoDB TSV
   if the user drops them in ``.scagent/knowledge/``.

Usage::

    from scagent.knowledge import MarkerDB
    db = MarkerDB()
    db.query("NK cells")
    db.query("CD8+ T cells", species="human", tissue="PBMC")
    db.validate_annotation(cluster_markers=["GNLY","NKG7","GZMB"], label="NK cells")
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

@dataclass
class MarkerHit:
    """One marker gene hit from a query."""
    gene: str
    cell_type: str
    source: str            # "canonical", "celltypist:<model>", "cellmarker2", "panglaodb"
    species: str = "human"
    tissue: str = ""
    rank: int = 0          # 1-based rank within source (lower = stronger)
    weight: float = 0.0    # model coefficient or evidence count

    def __repr__(self):
        src = self.source
        return f"MarkerHit({self.gene}, {self.cell_type}, src={src}, rank={self.rank})"


@dataclass
class ValidationResult:
    """Result of comparing cluster markers against known cell types."""
    proposed_label: str
    matched_markers: list[str] = field(default_factory=list)
    total_known: int = 0
    overlap_ratio: float = 0.0
    confidence: str = "low"    # low / medium / high
    evidence: str = ""
    alternative_labels: list[dict] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Built-in canonical markers (curated, always available)
# ---------------------------------------------------------------------------

# Structure: cell_type → {species → {tissue → [genes]}}
# Genes ordered by importance (most discriminative first).
# Sources: CellTypist Immune_All_Low, PanglaoDB, CellMarker 2.0, literature

_CANONICAL: dict[str, dict[str, dict[str, list[str]]]] = {
    # --- T cells ---
    "T cells": {
        "human": {"": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC", "TRBC1", "TRBC2"]},
        "mouse": {"": ["Cd3d", "Cd3e", "Cd3g", "Cd2", "Trac"]},
    },
    "CD4+ T cells": {
        "human": {"": ["CD4", "CD3D", "CD3E", "IL7R", "TRAC", "TCF7", "LEF1"]},
    },
    "CD8+ T cells": {
        "human": {"": ["CD8A", "CD8B", "CD3D", "CD3E", "GZMK", "GZMA"]},
    },
    "Naive CD4+ T cells": {
        "human": {"": ["CCR7", "SELL", "LEF1", "TCF7", "CD4", "CD3D", "IL7R"]},
    },
    "Memory CD4+ T cells": {
        "human": {"": ["IL7R", "S100A4", "CD4", "CD3D", "CD69", "GPR183"]},
    },
    "Regulatory T cells": {
        "human": {"": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "CD4"]},
    },
    "CD8+ effector T cells": {
        "human": {"": ["GZMB", "GNLY", "PRF1", "NKG7", "CD8A", "GZMH", "FGFBP2"]},
    },
    "Gamma-delta T cells": {
        "human": {"": ["TRDC", "TRGC1", "TRGC2", "TRDV2", "CD3D", "NKG7"]},
    },
    # --- NK cells ---
    "NK cells": {
        "human": {"": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "FCGR3A", "GZMB", "PRF1"]},
        "mouse": {"": ["Nkg7", "Gzma", "Klrb1c", "Ncr1", "Prf1"]},
    },
    "CD56bright NK cells": {
        "human": {"": ["NCAM1", "XCL1", "XCL2", "GZMK", "KLRC1", "SELL"]},
    },
    "CD56dim NK cells": {
        "human": {"": ["FCGR3A", "FGFBP2", "GZMB", "PRF1", "SPON2", "NKG7"]},
    },
    # --- B cells ---
    "B cells": {
        "human": {"": ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "PAX5", "CD74"]},
        "mouse": {"": ["Cd79a", "Cd79b", "Ms4a1", "Cd19", "Pax5"]},
    },
    "Naive B cells": {
        "human": {"": ["MS4A1", "CD79A", "TCL1A", "IGHD", "IGHM", "FCER2", "IL4R"]},
    },
    "Memory B cells": {
        "human": {"": ["CD27", "MS4A1", "CD79A", "TNFRSF13B", "AIM2"]},
    },
    "Plasma cells": {
        "human": {"": ["MZB1", "SDC1", "IGHG1", "JCHAIN", "XBP1", "PRDM1"]},
    },
    "Plasmablasts": {
        "human": {"": ["MZB1", "XBP1", "JCHAIN", "IGHG1", "SDC1", "IRF4"]},
    },
    # --- Myeloid ---
    "Monocytes": {
        "human": {"": ["CD14", "LYZ", "S100A9", "S100A8", "VCAN", "FCN1"]},
        "mouse": {"": ["Lyz2", "Ly6c2", "Ccr2", "Csf1r"]},
    },
    "Classical monocytes": {
        "human": {"": ["CD14", "LYZ", "S100A9", "S100A8", "VCAN", "FCN1", "MARC1"]},
    },
    "Non-classical monocytes": {
        "human": {"": ["FCGR3A", "MS4A7", "CX3CR1", "CDKN1C", "HES4", "LILRB2"]},
    },
    "Macrophages": {
        "human": {"": ["CD68", "CD163", "MRC1", "MSR1", "MARCO", "C1QA", "C1QB"]},
    },
    "Dendritic cells": {
        "human": {"": ["FCER1A", "CST3", "HLA-DQA1", "HLA-DPB1", "CLEC10A"]},
    },
    "cDC1": {
        "human": {"": ["CLEC9A", "XCR1", "BATF3", "IRF8", "CADM1", "WDFY4"]},
    },
    "cDC2": {
        "human": {"": ["CD1C", "FCER1A", "CLEC10A", "HLA-DQA1", "JAML"]},
    },
    "pDC": {
        "human": {"": ["LILRA4", "IL3RA", "JCHAIN", "GZMB", "IRF7", "ITM2C"]},
    },
    # --- Granulocytes ---
    "Neutrophils": {
        "human": {"": ["CSF3R", "FCGR3B", "CXCR2", "S100A8", "S100A9", "MMP9"]},
    },
    "Mast cells": {
        "human": {"": ["TPSAB1", "TPSB2", "CPA3", "HPGDS", "KIT", "GATA2"]},
    },
    # --- Other ---
    "Platelets": {
        "human": {"": ["PPBP", "PF4", "GP9", "ITGA2B", "TUBB1", "TREML1"]},
    },
    "Erythrocytes": {
        "human": {"": ["HBB", "HBA1", "HBA2", "HBD", "ALAS2", "SLC4A1"]},
    },
    "HSC": {
        "human": {"": ["CD34", "CRHBP", "AVP", "MLLT3", "SPINK2"]},
    },
    "Endothelial cells": {
        "human": {"": ["PECAM1", "CDH5", "VWF", "FLT1", "KDR", "ENG"]},
    },
    "Fibroblasts": {
        "human": {"": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "PDGFRA"]},
    },
    "Epithelial cells": {
        "human": {"": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "CLDN4"]},
    },
}

# Aliases → canonical name
_CELL_TYPE_ALIASES: dict[str, str] = {
    "t cells": "T cells", "t cell": "T cells",
    "cd4 t cells": "CD4+ T cells", "cd4+ t": "CD4+ T cells",
    "cd4 t": "CD4+ T cells", "helper t cells": "CD4+ T cells",
    "cd8 t cells": "CD8+ T cells", "cd8+ t": "CD8+ T cells",
    "cd8 t": "CD8+ T cells", "cytotoxic t cells": "CD8+ T cells",
    "tregs": "Regulatory T cells", "treg": "Regulatory T cells",
    "cd8 effector": "CD8+ effector T cells",
    "gd t cells": "Gamma-delta T cells", "gamma delta t": "Gamma-delta T cells",
    "nk cells": "NK cells", "nk": "NK cells", "natural killer": "NK cells",
    "b cells": "B cells", "b cell": "B cells",
    "plasma cells": "Plasma cells", "plasma": "Plasma cells",
    "monocytes": "Monocytes", "mono": "Monocytes",
    "cd14 monocytes": "Classical monocytes", "cd14+ monocytes": "Classical monocytes",
    "cd16 monocytes": "Non-classical monocytes", "cd16+ monocytes": "Non-classical monocytes",
    "classical mono": "Classical monocytes",
    "non-classical mono": "Non-classical monocytes", "ncm": "Non-classical monocytes",
    "macrophages": "Macrophages", "mac": "Macrophages",
    "dc": "Dendritic cells", "dendritic": "Dendritic cells",
    "cdc1": "cDC1", "cdc2": "cDC2", "pdc": "pDC",
    "conventional dc1": "cDC1", "conventional dc2": "cDC2",
    "plasmacytoid dc": "pDC", "plasmacytoid dendritic cells": "pDC",
    "neutrophils": "Neutrophils", "neut": "Neutrophils",
    "mast cells": "Mast cells", "mast": "Mast cells",
    "platelets": "Platelets", "plt": "Platelets", "megakaryocytes": "Platelets",
    "erythrocytes": "Erythrocytes", "rbc": "Erythrocytes", "red blood cells": "Erythrocytes",
    "hsc": "HSC", "stem cells": "HSC", "hematopoietic stem cells": "HSC",
    "endothelial": "Endothelial cells", "fibroblasts": "Fibroblasts",
    "epithelial": "Epithelial cells",
}


def _resolve_cell_type(name: str) -> str:
    """Normalize a cell type name to canonical form."""
    if name in _CANONICAL:
        return name
    return _CELL_TYPE_ALIASES.get(name.lower().strip(), name)


# ---------------------------------------------------------------------------
# MarkerDB
# ---------------------------------------------------------------------------

class MarkerDB:
    """Unified marker gene knowledge base.

    Parameters
    ----------
    knowledge_dir
        Directory for user-supplied marker databases. If CellMarker 2.0
        xlsx or PanglaoDB TSV files are placed here, they are loaded
        automatically.
    celltypist_models
        List of CellTypist model names to load for marker extraction.
        Default: ``["Immune_All_Low.pkl"]``.
    """

    def __init__(
        self,
        knowledge_dir: str | Path | None = None,
        celltypist_models: list[str] | None = None,
    ):
        self.knowledge_dir = Path(knowledge_dir) if knowledge_dir else None
        self._ct_models: dict[str, object] = {}  # model_name → Model
        self._ct_model_names = celltypist_models if celltypist_models is not None else ["Immune_All_Low.pkl"]
        self._external_markers: list[dict] = []  # from CellMarker2/PanglaoDB
        self._loaded_external = False

    # ------------------------------------------------------------------
    # CellTypist (lazy)
    # ------------------------------------------------------------------

    def _ensure_celltypist(self, model_name: str):
        if model_name not in self._ct_models:
            try:
                import celltypist
                self._ct_models[model_name] = celltypist.models.Model.load(
                    model=model_name,
                )
                logger.info("Loaded CellTypist model: %s", model_name)
            except Exception as e:
                logger.warning("Failed to load CellTypist model %s: %s", model_name, e)

    def _celltypist_markers(
        self, cell_type: str, model_name: str, top_n: int = 15,
    ) -> list[MarkerHit]:
        """Extract top markers from a CellTypist model."""
        self._ensure_celltypist(model_name)
        model = self._ct_models.get(model_name)
        if model is None:
            return []
        # Find best matching cell type in model
        ct_match = _best_match(cell_type, list(model.cell_types))
        if ct_match is None:
            return []
        try:
            genes = model.extract_top_markers(cell_type=ct_match, top_n=top_n)
            return [
                MarkerHit(
                    gene=g, cell_type=ct_match,
                    source=f"celltypist:{model_name}",
                    species="human", rank=i + 1,
                )
                for i, g in enumerate(genes)
            ]
        except Exception as e:
            logger.warning("extract_top_markers failed for %s: %s", ct_match, e)
            return []

    # ------------------------------------------------------------------
    # External databases (lazy)
    # ------------------------------------------------------------------

    def _load_external(self):
        """Load CellMarker 2.0 xlsx or PanglaoDB TSV if present."""
        if self._loaded_external:
            return
        self._loaded_external = True
        if self.knowledge_dir is None or not self.knowledge_dir.exists():
            return

        # CellMarker 2.0 — expects xlsx with columns:
        #   species, tissue_class, cell_name, cell_marker, ...
        for xlsx in self.knowledge_dir.glob("*[Cc]ell[Mm]arker*.xlsx"):
            self._load_cellmarker2(xlsx)

        # PanglaoDB — expects TSV with columns:
        #   species, official gene symbol, cell type, ...
        for tsv in self.knowledge_dir.glob("*[Pp]anglao*.tsv"):
            self._load_panglaodb(tsv)

    def _load_cellmarker2(self, path: Path):
        try:
            import pandas as pd
            df = pd.read_excel(path)
            # Normalize column names
            cols = {c.lower().replace(" ", "_"): c for c in df.columns}
            species_col = cols.get("species", cols.get("speceis", ""))
            cell_col = cols.get("cell_name", cols.get("cell_type", ""))
            marker_col = cols.get("cell_marker", cols.get("marker", ""))
            tissue_col = cols.get("tissue_class", cols.get("tissue_type", ""))
            if not (cell_col and marker_col):
                logger.warning("CellMarker2 xlsx missing expected columns: %s", list(df.columns))
                return
            for _, row in df.iterrows():
                markers_raw = str(row.get(marker_col, ""))
                for gene in markers_raw.replace(",", " ").replace(";", " ").split():
                    gene = gene.strip().strip("[]")
                    if gene and gene.upper() != "NA":
                        self._external_markers.append({
                            "gene": gene,
                            "cell_type": str(row.get(cell_col, "")),
                            "species": str(row.get(species_col, "Human")).lower(),
                            "tissue": str(row.get(tissue_col, "")),
                            "source": "cellmarker2",
                        })
            logger.info("Loaded CellMarker 2.0: %d entries from %s",
                        len(self._external_markers), path.name)
        except Exception as e:
            logger.warning("Failed to load CellMarker2 %s: %s", path, e)

    def _load_panglaodb(self, path: Path):
        try:
            import pandas as pd
            df = pd.read_csv(path, sep="\t")
            cols = {c.lower().replace(" ", "_"): c for c in df.columns}
            gene_col = cols.get("official_gene_symbol", cols.get("gene", ""))
            cell_col = cols.get("cell_type", "")
            species_col = cols.get("species", "")
            if not (gene_col and cell_col):
                logger.warning("PanglaoDB TSV missing expected columns: %s", list(df.columns))
                return
            n_before = len(self._external_markers)
            for _, row in df.iterrows():
                gene = str(row.get(gene_col, "")).strip()
                if gene and gene.upper() != "NA":
                    species_val = str(row.get(species_col, "")).lower()
                    self._external_markers.append({
                        "gene": gene,
                        "cell_type": str(row.get(cell_col, "")),
                        "species": species_val if species_val else "human",
                        "tissue": "",
                        "source": "panglaodb",
                    })
            logger.info("Loaded PanglaoDB: %d entries from %s",
                        len(self._external_markers) - n_before, path.name)
        except Exception as e:
            logger.warning("Failed to load PanglaoDB %s: %s", path, e)

    def _external_hits(
        self, cell_type: str, species: str, tissue: str,
    ) -> list[MarkerHit]:
        self._load_external()
        ct_lower = cell_type.lower()
        hits = []
        for entry in self._external_markers:
            if ct_lower not in entry["cell_type"].lower():
                continue
            if species and not _species_match(species, entry.get("species", "")):
                continue
            if tissue and tissue.lower() not in entry.get("tissue", "").lower():
                continue
            hits.append(MarkerHit(
                gene=entry["gene"], cell_type=entry["cell_type"],
                source=entry["source"], species=entry.get("species", "human"),
                tissue=entry.get("tissue", ""),
            ))
        # Deduplicate by gene and assign ranks
        seen = set()
        ranked = []
        for h in hits:
            if h.gene not in seen:
                seen.add(h.gene)
                h.rank = len(ranked) + 1
                ranked.append(h)
        return ranked

    # ------------------------------------------------------------------
    # Query
    # ------------------------------------------------------------------

    def query(
        self,
        cell_type: str,
        *,
        species: str = "human",
        tissue: str = "",
        top_n: int = 15,
        sources: list[str] | None = None,
    ) -> list[MarkerHit]:
        """Query markers for a cell type.

        Parameters
        ----------
        cell_type
            Cell type name (aliases resolved automatically).
        species
            ``"human"`` or ``"mouse"``.
        tissue
            Optional tissue filter (e.g. ``"PBMC"``, ``"lung"``).
        top_n
            Max markers per source.
        sources
            Restrict to specific sources. Default: all available.
            Options: ``"canonical"``, ``"celltypist"``, ``"cellmarker2"``,
            ``"panglaodb"``.

        Returns
        -------
        List of :class:`MarkerHit`, ordered by source priority then rank.
        """
        resolved = _resolve_cell_type(cell_type)
        all_hits: list[MarkerHit] = []
        use = set(sources) if sources else {"canonical", "celltypist", "cellmarker2", "panglaodb"}

        # 1. Built-in canonical
        if "canonical" in use:
            canon = _canonical_lookup(resolved, species, tissue)
            all_hits.extend(canon[:top_n])

        # 2. CellTypist
        if "celltypist" in use:
            for model_name in self._ct_model_names:
                ct_hits = self._celltypist_markers(resolved, model_name, top_n)
                all_hits.extend(ct_hits)

        # 3. External databases
        if use & {"cellmarker2", "panglaodb"}:
            ext = self._external_hits(resolved, species, tissue)
            if "cellmarker2" not in use:
                ext = [h for h in ext if h.source != "cellmarker2"]
            if "panglaodb" not in use:
                ext = [h for h in ext if h.source != "panglaodb"]
            all_hits.extend(ext[:top_n])

        return all_hits

    # ------------------------------------------------------------------
    # Validate annotation
    # ------------------------------------------------------------------

    def validate_annotation(
        self,
        cluster_markers: list[str],
        label: str,
        *,
        species: str = "human",
        top_n: int = 20,
    ) -> ValidationResult:
        """Check if DE markers support a proposed cell type label.

        Parameters
        ----------
        cluster_markers
            Top DE genes for the cluster (from ``rank_genes_groups``).
        label
            Proposed cell type label.
        species
            ``"human"`` or ``"mouse"``.
        top_n
            How many known markers to compare against.

        Returns
        -------
        :class:`ValidationResult` with overlap ratio, confidence, and
        alternative cell type suggestions.
        """
        resolved = _resolve_cell_type(label)
        known = self.query(resolved, species=species, top_n=top_n)
        known_genes = {h.gene.upper() for h in known}
        cluster_upper = {g.upper() for g in cluster_markers}

        matched = sorted(known_genes & cluster_upper)
        overlap = len(matched) / len(known_genes) if known_genes else 0.0

        if overlap >= 0.3:
            conf = "high"
        elif overlap >= 0.15:
            conf = "medium"
        else:
            conf = "low"

        evidence_parts = []
        if matched:
            evidence_parts.append(f"Matched markers: {', '.join(sorted(matched))}")
        if known_genes - cluster_upper:
            missing = sorted(known_genes - cluster_upper)[:5]
            evidence_parts.append(f"Missing expected: {', '.join(missing)}")

        # Check alternative labels
        alternatives = self._check_alternatives(cluster_upper, species, exclude=resolved)

        return ValidationResult(
            proposed_label=resolved,
            matched_markers=matched,
            total_known=len(known_genes),
            overlap_ratio=round(overlap, 3),
            confidence=conf,
            evidence="; ".join(evidence_parts),
            alternative_labels=alternatives,
        )

    def _check_alternatives(
        self, cluster_genes: set[str], species: str,
        exclude: str = "", top_k: int = 3,
    ) -> list[dict]:
        """Score all canonical cell types against cluster markers."""
        scores = []
        for ct_name in _CANONICAL:
            if ct_name == exclude:
                continue
            canon = _canonical_lookup(ct_name, species, "")
            known = {h.gene.upper() for h in canon}
            if not known:
                continue
            overlap = len(known & cluster_genes)
            if overlap > 0:
                scores.append({
                    "cell_type": ct_name,
                    "matched": overlap,
                    "total": len(known),
                    "ratio": round(overlap / len(known), 3),
                })
        scores.sort(key=lambda x: x["ratio"], reverse=True)
        return scores[:top_k]

    # ------------------------------------------------------------------
    # List available cell types
    # ------------------------------------------------------------------

    def list_cell_types(self, source: str = "canonical") -> list[str]:
        """List all queryable cell types from a source."""
        if source == "canonical":
            return sorted(_CANONICAL.keys())
        if source == "celltypist":
            types = set()
            for model_name in self._ct_model_names:
                self._ensure_celltypist(model_name)
                model = self._ct_models.get(model_name)
                if model:
                    types.update(model.cell_types)
            return sorted(types)
        return []

    def status(self) -> dict:
        """Summary of loaded sources."""
        self._load_external()
        ct_loaded = list(self._ct_models.keys())
        ext_sources = {}
        for e in self._external_markers:
            s = e["source"]
            ext_sources[s] = ext_sources.get(s, 0) + 1
        return {
            "canonical_cell_types": len(_CANONICAL),
            "celltypist_models_configured": self._ct_model_names,
            "celltypist_models_loaded": ct_loaded,
            "external_databases": ext_sources,
            "total_external_entries": len(self._external_markers),
        }


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _canonical_lookup(
    cell_type: str, species: str, tissue: str,
) -> list[MarkerHit]:
    """Look up canonical markers."""
    entry = _CANONICAL.get(cell_type)
    if not entry:
        return []
    sp = entry.get(species.lower(), entry.get("human", {}))
    # Try tissue-specific first, then generic
    genes = sp.get(tissue.lower(), sp.get("", []))
    return [
        MarkerHit(
            gene=g, cell_type=cell_type, source="canonical",
            species=species, tissue=tissue, rank=i + 1,
        )
        for i, g in enumerate(genes)
    ]


_SPECIES_ALIASES: dict[str, str] = {
    "human": "human", "homo sapiens": "human", "hs": "human", "h. sapiens": "human",
    "mouse": "mouse", "mus musculus": "mouse", "mm": "mouse", "m. musculus": "mouse",
}


def _species_match(query: str, entry_species: str) -> bool:
    """Check if a query species matches an entry species, handling aliases."""
    q = _SPECIES_ALIASES.get(query.lower().strip(), query.lower().strip())
    e = _SPECIES_ALIASES.get(entry_species.lower().strip(), entry_species.lower().strip())
    return q == e


def _best_match(query: str, candidates: list[str]) -> str | None:
    """Find best matching cell type name (case-insensitive substring)."""
    q = query.lower().strip()
    # Exact match
    for c in candidates:
        if c.lower() == q:
            return c
    # Substring match
    for c in candidates:
        if q in c.lower():
            return c
    # Reverse substring
    for c in candidates:
        if c.lower() in q:
            return c
    return None
