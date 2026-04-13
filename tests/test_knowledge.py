#!/usr/bin/env python3
"""Unit tests for scagent.knowledge (MarkerDB)."""

import pytest
from scagent.knowledge import (
    MarkerDB, MarkerHit, ValidationResult,
    _resolve_cell_type, _canonical_lookup, _best_match,
    _CANONICAL,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def db():
    """MarkerDB with canonical only (no external DBs)."""
    return MarkerDB(knowledge_dir=None, celltypist_models=[])


@pytest.fixture
def db_celltypist():
    """MarkerDB with CellTypist Immune_All_Low."""
    return MarkerDB(knowledge_dir=None, celltypist_models=["Immune_All_Low.pkl"])


# ---------------------------------------------------------------------------
# Canonical lookup
# ---------------------------------------------------------------------------

class TestCanonicalLookup:
    def test_nk_cells(self, db):
        hits = db.query("NK cells", sources=["canonical"])
        genes = [h.gene for h in hits]
        assert "NKG7" in genes
        assert "GNLY" in genes
        assert all(h.source == "canonical" for h in hits)

    def test_cd8_t_cells(self, db):
        hits = db.query("CD8+ T cells", sources=["canonical"])
        genes = [h.gene for h in hits]
        assert "CD8A" in genes
        assert "CD3D" in genes

    def test_b_cells(self, db):
        hits = db.query("B cells", sources=["canonical"])
        genes = [h.gene for h in hits]
        assert "CD79A" in genes
        assert "MS4A1" in genes

    def test_monocytes(self, db):
        hits = db.query("Monocytes", sources=["canonical"])
        genes = [h.gene for h in hits]
        assert "CD14" in genes
        assert "LYZ" in genes

    def test_classical_monocytes(self, db):
        hits = db.query("Classical monocytes", sources=["canonical"])
        genes = [h.gene for h in hits]
        assert "CD14" in genes
        assert "MARC1" in genes

    def test_tregs(self, db):
        hits = db.query("Regulatory T cells", sources=["canonical"])
        genes = [h.gene for h in hits]
        assert "FOXP3" in genes
        assert "IL2RA" in genes

    def test_plasma_cells(self, db):
        hits = db.query("Plasma cells", sources=["canonical"])
        genes = [h.gene for h in hits]
        assert "MZB1" in genes
        assert "XBP1" in genes

    def test_mouse_species(self, db):
        hits = db.query("NK cells", species="mouse", sources=["canonical"])
        genes = [h.gene for h in hits]
        assert "Nkg7" in genes  # mouse gene symbol
        assert all(h.species == "mouse" for h in hits)

    def test_unknown_cell_type_returns_empty(self, db):
        hits = db.query("Alien cells", sources=["canonical"])
        assert hits == []

    def test_rank_ordering(self, db):
        hits = db.query("NK cells", sources=["canonical"])
        ranks = [h.rank for h in hits]
        assert ranks == list(range(1, len(ranks) + 1))

    def test_top_n_limits(self, db):
        hits = db.query("NK cells", top_n=3, sources=["canonical"])
        assert len(hits) <= 3


# ---------------------------------------------------------------------------
# Cell type alias resolution
# ---------------------------------------------------------------------------

class TestAliases:
    @pytest.mark.parametrize("alias,expected", [
        ("NK", "NK cells"),
        ("nk cells", "NK cells"),
        ("natural killer", "NK cells"),
        ("Tregs", "Regulatory T cells"),
        ("treg", "Regulatory T cells"),
        ("CD14 monocytes", "Classical monocytes"),
        ("cd16 monocytes", "Non-classical monocytes"),
        ("pdc", "pDC"),
        ("helper t cells", "CD4+ T cells"),
        ("cytotoxic t cells", "CD8+ T cells"),
        ("plasma", "Plasma cells"),
        ("rbc", "Erythrocytes"),
        ("hsc", "HSC"),
    ])
    def test_alias_resolves(self, alias, expected):
        assert _resolve_cell_type(alias) == expected

    def test_canonical_name_unchanged(self):
        assert _resolve_cell_type("NK cells") == "NK cells"

    def test_unknown_returns_as_is(self):
        assert _resolve_cell_type("Alien cells") == "Alien cells"

    def test_query_via_alias(self, db):
        hits = db.query("nk", sources=["canonical"])
        assert len(hits) > 0
        assert hits[0].cell_type == "NK cells"


# ---------------------------------------------------------------------------
# CellTypist integration
# ---------------------------------------------------------------------------

class TestCellTypist:
    def test_celltypist_returns_hits(self, db_celltypist):
        hits = db_celltypist.query("NK cells", sources=["celltypist"])
        assert len(hits) > 0
        assert all("celltypist" in h.source for h in hits)

    def test_celltypist_nk_has_nkg7(self, db_celltypist):
        hits = db_celltypist.query("NK cells", sources=["celltypist"])
        genes = [h.gene for h in hits]
        assert "NKG7" in genes

    def test_celltypist_b_cells(self, db_celltypist):
        hits = db_celltypist.query("B cells", sources=["celltypist"])
        genes = [h.gene for h in hits]
        # CellTypist should find CD79A or MS4A1
        assert any(g in genes for g in ["CD79A", "MS4A1", "CD79B"])

    def test_celltypist_unknown_returns_empty(self, db_celltypist):
        hits = db_celltypist.query("Alien cells", sources=["celltypist"])
        assert hits == []

    def test_combined_canonical_and_celltypist(self, db_celltypist):
        hits = db_celltypist.query("NK cells")
        sources = {h.source for h in hits}
        assert any(s == "canonical" for s in sources)
        assert any("celltypist" in s for s in sources)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

class TestValidation:
    def test_high_confidence(self, db):
        result = db.validate_annotation(
            cluster_markers=["NKG7", "GNLY", "KLRD1", "FCGR3A", "GZMB", "PRF1", "KLRF1"],
            label="NK cells",
        )
        assert result.confidence == "high"
        assert result.overlap_ratio >= 0.3
        assert len(result.matched_markers) >= 3

    def test_low_confidence_wrong_label(self, db):
        result = db.validate_annotation(
            cluster_markers=["CD79A", "MS4A1", "CD19", "PAX5", "BANK1"],
            label="NK cells",
        )
        assert result.confidence == "low"
        assert result.overlap_ratio < 0.15

    def test_alternatives_suggested(self, db):
        result = db.validate_annotation(
            cluster_markers=["CD79A", "MS4A1", "CD19", "PAX5", "BANK1"],
            label="NK cells",
        )
        alt_types = [a["cell_type"] for a in result.alternative_labels]
        assert "B cells" in alt_types or "Naive B cells" in alt_types

    def test_correct_label_has_evidence(self, db):
        result = db.validate_annotation(
            cluster_markers=["FOXP3", "IL2RA", "CTLA4", "TIGIT", "CD4"],
            label="Regulatory T cells",
        )
        assert result.confidence in ("high", "medium")
        assert "FOXP3" in result.matched_markers

    def test_validation_via_alias(self, db):
        result = db.validate_annotation(
            cluster_markers=["NKG7", "GNLY", "KLRD1", "GZMB"],
            label="nk",
        )
        assert result.proposed_label == "NK cells"
        assert result.confidence in ("high", "medium")


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

class TestBestMatch:
    def test_exact(self):
        assert _best_match("NK cells", ["NK cells", "T cells"]) == "NK cells"

    def test_case_insensitive(self):
        assert _best_match("nk cells", ["NK cells", "T cells"]) == "NK cells"

    def test_substring(self):
        assert _best_match("NK", ["NK cells", "T cells"]) == "NK cells"

    def test_no_match(self):
        assert _best_match("Alien", ["NK cells", "T cells"]) is None


class TestListCellTypes:
    def test_canonical_list(self, db):
        types = db.list_cell_types("canonical")
        assert "NK cells" in types
        assert "T cells" in types
        assert "B cells" in types
        assert len(types) == len(_CANONICAL)

    def test_celltypist_list(self, db_celltypist):
        types = db_celltypist.list_cell_types("celltypist")
        assert len(types) > 30  # Immune_All_Low has 98 types


class TestStatus:
    def test_status_no_external(self, db):
        s = db.status()
        assert s["canonical_cell_types"] == len(_CANONICAL)
        assert s["total_external_entries"] == 0

    def test_status_with_celltypist(self, db_celltypist):
        # Force load
        db_celltypist.query("NK cells")
        s = db_celltypist.status()
        assert "Immune_All_Low.pkl" in s["celltypist_models_loaded"]


# ---------------------------------------------------------------------------
# Tool wrappers
# ---------------------------------------------------------------------------

class TestToolWrappers:
    def test_query_markers_tool(self):
        from scagent.tools.knowledge_tools import query_markers
        result = query_markers("NK cells", sources=["canonical"])
        assert result["cell_type"] == "NK cells"
        assert len(result["markers"]) > 0
        assert "canonical" in result["sources"]

    def test_validate_annotation_tool(self):
        from scagent.tools.knowledge_tools import validate_annotation
        result = validate_annotation(
            cluster_markers=["NKG7", "GNLY", "GZMB", "PRF1"],
            label="NK cells",
        )
        assert result["label"] == "NK cells"
        assert result["confidence"] in ("high", "medium", "low")
        assert "overlap_ratio" in result

    def test_validate_annotation_returns_alternatives(self):
        from scagent.tools.knowledge_tools import validate_annotation
        result = validate_annotation(
            cluster_markers=["CD79A", "MS4A1", "CD19"],
            label="T cells",
        )
        assert len(result["alternatives"]) > 0


# ---------------------------------------------------------------------------
# External DB loading (with fixtures)
# ---------------------------------------------------------------------------

class TestExternalDB:
    def test_no_dir_no_crash(self):
        db = MarkerDB(knowledge_dir="/nonexistent/path")
        hits = db.query("NK cells", sources=["cellmarker2", "panglaodb"])
        assert hits == []

    def test_empty_dir_no_crash(self, tmp_path):
        db = MarkerDB(knowledge_dir=tmp_path)
        hits = db.query("NK cells", sources=["cellmarker2", "panglaodb"])
        assert hits == []

    def test_panglaodb_tsv_loading(self, tmp_path):
        """Create a minimal PanglaoDB-style TSV and load it."""
        tsv = tmp_path / "PanglaoDB_markers.tsv"
        tsv.write_text(
            "species\tofficial gene symbol\tcell type\tsensitivity_human\n"
            "Hs\tNKG7\tNK cells\t1.0\n"
            "Hs\tGNLY\tNK cells\t0.95\n"
            "Hs\tCD3D\tT cells\t1.0\n"
        )
        db = MarkerDB(knowledge_dir=tmp_path, celltypist_models=[])
        hits = db.query("NK cells", sources=["panglaodb"])
        genes = [h.gene for h in hits]
        assert "NKG7" in genes
        assert "GNLY" in genes
        assert "CD3D" not in genes  # T cell marker, not NK
