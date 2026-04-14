"""Tests for immune repertoire tools."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import anndata as ad


def _has_scirpy() -> bool:
    try:
        import scirpy
        return True
    except ImportError:
        return False


class TestLoadVDJValidation:
    def test_import_error_message(self):
        if _has_scirpy():
            pytest.skip("scirpy is installed")

        from scagent.tools.repertoire import load_vdj

        adata = ad.AnnData(np.random.rand(10, 5).astype(np.float32))
        with pytest.raises(ImportError, match="scirpy"):
            load_vdj(adata, vdj_path="/nonexistent/path.csv")

    def test_file_not_found(self):
        if not _has_scirpy():
            pytest.skip("scirpy not installed")

        from scagent.tools.repertoire import load_vdj

        adata = ad.AnnData(np.random.rand(10, 5).astype(np.float32))
        with pytest.raises(FileNotFoundError):
            load_vdj(adata, vdj_path="/nonexistent/path.csv")


class TestRepertoireDAG:
    def test_immune_repertoire_dag_has_vdj_steps(self):
        from unittest.mock import MagicMock
        from scagent.dag import _immune_repertoire_steps

        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        steps = _immune_repertoire_steps(ctx)
        step_ids = [s.id for s in steps]

        assert "load_vdj" in step_ids
        assert "clonotype_analysis" in step_ids
        assert "repertoire_overlap" in step_ids

    def test_clonotype_after_vdj_and_annotation(self):
        from unittest.mock import MagicMock
        from scagent.dag import _immune_repertoire_steps

        ctx = MagicMock()
        ctx.needs_batch_correction.return_value = False

        steps = _immune_repertoire_steps(ctx)
        clono = next(s for s in steps if s.id == "clonotype_analysis")

        assert "load_vdj" in clono.depends_on
        assert "annotation" in clono.depends_on
