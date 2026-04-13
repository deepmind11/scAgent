"""Shared fixtures for integration tests that need real PBMC 10k data.

test_qc_pipeline.py needs raw-loaded adata (with QC metrics computed so
filter thresholds can reference pct_counts_mt).

test_clustering_pipeline.py needs fully preprocessed adata (post-QC,
normalized, ready for HVG → PCA → clustering).
"""

import pytest

DATA_FILE = "data/pbmc10k/filtered_feature_bc_matrix.h5"


@pytest.fixture(scope="module")
def adata(request):
    """Provide the right starting adata depending on which test module is running.

    - test_qc_pipeline       → raw loaded + QC metrics computed
    - test_clustering_pipeline → fully preprocessed (load → QC → normalize)
    """
    module_name = request.module.__name__

    if "clustering" in module_name:
        from scagent.tools.pipeline import run_preprocessing

        adata, _info = run_preprocessing(
            DATA_FILE, max_pct_mito=20.0, random_state=0
        )
        return adata
    else:
        # Default: raw loaded with QC metrics (needed for filter thresholds)
        from scagent.tools.loading import load_10x_h5
        from scagent.tools.qc import calculate_qc_metrics

        adata, _info = load_10x_h5(DATA_FILE)
        calculate_qc_metrics(adata)
        return adata
