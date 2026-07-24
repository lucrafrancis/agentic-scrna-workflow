"""Out-of-order tool calls must return an error summary (not raise), so the agent can react.

None of these import torch or celltypist: each tool's guardrail check runs before any heavy
import, so the suite stays fast and offline.
"""

from __future__ import annotations

import pytest

from agent import tools
from agent.session import SESSION


@pytest.fixture
def loaded(synthetic_h5ad):
    SESSION.load(synthetic_h5ad)


def test_compute_qc_needs_gene_check(loaded):
    assert tools.compute_qc()["error"] == "gene_identifiers_not_checked"


def test_recommend_needs_qc(loaded):
    tools.check_gene_identifiers()
    assert "error" in tools.recommend_qc_thresholds()


def test_filter_needs_qc(loaded):
    tools.check_gene_identifiers()
    assert "error" in tools.filter_cells_and_genes(min_genes=5, max_pct_mt=90, min_cells=1)


def test_detect_doublets_refuses_normalized_data(loaded):
    """Scrublet on log-normalized X only warns on stderr; the agent must see an error."""
    tools.check_gene_identifiers()
    tools.compute_qc()
    tools.filter_cells_and_genes(min_genes=5, max_pct_mt=90, min_cells=1)
    tools.normalize(n_top_genes=50)
    assert tools.detect_doublets()["error"] == "x_not_counts"


def test_filter_doublets_needs_detection(loaded):
    assert tools.filter_doublets(threshold=0.1)["error"] == "doublets_not_detected"


def test_pca_needs_normalize(loaded):
    assert tools.run_pca()["error"] == "not_normalized"


def test_scvi_needs_counts(loaded):
    assert tools.run_scvi()["error"] == "no_raw_counts"


def test_cluster_needs_representation(loaded):
    assert tools.cluster()["error"] == "no_representation"


def test_markers_need_clusters(loaded):
    assert tools.identify_markers()["error"] == "no_clusters"


def test_annotate_needs_clusters(loaded):
    assert tools.annotate_celltypes()["error"] == "no_clusters"
