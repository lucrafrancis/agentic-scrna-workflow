"""A fast end-to-end pass on synthetic data, exercising the tool contract.

Runs the pipeline through report generation (skipping scVI training and CellTypist, which
need heavy compute / network). Asserts every tool returns a JSON-serializable dict with no
error, and that the expected artifacts and session state result.
"""

from __future__ import annotations

import json

from agent import tools
from agent.session import SESSION


def test_pipeline_end_to_end(synthetic_h5ad):
    SESSION.load(synthetic_h5ad)

    summaries = [
        tools.inspect_dataset(),
        tools.check_gene_identifiers(),
        tools.compute_qc(),
        tools.recommend_qc_thresholds(),
        tools.filter_cells_and_genes(min_genes=5, max_pct_mt=90, min_cells=1),
        tools.detect_doublets(),
    ]
    summaries.append(tools.filter_doublets(threshold=summaries[-1]["recommended_threshold"]))
    summaries.append(tools.normalize(n_top_genes=50))
    summaries.append(tools.run_pca(n_comps=10))
    summaries.append(tools.cluster(resolution=1.0))
    summaries.append(tools.identify_markers(n_genes=5))
    summaries.append(tools.summarize_findings())
    summaries.append(tools.generate_report(report_markdown="# test report"))

    for summary in summaries:
        assert isinstance(summary, dict)
        assert "error" not in summary, summary
        json.dumps(summary)  # must be JSON-serializable to reach the LLM

    assert SESSION.representation == "X_pca"
    assert SESSION.paths.report.exists()
    assert SESSION.paths.annotated.exists()
    assert (SESSION.paths.figures / "umap.png").exists()


def test_mito_genes_found_in_qc(synthetic_h5ad):
    SESSION.load(synthetic_h5ad)
    tools.check_gene_identifiers()
    assert tools.compute_qc()["n_mito_genes_found"] == 5


def test_doublets_run_per_batch(batched_h5ad):
    """With a batch key, Scrublet's simulated-doublet model is built within each batch."""
    SESSION.load(batched_h5ad)
    summary = tools.detect_doublets()

    assert summary["batch_key"] == "batch"
    per_batch = summary["candidate_thresholds"]["scrublet_auto_per_batch"]
    assert set(per_batch) == {"A", "B"}
    json.dumps(summary)


def test_hvgs_selected_within_batch(batched_h5ad):
    """With a batch key, HVGs are ranked per batch so the batch effect doesn't drive them."""
    SESSION.load(batched_h5ad)
    tools.check_gene_identifiers()
    tools.compute_qc()
    tools.filter_cells_and_genes(min_genes=5, max_pct_mt=90, min_cells=1)
    summary = tools.normalize(n_top_genes=50)

    assert summary["hvg_batch_key"] == "batch"
    assert summary["n_hvgs_flagged"] > 0
    json.dumps(summary)


def test_normalize_records_target_sum(synthetic_h5ad):
    SESSION.load(synthetic_h5ad)
    tools.check_gene_identifiers()
    tools.compute_qc()
    tools.filter_cells_and_genes(min_genes=5, max_pct_mt=90, min_cells=1)
    tools.normalize(target_sum=1e6, n_top_genes=50)
    assert SESSION.normalize_target_sum == 1e6


def test_normalize_is_idempotent_guarded(synthetic_h5ad):
    SESSION.load(synthetic_h5ad)
    tools.check_gene_identifiers()
    tools.compute_qc()
    tools.filter_cells_and_genes(min_genes=5, max_pct_mt=90, min_cells=1)
    assert "checkpoint" in tools.normalize(n_top_genes=50)
    assert tools.normalize(n_top_genes=50)["error"] == "already_normalized"


NARRATIVE = """# Test report

## Overview
Short intro.

## Quality Control
QC text.

## Doublet Detection
Doublet text.

## Clustering
Cluster text.

## Conclusions
Done.
"""


def test_report_built_from_tool_log(logged_run):
    """The report's code-built parts (decisions table, figures, Methods) come from the log."""
    logged_run("summarize_findings")
    result = logged_run("generate_report", report_markdown=NARRATIVE)
    assert "error" not in result, result

    text = SESSION.paths.report.read_text()
    assert "## Key analysis decisions" in text
    assert "| Mitochondrial cutoff | 5% |" in text
    assert "**90%** (changed)" in text
    assert "## Methods" in text
    # The decisions table follows the overview; each figure sits in its own section.
    assert text.index("## Overview") < text.index("## Key analysis decisions") < text.index("## Quality Control")
    qc_section = text[text.index("## Quality Control") : text.index("## Doublet Detection")]
    assert "figures/qc_thresholds.png" in qc_section
    doublet_section = text[text.index("## Doublet Detection") : text.index("## Clustering")]
    assert "figures/doublet_scores.png" in doublet_section
    assert "figures/umap.png" in text[text.index("## Clustering") : text.index("## Conclusions")]
    for name in ("qc_thresholds.png", "doublet_scores.png", "umap.png"):
        assert (SESSION.paths.figures / name).exists()


def test_doublets_batch_key_can_be_chosen(batched_h5ad):
    """The agent can override the detected batch column, e.g. when donors shared one run."""
    SESSION.load(batched_h5ad)
    assert tools.detect_doublets(batch_key=None)["batch_key"] is None
    assert tools.detect_doublets(batch_key="nope")["error"] == "invalid_batch_key"


def test_check_markers_reports_rank_stats_and_absence(logged_run):
    """GENE0-29 are raised in population A only (see make_adata), so in A's cluster GENE0 ranks
    near the top and is widely expressed, and in B's cluster it is not enriched."""
    adata = SESSION.adata
    a_cells = adata.obs_names.str.replace("cell", "").astype(int) < 200
    cluster_a = adata.obs.loc[a_cells, "leiden"].astype(str).value_counts().index[0]
    cluster_b = adata.obs.loc[~a_cells, "leiden"].astype(str).value_counts().index[0]

    result = tools.check_markers(["GENE0", "NOT_A_GENE"], clusters=[cluster_a, cluster_b])

    a, b = result["genes"]["GENE0"][cluster_a], result["genes"]["GENE0"][cluster_b]
    assert int(a["rank"].split(" of ")[0]) <= 30 and a["log2FC"] > 1 and a["pct_in_cluster"] > 90
    assert b["log2FC"] < 0
    assert result["not_found"] == ["NOT_A_GENE"]
    assert tools.check_markers(["GENE0"], clusters=["99"])["error"] == "bad_clusters"
    json.dumps(result)
