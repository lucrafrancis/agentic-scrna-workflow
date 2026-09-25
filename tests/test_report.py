"""The report's placeholder rules: numbers come from code, never typed by the model."""

from __future__ import annotations

import pandas as pd

from agent import report
from agent.session import SESSION


def _render(markdown):
    return report.render(markdown, report.read_log(SESSION.paths.tool_log), SESSION.adata)


def _label_cells():
    """Stand-in for CellTypist (which needs a network download): two labelled cell types."""
    adata = SESSION.adata
    adata.obs["cell_type"] = pd.Categorical(
        ["Type A" if c in ("0", "1") else "Type B" for c in adata.obs["leiden"].astype(str)]
    )


def test_facts_fill_placeholders(logged_run):
    text, problems = _render(
        "# R\n\n## QC\nCutoff {{filter.max_pct_mt}} kept {{filter.cells_after}} of {{input.n_cells}} cells; "
        "doublets used {{doublet.rule}}."
    )
    assert problems == []
    assert "Cutoff 90%" in text
    assert "of 400 cells" in text
    assert "{{" not in text


def test_typed_numbers_are_rejected(logged_run):
    _, problems = _render("# R\n\n## QC\nWe removed 57 cells (2.1%) at a 5% cutoff, 3,291/4,000 at default.")
    flagged = " ".join(problems)
    for number in ("'57'", "'2.1%'", "'5%'", "'3,291'", "'4,000'"):
        assert number in flagged
    assert "'000'" not in flagged


def test_names_headings_and_clusters_are_allowed(logged_run):
    cluster_id = str(SESSION.adata.obs["leiden"].astype(str).iloc[0])
    text, problems = _render(
        f"# PBMC3k report\n\n## 2. Markers\nCD14, S100A8/9, MT-CO1, HLA-DRB1 and PC1 mark cluster {cluster_id} "
        "in the 10x data; doublets used median + 3×MAD.\n\n1. A numbered list item.\n"
    )
    assert problems == []


def test_unknown_placeholders_and_missing_clusters_are_rejected(logged_run):
    _, problems = _render("# R\n\n## QC\n{{no.such.fact}} and cluster 99.")
    assert any("unknown placeholder {{no.such.fact}}" in p for p in problems)
    assert any("cluster 99 does not exist" in p for p in problems)


def test_methods_section_is_rejected(logged_run):
    _, problems = _render("# R\n\n## Methods Summary\nScanpy.")
    assert any("Methods Summary" in p for p in problems)


def test_cell_type_placeholders_and_required_tables(logged_run):
    _label_cells()
    n_a = int((SESSION.adata.obs["cell_type"] == "Type A").sum())
    text, problems = _render("# R\n\n## Annotation\nType A: {{celltype:Type A}}. Both: {{celltypes:Type A|Type B}}.")
    assert f"Type A: {n_a:,} cells" in text
    assert f"Both: {SESSION.adata.n_obs:,} cells (100.0%)" in text
    text, _ = _render("# R\n\n## Annotation\nThere were {{celltype:Type A}} cells.")
    assert f"There were {n_a:,} cells (" in text and "%) cells" not in text
    assert "required table not placed: {{table:clusters}}" in problems
    assert "required table not placed: {{table:composition}}" in problems

    text, problems = _render("# R\n\n## Annotation\n{{table:clusters}}\n\n{{table:composition}}")
    assert problems == []
    assert "| Cluster | Cells | Top marker genes |" in text


def test_generate_report_rejects_then_writes_with_warnings(logged_run):
    bad = "# R\n\n## QC\nWe removed 57 cells."
    for _ in range(2):
        result = logged_run("generate_report", report_markdown=bad)
        assert result["error"] == "report_rejected"
        assert not SESSION.paths.report.exists()
        assert "filter.cells_removed" in result["available_facts"]

    result = logged_run("generate_report", report_markdown=bad)
    assert "warning" in result
    text = SESSION.paths.report.read_text()
    assert "failed automatic checks" in text


def test_report_gets_a_title_when_the_narrative_has_none(logged_run):
    logged_run("summarize_findings")
    logged_run("generate_report", report_markdown="## Overview\nNo title here.")
    text = SESSION.paths.report.read_text()
    assert text.startswith("# Single-cell RNA-seq analysis: synthetic\n")
    assert text.index("## Overview") < text.index("## Key analysis decisions")
