"""Condition comparisons: composition and pseudobulk DE on synthetic donors x conditions.

4 donors, each with a ctrl and a stim sample, two cell types. One gene is raised in stim
cells, so paired pseudobulk DE must find it; nothing else differs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from agent import report, tools
from agent.session import SESSION
from tests.conftest import make_adata

STIM_GENE = "GENE100"


@pytest.fixture
def two_condition_h5ad(tmp_path):
    adata = make_adata(n_cells=800)
    n = adata.n_obs
    adata.obs["donor"] = pd.Categorical([f"d{i % 4}" for i in range(n)])
    adata.obs["condition"] = pd.Categorical(["ctrl" if (i // 4) % 2 == 0 else "stim" for i in range(n)])
    adata.obs["cell_type"] = pd.Categorical(["Type A" if i < n // 2 else "Type B" for i in range(n)])
    # Raise one gene in stim cells only.
    X = adata.X.toarray()
    stim = (adata.obs["condition"] == "stim").to_numpy()
    X[stim, list(adata.var_names).index(STIM_GENE)] += np.random.default_rng(1).poisson(6.0, stim.sum())
    adata.X = sparse.csr_matrix(X)
    path = tmp_path / "two_condition.h5ad"
    adata.write_h5ad(path)
    return path


def test_compare_composition(two_condition_h5ad):
    SESSION.load(two_condition_h5ad)
    result = tools.compare_composition(["condition", "stim", "ctrl"], sample_key="donor")

    assert "error" not in result, result
    assert result["paired"] is True
    assert result["test"] == "Wilcoxon signed-rank"
    assert result["n_samples"] == {"stim": 4, "ctrl": 4}
    assert set(result["results"]) == {"Type A", "Type B"}
    assert result["n_significant"] == 0  # composition is identical by construction

    # Direction facts carry their own meaning (the model once read "0 of 8" the wrong way).
    text, _ = report.render("# R\n\n## Composition\n{{comp.type_a.direction}}.", [], SESSION.adata)
    assert "higher in stim in" in text and "donors, lower or equal in" in text

    # Guardrails: the condition can't be the replicate, and 2 donors is too few to test.
    assert tools.compare_composition(["condition", "stim", "ctrl"], sample_key="condition")["error"] == "bad_sample_key"
    SESSION.adata = SESSION.adata[SESSION.adata.obs["donor"].isin(["d0", "d1"])].copy()
    assert tools.compare_composition(["condition", "stim", "ctrl"], sample_key="donor")["error"] == "too_few_samples"


def test_run_pseudobulk_de(two_condition_h5ad):
    SESSION.load(two_condition_h5ad)
    SESSION.begin_run()
    result = tools.run_pseudobulk_de(["condition", "stim", "ctrl"], sample_key="donor", covariates=["donor"])

    assert "error" not in result, result
    assert result["design"] == "~donor + condition"
    for ct in ("Type A", "Type B"):
        r = result["results"][ct]
        assert r["n_samples"] == {"ctrl": 4, "stim": 4}
        assert r["top_up"][0]["gene"] == STIM_GENE
        assert r["n_significant"] <= 3  # essentially only the planted gene

    hit = tools.query_genes([STIM_GENE], ["Type A"])["genes"][STIM_GENE]["Type A"]
    assert hit["significant"] and hit["log2FC"] > 1

    # The report cites it through a placeholder, and draws the volcano plot.
    text, problems = report.render("# R\n\n## Differential expression\n{{gene:Type A:GENE100}}. "
                                   "{{table:de_summary}}", [], SESSION.adata)
    assert f"{STIM_GENE} (log2FC" in text
    assert "| Type A | 4 / 4 |" in text
    SESSION.paths.figures.mkdir(parents=True)
    assert report._volcano_figure(SESSION.paths.figures) is not None
    for name in ("volcano_type_a.png", "volcano_type_b.png"):  # one per cell type
        assert (SESSION.paths.figures / name).exists()

    # Guardrails: a covariate identical to the condition is confounded; the condition can't
    # be a covariate; unknown levels are refused.
    SESSION.adata.obs["arm"] = SESSION.adata.obs["condition"].map({"ctrl": "A", "stim": "B"}).astype("category")
    assert tools.run_pseudobulk_de(["condition", "stim", "ctrl"], "donor", covariates=["arm"])["error"] == "confounded_design"
    assert tools.run_pseudobulk_de(["condition", "stim", "ctrl"], "donor", covariates=["condition"])["error"] == "bad_covariates"
    assert tools.run_pseudobulk_de(["condition", "treated", "ctrl"], "donor")["error"] == "bad_level"


def test_scvi_batch_key_list_is_validated(two_condition_h5ad):
    """A list of columns means one batch per combination; unknown columns are refused
    before any training."""
    SESSION.load(two_condition_h5ad)
    tools.check_gene_identifiers()
    tools.compute_qc()
    tools.filter_cells_and_genes(min_genes=5, max_pct_mt=90, min_cells=1)
    tools.normalize(n_top_genes=50)
    result = tools.run_scvi(batch_key=["donor", "no_such_column"])
    assert result["error"] == "invalid_batch_key"
    assert "no_such_column" in result["message"]


def test_relabel_clusters(two_condition_h5ad):
    """Whole clusters get a new label with a reason; the CellTypist label is kept, stale
    comparison results are cleared, and the report shows both labels."""
    SESSION.load(two_condition_h5ad)
    SESSION.begin_run()
    adata = SESSION.adata
    adata.obs["leiden"] = pd.Categorical(np.where(adata.obs["cell_type"] == "Type A", "0", "1"))
    tools.run_pseudobulk_de(["condition", "stim", "ctrl"], sample_key="donor", covariates=["donor"])

    result = tools.relabel_clusters({"0": "CD14+ monocytes"}, reason="LYZ, CD14, S100A8 high")

    assert "error" not in result, result
    assert result["changes"] == [{"cluster": "0", "from": "Type A", "to": "CD14+ monocytes", "n_cells": 400}]
    assert result["cleared"] == ["run_pseudobulk_de"] and not SESSION.de_results
    assert set(adata.obs["cell_type"]) == {"CD14+ monocytes", "Type B"}
    assert set(adata.obs["cell_type_celltypist"]) == {"Type A", "Type B"}

    SESSION.paths.figures.mkdir(parents=True, exist_ok=True)
    assert "relabelled by the agent" in report._composition_figure(adata, SESSION.paths.figures).caption
    table = report.tables([], adata)["clusters"]
    assert "| **CD14+ monocytes** | Type A |" in table

    # Guardrails: unknown clusters, empty labels and missing reasons are refused.
    assert tools.relabel_clusters({"9": "X"}, reason="markers")["error"] == "bad_labels"
    assert tools.relabel_clusters({"1": " "}, reason="markers")["error"] == "bad_labels"
    assert tools.relabel_clusters({"1": "X"}, reason=" ")["error"] == "no_reason"
