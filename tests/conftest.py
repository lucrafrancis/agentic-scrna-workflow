"""Test fixtures.

Everything here runs offline: a tiny synthetic AnnData stands in for real data so the
whole pipeline can be exercised in seconds without downloads, the API, or scVI training.
Outputs are redirected to a temp dir and the SESSION singleton is reset between tests.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from agent import config
from agent.session import SESSION


def make_adata(n_cells=400, n_genes=300, with_batch=False, duplicate_barcodes=False, seed=0):
    """A small raw-counts AnnData with two clear populations and a few MT- genes.

    Two distinct marker blocks make clustering find >=2 groups; MT- genes let the mito path
    run. Values are Poisson counts so has_raw_counts is True.
    """
    rng = np.random.default_rng(seed)
    X = rng.poisson(1.0, size=(n_cells, n_genes)).astype(np.float32)
    half = n_cells // 2
    X[:half, 0:30] += rng.poisson(8.0, size=(half, 30))  # population A markers
    X[half:, 30:60] += rng.poisson(8.0, size=(n_cells - half, 30))  # population B markers

    var_names = [f"GENE{i}" for i in range(n_genes)]
    for i in range(5):
        var_names[60 + i] = f"MT-{i}"  # mitochondrial genes, outside the marker blocks

    if duplicate_barcodes:
        obs_names = [f"cell{i}" for i in range(half)] + [f"cell{i}" for i in range(n_cells - half)]
    else:
        obs_names = [f"cell{i}" for i in range(n_cells)]

    adata = ad.AnnData(
        X=sparse.csr_matrix(X),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    if with_batch:
        adata.obs["batch"] = pd.Categorical(["A"] * half + ["B"] * (n_cells - half))
    return adata


@pytest.fixture
def synthetic_h5ad(tmp_path):
    path = tmp_path / "synthetic.h5ad"
    make_adata().write_h5ad(path)
    return path


@pytest.fixture
def batched_h5ad(tmp_path):
    path = tmp_path / "batched.h5ad"
    make_adata(with_batch=True).write_h5ad(path)
    return path


@pytest.fixture
def dup_barcode_h5ad(tmp_path):
    path = tmp_path / "dup.h5ad"
    make_adata(duplicate_barcodes=True).write_h5ad(path)
    return path


@pytest.fixture(autouse=True)
def isolate(tmp_path, monkeypatch):
    """Redirect outputs to a temp dir and reset the SESSION singleton for each test."""
    monkeypatch.setattr(config, "OUTPUT_DIR", tmp_path / "outputs")
    SESSION.__init__()
    yield


@pytest.fixture
def logged_run(synthetic_h5ad):
    """Run the pipeline up to the report with every call dispatched and logged by the loop,
    as in a real run, so report code that reads tool_calls.jsonl sees a real log. Returns
    the `call` helper for further tool calls."""
    from agent.loop import _log_tool_call, _run_tool

    SESSION.load(synthetic_h5ad)
    SESSION.begin_run()

    def call(name, **args):
        summary = _run_tool(name, args)
        _log_tool_call(name, args, summary)
        return summary

    call("inspect_dataset")
    call("check_gene_identifiers")
    call("compute_qc")
    call("recommend_qc_thresholds")
    call("filter_cells_and_genes", min_genes=5, max_pct_mt=90, min_cells=1)
    threshold = call("detect_doublets")["recommended_threshold"]
    call("filter_doublets", threshold=threshold)
    call("normalize", n_top_genes=50)
    call("run_pca", n_comps=10)
    call("cluster", resolution=1.0)
    call("identify_markers", n_genes=5)
    return call
