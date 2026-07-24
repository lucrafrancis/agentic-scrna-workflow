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
