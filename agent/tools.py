"""Analysis tools.

THE CONTRACT (see CLAUDE.md):
  Every tool is a deterministic Python function that does real work AND returns a
  structured, JSON-serializable summary dict. That dict is the *only* thing the LLM
  sees — it decides the next step entirely from it. A tool that returns None is a bug.

Tools operate on a shared AnnData held by the running session (see loop.py). Each tool
validates the state it expects and returns a clear error summary (not a stack trace) the
agent can react to.

All functions below are stubs — signatures and return contracts only. Implementations
land in the next step.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import scanpy as sc
from scipy import sparse

from agent.session import SESSION

Summary = dict[str, Any]

# obs columns that commonly denote a batch/sample grouping, checked by inspect_dataset.
_BATCH_HINTS = {"batch", "sample", "donor", "patient", "subject", "condition", "dataset"}

# Quantiles reported for every QC metric distribution.
_QUANTILES = {"min": 0.0, "p25": 0.25, "median": 0.5, "p75": 0.75, "p95": 0.95, "p99": 0.99, "max": 1.0}


def _quantile_summary(values) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    return {name: round(float(np.quantile(arr, q)), 4) for name, q in _QUANTILES.items()}


def _looks_like_counts(X) -> bool:
    """Heuristic: True if X appears to be non-negative integer counts (a sample of values)."""
    data = X.data[:1000] if sparse.issparse(X) else np.asarray(X).ravel()[:1000]
    if data.size == 0:
        return True
    return bool(np.all(data >= 0) and np.allclose(data, np.round(data)))


def inspect_dataset() -> Summary:
    """n cells/genes, obs/var columns, candidate batch key, and count/layer state.

    Non-mutating; the agent calls this first to orient itself. Reports enough for the
    dimensionality-reduction decision (batch key) and the counts guardrail (has_raw_counts).
    """
    adata = SESSION.require_adata()
    obs_cols = list(adata.obs.columns)

    candidate_batch_key, n_batches = None, None
    for col in obs_cols:
        if col.lower() in _BATCH_HINTS and adata.obs[col].nunique() > 1:
            candidate_batch_key = col
            n_batches = int(adata.obs[col].nunique())
            break

    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": obs_cols,
        "var_columns": list(adata.var.columns),
        "candidate_batch_key": candidate_batch_key,
        "n_batches": n_batches,
        "has_raw_counts": _looks_like_counts(adata.X),
        "existing_layers": list(adata.layers.keys()),
        "gene_id_sample": [str(g) for g in adata.var_names[:5]],
    }


def check_gene_identifiers() -> Summary:
    """Detect the gene-identifier format and prepare for symbol-dependent steps.

    Runs before compute_qc: classifies var_names as symbol / ensembl / other, guesses
    organism and the mitochondrial prefix, and makes var_names unique. Stores the format
    and mito prefix on the session for compute_qc to use. Non-mutating (var cleanup only).
    """
    adata = SESSION.require_adata()
    adata.var_names_make_unique()
    names = adata.var_names.astype(str)
    n = len(names)

    ensembl = sum(bool(re.match(r"ENS[A-Z]*G\d+", s.upper())) for s in names)
    symbolic = sum(bool(re.match(r"[A-Za-z][A-Za-z0-9._-]*$", s)) for s in names)
    if n and ensembl / n > 0.5:
        gene_format = "ensembl"
    elif n and symbolic / n > 0.5:
        gene_format = "symbol"
    else:
        gene_format = "other"

    # Organism / mito prefix only meaningful when the index is symbols.
    organism, mito_prefix = None, None
    if gene_format == "symbol":
        if any(s.startswith("MT-") for s in names):
            organism, mito_prefix = "human", "MT-"
        elif any(s.startswith("mt-") for s in names):
            organism, mito_prefix = "mouse", "mt-"
        else:
            mito_prefix = "MT-"  # assume human symbols; compute_qc reports if 0 are found

    symbols_available = gene_format == "symbol"
    SESSION.gene_format = gene_format
    SESSION.mito_prefix = mito_prefix

    return {
        "gene_format": gene_format,
        "organism_guess": organism,
        "mito_prefix": mito_prefix,
        "symbols_available": symbols_available,
        "n_genes": int(n),
        "gene_id_sample": [str(g) for g in names[:5]],
        # Option A: flag rather than map. Symbol-dependent steps are degraded, not fixed.
        "note": None
        if symbols_available
        else "No gene symbols in var_names; mito detection and marker interpretation are degraded.",
    }


def compute_qc() -> Summary:
    """Compute per-cell QC (n_genes, total_counts, pct_mito) and return distributions.

    Requires check_gene_identifiers to have run (needs the mito prefix). Annotates obs/var
    with QC metrics; does not remove anything, so it does not checkpoint.
    """
    adata = SESSION.require_adata()
    if SESSION.gene_format is None:
        return {
            "error": "gene_identifiers_not_checked",
            "message": "Call check_gene_identifiers before compute_qc.",
        }

    prefix = SESSION.mito_prefix
    adata.var["mt"] = adata.var_names.str.startswith(prefix) if prefix else False
    n_mito = int(np.asarray(adata.var["mt"]).sum())
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    return {
        "n_mito_genes_found": n_mito,
        "n_genes_by_counts": _quantile_summary(adata.obs["n_genes_by_counts"]),
        "total_counts": _quantile_summary(adata.obs["total_counts"]),
        "pct_counts_mt": _quantile_summary(adata.obs["pct_counts_mt"]),
    }


def recommend_qc_thresholds() -> Summary:
    """Suggest tutorial-default thresholds and their projected impact on THIS dataset.

    Suggestion only, no mutation. The agent decides whether to accept or adjust before
    calling filter_cells_and_genes.
    """
    adata = SESSION.require_adata()
    if "n_genes_by_counts" not in adata.obs:
        return {"error": "qc_not_computed", "message": "Call compute_qc before recommend_qc_thresholds."}

    min_genes, min_cells, max_pct_mt = 200, 3, 5.0
    below_min_genes = int((adata.obs["n_genes_by_counts"] < min_genes).sum())
    above_max_mt = int((adata.obs["pct_counts_mt"] > max_pct_mt).sum())
    keep = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["pct_counts_mt"] <= max_pct_mt)
    genes_below_min_cells = int((adata.var["n_cells_by_counts"] < min_cells).sum())

    return {
        "recommended": {"min_genes": min_genes, "max_pct_mt": max_pct_mt, "min_cells": min_cells},
        "rationale": {
            "min_genes": "standard scanpy tutorial floor for viable cells",
            "max_pct_mt": "PBMC-typical mitochondrial cap; high mito % marks stressed/dying cells",
            "min_cells": "drop genes detected in too few cells (uninformative)",
        },
        "projected_impact": {
            "cells_total": int(adata.n_obs),
            "cells_removed_min_genes": below_min_genes,
            "cells_removed_max_pct_mt": above_max_mt,
            "cells_removed_total": int((~keep).sum()),
            "genes_total": int(adata.n_vars),
            "genes_removed_min_cells": genes_below_min_cells,
        },
    }


def filter_cells_and_genes(min_genes: int = 200, max_pct_mt: float = 5.0, min_cells: int = 3) -> Summary:
    """Apply cell and gene filters, then checkpoint. First mutating tool in the arc."""
    adata = SESSION.require_adata()
    if "n_genes_by_counts" not in adata.obs:
        return {"error": "qc_not_computed", "message": "Call compute_qc before filtering."}

    before = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mt].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    SESSION.adata = adata

    after = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}
    checkpoint = SESSION.checkpoint("after_filter")
    return {
        "applied": {"min_genes": min_genes, "max_pct_mt": max_pct_mt, "min_cells": min_cells},
        "before": before,
        "after": after,
        "removed": {
            "cells": before["n_cells"] - after["n_cells"],
            "genes": before["n_genes"] - after["n_genes"],
        },
        "checkpoint": checkpoint,
    }


def detect_doublets() -> Summary:
    """Scrublet via sc.pp.scrublet. Returns predicted doublet rate. Run before normalization."""
    raise NotImplementedError


def normalize() -> Summary:
    """Stash raw counts to layers['counts'], then log-normalize + HVGs. Guardrail: counts first."""
    raise NotImplementedError


def run_pca(n_comps: int = 50) -> Summary:
    """PCA on the normalized matrix. Correct choice for a single clean batch."""
    raise NotImplementedError


def run_scvi(batch_key: str | None = None) -> Summary:
    """scVI latent from raw counts (layers['counts']), with batch correction. Not for single batch."""
    raise NotImplementedError


def cluster(resolution: float = 1.0) -> Summary:
    """Neighbors + Leiden on the chosen representation (PCA or scVI). Returns n clusters, sizes."""
    raise NotImplementedError


def identify_markers() -> Summary:
    """Rank genes per cluster. Returns top markers per cluster."""
    raise NotImplementedError


def annotate_celltypes() -> Summary:
    """CellTypist annotation. Returns per-cluster majority cell type + confidence."""
    raise NotImplementedError


def summarize_findings() -> Summary:
    """Assemble the biological narrative from prior tool outputs."""
    raise NotImplementedError


def generate_report() -> Summary:
    """Write outputs/report.md (+ figures). Returns the report path."""
    raise NotImplementedError
