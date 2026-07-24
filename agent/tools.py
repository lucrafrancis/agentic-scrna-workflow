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

from typing import Any

Summary = dict[str, Any]


def inspect_dataset() -> Summary:
    """n cells/genes, organism, obs/var columns, candidate batch key, count layer state."""
    raise NotImplementedError


def compute_qc() -> Summary:
    """Compute per-cell QC (n_genes, total_counts, pct_mito, ...) and return distributions."""
    raise NotImplementedError


def recommend_qc_thresholds() -> Summary:
    """Suggest filtering thresholds from the QC distributions. Suggestion only, no mutation."""
    raise NotImplementedError


def filter_cells_and_genes(min_genes: int, max_pct_mito: float, min_cells: int) -> Summary:
    """Apply filtering. Returns before/after cell & gene counts."""
    raise NotImplementedError


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
