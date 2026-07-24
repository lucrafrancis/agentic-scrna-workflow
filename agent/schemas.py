"""JSON tool definitions sent to Claude, plus the name -> function dispatch table.

The schemas here MUST stay in sync with the signatures in tools.py: the `input_schema`
describes the arguments Claude is allowed to pass, and TOOL_FUNCTIONS maps the tool name
back to the Python function the loop actually calls.

Kept deliberately in its own module so the loop stays generic — it never hard-codes tool
names, it just iterates over whatever is registered here.
"""

from __future__ import annotations

from typing import Any, Callable

from agent import tools

# Anthropic tool schema list. We expose ONLY implemented tools so the agent cannot call a
# stub. Entries are added here as each tool is implemented.
TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "inspect_dataset",
        "description": "Inspect the AnnData: cell/gene counts, obs/var columns, candidate "
        "batch key, and whether raw counts are available. Call this first.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "check_gene_identifiers",
        "description": "Detect the gene-identifier format (symbol / ensembl / other), guess "
        "organism and the mitochondrial prefix, and make gene names unique. Run this before "
        "compute_qc, which relies on the detected mito prefix.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "compute_qc",
        "description": "Compute per-cell QC metrics (n_genes, total_counts, pct mitochondrial) "
        "and return their distributions. Requires check_gene_identifiers first.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "recommend_qc_thresholds",
        "description": "Suggest tutorial-default QC thresholds and their projected impact on "
        "this dataset. Suggestion only. Requires compute_qc first.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "filter_cells_and_genes",
        "description": "Filter cells (min genes per cell, max mitochondrial %) and genes "
        "(min cells per gene), then checkpoint. Requires compute_qc first.",
        "input_schema": {
            "type": "object",
            "properties": {
                "min_genes": {"type": "integer", "description": "Minimum genes per cell to keep it."},
                "max_pct_mt": {"type": "number", "description": "Maximum mitochondrial % to keep a cell."},
                "min_cells": {"type": "integer", "description": "Minimum cells a gene must appear in."},
            },
            "required": ["min_genes", "max_pct_mt", "min_cells"],
        },
    },
    {
        "name": "detect_doublets",
        "description": "Run Scrublet and report the doublet-score distribution plus candidate "
        "thresholds (bimodal valley, median+3*MAD, Scrublet automatic) and a recommendation. "
        "Non-mutating. Inspect the distribution, then choose a threshold for filter_doublets.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "filter_doublets",
        "description": "Remove cells whose doublet_score is >= the chosen threshold, then "
        "checkpoint. Requires detect_doublets first. Explain why you chose this threshold.",
        "input_schema": {
            "type": "object",
            "properties": {
                "threshold": {"type": "number", "description": "Doublet-score cutoff; cells at or above are removed."}
            },
            "required": ["threshold"],
        },
    },
    {
        "name": "normalize",
        "description": "Stash raw counts, then total-count normalize, log1p, and flag highly "
        "variable genes, then checkpoint. Run after doublet handling and before dimensionality "
        "reduction.",
        "input_schema": {
            "type": "object",
            "properties": {
                "target_sum": {"type": "number", "description": "Per-cell count target (default 1e4)."},
                "n_top_genes": {"type": "integer", "description": "Number of HVGs to flag (default 2000)."},
            },
            "required": [],
        },
    },
    {
        "name": "run_pca",
        "description": "PCA on the normalized HVG matrix; sets the representation used for "
        "clustering. The right choice for a single clean batch (no batch key). Requires "
        "normalize first.",
        "input_schema": {
            "type": "object",
            "properties": {
                "n_comps": {"type": "integer", "description": "Number of principal components (default 50)."}
            },
            "required": [],
        },
    },
    {
        "name": "run_scvi",
        "description": "Train scVI on raw counts (HVGs) to get a latent embedding with batch "
        "correction; sets the representation used for clustering. Prefer over PCA only when a "
        "batch key spans multiple batches. Requires normalize first.",
        "input_schema": {
            "type": "object",
            "properties": {
                "batch_key": {"type": "string", "description": "obs column identifying batches to correct."},
                "max_epochs": {"type": "integer", "description": "Training epochs (default: scVI auto)."},
            },
            "required": [],
        },
    },
]

# name -> callable. The loop dispatches through this; it never imports tools directly.
TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "inspect_dataset": tools.inspect_dataset,
    "check_gene_identifiers": tools.check_gene_identifiers,
    "compute_qc": tools.compute_qc,
    "recommend_qc_thresholds": tools.recommend_qc_thresholds,
    "filter_cells_and_genes": tools.filter_cells_and_genes,
    "detect_doublets": tools.detect_doublets,
    "filter_doublets": tools.filter_doublets,
    "normalize": tools.normalize,
    "run_pca": tools.run_pca,
    "run_scvi": tools.run_scvi,
    "cluster": tools.cluster,
    "identify_markers": tools.identify_markers,
    "annotate_celltypes": tools.annotate_celltypes,
    "summarize_findings": tools.summarize_findings,
    "generate_report": tools.generate_report,
}
