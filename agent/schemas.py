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

# Anthropic tool schema list: the complete set of actions available to the agent.
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
                "batch_key": {
                    "anyOf": [{"type": "string"}, {"type": "array", "items": {"type": "string"}}],
                    "description": "obs column identifying batches to correct, or a list of columns "
                    "combined into one batch per combination (e.g. ['donor', 'condition']).",
                },
                "max_epochs": {"type": "integer", "description": "Training epochs (default: scVI auto)."},
            },
            "required": [],
        },
    },
    {
        "name": "cluster",
        "description": "Build a neighbor graph on the chosen representation, run Leiden "
        "clustering, and compute UMAP coordinates, then checkpoint. Requires a dimensionality "
        "reduction (run_pca or run_scvi) first. Higher resolution yields more clusters.",
        "input_schema": {
            "type": "object",
            "properties": {
                "resolution": {"type": "number", "description": "Leiden resolution (default 1.0); higher = more clusters."}
            },
            "required": [],
        },
    },
    {
        "name": "identify_markers",
        "description": "Rank marker genes per cluster (Wilcoxon) and return the top genes per "
        "cluster, for interpreting cluster identity. Requires clustering first.",
        "input_schema": {
            "type": "object",
            "properties": {
                "n_genes": {"type": "integer", "description": "Top marker genes to return per cluster (default 10)."}
            },
            "required": [],
        },
    },
    {
        "name": "annotate_celltypes",
        "description": "Annotate cell types with CellTypist (majority voting over Leiden "
        "clusters). Returns per-cluster cell-type labels and overall counts. Requires "
        "clustering first. Default model suits immune/PBMC data.",
        "input_schema": {
            "type": "object",
            "properties": {
                "model": {"type": "string", "description": "CellTypist model (default Immune_All_Low.pkl)."}
            },
            "required": [],
        },
    },
    {
        "name": "compare_composition",
        "description": "Compare cell-type proportions between two conditions across replicate "
        "samples (e.g. donors). Proportions are computed per sample; the test is paired "
        "(Wilcoxon signed-rank) when every sample has both conditions, else Mann-Whitney U, with "
        "Benjamini-Hochberg correction. Refuses with too few samples per condition. Run after "
        "annotate_celltypes. Non-mutating.",
        "input_schema": {
            "type": "object",
            "properties": {
                "contrast": {
                    "type": "array", "items": {"type": "string"}, "minItems": 3, "maxItems": 3,
                    "description": "[condition_column, test_level, reference_level], e.g. ['condition', 'stim', 'ctrl'].",
                },
                "sample_key": {"type": "string", "description": "obs column identifying replicates, e.g. 'donor'."},
                "celltype_key": {"type": "string", "description": "obs column with cell types (default 'cell_type')."},
            },
            "required": ["contrast", "sample_key"],
        },
    },
    {
        "name": "run_pseudobulk_de",
        "description": "Pseudobulk differential expression per cell type with PyDESeq2 (Wald "
        "test). Raw counts are summed per sample (sample_key x condition) within each cell "
        "type; design = ~ covariates + factor. Put the replicate column (e.g. 'donor') in "
        "covariates for a paired design when each replicate has both conditions. Refuses "
        "confounded designs and covariates that vary within a sample; skips cell types with "
        "fewer than 2 samples per condition, and says why. Run after annotate_celltypes.",
        "input_schema": {
            "type": "object",
            "properties": {
                "contrast": {
                    "type": "array", "items": {"type": "string"}, "minItems": 3, "maxItems": 3,
                    "description": "[factor, test_level, reference_level], e.g. ['condition', 'stim', 'ctrl'].",
                },
                "sample_key": {"type": "string", "description": "obs column identifying replicates, e.g. 'donor'."},
                "covariates": {
                    "type": "array", "items": {"type": "string"},
                    "description": "Categorical columns to adjust for, e.g. ['donor'] (paired) or ['batch'].",
                },
                "celltype_key": {"type": "string", "description": "obs column with cell types (default 'cell_type')."},
                "cell_types": {
                    "type": "array", "items": {"type": "string"},
                    "description": "Cell types to test (default: all).",
                },
                "min_cells": {
                    "type": "integer",
                    "description": "Minimum cells for a sample to enter a cell type's pseudobulk (default 10).",
                },
            },
            "required": ["contrast", "sample_key"],
        },
    },
    {
        "name": "get_top_genes",
        "description": "Top significant DE genes for one tested cell type, ranked by padj then "
        "|log2FC|. Requires run_pseudobulk_de.",
        "input_schema": {
            "type": "object",
            "properties": {
                "cell_type": {"type": "string", "description": "A cell type tested by run_pseudobulk_de."},
                "n": {"type": "integer", "description": "Number of genes (default 20)."},
                "direction": {"type": "string", "enum": ["up", "down", "both"], "description": "Default 'both'."},
            },
            "required": ["cell_type"],
        },
    },
    {
        "name": "query_genes",
        "description": "log2FC and padj for specific genes in each tested cell type. Use it to "
        "check a gene before writing about it. Requires run_pseudobulk_de.",
        "input_schema": {
            "type": "object",
            "properties": {
                "genes": {"type": "array", "items": {"type": "string"}, "description": "Gene symbols."},
                "cell_types": {
                    "type": "array", "items": {"type": "string"},
                    "description": "Cell types to report (default: all tested).",
                },
            },
            "required": ["genes"],
        },
    },
    {
        "name": "summarize_findings",
        "description": "Consolidate the final analysis state (cell counts, cluster-to-"
        "cell-type mapping, top markers, cell-type counts) into one factual summary to write "
        "the report from. Also returns 'facts' (placeholder names with their current values, "
        "for {{name}} in the report) and the tables available as {{table:name}}. Non-mutating.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "generate_report",
        "description": "Write the final Markdown report. Never type a number: use {{fact}}, "
        "{{celltype:NAME}}, {{celltypes:A|B}} and {{table:name}} placeholders, which the tool "
        "fills in. Any problem rejects the whole report with a list of all issues; fix them "
        "all and resubmit. The tool adds the decisions table, captioned figures and Methods, "
        "and writes the annotated .h5ad. Call this last.",
        "input_schema": {
            "type": "object",
            "properties": {
                "report_markdown": {
                    "type": "string",
                    "description": "The full analysis report as Markdown, written by you from "
                    "the findings, with placeholders in place of every number.",
                }
            },
            "required": ["report_markdown"],
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
    "compare_composition": tools.compare_composition,
    "run_pseudobulk_de": tools.run_pseudobulk_de,
    "get_top_genes": tools.get_top_genes,
    "query_genes": tools.query_genes,
    "summarize_findings": tools.summarize_findings,
    "generate_report": tools.generate_report,
}
