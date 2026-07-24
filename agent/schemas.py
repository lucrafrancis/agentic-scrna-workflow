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

# Anthropic tool schema list. Stubbed with two representative entries; the rest are filled
# in alongside the tool implementations.
TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "inspect_dataset",
        "description": "Inspect the AnnData: cell/gene counts, organism, obs/var columns, "
        "candidate batch key, and whether raw counts are available. Call this first.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "run_scvi",
        "description": "Compute an scVI latent embedding from raw counts with batch "
        "correction. Prefer over PCA only when the data has >1 batch.",
        "input_schema": {
            "type": "object",
            "properties": {
                "batch_key": {
                    "type": "string",
                    "description": "obs column identifying batches.",
                }
            },
            "required": [],
        },
    },
    # ... remaining tool schemas added with their implementations.
]

# name -> callable. The loop dispatches through this; it never imports tools directly.
TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "inspect_dataset": tools.inspect_dataset,
    "compute_qc": tools.compute_qc,
    "recommend_qc_thresholds": tools.recommend_qc_thresholds,
    "filter_cells_and_genes": tools.filter_cells_and_genes,
    "detect_doublets": tools.detect_doublets,
    "normalize": tools.normalize,
    "run_pca": tools.run_pca,
    "run_scvi": tools.run_scvi,
    "cluster": tools.cluster,
    "identify_markers": tools.identify_markers,
    "annotate_celltypes": tools.annotate_celltypes,
    "summarize_findings": tools.summarize_findings,
    "generate_report": tools.generate_report,
}
