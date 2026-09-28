"""Replay of the agent's tool calls for pbmc3k, without the LLM.

Generated: 2026-09-28T12:45:33
Git commit: 52f5efd97eb1bf4a6f9a346598659f3e505c0dc9
Source log: tool_calls.jsonl (15 of 16 calls; failed calls omitted)

Re-runs the same tool functions with the same arguments and seed. Output goes to
outputs/pbmc3k_replay/, with its own tool log and report. Package versions are pinned by
uv.lock at the commit above: check it out first if the code has changed since.

Caveats: annotate_celltypes downloads CellTypist models, and scVI training can differ
slightly between hardware (CPU vs GPU), which can change clusters downstream.

  uv run python outputs/pbmc3k/replay.py
"""

import os
import sys
from pathlib import Path

# The repo this script sits in (it may have been copied, e.g. into benchmarks/runs/), else
# the repo it was generated in.
ROOT = next((p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists()),
            Path('/Users/lucfra/scratch/agentic-scrna-workflow'))
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)  # the dataset path below is relative to the repo root

from agent import config, tools
from agent.loop import _log_tool_call, _run_tool
from agent.session import SESSION


def _run(name, **args):
    print(f"-> {name}")
    summary = _run_tool(name, args)
    _log_tool_call(name, args, summary)
    if "error" in summary:
        sys.exit(f"   FAILED: {summary['error']}: {summary.get('message', '')}")
    return summary


config.set_global_seed()
SESSION.load('data/pbmc3k.h5ad')
SESSION.paths = config.RunPaths('pbmc3k_replay')
SESSION.begin_run()
print(f"Replay output: {SESSION.paths.dir}")

_run("inspect_dataset")

_run("check_gene_identifiers")

_run("compute_qc")

_run("recommend_qc_thresholds")

_run(
    "filter_cells_and_genes",
    min_genes=200,
    max_pct_mt=5.0,
    min_cells=3,
)

_run(
    "detect_doublets",
    batch_key=None,
)

_run(
    "filter_doublets",
    threshold=0.1,
)

_run("normalize")

_run("run_pca")

_run(
    "cluster",
    resolution=1.0,
)

_run(
    "identify_markers",
    n_genes=25,
)

_run(
    "annotate_celltypes",
    model='Immune_All_Low.pkl',
)

_run(
    "check_markers",
    genes=['CD3D', 'CD3E', 'CD8A', 'CD8B', 'CD4', 'IL7R', 'CCR7', 'GZMK', 'NKG7', 'GNLY',
     'FCGR3A', 'NCAM1', 'MS4A1', 'CD79A', 'CD19', 'LYZ', 'CD14', 'FCN1', 'MS4A7', 'FCER1A',
     'CD1C', 'CLEC9A', 'PPBP', 'PF4', 'S100A8'],
)

_run("summarize_findings")

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This dataset comprises PBMCs from a single healthy donor, processed on a single 10x '
     'Genomics run — a "clean" single-batch scenario with no technical replicate structure '
     'to correct for. The analysis pipeline proceeded from raw counts through QC, doublet '
     'removal, normalization, PCA-based clustering, and CellTypist-assisted annotation. '
     'Starting from {{input.n_cells}} cells and {{input.n_genes}} genes, the final '
     'annotated dataset contains {{final.n_cells}} cells and {{final.n_genes}} genes '
     '({{final.pct_removed}} of the original cells removed across QC and doublet filtering '
     'steps).\n'
     '\n'
     '## Quality control\n'
     '\n'
     'Gene identifiers were confirmed as human gene symbols, with the mitochondrial prefix '
     '"MT-" identifying {{qc.n_mito_genes}} mitochondrial genes. Per-cell QC metrics '
     'showed a healthy, relatively narrow distribution: median genes/cell of '
     '{{qc.genes.median}} (range {{qc.genes.min}}–{{qc.genes.max}}) and median total '
     'counts of {{qc.counts.median}}. Mitochondrial fraction was modest overall (median '
     '{{qc.mito.median}}, p95 {{qc.mito.p95}}) with a long tail up to {{qc.mito.max}}, '
     'consistent with a small population of stressed or dying cells rather than widespread '
     'ambient contamination.\n'
     '\n'
     'I used the standard tutorial-default thresholds ({{filter.min_genes}} min '
     'genes/cell, {{filter.max_pct_mt}} max mitochondrial fraction, {{filter.min_cells}} '
     "min cells/gene), since they matched this dataset's own percentiles well — the mito "
     'cutoff sits just above the p95 mitochondrial percentage, targeting genuinely '
     'high-mito outlier cells rather than cutting into the bulk distribution. This removed '
     '{{filter.cells_removed}} cells ({{filter.pct_cells_removed}}) and '
     '{{filter.genes_removed}} genes not detected in enough cells, leaving '
     '{{filter.cells_after}} cells and {{filter.genes_after}} genes.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'Scrublet was run on the full dataset as a single batch (one 10x run, no batches to '
     'split by). The doublet-score distribution was **not bimodal** — there was no clear '
     'valley separating singlets from doublets, which is common for a real (not '
     'synthetic-boosted) PBMC run at moderate loading density. Per the stated rule, I '
     'therefore used **median + 3×MAD** ({{doublet.threshold}}) rather than a fixed or '
     "bimodal-valley cutoff. This is more conservative than Scrublet's own automatic "
     'threshold ({{doublet.scrublet_auto}}, which would have flagged only '
     '{{doublet.scrublet_auto_cells_flagged}} cells), but appropriate here since there is '
     'no evidence of upstream doublet removal (no genotype demultiplexing or hashing was '
     'reported for this single-donor run). This removed {{doublet.cells_removed}} cells '
     '({{doublet.pct_removed}}), consistent with expected multiplet rates for a run of '
     'this size.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'With only one donor and one sequencing run, there is no batch structure to correct '
     'for. PCA on the normalized, highly-variable-gene matrix ({{norm.n_hvgs}} HVGs, '
     "{{pca.n_comps}} components) is the simpler and more appropriate choice over scVI's "
     'batch-correcting latent space, which would offer no benefit here and risks removing '
     'genuine biological variation. The top PC explains {{pca.pc1_variance}} of variance, '
     'consistent with the dominant lymphoid/myeloid axis expected in PBMCs.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering at resolution {{cluster.resolution}} on the PCA embedding yielded '
     '{{cluster.n}} clusters, ranging from {{cluster.smallest}} to {{cluster.largest}} '
     'cells, mapping onto the expected major PBMC lineages (T cell subsets, B cells, '
     'monocyte subsets, NK cells, dendritic cells and platelets).\n'
     '\n'
     '{{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, fine-grained immune model) assigned each cluster a '
     'label by majority vote, cross-checked against the coarser Immune_All_High.pkl model '
     'as a second opinion; the two models agreed at their respective resolutions with no '
     'contradictions. I verified every final cell type against canonical markers using '
     '`check_markers`:\n'
     '\n'
     '- **Tem/Trm cytotoxic T cells** (cluster 0): high CD3D/CD3E, CD8A/CD8B, and strong '
     'GZMK, consistent with an effector/memory CD8+ T cell phenotype.\n'
     '- **B cells** (cluster 1): near-universal CD79A and MS4A1 expression, with CD3D '
     'essentially absent, a clean B cell signature.\n'
     '- **Tcm/Naive helper T cells** (clusters 2 and 3): both express CD3D/CD3E and IL7R '
     'strongly; cluster 3 additionally shows strong CCR7 enrichment, marking it as the '
     'more naive-like subset, while cluster 2 is the more central-memory-like counterpart '
     '— both correctly grouped under the same CellTypist label since the model does not '
     'further split naive vs. central memory CD4 T cells.\n'
     '- **Classical monocytes** (cluster 4): near-universal LYZ, CD14, FCN1 and S100A8 '
     'expression, the canonical CD14+ monocyte signature.\n'
     '- **DC** (cluster 5): FCER1A and CD1C both strongly enriched, marking conventional '
     'dendritic cells, distinct from the monocyte clusters.\n'
     '- **CD16+ NK cells** (cluster 6): near-universal GNLY and NKG7 with strong FCGR3A '
     '(CD16), and CD3D/CD3E essentially absent, ruling out a T/NKT identity.\n'
     '- **Non-classical monocytes** (cluster 7): near-universal FCGR3A and MS4A7 with LYZ, '
     'but low CD14, the classic CD14(low)CD16+ non-classical monocyte pattern that '
     'distinguishes it from cluster 4.\n'
     '- **Megakaryocytes/platelets** (cluster 8): unanimous PPBP and PF4 expression, an '
     'unambiguous platelet signature in this small cluster.\n'
     '\n'
     "No relabeling was necessary — every cluster's canonical markers matched its assigned "
     'label and the second-opinion model was consistent.\n'
     '\n'
     '{{table:composition}}\n'
     '\n'
     'Final cell-type composition: {{celltype:Tcm/Naive helper T cells}} Tcm/Naive helper '
     'T cells, {{celltype:Classical monocytes}} classical monocytes, {{celltype:B cells}} '
     'B cells, {{celltype:Tem/Trm cytotoxic T cells}} cytotoxic T cells, '
     '{{celltype:Non-classical monocytes}} non-classical monocytes, {{celltype:CD16+ NK '
     'cells}} CD16+ NK cells, {{celltype:DC}} dendritic cells, and '
     '{{celltype:Megakaryocytes/platelets}} megakaryocytes/platelets — proportions broadly '
     'typical of a healthy PBMC sample.\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- This is a single donor processed in a single run, with no biological or technical '
     'replicates, so no composition or differential-expression comparison across '
     'conditions was performed (none was requested, and none would be statistically '
     'meaningful without replicate samples).\n'
     '- Clusters 2 and 3 are grouped under one CellTypist label (Tcm/Naive helper T cells) '
     'despite showing a CCR7 gradient suggestive of naive vs. central-memory substructure; '
     'a higher clustering resolution or manual sub-gating could resolve this further if '
     'finer T cell subsetting is of interest.\n'
     '- The megakaryocyte/platelet cluster is very small '
     '({{celltype:Megakaryocytes/platelets}} cells), so its marker statistics, while clear '
     '(unanimous PPBP/PF4), rest on limited cell numbers.\n'
     '- Median+3×MAD doublet filtering is a heuristic; a small number of true doublets may '
     'remain and a small number of borderline singlets may have been removed.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'Standard QC and doublet thresholds were well matched to this clean, single-donor '
     'PBMC dataset, requiring only modest filtering. PCA-based clustering (appropriate '
     'given the absence of batch structure) resolved the expected major PBMC lineages, and '
     'CellTypist annotation—confirmed against canonical markers for every cluster—produced '
     'a biologically coherent, unambiguous cell-type map spanning T cell subsets, B cells, '
     'monocyte subsets, NK cells, dendritic cells, and platelets.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
