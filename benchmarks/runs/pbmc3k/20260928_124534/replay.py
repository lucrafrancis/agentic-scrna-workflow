"""Replay of the agent's tool calls for pbmc3k, without the LLM.

Generated: 2026-09-28T12:47:39
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

_run("identify_markers")

_run("annotate_celltypes")

_run(
    "check_markers",
    genes=['CD3D', 'CD3E', 'CD8A', 'CD8B', 'GZMA', 'NKG7', 'CCL5', 'CD79A', 'MS4A1', 'CD79B',
     'IL7R', 'CD4', 'CD14', 'LYZ', 'S100A9', 'FCGR3A', 'FCER1A', 'CD1C', 'HLA-DRA', 'GNLY',
     'NCAM1', 'MS4A7', 'PPBP', 'PF4', 'ITGA2B'],
    clusters=['0', '1', '2', '3', '4', '5', '6', '7', '8'],
)

_run("summarize_findings")

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This report covers an end-to-end analysis of a single 10x Genomics PBMC sample from '
     'one healthy donor (starting with {{input.n_cells}} cells and {{input.n_genes}} '
     'genes, raw counts). Because this is a single, clean 10x run with no batch structure, '
     'the workflow used PCA (not batch-corrected scVI) for dimensionality reduction, and '
     'no cross-condition comparison was performed (there is only one donor and one '
     'condition).\n'
     '\n'
     '## Quality control\n'
     '\n'
     'Gene identifiers were confirmed as human gene symbols, letting mitochondrial genes '
     '({{qc.n_mito_genes}} genes with the `MT-` prefix) be identified correctly for QC. '
     'Per-cell QC distributions showed a median of {{qc.genes.median}} genes and '
     '{{qc.counts.median}} counts per cell, with mitochondrial fraction at a median of '
     '{{qc.mito.median}} (p95 {{qc.mito.p95}}, max {{qc.mito.max}}).\n'
     '\n'
     'Standard tutorial-default thresholds fit this dataset well and were applied '
     'unchanged: minimum {{filter.min_genes}} genes/cell (removed no cells — the '
     'population was already well above this floor), maximum {{filter.max_pct_mt}} '
     'mitochondrial content (removed {{default.mito_cells_removed}} cells, '
     '{{default.mito_pct_removed}} of the dataset, trimming the long tail of '
     'stressed/dying cells beyond the 99th percentile), and genes detected in at least '
     '{{filter.min_cells}} cells (dropping {{filter.genes_removed}} uninformative genes). '
     'After filtering, {{filter.cells_after}} cells and {{filter.genes_after}} genes '
     'remained.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'Scrublet was run on the whole sample as one batch (single 10x run, no '
     'genotype/hashing demultiplexing reported, so no light-touch exception applies). The '
     'doublet-score distribution was unimodal and decayed smoothly with no bimodal valley, '
     'so the median+3×MAD rule was used to set the threshold rather than the '
     '(inapplicable) bimodal-valley or the more permissive Scrublet-automatic option. This '
     'gave a threshold of {{doublet.threshold}} (median {{doublet.score_median}}, max '
     '{{doublet.score_max}}), flagging and removing {{doublet.cells_removed}} cells '
     '({{doublet.pct_removed}}) as likely doublets — a plausible rate for this loading '
     'density. {{final.n_cells}} cells and {{final.n_genes}} genes entered downstream '
     'analysis.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'No batch key was detected and the sample is a single clean 10x run from one donor, '
     'so PCA on the {{norm.n_hvgs}} highly variable genes (after total-count normalization '
     "to {{norm.target_sum}} and log1p) was the appropriate, simplest choice — scVI's "
     'batch correction would add complexity with no batches to correct. The leading '
     'principal components captured {{pca.pc1_variance}}, {{pca.pc2_variance}}, and '
     '{{pca.pc3_variance}} of variance respectively, consistent with a few dominant axes '
     'of lineage variation (myeloid vs lymphoid, then finer distinctions) typical of '
     'PBMCs.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering at resolution {{cluster.resolution}} on the PCA embedding yielded '
     '{{cluster.n}} clusters, ranging in size from {{cluster.smallest}} to '
     '{{cluster.largest}} cells — sizes and counts consistent with the expected mixture of '
     'abundant T-cell/monocyte populations and rarer populations (DCs, platelets).\n'
     '\n'
     '{{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, fine-grained) was used for initial labels, with '
     'Immune_All_High.pkl (coarse) as a second opinion. The two models agreed at the '
     'expected level of granularity (e.g., the fine model split lymphocytes into Tem/Trm '
     'cytotoxic T cells, Tcm/Naive helper T cells, and CD16+ NK cells, while the coarse '
     'model grouped the first two as "T cells" and the NK cluster as the broader "ILC" '
     'category that includes NK cells — not a contradiction).\n'
     '\n'
     'Every final cell type was checked against canonical markers before accepting the '
     'labels:\n'
     '- **Tem/Trm cytotoxic T cells** (cluster 0): high CD3D, CD3E, CD8A, CD8B, GZMA, NKG7 '
     'and CCL5, consistent with cytotoxic CD8 T cells.\n'
     '- **B cells** (cluster 1): high CD79A, MS4A1 and CD79B, with CD3D essentially '
     'absent.\n'
     '- **Tcm/Naive helper T cells** (clusters 2 and 3): CD3D/CD3E-positive with high IL7R '
     'and low CD8A, consistent with CD4 T cells; cluster 3 is dominated by '
     'ribosomal-protein genes, typical of quiescent naive lymphocytes.\n'
     '- **Classical monocytes** (cluster 4): near-universal CD14, LYZ and S100A9 '
     'expression.\n'
     '- **DC** (cluster 5): high HLA-DRA together with FCER1A and CD1C, both significantly '
     'enriched and largely absent elsewhere — a myeloid dendritic cell signature.\n'
     '- **CD16+ NK cells** (cluster 6): near-universal GNLY and NKG7, high FCGR3A, and '
     'CD3D essentially absent, ruling out a T/NKT identity.\n'
     '- **Non-classical monocytes** (cluster 7): high FCGR3A and MS4A7 with comparatively '
     'low CD14, the classical non-classical/CD14-dim monocyte profile.\n'
     '- **Megakaryocytes/platelets** (cluster 8): near-universal PPBP, PF4 and ITGA2B '
     '(CD41), a clean platelet signature.\n'
     '\n'
     "Because all clusters' markers matched their CellTypist labels (and the "
     "second-opinion model's differences were only granularity, not contradiction), no "
     'cluster relabeling was necessary.\n'
     '\n'
     '{{table:composition}}\n'
     '\n'
     'Final composition: {{celltype:Tcm/Naive helper T cells}} CD4 T cells, '
     '{{celltype:Classical monocytes}} classical monocytes, {{celltype:B cells}} B cells, '
     '{{celltype:Tem/Trm cytotoxic T cells}} cytotoxic T cells, {{celltype:Non-classical '
     'monocytes}} non-classical monocytes, {{celltype:CD16+ NK cells}} NK cells, '
     '{{celltype:DC}} dendritic cells, and {{celltype:Megakaryocytes/platelets}} '
     'megakaryocytes/platelets — proportions broadly in line with expectations for a '
     'healthy PBMC sample.\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- This is a single donor/single run, so no batch correction or condition comparison '
     "was needed or performed; results reflect this one individual's PBMC composition and "
     'cannot be generalized without replication.\n'
     '- The doublet threshold (median+3×MAD) is a heuristic; a small number of true '
     'rare-cell profiles that resemble doublets transcriptionally (e.g. some large '
     'cytotoxic cells) could in principle be lost, though the smooth, unimodal score '
     'distribution gives no indication of a distinct doublet population being missed or '
     'over-cut.\n'
     '- The megakaryocyte/platelet cluster is very small '
     '({{celltype:Megakaryocytes/platelets}} cells), so its marker statistics should be '
     'interpreted cautiously; individual DE p-values within this cluster are less well '
     'powered.\n'
     "- Cluster 3's marker list is dominated by ribosomal genes rather than a distinctive "
     'positive signature; its T-cell identity rests mainly on CD3D/CD3E/IL7R positivity '
     'confirmed via check_markers rather than on unique cluster-defining genes.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'Standard QC and doublet filtering removed a modest fraction of low-quality/multiplet '
     'cells ({{final.pct_removed}} total), leaving a clean dataset of {{final.n_cells}} '
     'cells. PCA-based clustering at resolution {{cluster.resolution}} resolved the '
     'expected major PBMC lineages, and CellTypist annotation — cross-validated against '
     'canonical markers for every cluster — produced confident, mutually consistent labels '
     'for all {{celltype.n}} cell types without requiring any manual relabeling.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
