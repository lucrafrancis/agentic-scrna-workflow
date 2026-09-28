"""Replay of the agent's tool calls for pbmc3k, without the LLM.

Generated: 2026-09-28T12:49:42
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
    genes=['CD3D', 'CD3E', 'CD8A', 'CD8B', 'IL7R', 'CD4', 'CD79A', 'MS4A1', 'CD14', 'LYZ',
     'S100A8', 'FCGR3A', 'MS4A7', 'FCER1A', 'CD1C', 'GNLY', 'NKG7', 'NCAM1', 'PF4', 'PPBP',
     'GZMK', 'CCL5'],
)

_run("summarize_findings")

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This analysis processed a 10x Genomics PBMC dataset from a single healthy donor '
     '({{input.n_cells}} cells, {{input.n_genes}} genes at input) through QC filtering, '
     'doublet removal, normalization, PCA-based dimensionality reduction, Leiden '
     'clustering, and CellTypist-assisted cell-type annotation. As this is a single clean '
     '10x run from one donor with no batch structure, no batch correction was required and '
     'standard PCA was used for the embedding.\n'
     '\n'
     '## Quality control\n'
     '\n'
     'Gene identifiers were confirmed as human gene symbols with a standard `MT-` '
     'mitochondrial prefix ({{qc.n_mito_genes}} mitochondrial genes detected). Per-cell QC '
     'distributions showed a healthy library: median genes/cell of {{qc.genes.median}} '
     '(range {{qc.genes.min}}–{{qc.genes.max}}), median total counts of '
     '{{qc.counts.median}}, and median mitochondrial fraction of {{qc.mito.median}} (p99 = '
     '{{qc.mito.p99}}, max = {{qc.mito.max}}).\n'
     '\n'
     'The tutorial-standard thresholds (min_genes={{filter.min_genes}}, '
     'max_pct_mt={{filter.max_pct_mt}}, min_cells={{filter.min_cells}}) fit this data '
     'well: the min_genes floor removed no cells (the dataset was already reasonably '
     'clean), while the mitochondrial cap removed {{filter.cells_removed}} cells '
     '({{filter.pct_cells_removed}}) that sat in the long high-mito tail beyond the '
     'p95–p99 range — consistent with stressed or dying cells rather than a genuine cell '
     'population. Gene filtering (present in at least {{filter.min_cells}} cells) dropped '
     '{{filter.genes_removed}} largely undetected genes, leaving {{filter.genes_after}} '
     'informative genes and {{filter.cells_after}} cells.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'Scrublet was run on the whole dataset as a single 10x run (no batch splitting '
     'needed). The doublet-score distribution was unimodal/right-skewed '
     '({{doublet.distribution}}), so per the standard rule I used median + 3×MAD '
     '({{doublet.median_3mad}}) rather than a bimodal valley (none was present) and rather '
     'than the lighter Scrublet-automatic threshold ({{doublet.scrublet_auto}}), which is '
     'reserved for datasets where doublets were already removed upstream (e.g. by '
     'demultiplexing) — not the case here since this is raw, unprocessed 10x output from '
     'one donor. At threshold {{doublet.threshold}}, {{doublet.cells_removed}} cells '
     '({{doublet.pct_removed}}) were flagged and removed, in line with expected 10x '
     'doublet rates for this cell loading.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'With a single donor and no batch key, PCA on the {{norm.n_hvgs}} highly variable '
     'genes (after total-count normalization to {{norm.target_sum}} and log1p) is the '
     "simpler and correct choice over scVI's batch-correcting latent space — there is no "
     'batch effect to correct here. The leading PCs captured a substantial share of '
     'variance ({{pca.pc1_variance}}, {{pca.pc2_variance}}, {{pca.pc3_variance}} for '
     'PC1–PC3 respectively), and {{pca.n_comps}} PCs were retained for downstream '
     'neighbor-graph construction.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering at resolution {{cluster.resolution}} produced {{cluster.n}} '
     'clusters ranging from {{cluster.smallest}} to {{cluster.largest}} cells, reflecting '
     'the expected diversity of major PBMC lineages plus a small platelet population.\n'
     '\n'
     '{{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, fine-grained model) assigned each cluster a majority '
     'label via majority voting, cross-checked against the coarser Immune_All_High.pkl '
     'model as a second opinion; both models agreed on cell identity at their respective '
     'resolutions, with no contradictions requiring relabeling.\n'
     '\n'
     'Canonical marker genes, checked with check_markers against the full per-cluster '
     'differential expression ranking, confirmed every assigned label:\n'
     '- **Tem/Trm cytotoxic T cells** (cluster 0): high CD3D co-expressed with CD8A, CD8B, '
     'GZMK and CCL5 (all strongly enriched with high rank and low padj, and largely absent '
     'from other clusters), consistent with cytotoxic/effector CD8 T cells.\n'
     '- **B cells** (cluster 1): near-exclusive CD79A and MS4A1 expression (both '
     'top-ranked, highly significant, and depleted elsewhere), with the T-cell marker CD3D '
     'strongly depleted in this cluster.\n'
     '- **Tcm/Naive helper T cells** (clusters 2 and 3): CD3D and IL7R both significantly '
     'enriched, while CD8A/CD8B are not enriched, consistent with CD4 T cells; cluster 3 '
     'is dominated by ribosomal-protein transcripts typical of a quiescent naive subset, '
     'which CellTypist grouped with cluster 2 under the same fine label.\n'
     '- **Classical monocytes** (cluster 4): strong, near-universal enrichment of CD14, '
     'LYZ and S100A8, essentially absent in other clusters.\n'
     '- **DC** (cluster 5): high FCER1A and CD1C enrichment, without the CD14 enrichment '
     'seen in classical monocytes, consistent with conventional dendritic cells.\n'
     '- **CD16+ NK cells** (cluster 6): GNLY, NKG7 and FCGR3A all strongly and '
     'significantly enriched, with CD3D essentially absent, ruling out a T/NKT identity.\n'
     '- **Non-classical monocytes** (cluster 7): FCGR3A and MS4A7 strongly enriched '
     'without the CD14/S100A8 signature that marks classical monocytes.\n'
     '- **Megakaryocytes/platelets** (cluster 8): near-universal PF4 and PPBP expression '
     '(top-ranked in the cluster, absent elsewhere), a small '
     '(n={{celltype:Megakaryocytes/platelets}}) but clearly distinct population.\n'
     '\n'
     'Final cell-type composition:\n'
     '\n'
     '{{table:composition}}\n'
     '\n'
     'The most abundant populations were {{celltype:Tcm/Naive helper T cells}} CD4 T '
     'cells, {{celltype:Classical monocytes}} classical monocytes, and {{celltype:B '
     'cells}} B cells, alongside smaller {{celltype:Tem/Trm cytotoxic T cells}} cytotoxic '
     'T cell, {{celltype:Non-classical monocytes}} non-classical monocyte, '
     '{{celltype:CD16+ NK cells}} NK cell, {{celltype:DC}} dendritic cell, and '
     '{{celltype:Megakaryocytes/platelets}} platelet populations — a composition '
     'consistent with expectations for healthy human PBMCs.\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- This is a single donor/single run dataset with no biological replicates or '
     'condition to compare, so no composition or differential-expression testing across '
     'conditions was performed.\n'
     '- Clusters 2 and 3 both received the same fine-grained "Tcm/Naive helper T cells" '
     "label; cluster 3's ribosomal-high signature suggests it may represent a "
     'naive/quiescent subset rather than a truly distinct cell type, and finer '
     'sub-clustering could resolve this if desired.\n'
     '- Doublet removal used a MAD-based statistical threshold rather than a clear bimodal '
     'valley, since the score distribution was unimodal; some borderline low-scoring '
     'doublets or high-scoring genuine (e.g. large, transcript-rich) cells may be '
     'imperfectly classified.\n'
     '- The mitochondrial % cutoff ({{filter.max_pct_mt}}) is a widely used PBMC heuristic '
     'but is somewhat arbitrary; a small number of legitimate high-mito cells (e.g. '
     'metabolically active cells) may have been excluded.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'After QC filtering and doublet removal ({{final.pct_removed}} of original cells '
     'removed in total across mitochondrial and doublet filters), normalization, and '
     'PCA-based clustering, {{final.n_cells}} high-quality cells resolved into '
     '{{cluster.n}} clusters mapping to {{celltype.n}} annotated cell types, all confirmed '
     'by canonical marker expression. The recovered populations — CD4 and CD8 T cells, B '
     'cells, classical and non-classical monocytes, NK cells, dendritic cells, and '
     'platelets — represent the expected major lineages of peripheral blood, indicating a '
     'technically sound and biologically sensible annotation of this dataset.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
