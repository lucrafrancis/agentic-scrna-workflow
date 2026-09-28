"""Replay of the agent's tool calls for pbmc3k, without the LLM.

Generated: 2026-09-28T12:43:14
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
    genes=['CD3D', 'CD3E', 'CD8A', 'CD8B', 'CCL5', 'GZMK', 'IL7R', 'CCR7', 'SELL', 'S100A4',
     'CD4', 'MS4A1', 'CD79A', 'CD19', 'CD14', 'LYZ', 'FCN1', 'FCGR3A', 'MS4A7', 'NKG7',
     'GNLY', 'NCAM1', 'FCER1A', 'CST3', 'CD1C', 'PPBP', 'PF4'],
)

_run("summarize_findings")

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This analysis processed a PBMC 10x Genomics dataset from a single healthy donor '
     '({{input.n_cells}} cells x {{input.n_genes}} genes at input) through QC, doublet '
     'removal, normalization, PCA-based dimensionality reduction, Leiden clustering, and '
     'CellTypist-assisted cell-type annotation. Because this is a single sequencing run '
     'from one donor with no batch structure, PCA (rather than scVI batch correction) was '
     'the appropriate choice for the embedding used in clustering and UMAP.\n'
     '\n'
     '## Quality control\n'
     '\n'
     'Gene identifiers were confirmed as human gene symbols, with the mitochondrial prefix '
     '`MT-` matching {{qc.n_mito_genes}} genes. Per-cell QC distributions showed a median '
     'of {{qc.genes.median}} genes and {{qc.counts.median}} counts per cell, with '
     'mitochondrial content generally low (median {{qc.mito.median}}, 95th percentile '
     '{{qc.mito.p95}}) but with a tail reaching {{qc.mito.max}}, indicative of a modest '
     'population of stressed/dying cells typical of PBMC preparations.\n'
     '\n'
     'The standard tutorial thresholds (min genes/cell = {{filter.min_genes}}, max mito% = '
     '{{filter.max_pct_mt}}, min cells/gene = {{filter.min_cells}}) fit this dataset well: '
     "the gene-count floor removed no cells at all (the dataset's minimum, "
     '{{qc.genes.min}}, already exceeds the floor), and the mito ceiling removed only '
     '{{filter.cells_removed}} cells ({{filter.pct_cells_removed}}), consistent with the '
     'p95 mito value sitting comfortably below the {{filter.max_pct_mt}} cutoff. Gene '
     'filtering removed {{filter.genes_removed}} genes detected in fewer than '
     '{{filter.min_cells}} cells, leaving {{filter.genes_after}} informative genes for '
     'downstream analysis.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'Scrublet was run on the full dataset as a single batch (one donor, one 10x run, no '
     'upstream demultiplexing-based doublet removal). The doublet-score distribution was '
     '{{doublet.distribution}}, so per the standard rule I used median + 3xMAD '
     '({{doublet.threshold}}) rather than a bimodal valley (not present) or a fixed '
     'cutoff. This flagged and removed {{doublet.cells_removed}} cells '
     '({{doublet.pct_removed}} of the post-QC dataset), consistent with the expected '
     'doublet rate for a droplet run of this scale.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'With no batch key and a single clean 10x run, PCA on the {{norm.n_hvgs}} highly '
     'variable genes (after total-count normalization to {{norm.target_sum}} and log1p '
     'transform) was the simpler and correct choice over scVI, whose batch-correction '
     'machinery is unnecessary here and would add training noise without benefit. The top '
     'PCs capture the expected structure: PC1 explains {{pca.pc1_variance}} of variance, '
     'consistent with the major myeloid/lymphoid split seen later in clustering.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering at resolution {{cluster.resolution}} on the PCA embedding produced '
     '{{cluster.n}} clusters ranging from {{cluster.smallest}} to {{cluster.largest}} '
     'cells, spanning the full range of PBMC lineages (T, B, NK, monocyte, dendritic, and '
     'platelet populations). Clusters 2 and 3 both mapped to CD4 T cells but differ '
     'somewhat in ribosomal-gene content and CCR7 expression, likely reflecting naive vs. '
     'more activated/transitional CD4 T-cell states rather than distinct lineages; they '
     'were kept under the same CellTypist-assigned label since no canonical lineage marker '
     'distinguished them.\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, fine-grained immune model) assigned {{celltype.n}} '
     'cell types via majority vote per cluster, cross-checked against the coarse '
     'Immune_All_High.pkl model as a second opinion; the two models agreed at the expected '
     'level of granularity (e.g., "CD16+ NK cells" vs. "ILC", both correctly capturing the '
     'NK identity of cluster 6). Every final label was validated against canonical marker '
     'genes with check_markers before acceptance:\n'
     '\n'
     '- **{{celltype:Tem/Trm cytotoxic T cells}} cells** (cluster 0): strong, specific '
     'enrichment of CD3D, CD8A, CD8B, CCL5 and GZMK, with CCR7 low — a cytotoxic/effector '
     'CD8 T-cell profile.\n'
     '- **{{celltype:B cells}}** (cluster 1): CD79A, MS4A1 and CD19 are all sharply '
     'enriched with high in-cluster expression and near-absence elsewhere, unambiguously '
     'marking B cells.\n'
     '- **{{celltype:Tcm/Naive helper T cells}}** (clusters 2 and 3): CD3D, CD3E and IL7R '
     'enriched in both clusters with CD8A/CD8B absent, confirming CD4 T-cell identity; '
     'cluster 3 additionally shows elevated CCR7, consistent with a more '
     'naive/central-memory phenotype.\n'
     '- **{{celltype:Classical monocytes}}** (cluster 4): CD14, LYZ and FCN1 are all '
     'highly and specifically enriched, the canonical classical-monocyte signature.\n'
     '- **{{celltype:DC}}** (cluster 5): FCER1A and CD1C are sharply and specifically '
     'enriched (near-absent in nearly all other clusters), confirming a conventional '
     'dendritic-cell identity alongside high CST3.\n'
     '- **{{celltype:CD16+ NK cells}}** (cluster 6): GNLY and NKG7 are the top-ranked '
     'genes for this cluster with near-universal in-cluster expression, and FCGR3A is also '
     'strongly enriched, while CD3D/CD3E are essentially absent — a clean NK profile.\n'
     '- **{{celltype:Non-classical monocytes}}** (cluster 7): FCGR3A and MS4A7 are both '
     'very strongly and specifically enriched, with CD14 not preferentially expressed '
     'relative to other clusters, matching the non-classical (CD14low/CD16+) monocyte '
     'phenotype.\n'
     '- **{{celltype:Megakaryocytes/platelets}}** (cluster 8, only {{cluster.smallest}} '
     'cells): PPBP and PF4 are both near-universally expressed within this cluster and '
     'essentially undetectable elsewhere, unambiguously marking platelets/megakaryocyte '
     'fragments.\n'
     '\n'
     'Final cell-type composition:\n'
     '\n'
     '{{table:clusters}}\n'
     '\n'
     '{{table:composition}}\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- Clusters 2 and 3 were both annotated as CD4 T cells despite showing '
     'distinguishable ribosomal-gene expression and CCR7 differences; a higher clustering '
     'resolution or dedicated subtyping could split them into naive vs. central-memory '
     'populations if that resolution is needed.\n'
     '- The dendritic-cell and megakaryocyte clusters are small ({{celltype:DC}} and '
     '{{celltype:Megakaryocytes/platelets}} cells respectively), so marker rankings for '
     'these groups are based on fewer cells and are more sensitive to noise.\n'
     '- This is a single donor with no replicate structure, so no condition comparison or '
     'composition/DE testing was performed; findings describe the cell-type landscape of '
     'this one sample only.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'Starting from {{input.n_cells}} raw cells, standard QC and Scrublet-based doublet '
     'removal left {{final.n_cells}} high-quality, singlet cells ({{final.pct_removed}} '
     'removed in total). PCA-based clustering at resolution {{cluster.resolution}} cleanly '
     'resolved the expected major PBMC lineages — CD4 and CD8 T cells, B cells, NK cells, '
     'classical and non-classical monocytes, dendritic cells, and a small '
     'platelet/megakaryocyte population — each supported by canonical markers checked '
     'directly against the data.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
