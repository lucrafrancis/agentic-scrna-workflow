"""Replay of the agent's tool calls for pbmc3k, without the LLM.

Generated: 2026-09-28T12:41:08
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

_run("annotate_celltypes")

_run(
    "check_markers",
    genes=['CD3D', 'CD3E', 'CD8A', 'CD8B', 'CD4', 'IL7R', 'CCR7', 'GZMK', 'CCL5', 'CD79A',
     'MS4A1', 'CD19', 'LYZ', 'CD14', 'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'CD1C', 'GNLY',
     'NKG7', 'GZMB', 'NCAM1', 'PPBP', 'PF4', 'ITGA2B'],
)

_run("summarize_findings")

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This report covers a standard single-sample PBMC analysis: one healthy donor, one '
     '10x Genomics run, {{input.n_cells}} cells and {{input.n_genes}} genes as input. No '
     'batch key was detected, consistent with the description of a single run, so no batch '
     'integration was needed at any stage. The workflow proceeded from raw counts through '
     'QC filtering, doublet removal, normalization, PCA-based dimensionality reduction, '
     'Leiden clustering, and CellTypist-assisted cell-type annotation, with every '
     'canonical marker checked before finalizing labels.\n'
     '\n'
     '## Quality control\n'
     '\n'
     'Gene identifiers were confirmed as human gene symbols, with the mitochondrial gene '
     'prefix `MT-` detected ({{qc.n_mito_genes}} mitochondrial genes found). Per-cell QC '
     'metrics were unremarkable for healthy PBMCs: median genes/cell {{qc.genes.median}}, '
     'median UMI counts {{qc.counts.median}}, median mitochondrial fraction '
     '{{qc.mito.median}} (up to {{qc.mito.max}} in a tail of outlier cells).\n'
     '\n'
     'I applied the tutorial-standard thresholds recommended for this dataset: minimum '
     '{{filter.min_genes}} genes/cell, maximum {{filter.max_pct_mt}} mitochondrial '
     'content, and genes present in at least {{filter.min_cells}} cells. These are well '
     'suited here because the gene-count distribution had essentially no low-count tail '
     '(min was already {{qc.genes.min}} genes), so the min_genes cutoff removed nothing; '
     'the mitochondrial cap did the actual work, removing {{filter.cells_removed}} likely '
     'dying/stressed cells ({{filter.pct_cells_removed}} of the input). Gene filtering '
     'dropped {{filter.genes_removed}} genes detected in too few cells, mostly '
     'uninformative low-expression transcripts. After this step: {{filter.cells_after}} '
     'cells and {{filter.genes_after}} genes remained.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'Scrublet was run on this single 10x run as one batch (no per-run splitting needed). '
     'The doublet-score distribution was unimodal with a long right tail rather than two '
     'clearly separated peaks ({{doublet.distribution}}), so per the decision rule I used '
     'median + 3×MAD ({{doublet.threshold}}) rather than a bimodal valley. There is no '
     'indication doublets were removed upstream (e.g. via genotype demultiplexing), so the '
     'standard rule — not the lighter-touch Scrublet-automatic threshold — was the '
     'appropriate choice. This flagged and removed {{doublet.cells_removed}} cells '
     '({{doublet.pct_removed}}), a plausible heterotypic doublet rate for this loading '
     'density. After doublet removal: {{doublet.cells_after}} cells remained for '
     'downstream analysis.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'Since this dataset is a single clean batch from one donor with no technical grouping '
     'variable, PCA on the normalized, log-transformed HVG matrix (top {{norm.n_hvgs}} '
     "genes, target sum {{norm.target_sum}}) is the correct and simplest choice — scVI's "
     'batch-correction machinery is unnecessary overhead here and could over-smooth '
     'genuine biological structure in a single sample. PCA was run with {{pca.n_comps}} '
     'components; the leading PCs (PC1 {{pca.pc1_variance}}, PC2 {{pca.pc2_variance}}, PC3 '
     '{{pca.pc3_variance}} of variance) show a clear elbow typical of well-structured '
     'immune cell data.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering (resolution {{cluster.resolution}}) on the PCA representation '
     'yielded {{cluster.n}} clusters, ranging in size from {{cluster.smallest}} to '
     '{{cluster.largest}} cells. Cluster sizes and UMAP topology are consistent with the '
     'expected composition of PBMCs: a few large lymphocyte and monocyte populations plus '
     'small, distinct populations (dendritic cells, platelets).\n'
     '\n'
     '{{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, fine-grained) was used for majority-vote annotation '
     'per cluster, cross-checked against the coarser Immune_All_High.pkl model as a second '
     'opinion. The two models agreed at the appropriate granularity for every cluster '
     '(e.g. cluster 6 was called "CD16+ NK cells" by the fine model and the broader "ILC" '
     'category by the coarse model — not a contradiction, since NK cells are ILCs).\n'
     '\n'
     'Before accepting labels, canonical markers were checked directly against the full '
     'marker ranking for every final cell type (not just ambiguous ones):\n'
     '\n'
     '- **Cytotoxic T cells** (cluster 0): CCL5 was the single top-ranked marker gene for '
     'this cluster, expressed in nearly all of its cells and strongly depleted elsewhere, '
     'together with high CD8A/CD8B and GZMK, and CD3D/CD3E positivity confirming a '
     'T-lineage identity co-expressed with cytotoxic granzymes.\n'
     '- **B cells** (cluster 1): CD79A and MS4A1 were both highly specific — expressed in '
     'the large majority of cluster cells and rare outside it — with CD3D strongly '
     'depleted — a clean B-cell signature.\n'
     '- **CD4 T cells / Tcm-Naive helper T cells** (clusters 2 and 3): both show strong '
     'CD3D/CD3E/IL7R positivity and CD8A/CD79A absence. Cluster 3 shows markedly higher '
     'CCR7 than cluster 2, suggesting a more naive-like subpopulation within the same '
     "broad CellTypist label; both were kept under the same annotation since CellTypist's "
     'own majority vote did not distinguish them and the marker signal is a difference of '
     'degree, not identity.\n'
     '- **Classical monocytes** (cluster 4): LYZ was expressed in essentially every cell '
     'of this cluster and CD14 in the majority, both strongly enriched versus every other '
     'cluster — a textbook classical monocyte signature.\n'
     '- **Dendritic cells** (cluster 5): FCER1A and CD1C were both highly specific to this '
     'cluster and nearly absent elsewhere, confirming a conventional DC identity distinct '
     'from monocytes despite shared CD74/HLA-DR expression.\n'
     '- **CD16+ NK cells** (cluster 6): GNLY, NKG7 and GZMB were the top-ranked genes for '
     'this cluster, expressed in nearly all its cells; FCGR3A was also strongly enriched, '
     'and CD3D was essentially absent — ruling out a T-cell identity and confirming NK '
     'cells.\n'
     '- **Non-classical monocytes** (cluster 7): FCGR3A and MS4A7 were both highly '
     'specific and strongly enriched, while CD14 was low and not significantly enriched — '
     'the expected classical-vs-non-classical monocyte split confirmed by reciprocal '
     'CD14/FCGR3A patterns between clusters 4 and 7.\n'
     '- **Megakaryocytes/platelets** (cluster 8): PPBP, PF4 and ITGA2B were essentially '
     'exclusive to this small cluster, confirming the platelet/megakaryocyte identity '
     'despite the very small cluster size ({{celltype:Megakaryocytes/platelets}} cells).\n'
     '\n'
     "No cluster's marker evidence contradicted its CellTypist label, so no relabeling was "
     'necessary.\n'
     '\n'
     '{{table:composition}}\n'
     '\n'
     'Final annotated dataset: {{final.n_cells}} cells across {{celltype.n}} cell types, '
     'following removal of {{final.cells_removed}} cells ({{final.pct_removed}} of the '
     'original input) through QC and doublet filtering combined.\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- This is a single donor/single run, so no batch-correction or condition comparison '
     "was performed or appropriate; all conclusions describe this one sample's composition "
     'and cannot be generalized without replicates.\n'
     '- Clusters 2 and 3 (both CD4 T cells) likely represent a naive/memory continuum '
     'rather than two distinct cell types; the marker difference (CCR7) is one of degree '
     'and would benefit from a targeted resolution increase or additional markers (e.g. '
     'CCR7/SELL/CD27 panels) if finer T-cell substates are needed.\n'
     '- The dendritic cell and megakaryocyte/platelet clusters are small ({{celltype:DC}} '
     'and {{celltype:Megakaryocytes/platelets}} cells respectively), so their marker '
     'statistics and any downstream inference about them carry wider uncertainty than the '
     'larger populations.\n'
     '- Doublet removal used an unsupervised score threshold (median + 3×MAD); while '
     'principled, some borderline transcriptionally-mixed cells (e.g. platelet-adhered '
     'lymphocytes, a known PBMC artifact) may still remain or have been removed '
     'unnecessarily.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'Standard QC and doublet-filtering thresholds were appropriate for this clean, '
     'single-sample PBMC dataset, removing a modest fraction of low-quality and doublet '
     'cells ({{final.pct_removed}} total) without needing batch correction. PCA-based '
     'clustering at resolution {{cluster.resolution}} resolved the expected major PBMC '
     'lineages — cytotoxic and helper T cells, B cells, classical and non-classical '
     'monocytes, dendritic cells, NK cells, and a small platelet/megakaryocyte population '
     '— each supported by canonical markers confirmed directly against the data, giving '
     'confidence in the final annotation.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
