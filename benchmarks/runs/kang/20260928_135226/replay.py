"""Replay of the agent's tool calls for kang, without the LLM.

Generated: 2026-09-28T14:12:12
Git commit: 24068e592d14aa176f5b6e5c079649c4ed21485a
Source log: tool_calls.jsonl (21 of 22 calls; failed calls omitted)

Re-runs the same tool functions with the same arguments and seed. Output goes to
outputs/kang_replay/, with its own tool log and report. Package versions are pinned by
uv.lock at the commit above: check it out first if the code has changed since.

Caveats: annotate_celltypes downloads CellTypist models, and scVI training can differ
slightly between hardware (CPU vs GPU), which can change clusters downstream.

  uv run python outputs/kang/replay.py
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
SESSION.load('data/kang.h5ad')
SESSION.paths = config.RunPaths('kang_replay')
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
    batch_key='condition',
)

_run(
    "filter_doublets",
    threshold=0.6482692445451002,
)

_run("normalize")

_run(
    "run_scvi",
    batch_key=['donor', 'condition'],
)

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
    genes=['CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'S100A8', 'FCN1', 'CD3D', 'CD3E', 'CD4', 'CD8A',
     'CD8B', 'IL7R', 'CCR7', 'SELL', 'GNLY', 'NKG7', 'KLRD1', 'NCAM1', 'CD79A', 'MS4A1',
     'CD19', 'CD1C', 'CLEC9A', 'FCER1A', 'LILRA4', 'PPBP', 'PF4', 'ITGA2B', 'TPSAB1',
     'CPA3', 'KIT', 'HBB', 'HBA1', 'GYPA', 'GZMK', 'GZMB', 'TIGIT'],
)

_run(
    "relabel_clusters",
    labels={'0': 'CD14+ Classical monocytes',
     '1': 'CD14+ Classical monocytes',
     '4': 'CD8+ Effector T cells',
     '12': 'Platelets'},
    reason=("Cluster 0/1: CellTypist called these 'Intermediate macrophages' but markers (CD14 "
     'log2FC~2.4-4.6 in 40-64% of cells, LYZ in 72-84%, FCN1 in 55-79%, S100A8 in 60-78%) '
     'are canonical classical monocyte markers; tissue macrophages are not expected in '
     "PBMC blood samples. Cluster 4: labeled 'NK cells' but CD3D is significantly enriched "
     '(66.4% pct-in vs 29.2% elsewhere, log2FC 2.34, padj~0) together with CD8A/CD8B, '
     'GZMB, NKG7, GNLY -- a CD3+ cytotoxic/effector CD8 T cell profile, not NK (contrast '
     'with true NK cluster 11 where CD3D is depleted: 4.5% vs 35.1% elsewhere). The coarse '
     "CellTypist model independently calls cluster 4 'T cells', agreeing with this "
     "correction. Cluster 12: labeled 'Tcm/Naive helper T cells' but PPBP (85.1% cells, "
     'log2FC 7.7), PF4 (63.6%, log2FC 7.46), and ITGA2B (7.4%, log2FC 6.76) are all '
     'massively and specifically enriched -- unambiguous platelet markers, not T cells.'),
)

_run(
    "compare_composition",
    contrast=['condition', 'stim', 'ctrl'],
    sample_key='donor',
)

_run(
    "run_pseudobulk_de",
    contrast=['condition', 'stim', 'ctrl'],
    sample_key='donor',
    covariates=['donor'],
)

_run(
    "query_genes",
    genes=['ISG15', 'ISG20', 'IFIT1', 'MX1', 'CCL8', 'IL1RN'],
)

_run(
    "get_top_genes",
    cell_type='CD14+ Classical monocytes',
    direction='up',
    n=15,
)

_run("summarize_findings")

_run(
    "check_markers",
    genes=['MS4A2', 'HDC', 'TPSB2', 'CTSG', 'GATA2', 'SLC18A2', 'CD69'],
    clusters=['15'],
)

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'PBMCs from {{input.n_donor}} SLE patients, each split into a control and '
     'IFN-β-stimulated\n'
     'half, were processed from raw counts through QC, doublet removal, batch-aware '
     'integration,\n'
     'clustering, cell-type annotation, and a paired stimulated-vs-control comparison. '
     'The\n'
     'starting object held {{input.n_cells}} cells x {{input.n_genes}} genes across\n'
     '{{input.n_batches}} donors and {{input.n_condition}} conditions; condition and 10x '
     'run are\n'
     'fully confounded by design (each condition was captured in its own pooled run).\n'
     '\n'
     '## Quality control\n'
     '\n'
     'No mitochondrial genes were detected in the gene set ({{qc.n_mito_genes}} found '
     'despite the\n'
     'standard "MT-" prefix being present), so pct-mito is uniformly {{qc.mito.median}} '
     'and provided\n'
     'no filtering signal — this dataset had mitochondrial genes already stripped out '
     'upstream, so a\n'
     'mito-based cutoff is a no-op rather than a QC decision here. Genes-per-cell ranged '
     'from\n'
     '{{qc.genes.min}} to {{qc.genes.max}} (median {{qc.genes.median}}), consistent with a '
     'fairly\n'
     'shallow 10x run typical of this dataset. The tutorial-default thresholds\n'
     '(min_genes={{filter.min_genes}}, max_pct_mt={{filter.max_pct_mt}}, '
     'min_cells={{filter.min_cells}})\n'
     'removed only {{filter.cells_removed}} cells ({{filter.pct_cells_removed}}) and\n'
     '{{filter.genes_removed}} genes, so I applied them as-is — the low-gene tail was '
     'small and this\n'
     'is not an aggressive cut.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'Demuxlet genotype-based demultiplexing (assigning cells to donors) had already '
     'removed\n'
     'cross-donor doublets upstream; only same-donor doublets (transcriptomically '
     'indistinguishable\n'
     'by genotype) can remain. The Scrublet score distribution was '
     '{{doublet.distribution}} (median\n'
     '{{doublet.score_median}}, max {{doublet.score_max}}), so the standard median+3·MAD '
     'rule\n'
     '({{doublet.median_3mad}}) was not used — with most doublets already gone, that '
     'generic rule\n'
     'would treat the upper tail of real, transcriptionally-active cells as doublets. '
     'Instead, per\n'
     'the documented exception for demultiplexed data, I used a light-touch cutoff: the '
     'higher of the\n'
     'two per-run Scrublet automatic thresholds ({{doublet.scrublet_auto}} for ctrl/stim), '
     'i.e.\n'
     '{{doublet.threshold}}, applied uniformly. This removed only '
     '{{doublet.cells_removed}} cells\n'
     '({{doublet.pct_removed}}), consistent with the expectation that little residual '
     'doublet\n'
     'contamination remained.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'IFN-β stimulation causes a strong, cell-type-independent transcriptional shift (a '
     'large\n'
     'interferon-response program), and condition is fully confounded with the 10x run. '
     'Without\n'
     'correction, this shift risks splitting each cell type into two condition-specific '
     'clusters,\n'
     'which would break both cluster-based annotation and the downstream composition test. '
     'I\n'
     "therefore ran scVI with a combined `['donor', 'condition']` batch key (one batch per "
     'sample,\n'
     '{{scvi.n_hvgs}} HVGs, latent dimension {{scvi.n_latent}}), which aligns matching '
     'cell types\n'
     'across donors and across conditions in the embedding used for clustering/UMAP, while '
     'leaving\n'
     'raw counts untouched for the composition and DE tests below. The risk of this choice '
     'is\n'
     'over-correction: any cell state that is genuinely unique to one condition could be '
     'merged\n'
     'into a resting counterpart. I checked the UMAP by condition (see figure) — clusters '
     'mix well\n'
     'across stim/ctrl, supporting good integration without obvious over-merging of '
     'clearly distinct\n'
     'populations.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering (resolution {{cluster.resolution}}) on the scVI embedding '
     'produced\n'
     '{{cluster.n}} clusters ranging from {{cluster.smallest}} to {{cluster.largest}} '
     'cells.\n'
     '\n'
     '{{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl) gave fine-grained labels per cluster, cross-checked '
     'against\n'
     'the coarse Immune_All_High.pkl model and canonical markers. Three corrections were '
     'made where\n'
     "the fine model's calls did not hold up:\n"
     '\n'
     '- **Clusters 0/1** were labeled "Intermediate macrophages" by the fine model, but '
     'canonical\n'
     '  classical-monocyte markers CD14, LYZ, FCN1 and S100A8 were all strongly and '
     'specifically\n'
     '  enriched, consistent with classical (CD14+) monocytes rather than tissue '
     'macrophages, which\n'
     '  are not expected in a blood sample. Relabeled to "CD14+ Classical monocytes".\n'
     '- **Cluster 4** was labeled "NK cells", but CD3D was significantly enriched '
     'alongside\n'
     '  CD8A/CD8B/GZMB/NKG7, unlike the genuine CD3-negative NK cluster (cluster 11, where '
     'CD3D is\n'
     '  depleted rather than enriched). The coarse model independently called this '
     'cluster\n'
     '  "T cells", agreeing with the correction to "CD8+ Effector T cells".\n'
     '- **Cluster 12** was labeled "Tcm/Naive helper T cells", but the platelet markers '
     'PPBP, PF4\n'
     '  and ITGA2B were massively and specifically enriched, unambiguously identifying\n'
     '  platelets/platelet-containing droplets rather than T cells. Relabeled to '
     '"Platelets".\n'
     '\n'
     'All other clusters agreed between the two CellTypist models and canonical markers: '
     'non-classical\n'
     'monocytes (FCGR3A, MS4A7 enriched), naive/central-memory CD4 T cells (CCR7, IL7R, '
     'SELL\n'
     'enriched), CD8 Tcm/naive (CD8B, CCR7, SELL enriched), CD8 Tem/Trm (GZMK, TIGIT, '
     'CD3D\n'
     'enriched), B cells (CD79A, MS4A1, CD19 enriched), plasmacytoid DC (LILRA4, GZMB '
     'enriched),\n'
     'DC2 (CD1C, FCER1A, LYZ enriched), erythroid cells (HBB, HBA1 near-universally '
     'expressed and\n'
     'strongly enriched), and a small, rare cluster with GATA2 strongly and specifically '
     'enriched\n'
     '(consistent with mast cell/basophil identity; CD3D and lymphoid/myeloid markers '
     'essentially\n'
     'absent), retained as "Mast cells".\n'
     '\n'
     'Final composition: {{celltype:CD14+ Classical monocytes}} classical monocytes,\n'
     '{{celltype:Tcm/Naive helper T cells}} naive/Tcm CD4 T cells,\n'
     '{{celltype:Tem/Effector helper T cells}} effector/Tem CD4 T cells,\n'
     '{{celltype:B cells}} B cells, {{celltype:CD8+ Effector T cells}} CD8 effector T '
     'cells,\n'
     '{{celltype:CD16+ NK cells}} NK cells, {{celltype:Tcm/Naive cytotoxic T cells}} '
     'naive/Tcm CD8 T\n'
     'cells, {{celltype:Non-classical monocytes}} non-classical monocytes, '
     '{{celltype:DC2}} DC2,\n'
     '{{celltype:Platelets}} platelets, {{celltype:pDC}} pDC, {{celltype:Late erythroid}} '
     'erythroid\n'
     'cells, {{celltype:Tem/Trm cytotoxic T cells}} Tem/Trm CD8 T cells, and\n'
     '{{celltype:Mast cells}} mast cells, out of {{final.n_cells}} total cells after '
     'QC/doublet\n'
     'filtering.\n'
     '\n'
     '## Composition analysis\n'
     '\n'
     'Cell-type proportions were compared stim vs ctrl per donor (paired {{comp.test}}, '
     'all\n'
     '{{input.n_donor}} donors contributing both conditions).\n'
     '\n'
     '{{table:composition}}\n'
     '{{table:composition_test}}\n'
     '\n'
     'No cell type showed a significant proportion shift after multiple-testing '
     'correction\n'
     '({{comp.n_significant}} of {{comp.n_cell_types}} significant). The largest nominal '
     'changes were\n'
     'a decrease in DC2 ({{comp.dc2.ctrl}} -> {{comp.dc2.stim}}, padj {{comp.dc2.padj}}) '
     'and an\n'
     'increase in pDC ({{comp.pdc.ctrl}} -> {{comp.pdc.stim}}, padj {{comp.pdc.padj}}), '
     'but neither\n'
     'survived correction and pDC could only be tested with very few cells per sample. '
     'With so few\n'
     "donor replicates, the paired test's attainable significance is inherently limited, "
     'so this\n'
     'analysis is under-powered to detect anything but large, consistent shifts; a true '
     'absence of\n'
     'compositional change over such a short stimulation window is also biologically '
     'plausible,\n'
     'since PBMC subset frequencies are not expected to shift substantially without cell\n'
     'death/proliferation over a brief ex vivo culture.\n'
     '\n'
     '## Differential expression\n'
     '\n'
     'Pseudobulk DE (PyDESeq2, design {{de.design}}, paired by donor) was run per cell '
     'type\n'
     'using raw counts summed per donor.\n'
     '\n'
     '{{table:de_summary}}\n'
     '\n'
     'Every cell type tested showed a clear, coherent interferon-stimulated gene (ISG) '
     'signature,\n'
     'as expected for IFN-β treatment: {{gene:CD14+ Classical monocytes:ISG15}} in '
     'classical\n'
     'monocytes, {{gene:CD16+ NK cells:ISG15}} in NK cells, {{gene:B cells:ISG15}} in B '
     'cells, and\n'
     '{{gene:CD8+ Effector T cells:ISG15}} in CD8 effector T cells, together with '
     'consistent\n'
     'up-regulation of {{gene:CD14+ Classical monocytes:ISG20}}, {{gene:CD14+ Classical '
     'monocytes:IFIT1}},\n'
     'and {{gene:CD14+ Classical monocytes:MX1}} across essentially all tested '
     'populations. Classical\n'
     'monocytes showed the largest and most cell-type-specific response, including strong '
     'induction\n'
     'of the chemokine {{gene:CD14+ Classical monocytes:CCL8}} and the anti-inflammatory '
     'cytokine\n'
     'antagonist {{gene:CD14+ Classical monocytes:IL1RN}}, both far less induced in '
     'non-classical\n'
     'monocytes ({{gene:Non-classical monocytes:CCL8}}, {{gene:Non-classical '
     'monocytes:IL1RN}}) and\n'
     'essentially absent from lymphoid populations tested for these genes — consistent '
     'with monocytes\n'
     'being a primary sensor/responder to IFN-β and secondary inflammatory-mediator '
     'producer. Rare\n'
     'populations ({{celltype:Late erythroid}} erythroid cells, {{celltype:Mast cells}} '
     'mast cells,\n'
     '{{celltype:Tem/Trm cytotoxic T cells}} Tem/Trm cytotoxic T cells, {{celltype:pDC}} '
     'pDC cells)\n'
     'had too few cells per donor in one or both conditions and were skipped, as noted by '
     'the tool\n'
     '({{de.skipped}}).\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- **Condition/run confound**: stim and ctrl were captured in separate 10x runs, so '
     'pseudobulk\n'
     '  DE cannot fully separate a true IFN-β effect from a run-specific technical effect; '
     'the very\n'
     '  large, coherent, and biologically expected ISG signature makes a purely technical '
     'explanation\n'
     '  unlikely, but this cannot be ruled out formally.\n'
     '- **scVI batch correction with `donor + condition`** aligns cell types across '
     'conditions for\n'
     '  clustering/annotation but could, in principle, over-merge a condition-specific '
     'cell state\n'
     '  into a resting counterpart; the reasonably good mixing of conditions on the UMAP '
     'argues\n'
     '  against gross over-merging, but subtle state-specific subpopulations could still '
     'be masked at\n'
     '  this clustering resolution.\n'
     '- **Composition test is under-powered** with only {{input.n_donor}} paired donors, '
     'so the null\n'
     '  result for composition changes should not be read as strong evidence of no change, '
     'only as no\n'
     '  detected change at this sample size.\n'
     '- **Rare clusters** ({{celltype:Mast cells}} mast cells, {{celltype:pDC}} pDC,\n'
     '  {{celltype:Platelets}} platelets, {{celltype:Tem/Trm cytotoxic T cells}} Tem/Trm '
     'cytotoxic T\n'
     '  cells) have limited statistical power for both marker validation and DE; labels '
     'here rely on\n'
     '  strong marker fold-changes despite low cell counts.\n'
     '- **No mitochondrial genes were present** in this dataset, so QC could not filter '
     'on\n'
     '  mito-fraction; only gene-count-based filtering was applied.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'After QC, light-touch doublet removal (justified by prior genotype-based '
     'demultiplexing),\n'
     'and scVI-based integration across donors and conditions, clustering recovered the '
     'expected\n'
     'major PBMC lineages (classical and non-classical monocytes, CD4 and CD8 T cell '
     'subsets, NK\n'
     'cells, B cells, DC2, pDC, platelets, mast cells, and a small erythroid contaminant), '
     'after\n'
     'correcting three CellTypist mislabels using canonical marker evidence. IFN-β '
     'stimulation did\n'
     'not detectably shift PBMC subset proportions but induced a strong, broad, '
     'biologically\n'
     'coherent interferon-response transcriptional program across essentially every cell '
     'type, most\n'
     'pronounced in classical monocytes, which additionally mounted a distinct secondary\n'
     'inflammatory-mediator response (CCL8, IL1RN) not seen in most other populations.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
