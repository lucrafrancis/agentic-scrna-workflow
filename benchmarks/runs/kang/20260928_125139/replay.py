"""Replay of the agent's tool calls for kang, without the LLM.

Generated: 2026-09-28T13:12:44
Git commit: 24068e592d14aa176f5b6e5c079649c4ed21485a
Source log: tool_calls.jsonl (22 of 24 calls; failed calls omitted)

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

_run("identify_markers")

_run("annotate_celltypes")

_run(
    "check_markers",
    genes=['CD3D', 'CD3E', 'CD8A', 'CD8B', 'CD4', 'IL7R', 'CCR7', 'SELL', 'GNLY', 'NKG7', 'KLRD1',
     'NCAM1', 'FCGR3A', 'CD14', 'LYZ', 'FCN1', 'MS4A7', 'CD1C', 'FCER1A', 'CLEC10A',
     'CLEC4C', 'IL3RA', 'TCF4', 'MS4A1', 'CD79A', 'CD19', 'PPBP', 'PF4', 'ITGA2B', 'TPSAB1',
     'CPA3', 'KIT', 'HBB', 'GZMK', 'GZMA', 'GZMB'],
)

_run(
    "relabel_clusters",
    labels={'0': 'CD14+ Monocytes', '1': 'CD14+ Monocytes', '4': 'CD8+ T cells', '12': 'Platelets'},
    reason=("Cluster 0/1: CellTypist called these 'Intermediate macrophages', but PBMC blood "
     "shouldn't contain tissue macrophages; markers show canonical classical monocyte "
     'identity (CD14 log2FC 2.35/4.58, pct 39.6%/64.0%; LYZ pct 84.1%/72.3%; FCN1 pct '
     "78.9%/54.3%), so relabeled CD14+ Monocytes. Cluster 4: labeled 'NK cells' but CD3D "
     '(pct 66.4% vs 29.2% elsewhere, log2FC 2.34, padj~0) and CD8A/CD8B (pct 40.1%/30.7%, '
     'log2FC 4.46/3.74) are clearly enriched, which is incompatible with NK identity (real '
     'NK cluster 11 has CD3D pct only 4.5%, strongly depleted); the GNLY/NKG7/GZMB '
     'signature reflects a cytotoxic CD8+ effector T-cell state, so relabeled CD8+ T '
     'cells. Cluster 12: labeled a helper-T subset but dominated by platelet markers PPBP '
     '(pct 85.1%, log2FC 7.7), PF4 (pct 63.6%, log2FC 7.46), and ITGA2B (log2FC 6.76), so '
     'relabeled Platelets.'),
)

_run("identify_markers")

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
    "get_top_genes",
    cell_type='CD14+ Monocytes',
    n=15,
)

_run(
    "query_genes",
    genes=['ISG15', 'ISG20', 'IFI6', 'MX1', 'IFIT1', 'CXCL10', 'CCL8', 'IL1RN', 'STAT1'],
)

_run("summarize_findings")

_run(
    "check_markers",
    genes=['MS4A2', 'HDC', 'GATA2', 'CTSG', 'VWA5A', 'CMA1', 'ENPP3', 'PRSS33', 'SLC18A2'],
    clusters=['15'],
)

_run(
    "generate_report",
    report_markdown=('\n'
     '## Overview\n'
     '\n'
     'This analysis processes PBMCs from {{input.n_donor}} SLE patients (Kang et al., Nat '
     'Biotechnol), each split into a paired IFN-β-stimulated and unstimulated (control) '
     'sample, to identify cell types and test how stimulation shifts cell-type composition '
     'and gene expression. Starting from {{input.n_cells}} cells x {{input.n_genes}} '
     'genes, the pipeline went through QC filtering, doublet removal, normalization, '
     'batch-corrected embedding, clustering, CellTypist annotation (with manual correction '
     'where markers disagreed), and paired condition comparisons, ending with '
     '{{final.n_cells}} cells across {{celltype.n}} annotated cell types.\n'
     '\n'
     'A key structural fact drives several downstream choices: the two conditions were '
     'captured in **two separate 10x runs** (stim pooled, ctrl pooled), so `condition` is '
     'fully confounded with the technical (10x run) batch. Donors were later assigned by '
     'genotype demultiplexing, which also removed cross-donor doublets.\n'
     '\n'
     '## Quality control\n'
     '\n'
     'No mitochondrial genes were detected in this dataset ({{qc.n_mito_genes}} found; '
     '{{qc.mito.median}} pct_mt at every percentile), meaning it had already been stripped '
     'of MT genes upstream — the mito filter is effectively a no-op here rather than a '
     'real QC signal, so I relied on the gene-count floor instead. I used the standard '
     'tutorial-default thresholds (min_genes={{filter.min_genes}}, '
     'max_pct_mt={{filter.max_pct_mt}}, min_cells={{filter.min_cells}}), which is '
     'reasonable given the counts distribution (median genes/cell = {{qc.genes.median}}, '
     'median total counts = {{qc.counts.median}}) shows no bimodal low-quality tail beyond '
     'the standard floor. This removed only {{filter.cells_removed}} cells '
     '({{filter.pct_cells_removed}}) and {{filter.genes_removed}} genes.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'The dataset description states doublets were largely already removed by '
     'genotype-based demultiplexing (demuxlet), except doublets formed from two cells of '
     "the *same* donor, which don't produce a mixed-genotype signal. This is the "
     'documented exception case: a strict median+3×MAD cutoff ({{doublet.median_3mad}}) '
     "would over-remove genuine cells from a population that's already mostly clean. "
     'Scrublet was run separately per 10x run (`condition`, since donors were pooled '
     'within each run rather than being separate runs), giving per-run automatic '
     'thresholds of {{doublet.scrublet_auto}} (ctrl / stim). I used the higher of the two '
     '({{doublet.threshold}}, from the stim run) as a light-touch cutoff, which removed '
     'only {{doublet.cells_removed}} cells ({{doublet.pct_removed}}) — consistent with the '
     'expectation that few residual doublets remain.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'Because `condition` is confounded with the 10x run, and IFN-β stimulation is known '
     'to cause a strong, cell-type-wide transcriptional shift, clustering directly on '
     'uncorrected PCA risks splitting each cell type into stim/ctrl-specific clusters, '
     'which would make annotation inconsistent between conditions and any composition test '
     "meaningless. I used **scVI** with batch_key = `['donor', 'condition']` (one batch "
     'per sample), correcting for both donor-to-donor variation and the condition/run '
     'effect in the embedding used for clustering and UMAP, while leaving raw counts '
     'untouched for the downstream DE and composition analyses (which is exactly where the '
     'biological stim-vs-ctrl signal needs to be preserved). This is a deliberate '
     'trade-off: over-correction could in principle merge a stimulation-specific state '
     'into another type, so I checked whether clusters mixed by condition (see Clustering '
     'below) and rely on the pseudobulk DE — run on raw counts — to recover the real '
     'condition effect regardless of how the embedding was corrected.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering (resolution {{cluster.resolution}}) on the scVI embedding produced '
     '{{cluster.n}} clusters ranging from {{cluster.smallest}} to {{cluster.largest}} '
     'cells. Clusters correspond well to expected PBMC lineages (monocyte, T-cell, NK, '
     'B-cell, DC, platelet, erythroid and mast-cell clusters), and — reassuringly for the '
     'batch-correction choice — clusters group by cell identity rather than by condition, '
     'indicating the scVI correction achieved cross-condition alignment without obviously '
     'collapsing distinct states (see UMAP by condition in the figures).\n'
     '\n'
     '{{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, with Immune_All_High.pkl as a second opinion) gave '
     'an initial label for each cluster. Before accepting these, I checked canonical '
     'markers for every final cell type with `check_markers`, which is what caught three '
     'clear mislabels:\n'
     '\n'
     '- **Clusters 0 and 1** were called "Intermediate macrophages" by CellTypist, but '
     'PBMC blood samples do not contain tissue-resident macrophages. Both clusters show '
     'canonical classical-monocyte markers instead ({{gene:CD14+ Monocytes:CD14}} and '
     'strong LYZ/FCN1 expression), so they were relabeled **CD14+ Monocytes**.\n'
     '- **Cluster 4** was called "NK cells", but {{gene:CD8+ T cells:CD3D}} and '
     '{{gene:CD8+ T cells:CD8A}}/CD8B are all clearly enriched in this cluster (unlike the '
     'genuine NK cluster 11, where CD3D is strongly depleted) — real NK cells are '
     'CD3-negative. This cluster is a cytotoxic/effector **CD8+ T cell** population (its '
     'high GNLY/NKG7/GZMB reflects an effector phenotype, not NK identity), relabeled '
     'accordingly.\n'
     '- **Cluster 12** was called a helper-T subset, but it is dominated by '
     '{{gene:Platelets:PPBP}} and {{gene:Platelets:PF4}}, unambiguous platelet markers, so '
     'it was relabeled **Platelets**.\n'
     '\n'
     "All other clusters' CellTypist labels were confirmed by canonical markers and left "
     'unchanged (e.g. {{gene:B cells:MS4A1}} and CD79A for B cells; {{gene:CD16+ NK '
     'cells:GNLY}} and NKG7/KLRD1 for NK cells; FCGR3A/MS4A7 for non-classical monocytes; '
     'TCF4/IL3RA/CLEC4C for pDC; FCER1A/CLEC10A for DC2; CCR7/SELL for the '
     'naive/central-memory T-cell clusters; GZMK/GZMA for the Tem/Trm cytotoxic cluster; '
     'HBB for the erythroid cluster). The small Mast cell cluster is supported by GATA2, a '
     'master mast-cell/basophil transcription factor that was strongly and significantly '
     'enriched (essentially absent outside the cluster); its tryptase gene TPSAB1 trends '
     'in the same direction but does not reach significance given the very small size of '
     'this cluster ({{celltype:Mast cells}} cells).\n'
     '\n'
     'Final cell-type composition:\n'
     '\n'
     '{{table:composition}}\n'
     '\n'
     '## Composition analysis\n'
     '\n'
     'Cell-type proportions were compared per donor (paired, since every donor contributed '
     'both conditions) using a Wilcoxon signed-rank test. **No cell type showed a '
     'significant proportion shift** after correction ({{comp.n_significant}} of '
     '{{comp.n_cell_types}} significant). The largest (non-significant) trends were a '
     'decrease in DC2 proportion in stim ({{comp.dc2.ctrl}} → {{comp.dc2.stim}}, log2 '
     'ratio {{comp.dc2.log2_ratio}}, padj {{comp.dc2.padj}}) and an increase in pDC '
     'proportion ({{comp.pdc.ctrl}} → {{comp.pdc.stim}}, log2 ratio '
     '{{comp.pdc.log2_ratio}}, padj {{comp.pdc.padj}}), but with only {{input.n_donor}} '
     'paired donors the test cannot detect anything short of a large, consistent effect. '
     'This null result is consistent with the short (6h) stimulation window being long '
     'enough to reprogram gene expression but not to substantially alter the mix of '
     'circulating cell types.\n'
     '\n'
     '{{table:composition_test}}\n'
     '\n'
     '## Differential expression\n'
     '\n'
     'Pseudobulk DE (PyDESeq2, design {{de.design}}, paired by donor) was run per cell '
     'type on raw counts. Every cell type with sufficient samples showed a robust '
     'interferon-stimulated gene (ISG) signature: {{gene:CD14+ Monocytes:ISG15}}, '
     '{{gene:CD14+ Monocytes:ISG20}}, and other canonical ISGs ({{gene:CD8+ T cells:IFI6}} '
     'etc.) were significantly upregulated in every tested cell type, essentially '
     'independent of lineage — the expected, textbook response to IFN-β. CD14+ Monocytes '
     'had by far the largest and most complex response (n_significant = '
     '{{de.cd14_monocytes.n_significant}} of {{de.cd14_monocytes.genes_tested}} genes '
     'tested, roughly balanced up/down), consistent with monocytes being major '
     'IFN-responsive antigen-presenting cells; they also induced strong chemokine '
     'responses such as {{gene:CD14+ Monocytes:CCL8}} and IL1RN, not seen as prominently '
     'in lymphocytes. Non-classical monocytes and DC2 also showed large, significant '
     'responses ({{de.non_classical_monocytes.n_significant}} and {{de.dc2.n_significant}} '
     'significant genes respectively). T- and B-lymphocyte subsets showed smaller but '
     'still highly significant ISG inductions (e.g. B cells: {{de.b_cells.n_significant}} '
     'significant genes; CD8+ T cells: {{de.cd8_t_cells.n_significant}}). Downregulated '
     'genes were dominated by ribosomal protein genes (RPL/RPS) and general '
     'translation/iron-storage genes (e.g. FTH1), a common pattern in activated-cell '
     'pseudobulk contrasts and likely reflecting a shift in cellular metabolic state '
     'rather than a specific IFN-driven repression program.\n'
     '\n'
     'Late erythroid, Mast cells, Tem/Trm cytotoxic T cells and pDC were skipped for DE '
     'because too few donors had at least {{de.min_cells}} cells of that type in one or '
     'both conditions — these are the rarest clusters in the dataset and the '
     'composition/DE tests are underpowered for them regardless of biology.\n'
     '\n'
     '{{table:de_summary}}\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- **Condition is fully confounded with the 10x run** (stim and ctrl were captured in '
     'separate captures). Any batch-specific technical effect (e.g. run-to-run efficiency '
     'differences) cannot be distinguished from the true condition effect in this design; '
     'the pseudobulk DE and composition results should be read with this in mind, although '
     'the consistency of the signature with the well-established IFN-β biology, and its '
     'presence across every lineage, makes a purely technical explanation unlikely.\n'
     '- The scVI batch correction used to build the clustering/UMAP embedding included '
     '`condition` in the batch key. This choice helps align cell types across conditions '
     'for annotation and composition testing, but risks merging a condition-specific cell '
     'state into an existing cluster; DE analyses use raw counts and are unaffected, but '
     'any subtler stimulation-induced subpopulation could be underrepresented in the '
     'cluster structure itself.\n'
     '- The mitochondrial-percentage QC metric was uninformative here (no MT genes present '
     'in the object), so a standard cause of low-quality-cell removal was unavailable; '
     'only the gene-count floor was used to catch low-complexity cells.\n'
     '- Very rare clusters (Mast cells, n={{celltype:Mast cells}}; Tem/Trm cytotoxic T '
     'cells, n={{celltype:Tem/Trm cytotoxic T cells}}; Late erythroid, n={{celltype:Late '
     'erythroid}}) rest on small numbers of cells and were excluded from DE/composition '
     'testing; their labels, while marker-supported, should be treated as lower-confidence '
     'than the major lineages.\n'
     '- Composition testing with only {{input.n_donor}} paired donors has limited power; a '
     'true but modest compositional shift could easily go undetected.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'Fourteen PBMC cell types were identified and annotated, with three CellTypist labels '
     'corrected after marker inspection (two monocyte clusters mislabeled as macrophages, '
     'one CD8+ T-cell cluster mislabeled as NK cells, and one platelet cluster mislabeled '
     'as helper T cells). Six-hour IFN-β stimulation did not significantly change the '
     'proportions of any PBMC cell type, but it induced a strong, highly significant, and '
     'pervasive interferon-stimulated gene signature across essentially every cell type, '
     'with the largest transcriptional response in classical (CD14+) monocytes, '
     "non-classical monocytes and DC2 — consistent with these myeloid populations' central "
     'role in mounting and amplifying an innate antiviral response.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
