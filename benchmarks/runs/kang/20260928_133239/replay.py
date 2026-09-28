"""Replay of the agent's tool calls for kang, without the LLM.

Generated: 2026-09-28T13:52:24
Git commit: 24068e592d14aa176f5b6e5c079649c4ed21485a
Source log: tool_calls.jsonl (23 of 24 calls; failed calls omitted)

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
    genes=['CD14', 'FCN1', 'S100A8', 'S100A9', 'LYZ', 'FCGR3A', 'MS4A7', 'CD3D', 'CD3E', 'CD4',
     'CD8A', 'IL7R', 'CCR7', 'SELL', 'NCAM1', 'GNLY', 'NKG7', 'KLRD1', 'PPBP', 'PF4',
     'ITGA2B', 'TPSAB1', 'KIT', 'CPA3', 'CD79A', 'MS4A1', 'LILRA4', 'IL3RA', 'CD1C',
     'FCER1A', 'GZMK'],
    clusters=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '15'],
)

_run(
    "relabel_clusters",
    labels={'0': 'CD14+ Monocytes',
     '1': 'CD14+ Monocytes',
     '4': 'CD8+ effector T cells',
     '12': 'Platelets'},
    reason=('Cluster 0/1: top markers CD14, FCN1, S100A8/S100A9, LYZ are classical (CD14+) '
     'monocyte markers, not a macrophage signature (macrophages are not expected in '
     "blood/PBMC). Cluster 4: CellTypist called it 'NK cells' but CD3D is strongly and "
     'significantly enriched (66.4% of cells express it, log2FC=2.34) along with CD8A '
     '(40.1%, log2FC=4.46), which is incompatible with NK identity (CD3-negative by '
     'definition); the GZMB/NKG7/GNLY signature instead reflects a cytotoxic/effector CD8 '
     'T-cell state. Cluster 12: PPBP (85.1% of cells, log2FC=7.7) and PF4 (63.6%, '
     'log2FC=7.46) are canonical, highly specific platelet markers, clearly contradicting '
     "the 'T helper cell' label."),
)

_run(
    "identify_markers",
    n_genes=25,
)

_run(
    "check_markers",
    genes=['HBB', 'HBA1', 'ALAS2'],
    clusters=['14'],
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
    genes=['ISG15', 'ISG20', 'IFI6', 'IFIT1', 'MX1', 'CXCL10', 'CCL8', 'IL1RN'],
)

_run(
    "get_top_genes",
    cell_type='CD14+ Monocytes',
    direction='up',
    n=15,
)

_run("summarize_findings")

_run(
    "check_markers",
    genes=['MS4A2', 'HDC', 'GATA2', 'CTSG', 'CMA1', 'ENPP3', 'SERPINB1'],
    clusters=['15'],
)

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This dataset (Kang et al., PBMCs from lupus patients) contains cells from '
     '{{input.n_donor}} donors, each split into a control half and an IFN-β-stimulated '
     'half, then pooled and captured as two separate 10x runs — one per condition — before '
     'genotype-based demultiplexing (demuxlet) assigned cells back to donors. Because '
     'condition and 10x run coincide, this confound is carried through the whole analysis '
     '(see Caveats). Starting from {{input.n_cells}} cells x {{input.n_genes}} genes, the '
     'workflow below performs QC, doublet removal, batch-aware integration, clustering, '
     'annotation, and a paired stim-vs-ctrl comparison of composition and expression.\n'
     '\n'
     '## Quality control\n'
     '\n'
     'Gene symbols were detected with a standard `MT-` mitochondrial prefix, but the '
     'search found {{qc.n_mito_genes}} mitochondrial genes in the object — mitochondrial '
     'genes appear to have already been stripped upstream, so `pct_counts_mt` is uniformly '
     '{{qc.mito.max}} for every cell and could not be used as a QC axis. Genes-per-cell '
     'was otherwise well-behaved (median {{qc.genes.median}}, p99 {{qc.genes.p99}}) and '
     'total counts were modest but typical for a short-stimulation dataset (median '
     '{{qc.counts.median}}). The tutorial-default thresholds '
     '(min_genes={{default.min_genes}}, max_pct_mt={{default.max_pct_mt}}, '
     'min_cells={{default.min_cells}}) removed only a small low-quality tail — '
     '{{filter.cells_removed}} cells ({{filter.pct_cells_removed}}) and '
     '{{filter.genes_removed}} rarely-detected genes — so I kept them as-is rather than '
     'tightening further.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     "The user's description states that demuxlet genotype demultiplexing already removed "
     'most cross-donor doublets, leaving only doublets formed from two cells of the *same* '
     'donor (which genotype cannot distinguish). This matches the documented exception to '
     'the standard doublet workflow: a strict median+3×MAD cut would treat this '
     'already-cleaned distribution as if it still contained the full doublet load and '
     'remove real cells from the upper score tail. Scrublet was run separately per 10x run '
     '(`condition`, the true technical grouping — donors were pooled *within* each run, so '
     '`donor` is not a separate-run marker). The score distribution was not bimodal '
     '({{doublet.distribution}}), consistent with most obvious doublets already being '
     'gone. I used the higher of the two per-run automatic thresholds '
     '({{doublet.threshold}}, from the stim run) as a light-touch cutoff, removing only '
     '{{doublet.cells_removed}} cells ({{doublet.pct_removed}}) — the residual same-donor '
     'doublets are a small, largely irreducible fraction, not worth aggressive filtering.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'Condition is fully confounded with 10x run, and IFN-β stimulation is known to '
     'produce a strong, genome-wide transcriptional shift (interferon-stimulated genes). '
     'Left uncorrected, this can split each cell type into two condition-specific '
     'clusters, undermining both annotation and the composition test. I therefore ran scVI '
     'with batch_key = [`donor`, `condition`] (one batch per donor-condition combination), '
     'so the latent space and downstream clustering/UMAP align matching cell types across '
     'stim and ctrl while raw counts — used later for pseudobulk DE and composition — are '
     'left untouched. The risk of this choice is over-correction: a cell state that only '
     'exists under stimulation could be merged into its resting counterpart. I checked for '
     'this by confirming that ISG expression differences are still detected as strong, '
     'significant, and cell-type-specific in the DE step below, i.e. the correction '
     'integrated cell identity without erasing the condition effect.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering at resolution {{cluster.resolution}} on the scVI embedding '
     'produced {{cluster.n}} clusters, ranging from {{cluster.smallest}} to '
     '{{cluster.largest}} cells (see {{table:clusters}}).\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, majority vote per cluster) gave an initial '
     'fine-grained annotation, cross-checked against the coarser Immune_All_High.pkl model '
     'and canonical markers for every cluster. Three annotations were overridden:\n'
     '\n'
     '- **Clusters 0 and 1** were called "Intermediate macrophages" by CellTypist, but '
     'their top markers (CD14, FCN1, S100A8/S100A9, LYZ, all strongly and specifically '
     'enriched) are a classical-monocyte signature, not a tissue-macrophage one, which is '
     'not an expected PBMC population. Relabeled to {{celltype:CD14+ Monocytes}}.\n'
     '- **Cluster 4** was called "NK cells", but CD3D is strongly and specifically '
     'enriched in this cluster (checked via check_markers) together with CD8A — NK cells '
     'are CD3-negative by definition, so this is a {{celltype:CD8+ effector T cells}} '
     'population instead, with GZMB/NKG7/GNLY reflecting cytotoxic activity rather than NK '
     'identity.\n'
     '- **Cluster 12** was called a naive/Tcm helper T-cell cluster, but its top markers '
     'are PPBP and PF4 ({{gene:Platelets:PPBP}} and {{gene:Platelets:PF4}}), textbook, '
     'highly specific platelet markers. Relabeled to {{celltype:Platelets}}.\n'
     '\n'
     'All remaining clusters agreed with CellTypist and canonical markers: classical '
     'monocyte markers absent, FCGR3A/MS4A7 near-universal in {{celltype:Non-classical '
     'monocytes}}; CD79A/MS4A1 enriched in {{celltype:B cells}}; GNLY/NKG7 near-universal '
     'with CD3D absent in {{celltype:CD16+ NK cells}}; LILRA4/IL3RA enriched in '
     '{{celltype:pDC}}; CD1C/FCER1A/HLA-DR enriched in {{celltype:DC2}}; hemoglobin genes '
     'HBB and HBA1 (near-ubiquitous detection, large fold-change) in {{celltype:Late '
     'erythroid}}; CCR7/SELL/IL7R enrichment distinguishing naive/central-memory '
     '({{celltypes:Tcm/Naive helper T cells|Tcm/Naive cytotoxic T cells}}) from '
     'effector/cytotoxic subsets ({{celltypes:Tem/Effector helper T cells|Tem/Trm '
     'cytotoxic T cells}}), the latter also carrying GZMK/GZMA/TIGIT/LAG3; and the '
     'transcription factor GATA2, a canonical mast-cell/basophil-lineage marker, is '
     'significantly and strongly enriched in the small {{celltype:Mast cells}} cluster, '
     'consistent with (though not as specific as) its co-enrichment of the tryptase gene '
     'TPSAB1. Final cell-type composition is in {{table:composition}}, totalling '
     '{{final.n_cells}} cells across {{celltype.n}} types after QC and doublet filtering '
     '({{final.pct_removed}} of input cells removed overall).\n'
     '\n'
     '## Composition analysis (stim vs ctrl)\n'
     '\n'
     'Cell-type proportions were compared per donor (paired Wilcoxon signed-rank test, all '
     '{{input.n_donor}} donors have both conditions) with Benjamini-Hochberg correction; '
     'see {{table:composition}} and {{table:composition_test}}. No cell type showed a '
     'significant proportion shift ({{comp.n_significant}} of {{comp.n_cell_types}} '
     'significant after correction). The largest (non-significant) trend was a decrease in '
     '{{celltype:DC2}} under stimulation (ctrl {{comp.dc2.ctrl}} vs stim '
     '{{comp.dc2.stim}}, padj {{comp.dc2.padj}}), and a non-significant increase in '
     '{{celltype:pDC}} (ctrl {{comp.pdc.ctrl}} vs stim {{comp.pdc.stim}}, padj '
     "{{comp.pdc.padj}}). With only {{input.n_donor}} paired donors the test's minimum "
     'attainable p-value is limited, so this is more consistent with "no detectable '
     'compositional shift over the short stimulation window" than definitive proof of no '
     'change — a longer stimulation might reveal shifts that this protocol does not.\n'
     '\n'
     '## Differential expression\n'
     '\n'
     'Pseudobulk DE (PyDESeq2, design `~donor + condition`, donor as a paired covariate) '
     'was run per cell type on raw counts; {{de.n_cell_types_tested}} cell types had '
     'enough samples per condition to test (skipped: {{de.skipped}}, all too rare or too '
     'unevenly captured for pseudobulk after the min-cells-per-sample filter). Results are '
     'summarized in {{table:de_summary}}.\n'
     '\n'
     'The dominant signal in every tested cell type is a canonical, coherent type-I '
     'interferon response: ISG15, ISG20, IFI6, IFIT1 and MX1 are significantly upregulated '
     'with large, consistent log2 fold-changes across essentially all cell types, e.g. '
     '{{gene:CD14+ Monocytes:ISG15}} in monocytes, {{gene:CD16+ NK cells:ISG15}} in NK '
     'cells, {{gene:B cells:ISG15}} in B cells, and {{gene:Tcm/Naive helper T '
     'cells:ISG15}} in naive/Tcm CD4 T cells, alongside {{gene:CD14+ Monocytes:IFIT1}} and '
     '{{gene:CD14+ Monocytes:MX1}}. This is the expected biology for IFN-β stimulation and '
     "confirms the ISG signal survived scVI's batch correction rather than being "
     'integrated away.\n'
     '\n'
     'Beyond the shared ISG core, {{celltype:CD14+ Monocytes}} show the largest and most '
     'cell-type-specific response ({{de.cd14_monocytes.n_significant}} of '
     '{{de.cd14_monocytes.genes_tested}} genes tested significant), including strong '
     'induction of the inflammatory chemokine {{gene:CD14+ Monocytes:CCL8}} and the '
     'anti-inflammatory decoy receptor antagonist {{gene:CD14+ Monocytes:IL1RN}}, both far '
     'beyond the generic ISG module — consistent with monocytes acting as major '
     'amplifiers/responders of the IFN-β response. {{celltype:Platelets}}, by contrast, '
     'show comparatively few significant genes ({{de.platelets.n_significant}} of '
     '{{de.platelets.genes_tested}}), predominantly the core ISGs, reflecting both a '
     'genuine muted response and reduced power from a small pseudobulk sample.\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- **Condition/run confound**: stim and ctrl cells were captured in two separate 10x '
     'runs, so condition is fully confounded with any run-level technical batch effect. '
     'scVI integration (batch_key = donor + condition) aligns cell identities across the '
     'two runs for clustering/annotation, but it cannot separate a true biological IFN-β '
     'effect from a run-specific technical shift in the pseudobulk DE test; the paired '
     'donor design (each donor contributes both conditions) is the main mitigation, since '
     'a donor-specific batch artifact would need to systematically covary with condition '
     'to create a false signal.\n'
     '- **Rare types excluded from DE**: {{celltype:Mast cells}}, {{celltype:pDC}}, '
     '{{celltype:Tem/Trm cytotoxic T cells}}, and {{celltype:Late erythroid}} either '
     'lacked enough cells per donor per condition for pseudobulk, or (mast cells) had too '
     'few cells overall for fully confident marker-based annotation in the first place — '
     'treat these labels as provisional.\n'
     '- **Ambiguous transcriptional state cluster**: cluster 3 (part of '
     '{{celltype:Tcm/Naive helper T cells}}) is defined more by an '
     'immediate-early/activation signature (CD69, BTG1, immediate-early transcripts) than '
     'by strong CD3D enrichment specifically, though IL7R/CCR7/SELL enrichment and absence '
     'of other lineage markers support a CD4 T-cell identity; it may represent an '
     'activated/stressed T-cell state rather than a distinct resting subset.\n'
     '- **No mitochondrial QC**: mitochondrial genes were absent from the gene list, so '
     'the usual pct-mt-based dying-cell filter could not be applied; QC relied on '
     'gene/count count filters only.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'After QC, light doublet removal, and scVI-based integration across donors and the '
     'condition-confounded 10x runs, {{final.n_cells}} cells resolved into {{celltype.n}} '
     'annotated PBMC types spanning monocytes, DCs, B, T, NK, platelet and erythroid '
     'populations. The stimulation did not measurably shift cell-type proportions in this '
     'paired cohort, but it induced a strong, broadly shared interferon-stimulated gene '
     'program across virtually every cell type, with monocytes mounting the largest and '
     'most distinct additional inflammatory response (chemokines, IL1RN) on top of that '
     'shared core.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
