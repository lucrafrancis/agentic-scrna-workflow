"""Replay of the agent's tool calls for kang, without the LLM.

Generated: 2026-09-28T14:30:35
Git commit: 24068e592d14aa176f5b6e5c079649c4ed21485a
Source log: tool_calls.jsonl (20 of 21 calls; failed calls omitted)

Re-runs the same tool functions with the same arguments and seed. Output goes to
outputs/kang_replay/, with its own tool log and report. Package versions are pinned by
uv.lock at the commit above: check it out first if the code has changed since.

Caveats: annotate_celltypes downloads CellTypist models, and scVI training can differ
slightly between hardware (CPU vs GPU), which can change clusters downstream.

  uv run python examples/kang/replay.py
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
    genes=['CD14', 'FCN1', 'LYZ', 'S100A8', 'CD68', 'CD163', 'FCGR3A', 'MS4A7', 'CD3D', 'CD3E',
     'CD8A', 'CD8B', 'CD4', 'NKG7', 'GNLY', 'KLRD1', 'FGFBP2', 'PPBP', 'PF4', 'GNG11',
     'TUBB1', 'TPSAB1', 'CPA3', 'KIT', 'MS4A1', 'CD79A', 'IL3RA', 'GZMB', 'TCF4', 'HBB'],
    clusters=['0', '1', '2', '4', '6', '9', '10', '11', '12', '13', '14', '15'],
)

_run(
    "check_markers",
    genes=['CD3D', 'CD3E', 'CD4', 'IL7R', 'CCR7', 'SELL', 'CD8A', 'FOXP3', 'MS4A2', 'FCER1A',
     'GATA2'],
    clusters=['3', '5', '7', '8', '15'],
)

_run(
    "relabel_clusters",
    labels={'0': 'CD14+ Monocytes',
     '1': 'CD14+ Monocytes',
     '4': 'CD8+ Effector T cells',
     '12': 'Platelets'},
    reason=("Cluster 0/1: CellTypist called these 'Intermediate macrophages', but macrophages are "
     'tissue-resident and not expected in PBMC blood samples; markers show classical '
     'monocyte identity (CD14 39.6%/64.0%, FCN1 78.9%/54.3%, LYZ 84.1%/72.3%, S100A8 '
     "60.6%/77.5%, all highly significant). Cluster 4: labeled 'NK cells' but strongly "
     'expresses CD3D (66.4%, padj~0), CD3E, CD8A (40.1%) and CD8B (30.7%), i.e. a bona '
     'fide T-cell receptor complex, unlike the true NK cluster 11 (CD3D only 4.5%); this '
     'is a CD3+CD8+ cytotoxic/effector T-cell population despite high GZMB/NKG7/GNLY. '
     'Cluster 12: labeled a T-helper subset but its top markers are PPBP (85.1%), PF4 '
     '(63.6%), GNG11 (61.6%), TUBB1 (36%) - an unambiguous platelet/megakaryocyte '
     'signature.'),
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
    genes=['ISG15', 'ISG20', 'IFIT1', 'MX1', 'CCL8', 'IL1RN', 'APOBEC3A', 'B2M', 'LY6E', 'IFI6'],
)

_run("summarize_findings")

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This analysis processed PBMCs from {{input.n_donor}} SLE patients, each split into '
     'interferon-beta–stimulated and control halves, captured as two separate 10x runs '
     '(`condition` = ctrl/stim, one run per condition, {{input.n_batches}} donors pooled '
     'per run). Starting from {{input.n_cells}} cells x {{input.n_genes}} genes, the '
     'workflow proceeded through QC filtering, light-touch doublet removal, scVI-based '
     'batch-corrected clustering, CellTypist annotation with marker-based correction, and '
     'paired composition/differential-expression testing between stimulated and control '
     'samples across donors.\n'
     '\n'
     '## Quality control\n'
     '\n'
     'Gene-identifier detection found standard human gene symbols, but no mitochondrial '
     'genes (`MT-` prefix) were present in the object at all ({{qc.n_mito_genes}} found), '
     'so pct-mitochondrial is {{qc.mito.median}} for every cell — mitochondrial genes were '
     'evidently removed upstream, before this dataset was provided. This makes the '
     'mitochondrial filter a no-op safeguard rather than a real QC lever here.\n'
     '\n'
     'Genes-per-cell ranged from a minimum of {{qc.genes.min}} to a maximum of '
     '{{qc.genes.max}} (median {{qc.genes.median}}), and total counts per cell from '
     '{{qc.counts.min}} to {{qc.counts.max}} (median {{qc.counts.median}}). The '
     'recommended tutorial-default thresholds (min_genes={{filter.min_genes}}, '
     'max_pct_mt={{filter.max_pct_mt}}, min_cells={{filter.min_cells}}) were adopted '
     'as-is: they removed only {{filter.cells_removed}} cells '
     '({{filter.pct_cells_removed}}) and {{filter.genes_removed}} genes, a light and '
     'appropriate trim for an already reasonably clean dataset.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'Demuxlet genotype-based demultiplexing was used upstream to assign cells to donors, '
     'which removes doublets formed by two cells from different patients — but '
     'same-patient doublets survive undetected by genotype. Scrublet was run separately '
     'per 10x run (`condition`, since stim and ctrl are the two distinct captures), giving '
     'per-batch automatic thresholds of {{doublet.scrublet_auto}} (ctrl / stim). The score '
     'distribution was not bimodal ({{doublet.distribution}}), and because most doublets '
     'were already removed upstream, the standard median+3*MAD rule '
     '({{doublet.median_3mad}}) would flag a much larger fraction of cells and '
     'over-aggressively cut real cells from the upper tail of this already-cleaned '
     'distribution. Following the light-touch guidance for pre-demultiplexed data, the '
     'higher of the two per-batch Scrublet automatic thresholds ({{doublet.threshold}}) '
     'was applied dataset-wide as a conservative cutoff, removing only '
     '{{doublet.cells_removed}} cells ({{doublet.pct_removed}}) — consistent with the '
     'expectation that few doublets remained.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'Condition and 10x run are the same technical factor in this dataset, and IFN-β '
     'stimulation is known to induce a strong, broad transcriptional shift (the '
     'interferon-stimulated gene program) that risks splitting each cell type into '
     'separate stim/ctrl clusters if left uncorrected. To align shared cell types across '
     'conditions for consistent annotation while leaving raw counts untouched for the '
     'downstream DE and composition tests, scVI was trained with a combined batch key of '
     'donor and condition (one batch per sample, {{scvi.n_latent}}-dimensional latent '
     'space, {{scvi.epochs}} epochs) rather than plain PCA. The resulting embedding was '
     'used for neighbor graph, clustering and UMAP.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering at resolution {{cluster.resolution}} on the scVI embedding yielded '
     '{{cluster.n}} clusters, ranging from {{cluster.smallest}} to {{cluster.largest}} '
     'cells. {{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, majority vote per cluster) provided initial '
     'fine-grained labels, cross-checked against the coarser Immune_All_High.pkl model and '
     'canonical markers for every final cell type.\n'
     '\n'
     'Three systematic issues were corrected:\n'
     '\n'
     '- **Clusters 0 and 1** were called "Intermediate macrophages" by both CellTypist '
     'models, but macrophages are tissue-resident and not expected in blood. Markers '
     'instead showed a clear classical-monocyte signature (CD14, FCN1, LYZ, S100A8 all '
     'strongly enriched), so these were relabeled to CD14+ Monocytes, giving '
     '{{celltype:CD14+ Monocytes}} cells in this type.\n'
     '- **Cluster 4** was labeled "NK cells," but it strongly expressed the T-cell '
     'receptor complex genes CD3D and CD3E along with CD8A/CD8B, unlike the true NK '
     'cluster (cluster 11), which lacks CD3D. It was relabeled to CD8+ Effector T cells '
     '(retaining high GZMB/NKG7 as an activated/cytotoxic phenotype), giving '
     '{{celltype:CD8+ Effector T cells}} cells.\n'
     '- **Cluster 12** was labeled a naive T-helper subset, but its dominant markers '
     '(PPBP, PF4, GNG11, TUBB1) are an unambiguous platelet/megakaryocyte signature; '
     'relabeled to Platelets, giving {{celltype:Platelets}} cells.\n'
     '\n'
     'All other clusters agreed with CellTypist and canonical markers: Non-classical '
     'monocytes (FCGR3A, MS4A7), B cells (CD79A, MS4A1), CD16+ NK cells (GNLY, NKG7, '
     'KLRD1, FGFBP2, CD3D-negative), pDC (IL3RA, TCF4, GZMB), DC2 (LYZ, HLA-DR), '
     'naive/memory helper and cytotoxic T-cell subsets (distinguished by CCR7/SELL/IL7R vs '
     'CD8A), and Late erythroid (HBB detected in essentially every cell of the cluster). '
     'Cluster 15 ("Mast cells", only {{celltype:Mast cells}} cells) showed a significant '
     'GATA2 enrichment and a TPSAB1/FCER1A trend, while CD3D/CD3E were undetectable, '
     'ruling out the alternative T-cell label from the second model — this identification '
     'should be treated cautiously given the very small cluster size.\n'
     '\n'
     'Final cell-type composition: {{table:composition}}\n'
     '\n'
     '## Composition analysis\n'
     '\n'
     'Cell-type proportions were compared between stim and ctrl per donor (paired Wilcoxon '
     'signed-rank test, {{comp.test}}, {{input.n_donor}} donors per condition, '
     'Benjamini-Hochberg corrected). {{comp.n_significant}} of {{comp.n_cell_types}} cell '
     'types reached significance after correction. {{table:composition_test}}\n'
     '\n'
     'This null result is consistent with the short stimulation window used in this '
     'experiment: IFN-β reprograms gene expression rapidly but is not expected to alter '
     'the balance of major PBMC lineages so quickly, and the compositional nature of these '
     'proportions (a shift in one type mechanically moves the others) combined with only '
     '{{input.n_donor}} paired donors limits power to detect small shifts.\n'
     '\n'
     '## Differential expression\n'
     '\n'
     'Pseudobulk DE (PyDESeq2, design {{de.design}}, donor included as a paired covariate) '
     'was run per cell type on raw counts. {{de.skipped}} were skipped for having too few '
     'samples with sufficient cells per condition. {{table:de_summary}}\n'
     '\n'
     'Across every tested cell type, the top upregulated genes in stim vs ctrl were '
     'canonical type-I interferon-stimulated genes: {{gene:CD14+ Monocytes:ISG15}} in '
     'CD14+ Monocytes, {{gene:B cells:ISG15}} in B cells, {{gene:CD16+ NK cells:ISG15}} in '
     'CD16+ NK cells, and {{gene:CD8+ Effector T cells:ISG15}} in CD8+ Effector T cells, '
     'alongside consistent induction of {{gene:CD14+ Monocytes:ISG20}}, {{gene:CD14+ '
     'Monocytes:IFIT1}}, {{gene:CD14+ Monocytes:MX1}}, {{gene:CD14+ Monocytes:LY6E}} and '
     '{{gene:CD14+ Monocytes:IFI6}} in monocytes, and matching induction of the same genes '
     'in every other tested type (e.g. {{gene:Tcm/Naive helper T cells:ISG15}} in '
     'Tcm/Naive helper T cells, {{gene:DC2:ISG15}} in DC2, {{gene:Platelets:ISG15}} in '
     'Platelets). MHC-I component {{gene:CD14+ Monocytes:B2M}} was also consistently '
     'upregulated. CD14+ Monocytes additionally showed strong, more monocyte-specific '
     'induction of the chemokine {{gene:CD14+ Monocytes:CCL8}}, the anti-inflammatory '
     'cytokine antagonist {{gene:CD14+ Monocytes:IL1RN}}, and the antiviral deaminase '
     '{{gene:CD14+ Monocytes:APOBEC3A}} (also induced in Non-classical monocytes: '
     '{{gene:Non-classical monocytes:APOBEC3A}}), reflecting the well-documented '
     'sensitivity of monocytes to IFN-β. CD14+ Monocytes showed by far the largest number '
     'of significant genes ({{de.cd14_monocytes.n_significant}} of '
     '{{de.cd14_monocytes.genes_tested}} tested), consistent with monocytes being a '
     'principal IFN-responsive population, while smaller lymphocyte and platelet '
     'populations showed proportionally fewer significant genes, partly reflecting lower '
     'statistical power from fewer cells/pseudobulk counts.\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- **Condition is fully confounded with the 10x run**: stim and ctrl PBMCs were '
     'captured as two separate captures, so any residual run-specific technical effect '
     '(loading, reagent lot, capture efficiency) cannot be distinguished from the true '
     'biological IFN-β response in this DE analysis. The paired donor design and '
     'consistent, biologically sensible ISG signature across all cell types make a purely '
     'technical explanation unlikely, but this cannot be formally ruled out.\n'
     '- Mitochondrial-based QC was uninformative here because mitochondrial genes were '
     'absent from the input matrix, so this dataset lacks that avenue for identifying '
     'stressed/dying cells.\n'
     '- The scVI batch correction (donor x condition) was chosen to align shared cell '
     'types across stim/ctrl, which is necessary given the strong, ubiquitous ISG program, '
     'but a joint donor+condition batch key risks over-merging a genuinely '
     'condition-specific cell state into an existing type; the UMAP colored by condition '
     'should be inspected to confirm mixing without erasing a true stimulated-only '
     'subpopulation.\n'
     '- Very small clusters (Platelets, pDC, Late erythroid, Tem/Trm cytotoxic T cells, '
     'Mast cells) have limited statistical power for both composition and DE testing, and '
     'several were skipped from DE entirely for having too few samples with sufficient '
     'cells per condition.\n'
     '- The Mast cell label rests on a single significant marker (GATA2) in a cluster of '
     'only {{celltype:Mast cells}} cells, so this identity should be treated as '
     'tentative.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'Clustering on a donor+condition-corrected scVI embedding resolved {{final.n_cells}} '
     'cells into {{cluster.n}} clusters spanning the expected major PBMC lineages (CD14+ '
     'and non-classical monocytes, DCs, B cells, NK cells, several CD4/CD8 T-cell memory '
     'subsets, platelets, and a small erythroid/mast contamination), with three clusters '
     'requiring correction of biologically implausible or contradicted CellTypist labels '
     'before analysis. IFN-β stimulation did not significantly shift cell-type proportions '
     'in any lineage, but it strongly and consistently reprogrammed gene expression across '
     'every cell type via the same core ISG module, with monocytes showing the broadest '
     'and most quantitatively pronounced response and additional monocyte-specific '
     'chemokine/cytokine induction — a pattern well aligned with the known systemic action '
     'of type-I interferon on peripheral immune cells.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
