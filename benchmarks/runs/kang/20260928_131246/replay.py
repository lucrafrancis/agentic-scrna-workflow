"""Replay of the agent's tool calls for kang, without the LLM.

Generated: 2026-09-28T13:32:37
Git commit: 24068e592d14aa176f5b6e5c079649c4ed21485a
Source log: tool_calls.jsonl (23 of 25 calls; failed calls omitted)

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
    genes=['CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'IL7R', 'CCR7', 'SELL', 'NKG7', 'GNLY', 'KLRD1',
     'CD14', 'LYZ', 'FCGR3A', 'MS4A7', 'MS4A1', 'CD79A', 'CD19', 'PPBP', 'PF4', 'GP9',
     'TPSAB1', 'CPA3', 'KIT', 'HBB', 'GATA1', 'CD74', 'HLA-DRA', 'FCER1A', 'CLEC9A',
     'IL3RA', 'GZMB', 'GZMK', 'TCF7', 'LEF1', 'S100A8', 'S100A9', 'ISG15', 'IFI6', 'MX1'],
)

_run(
    "check_markers",
    genes=['TRAC', 'TRBC1', 'CD3G', 'NCAM1', 'FCER1G', 'KLRB1', 'CD2'],
    clusters=['4', '6', '11'],
)

_run(
    "relabel_clusters",
    labels={'0': 'CD14+ Monocytes',
     '1': 'CD14+ Monocytes',
     '4': 'CD8+ Effector T cells',
     '12': 'Platelets'},
    reason=("Cluster 0/1: CellTypist called these 'Intermediate macrophages', but "
     'CD14/LYZ/S100A8/S100A9 enrichment with high FCGR3A/MS4A7 marks classical blood '
     'monocytes, not tissue macrophages; the two clusters differ mainly in IFN-response '
     'genes (ISG15/IFI6/MX1), reflecting stim vs ctrl states of the same cell type rather '
     "than distinct types. Cluster 4: CellTypist called it 'NK cells', but "
     'CD3D/CD3E/CD3G/CD2 are all significantly enriched (up to 66% of cells) while FCER1G '
     '(an NK/myeloid marker) is significantly depleted (2.9% vs 37.5% elsewhere); with '
     'CD8A/CD8B enrichment this is cytotoxic/effector CD8+ T cells, matching the '
     "Immune_All_High second opinion of 'T cells'. Cluster 12: CellTypist called it a "
     'T-cell subtype, but PPBP, PF4, and GP9 are massively and specifically enriched (up '
     'to 85% of cells, log2FC>7) with near-zero expression elsewhere, the canonical '
     'platelet/megakaryocyte signature.'),
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
    genes=['ISG15', 'ISG20', 'IFI6', 'MX1', 'IFIT1', 'CCL8', 'IL1RN', 'CXCL10', 'STAT1'],
)

_run(
    "get_top_genes",
    cell_type='CD14+ Monocytes',
    direction='up',
    n=10,
)

_run(
    "get_top_genes",
    cell_type='B cells',
    direction='up',
    n=10,
)

_run("summarize_findings")

_run(
    "check_markers",
    genes=['MS4A2', 'GATA2', 'HDC', 'TPSB2', 'CTSG', 'SLC18A2', 'VWA5A', 'PRG2'],
    clusters=['15'],
)

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This analysis processes PBMCs from SLE patients ({{input.n_donor}} donors), each '
     'split into a control and an IFN-β-stimulated (6h) sample, starting from '
     '{{input.n_cells}} raw cells / {{input.n_genes}} genes and ending with '
     '{{final.n_cells}} cells / {{final.n_genes}} genes assigned to {{celltype.n}} cell '
     'types. Because the two conditions were captured in separate 10x runs, `condition` is '
     'completely confounded with the technical batch — this shapes several downstream '
     'choices (dimensionality reduction, and a caveat for DE) discussed below.\n'
     '\n'
     '## Quality control\n'
     '\n'
     'No mitochondrial genes ({{qc.n_mito_genes}} found) survived in this matrix, so the '
     '%-mito metric is uninformative here and could not be used to flag stressed cells; '
     'the mito filter (max_pct_mt={{filter.max_pct_mt}}) was left at the standard tutorial '
     'default purely as a no-op safety net. Filtering on the standard tutorial floor '
     '(min_genes={{filter.min_genes}}, min_cells={{filter.min_cells}}) removed only '
     '{{filter.cells_removed}} cells ({{filter.pct_cells_removed}}) and '
     '{{filter.genes_removed}} genes — this dataset was already fairly clean, consistent '
     'with genotype-based (demuxlet) cell assignment having already discarded '
     'low-quality/ambiguous droplets upstream.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     "The user's protocol notes that demuxlet genotyping already removed cross-donor "
     'doublets; only doublets formed by two cells of the *same* donor (transcriptionally '
     'indistinguishable from singlets to Scrublet) can remain. Consistent with this, the '
     'doublet-score distribution was unimodal ({{doublet.distribution}}), so the standard '
     'median+3×MAD rule ({{doublet.median_3mad}}) would over-flag genuine cells. I instead '
     'used the light-touch approach: Scrublet run separately per 10x run (= per '
     'condition), taking the higher of the two per-run automatic thresholds '
     '({{doublet.scrublet_auto}}) as the cutoff ({{doublet.threshold}}). This removed only '
     '{{doublet.cells_removed}} cells ({{doublet.pct_removed}}), appropriate given most '
     'doublets were already excised upstream.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     '`condition` and 10x run are the same variable here, and IFN-β is known to induce a '
     'strong, cell-type-spanning transcriptional shift that can otherwise split each true '
     'cell type into per-condition clusters. To align cell types across conditions while '
     'still allowing donor-level batch structure to be corrected, I ran scVI with batch '
     'key = donor × condition (one batch per sample, {{scvi.n_latent}}-dimensional latent '
     'space) rather than plain PCA. Raw counts were untouched, so pseudobulk DE and '
     'composition below still reflect real biology, not the correction.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering (resolution {{cluster.resolution}}) on the scVI embedding gave '
     '{{cluster.n}} clusters, sized from {{cluster.smallest}} to {{cluster.largest}} cells '
     '(full table below).\n'
     '\n'
     '{{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low, cross-checked against Immune_All_High) gave an initial '
     'per-cluster majority call. Checking canonical markers against both labels surfaced '
     'three clear contradictions, which were corrected:\n'
     '\n'
     '- **Clusters 0 and 1** were called "Intermediate macrophages", but they are '
     'CD14-high, LYZ-high, S100A8/S100A9-high classical **monocytes** ({{gene:CD14+ '
     'Monocytes:CD14}} confirms strong CD14 enrichment) — "macrophage" is a tissue label '
     'misapplied to blood monocytes. The two clusters differ mainly in '
     'interferon-stimulated genes rather than cell identity (a residual condition effect '
     'on the strongly IFN-responsive monocyte compartment, discussed in Caveats), so both '
     'were relabeled **CD14+ Monocytes**.\n'
     '- **Cluster 4** was called "NK cells", but CD3D/CD3E/CD3G/CD2 are all significantly '
     'enriched (expressed in the clear majority of cells) while FCER1G, an NK/myeloid '
     'marker, is significantly *depleted* relative to the rest of the data — the opposite '
     'of what an NK population should show. Combined with CD8A/CD8B enrichment, this is a '
     'cytotoxic/effector **CD8+ T cell** cluster (matching the Immune_All_High second '
     'opinion of "T cells"), so it was relabeled accordingly.\n'
     '- **Cluster 12** was called a helper-T-cell subtype, but PPBP and PF4 are massively '
     'and specifically enriched (expressed in most cells of the cluster, essentially '
     'absent elsewhere) — the canonical platelet signature — so it was relabeled '
     '**Platelets**.\n'
     '\n'
     'All final cell types were confirmed with canonical markers checked directly against '
     'the ranked marker tables, requiring significant enrichment versus the rest of the '
     'data: classical monocytes by {{gene:CD14+ Monocytes:LYZ}}; non-classical monocytes '
     'by FCGR3A/MS4A7 enrichment; genuine NK cells (cluster 11, unaffected by the '
     'relabeling) by GNLY/NKG7/KLRD1 enrichment together with CD3D absence and FCER1G '
     'presence; B cells by CD79A/MS4A1/CD19 enrichment; pDCs by IL3RA/HLA-DRA/GZMB; DC2 by '
     'HLA-DRA/CD74/FCER1A; platelets by PPBP/PF4/GP9; erythroid cells by near-universal '
     'HBB expression; mast cells by the transcription factor GATA2, which is essentially '
     'restricted to this tiny cluster and significantly enriched there (with T/B/monocyte '
     'lineage markers all absent), consistent with a mast cell/basophil-lineage identity; '
     'and naive/memory T-cell subsets by the expected CCR7/SELL/LEF1 (naive-like) or '
     'GZMK/LAG3/TIGIT (effector-like) combined with CD8A/CD8B presence or absence to '
     'separate CD4 from CD8 lineages.\n'
     '\n'
     'Final cell-type sizes: {{celltype:CD14+ Monocytes}} CD14+ monocytes, '
     '{{celltype:Tcm/Naive helper T cells}} naive/central-memory CD4 T cells, '
     '{{celltype:Tem/Effector helper T cells}} effector/memory CD4 T cells, {{celltype:B '
     'cells}} B cells, {{celltype:CD8+ Effector T cells}} effector CD8 T cells, '
     '{{celltype:CD16+ NK cells}} NK cells, {{celltype:Tcm/Naive cytotoxic T cells}} naive '
     'CD8 T cells, {{celltype:Non-classical monocytes}} non-classical monocytes, '
     '{{celltype:DC2}} conventional DCs, {{celltype:Platelets}} platelets, '
     '{{celltype:pDC}} plasmacytoid DCs, {{celltype:Late erythroid}} erythroid cells, '
     '{{celltype:Tem/Trm cytotoxic T cells}} tissue-resident/effector-memory CD8 T cells, '
     'and {{celltype:Mast cells}} mast cells.\n'
     '\n'
     '## Composition analysis\n'
     '\n'
     'Comparing stim vs ctrl proportions per donor ({{comp.test}}) found '
     '{{comp.n_significant}} significant shifts after multiple-testing correction, out of '
     '{{comp.n_cell_types}} cell types tested.\n'
     '\n'
     '{{table:composition_test}}\n'
     '\n'
     '{{table:composition}}\n'
     '\n'
     "No cell type's proportion changed significantly between conditions. This matches the "
     'expected biology of a short (6h) cytokine challenge: IFN-β reprograms '
     'transcriptional state broadly across cell types without driving proliferation, '
     'death, or migration on this timescale. The largest (non-significant) trends were a '
     'relative increase in pDC abundance and a decrease in DC2 abundance under stimulation '
     "— plausible given IFN's known effects on dendritic cell subsets — but with few "
     'donors and compositional (non-independent) proportions, this analysis is '
     "underpowered for effects smaller than a large fold-change (note the paired test's "
     'floor p-value).\n'
     '\n'
     '## Differential expression\n'
     '\n'
     'Pseudobulk DE (PyDESeq2, design {{de.design}}, paired by donor) was run per cell '
     'type; a few rare cell types were skipped for having too few cells/samples '
     '({{de.skipped}}). Across all {{de.n_cell_types_tested}} tested cell types, a '
     'strikingly consistent interferon-stimulated gene (ISG) signature dominates the '
     'stim-vs-ctrl response: {{gene:CD14+ Monocytes:ISG15}} in monocytes, {{gene:B '
     'cells:ISG15}} in B cells, {{gene:CD16+ NK cells:ISG15}} in NK cells, and {{gene:CD8+ '
     'Effector T cells:ISG15}} in CD8 T cells — all strongly and significantly '
     'upregulated, alongside ISG20, IFI6, MX1, IFIT1 and LY6E in essentially every '
     'lineage. This is exactly the expected direct transcriptional consequence of IFN-β '
     'signaling through the JAK-STAT/ISGF3 pathway, and its reproducibility across '
     'independent cell types is a strong internal consistency check on the stimulation '
     'itself.\n'
     '\n'
     '{{table:de_summary}}\n'
     '\n'
     'CD14+ monocytes show by far the largest response '
     '({{de.cd14_monocytes.n_significant}} significant genes out of '
     '{{de.cd14_monocytes.genes_tested}} tested), consistent with monocytes being primary '
     'IFN-β responders and antigen-presenting/inflammatory relay cells; top induced genes '
     'include chemokine {{gene:CD14+ Monocytes:CCL8}} and the anti-inflammatory decoy '
     'receptor antagonist {{gene:CD14+ Monocytes:IL1RN}}, alongside the shared ISG module. '
     'B cells also show a robust response ({{gene:B cells:ISG20}}), with induction of '
     'antiviral effectors such as {{gene:B cells:IFI6}}. Downregulated genes are dominated '
     'by ribosomal-protein transcripts across most cell types, likely reflecting a shift '
     'away from steady-state protein synthesis under acute cytokine stress rather than a '
     'specific IFN target program.\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- **Condition is fully confounded with 10x run.** Stimulated and control cells were '
     'each captured in a single, separate 10x lane, so pseudobulk DE cannot formally '
     'separate a stimulation effect from a run-specific technical effect. The consistency '
     'of the induced ISG module across essentially all tested cell types, and its match to '
     'known IFN-β biology, makes a purely technical explanation unlikely, but this '
     'confound cannot be excluded by the data alone.\n'
     '- **CD14+ monocytes still partially split by condition** (clusters 0 and 1) even '
     'after scVI batch correction on donor×condition. This is expected — IFN stimulation '
     'is a real, strong biological state change in monocytes, not a technical batch '
     'effect, so scVI correctly did not merge it away. Both clusters were labeled '
     'identically since they are the same cell type in different activation states, but '
     'this means the "cluster" and "cell type" granularity diverge for this population.\n'
     '- **Composition test is weakly powered.** With only a handful of paired donors and '
     "compositional (sum-to-one) proportions, the paired Wilcoxon test's minimum "
     'attainable p-value limits sensitivity to modest but real proportion shifts.\n'
     '- **Rare cell types** (mast cells, erythroid cells, Tem/Trm cytotoxic T cells, pDCs) '
     'had too few cells per donor for reliable pseudobulk DE and were excluded from that '
     'analysis; their annotations also rest on smaller marker-gene evidence than the major '
     'populations.\n'
     '- **No mitochondrial genes were present** in this gene set, so the standard QC '
     'mitochondrial-fraction filter could not be used to catch stressed/dying cells; QC '
     'relied on gene-count filtering alone.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'Starting from raw counts, this analysis produced a clean, batch-corrected PBMC atlas '
     'of {{celltype.n}} annotated cell types, dominated in abundance by CD14+ monocytes, '
     'naive/memory CD4 T cells, B cells and CD8 T cells. IFN-β stimulation for 6h produced '
     'no significant change in cell-type composition but a broad, highly reproducible '
     'interferon-stimulated gene program across essentially every lineage, most pronounced '
     'in monocytes — consistent with monocytes acting as primary responders and '
     'antigen-presenting relay cells in this system, and validating the biological premise '
     'of the original Kang et al. experiment.'),
)

print(f"Done. Report: {SESSION.paths.report}")
