"""Replay of the agent's tool calls for pbmc_multibatch, without the LLM.

Generated: 2026-09-28T14:36:06
Git commit: 24068e592d14aa176f5b6e5c079649c4ed21485a
Source log: tool_calls.jsonl (19 of 20 calls; failed calls omitted)

Re-runs the same tool functions with the same arguments and seed. Output goes to
outputs/pbmc_multibatch_replay/, with its own tool log and report. Package versions are pinned by
uv.lock at the commit above: check it out first if the code has changed since.

Caveats: annotate_celltypes downloads CellTypist models, and scVI training can differ
slightly between hardware (CPU vs GPU), which can change clusters downstream.

  uv run python examples/pbmc_multibatch/replay.py
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
SESSION.load('data/pbmc_multibatch.h5ad')
SESSION.paths = config.RunPaths('pbmc_multibatch_replay')
SESSION.begin_run()
print(f"Replay output: {SESSION.paths.dir}")

_run("inspect_dataset")

_run("check_gene_identifiers")

_run("compute_qc")

_run("recommend_qc_thresholds")

_run(
    "filter_cells_and_genes",
    min_genes=200,
    max_pct_mt=15,
    min_cells=3,
)

_run(
    "detect_doublets",
    batch_key='batch',
)

_run(
    "filter_doublets",
    threshold=0.09,
)

_run("normalize")

_run(
    "run_scvi",
    batch_key='batch',
)

_run(
    "cluster",
    resolution=1.0,
)

_run("identify_markers")

_run("annotate_celltypes")

_run(
    "check_markers",
    genes=['CD3D', 'CD3E', 'IL7R', 'CD4', 'CD8A', 'CD8B', 'CCR7', 'SELL', 'CD14', 'FCGR3A', 'LYZ',
     'S100A8', 'MS4A1', 'NKG7', 'MT-CO1', 'MALAT1'],
    clusters=['0', '1', '4', '5', '7', '8'],
)

_run(
    "check_markers",
    genes=['GNLY', 'NKG7', 'GZMB', 'CD79A', 'MS4A1', 'FCER1A', 'CD1C', 'GZMK', 'CCL5', 'KLRB1',
     'SLC4A10', 'FCGR3A', 'MS4A7', 'LILRA4', 'IL3RA', 'CD14'],
    clusters=['2', '3', '6', '9', '10', '11', '12'],
)

_run(
    "relabel_clusters",
    labels={'8': 'Low-quality/stressed cells'},
    reason=("Cluster 8 (24 cells) was labeled 'Tcm/Naive helper T cells' by CellTypist, but "
     'canonical T-cell markers are not enriched: CD3D (padj=1.0, 45.8% in vs 45.7% '
     'elsewhere), CD3E (padj=1.0), IL7R (padj=1.0), CCR7 (padj=1.0), and SELL is actually '
     'depleted (log2FC=-2.44, padj=0.016). Instead its top-ranked genes are MALAT1, JUN, '
     'and a long list of mitochondrial transcripts (MT-CO1, MT-ATP6, MT-ND5, MT-CYB, '
     'MT-ND4, MT-CO3, MT-CO2, MT-ND1) plus MTRNR2L12, a classic signature of '
     'stressed/damaged cells rather than a genuine lineage. No lineage marker (CD14, CD3, '
     'CD79A/MS4A1, GNLY) is enriched here.'),
)

_run("identify_markers")

_run(
    "check_markers",
    genes=['CD14', 'LYZ', 'S100A8', 'FCGR3A', 'MS4A7', 'CD3D', 'CD3E', 'IL7R', 'CCR7', 'SELL',
     'CD8A', 'CD8B', 'GZMK', 'CCL5', 'KLRB1', 'SLC4A10', 'GNLY', 'NKG7', 'GZMB', 'CD79A',
     'MS4A1', 'FCER1A', 'CD1C', 'IL3RA', 'LILRA4'],
    clusters=['0', '1', '2', '3', '4', '5', '6', '7', '9', '10', '11', '12'],
)

_run("summarize_findings")

_run(
    "generate_report",
    report_markdown=('## Overview\n'
     '\n'
     'This dataset combines PBMCs from two separate 10x runs (PBMC5k and PBMC10k) into one '
     'object of {{input.n_cells}} cells x {{input.n_genes}} genes, provided as raw counts '
     'with human gene symbols. Because the two runs are a genuine technical batch factor '
     '(different chip lanes/reagent lots run at different times), batch identity was '
     'treated as a covariate to correct for during integration rather than ignored.\n'
     '\n'
     '## Quality control\n'
     '\n'
     'Mitochondrial genes were identified via the `MT-` prefix ({{qc.n_mito_genes}} '
     'genes). The per-cell QC distributions showed a notably right-shifted mitochondrial '
     'fraction for this dataset (median {{qc.mito.median}}, p95 {{qc.mito.p95}}, p99 '
     '{{qc.mito.p99}}, max {{qc.mito.max}}) compared to the tutorial-default cutoff of '
     '{{default.max_pct_mt}}. Applying that default would have discarded '
     '{{default.mito_cells_removed}} cells ({{default.mito_pct_removed}} of the dataset) — '
     'clearly disproportionate, and consistent with a global shift in mito content rather '
     'than a distinct population of dying cells.\n'
     '\n'
     'Instead, a max mitochondrial threshold of {{filter.max_pct_mt}} was used, trimming '
     'only the extreme tail while retaining the bulk of the distribution. Combined with '
     'the standard minimum-genes-per-cell floor ({{filter.min_genes}}) — which removed no '
     'additional cells since the lowest cell already had {{qc.genes.min}} genes detected — '
     'and a minimum-cells-per-gene filter ({{filter.min_cells}}), this removed '
     '{{filter.cells_removed}} cells ({{filter.pct_cells_removed}}) and '
     '{{filter.genes_removed}} lowly-detected genes, leaving {{filter.cells_after}} cells '
     'x {{filter.genes_after}} genes.\n'
     '\n'
     '## Doublet detection\n'
     '\n'
     'Scrublet was run separately per 10x run (batch key `batch`) since doublet rates and '
     'score scales are run-specific. The combined score distribution was not bimodal, so '
     'the standard median + 3xMAD rule was used rather than a bimodal valley, giving a '
     'threshold of {{doublet.threshold}} (median {{doublet.score_median}}, max observed '
     '{{doublet.score_max}}). This flagged and removed {{doublet.cells_removed}} cells '
     '({{doublet.pct_removed}}), leaving {{doublet.cells_after}} cells. There was no '
     'indication in the task description that doublets had already been removed upstream '
     '(e.g. by hashing or genotype demultiplexing), so the standard rule — rather than the '
     'lighter-touch Scrublet-automatic threshold — was appropriate here.\n'
     '\n'
     '## Dimensionality reduction\n'
     '\n'
     'The dataset has {{input.n_batches}} batches from two independent 10x runs, a '
     'technical grouping rather than a biological condition. Since a real batch effect '
     'between separate sequencing runs is expected (loading, capture efficiency, ambient '
     'RNA differences), scVI was used with `batch` as the correction key to learn a joint '
     'latent embedding ({{scvi.n_latent}} latent dimensions, {{scvi.epochs}} training '
     'epochs on {{scvi.n_hvgs}} HVGs) rather than plain PCA, which would leave run-driven '
     'variation uncorrected in the neighbor graph and could split shared cell types by '
     'batch. Raw counts were preserved throughout for downstream marker testing.\n'
     '\n'
     '## Clustering\n'
     '\n'
     'Leiden clustering on the scVI latent space (resolution {{cluster.resolution}}) '
     'produced {{cluster.n}} clusters ranging from {{cluster.smallest}} to '
     '{{cluster.largest}} cells. {{table:clusters}}\n'
     '\n'
     '## Cell-type annotation\n'
     '\n'
     'CellTypist (Immune_All_Low.pkl, majority vote per cluster) provided initial labels, '
     'cross-checked against the coarser Immune_All_High.pkl model and against canonical '
     'marker genes for every final cell type, using check_markers on the full DE ranking.\n'
     '\n'
     '- **Classical monocytes** (clusters 0 and 1): CD14 and LYZ strongly and specifically '
     'enriched, with S100A8/S100A9 additionally enriched in cluster 0, and no FCGR3A/MS4A7 '
     'enrichment ruling out non-classical identity.\n'
     '- **Non-classical monocytes**: FCGR3A and MS4A7 both strongly enriched, CD14 not '
     'distinctly enriched relative to other clusters, consistent with the CD14dim/CD16+ '
     'non-classical subset.\n'
     '- **CD16+ NK cells**: GNLY, NKG7 and GZMB all highly and near-uniformly enriched, '
     'FCGR3A co-enriched, with CD3D/CD3E largely absent, ruling out a T/NK doublet '
     'population.\n'
     '- **Naive B cells**: CD79A and MS4A1 both sharply enriched and essentially absent '
     'elsewhere.\n'
     '- **DC2**: FCER1A and CD1C sharply and specifically enriched, distinguishing this '
     'small cluster from monocytes and pDCs.\n'
     '- **pDC**: IL3RA, LILRA4 and GZMB all enriched to near-ubiquitous expression in this '
     'small cluster, a canonical pDC signature.\n'
     '- **Tcm/Naive helper T cells** (cluster 4): CD3D/CD3E/IL7R enriched together with '
     'high CCR7 and SELL, and CD8A/CD8B essentially absent — a CD4 naive/central-memory '
     'phenotype (CD4 mRNA itself is not a reliable discriminator, as it is also expressed '
     'by monocytes).\n'
     '- **Tcm/Naive cytotoxic T cells** (cluster 5): CD3D/CD3E enriched together with '
     'CD8A/CD8B and elevated CCR7/SELL — CD8 naive/central memory.\n'
     '- **Tem/Effector helper T cells** (cluster 7): CD3D/CD3E/IL7R strongly enriched with '
     'low CCR7/SELL, and CD8A/CD8B absent — an effector/memory CD4 phenotype.\n'
     '- **Tem/Trm cytotoxic T cells** (cluster 9): CCL5, GZMA/GZMK and CD8A/CD8B enriched '
     'together with CD3D/CD3E — a cytotoxic CD8 effector/memory phenotype.\n'
     '- **MAIT cells**: KLRB1, GZMK and the MAIT-specific SLC4A10 all enriched, together '
     'with CD3E and CCL5.\n'
     '- **Low-quality/stressed cells** (cluster 8, {{celltype:Low-quality/stressed cells}} '
     'cells): CellTypist called this cluster "Tcm/Naive helper T cells", but none of CD3D, '
     'CD3E, IL7R or CCR7 were significantly enriched, and SELL was actually depleted '
     'relative to the rest of the dataset. Its defining genes were instead MALAT1, JUN and '
     'a long run of mitochondrial transcripts (MT-CO1, MT-ATP6, MT-ND5, MT-CYB, MT-ND4, '
     'MT-CO3, MT-CO2, MT-ND1) plus MTRNR2L12 — a stress/damage signature with no lineage '
     'marker enriched. It was relabeled accordingly rather than reported as a T-cell '
     'subset.\n'
     '\n'
     'Final cell-type composition: {{table:composition}}\n'
     '\n'
     '## Caveats\n'
     '\n'
     '- Cluster 8 ({{celltype:Low-quality/stressed cells}} cells) reflects '
     'damaged/stressed cells rather than a discrete biological population and should be '
     'excluded from downstream biological interpretation; it is retained in the annotated '
     'object for transparency.\n'
     '- The two clusters both labeled "Classical monocytes" (clusters 0 and 1) differ '
     'somewhat in markers (cluster 1 shows higher HLA-DR/CST3/CD68), which may reflect an '
     'activation or maturity gradient within classical monocytes rather than a distinct '
     'subset — worth revisiting with finer-resolution clustering if a monocyte-focused '
     'question arises.\n'
     '- No biological condition was present in this dataset (only the two technical 10x '
     'runs), so no composition or differential-expression comparison between conditions '
     'was performed; this analysis focused on quality control, batch-aware integration and '
     'cell-type annotation only.\n'
     '- The elevated mitochondrial fraction across this dataset (well above the '
     'tutorial-default cutoff of {{default.max_pct_mt}}) means a permissive filtering '
     'threshold was used; cells at the upper end of the retained range should be '
     'interpreted with some caution as possibly reflecting somewhat lower RNA quality '
     'rather than a wholly clean cell population.\n'
     '\n'
     '## Conclusions\n'
     '\n'
     'Starting from {{input.n_cells}} raw cells across two PBMC 10x runs, QC and doublet '
     "filtering tuned to this dataset's actual distributions (rather than default tutorial "
     'cutoffs) retained {{final.n_cells}} high-quality, singlet cells. Batch-aware '
     'integration with scVI resolved the two runs into a shared embedding, yielding '
     '{{cluster.n}} Leiden clusters that mapped onto {{celltype.n}} expected PBMC cell '
     'types (monocyte subsets, B cells, NK cells, CD4/CD8 T cell subsets including MAIT '
     'cells, DC2, and pDCs), each supported by canonical marker evidence, with one small '
     'artefactual cluster of stressed cells correctly separated out from genuine '
     'lineages.\n'),
)

print(f"Done. Report: {SESSION.paths.report}")
