# Single-cell RNA-seq analysis: kang

**kang** · 24,673 → 24,558 cells · 16 clusters · 14 cell types · scVI embedding

## Overview

PBMCs from 8 SLE patients, each split into a control and IFN-β-stimulated
half, were processed from raw counts through QC, doublet removal, batch-aware integration,
clustering, cell-type annotation, and a paired stimulated-vs-control comparison. The
starting object held 24,673 cells x 15,706 genes across
8 donors and 2 conditions; condition and 10x run are
fully confounded by design (each condition was captured in its own pooled run).

## Key analysis decisions

Each row compares a standard default with what the agent chose after reading the data. Built from the tool log (`tool_calls.jsonl`), not written by the model, except the relabelling reasons, which quote the agent's tool call.

| Step | Standard default | What the data showed | Agent's choice | Effect |
|---|---|---|---|---|
| Mitochondrial cutoff | 5% | median 0%; default would remove 0 of 24,673 cells (0.0%) | 5% (default kept) | 111 cells removed by QC (0.4%) |
| Cell and gene floors | ≥200 genes/cell, genes in ≥3 cells | fewest genes in a cell: 12 | ≥200 genes/cell, genes in ≥3 cells (default kept) | 5 genes removed |
| Doublet threshold | Scrublet automatic (0.276 / 0.648 per batch) | scores not bimodal; median + 3×MAD 0.102 | **0.648** (Scrublet auto, stim run) | 4 cells removed (0.02%) |
| Normalisation | 10,000 counts/cell, 2,000 HVGs | — | 10,000 counts/cell, 2,000 HVGs (default kept) | HVGs ranked within each batch ('donor') |
| Embedding | PCA | 8 batches in `obs['donor']` | **scVI**, correcting for donor × condition | 10 latent dimensions, 16 batches |
| Leiden resolution | 1.0 | — | 1.0 (default kept) | 16 clusters |
| Annotation model | Immune_All_Low | — | Immune_All_Low (default kept) | 13 cell types |
| Cluster labels | CellTypist majority vote | Cluster 0/1: CellTypist called these 'Intermediate macrophages' but markers (CD14 log2FC~2.4-4.6 in 40-64% of cells, LYZ in 72-84%, FCN1 in 55-79%, S100A8 in 60-78%) are canonical classical monocyte markers; tissue macrophages are not expected in PBMC blood samples. Cluster 4: labeled 'NK cells' but CD3D is significantly enriched (66.4% pct-in vs 29.2% elsewhere, log2FC 2.34, padj~0) together with CD8A/CD8B, GZMB, NKG7, GNLY -- a CD3+ cytotoxic/effector CD8 T cell profile, not NK (contrast with true NK cluster 11 where CD3D is depleted: 4.5% vs 35.1% elsewhere). The coarse CellTypist model independently calls cluster 4 'T cells', agreeing with this correction. Cluster 12: labeled 'Tcm/Naive helper T cells' but PPBP (85.1% cells, log2FC 7.7), PF4 (63.6%, log2FC 7.46), and ITGA2B (7.4%, log2FC 6.76) are all massively and specifically enriched -- unambiguous platelet markers, not T cells. | cluster 0: Intermediate macrophages → **CD14+ Classical monocytes**; cluster 1: Intermediate macrophages → **CD14+ Classical monocytes**; cluster 4: NK cells → **CD8+ Effector T cells**; cluster 12: Tcm/Naive helper T cells → **Platelets** | 8,381 cells relabelled |
| Composition test | — | 8 stim / 8 ctrl samples | paired Wilcoxon signed-rank by 'donor' | 0 of 14 cell types changed (padj < 0.05) |
| DE design | ~ condition | replicates in 'donor' | **~ donor + condition** (paired) | 10 cell types tested, 6,428 significant genes in total |
| Min cells per pseudobulk sample | 10 | — | 10 (default kept) | 4 cell types skipped |

## Quality control

No mitochondrial genes were detected in the gene set (0 found despite the
standard "MT-" prefix being present), so pct-mito is uniformly 0% and provided
no filtering signal — this dataset had mitochondrial genes already stripped out upstream, so a
mito-based cutoff is a no-op rather than a QC decision here. Genes-per-cell ranged from
12 to 2,757 (median 519), consistent with a fairly
shallow 10x run typical of this dataset. The tutorial-default thresholds
(min_genes=200, max_pct_mt=5%, min_cells=3)
removed only 111 cells (0.4%) and
5 genes, so I applied them as-is — the low-gene tail was small and this
is not an aggressive cut.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (24,673 cells). The data contains no mitochondrial reads (no mitochondrial genes), so the mitochondrial filter removes nothing. Minimum genes per cell: 200, which removes 111 cells (fewest observed: 12).*

## Doublet detection

Demuxlet genotype-based demultiplexing (assigning cells to donors) had already removed
cross-donor doublets upstream; only same-donor doublets (transcriptomically indistinguishable
by genotype) can remain. The Scrublet score distribution was not bimodal (median
0.046, max 0.685), so the standard median+3·MAD rule
(0.102) was not used — with most doublets already gone, that generic rule
would treat the upper tail of real, transcriptionally-active cells as doublets. Instead, per
the documented exception for demultiplexed data, I used a light-touch cutoff: the higher of the
two per-run Scrublet automatic thresholds (0.276 / 0.648 for ctrl/stim), i.e.
0.648, applied uniformly. This removed only 4 cells
(0.02%), consistent with the expectation that little residual doublet
contamination remained.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 24,562 cells; the distribution is not bimodal. Chosen threshold 0.648 (Scrublet auto, stim run) flags 4 cells (0.0%). The mean of Scrublet's per-batch automatic thresholds (dashed, 0.462) would flag 47.*

## Dimensionality reduction

IFN-β stimulation causes a strong, cell-type-independent transcriptional shift (a large
interferon-response program), and condition is fully confounded with the 10x run. Without
correction, this shift risks splitting each cell type into two condition-specific clusters,
which would break both cluster-based annotation and the downstream composition test. I
therefore ran scVI with a combined `['donor', 'condition']` batch key (one batch per sample,
2,000 HVGs, latent dimension 10), which aligns matching cell types
across donors and across conditions in the embedding used for clustering/UMAP, while leaving
raw counts untouched for the composition and DE tests below. The risk of this choice is
over-correction: any cell state that is genuinely unique to one condition could be merged
into a resting counterpart. I checked the UMAP by condition (see figure) — clusters mix well
across stim/ctrl, supporting good integration without obvious over-merging of clearly distinct
populations.

## Clustering

Leiden clustering (resolution 1.0) on the scVI embedding produced
16 clusters ranging from 21 to 3,982 cells.

| Cluster | Cells | Top marker genes | Cell type | CellTypist label |
|---|---|---|---|---|
| 0 | 3,647 | FTL, CD63, SOD2, LGALS3, ANXA5 | **CD14+ Classical monocytes** | Intermediate macrophages |
| 1 | 2,043 | TIMP1, FTH1, TYROBP, C15orf48, FTL | **CD14+ Classical monocytes** | Intermediate macrophages |
| 2 | 1,097 | FCGR3A, MS4A7, CXCL16, LST1, S100A11 | Non-classical monocytes | Non-classical monocytes |
| 3 | 1,348 | EIF1, BTG1, UBC, UBB, H3F3B | Tcm/Naive helper T cells | Tcm/Naive helper T cells |
| 4 | 2,449 | CCL5, B2M, HLA-A, NKG7, TMSB4X | **CD8+ Effector T cells** | NK cells |
| 5 | 3,231 | PABPC1, TMSB4X, RPL10, RPS4X, RPS2 | Tem/Effector helper T cells | Tem/Effector helper T cells |
| 6 | 30 | CD3D, CCL5, GAPDH, HLA-A, RARRES3 | Tem/Trm cytotoxic T cells | Tem/Trm cytotoxic T cells |
| 7 | 1,444 | RPS6, RPL32, RPS18, RPL13, RPS14 | Tcm/Naive cytotoxic T cells | Tcm/Naive cytotoxic T cells |
| 8 | 3,982 | RPL32, RPS6, RPS14, RPL13, RPS18 | Tcm/Naive helper T cells | Tcm/Naive helper T cells |
| 9 | 115 | SEC61B, TSPAN13, TXN, CD74, HERPUD1 | pDC | pDC |
| 10 | 408 | HLA-DRA, HLA-DPB1, HLA-DPA1, HLA-DRB1, CD74 | DC2 | DC2 |
| 11 | 1,811 | GNLY, NKG7, GZMB, HLA-A, APOBEC3G | CD16+ NK cells | CD16+ NK cells |
| 12 | 242 | PPBP, PF4, GNG11, SDPR, TMSB4X | **Platelets** | Tcm/Naive helper T cells |
| 13 | 2,626 | CD74, HLA-DRA, CD79A, HLA-DRB1, HLA-DPA1 | B cells | B cells |
| 14 | 64 | HBA2, HBB, HBA1, ALAS2, SNCA | Late erythroid | Late erythroid |
| 15 | 21 | RPS4X, SERPINB1, GNB2L1, RPS2, RPS24 | Mast cells | Mast cells |

![UMAP](figures/umap.png)

*UMAP of 24,558 cells computed on the scVI latent space, coloured by cell type, Leiden cluster, batch (donor) and condition (condition). scVI corrected for donor × condition. Batches that overlap within each cell type indicate the integration worked. Conditions that overlap within each cell type mean the labels are comparable between conditions; a region with only one condition may be a condition-specific state.*

## Cell-type annotation

CellTypist (Immune_All_Low.pkl) gave fine-grained labels per cluster, cross-checked against
the coarse Immune_All_High.pkl model and canonical markers. Three corrections were made where
the fine model's calls did not hold up:

- **Clusters 0/1** were labeled "Intermediate macrophages" by the fine model, but canonical
  classical-monocyte markers CD14, LYZ, FCN1 and S100A8 were all strongly and specifically
  enriched, consistent with classical (CD14+) monocytes rather than tissue macrophages, which
  are not expected in a blood sample. Relabeled to "CD14+ Classical monocytes".
- **Cluster 4** was labeled "NK cells", but CD3D was significantly enriched alongside
  CD8A/CD8B/GZMB/NKG7, unlike the genuine CD3-negative NK cluster (cluster 11, where CD3D is
  depleted rather than enriched). The coarse model independently called this cluster
  "T cells", agreeing with the correction to "CD8+ Effector T cells".
- **Cluster 12** was labeled "Tcm/Naive helper T cells", but the platelet markers PPBP, PF4
  and ITGA2B were massively and specifically enriched, unambiguously identifying
  platelets/platelet-containing droplets rather than T cells. Relabeled to "Platelets".

All other clusters agreed between the two CellTypist models and canonical markers: non-classical
monocytes (FCGR3A, MS4A7 enriched), naive/central-memory CD4 T cells (CCR7, IL7R, SELL
enriched), CD8 Tcm/naive (CD8B, CCR7, SELL enriched), CD8 Tem/Trm (GZMK, TIGIT, CD3D
enriched), B cells (CD79A, MS4A1, CD19 enriched), plasmacytoid DC (LILRA4, GZMB enriched),
DC2 (CD1C, FCER1A, LYZ enriched), erythroid cells (HBB, HBA1 near-universally expressed and
strongly enriched), and a small, rare cluster with GATA2 strongly and specifically enriched
(consistent with mast cell/basophil identity; CD3D and lymphoid/myeloid markers essentially
absent), retained as "Mast cells".

Final composition: 5,690 cells (23.2%) classical monocytes,
5,330 cells (21.7%) naive/Tcm CD4 T cells,
3,231 cells (13.2%) effector/Tem CD4 T cells,
2,626 cells (10.7%) B cells, 2,449 cells (10.0%) CD8 effector T cells,
1,811 cells (7.4%) NK cells, 1,444 cells (5.9%) naive/Tcm CD8 T
cells, 1,097 cells (4.5%) non-classical monocytes, 408 cells (1.7%) DC2,
242 cells (1.0%) platelets, 115 cells (0.5%) pDC, 64 cells (0.3%) erythroid
cells, 30 cells (0.1%) Tem/Trm CD8 T cells, and
21 cells (0.1%) mast cells, out of 24,558 total cells after QC/doublet
filtering.

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (28 of 39 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: CD4, NCAM1, CD19, CD1C, CLEC9A, FCER1A, ITGA2B, TPSAB1, KIT, HDC, SLC18A2.*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Composition analysis

Cell-type proportions were compared stim vs ctrl per donor (paired paired Wilcoxon signed-rank, all
8 donors contributing both conditions).

| Cell type | Cells | Share |
|---|---|---|
| CD14+ Classical monocytes | 5,690 | 23.2% |
| Tcm/Naive helper T cells | 5,330 | 21.7% |
| Tem/Effector helper T cells | 3,231 | 13.2% |
| B cells | 2,626 | 10.7% |
| CD8+ Effector T cells | 2,449 | 10.0% |
| CD16+ NK cells | 1,811 | 7.4% |
| Tcm/Naive cytotoxic T cells | 1,444 | 5.9% |
| Non-classical monocytes | 1,097 | 4.5% |
| DC2 | 408 | 1.7% |
| Platelets | 242 | 1.0% |
| pDC | 115 | 0.5% |
| Late erythroid | 64 | 0.3% |
| Tem/Trm cytotoxic T cells | 30 | 0.1% |
| Mast cells | 21 | 0.1% |
| **Total** | **24,558** | |

| Cell type | Mean share, ctrl | Mean share, stim | log2 ratio | padj |
|---|---|---|---|---|
| DC2 | 2.2% | 1.5% | -0.52 | 0.33 |
| CD14+ Classical monocytes | 25.4% | 24.4% | -0.06 | 0.68 |
| CD8+ Effector T cells | 9.9% | 9.4% | -0.08 | 0.68 |
| pDC | 0.3% | 0.8% | 1.27 | 0.68 |
| B cells | 9.8% | 10.1% | 0.05 | 0.91 |
| CD16+ NK cells | 7.1% | 7.4% | 0.08 | 0.91 |
| Late erythroid | 0.2% | 0.3% | 0.26 | 0.91 |
| Non-classical monocytes | 4.7% | 4.6% | -0.03 | 0.91 |
| Platelets | 0.9% | 1.1% | 0.23 | 0.91 |
| Tcm/Naive cytotoxic T cells | 5.4% | 5.6% | 0.07 | 0.91 |
| Tcm/Naive helper T cells | 19.9% | 20.5% | 0.04 | 0.91 |
| Tem/Effector helper T cells | 14.1% | 14.1% | 0.00 | 0.91 |
| Tem/Trm cytotoxic T cells | 0.1% | 0.1% | 0.45 | 0.91 |
| Mast cells | 0.1% | 0.1% | -0.11 | 0.95 |

No cell type showed a significant proportion shift after multiple-testing correction
(0 of 14 significant). The largest nominal changes were
a decrease in DC2 (2.2% -> 1.5%, padj 0.33) and an
increase in pDC (0.3% -> 0.8%, padj 0.68), but neither
survived correction and pDC could only be tested with very few cells per sample. With so few
donor replicates, the paired test's attainable significance is inherently limited, so this
analysis is under-powered to detect anything but large, consistent shifts; a true absence of
compositional change over such a short stimulation window is also biologically plausible,
since PBMC subset frequencies are not expected to shift substantially without cell
death/proliferation over a brief ex vivo culture.

![Cell-type composition](figures/composition.png)

*Cells per annotated type (24,558 cells, 14 types; CellTypist majority vote over Leiden clusters, with clusters relabelled by the agent from their markers (see the decisions table)).*

![Cell-type proportions by condition](figures/composition_change.png)

*Share of each cell type per sample, ctrl (blue) vs stim (orange); lines join each donor's two samples. Wilcoxon signed-rank test with Benjamini-Hochberg correction: 0 of 14 cell types differ at padj < 0.05. Log scale above 0.1%, linear below, so samples with none of a type sit at 0.*

## Differential expression

Pseudobulk DE (PyDESeq2, design `~donor + condition`, paired by donor) was run per cell type
using raw counts summed per donor.

| Cell type | Samples (test / ref) | Genes tested | Up | Down | Top up-regulated genes |
|---|---|---|---|---|---|
| B cells | 8 / 8 | 2,101 | 347 | 273 | ISG20, B2M, ISG15, LY6E, UBE2L6 |
| CD14+ Classical monocytes | 8 / 8 | 3,780 | 1,221 | 1,198 | SSB, CCL8, IL1RN, DYNLT1, ISG20 |
| CD16+ NK cells | 8 / 8 | 1,450 | 220 | 169 | ISG15, ISG20, TNFSF10, IFI6, IFIT1 |
| CD8+ Effector T cells | 8 / 8 | 1,248 | 165 | 116 | ISG20, ISG15, IFI6, SAT1, IFIT1 |
| DC2 | 7 / 8 | 1,085 | 296 | 229 | ISG15, ISG20, HSP90AA1, LY6E, MX1 |
| Non-classical monocytes | 8 / 8 | 1,612 | 425 | 392 | APOBEC3A, TNFSF10, MYL12A, ISG20, ISG15 |
| Platelets | 5 / 4 | 339 | 30 | 3 | ISG15, ISG20, IFI6, B2M, LY6E |
| Tcm/Naive cytotoxic T cells | 7 / 6 | 1,214 | 130 | 117 | ISG15, ISG20, LY6E, MX1, IFIT1 |
| Tcm/Naive helper T cells | 8 / 8 | 3,144 | 369 | 293 | ISG20, PSMB9, IFI16, LY6E, B2M |
| Tem/Effector helper T cells | 8 / 8 | 2,340 | 283 | 152 | ISG20, SAT1, MT2A, TMSB10, TNFSF10 |

Every cell type tested showed a clear, coherent interferon-stimulated gene (ISG) signature,
as expected for IFN-β treatment: ISG15 (log2FC 7.66, padj 1.1e-74) in classical
monocytes, ISG15 (log2FC 4.69, padj 6.3e-175) in NK cells, ISG15 (log2FC 5.54, padj 9.4e-140) in B cells, and
ISG15 (log2FC 4.77, padj 1e-152) in CD8 effector T cells, together with consistent
up-regulation of ISG20 (log2FC 5.72, padj 1.1e-110), IFIT1 (log2FC 8.65, padj 4.4e-59),
and MX1 (log2FC 6.30, padj 3.4e-49) across essentially all tested populations. Classical
monocytes showed the largest and most cell-type-specific response, including strong induction
of the chemokine CCL8 (log2FC 10.00, padj 5.7e-127) and the anti-inflammatory cytokine
antagonist IL1RN (log2FC 6.66, padj 3.4e-116), both far less induced in non-classical
monocytes (CCL8 (log2FC 8.23, padj 3.4e-09), IL1RN (log2FC 3.94, padj 4.7e-15)) and
essentially absent from lymphoid populations tested for these genes — consistent with monocytes
being a primary sensor/responder to IFN-β and secondary inflammatory-mediator producer. Rare
populations (64 cells (0.3%) erythroid cells, 21 cells (0.1%) mast cells,
30 cells (0.1%) Tem/Trm cytotoxic T cells, 115 cells (0.5%) pDC cells)
had too few cells per donor in one or both conditions and were skipped, as noted by the tool
(Late erythroid, Mast cells, Tem/Trm cytotoxic T cells, pDC).

![Volcano plots per cell type](figures/de_volcano.png)

*Pseudobulk differential expression, stim vs ctrl, per cell type (PyDESeq2, design `~donor + condition`, Wald test). Red: up in stim; blue: down (padj < 0.05); grey: not significant. The three most significant up-regulated genes are labelled. Not tested (too few samples): Late erythroid, Mast cells, Tem/Trm cytotoxic T cells, pDC.*

## Caveats

- **Condition/run confound**: stim and ctrl were captured in separate 10x runs, so pseudobulk
  DE cannot fully separate a true IFN-β effect from a run-specific technical effect; the very
  large, coherent, and biologically expected ISG signature makes a purely technical explanation
  unlikely, but this cannot be ruled out formally.
- **scVI batch correction with `donor + condition`** aligns cell types across conditions for
  clustering/annotation but could, in principle, over-merge a condition-specific cell state
  into a resting counterpart; the reasonably good mixing of conditions on the UMAP argues
  against gross over-merging, but subtle state-specific subpopulations could still be masked at
  this clustering resolution.
- **Composition test is under-powered** with only 8 paired donors, so the null
  result for composition changes should not be read as strong evidence of no change, only as no
  detected change at this sample size.
- **Rare clusters** (21 cells (0.1%) mast cells, 115 cells (0.5%) pDC,
  242 cells (1.0%) platelets, 30 cells (0.1%) Tem/Trm cytotoxic T
  cells) have limited statistical power for both marker validation and DE; labels here rely on
  strong marker fold-changes despite low cell counts.
- **No mitochondrial genes were present** in this dataset, so QC could not filter on
  mito-fraction; only gene-count-based filtering was applied.

## Conclusions

After QC, light-touch doublet removal (justified by prior genotype-based demultiplexing),
and scVI-based integration across donors and conditions, clustering recovered the expected
major PBMC lineages (classical and non-classical monocytes, CD4 and CD8 T cell subsets, NK
cells, B cells, DC2, pDC, platelets, mast cells, and a small erythroid contaminant), after
correcting three CellTypist mislabels using canonical marker evidence. IFN-β stimulation did
not detectably shift PBMC subset proportions but induced a strong, broad, biologically
coherent interferon-response transcriptional program across essentially every cell type, most
pronounced in classical monocytes, which additionally mounted a distinct secondary
inflammatory-mediator response (CCL8, IL1RN) not seen in most other populations.

## Methods

*Generated from the tool log: these are the steps and parameters that actually ran.*

**Quality control.** Per-cell metrics were computed with scanpy `calculate_qc_metrics`; mitochondrial genes were those prefixed `MT-` (0 found).

**Filtering.** Cells with fewer than 200 detected genes or more than 5% mitochondrial reads were removed, as were genes detected in fewer than 3 cells (24,673 → 24,562 cells, 15,706 → 15,701 genes).

**Doublets.** Doublet scores were computed with Scrublet (scanpy `pp.scrublet`), run separately within each batch ('condition'), on raw counts; cells scoring ≥ 0.648 were removed (4 cells, 0.02%).

**Normalisation.** Raw counts were kept in `layers['counts']`; expression was scaled to 10,000 counts per cell and log1p-transformed. The top 2,000 highly variable genes were flagged, ranked within each batch ('donor').

**Embedding.** An scVI model (10 latent dimensions, 326 epochs) was trained on raw counts of 2,000 highly variable genes, correcting for donor × condition; its latent space was used for clustering, UMAP and annotation (raw counts were unchanged).

**Clustering.** A k-nearest-neighbour graph on that embedding was clustered with Leiden (resolution 1.0; 16 clusters) and embedded with UMAP.

**Markers and annotation.** Marker genes per cluster were ranked with a Wilcoxon rank-sum test. Cell types were assigned with CellTypist (model `Immune_All_Low.pkl`), using majority voting over the Leiden clusters.

**Label curation.** Where marker genes contradicted CellTypist, the agent relabelled whole clusters (4 clusters; listed in the decisions table with the evidence). The original CellTypist labels are kept in `obs['cell_type_celltypist']`.

**Composition.** For each sample ('donor' × condition), the share of cells of each type ('cell_type') was computed and compared between stim and ctrl with a paired Wilcoxon signed-rank test; p-values were adjusted across cell types with the Benjamini-Hochberg method.

**Differential expression.** Raw counts were summed per sample ('donor' × condition) within each cell type ('cell_type'); samples with fewer than 10 cells were dropped, and genes were kept with at least 10 counts in at least as many samples as the smaller condition. PyDESeq2 was fitted with design `~donor + condition` and stim was compared with ctrl using a Wald test; p-values were adjusted with Benjamini-Hochberg (padj < 0.05 called significant). Log2 fold changes are unshrunken.

### Software

| Software | Version |
|---|---|
| Python | 3.11.14 |
| scanpy | 1.11.5 |
| anndata | 0.12.19 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| scvi-tools | 1.4.2 |
| celltypist | 1.7.1 |
| statsmodels | 0.14.6 |
| pydeseq2 | 0.5.4 |
| LLM (analysis decisions and narrative) | claude-sonnet-5 |

### Reproducibility

Random seed 0 for all stochastic steps. Every tool call, with the arguments the agent chose and the summary it read back, is in `tool_calls.jsonl`. `replay.py` re-runs those calls without the model, reproducing this analysis; running the agent again samples new choices from the model, so it can take different decisions.

---

*This report was produced by an AI agent (claude-sonnet-5). The run summary, decisions table, figures, captions and Methods are generated by code from the tool log; the narrative is written by the model and should be checked against them before use.*
