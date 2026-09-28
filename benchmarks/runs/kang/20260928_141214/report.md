# Single-cell RNA-seq analysis: kang

**kang** · 24,673 → 24,558 cells · 16 clusters · 14 cell types · scVI embedding

## Overview

This analysis processed PBMCs from 8 SLE patients, each split into interferon-beta–stimulated and control halves, captured as two separate 10x runs (`condition` = ctrl/stim, one run per condition, 8 donors pooled per run). Starting from 24,673 cells x 15,706 genes, the workflow proceeded through QC filtering, light-touch doublet removal, scVI-based batch-corrected clustering, CellTypist annotation with marker-based correction, and paired composition/differential-expression testing between stimulated and control samples across donors.

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
| Cluster labels | CellTypist majority vote | Cluster 0/1: CellTypist called these 'Intermediate macrophages', but macrophages are tissue-resident and not expected in PBMC blood samples; markers show classical monocyte identity (CD14 39.6%/64.0%, FCN1 78.9%/54.3%, LYZ 84.1%/72.3%, S100A8 60.6%/77.5%, all highly significant). Cluster 4: labeled 'NK cells' but strongly expresses CD3D (66.4%, padj~0), CD3E, CD8A (40.1%) and CD8B (30.7%), i.e. a bona fide T-cell receptor complex, unlike the true NK cluster 11 (CD3D only 4.5%); this is a CD3+CD8+ cytotoxic/effector T-cell population despite high GZMB/NKG7/GNLY. Cluster 12: labeled a T-helper subset but its top markers are PPBP (85.1%), PF4 (63.6%), GNG11 (61.6%), TUBB1 (36%) - an unambiguous platelet/megakaryocyte signature. | cluster 0: Intermediate macrophages → **CD14+ Monocytes**; cluster 1: Intermediate macrophages → **CD14+ Monocytes**; cluster 4: NK cells → **CD8+ Effector T cells**; cluster 12: Tcm/Naive helper T cells → **Platelets** | 8,381 cells relabelled |
| Composition test | — | 8 stim / 8 ctrl samples | paired Wilcoxon signed-rank by 'donor' | 0 of 14 cell types changed (padj < 0.05) |
| DE design | ~ condition | replicates in 'donor' | **~ donor + condition** (paired) | 10 cell types tested, 6,428 significant genes in total |
| Min cells per pseudobulk sample | 10 | — | 10 (default kept) | 4 cell types skipped |

## Quality control

Gene-identifier detection found standard human gene symbols, but no mitochondrial genes (`MT-` prefix) were present in the object at all (0 found), so pct-mitochondrial is 0% for every cell — mitochondrial genes were evidently removed upstream, before this dataset was provided. This makes the mitochondrial filter a no-op safeguard rather than a real QC lever here.

Genes-per-cell ranged from a minimum of 12 to a maximum of 2,757 (median 519), and total counts per cell from 562 to 11,816 (median 1,246). The recommended tutorial-default thresholds (min_genes=200, max_pct_mt=5%, min_cells=3) were adopted as-is: they removed only 111 cells (0.4%) and 5 genes, a light and appropriate trim for an already reasonably clean dataset.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (24,673 cells). The data contains no mitochondrial reads (no mitochondrial genes), so the mitochondrial filter removes nothing. Minimum genes per cell: 200, which removes 111 cells (fewest observed: 12).*

## Doublet detection

Demuxlet genotype-based demultiplexing was used upstream to assign cells to donors, which removes doublets formed by two cells from different patients — but same-patient doublets survive undetected by genotype. Scrublet was run separately per 10x run (`condition`, since stim and ctrl are the two distinct captures), giving per-batch automatic thresholds of 0.276 / 0.648 (ctrl / stim). The score distribution was not bimodal (not bimodal), and because most doublets were already removed upstream, the standard median+3*MAD rule (0.102) would flag a much larger fraction of cells and over-aggressively cut real cells from the upper tail of this already-cleaned distribution. Following the light-touch guidance for pre-demultiplexed data, the higher of the two per-batch Scrublet automatic thresholds (0.648) was applied dataset-wide as a conservative cutoff, removing only 4 cells (0.02%) — consistent with the expectation that few doublets remained.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 24,562 cells; the distribution is not bimodal. Chosen threshold 0.648 (Scrublet auto, stim run) flags 4 cells (0.0%). The mean of Scrublet's per-batch automatic thresholds (dashed, 0.462) would flag 47.*

## Dimensionality reduction

Condition and 10x run are the same technical factor in this dataset, and IFN-β stimulation is known to induce a strong, broad transcriptional shift (the interferon-stimulated gene program) that risks splitting each cell type into separate stim/ctrl clusters if left uncorrected. To align shared cell types across conditions for consistent annotation while leaving raw counts untouched for the downstream DE and composition tests, scVI was trained with a combined batch key of donor and condition (one batch per sample, 10-dimensional latent space, 326 epochs) rather than plain PCA. The resulting embedding was used for neighbor graph, clustering and UMAP.

## Clustering

Leiden clustering at resolution 1.0 on the scVI embedding yielded 16 clusters, ranging from 21 to 3,982 cells. 

| Cluster | Cells | Top marker genes | Cell type | CellTypist label |
|---|---|---|---|---|
| 0 | 3,647 | FTL, CD63, SOD2, LGALS3, ANXA5 | **CD14+ Monocytes** | Intermediate macrophages |
| 1 | 2,043 | TIMP1, FTH1, TYROBP, C15orf48, FTL | **CD14+ Monocytes** | Intermediate macrophages |
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

CellTypist (Immune_All_Low.pkl, majority vote per cluster) provided initial fine-grained labels, cross-checked against the coarser Immune_All_High.pkl model and canonical markers for every final cell type.

Three systematic issues were corrected:

- **Clusters 0 and 1** were called "Intermediate macrophages" by both CellTypist models, but macrophages are tissue-resident and not expected in blood. Markers instead showed a clear classical-monocyte signature (CD14, FCN1, LYZ, S100A8 all strongly enriched), so these were relabeled to CD14+ Monocytes, giving 5,690 cells (23.2%) in this type.
- **Cluster 4** was labeled "NK cells," but it strongly expressed the T-cell receptor complex genes CD3D and CD3E along with CD8A/CD8B, unlike the true NK cluster (cluster 11), which lacks CD3D. It was relabeled to CD8+ Effector T cells (retaining high GZMB/NKG7 as an activated/cytotoxic phenotype), giving 2,449 cells (10.0%).
- **Cluster 12** was labeled a naive T-helper subset, but its dominant markers (PPBP, PF4, GNG11, TUBB1) are an unambiguous platelet/megakaryocyte signature; relabeled to Platelets, giving 242 cells (1.0%).

All other clusters agreed with CellTypist and canonical markers: Non-classical monocytes (FCGR3A, MS4A7), B cells (CD79A, MS4A1), CD16+ NK cells (GNLY, NKG7, KLRD1, FGFBP2, CD3D-negative), pDC (IL3RA, TCF4, GZMB), DC2 (LYZ, HLA-DR), naive/memory helper and cytotoxic T-cell subsets (distinguished by CCR7/SELL/IL7R vs CD8A), and Late erythroid (HBB detected in essentially every cell of the cluster). Cluster 15 ("Mast cells", only 21 cells (0.1%)) showed a significant GATA2 enrichment and a TPSAB1/FCER1A trend, while CD3D/CD3E were undetectable, ruling out the alternative T-cell label from the second model — this identification should be treated cautiously given the very small cluster size.

Final cell-type composition: 

| Cell type | Cells | Share |
|---|---|---|
| CD14+ Monocytes | 5,690 | 23.2% |
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

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (29 of 35 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: CD163, CD4, TPSAB1, KIT, FOXP3, FCER1A.*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Composition analysis

Cell-type proportions were compared between stim and ctrl per donor (paired Wilcoxon signed-rank test, paired Wilcoxon signed-rank, 8 donors per condition, Benjamini-Hochberg corrected). 0 of 14 cell types reached significance after correction. 

| Cell type | Mean share, ctrl | Mean share, stim | log2 ratio | padj |
|---|---|---|---|---|
| DC2 | 2.2% | 1.5% | -0.52 | 0.33 |
| CD14+ Monocytes | 25.4% | 24.4% | -0.06 | 0.68 |
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

This null result is consistent with the short stimulation window used in this experiment: IFN-β reprograms gene expression rapidly but is not expected to alter the balance of major PBMC lineages so quickly, and the compositional nature of these proportions (a shift in one type mechanically moves the others) combined with only 8 paired donors limits power to detect small shifts.

![Cell-type composition](figures/composition.png)

*Cells per annotated type (24,558 cells, 14 types; CellTypist majority vote over Leiden clusters, with clusters relabelled by the agent from their markers (see the decisions table)).*

![Cell-type proportions by condition](figures/composition_change.png)

*Share of each cell type per sample, ctrl (blue) vs stim (orange); lines join each donor's two samples. Wilcoxon signed-rank test with Benjamini-Hochberg correction: 0 of 14 cell types differ at padj < 0.05. Log scale above 0.1%, linear below, so samples with none of a type sit at 0.*

## Differential expression

Pseudobulk DE (PyDESeq2, design `~donor + condition`, donor included as a paired covariate) was run per cell type on raw counts. Late erythroid, Mast cells, Tem/Trm cytotoxic T cells, pDC were skipped for having too few samples with sufficient cells per condition. 

| Cell type | Samples (test / ref) | Genes tested | Up | Down | Top up-regulated genes |
|---|---|---|---|---|---|
| B cells | 8 / 8 | 2,101 | 347 | 273 | ISG20, B2M, ISG15, LY6E, UBE2L6 |
| CD14+ Monocytes | 8 / 8 | 3,780 | 1,221 | 1,198 | SSB, CCL8, IL1RN, DYNLT1, ISG20 |
| CD16+ NK cells | 8 / 8 | 1,450 | 220 | 169 | ISG15, ISG20, TNFSF10, IFI6, IFIT1 |
| CD8+ Effector T cells | 8 / 8 | 1,248 | 165 | 116 | ISG20, ISG15, IFI6, SAT1, IFIT1 |
| DC2 | 7 / 8 | 1,085 | 296 | 229 | ISG15, ISG20, HSP90AA1, LY6E, MX1 |
| Non-classical monocytes | 8 / 8 | 1,612 | 425 | 392 | APOBEC3A, TNFSF10, MYL12A, ISG20, ISG15 |
| Platelets | 5 / 4 | 339 | 30 | 3 | ISG15, ISG20, IFI6, B2M, LY6E |
| Tcm/Naive cytotoxic T cells | 7 / 6 | 1,214 | 130 | 117 | ISG15, ISG20, LY6E, MX1, IFIT1 |
| Tcm/Naive helper T cells | 8 / 8 | 3,144 | 369 | 293 | ISG20, PSMB9, IFI16, LY6E, B2M |
| Tem/Effector helper T cells | 8 / 8 | 2,340 | 283 | 152 | ISG20, SAT1, MT2A, TMSB10, TNFSF10 |

Across every tested cell type, the top upregulated genes in stim vs ctrl were canonical type-I interferon-stimulated genes: ISG15 (log2FC 7.66, padj 1.1e-74) in CD14+ Monocytes, ISG15 (log2FC 5.54, padj 9.4e-140) in B cells, ISG15 (log2FC 4.69, padj 6.3e-175) in CD16+ NK cells, and ISG15 (log2FC 4.77, padj 1e-152) in CD8+ Effector T cells, alongside consistent induction of ISG20 (log2FC 5.72, padj 1.1e-110), IFIT1 (log2FC 8.65, padj 4.4e-59), MX1 (log2FC 6.30, padj 3.4e-49), LY6E (log2FC 5.96, padj 1.6e-38) and IFI6 (log2FC 4.95, padj 1.2e-30) in monocytes, and matching induction of the same genes in every other tested type (e.g. ISG15 (log2FC 4.94, padj 2.3e-77) in Tcm/Naive helper T cells, ISG15 (log2FC 6.90, padj 5.9e-118) in DC2, ISG15 (log2FC 4.89, padj 3.5e-30) in Platelets). MHC-I component B2M (log2FC 0.79, padj 8.3e-30) was also consistently upregulated. CD14+ Monocytes additionally showed strong, more monocyte-specific induction of the chemokine CCL8 (log2FC 10.00, padj 5.7e-127), the anti-inflammatory cytokine antagonist IL1RN (log2FC 6.66, padj 3.4e-116), and the antiviral deaminase APOBEC3A (log2FC 6.06, padj 7e-64) (also induced in Non-classical monocytes: APOBEC3A (log2FC 4.19, padj 4.8e-100)), reflecting the well-documented sensitivity of monocytes to IFN-β. CD14+ Monocytes showed by far the largest number of significant genes (2,419 of 3,780 tested), consistent with monocytes being a principal IFN-responsive population, while smaller lymphocyte and platelet populations showed proportionally fewer significant genes, partly reflecting lower statistical power from fewer cells/pseudobulk counts.

![Volcano plots per cell type](figures/de_volcano.png)

*Pseudobulk differential expression, stim vs ctrl, per cell type (PyDESeq2, design `~donor + condition`, Wald test). Red: up in stim; blue: down (padj < 0.05); grey: not significant. The three most significant up-regulated genes are labelled. Not tested (too few samples): Late erythroid, Mast cells, Tem/Trm cytotoxic T cells, pDC.*

## Caveats

- **Condition is fully confounded with the 10x run**: stim and ctrl PBMCs were captured as two separate captures, so any residual run-specific technical effect (loading, reagent lot, capture efficiency) cannot be distinguished from the true biological IFN-β response in this DE analysis. The paired donor design and consistent, biologically sensible ISG signature across all cell types make a purely technical explanation unlikely, but this cannot be formally ruled out.
- Mitochondrial-based QC was uninformative here because mitochondrial genes were absent from the input matrix, so this dataset lacks that avenue for identifying stressed/dying cells.
- The scVI batch correction (donor x condition) was chosen to align shared cell types across stim/ctrl, which is necessary given the strong, ubiquitous ISG program, but a joint donor+condition batch key risks over-merging a genuinely condition-specific cell state into an existing type; the UMAP colored by condition should be inspected to confirm mixing without erasing a true stimulated-only subpopulation.
- Very small clusters (Platelets, pDC, Late erythroid, Tem/Trm cytotoxic T cells, Mast cells) have limited statistical power for both composition and DE testing, and several were skipped from DE entirely for having too few samples with sufficient cells per condition.
- The Mast cell label rests on a single significant marker (GATA2) in a cluster of only 21 cells (0.1%), so this identity should be treated as tentative.

## Conclusions

Clustering on a donor+condition-corrected scVI embedding resolved 24,558 cells into 16 clusters spanning the expected major PBMC lineages (CD14+ and non-classical monocytes, DCs, B cells, NK cells, several CD4/CD8 T-cell memory subsets, platelets, and a small erythroid/mast contamination), with three clusters requiring correction of biologically implausible or contradicted CellTypist labels before analysis. IFN-β stimulation did not significantly shift cell-type proportions in any lineage, but it strongly and consistently reprogrammed gene expression across every cell type via the same core ISG module, with monocytes showing the broadest and most quantitatively pronounced response and additional monocyte-specific chemokine/cytokine induction — a pattern well aligned with the known systemic action of type-I interferon on peripheral immune cells.

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
