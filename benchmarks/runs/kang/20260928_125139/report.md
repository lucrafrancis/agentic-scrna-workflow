# Single-cell RNA-seq analysis: kang

**kang** · 24,673 → 24,558 cells · 16 clusters · 14 cell types · scVI embedding

## Overview

This analysis processes PBMCs from 8 SLE patients (Kang et al., Nat Biotechnol), each split into a paired IFN-β-stimulated and unstimulated (control) sample, to identify cell types and test how stimulation shifts cell-type composition and gene expression. Starting from 24,673 cells x 15,706 genes, the pipeline went through QC filtering, doublet removal, normalization, batch-corrected embedding, clustering, CellTypist annotation (with manual correction where markers disagreed), and paired condition comparisons, ending with 24,558 cells across 14 annotated cell types.

A key structural fact drives several downstream choices: the two conditions were captured in **two separate 10x runs** (stim pooled, ctrl pooled), so `condition` is fully confounded with the technical (10x run) batch. Donors were later assigned by genotype demultiplexing, which also removed cross-donor doublets.

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
| Cluster labels | CellTypist majority vote | Cluster 0/1: CellTypist called these 'Intermediate macrophages', but PBMC blood shouldn't contain tissue macrophages; markers show canonical classical monocyte identity (CD14 log2FC 2.35/4.58, pct 39.6%/64.0%; LYZ pct 84.1%/72.3%; FCN1 pct 78.9%/54.3%), so relabeled CD14+ Monocytes. Cluster 4: labeled 'NK cells' but CD3D (pct 66.4% vs 29.2% elsewhere, log2FC 2.34, padj~0) and CD8A/CD8B (pct 40.1%/30.7%, log2FC 4.46/3.74) are clearly enriched, which is incompatible with NK identity (real NK cluster 11 has CD3D pct only 4.5%, strongly depleted); the GNLY/NKG7/GZMB signature reflects a cytotoxic CD8+ effector T-cell state, so relabeled CD8+ T cells. Cluster 12: labeled a helper-T subset but dominated by platelet markers PPBP (pct 85.1%, log2FC 7.7), PF4 (pct 63.6%, log2FC 7.46), and ITGA2B (log2FC 6.76), so relabeled Platelets. | cluster 0: Intermediate macrophages → **CD14+ Monocytes**; cluster 1: Intermediate macrophages → **CD14+ Monocytes**; cluster 4: NK cells → **CD8+ T cells**; cluster 12: Tcm/Naive helper T cells → **Platelets** | 8,381 cells relabelled |
| Composition test | — | 8 stim / 8 ctrl samples | paired Wilcoxon signed-rank by 'donor' | 0 of 14 cell types changed (padj < 0.05) |
| DE design | ~ condition | replicates in 'donor' | **~ donor + condition** (paired) | 10 cell types tested, 6,428 significant genes in total |
| Min cells per pseudobulk sample | 10 | — | 10 (default kept) | 4 cell types skipped |

## Quality control

No mitochondrial genes were detected in this dataset (0 found; 0% pct_mt at every percentile), meaning it had already been stripped of MT genes upstream — the mito filter is effectively a no-op here rather than a real QC signal, so I relied on the gene-count floor instead. I used the standard tutorial-default thresholds (min_genes=200, max_pct_mt=5%, min_cells=3), which is reasonable given the counts distribution (median genes/cell = 519, median total counts = 1,246) shows no bimodal low-quality tail beyond the standard floor. This removed only 111 cells (0.4%) and 5 genes.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (24,673 cells). The data contains no mitochondrial reads (no mitochondrial genes), so the mitochondrial filter removes nothing. Minimum genes per cell: 200, which removes 111 cells (fewest observed: 12).*

## Doublet detection

The dataset description states doublets were largely already removed by genotype-based demultiplexing (demuxlet), except doublets formed from two cells of the *same* donor, which don't produce a mixed-genotype signal. This is the documented exception case: a strict median+3×MAD cutoff (0.102) would over-remove genuine cells from a population that's already mostly clean. Scrublet was run separately per 10x run (`condition`, since donors were pooled within each run rather than being separate runs), giving per-run automatic thresholds of 0.276 / 0.648 (ctrl / stim). I used the higher of the two (0.648, from the stim run) as a light-touch cutoff, which removed only 4 cells (0.02%) — consistent with the expectation that few residual doublets remain.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 24,562 cells; the distribution is not bimodal. Chosen threshold 0.648 (Scrublet auto, stim run) flags 4 cells (0.0%). The mean of Scrublet's per-batch automatic thresholds (dashed, 0.462) would flag 47.*

## Dimensionality reduction

Because `condition` is confounded with the 10x run, and IFN-β stimulation is known to cause a strong, cell-type-wide transcriptional shift, clustering directly on uncorrected PCA risks splitting each cell type into stim/ctrl-specific clusters, which would make annotation inconsistent between conditions and any composition test meaningless. I used **scVI** with batch_key = `['donor', 'condition']` (one batch per sample), correcting for both donor-to-donor variation and the condition/run effect in the embedding used for clustering and UMAP, while leaving raw counts untouched for the downstream DE and composition analyses (which is exactly where the biological stim-vs-ctrl signal needs to be preserved). This is a deliberate trade-off: over-correction could in principle merge a stimulation-specific state into another type, so I checked whether clusters mixed by condition (see Clustering below) and rely on the pseudobulk DE — run on raw counts — to recover the real condition effect regardless of how the embedding was corrected.

## Clustering

Leiden clustering (resolution 1.0) on the scVI embedding produced 16 clusters ranging from 21 to 3,982 cells. Clusters correspond well to expected PBMC lineages (monocyte, T-cell, NK, B-cell, DC, platelet, erythroid and mast-cell clusters), and — reassuringly for the batch-correction choice — clusters group by cell identity rather than by condition, indicating the scVI correction achieved cross-condition alignment without obviously collapsing distinct states (see UMAP by condition in the figures).

| Cluster | Cells | Top marker genes | Cell type | CellTypist label |
|---|---|---|---|---|
| 0 | 3,647 | FTL, CD63, SOD2, LGALS3, ANXA5 | **CD14+ Monocytes** | Intermediate macrophages |
| 1 | 2,043 | TIMP1, FTH1, TYROBP, C15orf48, FTL | **CD14+ Monocytes** | Intermediate macrophages |
| 2 | 1,097 | FCGR3A, MS4A7, CXCL16, LST1, S100A11 | Non-classical monocytes | Non-classical monocytes |
| 3 | 1,348 | EIF1, BTG1, UBC, UBB, H3F3B | Tcm/Naive helper T cells | Tcm/Naive helper T cells |
| 4 | 2,449 | CCL5, B2M, HLA-A, NKG7, TMSB4X | **CD8+ T cells** | NK cells |
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

CellTypist (Immune_All_Low.pkl, with Immune_All_High.pkl as a second opinion) gave an initial label for each cluster. Before accepting these, I checked canonical markers for every final cell type with `check_markers`, which is what caught three clear mislabels:

- **Clusters 0 and 1** were called "Intermediate macrophages" by CellTypist, but PBMC blood samples do not contain tissue-resident macrophages. Both clusters show canonical classical-monocyte markers instead (CD14 (log2FC -2.59, padj 2.3e-43) and strong LYZ/FCN1 expression), so they were relabeled **CD14+ Monocytes**.
- **Cluster 4** was called "NK cells", but CD3D (log2FC -0.35, padj 0.0072) and CD8A (log2FC -0.03, padj 0.9)/CD8B are all clearly enriched in this cluster (unlike the genuine NK cluster 11, where CD3D is strongly depleted) — real NK cells are CD3-negative. This cluster is a cytotoxic/effector **CD8+ T cell** population (its high GNLY/NKG7/GZMB reflects an effector phenotype, not NK identity), relabeled accordingly.
- **Cluster 12** was called a helper-T subset, but it is dominated by PPBP (log2FC 0.68, padj 0.51) and PF4 (log2FC 0.33, padj 0.96), unambiguous platelet markers, so it was relabeled **Platelets**.

All other clusters' CellTypist labels were confirmed by canonical markers and left unchanged (e.g. MS4A1 (log2FC -0.43, padj 0.055) and CD79A for B cells; GNLY (log2FC 0.35, padj 0.012) and NKG7/KLRD1 for NK cells; FCGR3A/MS4A7 for non-classical monocytes; TCF4/IL3RA/CLEC4C for pDC; FCER1A/CLEC10A for DC2; CCR7/SELL for the naive/central-memory T-cell clusters; GZMK/GZMA for the Tem/Trm cytotoxic cluster; HBB for the erythroid cluster). The small Mast cell cluster is supported by GATA2, a master mast-cell/basophil transcription factor that was strongly and significantly enriched (essentially absent outside the cluster); its tryptase gene TPSAB1 trends in the same direction but does not reach significance given the very small size of this cluster (21 cells (0.1%)).

Final cell-type composition:

| Cell type | Cells | Share |
|---|---|---|
| CD14+ Monocytes | 5,690 | 23.2% |
| Tcm/Naive helper T cells | 5,330 | 21.7% |
| Tem/Effector helper T cells | 3,231 | 13.2% |
| B cells | 2,626 | 10.7% |
| CD8+ T cells | 2,449 | 10.0% |
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

*Canonical markers checked by the agent before accepting or changing labels (27 of 39 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: CD4, NCAM1, CD1C, FCER1A, CLEC10A, CD19, ITGA2B, TPSAB1, KIT, HDC, VWA5A, SLC18A2.*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Composition analysis

Cell-type proportions were compared per donor (paired, since every donor contributed both conditions) using a Wilcoxon signed-rank test. **No cell type showed a significant proportion shift** after correction (0 of 14 significant). The largest (non-significant) trends were a decrease in DC2 proportion in stim (2.2% → 1.5%, log2 ratio -0.52, padj 0.33) and an increase in pDC proportion (0.3% → 0.8%, log2 ratio 1.27, padj 0.68), but with only 8 paired donors the test cannot detect anything short of a large, consistent effect. This null result is consistent with the short (6h) stimulation window being long enough to reprogram gene expression but not to substantially alter the mix of circulating cell types.

| Cell type | Mean share, ctrl | Mean share, stim | log2 ratio | padj |
|---|---|---|---|---|
| DC2 | 2.2% | 1.5% | -0.52 | 0.33 |
| CD14+ Monocytes | 25.4% | 24.4% | -0.06 | 0.68 |
| CD8+ T cells | 9.9% | 9.4% | -0.08 | 0.68 |
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

![Cell-type composition](figures/composition.png)

*Cells per annotated type (24,558 cells, 14 types; CellTypist majority vote over Leiden clusters, with clusters relabelled by the agent from their markers (see the decisions table)).*

![Cell-type proportions by condition](figures/composition_change.png)

*Share of each cell type per sample, ctrl (blue) vs stim (orange); lines join each donor's two samples. Wilcoxon signed-rank test with Benjamini-Hochberg correction: 0 of 14 cell types differ at padj < 0.05. Log scale above 0.1%, linear below, so samples with none of a type sit at 0.*

## Differential expression

Pseudobulk DE (PyDESeq2, design `~donor + condition`, paired by donor) was run per cell type on raw counts. Every cell type with sufficient samples showed a robust interferon-stimulated gene (ISG) signature: ISG15 (log2FC 7.66, padj 1.1e-74), ISG20 (log2FC 5.72, padj 1.1e-110), and other canonical ISGs (IFI6 (log2FC 5.02, padj 2e-102) etc.) were significantly upregulated in every tested cell type, essentially independent of lineage — the expected, textbook response to IFN-β. CD14+ Monocytes had by far the largest and most complex response (n_significant = 2,419 of 3,780 genes tested, roughly balanced up/down), consistent with monocytes being major IFN-responsive antigen-presenting cells; they also induced strong chemokine responses such as CCL8 (log2FC 10.00, padj 5.7e-127) and IL1RN, not seen as prominently in lymphocytes. Non-classical monocytes and DC2 also showed large, significant responses (817 and 525 significant genes respectively). T- and B-lymphocyte subsets showed smaller but still highly significant ISG inductions (e.g. B cells: 620 significant genes; CD8+ T cells: 281). Downregulated genes were dominated by ribosomal protein genes (RPL/RPS) and general translation/iron-storage genes (e.g. FTH1), a common pattern in activated-cell pseudobulk contrasts and likely reflecting a shift in cellular metabolic state rather than a specific IFN-driven repression program.

Late erythroid, Mast cells, Tem/Trm cytotoxic T cells and pDC were skipped for DE because too few donors had at least 10 cells of that type in one or both conditions — these are the rarest clusters in the dataset and the composition/DE tests are underpowered for them regardless of biology.

| Cell type | Samples (test / ref) | Genes tested | Up | Down | Top up-regulated genes |
|---|---|---|---|---|---|
| B cells | 8 / 8 | 2,101 | 347 | 273 | ISG20, B2M, ISG15, LY6E, UBE2L6 |
| CD14+ Monocytes | 8 / 8 | 3,780 | 1,221 | 1,198 | SSB, CCL8, IL1RN, DYNLT1, ISG20 |
| CD16+ NK cells | 8 / 8 | 1,450 | 220 | 169 | ISG15, ISG20, TNFSF10, IFI6, IFIT1 |
| CD8+ T cells | 8 / 8 | 1,248 | 165 | 116 | ISG20, ISG15, IFI6, SAT1, IFIT1 |
| DC2 | 7 / 8 | 1,085 | 296 | 229 | ISG15, ISG20, HSP90AA1, LY6E, MX1 |
| Non-classical monocytes | 8 / 8 | 1,612 | 425 | 392 | APOBEC3A, TNFSF10, MYL12A, ISG20, ISG15 |
| Platelets | 5 / 4 | 339 | 30 | 3 | ISG15, ISG20, IFI6, B2M, LY6E |
| Tcm/Naive cytotoxic T cells | 7 / 6 | 1,214 | 130 | 117 | ISG15, ISG20, LY6E, MX1, IFIT1 |
| Tcm/Naive helper T cells | 8 / 8 | 3,144 | 369 | 293 | ISG20, PSMB9, IFI16, LY6E, B2M |
| Tem/Effector helper T cells | 8 / 8 | 2,340 | 283 | 152 | ISG20, SAT1, MT2A, TMSB10, TNFSF10 |

![Volcano plots per cell type](figures/de_volcano.png)

*Pseudobulk differential expression, stim vs ctrl, per cell type (PyDESeq2, design `~donor + condition`, Wald test). Red: up in stim; blue: down (padj < 0.05); grey: not significant. The three most significant up-regulated genes are labelled. Not tested (too few samples): Late erythroid, Mast cells, Tem/Trm cytotoxic T cells, pDC.*

## Caveats

- **Condition is fully confounded with the 10x run** (stim and ctrl were captured in separate captures). Any batch-specific technical effect (e.g. run-to-run efficiency differences) cannot be distinguished from the true condition effect in this design; the pseudobulk DE and composition results should be read with this in mind, although the consistency of the signature with the well-established IFN-β biology, and its presence across every lineage, makes a purely technical explanation unlikely.
- The scVI batch correction used to build the clustering/UMAP embedding included `condition` in the batch key. This choice helps align cell types across conditions for annotation and composition testing, but risks merging a condition-specific cell state into an existing cluster; DE analyses use raw counts and are unaffected, but any subtler stimulation-induced subpopulation could be underrepresented in the cluster structure itself.
- The mitochondrial-percentage QC metric was uninformative here (no MT genes present in the object), so a standard cause of low-quality-cell removal was unavailable; only the gene-count floor was used to catch low-complexity cells.
- Very rare clusters (Mast cells, n=21 cells (0.1%); Tem/Trm cytotoxic T cells, n=30 cells (0.1%); Late erythroid, n=64 cells (0.3%)) rest on small numbers of cells and were excluded from DE/composition testing; their labels, while marker-supported, should be treated as lower-confidence than the major lineages.
- Composition testing with only 8 paired donors has limited power; a true but modest compositional shift could easily go undetected.

## Conclusions

Fourteen PBMC cell types were identified and annotated, with three CellTypist labels corrected after marker inspection (two monocyte clusters mislabeled as macrophages, one CD8+ T-cell cluster mislabeled as NK cells, and one platelet cluster mislabeled as helper T cells). Six-hour IFN-β stimulation did not significantly change the proportions of any PBMC cell type, but it induced a strong, highly significant, and pervasive interferon-stimulated gene signature across essentially every cell type, with the largest transcriptional response in classical (CD14+) monocytes, non-classical monocytes and DC2 — consistent with these myeloid populations' central role in mounting and amplifying an innate antiviral response.

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
