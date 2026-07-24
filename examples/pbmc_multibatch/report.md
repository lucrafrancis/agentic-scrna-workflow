# PBMC Multi-Batch scRNA-seq Analysis Report

## Overview
This report describes the end-to-end analysis of `data/pbmc_multibatch.h5ad`, a peripheral blood mononuclear cell (PBMC) dataset comprising two batches (technical replicates/samples, `PBMC10k` and `PBMC5k`). The analysis proceeded from raw counts through QC filtering, doublet removal, batch-corrected dimensionality reduction, clustering, and cell-type annotation.

## 1. Input Data
- **Starting size:** 4,000 cells × 15,792 genes
- **Gene identifiers:** Human gene symbols (mitochondrial prefix `MT-`, 13 mito genes detected)
- **Batch structure:** 2 batches present in `obs['batch']`, treated as a technical grouping to be corrected for (not an experimental variable of interest)
- Raw counts were available and preserved throughout for later use by scVI.

## 2. Quality Control
Per-cell QC metrics were computed (genes/cell, total counts/cell, % mitochondrial reads). Key distributional stats:
- `n_genes_by_counts`: median 1568 (range 264–4360)
- `total_counts`: median 4783 (range 516–19,945)
- `pct_counts_mt`: median 6.9% (range 1.3%–20.0%)

**Threshold decision:** The standard tutorial-default mitochondrial cutoff (5%) was evaluated but rejected — it would have discarded 3,291/4,000 cells (82%), since the dataset's median mito% (6.9%) already exceeds that threshold. This indicated the tutorial default was miscalibrated for this dataset's biology/chemistry. Instead, a data-driven cutoff of **15% mito** (between the 95th and 99th percentiles) was used, removing only genuine high-mito outliers.

**Filters applied:**
- `min_genes = 200` (floor; removed 0 cells, as the minimum was already 264)
- `max_pct_mt = 15%`
- `min_cells = 3` (gene-level filter)

**Result:** 4,000 → 3,927 cells (73 removed); 15,792 → 15,702 genes (90 removed).

## 3. Doublet Detection
Scrublet was run on the filtered data. The doublet-score distribution was **not bimodal**, so per protocol the **median + 3×MAD** rule was used rather than a bimodal valley or the Scrublet automatic threshold (the latter, computed per-batch at ~0.40, appeared miscalibrated — far out in the sparse tail and would have flagged almost no cells despite a clearly heavy-tailed score distribution).

- **Threshold chosen:** 0.09 (median 0.038 + 3×MAD)
- **Cells flagged and removed:** 364 (9.3%)
- **Result:** 3,927 → 3,563 cells

## 4. Normalization
Raw counts were stashed (for scVI), then data were total-count normalized (target sum 1e4), log1p-transformed, and the top 2,000 highly variable genes (HVGs) were flagged (computed per-batch to avoid batch-driven HVG selection bias).

## 5. Dimensionality Reduction
Since the dataset contains 2 batches representing a **technical grouping** (separate samples/sequencing runs) rather than a biological variable of interest, **scVI** was used instead of plain PCA to learn a batch-corrected latent embedding. scVI was trained on raw counts restricted to the HVGs (10 latent dimensions, 400 epochs), producing the `X_scVI` representation used for all downstream clustering/visualization.

## 6. Clustering
Leiden clustering (resolution = 1.0) on the scVI latent space yielded **13 clusters**, ranging from 23 to 712 cells, with UMAP coordinates computed for visualization.

## 7. Marker Genes & Cell-Type Annotation
Wilcoxon rank-sum marker analysis and CellTypist (`Immune_All_Low.pkl`, majority-vote per cluster) were used to annotate clusters. Marker genes were highly concordant with the assigned identities:

| Cluster | Top Markers | Annotated Cell Type | n Cells |
|---|---|---|---|
| 0 | S100A9, S100A8, S100A12, VCAN, CD14 | Classical monocytes | 658 |
| 1 | CPVL, NEAT1, PSAP, CST3, FGL2 | Classical monocytes | 385 |
| 2 | GNLY, NKG7, PRF1, KLRD1, GZMB | CD16+ NK cells | 399 |
| 3 | CD79A, MS4A1, CD37, CD79B | Naive B cells | 288 |
| 4 | RPS3A, RPL30, RPL32 (ribosomal) | Tcm/Naive helper T cells | 712* |
| 5 | CD8B, ribosomal genes | Tcm/Naive cytotoxic T cells | 172 |
| 6 | HLA-DRB1/DPB1/DRA, FCER1A, CLEC10A | DC2 | 64 |
| 7 | IL32, TRAC, IL7R, CD3D | Tem/Effector helper T cells | 445 |
| 8 | MALAT1, JUN, MT-CO1/ATP6/ND5 (mito-high) | Tcm/Naive helper T cells (low-confidence) | 24 |
| 9 | CCL5, GZMA, CST7, NKG7, CD3D | Tem/Trm cytotoxic T cells | 192 |
| 10 | KLRB1, GZMK, IL7R, NKG7 | MAIT cells | 132 |
| 11 | LST1, AIF1, FCGR3A, MS4A7, CSF1R | Non-classical monocytes | 69 |
| 12 | IL3RA, PLD4, IRF8, TCF4, MZB1 | pDC | 23 |

*Clusters 4 and 8 were both majority-called "Tcm/Naive helper T cells"; combined count shown in totals below.

**Overall cell-type composition (n=3,563 cells):**
| Cell Type | Count | % |
|---|---|---|
| Classical monocytes | 1,043 | 29.3% |
| Tcm/Naive helper T cells | 736 | 20.7% |
| Tem/Effector helper T cells | 445 | 12.5% |
| CD16+ NK cells | 399 | 11.2% |
| Naive B cells | 288 | 8.1% |
| Tem/Trm cytotoxic T cells | 192 | 5.4% |
| Tcm/Naive cytotoxic T cells | 172 | 4.8% |
| MAIT cells | 132 | 3.7% |
| Non-classical monocytes | 69 | 1.9% |
| DC2 | 64 | 1.8% |
| pDC | 23 | 0.6% |

## 8. Caveats & Notes
- **Cluster 8** (24 cells) shows a marker profile dominated by MALAT1, JUN, and multiple mitochondrial genes (MT-CO1, MT-ATP6, MT-ND5, MT-CYB, MT-ND4, MT-CO3) rather than canonical lineage markers. This is a hallmark of stressed or partially degraded cells rather than a distinct biological population, despite being majority-labeled "Tcm/Naive helper T cells" by CellTypist. This small cluster should be interpreted with caution and could be excluded in follow-up analyses.
- The resulting composition (monocyte-dominant, with the expected T/B/NK/DC/pDC compartments) is consistent with a typical healthy PBMC sample.
- Batch correction via scVI was applied given the 2-batch technical structure; no batch-driven clustering artifacts were apparent among the major lineage clusters.

## 9. Pipeline Summary
| Step | Cells | Genes |
|---|---|---|
| Raw input | 4,000 | 15,792 |
| After QC filtering (min_genes=200, max_pct_mt=15%, min_cells=3) | 3,927 | 15,702 |
| After doublet removal (threshold=0.09) | 3,563 | 15,702 |
| Final (post-normalization/HVG/clustering) | 3,563 | 15,702 (2,000 HVGs used for scVI/PCA) |

Final annotated object, figures (UMAP, QC plots), and this report have been saved as part of the pipeline output.

## Figures

![umap](figures/umap.png)

![marker_dotplot_clusters](figures/marker_dotplot_clusters.png)

![marker_dotplot_celltype](figures/marker_dotplot_celltype.png)

![qc_violin](figures/qc_violin.png)
