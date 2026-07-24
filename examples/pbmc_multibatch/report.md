
# PBMC Multi-Batch scRNA-seq Analysis Report

## Overview
This report summarizes the end-to-end analysis of `data/pbmc_multibatch.h5ad`, a peripheral blood mononuclear cell (PBMC) dataset spanning **2 batches**. The workflow proceeded from raw counts through QC, doublet removal, batch-corrected dimensionality reduction, clustering, and cell-type annotation.

## 1. Data Overview
- Starting dataset: **4,000 cells × 15,792 genes**, raw counts available, 2 batches present in `obs['batch']`.
- Gene identifiers confirmed as **human gene symbols**, with mitochondrial genes identifiable via the `MT-` prefix (13 mito genes found).

## 2. Quality Control
Per-cell QC metrics were computed (genes/cell, total counts, %mitochondrial reads):
- Median genes/cell: 1,568 (range 264–4,360)
- Median total counts: 4,782 (range 516–19,945)
- Median %mito: 6.9% (range 1.3–20.0%)

**Threshold selection:** The tool-suggested tutorial defaults (`min_genes=200`, `max_pct_mt=5.0`, `min_cells=3`) were evaluated first. The default mito cap of 5% would have discarded **3,291 of 4,000 cells (82%)** — clearly inappropriate, since the dataset's own median mito% (6.9%) already exceeds that cutoff, indicating a naturally higher baseline mitochondrial content than the tutorial assumption. Instead, thresholds were set from the dataset's own distribution: `min_genes=200` (removes 0 cells, safe floor), `max_pct_mt=15%` (trims only the extreme high-mito tail beyond the 99th percentile region), and `min_cells=3` for gene filtering.

**Result:** 73 cells and 90 genes removed → **3,927 cells × 15,702 genes** retained.

## 3. Doublet Detection
Scrublet was run on the filtered data. The doublet-score distribution was **not bimodal** (no clear valley), so per protocol the **median + 3×MAD** rule was applied: threshold = **0.0551**. This flagged 509 cells (12.96%) as likely doublets, which were removed.

**Result:** 3,418 cells retained after doublet removal.

## 4. Normalization
Raw counts were stashed, followed by total-count normalization (target sum 1e4), log1p transformation, and highly-variable gene selection (2,000 HVGs flagged) — standard preprocessing prior to dimensionality reduction.

## 5. Dimensionality Reduction
Since the dataset contains a genuine `batch` covariate with 2 distinct batches, **scVI** was used (rather than plain PCA) to learn a batch-corrected 10-dimensional latent embedding from raw counts restricted to the HVG set (400 training epochs). This representation (`X_scVI`) was used for all downstream neighbor graph, clustering, and UMAP computation.

## 6. Clustering
Leiden clustering (resolution = 1.0) on the scVI latent space yielded **13 clusters**, ranging in size from 27 to 707 cells.

## 7. Marker Genes & Cell-Type Annotation
Wilcoxon rank-sum marker identification followed by CellTypist (Immune_All_Low model, majority voting per cluster) produced the following annotation:

| Cluster | Size | Top Markers | Cell Type |
|---|---|---|---|
| 0 | 369 | S100A9, S100A8, S100A12, VCAN, LYZ | Classical monocytes |
| 1 | 609 | FTH1, CTSS, NEAT1, CST3, FGL2 | Classical monocytes |
| 2 | 388 | GNLY, NKG7, PRF1, KLRD1, CST7 | CD16+ NK cells |
| 3 | 115 | MS4A1, CD79A, HLA-DQA1, CD79B, BANK1 | Memory B cells |
| 4 | 707 | RPS3A, RPL30, RPL32, RPL11, RPS13 | Tcm/Naive helper T cells |
| 5 | 155 | CD8B, RPS12, RPS3A, RPS8, RPL32 | Tcm/Naive cytotoxic T cells |
| 6 | 465 | TRAC, IL7R, IL32, LTB, LDHB | Tem/Effector helper T cells |
| 7 | 149 | GZMK, KLRB1, IL7R, DUSP2, LYAR | MAIT cells |
| 8 | 27 | MALAT1, MTRNR2L12, MT-CO1, MT-ATP6, MT-ND5 | Tcm/Naive helper T cells* |
| 9 | 31 | FCER1A, HLA-DRB1, HLA-DPB1, CD1C, HLA-DPA1 | DC2 |
| 10 | 167 | IGHM, CD79A, IGHD, CD37, CD74 | Naive B cells |
| 11 | 182 | CCL5, CST7, NKG7, GZMA, CD3D | Tem/Trm cytotoxic T cells |
| 12 | 54 | LST1, AIF1, COTL1, MS4A7, FCGR3A | Non-classical monocytes |

*Cluster 8 caveat: its top markers are dominated by mitochondrial genes (MT-CO1/ATP6/ND5) and MALAT1 rather than clear lineage markers, and it is the smallest cluster (27 cells). This signature is characteristic of stressed/low-quality cells that passed the mito filter rather than a distinct biological population; its CellTypist label (Tcm/Naive helper T cells) should be treated with low confidence.

### Overall cell-type composition (n = 3,418 cells)
| Cell Type | Count |
|---|---|
| Classical monocytes | 978 |
| Tcm/Naive helper T cells | 734 |
| Tem/Effector helper T cells | 465 |
| CD16+ NK cells | 388 |
| Tem/Trm cytotoxic T cells | 182 |
| Naive B cells | 167 |
| Tcm/Naive cytotoxic T cells | 155 |
| MAIT cells | 149 |
| Memory B cells | 115 |
| Non-classical monocytes | 54 |
| DC2 | 31 |

## 8. Conclusions
The dataset represents a typical human PBMC sample with all expected major lineages present: monocytes (classical + non-classical), T-cell subsets (CD4 naive/effector, CD8 naive/cytotoxic, MAIT), B cells (naive + memory), NK cells, and dendritic cells (DC2). Batch effects between the two input batches were corrected using scVI prior to clustering, and marker-gene profiles are consistent with canonical PBMC identities. One small cluster (cluster 8, 27 cells) appears to reflect low-quality/stressed cells rather than a true cell type and should be interpreted cautiously or excluded in downstream analyses.

Processed data, UMAP embeddings, and QC figures accompany this report, and the final annotated `.h5ad` object has been saved.

## Figures

![umap](figures/umap.png)

![qc_violin](figures/qc_violin.png)
