
# PBMC3k Single-Cell RNA-seq Analysis Report

## Overview
This report summarizes the end-to-end analysis of the PBMC3k dataset (peripheral blood mononuclear cells), from raw counts through quality control, doublet removal, normalization, clustering, and cell-type annotation.

## 1. Dataset & Gene Identifiers
- **Input**: 2,700 cells × 32,738 genes, raw counts available, no batch key present (single, clean batch).
- **Gene identifiers**: Human gene symbols confirmed. Mitochondrial genes identified via the `MT-` prefix (13 mitochondrial genes found).

## 2. Quality Control
Per-cell QC metrics were computed (genes per cell, total counts, % mitochondrial reads):
- Genes/cell: median 817 (range 212–3422)
- Total counts/cell: median 2197 (range 548–15844)
- % mitochondrial: median 2.0% (p95 = 4.0%, p99 = 5.9%, max = 22.6%)

**Thresholds applied** (standard tutorial defaults, well-matched to this dataset's distributions):
- Minimum genes/cell: 200 (removed 0 cells — all cells already exceeded this floor)
- Maximum mitochondrial %: 5.0% (removed 57 high-mito, likely stressed/dying cells)
- Minimum cells per gene: 3 (removed 19,041 genes detected too sparsely to be informative)

**Result**: 2,643 cells × 13,697 genes retained after filtering.

## 3. Doublet Detection
Scrublet was run to score putative doublets. The score distribution was **unimodal** (not bimodal), so per protocol the **median + 3×MAD** rule (threshold = 0.10) was used rather than a bimodal-valley cutoff. This flagged 222 cells (8.4%) sitting in the long right tail of the score distribution (scores above ~0.11, up to a max of 0.51), consistent with heterotypic doublets.

**Result**: 2,421 cells retained after doublet removal.

## 4. Normalization
Raw counts were stashed for later use (e.g., scVI-style tools), then data were total-count normalized (target sum 1e4), log1p-transformed, and 2,000 highly variable genes flagged for downstream dimensionality reduction.

## 5. Dimensionality Reduction
No batch key was present in this dataset (single clean batch), so **PCA** was chosen over scVI as the simpler, appropriate method — batch correction is unnecessary here. PCA was run on the HVG matrix (50 components); the top PCs captured the dominant structure (PC1 = 10.3%, PC2 = 3.5%, PC3 = 2.5% variance), typical for PBMC data with well-separated major lineages.

## 6. Clustering
Leiden clustering (resolution = 1.0) on the PCA representation identified **9 clusters**, ranging in size from 550 to 7 cells. UMAP was computed for visualization.

## 7. Marker Genes & Cell-Type Annotation
Wilcoxon rank-sum marker gene analysis per cluster, combined with CellTypist (Immune_All_Low model, majority vote per cluster), gave clear, mutually consistent identities:

| Cluster | Size | Top Markers | Annotated Cell Type |
|---|---|---|---|
| 0 | 260 | CCL5, NKG7, B2M, GZMA, CST7 | Tem/Trm cytotoxic T cells |
| 1 | 315 | CD74, CD79A, HLA-DRA, CD79B, HLA-DPB1 | B cells |
| 2 | 515 | LTB, IL32, LDHB, IL7R, CD3D | Tcm/Naive helper T cells |
| 3 | 550 | RPS12, RPS6, RPL32, RPS27, RPS3A | Tcm/Naive helper T cells |
| 4 | 446 | LYZ, S100A9, S100A8, TYROBP, FCN1 | Classical monocytes |
| 5 | 41 | CD74, HLA-DRB1, HLA-DRB5, GAPDH, ACTB | DC |
| 6 | 140 | GNLY, NKG7, GZMB, PRF1, CTSW | CD16+ NK cells |
| 7 | 147 | LST1, FCER1G, AIF1, COTL1, FCGR3A | Non-classical monocytes |
| 8 | 7 | PF4, GNG11, PPBP, SDPR, SPARC | Megakaryocytes/platelets |

Marker genes are biologically coherent with the assigned identities:
- **Cytotoxic T / NK clusters** (0, 6) show canonical cytotoxic granule genes (GZMA/GZMB, PRF1, NKG7, GNLY).
- **T helper clusters** (2, 3) show CD3D/IL7R/IL32 (cluster 2) and a ribosomal-gene-high signature (cluster 3), consistent with quiescent/naive T cells.
- **B cells** (1) show CD79A/CD79B/MS4A1/HLA-DR genes.
- **Monocyte subsets** (4, 7) are distinguished by classical (LYZ, S100A8/9, FCN1) vs. non-classical (FCGR3A, LST1, AIF1) monocyte markers.
- **DCs** (5) show high HLA-DR/CD74 antigen-presentation genes.
- **Megakaryocytes/platelets** (8), though a very small cluster (n=7), are unambiguously marked by PF4, PPBP, GNG11, SDPR — classic platelet genes.

## 8. Final Cell-Type Composition

| Cell Type | Count | % of Total |
|---|---|---|
| Tcm/Naive helper T cells | 1,065 | 44.0% |
| Classical monocytes | 446 | 18.4% |
| B cells | 315 | 13.0% |
| Tem/Trm cytotoxic T cells | 260 | 10.7% |
| Non-classical monocytes | 147 | 6.1% |
| CD16+ NK cells | 140 | 5.8% |
| DC | 41 | 1.7% |
| Megakaryocytes/platelets | 7 | 0.3% |

**Total cells in final annotated dataset: 2,421**

## 9. Conclusions
The analysis recovered the expected major PBMC lineages (T cells, B cells, monocytes, NK cells, dendritic cells, and a small platelet contaminant population) with clear, canonical marker gene support and consistent CellTypist annotation. This is in line with expectations for the well-characterized PBMC3k reference dataset. QC and doublet filtering removed a modest fraction of low-quality/doublet cells (2,700 → 2,421, ~10.3% total removed) without distorting the expected cell-type proportions. No batch correction was needed, as the dataset comprises a single sequencing batch.

## Figures

![umap](figures/umap.png)

![marker_dotplot_clusters](figures/marker_dotplot_clusters.png)

![marker_dotplot_celltype](figures/marker_dotplot_celltype.png)

![qc_violin](figures/qc_violin.png)
