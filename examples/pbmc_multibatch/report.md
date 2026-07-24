# PBMC Multi-Batch scRNA-seq Analysis Report

## Overview
Dataset: `data/pbmc_multibatch.h5ad` — 4,000 cells × 15,792 genes, human PBMCs pooled across **2 batches** (`batch` column: PBMC10k, PBMC5k), with raw counts available.

## Gene Identifiers & QC
- Gene symbols confirmed (human), mitochondrial genes identified via **MT-** prefix (13 genes found).
- QC metrics computed per cell: `n_genes`, `total_counts`, `pct_counts_mt`.
- **Mitochondrial content was unusually high genome-wide** (median 6.9%, p75 8.6%, p95 12.2%), well above the typical PBMC tutorial baseline. The tool's default recommendation (`max_pct_mt=5%`) would have discarded **82% of cells** — clearly inappropriate for this dataset's actual biology/library prep. Instead, a threshold of **15% mito** was chosen from the observed distribution, targeting only the extreme tail (~p97+) of stressed/dying cells.
- Filtering applied: `min_genes=200` (removed 0 cells — distribution was already healthy), `max_pct_mt=15`, `min_cells=3` (gene-level).
- Result: **4000 → 3927 cells** (73 removed), 15792 → 15702 genes (90 removed).

## Doublet Detection
- Scrublet doublet scores were computed; the distribution was **not bimodal** (no clear valley), so per protocol the **median + 3×MAD** rule was used: threshold = **0.09**.
- This flagged 364 cells (9.3%) as doublets — a plausible rate for a pooled/multiplexed multi-batch PBMC dataset.
- Result: **3927 → 3563 cells** after removal.

## Normalization & Batch Correction
- Raw counts stashed prior to normalization (required for scVI).
- Total-count normalization (target sum 1e4), log1p transform, and HVG flagging (2000 genes) performed.
- Since the dataset contains **2 real batches**, **scVI** was used (rather than plain PCA) to learn a batch-corrected 10-dimensional latent embedding (`X_scVI`), trained for 400 epochs on raw counts restricted to HVGs. This is the appropriate choice per protocol when a genuine multi-batch structure is present.

## Clustering & Markers
- Leiden clustering (resolution 1.0) on the scVI latent space + UMAP embedding yielded **14 clusters**, ranging from 23 to 703 cells.
- Wilcoxon marker-gene ranking per cluster identified clear, biologically coherent signatures (e.g., S100A8/9/CD14/LYZ for classical monocytes, CD3D/TRAC for T cells, CD79A/MS4A1 for B cells, NKG7/GNLY for NK cells, IL3RA/IRF8 for pDCs).

## Cell Type Annotation (CellTypist, Immune_All_Low model)
Majority-vote annotation per cluster produced 11 immune cell types, matching the manual marker interpretation closely:

| Cell Type | Count | % of cells |
|---|---|---|
| Classical monocytes | 1017 | 28.5% |
| Tcm/Naive helper T cells | 728 | 20.4% |
| CD16+ NK cells | 397 | 11.1% |
| Tem/Effector helper T cells | 318 | 8.9% |
| Naive B cells | 288 | 8.1% |
| MAIT cells | 279 | 7.8% |
| Tem/Trm cytotoxic T cells | 194 | 5.4% |
| Tcm/Naive cytotoxic T cells | 164 | 4.6% |
| Non-classical monocytes | 83 | 2.3% |
| DC2 | 72 | 2.0% |
| pDC | 23 | 0.6% |

### Cluster-to-cell-type mapping
| Cluster | Size | Cell Type | Top markers |
|---|---|---|---|
| 0 | 413 | Classical monocytes | S100A9, S100A8, S100A12, VCAN, S100A6 |
| 1 | 318 | Tem/Effector helper T cells | TRAC, IL32, LDHB, LTB, CD3D |
| 2 | 25 | Tcm/Naive helper T cells* | MALAT1, JUN, MT-CO1, MTRNR2L12, MT-ATP6 |
| 3 | 604 | Classical monocytes | NEAT1, FGL2, CTSS, PSAP, FTH1 |
| 4 | 282 | CD16+ NK cells | GNLY, GZMH, NKG7, CCL5, FGFBP2 |
| 5 | 288 | Naive B cells | CD79A, MS4A1, CD37, CD79B, HLA-DQA1 |
| 6 | 703 | Tcm/Naive helper T cells | RPS3A, RPL30, RPL32, RPL11, RPS27A |
| 7 | 164 | Tcm/Naive cytotoxic T cells | CD8B, RPS12, RPS3A, RPL32, RPS8 |
| 8 | 72 | DC2 | HLA-DRB1, HLA-DPB1, FCER1A, HLA-DPA1, HLA-DRA |
| 9 | 279 | MAIT cells | IL7R, KLRB1, GZMK, IL32, LYAR |
| 10 | 115 | CD16+ NK cells | NKG7, KLRD1, GNLY, KLRF1, CLIC3 |
| 11 | 194 | Tem/Trm cytotoxic T cells | CCL5, GZMA, CST7, NKG7, CD3D |
| 12 | 83 | Non-classical monocytes | LST1, AIF1, COTL1, FTL, FCGR3A |
| 13 | 23 | pDC | IL3RA, CCDC50, PLD4, IRF8, TCF4 |

\* Cluster 2 (25 cells, 0.7% of cells) shows a stress/low-quality signature dominated by MALAT1, JUN, and multiple mitochondrial transcripts rather than canonical T-cell markers. CellTypist assigned it a "Tcm/Naive helper T cell" label by majority vote, but this cluster likely represents residual stressed or partially degraded cells rather than a distinct biological T-cell state. It is small and does not materially affect downstream interpretation, but should be interpreted with caution (or filtered) in follow-up analyses.

## Summary
This PBMC dataset was processed from raw counts through QC, doublet removal, batch-corrected scVI embedding, clustering, and cell-type annotation. Final dataset: **3,563 high-quality cells × 15,702 genes**, organized into 14 transcriptionally distinct clusters spanning the expected major PBMC lineages: monocytes (classical & non-classical), T-cell subsets (CD4 helper, CD8 cytotoxic, MAIT), B cells, NK cells, and dendritic cells (DC2, pDC). The two input batches were successfully integrated via scVI, and cluster identities are well supported by canonical marker genes and independent CellTypist annotation. Key analytical judgment calls — the mito-percentage threshold (adapted from the dataset's own elevated baseline rather than a generic default) and the doublet threshold (median+3×MAD, given a non-bimodal score distribution) — are documented above with rationale.

## Figures

![umap](figures/umap.png)

![marker_dotplot_clusters](figures/marker_dotplot_clusters.png)

![marker_dotplot_celltype](figures/marker_dotplot_celltype.png)

![qc_violin](figures/qc_violin.png)
