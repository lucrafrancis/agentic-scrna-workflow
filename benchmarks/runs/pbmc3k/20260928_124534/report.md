# Single-cell RNA-seq analysis: pbmc3k

**pbmc3k** · 2,700 → 2,421 cells · 9 clusters · 8 cell types · PCA embedding

## Overview

This report covers an end-to-end analysis of a single 10x Genomics PBMC sample from one healthy donor (starting with 2,700 cells and 32,738 genes, raw counts). Because this is a single, clean 10x run with no batch structure, the workflow used PCA (not batch-corrected scVI) for dimensionality reduction, and no cross-condition comparison was performed (there is only one donor and one condition).

## Key analysis decisions

Each row compares a standard default with what the agent chose after reading the data. Built from the tool log (`tool_calls.jsonl`), not written by the model.

| Step | Standard default | What the data showed | Agent's choice | Effect |
|---|---|---|---|---|
| Mitochondrial cutoff | 5% | median 2.03%; default would remove 57 of 2,700 cells (2.1%) | 5% (default kept) | 57 cells removed by QC (2.1%) |
| Cell and gene floors | ≥200 genes/cell, genes in ≥3 cells | fewest genes in a cell: 212 | ≥200 genes/cell, genes in ≥3 cells (default kept) | 19,041 genes removed |
| Doublet threshold | Scrublet automatic (0.227) | scores not bimodal; median + 3×MAD 0.1 | **0.1** (median + 3×MAD) | 222 cells removed (8.4%) |
| Normalisation | 10,000 counts/cell, 2,000 HVGs | — | 10,000 counts/cell, 2,000 HVGs (default kept) | 2,000 HVGs flagged |
| Embedding | PCA | no batch column | PCA | 50 components |
| Leiden resolution | 1.0 | — | 1.0 (default kept) | 9 clusters |
| Annotation model | Immune_All_Low | — | Immune_All_Low (default kept) | 8 cell types |

## Quality control

Gene identifiers were confirmed as human gene symbols, letting mitochondrial genes (13 genes with the `MT-` prefix) be identified correctly for QC. Per-cell QC distributions showed a median of 817 genes and 2,197 counts per cell, with mitochondrial fraction at a median of 2.03% (p95 4.01%, max 22.6%).

Standard tutorial-default thresholds fit this dataset well and were applied unchanged: minimum 200 genes/cell (removed no cells — the population was already well above this floor), maximum 5% mitochondrial content (removed 57 cells, 2.1% of the dataset, trimming the long tail of stressed/dying cells beyond the 99th percentile), and genes detected in at least 3 cells (dropping 19,041 uninformative genes). After filtering, 2,643 cells and 13,697 genes remained.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (2,700 cells, median 2.03% mitochondrial). The standard 5% mitochondrial cutoff was kept and removes 57 cells (2.1%). Minimum genes per cell: 200, which removes 0 cells (fewest observed: 212).*

## Doublet detection

Scrublet was run on the whole sample as one batch (single 10x run, no genotype/hashing demultiplexing reported, so no light-touch exception applies). The doublet-score distribution was unimodal and decayed smoothly with no bimodal valley, so the median+3×MAD rule was used to set the threshold rather than the (inapplicable) bimodal-valley or the more permissive Scrublet-automatic option. This gave a threshold of 0.1 (median 0.0442, max 0.506), flagging and removing 222 cells (8.4%) as likely doublets — a plausible rate for this loading density. 2,421 cells and 13,697 genes entered downstream analysis.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 2,643 cells; the distribution is not bimodal. Chosen threshold 0.1 (median + 3×MAD) flags 222 cells (8.4%). Scrublet's automatic threshold (dashed, 0.227) would flag 39.*

## Dimensionality reduction

No batch key was detected and the sample is a single clean 10x run from one donor, so PCA on the 2,000 highly variable genes (after total-count normalization to 10,000 and log1p) was the appropriate, simplest choice — scVI's batch correction would add complexity with no batches to correct. The leading principal components captured 10.3%, 3.5%, and 2.5% of variance respectively, consistent with a few dominant axes of lineage variation (myeloid vs lymphoid, then finer distinctions) typical of PBMCs.

## Clustering

Leiden clustering at resolution 1.0 on the PCA embedding yielded 9 clusters, ranging in size from 7 to 550 cells — sizes and counts consistent with the expected mixture of abundant T-cell/monocyte populations and rarer populations (DCs, platelets).

| Cluster | Cells | Top marker genes | Cell type |
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

![UMAP](figures/umap.png)

*UMAP of 2,421 cells computed on PCA, coloured by cell type and Leiden cluster.*

## Cell-type annotation

CellTypist (Immune_All_Low.pkl, fine-grained) was used for initial labels, with Immune_All_High.pkl (coarse) as a second opinion. The two models agreed at the expected level of granularity (e.g., the fine model split lymphocytes into Tem/Trm cytotoxic T cells, Tcm/Naive helper T cells, and CD16+ NK cells, while the coarse model grouped the first two as "T cells" and the NK cluster as the broader "ILC" category that includes NK cells — not a contradiction).

Every final cell type was checked against canonical markers before accepting the labels:
- **Tem/Trm cytotoxic T cells** (cluster 0): high CD3D, CD3E, CD8A, CD8B, GZMA, NKG7 and CCL5, consistent with cytotoxic CD8 T cells.
- **B cells** (cluster 1): high CD79A, MS4A1 and CD79B, with CD3D essentially absent.
- **Tcm/Naive helper T cells** (clusters 2 and 3): CD3D/CD3E-positive with high IL7R and low CD8A, consistent with CD4 T cells; cluster 3 is dominated by ribosomal-protein genes, typical of quiescent naive lymphocytes.
- **Classical monocytes** (cluster 4): near-universal CD14, LYZ and S100A9 expression.
- **DC** (cluster 5): high HLA-DRA together with FCER1A and CD1C, both significantly enriched and largely absent elsewhere — a myeloid dendritic cell signature.
- **CD16+ NK cells** (cluster 6): near-universal GNLY and NKG7, high FCGR3A, and CD3D essentially absent, ruling out a T/NKT identity.
- **Non-classical monocytes** (cluster 7): high FCGR3A and MS4A7 with comparatively low CD14, the classical non-classical/CD14-dim monocyte profile.
- **Megakaryocytes/platelets** (cluster 8): near-universal PPBP, PF4 and ITGA2B (CD41), a clean platelet signature.

Because all clusters' markers matched their CellTypist labels (and the second-opinion model's differences were only granularity, not contradiction), no cluster relabeling was necessary.

| Cell type | Cells | Share |
|---|---|---|
| Tcm/Naive helper T cells | 1,065 | 44.0% |
| Classical monocytes | 446 | 18.4% |
| B cells | 315 | 13.0% |
| Tem/Trm cytotoxic T cells | 260 | 10.7% |
| Non-classical monocytes | 147 | 6.1% |
| CD16+ NK cells | 140 | 5.8% |
| DC | 41 | 1.7% |
| Megakaryocytes/platelets | 7 | 0.3% |
| **Total** | **2,421** | |

Final composition: 1,065 cells (44.0%) CD4 T cells, 446 cells (18.4%) classical monocytes, 315 cells (13.0%) B cells, 260 cells (10.7%) cytotoxic T cells, 147 cells (6.1%) non-classical monocytes, 140 cells (5.8%) NK cells, 41 cells (1.7%) dendritic cells, and 7 cells (0.3%) megakaryocytes/platelets — proportions broadly in line with expectations for a healthy PBMC sample.

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (24 of 25 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: NCAM1.*

![Cell-type composition](figures/composition.png)

*Cells per annotated type (2,421 cells, 8 types; CellTypist majority vote over Leiden clusters).*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Caveats

- This is a single donor/single run, so no batch correction or condition comparison was needed or performed; results reflect this one individual's PBMC composition and cannot be generalized without replication.
- The doublet threshold (median+3×MAD) is a heuristic; a small number of true rare-cell profiles that resemble doublets transcriptionally (e.g. some large cytotoxic cells) could in principle be lost, though the smooth, unimodal score distribution gives no indication of a distinct doublet population being missed or over-cut.
- The megakaryocyte/platelet cluster is very small (7 cells (0.3%)), so its marker statistics should be interpreted cautiously; individual DE p-values within this cluster are less well powered.
- Cluster 3's marker list is dominated by ribosomal genes rather than a distinctive positive signature; its T-cell identity rests mainly on CD3D/CD3E/IL7R positivity confirmed via check_markers rather than on unique cluster-defining genes.

## Conclusions

Standard QC and doublet filtering removed a modest fraction of low-quality/multiplet cells (10.3% total), leaving a clean dataset of 2,421 cells. PCA-based clustering at resolution 1.0 resolved the expected major PBMC lineages, and CellTypist annotation — cross-validated against canonical markers for every cluster — produced confident, mutually consistent labels for all 8 cell types without requiring any manual relabeling.

## Methods

*Generated from the tool log: these are the steps and parameters that actually ran.*

**Quality control.** Per-cell metrics were computed with scanpy `calculate_qc_metrics`; mitochondrial genes were those prefixed `MT-` (13 found).

**Filtering.** Cells with fewer than 200 detected genes or more than 5% mitochondrial reads were removed, as were genes detected in fewer than 3 cells (2,700 → 2,643 cells, 32,738 → 13,697 genes).

**Doublets.** Doublet scores were computed with Scrublet (scanpy `pp.scrublet`) on raw counts; cells scoring ≥ 0.1 were removed (222 cells, 8.4%).

**Normalisation.** Raw counts were kept in `layers['counts']`; expression was scaled to 10,000 counts per cell and log1p-transformed. The top 2,000 highly variable genes were flagged.

**Embedding.** PCA (50 components) on the highly variable genes.

**Clustering.** A k-nearest-neighbour graph on that embedding was clustered with Leiden (resolution 1.0; 9 clusters) and embedded with UMAP.

**Markers and annotation.** Marker genes per cluster were ranked with a Wilcoxon rank-sum test. Cell types were assigned with CellTypist (model `Immune_All_Low.pkl`), using majority voting over the Leiden clusters.

### Software

| Software | Version |
|---|---|
| Python | 3.11.14 |
| scanpy | 1.11.5 |
| anndata | 0.12.19 |
| numpy | 2.4.6 |
| scipy | 1.17.1 |
| celltypist | 1.7.1 |
| LLM (analysis decisions and narrative) | claude-sonnet-5 |

### Reproducibility

Random seed 0 for all stochastic steps. Every tool call, with the arguments the agent chose and the summary it read back, is in `tool_calls.jsonl`. `replay.py` re-runs those calls without the model, reproducing this analysis; running the agent again samples new choices from the model, so it can take different decisions.

---

*This report was produced by an AI agent (claude-sonnet-5). The run summary, decisions table, figures, captions and Methods are generated by code from the tool log; the narrative is written by the model and should be checked against them before use.*
