# Single-cell RNA-seq analysis: pbmc3k

**pbmc3k** · 2,700 → 2,421 cells · 9 clusters · 8 cell types · PCA embedding

## Overview

This dataset comprises 2,700 droplets from a single healthy-donor PBMC sample (10x Genomics, one run). Since there is only one biological sample and no batch structure, the analysis follows the standard single-sample workflow: QC filtering, doublet removal, normalization, PCA-based dimensionality reduction, Leiden clustering, and CellTypist-based annotation. No condition comparison was requested or possible (single donor, single condition).

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

Gene symbols were confirmed as human with the standard `MT-` mitochondrial prefix (13 mitochondrial genes detected). Per-cell distributions were unremarkable for healthy PBMCs: median 817 genes/cell and median 2,197 total counts/cell, with mitochondrial fraction generally low (median 2.03%, 95th percentile 4.01%) but with a tail reaching 22.6%.

The recommended tutorial-standard thresholds (min 200 genes/cell, max 5% mitochondrial, genes present in ≥3 cells) were adopted as-is, since the projected impact was modest and consistent with removing genuinely stressed/dying cells rather than a broad swath of the data. This removed 57 cells (2.1%), essentially all for high mitochondrial content (no cells failed the gene-count floor), leaving 2,643 cells and 13,697 genes (many genes dropped were simply undetected outside the ≥3-cell floor).

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (2,700 cells, median 2.03% mitochondrial). The standard 5% mitochondrial cutoff was kept and removes 57 cells (2.1%). Minimum genes per cell: 200, which removes 0 cells (fewest observed: 212).*

## Doublet detection

Scrublet was run on the full dataset as a single 10x run (no batch key). The doublet-score distribution was unimodal (not bimodal), so per the standard rule I used median + 3×MAD (0.1) rather than a bimodal valley. This is a plain 10x PBMC run with no mention of upstream doublet removal (e.g. hashing or genotype demultiplexing), so the standard, stricter median+3×MAD rule was appropriate rather than the lighter-touch Scrublet-automatic threshold (0.227), which would have flagged far fewer cells (39) and likely left true doublets in place. Filtering at 0.1 removed 222 cells (8.4%), leaving 2,421 cells for downstream analysis.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 2,643 cells; the distribution is not bimodal. Chosen threshold 0.1 (median + 3×MAD) flags 222 cells (8.4%). Scrublet's automatic threshold (dashed, 0.227) would flag 39.*

## Dimensionality reduction

With a single clean batch (one donor, one run, no batch key), PCA on the 2,000 highly variable genes is the correct and simplest choice — scVI's batch-correction machinery is unnecessary here and would add noise without benefit. The top principal components captured a sensible amount of structure (PC1 10.3%, PC2 3.5%, PC3 2.5% of variance), consistent with a handful of major PBMC lineages driving most of the variation.

## Clustering

Leiden clustering (resolution 1.0) on the PCA embedding yielded 9 clusters, ranging from 7 to 550 cells. Cluster sizes and marker profiles map cleanly onto expected PBMC lineages (T cell subsets, B cells, monocyte subsets, NK cells, dendritic cells, and a small platelet/megakaryocyte cluster), as detailed below.

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

CellTypist (Immune_All_Low.pkl, fine-grained) assigned 8 cell types by per-cluster majority vote, cross-checked against the coarser Immune_All_High.pkl model and against canonical marker genes:

- **Cluster 0 — Tem/Trm cytotoxic T cells**: top markers CCL5, NKG7, GZMA, GZMK, CD8A/CD8B — classic cytotoxic CD8+ T cell signature.
- **Cluster 1 — B cells**: CD79A, CD79B, MS4A1, HLA-DR genes — canonical B cell markers, unambiguous.
- **Cluster 2 — Tcm/Naive helper T cells**: CD3D/CD3E and IL7R strongly enriched, CD8A absent — CD4+ T cell compartment.
- **Cluster 3 — Tcm/Naive helper T cells**: also CD3D/CD3E+, IL7R+, but with markedly higher CCR7 and SELL than cluster 2 and a ribosomal-gene-dominated profile — the more naive/quiescent end of the same CD4 T cell lineage, correctly grouped with cluster 2 under the same fine label.
- **Cluster 4 — Classical monocytes**: LYZ, S100A8/S100A9, FCN1, CST3 — canonical CD14+ monocyte signature.
- **Cluster 5 — DC**: CD74, HLA-DR series, FCER1A — dendritic cell profile, distinct from classical monocytes by lower S100A8/9.
- **Cluster 6 — CD16+ NK cells**: GNLY and NKG7 are the top two ranked genes genome-wide for this cluster, plus GZMB, PRF1, FCGR3A; CD3D, CD3E and CD79A are all absent (low % expressing, negative log2FC). The coarse model's alternative label "ILC" is a broader innate-lymphoid superset that includes NK cells rather than a genuine contradiction, so the fine-grained NK label was kept.
- **Cluster 7 — Non-classical monocytes**: LST1, FCER1G, AIF1, FCGR3A, MS4A7 — canonical CD16+ non-classical monocyte markers, distinct from cluster 4's classical monocyte profile.
- **Cluster 8 — Megakaryocytes/platelets**: PF4, PPBP, GNG11, SDPR — unambiguous platelet/megakaryocyte signature in this small cluster.

All labels were consistent with marker evidence, so no cluster relabeling was necessary. Final cell-type composition:

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

Notable populations: 1,065 cells (44.0%) CD4 T cells, 446 cells (18.4%) classical monocytes, 315 cells (13.0%) B cells, 260 cells (10.7%) cytotoxic T cells, 147 cells (6.1%) non-classical monocytes, 140 cells (5.8%) NK cells, 41 cells (1.7%) dendritic cells, and a small 7 cells (0.3%)-cell megakaryocyte/platelet population.

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (10 of 13 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: SELL, NCAM1, KLRB1.*

![Cell-type composition](figures/composition.png)

*Cells per annotated type (2,421 cells, 8 types; CellTypist majority vote over Leiden clusters).*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Caveats

- This is a single donor and single condition — no composition or differential-expression comparison across conditions was performed or appropriate here.
- Clusters 2 and 3 represent a continuum of CD4+ T cell states (naive through central-memory-like) rather than sharply distinct types; CellTypist's fine model assigns both the same label, and higher-resolution clustering could split them further if finer subtyping were needed.
- Doublet filtering used a standard, relatively strict statistical threshold (median+3×MAD); a small number of genuine transitional or high-RNA-content cells could in principle be lost, though the flagged fraction (8.4%) is in the expected range for this chemistry.
- Wilcoxon marker test p-values treat cells as independent observations and are not a substitute for biological replication; they were used here only for ranking and describing markers within one sample.

## Conclusions

Starting from 2,700 raw droplets, standard QC and doublet filtering retained 2,421 high-quality single cells (10.3% total removed) with 13,697 genes. PCA-based clustering at resolution 1.0 resolved 9 transcriptionally distinct clusters that map cleanly onto the expected major human PBMC lineages — CD4 and CD8 T cells, B cells, classical and non-classical monocytes, NK cells, dendritic cells, and megakaryocytes/platelets — each supported by canonical marker genes and largely concordant CellTypist calls across two reference granularities.

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

Random seed 0 for all stochastic steps. Every tool call, with the arguments the agent chose and the summary it read back, is in `tool_calls.jsonl`. The agent's choices are sampled from the model, so a re-run can take different decisions.

---

*This report was produced by an AI agent (claude-sonnet-5). The run summary, decisions table, figures, captions and Methods are generated by code from the tool log; the narrative is written by the model and should be checked against them before use.*
