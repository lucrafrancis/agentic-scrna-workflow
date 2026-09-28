# Single-cell RNA-seq analysis: pbmc3k

**pbmc3k** · 2,700 → 2,421 cells · 9 clusters · 8 cell types · PCA embedding

## Overview

This report covers a standard single-sample PBMC analysis: one healthy donor, one 10x Genomics run, 2,700 cells and 32,738 genes as input. No batch key was detected, consistent with the description of a single run, so no batch integration was needed at any stage. The workflow proceeded from raw counts through QC filtering, doublet removal, normalization, PCA-based dimensionality reduction, Leiden clustering, and CellTypist-assisted cell-type annotation, with every canonical marker checked before finalizing labels.

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

Gene identifiers were confirmed as human gene symbols, with the mitochondrial gene prefix `MT-` detected (13 mitochondrial genes found). Per-cell QC metrics were unremarkable for healthy PBMCs: median genes/cell 817, median UMI counts 2,197, median mitochondrial fraction 2.03% (up to 22.6% in a tail of outlier cells).

I applied the tutorial-standard thresholds recommended for this dataset: minimum 200 genes/cell, maximum 5% mitochondrial content, and genes present in at least 3 cells. These are well suited here because the gene-count distribution had essentially no low-count tail (min was already 212 genes), so the min_genes cutoff removed nothing; the mitochondrial cap did the actual work, removing 57 likely dying/stressed cells (2.1% of the input). Gene filtering dropped 19,041 genes detected in too few cells, mostly uninformative low-expression transcripts. After this step: 2,643 cells and 13,697 genes remained.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (2,700 cells, median 2.03% mitochondrial). The standard 5% mitochondrial cutoff was kept and removes 57 cells (2.1%). Minimum genes per cell: 200, which removes 0 cells (fewest observed: 212).*

## Doublet detection

Scrublet was run on this single 10x run as one batch (no per-run splitting needed). The doublet-score distribution was unimodal with a long right tail rather than two clearly separated peaks (not bimodal), so per the decision rule I used median + 3×MAD (0.1) rather than a bimodal valley. There is no indication doublets were removed upstream (e.g. via genotype demultiplexing), so the standard rule — not the lighter-touch Scrublet-automatic threshold — was the appropriate choice. This flagged and removed 222 cells (8.4%), a plausible heterotypic doublet rate for this loading density. After doublet removal: 2,421 cells remained for downstream analysis.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 2,643 cells; the distribution is not bimodal. Chosen threshold 0.1 (median + 3×MAD) flags 222 cells (8.4%). Scrublet's automatic threshold (dashed, 0.227) would flag 39.*

## Dimensionality reduction

Since this dataset is a single clean batch from one donor with no technical grouping variable, PCA on the normalized, log-transformed HVG matrix (top 2,000 genes, target sum 10,000) is the correct and simplest choice — scVI's batch-correction machinery is unnecessary overhead here and could over-smooth genuine biological structure in a single sample. PCA was run with 50 components; the leading PCs (PC1 10.3%, PC2 3.5%, PC3 2.5% of variance) show a clear elbow typical of well-structured immune cell data.

## Clustering

Leiden clustering (resolution 1.0) on the PCA representation yielded 9 clusters, ranging in size from 7 to 550 cells. Cluster sizes and UMAP topology are consistent with the expected composition of PBMCs: a few large lymphocyte and monocyte populations plus small, distinct populations (dendritic cells, platelets).

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

CellTypist (Immune_All_Low.pkl, fine-grained) was used for majority-vote annotation per cluster, cross-checked against the coarser Immune_All_High.pkl model as a second opinion. The two models agreed at the appropriate granularity for every cluster (e.g. cluster 6 was called "CD16+ NK cells" by the fine model and the broader "ILC" category by the coarse model — not a contradiction, since NK cells are ILCs).

Before accepting labels, canonical markers were checked directly against the full marker ranking for every final cell type (not just ambiguous ones):

- **Cytotoxic T cells** (cluster 0): CCL5 was the single top-ranked marker gene for this cluster, expressed in nearly all of its cells and strongly depleted elsewhere, together with high CD8A/CD8B and GZMK, and CD3D/CD3E positivity confirming a T-lineage identity co-expressed with cytotoxic granzymes.
- **B cells** (cluster 1): CD79A and MS4A1 were both highly specific — expressed in the large majority of cluster cells and rare outside it — with CD3D strongly depleted — a clean B-cell signature.
- **CD4 T cells / Tcm-Naive helper T cells** (clusters 2 and 3): both show strong CD3D/CD3E/IL7R positivity and CD8A/CD79A absence. Cluster 3 shows markedly higher CCR7 than cluster 2, suggesting a more naive-like subpopulation within the same broad CellTypist label; both were kept under the same annotation since CellTypist's own majority vote did not distinguish them and the marker signal is a difference of degree, not identity.
- **Classical monocytes** (cluster 4): LYZ was expressed in essentially every cell of this cluster and CD14 in the majority, both strongly enriched versus every other cluster — a textbook classical monocyte signature.
- **Dendritic cells** (cluster 5): FCER1A and CD1C were both highly specific to this cluster and nearly absent elsewhere, confirming a conventional DC identity distinct from monocytes despite shared CD74/HLA-DR expression.
- **CD16+ NK cells** (cluster 6): GNLY, NKG7 and GZMB were the top-ranked genes for this cluster, expressed in nearly all its cells; FCGR3A was also strongly enriched, and CD3D was essentially absent — ruling out a T-cell identity and confirming NK cells.
- **Non-classical monocytes** (cluster 7): FCGR3A and MS4A7 were both highly specific and strongly enriched, while CD14 was low and not significantly enriched — the expected classical-vs-non-classical monocyte split confirmed by reciprocal CD14/FCGR3A patterns between clusters 4 and 7.
- **Megakaryocytes/platelets** (cluster 8): PPBP, PF4 and ITGA2B were essentially exclusive to this small cluster, confirming the platelet/megakaryocyte identity despite the very small cluster size (7 cells (0.3%)).

No cluster's marker evidence contradicted its CellTypist label, so no relabeling was necessary.

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

Final annotated dataset: 2,421 cells across 8 cell types, following removal of 279 cells (10.3% of the original input) through QC and doublet filtering combined.

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (24 of 26 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: CD19, NCAM1.*

![Cell-type composition](figures/composition.png)

*Cells per annotated type (2,421 cells, 8 types; CellTypist majority vote over Leiden clusters).*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Caveats

- This is a single donor/single run, so no batch-correction or condition comparison was performed or appropriate; all conclusions describe this one sample's composition and cannot be generalized without replicates.
- Clusters 2 and 3 (both CD4 T cells) likely represent a naive/memory continuum rather than two distinct cell types; the marker difference (CCR7) is one of degree and would benefit from a targeted resolution increase or additional markers (e.g. CCR7/SELL/CD27 panels) if finer T-cell substates are needed.
- The dendritic cell and megakaryocyte/platelet clusters are small (41 cells (1.7%) and 7 cells (0.3%) respectively), so their marker statistics and any downstream inference about them carry wider uncertainty than the larger populations.
- Doublet removal used an unsupervised score threshold (median + 3×MAD); while principled, some borderline transcriptionally-mixed cells (e.g. platelet-adhered lymphocytes, a known PBMC artifact) may still remain or have been removed unnecessarily.

## Conclusions

Standard QC and doublet-filtering thresholds were appropriate for this clean, single-sample PBMC dataset, removing a modest fraction of low-quality and doublet cells (10.3% total) without needing batch correction. PCA-based clustering at resolution 1.0 resolved the expected major PBMC lineages — cytotoxic and helper T cells, B cells, classical and non-classical monocytes, dendritic cells, NK cells, and a small platelet/megakaryocyte population — each supported by canonical markers confirmed directly against the data, giving confidence in the final annotation.

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
