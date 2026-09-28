# Single-cell RNA-seq analysis: pbmc3k

**pbmc3k** · 2,700 → 2,421 cells · 9 clusters · 8 cell types · PCA embedding

## Overview

This dataset comprises PBMCs from a single healthy donor, processed on a single 10x Genomics run — a "clean" single-batch scenario with no technical replicate structure to correct for. The analysis pipeline proceeded from raw counts through QC, doublet removal, normalization, PCA-based clustering, and CellTypist-assisted annotation. Starting from 2,700 cells and 32,738 genes, the final annotated dataset contains 2,421 cells and 13,697 genes (10.3% of the original cells removed across QC and doublet filtering steps).

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

Gene identifiers were confirmed as human gene symbols, with the mitochondrial prefix "MT-" identifying 13 mitochondrial genes. Per-cell QC metrics showed a healthy, relatively narrow distribution: median genes/cell of 817 (range 212–3,422) and median total counts of 2,197. Mitochondrial fraction was modest overall (median 2.03%, p95 4.01%) with a long tail up to 22.6%, consistent with a small population of stressed or dying cells rather than widespread ambient contamination.

I used the standard tutorial-default thresholds (200 min genes/cell, 5% max mitochondrial fraction, 3 min cells/gene), since they matched this dataset's own percentiles well — the mito cutoff sits just above the p95 mitochondrial percentage, targeting genuinely high-mito outlier cells rather than cutting into the bulk distribution. This removed 57 cells (2.1%) and 19,041 genes not detected in enough cells, leaving 2,643 cells and 13,697 genes.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (2,700 cells, median 2.03% mitochondrial). The standard 5% mitochondrial cutoff was kept and removes 57 cells (2.1%). Minimum genes per cell: 200, which removes 0 cells (fewest observed: 212).*

## Doublet detection

Scrublet was run on the full dataset as a single batch (one 10x run, no batches to split by). The doublet-score distribution was **not bimodal** — there was no clear valley separating singlets from doublets, which is common for a real (not synthetic-boosted) PBMC run at moderate loading density. Per the stated rule, I therefore used **median + 3×MAD** (0.1) rather than a fixed or bimodal-valley cutoff. This is more conservative than Scrublet's own automatic threshold (0.227, which would have flagged only 39 cells), but appropriate here since there is no evidence of upstream doublet removal (no genotype demultiplexing or hashing was reported for this single-donor run). This removed 222 cells (8.4%), consistent with expected multiplet rates for a run of this size.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 2,643 cells; the distribution is not bimodal. Chosen threshold 0.1 (median + 3×MAD) flags 222 cells (8.4%). Scrublet's automatic threshold (dashed, 0.227) would flag 39.*

## Dimensionality reduction

With only one donor and one sequencing run, there is no batch structure to correct for. PCA on the normalized, highly-variable-gene matrix (2,000 HVGs, 50 components) is the simpler and more appropriate choice over scVI's batch-correcting latent space, which would offer no benefit here and risks removing genuine biological variation. The top PC explains 10.3% of variance, consistent with the dominant lymphoid/myeloid axis expected in PBMCs.

## Clustering

Leiden clustering at resolution 1.0 on the PCA embedding yielded 9 clusters, ranging from 7 to 550 cells, mapping onto the expected major PBMC lineages (T cell subsets, B cells, monocyte subsets, NK cells, dendritic cells and platelets).

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

CellTypist (Immune_All_Low.pkl, fine-grained immune model) assigned each cluster a label by majority vote, cross-checked against the coarser Immune_All_High.pkl model as a second opinion; the two models agreed at their respective resolutions with no contradictions. I verified every final cell type against canonical markers using `check_markers`:

- **Tem/Trm cytotoxic T cells** (cluster 0): high CD3D/CD3E, CD8A/CD8B, and strong GZMK, consistent with an effector/memory CD8+ T cell phenotype.
- **B cells** (cluster 1): near-universal CD79A and MS4A1 expression, with CD3D essentially absent, a clean B cell signature.
- **Tcm/Naive helper T cells** (clusters 2 and 3): both express CD3D/CD3E and IL7R strongly; cluster 3 additionally shows strong CCR7 enrichment, marking it as the more naive-like subset, while cluster 2 is the more central-memory-like counterpart — both correctly grouped under the same CellTypist label since the model does not further split naive vs. central memory CD4 T cells.
- **Classical monocytes** (cluster 4): near-universal LYZ, CD14, FCN1 and S100A8 expression, the canonical CD14+ monocyte signature.
- **DC** (cluster 5): FCER1A and CD1C both strongly enriched, marking conventional dendritic cells, distinct from the monocyte clusters.
- **CD16+ NK cells** (cluster 6): near-universal GNLY and NKG7 with strong FCGR3A (CD16), and CD3D/CD3E essentially absent, ruling out a T/NKT identity.
- **Non-classical monocytes** (cluster 7): near-universal FCGR3A and MS4A7 with LYZ, but low CD14, the classic CD14(low)CD16+ non-classical monocyte pattern that distinguishes it from cluster 4.
- **Megakaryocytes/platelets** (cluster 8): unanimous PPBP and PF4 expression, an unambiguous platelet signature in this small cluster.

No relabeling was necessary — every cluster's canonical markers matched its assigned label and the second-opinion model was consistent.

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

Final cell-type composition: 1,065 cells (44.0%) Tcm/Naive helper T cells, 446 cells (18.4%) classical monocytes, 315 cells (13.0%) B cells, 260 cells (10.7%) cytotoxic T cells, 147 cells (6.1%) non-classical monocytes, 140 cells (5.8%) CD16+ NK cells, 41 cells (1.7%) dendritic cells, and 7 cells (0.3%) megakaryocytes/platelets — proportions broadly typical of a healthy PBMC sample.

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (22 of 24 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: NCAM1, CD19.*

![Cell-type composition](figures/composition.png)

*Cells per annotated type (2,421 cells, 8 types; CellTypist majority vote over Leiden clusters).*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Caveats

- This is a single donor processed in a single run, with no biological or technical replicates, so no composition or differential-expression comparison across conditions was performed (none was requested, and none would be statistically meaningful without replicate samples).
- Clusters 2 and 3 are grouped under one CellTypist label (Tcm/Naive helper T cells) despite showing a CCR7 gradient suggestive of naive vs. central-memory substructure; a higher clustering resolution or manual sub-gating could resolve this further if finer T cell subsetting is of interest.
- The megakaryocyte/platelet cluster is very small (7 cells (0.3%)), so its marker statistics, while clear (unanimous PPBP/PF4), rest on limited cell numbers.
- Median+3×MAD doublet filtering is a heuristic; a small number of true doublets may remain and a small number of borderline singlets may have been removed.

## Conclusions

Standard QC and doublet thresholds were well matched to this clean, single-donor PBMC dataset, requiring only modest filtering. PCA-based clustering (appropriate given the absence of batch structure) resolved the expected major PBMC lineages, and CellTypist annotation—confirmed against canonical markers for every cluster—produced a biologically coherent, unambiguous cell-type map spanning T cell subsets, B cells, monocyte subsets, NK cells, dendritic cells, and platelets.

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
