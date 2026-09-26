# Single-cell RNA-seq analysis: pbmc3k

**pbmc3k** · 2,700 → 2,421 cells · 9 clusters · 8 cell types · PCA embedding

## Overview

This analysis takes a single healthy-donor PBMC 10x dataset (2,700 cells, 32,738 genes) from raw counts to an annotated cell atlas. Because this is one donor processed in a single 10x run, no batch correction was required — dimensionality reduction and clustering could operate directly on PCA of the normalized, log-transformed data.

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

Gene identifiers were confirmed as human gene symbols with the standard "MT-" mitochondrial prefix (13 mitochondrial genes detected), enabling accurate mitochondrial-content QC. Per-cell distributions were unremarkable for healthy PBMCs: genes/cell ranged from 212 to 3,422 (median 817), and mitochondrial fraction was low for most cells (median 2.03%) with a long tail up to 22.6%, marking a subset of stressed/dying cells.

The tutorial-default thresholds (min_genes=200, max_pct_mt=5%, min_cells=3) fit this dataset well: since the lowest cell already had 212 genes, the min_genes floor removed no cells, while the mitochondrial cap removed 57 cells (2.1%) from the high-mito tail — exactly the stressed/dying cells that threshold is meant to catch. After filtering, 2,643 cells and 13,697 genes remained (57 cells and 19,041 low-count genes dropped).

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (2,700 cells, median 2.03% mitochondrial). The standard 5% mitochondrial cutoff was kept and removes 57 cells (2.1%). Minimum genes per cell: 200, which removes 0 cells (fewest observed: 212).*

## Doublet detection

Scrublet was run on the whole dataset as a single 10x run (no batch splitting needed). The doublet-score distribution was **not bimodal** (not bimodal, median 0.0442, max 0.506), so per the stated rule I used **median + 3×MAD** (0.1) rather than a visual valley. There is no indication in the dataset description that doublets were already removed upstream (e.g. via genotype demultiplexing or hashing), so the standard median+3MAD threshold — not a lighter touch — was appropriate. This removed 222 cells (8.4%), leaving 2,421 cells for downstream analysis. This rate is in the typical range for a standard, unhashed 10x PBMC run.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 2,643 cells; the distribution is not bimodal. Chosen threshold 0.1 (median + 3×MAD) flags 222 cells (8.4%). Scrublet's automatic threshold (dashed, 0.227) would flag 39.*

## Dimensionality reduction

With a single clean batch and no biologically or technically distinct sample groups, standard **PCA** on the normalized, HVG-selected matrix (2,000 HVGs, target_sum=10,000) was the correct and simplest choice — scVI's batch-correction machinery would add complexity without benefit here. The top PCs captured the expected structure for PBMCs (PC1 10.3% of variance, PC2 3.5%, PC3 2.5%), consistent with the major axes being lymphoid-vs-myeloid identity and finer lineage distinctions.

## Clustering

Leiden clustering at resolution 1.0 on the PCA embedding produced 9 clusters, ranging from a large T-cell cluster (550 cells) down to a small, well-separated platelet cluster (7 cells) — matching the expected composition of PBMCs, where lymphocytes dominate and megakaryocyte/platelet contamination is a small but recognizable minority population.

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

CellTypist (Immune_All_Low.pkl, fine-grained model) assigned 8 cell types via majority vote per cluster, and I cross-checked every cluster's markers (both the ranked marker list and explicit lookups of canonical genes) against this label and against the coarse-grained Immune_All_High.pkl model run as a second opinion:

- **Cluster 0 → 260 cells (10.7%)**: CD3D/CD3E positive together with CD8A/CD8B, high NKG7/GZMA/GZMK, and low CCR7/SELL — an activated/effector, not naive, CD8 T-cell profile. Confirmed.
- **Cluster 1 → 315 cells (13.0%)**: CD79A strongly enriched and essentially absent elsewhere, together with MS4A1 and HLA-DR genes, while CD3D is negative. Confirmed.
- **Clusters 2 & 3 → 1,065 cells (44.0%)**: both clusters are CD3D/CD3E/IL7R positive and CD8-negative — genuine CD4 T cells. Cluster 3 shows relatively higher CCR7/SELL (more naive-like, and its ribosomal-gene-dominated marker list is consistent with the low transcriptional/metabolic activity typical of naive T cells) while cluster 2 shows relatively higher IL7R (more central-memory-like). CellTypist's low-resolution model does not separate these substates further, and the underlying markers do not contradict grouping them together as CD4 T helper cells.
- **Cluster 4 → 446 cells (18.4%)**: CD14 strongly enriched in this cluster and nearly absent elsewhere, with LYZ/S100A8/S100A9/FCN1 all top-ranked markers. Confirmed.
- **Cluster 5 → 41 cells (1.7%)**: FCER1A strongly enriched in this cluster and nearly absent elsewhere, together with high CST3 and HLA-DR genes — a classic myeloid dendritic cell signature. Confirmed.
- **Cluster 6 → 140 cells (5.8%)**: GNLY and NKG7 both essentially ubiquitous in this cluster, KLRD1 clearly enriched, and CD3D/CD3E negative. The coarse model's alternative label "ILC" is a broader lineage umbrella that includes NK cells, not a contradiction. Confirmed.
- **Cluster 7 → 147 cells (6.1%)**: FCGR3A and MS4A7 both strongly enriched while CD14 stays low — the canonical CD14-low/CD16-high non-classical monocyte pattern. Confirmed.
- **Cluster 8 → 7 cells (0.3%)**: PPBP and PF4 both essentially ubiquitous in this cluster and essentially absent elsewhere — an unambiguous, small platelet/megakaryocyte-derived population.

No cluster required relabeling; CellTypist's calls were consistently corroborated by canonical marker expression and absence.

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

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (19 of 20 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: SELL.*

![Cell-type composition](figures/composition.png)

*Cells per annotated type (2,421 cells, 8 types; CellTypist majority vote over Leiden clusters).*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Caveats

- This is a single donor and a single 10x run, so results describe this individual's PBMC composition and cannot be generalized without replication; no batch effects needed correcting, but also no biological replicate structure is available to assess donor-to-donor variability.
- Clusters 2 and 3 likely represent a naive/memory CD4 T-cell continuum rather than two sharply distinct populations; the boundary drawn by Leiden at this resolution should be interpreted as a coarse split, not a hard biological boundary.
- The megakaryocyte/platelet cluster is very small (7 cells (0.3%)); it is a well-supported population by markers but statistical power for anything beyond identity confirmation is limited.
- Doublet removal used an unsupervised score threshold (median+3×MAD); as with any such cutoff, some true doublets near the boundary may remain and a few unusual-but-real cells may have been removed.

## Conclusions

Starting from 2,700 raw cells, standard QC and doublet filtering yielded a clean set of 2,421 cells (10.3% removed overall) that resolved into 9 Leiden clusters representing all major expected PBMC lineages: cytotoxic and helper T cells, B cells, classical and non-classical monocytes, dendritic cells, NK cells, and a small platelet/megakaryocyte population. Cell-type calls were confirmed against canonical marker genes and cross-validated between two CellTypist model resolutions, giving a confident, fully annotated PBMC reference dataset ready for downstream use.

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
