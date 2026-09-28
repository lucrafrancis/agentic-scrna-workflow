# Single-cell RNA-seq analysis: pbmc_multibatch

**pbmc_multibatch** · 4,000 → 3,563 cells · 13 clusters · 12 cell types · scVI embedding

## Overview

This dataset combines PBMCs from two separate 10x runs (PBMC5k and PBMC10k) into one object of 4,000 cells x 15,792 genes, provided as raw counts with human gene symbols. Because the two runs are a genuine technical batch factor (different chip lanes/reagent lots run at different times), batch identity was treated as a covariate to correct for during integration rather than ignored.

## Key analysis decisions

Each row compares a standard default with what the agent chose after reading the data. Built from the tool log (`tool_calls.jsonl`), not written by the model, except the relabelling reasons, which quote the agent's tool call.

| Step | Standard default | What the data showed | Agent's choice | Effect |
|---|---|---|---|---|
| Mitochondrial cutoff | 5% | median 6.93%; default would remove 3,291 of 4,000 cells (82.3%) | **15%** (changed) | 73 cells removed by QC (1.8%) |
| Cell and gene floors | ≥200 genes/cell, genes in ≥3 cells | fewest genes in a cell: 264 | ≥200 genes/cell, genes in ≥3 cells (default kept) | 90 genes removed |
| Doublet threshold | Scrublet automatic (0.404 / 0.397 per batch) | scores not bimodal; median + 3×MAD 0.09 | **0.09** (median + 3×MAD) | 364 cells removed (9.27%) |
| Normalisation | 10,000 counts/cell, 2,000 HVGs | — | 10,000 counts/cell, 2,000 HVGs (default kept) | HVGs ranked within each batch ('batch') |
| Embedding | PCA | 2 batches in `obs['batch']` | **scVI**, correcting for batch | 10 latent dimensions, 2 batches |
| Leiden resolution | 1.0 | — | 1.0 (default kept) | 13 clusters |
| Annotation model | Immune_All_Low | — | Immune_All_Low (default kept) | 11 cell types |
| Cluster labels | CellTypist majority vote | Cluster 8 (24 cells) was labeled 'Tcm/Naive helper T cells' by CellTypist, but canonical T-cell markers are not enriched: CD3D (padj=1.0, 45.8% in vs 45.7% elsewhere), CD3E (padj=1.0), IL7R (padj=1.0), CCR7 (padj=1.0), and SELL is actually depleted (log2FC=-2.44, padj=0.016). Instead its top-ranked genes are MALAT1, JUN, and a long list of mitochondrial transcripts (MT-CO1, MT-ATP6, MT-ND5, MT-CYB, MT-ND4, MT-CO3, MT-CO2, MT-ND1) plus MTRNR2L12, a classic signature of stressed/damaged cells rather than a genuine lineage. No lineage marker (CD14, CD3, CD79A/MS4A1, GNLY) is enriched here. | cluster 8: Tcm/Naive helper T cells → **Low-quality/stressed cells** | 24 cells relabelled |

## Quality control

Mitochondrial genes were identified via the `MT-` prefix (13 genes). The per-cell QC distributions showed a notably right-shifted mitochondrial fraction for this dataset (median 6.93%, p95 12.2%, p99 16.8%, max 20%) compared to the tutorial-default cutoff of 5%. Applying that default would have discarded 3,291 cells (82.3% of the dataset) — clearly disproportionate, and consistent with a global shift in mito content rather than a distinct population of dying cells.

Instead, a max mitochondrial threshold of 15% was used, trimming only the extreme tail while retaining the bulk of the distribution. Combined with the standard minimum-genes-per-cell floor (200) — which removed no additional cells since the lowest cell already had 264 genes detected — and a minimum-cells-per-gene filter (3), this removed 73 cells (1.8%) and 90 lowly-detected genes, leaving 3,927 cells x 15,702 genes.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (4,000 cells, median 6.93% mitochondrial). The standard 5% mitochondrial cutoff (dashed) would remove 3,291 cells (82.3%); the chosen 15% cutoff removes 73 (1.8%). Minimum genes per cell: 200, which removes 0 cells (fewest observed: 264).*

## Doublet detection

Scrublet was run separately per 10x run (batch key `batch`) since doublet rates and score scales are run-specific. The combined score distribution was not bimodal, so the standard median + 3xMAD rule was used rather than a bimodal valley, giving a threshold of 0.09 (median 0.0384, max observed 0.465). This flagged and removed 364 cells (9.27%), leaving 3,563 cells. There was no indication in the task description that doublets had already been removed upstream (e.g. by hashing or genotype demultiplexing), so the standard rule — rather than the lighter-touch Scrublet-automatic threshold — was appropriate here.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 3,927 cells; the distribution is not bimodal. Chosen threshold 0.09 (median + 3×MAD) flags 364 cells (9.3%). The mean of Scrublet's per-batch automatic thresholds (dashed, 0.401) would flag 5.*

## Dimensionality reduction

The dataset has 2 batches from two independent 10x runs, a technical grouping rather than a biological condition. Since a real batch effect between separate sequencing runs is expected (loading, capture efficiency, ambient RNA differences), scVI was used with `batch` as the correction key to learn a joint latent embedding (10 latent dimensions, 400 training epochs on 2,000 HVGs) rather than plain PCA, which would leave run-driven variation uncorrected in the neighbor graph and could split shared cell types by batch. Raw counts were preserved throughout for downstream marker testing.

## Clustering

Leiden clustering on the scVI latent space (resolution 1.0) produced 13 clusters ranging from 23 to 712 cells. 

| Cluster | Cells | Top marker genes | Cell type | CellTypist label |
|---|---|---|---|---|
| 0 | 658 | S100A9, S100A8, S100A12, VCAN, S100A6 | Classical monocytes | Classical monocytes |
| 1 | 385 | CPVL, NEAT1, PSAP, CST3, FGL2 | Classical monocytes | Classical monocytes |
| 2 | 399 | GNLY, NKG7, PRF1, KLRD1, CST7 | CD16+ NK cells | CD16+ NK cells |
| 3 | 288 | CD79A, MS4A1, CD37, CD79B, HLA-DQA1 | Naive B cells | Naive B cells |
| 4 | 712 | RPS3A, RPL30, RPL32, RPL11, RPS27 | Tcm/Naive helper T cells | Tcm/Naive helper T cells |
| 5 | 172 | CD8B, RPS3A, RPS12, RPL32, RPS8 | Tcm/Naive cytotoxic T cells | Tcm/Naive cytotoxic T cells |
| 6 | 64 | HLA-DRB1, HLA-DPB1, HLA-DRA, HLA-DPA1, CST3 | DC2 | DC2 |
| 7 | 445 | IL32, TRAC, IL7R, LTB, LDHB | Tem/Effector helper T cells | Tem/Effector helper T cells |
| 8 | 24 | MALAT1, JUN, MTRNR2L12, MT-CO1, MT-ATP6 | **Low-quality/stressed cells** | Tcm/Naive helper T cells |
| 9 | 192 | CCL5, GZMA, CST7, NKG7, CD3D | Tem/Trm cytotoxic T cells | Tem/Trm cytotoxic T cells |
| 10 | 132 | KLRB1, GZMK, IL7R, DUSP2, NKG7 | MAIT cells | MAIT cells |
| 11 | 69 | LST1, AIF1, COTL1, FCGR3A, MS4A7 | Non-classical monocytes | Non-classical monocytes |
| 12 | 23 | IL3RA, CCDC50, PLD4, IRF8, TCF4 | pDC | pDC |

![UMAP](figures/umap.png)

*UMAP of 3,563 cells computed on the scVI latent space, coloured by cell type, Leiden cluster and batch (batch). scVI corrected for batch. Batches that overlap within each cell type indicate the integration worked.*

## Cell-type annotation

CellTypist (Immune_All_Low.pkl, majority vote per cluster) provided initial labels, cross-checked against the coarser Immune_All_High.pkl model and against canonical marker genes for every final cell type, using check_markers on the full DE ranking.

- **Classical monocytes** (clusters 0 and 1): CD14 and LYZ strongly and specifically enriched, with S100A8/S100A9 additionally enriched in cluster 0, and no FCGR3A/MS4A7 enrichment ruling out non-classical identity.
- **Non-classical monocytes**: FCGR3A and MS4A7 both strongly enriched, CD14 not distinctly enriched relative to other clusters, consistent with the CD14dim/CD16+ non-classical subset.
- **CD16+ NK cells**: GNLY, NKG7 and GZMB all highly and near-uniformly enriched, FCGR3A co-enriched, with CD3D/CD3E largely absent, ruling out a T/NK doublet population.
- **Naive B cells**: CD79A and MS4A1 both sharply enriched and essentially absent elsewhere.
- **DC2**: FCER1A and CD1C sharply and specifically enriched, distinguishing this small cluster from monocytes and pDCs.
- **pDC**: IL3RA, LILRA4 and GZMB all enriched to near-ubiquitous expression in this small cluster, a canonical pDC signature.
- **Tcm/Naive helper T cells** (cluster 4): CD3D/CD3E/IL7R enriched together with high CCR7 and SELL, and CD8A/CD8B essentially absent — a CD4 naive/central-memory phenotype (CD4 mRNA itself is not a reliable discriminator, as it is also expressed by monocytes).
- **Tcm/Naive cytotoxic T cells** (cluster 5): CD3D/CD3E enriched together with CD8A/CD8B and elevated CCR7/SELL — CD8 naive/central memory.
- **Tem/Effector helper T cells** (cluster 7): CD3D/CD3E/IL7R strongly enriched with low CCR7/SELL, and CD8A/CD8B absent — an effector/memory CD4 phenotype.
- **Tem/Trm cytotoxic T cells** (cluster 9): CCL5, GZMA/GZMK and CD8A/CD8B enriched together with CD3D/CD3E — a cytotoxic CD8 effector/memory phenotype.
- **MAIT cells**: KLRB1, GZMK and the MAIT-specific SLC4A10 all enriched, together with CD3E and CCL5.
- **Low-quality/stressed cells** (cluster 8, 24 cells (0.7%)): CellTypist called this cluster "Tcm/Naive helper T cells", but none of CD3D, CD3E, IL7R or CCR7 were significantly enriched, and SELL was actually depleted relative to the rest of the dataset. Its defining genes were instead MALAT1, JUN and a long run of mitochondrial transcripts (MT-CO1, MT-ATP6, MT-ND5, MT-CYB, MT-ND4, MT-CO3, MT-CO2, MT-ND1) plus MTRNR2L12 — a stress/damage signature with no lineage marker enriched. It was relabeled accordingly rather than reported as a T-cell subset.

Final cell-type composition: 

| Cell type | Cells | Share |
|---|---|---|
| Classical monocytes | 1,043 | 29.3% |
| Tcm/Naive helper T cells | 712 | 20.0% |
| Tem/Effector helper T cells | 445 | 12.5% |
| CD16+ NK cells | 399 | 11.2% |
| Naive B cells | 288 | 8.1% |
| Tem/Trm cytotoxic T cells | 192 | 5.4% |
| Tcm/Naive cytotoxic T cells | 172 | 4.8% |
| MAIT cells | 132 | 3.7% |
| Non-classical monocytes | 69 | 1.9% |
| DC2 | 64 | 1.8% |
| Low-quality/stressed cells | 24 | 0.7% |
| pDC | 23 | 0.6% |
| **Total** | **3,563** | |

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (27 of 28 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: MT-CO1.*

![Cell-type composition](figures/composition.png)

*Cells per annotated type (3,563 cells, 12 types; CellTypist majority vote over Leiden clusters, with clusters relabelled by the agent from their markers (see the decisions table)).*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Caveats

- Cluster 8 (24 cells (0.7%)) reflects damaged/stressed cells rather than a discrete biological population and should be excluded from downstream biological interpretation; it is retained in the annotated object for transparency.
- The two clusters both labeled "Classical monocytes" (clusters 0 and 1) differ somewhat in markers (cluster 1 shows higher HLA-DR/CST3/CD68), which may reflect an activation or maturity gradient within classical monocytes rather than a distinct subset — worth revisiting with finer-resolution clustering if a monocyte-focused question arises.
- No biological condition was present in this dataset (only the two technical 10x runs), so no composition or differential-expression comparison between conditions was performed; this analysis focused on quality control, batch-aware integration and cell-type annotation only.
- The elevated mitochondrial fraction across this dataset (well above the tutorial-default cutoff of 5%) means a permissive filtering threshold was used; cells at the upper end of the retained range should be interpreted with some caution as possibly reflecting somewhat lower RNA quality rather than a wholly clean cell population.

## Conclusions

Starting from 4,000 raw cells across two PBMC 10x runs, QC and doublet filtering tuned to this dataset's actual distributions (rather than default tutorial cutoffs) retained 3,563 high-quality, singlet cells. Batch-aware integration with scVI resolved the two runs into a shared embedding, yielding 13 Leiden clusters that mapped onto 12 expected PBMC cell types (monocyte subsets, B cells, NK cells, CD4/CD8 T cell subsets including MAIT cells, DC2, and pDCs), each supported by canonical marker evidence, with one small artefactual cluster of stressed cells correctly separated out from genuine lineages.

## Methods

*Generated from the tool log: these are the steps and parameters that actually ran.*

**Quality control.** Per-cell metrics were computed with scanpy `calculate_qc_metrics`; mitochondrial genes were those prefixed `MT-` (13 found).

**Filtering.** Cells with fewer than 200 detected genes or more than 15% mitochondrial reads were removed, as were genes detected in fewer than 3 cells (4,000 → 3,927 cells, 15,792 → 15,702 genes).

**Doublets.** Doublet scores were computed with Scrublet (scanpy `pp.scrublet`), run separately within each batch ('batch'), on raw counts; cells scoring ≥ 0.09 were removed (364 cells, 9.27%).

**Normalisation.** Raw counts were kept in `layers['counts']`; expression was scaled to 10,000 counts per cell and log1p-transformed. The top 2,000 highly variable genes were flagged, ranked within each batch ('batch').

**Embedding.** An scVI model (10 latent dimensions, 400 epochs) was trained on raw counts of 2,000 highly variable genes, correcting for batch; its latent space was used for clustering, UMAP and annotation (raw counts were unchanged).

**Clustering.** A k-nearest-neighbour graph on that embedding was clustered with Leiden (resolution 1.0; 13 clusters) and embedded with UMAP.

**Markers and annotation.** Marker genes per cluster were ranked with a Wilcoxon rank-sum test. Cell types were assigned with CellTypist (model `Immune_All_Low.pkl`), using majority voting over the Leiden clusters.

**Label curation.** Where marker genes contradicted CellTypist, the agent relabelled whole clusters (1 clusters; listed in the decisions table with the evidence). The original CellTypist labels are kept in `obs['cell_type_celltypist']`.

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
| LLM (analysis decisions and narrative) | claude-sonnet-5 |

### Reproducibility

Random seed 0 for all stochastic steps. Every tool call, with the arguments the agent chose and the summary it read back, is in `tool_calls.jsonl`. `replay.py` re-runs those calls without the model, reproducing this analysis; running the agent again samples new choices from the model, so it can take different decisions.

---

*This report was produced by an AI agent (claude-sonnet-5). The run summary, decisions table, figures, captions and Methods are generated by code from the tool log; the narrative is written by the model and should be checked against them before use.*
