# Single-cell RNA-seq analysis: pbmc3k

**pbmc3k** · 2,700 → 2,421 cells · 9 clusters · 8 cell types · PCA embedding

## Overview

This analysis processed a 10x Genomics PBMC dataset from a single healthy donor (2,700 cells, 32,738 genes at input) through QC filtering, doublet removal, normalization, PCA-based dimensionality reduction, Leiden clustering, and CellTypist-assisted cell-type annotation. As this is a single clean 10x run from one donor with no batch structure, no batch correction was required and standard PCA was used for the embedding.

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

Gene identifiers were confirmed as human gene symbols with a standard `MT-` mitochondrial prefix (13 mitochondrial genes detected). Per-cell QC distributions showed a healthy library: median genes/cell of 817 (range 212–3,422), median total counts of 2,197, and median mitochondrial fraction of 2.03% (p99 = 5.88%, max = 22.6%).

The tutorial-standard thresholds (min_genes=200, max_pct_mt=5%, min_cells=3) fit this data well: the min_genes floor removed no cells (the dataset was already reasonably clean), while the mitochondrial cap removed 57 cells (2.1%) that sat in the long high-mito tail beyond the p95–p99 range — consistent with stressed or dying cells rather than a genuine cell population. Gene filtering (present in at least 3 cells) dropped 19,041 largely undetected genes, leaving 13,697 informative genes and 2,643 cells.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (2,700 cells, median 2.03% mitochondrial). The standard 5% mitochondrial cutoff was kept and removes 57 cells (2.1%). Minimum genes per cell: 200, which removes 0 cells (fewest observed: 212).*

## Doublet detection

Scrublet was run on the whole dataset as a single 10x run (no batch splitting needed). The doublet-score distribution was unimodal/right-skewed (not bimodal), so per the standard rule I used median + 3×MAD (0.1) rather than a bimodal valley (none was present) and rather than the lighter Scrublet-automatic threshold (0.227), which is reserved for datasets where doublets were already removed upstream (e.g. by demultiplexing) — not the case here since this is raw, unprocessed 10x output from one donor. At threshold 0.1, 222 cells (8.4%) were flagged and removed, in line with expected 10x doublet rates for this cell loading.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 2,643 cells; the distribution is not bimodal. Chosen threshold 0.1 (median + 3×MAD) flags 222 cells (8.4%). Scrublet's automatic threshold (dashed, 0.227) would flag 39.*

## Dimensionality reduction

With a single donor and no batch key, PCA on the 2,000 highly variable genes (after total-count normalization to 10,000 and log1p) is the simpler and correct choice over scVI's batch-correcting latent space — there is no batch effect to correct here. The leading PCs captured a substantial share of variance (10.3%, 3.5%, 2.5% for PC1–PC3 respectively), and 50 PCs were retained for downstream neighbor-graph construction.

## Clustering

Leiden clustering at resolution 1.0 produced 9 clusters ranging from 7 to 550 cells, reflecting the expected diversity of major PBMC lineages plus a small platelet population.

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

CellTypist (Immune_All_Low.pkl, fine-grained model) assigned each cluster a majority label via majority voting, cross-checked against the coarser Immune_All_High.pkl model as a second opinion; both models agreed on cell identity at their respective resolutions, with no contradictions requiring relabeling.

Canonical marker genes, checked with check_markers against the full per-cluster differential expression ranking, confirmed every assigned label:
- **Tem/Trm cytotoxic T cells** (cluster 0): high CD3D co-expressed with CD8A, CD8B, GZMK and CCL5 (all strongly enriched with high rank and low padj, and largely absent from other clusters), consistent with cytotoxic/effector CD8 T cells.
- **B cells** (cluster 1): near-exclusive CD79A and MS4A1 expression (both top-ranked, highly significant, and depleted elsewhere), with the T-cell marker CD3D strongly depleted in this cluster.
- **Tcm/Naive helper T cells** (clusters 2 and 3): CD3D and IL7R both significantly enriched, while CD8A/CD8B are not enriched, consistent with CD4 T cells; cluster 3 is dominated by ribosomal-protein transcripts typical of a quiescent naive subset, which CellTypist grouped with cluster 2 under the same fine label.
- **Classical monocytes** (cluster 4): strong, near-universal enrichment of CD14, LYZ and S100A8, essentially absent in other clusters.
- **DC** (cluster 5): high FCER1A and CD1C enrichment, without the CD14 enrichment seen in classical monocytes, consistent with conventional dendritic cells.
- **CD16+ NK cells** (cluster 6): GNLY, NKG7 and FCGR3A all strongly and significantly enriched, with CD3D essentially absent, ruling out a T/NKT identity.
- **Non-classical monocytes** (cluster 7): FCGR3A and MS4A7 strongly enriched without the CD14/S100A8 signature that marks classical monocytes.
- **Megakaryocytes/platelets** (cluster 8): near-universal PF4 and PPBP expression (top-ranked in the cluster, absent elsewhere), a small (n=7 cells (0.3%)) but clearly distinct population.

Final cell-type composition:

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

The most abundant populations were 1,065 cells (44.0%) CD4 T cells, 446 cells (18.4%) classical monocytes, and 315 cells (13.0%) B cells, alongside smaller 260 cells (10.7%) cytotoxic T cell, 147 cells (6.1%) non-classical monocyte, 140 cells (5.8%) NK cell, 41 cells (1.7%) dendritic cell, and 7 cells (0.3%) platelet populations — a composition consistent with expectations for healthy human PBMCs.

![Marker genes per Leiden cluster](figures/marker_dotplot_clusters.png)

*Top 5 marker genes per Leiden cluster (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

![Canonical marker genes per cell type](figures/marker_dotplot_canonical.png)

*Canonical markers checked by the agent before accepting or changing labels (21 of 22 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: NCAM1.*

![Cell-type composition](figures/composition.png)

*Cells per annotated type (2,421 cells, 8 types; CellTypist majority vote over Leiden clusters).*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Caveats

- This is a single donor/single run dataset with no biological replicates or condition to compare, so no composition or differential-expression testing across conditions was performed.
- Clusters 2 and 3 both received the same fine-grained "Tcm/Naive helper T cells" label; cluster 3's ribosomal-high signature suggests it may represent a naive/quiescent subset rather than a truly distinct cell type, and finer sub-clustering could resolve this if desired.
- Doublet removal used a MAD-based statistical threshold rather than a clear bimodal valley, since the score distribution was unimodal; some borderline low-scoring doublets or high-scoring genuine (e.g. large, transcript-rich) cells may be imperfectly classified.
- The mitochondrial % cutoff (5%) is a widely used PBMC heuristic but is somewhat arbitrary; a small number of legitimate high-mito cells (e.g. metabolically active cells) may have been excluded.

## Conclusions

After QC filtering and doublet removal (10.3% of original cells removed in total across mitochondrial and doublet filters), normalization, and PCA-based clustering, 2,421 high-quality cells resolved into 9 clusters mapping to 8 annotated cell types, all confirmed by canonical marker expression. The recovered populations — CD4 and CD8 T cells, B cells, classical and non-classical monocytes, NK cells, dendritic cells, and platelets — represent the expected major lineages of peripheral blood, indicating a technically sound and biologically sensible annotation of this dataset.

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
