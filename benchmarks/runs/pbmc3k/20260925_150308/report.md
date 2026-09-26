# Single-cell RNA-seq analysis: pbmc3k

**pbmc3k** · 2,700 → 2,421 cells · 9 clusters · 8 cell types · PCA embedding

## Overview

This report covers the standard processing of a single-donor, single-run 10x PBMC dataset (`pbmc3k`) from raw counts through quality control, doublet removal, clustering, and cell-type annotation. The dataset started with 2,700 cells and 32,738 genes and ended with 2,421 cells and 13,697 genes annotated into 8 cell types. Because this is one donor profiled in one run with no batch structure and no experimental conditions to contrast, the analysis is a straightforward characterization rather than a comparison.

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

Gene identifiers were confirmed as human gene symbols with the standard `MT-` mitochondrial prefix (13 mitochondrial genes detected), so mitochondrial-fraction QC could proceed reliably. Per-cell metrics were unremarkable for healthy PBMCs: median genes/cell 817 (IQR 690–953), median UMIs/cell 2,197, and median mitochondrial fraction 2.03% (95th percentile 4.01%, max 22.6%).

I applied the tutorial-standard thresholds recommended for this data — minimum 200 genes/cell, maximum 5% mitochondrial content, and genes present in at least 3 cells — since the distributions showed no need to deviate from these defaults (essentially no cells failed the gene-count floor; only the high-mito tail, 2.1% of cells, was trimmed as likely dying/stressed cells). This removed 57 cells and 19,041 rarely-detected genes, leaving 2,643 cells and 13,697 genes for downstream analysis.

![QC distributions with cutoffs](figures/qc_thresholds.png)

*Per-cell QC before filtering (2,700 cells, median 2.03% mitochondrial). The standard 5% mitochondrial cutoff was kept and removes 57 cells (2.1%). Minimum genes per cell: 200, which removes 0 cells (fewest observed: 212).*

## Doublet detection

Scrublet was run on the full dataset as a single 10x run (no batch structure to split on). The doublet-score distribution was **not bimodal** (median 0.0442, max 0.506), so per the standard rule I used **median + 3×MAD** (0.1) rather than a bimodal valley. This is a first-pass 10x PBMC sample with no stated upstream genotype/hashing-based doublet removal, so the standard (non-permissive) rule was appropriate rather than the lighter-touch Scrublet-automatic threshold (0.227, which would have flagged far fewer cells). Filtering at 0.1 removed 222 cells (8.4% of the QC-passed population), consistent with typical doublet rates for this loading density.

![Doublet score distribution](figures/doublet_scores.png)

*Scrublet doublet scores for 2,643 cells; the distribution is not bimodal. Chosen threshold 0.1 (median + 3×MAD) flags 222 cells (8.4%). Scrublet's automatic threshold (dashed, 0.227) would flag 39.*

## Dimensionality reduction

This is a single clean sample from one donor and one 10x run, with no batch key and no condition to correct for. Standard PCA on the normalized, log-transformed, HVG-selected matrix (2,000 genes) is therefore the correct and simplest choice — batch-correction methods like scVI are unnecessary here and could only add noise. The first three PCs captured 10.3%, 3.5%, and 2.5% of variance respectively, consistent with a few dominant axes of immune-lineage variation (T/NK vs B vs myeloid) followed by more subtle sub-structure.

## Clustering

Leiden clustering at resolution 1.0 on the PCA embedding yielded 9 clusters, ranging from 7 to 550 cells. The size range spans from major lineages (T cells, monocytes) down to a very small, sharply distinct cluster of platelets — expected in PBMC preparations. 

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

CellTypist (fine-grained `Immune_All_Low.pkl` model, majority vote per cluster) assigned each cluster to a canonical PBMC cell type, and a coarse-grained second opinion (`Immune_All_High.pkl`) was used as a cross-check.

- **Cluster 0** — CCL5, NKG7, GZMA, CD8A/CD8B high → labeled 260 cells (10.7%), consistent with the top markers (cytotoxic CD8+ T cells).
- **Cluster 1** — CD79A, CD79B, MS4A1, HLA-DR genes high → 315 cells (13.0%), a clean, unambiguous call.
- **Clusters 2 and 3** — both CD3D+/CD3E+/IL7R+ and MS4A1/CD79A-negative on direct marker lookup, confirming genuine CD4 T cells rather than B cells; cluster 3 additionally shows elevated CCR7, consistent with a more naive phenotype within the same broad 1,065 cells (44.0%) category that CellTypist assigned to both clusters.
- **Cluster 4** — LYZ, S100A8/S100A9, FCN1 → 446 cells (18.4%).
- **Cluster 5** — CD74 and HLA-DR genes dominant, plus FCER1A → 41 cells (1.7%) (dendritic cells).
- **Cluster 6** — GNLY, NKG7, GZMB, KLRD1, FCGR3A strongly positive while CD3D/CD3E are absent (checked directly) → 140 cells (5.8%). This directly contradicts the coarse model's "ILC" call; since GNLY/NKG7/KLRD1/FCGR3A are all classic NK markers and T-cell receptor genes are absent, the fine-grained NK label was kept without relabeling.
- **Cluster 7** — LST1, FCGR3A, MS4A7, CD68 → 147 cells (6.1%).
- **Cluster 8** — PF4, PPBP, GNG11, SDPR → 7 cells (0.3%), a small but very sharply defined population.

All CellTypist calls were checked against cluster markers with `check_markers`, including genes that should be absent (e.g., MS4A1/CD79A in T-cell clusters, CD3D/CD3E in the NK cluster). Every label was well supported, so no cluster required relabeling. Final composition: 

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

*Canonical markers checked by the agent before accepting or changing labels (11 of 12 shown), ordered by the cell type each marks most, so the enriched dots run down the diagonal. A gene is shown if, in at least one type, padj < 0.05, log2FC > 1 and it is expressed in at least 25% of cells. Dot size: fraction of cells expressing; colour: mean expression scaled per gene. Not enriched in any type: NCAM1.*

![Cell-type composition](figures/composition.png)

*Cells per annotated type (2,421 cells, 8 types; CellTypist majority vote over Leiden clusters).*

![Marker genes per cell type](figures/marker_dotplot_celltype.png)

*Top 5 marker genes per cell type (Wilcoxon rank-sum). Dot size: fraction of cells expressing; colour: mean expression scaled per gene.*

## Caveats

- This is a single donor and a single 10x run, so all conclusions describe this individual's PBMC composition at one time point and should not be generalized without replication.
- Clusters 2 and 3 both received the same CellTypist label ("Tcm/Naive helper T cells") despite forming separate Leiden clusters; the marker evidence (CCR7 gradient) suggests they represent naive vs. more activated/memory CD4 T-cell states rather than two independent cell types — a finer distinction than the annotation model draws.
- The platelet cluster is very small (7 cells) and, while its markers (PF4, PPBP, GNG11) are unambiguous, cluster-level statistics for such a small group should be interpreted cautiously.
- Doublet filtering used a fixed statistical rule (median + 3×MAD) rather than a visually confirmed bimodal valley, since the score distribution did not show a clear valley; a small number of borderline calls near the threshold are possible in either direction.

## Conclusions

Standard QC and doublet filtering removed 10.3% of the original cells, yielding a clean dataset of 2,421 cells that resolved into 9 transcriptionally distinct clusters. These map onto the expected PBMC lineages — CD4 and CD8 T cells, B cells, CD16+ NK cells, classical and non-classical monocytes, dendritic cells, and a small platelet population — with marker expression matching canonical identities in every case, giving good confidence in the final annotation.

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
