"""Analysis tools.

THE CONTRACT (see CLAUDE.md):
  Every tool is a deterministic Python function that does real work AND returns a
  structured, JSON-serializable summary dict. That dict is the *only* thing the LLM
  sees — it decides the next step entirely from it. A tool that returns None is a bug.

Tools operate on a shared AnnData held by the running session (see loop.py). Each tool
validates the state it expects and returns a clear error summary (not a stack trace) the
agent can react to.

All functions below are stubs — signatures and return contracts only. Implementations
land in the next step.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import scanpy as sc
from scipy import sparse
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde

from agent import config
from agent.session import SESSION

Summary = dict[str, Any]

# obs columns that commonly denote a batch/sample grouping, checked by inspect_dataset.
_BATCH_HINTS = {"batch", "sample", "donor", "patient", "subject", "condition", "dataset"}

# Quantiles reported for every QC metric distribution.
_QUANTILES = {"min": 0.0, "p25": 0.25, "median": 0.5, "p75": 0.75, "p95": 0.95, "p99": 0.99, "max": 1.0}


def _quantile_summary(values) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    return {name: round(float(np.quantile(arr, q)), 4) for name, q in _QUANTILES.items()}


def _looks_like_counts(X) -> bool:
    """Heuristic: True if X appears to be non-negative integer counts (a sample of values)."""
    data = X.data[:1000] if sparse.issparse(X) else np.asarray(X).ravel()[:1000]
    if data.size == 0:
        return True
    return bool(np.all(data >= 0) and np.allclose(data, np.round(data)))


def inspect_dataset() -> Summary:
    """n cells/genes, obs/var columns, candidate batch key, and count/layer state.

    Non-mutating; the agent calls this first to orient itself. Reports enough for the
    dimensionality-reduction decision (batch key) and the counts guardrail (has_raw_counts).
    """
    adata = SESSION.require_adata()
    obs_cols = list(adata.obs.columns)

    candidate_batch_key, n_batches = None, None
    for col in obs_cols:
        if col.lower() in _BATCH_HINTS and adata.obs[col].nunique() > 1:
            candidate_batch_key = col
            n_batches = int(adata.obs[col].nunique())
            break

    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": obs_cols,
        "var_columns": list(adata.var.columns),
        "candidate_batch_key": candidate_batch_key,
        "n_batches": n_batches,
        "has_raw_counts": _looks_like_counts(adata.X),
        "existing_layers": list(adata.layers.keys()),
        "gene_id_sample": [str(g) for g in adata.var_names[:5]],
        "duplicate_barcodes_fixed": SESSION.n_duplicate_barcodes,
    }


def check_gene_identifiers() -> Summary:
    """Detect the gene-identifier format and prepare for symbol-dependent steps.

    Runs before compute_qc: classifies var_names as symbol / ensembl / other, guesses
    organism and the mitochondrial prefix, and makes var_names unique. Stores the format
    and mito prefix on the session for compute_qc to use. Non-mutating (var cleanup only).
    """
    adata = SESSION.require_adata()
    adata.var_names_make_unique()
    names = adata.var_names.astype(str)
    n = len(names)

    ensembl = sum(bool(re.match(r"ENS[A-Z]*G\d+", s.upper())) for s in names)
    symbolic = sum(bool(re.match(r"[A-Za-z][A-Za-z0-9._-]*$", s)) for s in names)
    if n and ensembl / n > 0.5:
        gene_format = "ensembl"
    elif n and symbolic / n > 0.5:
        gene_format = "symbol"
    else:
        gene_format = "other"

    # Organism / mito prefix only meaningful when the index is symbols.
    organism, mito_prefix = None, None
    if gene_format == "symbol":
        if any(s.startswith("MT-") for s in names):
            organism, mito_prefix = "human", "MT-"
        elif any(s.startswith("mt-") for s in names):
            organism, mito_prefix = "mouse", "mt-"
        else:
            mito_prefix = "MT-"  # assume human symbols; compute_qc reports if 0 are found

    symbols_available = gene_format == "symbol"
    SESSION.gene_format = gene_format
    SESSION.mito_prefix = mito_prefix

    return {
        "gene_format": gene_format,
        "organism_guess": organism,
        "mito_prefix": mito_prefix,
        "symbols_available": symbols_available,
        "n_genes": int(n),
        "gene_id_sample": [str(g) for g in names[:5]],
        # Option A: flag rather than map. Symbol-dependent steps are degraded, not fixed.
        "note": None
        if symbols_available
        else "No gene symbols in var_names; mito detection and marker interpretation are degraded.",
    }


def compute_qc() -> Summary:
    """Compute per-cell QC (n_genes, total_counts, pct_mito) and return distributions.

    Requires check_gene_identifiers to have run (needs the mito prefix). Annotates obs/var
    with QC metrics; does not remove anything, so it does not checkpoint.
    """
    adata = SESSION.require_adata()
    if SESSION.gene_format is None:
        return {
            "error": "gene_identifiers_not_checked",
            "message": "Call check_gene_identifiers before compute_qc.",
        }

    prefix = SESSION.mito_prefix
    adata.var["mt"] = adata.var_names.str.startswith(prefix) if prefix else False
    n_mito = int(np.asarray(adata.var["mt"]).sum())
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    return {
        "n_mito_genes_found": n_mito,
        "n_genes_by_counts": _quantile_summary(adata.obs["n_genes_by_counts"]),
        "total_counts": _quantile_summary(adata.obs["total_counts"]),
        "pct_counts_mt": _quantile_summary(adata.obs["pct_counts_mt"]),
    }


def recommend_qc_thresholds() -> Summary:
    """Suggest tutorial-default thresholds and their projected impact on THIS dataset.

    Suggestion only, no mutation. The agent decides whether to accept or adjust before
    calling filter_cells_and_genes.
    """
    adata = SESSION.require_adata()
    if "n_genes_by_counts" not in adata.obs:
        return {"error": "qc_not_computed", "message": "Call compute_qc before recommend_qc_thresholds."}

    min_genes, min_cells, max_pct_mt = 200, 3, 5.0
    below_min_genes = int((adata.obs["n_genes_by_counts"] < min_genes).sum())
    above_max_mt = int((adata.obs["pct_counts_mt"] > max_pct_mt).sum())
    keep = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["pct_counts_mt"] <= max_pct_mt)
    genes_below_min_cells = int((adata.var["n_cells_by_counts"] < min_cells).sum())

    return {
        "recommended": {"min_genes": min_genes, "max_pct_mt": max_pct_mt, "min_cells": min_cells},
        "rationale": {
            "min_genes": "standard scanpy tutorial floor for viable cells",
            "max_pct_mt": "PBMC-typical mitochondrial cap; high mito % marks stressed/dying cells",
            "min_cells": "drop genes detected in too few cells (uninformative)",
        },
        "projected_impact": {
            "cells_total": int(adata.n_obs),
            "cells_removed_min_genes": below_min_genes,
            "cells_removed_max_pct_mt": above_max_mt,
            "cells_removed_total": int((~keep).sum()),
            "genes_total": int(adata.n_vars),
            "genes_removed_min_cells": genes_below_min_cells,
        },
    }


def filter_cells_and_genes(min_genes: int = 200, max_pct_mt: float = 5.0, min_cells: int = 3) -> Summary:
    """Apply cell and gene filters, then checkpoint. First mutating tool in the arc."""
    adata = SESSION.require_adata()
    if "n_genes_by_counts" not in adata.obs:
        return {"error": "qc_not_computed", "message": "Call compute_qc before filtering."}

    before = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mt].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    SESSION.adata = adata

    after = {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}
    checkpoint = SESSION.checkpoint("after_filter")
    return {
        "applied": {"min_genes": min_genes, "max_pct_mt": max_pct_mt, "min_cells": min_cells},
        "before": before,
        "after": after,
        "removed": {
            "cells": before["n_cells"] - after["n_cells"],
            "genes": before["n_genes"] - after["n_genes"],
        },
        "checkpoint": checkpoint,
    }


def _valley_threshold(scores: np.ndarray) -> tuple[float | None, bool]:
    """If the score distribution is bimodal, return the density valley between the two
    dominant modes and True; otherwise (None, False). Used as the preferred doublet cutoff.
    """
    if len(np.unique(scores)) < 10:
        return None, False
    grid = np.linspace(float(scores.min()), float(scores.max()), 256)
    density = gaussian_kde(scores)(grid)
    peaks, _ = find_peaks(density, height=density.max() * 0.05)
    if len(peaks) < 2:
        return None, False
    lo, hi = sorted(peaks[np.argsort(density[peaks])[-2:]])  # the two tallest peaks
    valley = lo + int(np.argmin(density[lo : hi + 1]))
    return float(grid[valley]), True


def detect_doublets() -> Summary:
    """Run Scrublet and report the score distribution plus candidate thresholds.

    Non-mutating: stores per-cell doublet scores in obs but removes nothing. The agent
    inspects the distribution, chooses a threshold, and calls filter_doublets. Recommends
    the bimodal-valley cutoff when the distribution is bimodal, else median + 3*MAD.
    """
    adata = SESSION.require_adata()
    sc.pp.scrublet(adata, random_state=config.SEED)
    scores = np.asarray(adata.obs["doublet_score"], dtype=float)

    median = float(np.median(scores))
    mad = float(np.median(np.abs(scores - median)))  # raw (unscaled) MAD, per project spec
    mad_threshold = median + 3 * mad
    valley_threshold, is_bimodal = _valley_threshold(scores)
    scrublet_auto = adata.uns.get("scrublet", {}).get("threshold")
    scrublet_auto = float(scrublet_auto) if scrublet_auto is not None else None

    recommended = valley_threshold if is_bimodal else mad_threshold
    rule = "bimodal_valley" if is_bimodal else "median+3*MAD"
    n_flagged = int((scores >= recommended).sum())
    counts, edges = np.histogram(scores, bins=20)

    return {
        "is_bimodal": is_bimodal,
        "score_distribution": _quantile_summary(scores),
        "histogram": {"counts": counts.tolist(), "bin_edges": [round(e, 4) for e in edges]},
        "candidate_thresholds": {
            "bimodal_valley": valley_threshold,
            "median_3mad": round(mad_threshold, 4),
            "scrublet_auto": scrublet_auto,
        },
        "recommended_threshold": round(float(recommended), 4),
        "recommended_rule": rule,
        "n_flagged_at_recommended": n_flagged,
        "pct_flagged_at_recommended": round(100 * n_flagged / len(scores), 2),
    }


def filter_doublets(threshold: float) -> Summary:
    """Remove cells with doublet_score >= threshold, then checkpoint. Mutating.

    The agent supplies the threshold it chose from detect_doublets' candidates.
    """
    adata = SESSION.require_adata()
    if "doublet_score" not in adata.obs:
        return {"error": "doublets_not_detected", "message": "Call detect_doublets first."}

    before = int(adata.n_obs)
    adata = adata[adata.obs["doublet_score"] < threshold].copy()
    SESSION.adata = adata
    after = int(adata.n_obs)
    checkpoint = SESSION.checkpoint("after_doublet_filter")
    return {
        "threshold": threshold,
        "before_cells": before,
        "after_cells": after,
        "removed": before - after,
        "pct_removed": round(100 * (before - after) / before, 2),
        "checkpoint": checkpoint,
    }


def normalize(target_sum: float = 1e4, n_top_genes: int = 2000) -> Summary:
    """Stash raw counts, log-normalize, and flag HVGs, then checkpoint. Mutating.

    Guardrails: raw counts are copied to layers['counts'] BEFORE normalizing (scVI needs
    them); refuses to run if X is not counts or if counts are already stashed (idempotency).
    HVGs are flagged, not subset — later steps (markers) need all genes.
    """
    adata = SESSION.require_adata()
    if "counts" in adata.layers:
        return {"error": "already_normalized", "message": "layers['counts'] exists; normalize already ran."}
    if not _looks_like_counts(adata.X):
        return {"error": "x_not_counts", "message": "X does not look like raw counts; refusing to normalize."}

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    n_hvgs = int(adata.var["highly_variable"].sum())

    checkpoint = SESSION.checkpoint("after_normalize")
    return {
        "target_sum": target_sum,
        "n_top_genes_requested": n_top_genes,
        "n_hvgs_flagged": n_hvgs,
        "raw_counts_stashed": True,
        "checkpoint": checkpoint,
    }


def run_pca(n_comps: int = 50) -> Summary:
    """PCA on the normalized matrix (HVGs), then checkpoint. Sets the session representation.

    Guardrail: requires normalize to have run (reads log-normalized X, not raw counts).
    Does not overwrite X (PCA zero-centers internally), so markers/annotation still see the
    log-normalized values.
    """
    adata = SESSION.require_adata()
    if "highly_variable" not in adata.var:
        return {"error": "not_normalized", "message": "Call normalize before run_pca."}

    n_hvgs = int(adata.var["highly_variable"].sum())
    n_comps = int(min(n_comps, adata.n_obs - 1, n_hvgs - 1))
    sc.pp.pca(adata, n_comps=n_comps, use_highly_variable=True)
    SESSION.representation = "X_pca"

    var_ratio = np.asarray(adata.uns["pca"]["variance_ratio"], dtype=float)
    checkpoint = SESSION.checkpoint("after_pca")
    return {
        "representation": "X_pca",
        "n_comps": n_comps,
        "variance_ratio_top10": [round(float(v), 4) for v in var_ratio[:10]],
        "cumulative_variance": round(float(var_ratio.sum()), 4),
        "checkpoint": checkpoint,
    }


def run_scvi(batch_key: str | None = None, max_epochs: int | None = None) -> Summary:
    """scVI latent from raw counts (layers['counts'], HVGs), then checkpoint.

    Guardrail: requires layers['counts'] (stashed by normalize) — scVI models counts, not
    log-normalized data. Trains on the HVG subset. Prefer PCA for a single clean batch;
    scVI's value is batch correction when batch_key spans multiple batches.
    """
    adata = SESSION.require_adata()
    if "counts" not in adata.layers:
        return {"error": "no_raw_counts", "message": "scVI needs raw counts; call normalize first."}
    if batch_key is not None and batch_key not in adata.obs:
        return {"error": "invalid_batch_key", "message": f"'{batch_key}' is not an obs column."}

    import scvi  # heavy (torch); import lazily so the rest of the toolset stays light

    scvi.settings.seed = config.SEED
    use_hvg = "highly_variable" in adata.var
    train = adata[:, adata.var["highly_variable"]].copy() if use_hvg else adata.copy()
    scvi.model.SCVI.setup_anndata(train, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(train)
    model.train(max_epochs=max_epochs)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    SESSION.representation = "X_scVI"
    try:
        epochs_trained = len(next(iter(model.history.values())))
    except (StopIteration, AttributeError):
        epochs_trained = None

    checkpoint = SESSION.checkpoint("after_scvi")
    return {
        "representation": "X_scVI",
        "n_latent": int(adata.obsm["X_scVI"].shape[1]),
        "batch_key": batch_key,
        "n_hvgs_used": int(train.n_vars),
        "epochs_trained": epochs_trained,
        "checkpoint": checkpoint,
    }


def cluster(resolution: float = 1.0) -> Summary:
    """Neighbor graph + Leiden clustering on the chosen representation, plus UMAP, then
    checkpoint.

    Reads SESSION.representation (X_pca or X_scVI, whichever the DR tool set) so clustering
    follows the agent's earlier choice. Guardrail: a representation must exist.
    """
    adata = SESSION.require_adata()
    rep = SESSION.representation
    if rep is None or rep not in adata.obsm:
        return {"error": "no_representation", "message": "Run run_pca or run_scvi before clustering."}

    sc.pp.neighbors(adata, use_rep=rep, random_state=config.SEED)
    sc.tl.leiden(
        adata, resolution=resolution, flavor="igraph", n_iterations=2,
        directed=False, random_state=config.SEED,
    )
    sc.tl.umap(adata, random_state=config.SEED)

    sizes = adata.obs["leiden"].value_counts().sort_index()
    checkpoint = SESSION.checkpoint("after_cluster")
    return {
        "representation_used": rep,
        "resolution": resolution,
        "n_clusters": int(sizes.shape[0]),
        "cluster_sizes": {str(k): int(v) for k, v in sizes.items()},
        "umap_computed": True,
        "checkpoint": checkpoint,
    }


def identify_markers(n_genes: int = 10) -> Summary:
    """Rank marker genes per Leiden cluster (Wilcoxon) and return the top genes per cluster.

    Reads the log-normalized X (preserved through PCA). Non-mutating (writes ranking to uns),
    so it does not checkpoint. Guardrail: requires clustering first.
    """
    adata = SESSION.require_adata()
    if "leiden" not in adata.obs:
        return {"error": "no_clusters", "message": "Run cluster before identify_markers."}

    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    names = adata.uns["rank_genes_groups"]["names"]
    groups = list(names.dtype.names)
    top = {g: [str(names[g][i]) for i in range(min(n_genes, len(names[g])))] for g in groups}

    return {
        "method": "wilcoxon",
        "n_clusters": len(groups),
        "n_genes_per_cluster": n_genes,
        "top_markers_per_cluster": top,
    }


def annotate_celltypes(model: str = "Immune_All_Low.pkl") -> Summary:
    """Annotate cell types with CellTypist (majority voting over Leiden clusters), then
    checkpoint.

    Our normalize output (1e4 counts + log1p) is exactly CellTypist's expected input. Uses
    the Leiden clusters for majority voting so labels align with clustering. Guardrail:
    requires clustering first. The resulting checkpoint is the annotated deliverable.
    """
    adata = SESSION.require_adata()
    if "leiden" not in adata.obs:
        return {"error": "no_clusters", "message": "Run cluster before annotate_celltypes."}

    import celltypist
    from celltypist import models

    models.download_models(model=[model], force_update=False)
    predictions = celltypist.annotate(
        adata, model=model, majority_voting=True, over_clustering="leiden"
    )
    labels = predictions.predicted_labels.loc[adata.obs_names, "majority_voting"].astype(str)
    adata.obs["cell_type"] = labels

    per_cluster = adata.obs.groupby("leiden", observed=True)["cell_type"].agg(
        lambda s: s.value_counts().index[0]
    )
    counts = adata.obs["cell_type"].value_counts()

    checkpoint = SESSION.checkpoint("after_annotate")
    return {
        "model": model,
        "n_cell_types": int(counts.shape[0]),
        "cell_type_counts": {str(k): int(v) for k, v in counts.items()},
        "per_cluster_majority": {str(k): str(v) for k, v in per_cluster.items()},
        "checkpoint": checkpoint,
    }


def summarize_findings() -> Summary:
    """Consolidate the final analysis state into one factual summary for the write-up.

    Non-mutating: reads the current adata (cell counts, cluster->cell-type mapping, top
    markers, cell-type counts) so the agent can compose the report narrative from a single
    joined view rather than re-reading each earlier tool result.
    """
    adata = SESSION.require_adata()
    out: Summary = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "representation": SESSION.representation,
    }
    if "leiden" in adata.obs:
        out["n_clusters"] = int(adata.obs["leiden"].nunique())
    if "cell_type" in adata.obs:
        out["cell_type_counts"] = {str(k): int(v) for k, v in adata.obs["cell_type"].value_counts().items()}
        per_cluster = adata.obs.groupby("leiden", observed=True)["cell_type"].agg(
            lambda s: s.value_counts().index[0]
        )
        out["cluster_to_cell_type"] = {str(k): str(v) for k, v in per_cluster.items()}
    if "rank_genes_groups" in adata.uns:
        names = adata.uns["rank_genes_groups"]["names"]
        out["top_markers_per_cluster"] = {
            g: [str(names[g][i]) for i in range(min(5, len(names[g])))] for g in names.dtype.names
        }
    return out


def generate_report(report_markdown: str) -> Summary:
    """Render figures, write the annotated .h5ad, and assemble outputs/report.md.

    The agent supplies the narrative (report_markdown); this tool handles the deterministic
    artifacts: UMAP + QC figures, the annotated deliverable, and stitching them into the
    report. Returns the paths written.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    adata = SESSION.require_adata()
    paths = SESSION.paths
    paths.figures.mkdir(parents=True, exist_ok=True)
    figures: list = []

    if "X_umap" in adata.obsm:
        colors = [c for c in ("cell_type", "leiden") if c in adata.obs]
        sc.pl.umap(adata, color=colors, show=False, wspace=0.4)
        path = paths.figures / "umap.png"
        plt.savefig(path, dpi=120, bbox_inches="tight")
        plt.close()
        figures.append(path)

    qc_cols = [c for c in ("n_genes_by_counts", "total_counts", "pct_counts_mt") if c in adata.obs]
    if qc_cols:
        sc.pl.violin(adata, qc_cols, multi_panel=True, show=False)
        path = paths.figures / "qc_violin.png"
        plt.savefig(path, dpi=120, bbox_inches="tight")
        plt.close()
        figures.append(path)

    adata.write_h5ad(paths.annotated)

    figures_md = ""
    if figures:
        figures_md = "\n\n## Figures\n\n" + "\n\n".join(
            f"![{p.stem}](figures/{p.name})" for p in figures
        )
    paths.report.write_text(report_markdown.rstrip() + figures_md + "\n")

    return {
        "report_path": str(paths.report),
        "figures": [str(p) for p in figures],
        "annotated_h5ad": str(paths.annotated),
        "report_chars": len(report_markdown),
    }
