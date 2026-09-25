"""Analysis tools.

THE CONTRACT:
  Every tool is a deterministic Python function that does real work AND returns a
  structured, JSON-serializable summary dict. That dict is the *only* thing the LLM
  sees — it decides the next step entirely from it. A tool that returns None is a bug.

Tools operate on a shared AnnData held by the running session (see loop.py). Each tool
validates the state it expects and returns a clear error summary (not a stack trace) the
agent can react to.
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde

from agent import config, report
from agent.session import SESSION

Summary = dict[str, Any]

# obs columns that commonly denote a technical batch/sample grouping, checked by
# inspect_dataset. Deliberately excludes experimental-design columns such as `condition` or
# `treatment`: correcting those away would remove the biology the analysis is about.
_BATCH_HINTS = {"batch", "sample", "donor", "patient", "subject", "dataset"}

# inspect_dataset lists the values of obs columns with at most this many distinct values.
_MAX_LEVELS_SHOWN = 20

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


def _as_float(value) -> float | None:
    return float(value) if value is not None else None


def _detect_batch_key(adata) -> str | None:
    """First obs column matching a batch hint with more than one value, else None."""
    for col in adata.obs.columns:
        if col.lower() in _BATCH_HINTS and adata.obs[col].nunique() > 1:
            return col
    return None


def inspect_dataset() -> Summary:
    """n cells/genes, obs/var columns, candidate batch key, and count/layer state.

    Non-mutating; the agent calls this first to orient itself. Reports enough for the
    dimensionality-reduction decision (batch key) and the counts guardrail (has_raw_counts).
    """
    adata = SESSION.require_adata()
    obs_cols = list(adata.obs.columns)

    candidate_batch_key = _detect_batch_key(adata)
    n_batches = int(adata.obs[candidate_batch_key].nunique()) if candidate_batch_key else None

    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": obs_cols,
        # Levels of the low-cardinality columns, so the agent can tell a condition to keep
        # (e.g. ctrl/stim) from a replicate or batch grouping (e.g. donors).
        "obs_levels": {
            c: sorted(map(str, adata.obs[c].unique()))
            for c in obs_cols
            if adata.obs[c].nunique() <= _MAX_LEVELS_SHOWN
            and not pd.api.types.is_float_dtype(adata.obs[c])
        },
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
    if SESSION.qc_before_filter is None:
        SESSION.qc_before_filter = adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
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


def detect_doublets(batch_key: str | None = "auto") -> Summary:
    """Run Scrublet and report the score distribution plus candidate thresholds.

    Non-mutating: stores per-cell doublet scores in obs but removes nothing. The agent
    inspects the distribution, chooses a threshold, and calls filter_doublets. Recommends
    the bimodal-valley cutoff when the distribution is bimodal, else median + 3*MAD.

    Guardrails: Scrublet models raw counts, so this refuses to run after normalize (scanpy
    only warns on stderr, which the agent never sees). Scrublet runs per batch, so the
    simulated-doublet model is built within a single 10x run rather than across pooled
    runs. batch_key="auto" uses the detected batch column; the agent can instead name the
    column that marks 10x runs (which may differ from donors when donors were pooled into
    one run), or pass None.
    """
    adata = SESSION.require_adata()
    if "counts" in adata.layers or not _looks_like_counts(adata.X):
        return {
            "error": "x_not_counts",
            "message": "Scrublet needs raw counts; run detect_doublets before normalize.",
        }

    # Scrublet should run within each physical 10x run. "auto" uses the detected batch
    # column; the agent can name the column that marks runs, or None for one pooled run.
    if batch_key == "auto":
        batch_key = _detect_batch_key(adata)
    elif batch_key is not None and batch_key not in adata.obs:
        return {"error": "invalid_batch_key", "message": f"'{batch_key}' is not an obs column."}
    sc.pp.scrublet(adata, batch_key=batch_key, random_state=config.SEED)
    scores = np.asarray(adata.obs["doublet_score"], dtype=float)
    SESSION.doublet_scores = scores.copy()

    median = float(np.median(scores))
    mad = float(np.median(np.abs(scores - median)))  # raw (unscaled) MAD, per project spec
    mad_threshold = median + 3 * mad
    valley_threshold, is_bimodal = _valley_threshold(scores)
    # Batched runs get one automatic threshold per batch, not a single global one.
    scrublet_uns = adata.uns.get("scrublet", {})
    scrublet_auto_per_batch = (
        {str(b): _as_float(u.get("threshold")) for b, u in scrublet_uns.get("batches", {}).items()}
        if batch_key
        else None
    )
    scrublet_auto = None if batch_key else _as_float(scrublet_uns.get("threshold"))

    recommended = valley_threshold if is_bimodal else mad_threshold
    rule = "bimodal_valley" if is_bimodal else "median+3*MAD"
    n_flagged = int((scores >= recommended).sum())
    counts, edges = np.histogram(scores, bins=20)

    return {
        "is_bimodal": is_bimodal,
        "batch_key": batch_key,
        "score_distribution": _quantile_summary(scores),
        "histogram": {"counts": counts.tolist(), "bin_edges": [round(e, 4) for e in edges]},
        "candidate_thresholds": {
            "bimodal_valley": valley_threshold,
            "median_3mad": round(mad_threshold, 4),
            "scrublet_auto": scrublet_auto,
            "scrublet_auto_per_batch": scrublet_auto_per_batch,
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
    HVGs are flagged, not subset — later steps (markers) need all genes. When the data has a
    batch key, HVGs are ranked within each batch and combined, so the selection is not driven
    by the batch effect that run_scvi is there to correct.
    """
    adata = SESSION.require_adata()
    if "counts" in adata.layers:
        return {"error": "already_normalized", "message": "layers['counts'] exists; normalize already ran."}
    if not _looks_like_counts(adata.X):
        return {"error": "x_not_counts", "message": "X does not look like raw counts; refusing to normalize."}

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    batch_key = _detect_batch_key(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, batch_key=batch_key)
    n_hvgs = int(adata.var["highly_variable"].sum())
    # annotate_celltypes needs to know the scale X is on: CellTypist expects 1e4 + log1p.
    SESSION.normalize_target_sum = float(target_sum)

    checkpoint = SESSION.checkpoint("after_normalize")
    return {
        "target_sum": target_sum,
        "n_top_genes_requested": n_top_genes,
        "n_hvgs_flagged": n_hvgs,
        "hvg_batch_key": batch_key,
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


def run_scvi(batch_key: str | list[str] | None = None, max_epochs: int | None = None) -> Summary:
    """scVI latent from raw counts (layers['counts'], HVGs), then checkpoint.

    Guardrail: requires layers['counts'] (stashed by normalize) — scVI models counts, not
    log-normalized data. Trains on the HVG subset. Prefer PCA for a single clean batch;
    scVI's value is batch correction when batch_key spans multiple batches. A list of
    columns is combined into one batch per combination (e.g. ['donor', 'condition'] gives
    one batch per sample), stored as obs['scvi_batch'].
    """
    adata = SESSION.require_adata()
    if "counts" not in adata.layers:
        return {"error": "no_raw_counts", "message": "scVI needs raw counts; call normalize first."}
    keys = [batch_key] if isinstance(batch_key, str) else list(batch_key or [])
    missing = [k for k in keys if k not in adata.obs]
    if missing:
        return {"error": "invalid_batch_key", "message": f"Not obs columns: {missing}."}
    if len(keys) > 1:
        adata.obs["scvi_batch"] = pd.Categorical(adata.obs[keys].astype(str).agg("_".join, axis=1))
        batch_key = "scvi_batch"
    elif keys:
        batch_key = keys[0]

    import scvi  # heavy (torch); import lazily so the rest of the toolset stays light

    scvi.settings.seed = config.SEED
    use_hvg = "highly_variable" in adata.var
    train = adata[:, adata.var["highly_variable"]].copy() if use_hvg else adata.copy()
    scvi.model.SCVI.setup_anndata(train, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(train)
    model.train(max_epochs=max_epochs)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    SESSION.representation = "X_scVI"
    SESSION.scvi_batch_columns = keys
    try:
        epochs_trained = len(next(iter(model.history.values())))
    except (StopIteration, AttributeError):
        epochs_trained = None

    checkpoint = SESSION.checkpoint("after_scvi")
    return {
        "representation": "X_scVI",
        "n_latent": int(adata.obsm["X_scVI"].shape[1]),
        "batch_key": batch_key,
        "batch_columns": keys,
        "n_batches": int(adata.obs[batch_key].nunique()) if batch_key else 1,
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


def identify_markers(n_genes: int = 25) -> Summary:
    """Rank marker genes per Leiden cluster (Wilcoxon, each cluster vs all other cells) and
    return the top genes per cluster.

    Every gene is ranked, not just the top n_genes, and the full table (rank, log2FC,
    adjusted p, % of cells expressing in and outside the cluster) is kept on the session so
    check_markers can look up any gene. Reads the log-normalized X (preserved through PCA).
    Non-mutating (writes the ranking to uns), so it does not checkpoint. Guardrail: requires
    clustering first.
    """
    adata = SESSION.require_adata()
    if "leiden" not in adata.obs:
        return {"error": "no_clusters", "message": "Run cluster before identify_markers."}

    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", n_genes=adata.n_vars, pts=True)
    table = sc.get.rank_genes_groups_df(adata, group=None)
    if "group" not in table:  # a single cluster comes back without a group column
        table.insert(0, "group", adata.obs["leiden"].astype(str).iloc[0])
    table["group"] = table["group"].astype(str)
    table["rank"] = table.groupby("group").cumcount() + 1
    SESSION.marker_table = table.set_index(["group", "names"])

    names = adata.uns["rank_genes_groups"]["names"]
    groups = list(names.dtype.names)
    top = {g: [str(names[g][i]) for i in range(min(n_genes, len(names[g])))] for g in groups}
    return {
        "method": "wilcoxon",
        "n_clusters": len(groups),
        "n_genes_per_cluster": n_genes,
        "n_genes_ranked": int(adata.n_vars),
        "top_markers_per_cluster": top,
        "note": "Every gene is ranked; use check_markers to look up specific genes in any cluster.",
    }


def check_markers(genes: list[str], clusters: list[str] | None = None) -> Summary:
    """Marker statistics for specific genes in each cluster, from identify_markers' full
    ranking: the gene's rank among all genes for that cluster, log2 fold change and adjusted
    p (cluster vs all other cells), and the % of cells expressing it in and outside the
    cluster. The % expressing is what shows a marker is absent. Non-mutating.

    The p-values treat cells as independent, so they are far too small for inference; use
    them to rank and describe markers, not as proof.
    """
    table = SESSION.marker_table
    if table is None:
        return {"error": "no_markers", "message": "Run identify_markers before check_markers."}
    all_clusters = sorted(table.index.get_level_values("group").unique(), key=lambda c: int(c) if c.isdigit() else c)
    clusters = [str(c) for c in (clusters or all_clusters)]
    unknown_clusters = [c for c in clusters if c not in all_clusters]
    if unknown_clusters:
        return {"error": "bad_clusters", "message": f"Unknown clusters {unknown_clusters}. Clusters: {all_clusters}"}

    n_ranked = int(table.loc[all_clusters[0]].shape[0])
    genes_known = set(table.index.get_level_values("names"))
    out: dict[str, dict] = {}
    not_found = []
    for gene in genes:
        if gene not in genes_known:
            not_found.append(gene)
            continue
        out[gene] = {}
        for cl in clusters:
            r = table.loc[(cl, gene)]
            out[gene][cl] = {
                "rank": f"{int(r['rank'])} of {n_ranked}",
                "log2FC": round(float(r["logfoldchanges"]), 2),
                "padj": float(f"{r['pvals_adj']:.2g}"),
                "pct_in_cluster": round(100 * float(r["pct_nz_group"]), 1),
                "pct_elsewhere": round(100 * float(r["pct_nz_reference"]), 1),
            }
    result: Summary = {
        "genes": out,
        "note": "Cluster vs all other cells (Wilcoxon). p-values treat cells as independent: use them to "
        "rank and describe markers, not as proof. Low pct_in_cluster means the marker is absent.",
    }
    if not_found:
        result["not_found"] = not_found
    return result


_FINE_MODEL, _COARSE_MODEL = "Immune_All_Low.pkl", "Immune_All_High.pkl"


def annotate_celltypes(model: str = _FINE_MODEL) -> Summary:
    """Annotate cell types with CellTypist (majority voting over Leiden clusters), then
    checkpoint.

    CellTypist expects 1e4-normalized, log1p data. That is what normalize produces by
    default, but target_sum is agent-controllable: if this run normalized to anything else
    (e.g. 1e6 / CPM), CellTypist would only warn and return degraded labels, so we rebuild
    its input from the stashed raw counts instead. Uses the Leiden clusters for majority
    voting so labels align with clustering. Guardrail: requires clustering first. The
    resulting checkpoint is the annotated deliverable. Also runs the other CellTypist immune
    model (coarse if fine was chosen, and vice versa) and reports its per-cluster labels as a
    second opinion, stored in obs['cell_type_alt'].
    """
    adata = SESSION.require_adata()
    if "leiden" not in adata.obs:
        return {"error": "no_clusters", "message": "Run cluster before annotate_celltypes."}

    target_sum = SESSION.normalize_target_sum
    rescaled = False
    ct_input = adata
    if target_sum is not None and not np.isclose(target_sum, 1e4):
        if "counts" not in adata.layers:
            return {
                "error": "wrong_normalization",
                "message": f"CellTypist needs target_sum=1e4; this run used {target_sum:g} and "
                "layers['counts'] is missing, so the input cannot be rescaled.",
            }
        ct_input = adata.copy()
        ct_input.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(ct_input, target_sum=1e4)
        sc.pp.log1p(ct_input)
        rescaled = True

    import celltypist
    from celltypist import models

    def _majority_labels(model_name: str) -> pd.Series:
        models.download_models(model=[model_name], force_update=False)
        predictions = celltypist.annotate(
            ct_input, model=model_name, majority_voting=True, over_clustering="leiden"
        )
        return predictions.predicted_labels.loc[adata.obs_names, "majority_voting"].astype(str)

    def _per_cluster(column: str) -> dict[str, str]:
        top = adata.obs.groupby("leiden", observed=True)[column].agg(lambda s: s.value_counts().index[0])
        return {str(k): str(v) for k, v in top.items()}

    adata.obs["cell_type"] = _majority_labels(model)
    counts = adata.obs["cell_type"].value_counts()
    out: Summary = {
        "model": model,
        "normalize_target_sum": target_sum,
        "rescaled_to_1e4_for_celltypist": rescaled,
        "n_cell_types": int(counts.shape[0]),
        "cell_type_counts": {str(k): int(v) for k, v in counts.items()},
        "per_cluster_majority": _per_cluster("cell_type"),
    }
    # A second opinion at the other granularity: fine labels can be wrong for cells unlike the
    # reference (e.g. cultured monocytes called macrophages), where the coarse label is right.
    other = _COARSE_MODEL if model != _COARSE_MODEL else _FINE_MODEL
    adata.obs["cell_type_alt"] = _majority_labels(other)
    out["second_opinion"] = {
        "model": other,
        "per_cluster_majority": _per_cluster("cell_type_alt"),
        "note": "Where the two models disagree beyond granularity, check the cluster's markers. To use "
        f"the other model's labels, call annotate_celltypes(model='{other}').",
    }
    out["checkpoint"] = SESSION.checkpoint("after_annotate")
    return out


def relabel_clusters(labels: dict[str, str], reason: str) -> Summary:
    """Override the cell-type label of whole Leiden clusters, with a stated reason.

    For when marker genes contradict CellTypist (e.g. cells unlike its reference). The
    original CellTypist labels are kept in obs['cell_type_celltypist'] and every relabelling
    is recorded for the report. Composition and DE results computed with the old labels are
    cleared, so they must be rerun. Guardrails: clusters must exist, labels must be
    non-empty, and a reason is required.
    """
    adata = SESSION.require_adata()
    if "cell_type" not in adata.obs or "leiden" not in adata.obs:
        return {"error": "not_annotated", "message": "Run annotate_celltypes before relabel_clusters."}
    if not reason.strip():
        return {"error": "no_reason", "message": "Give the marker evidence for the new labels."}
    clusters = set(adata.obs["leiden"].astype(str))
    unknown = [c for c in labels if str(c) not in clusters]
    empty = [c for c, v in labels.items() if not str(v).strip()]
    if unknown or empty:
        return {"error": "bad_labels", "message": f"Unknown clusters {unknown}; empty labels for {empty}. "
                f"Clusters: {sorted(clusters, key=lambda c: int(c) if c.isdigit() else c)}"}

    if "cell_type_celltypist" not in adata.obs:
        adata.obs["cell_type_celltypist"] = adata.obs["cell_type"].astype(str)
    leiden = adata.obs["leiden"].astype(str)
    current = adata.obs["cell_type"].astype(str)
    changes = []
    for cluster, label in labels.items():
        mask = leiden == str(cluster)
        before = current[mask].value_counts().index[0]
        current = current.where(~mask, label.strip())
        changes.append({"cluster": str(cluster), "from": before, "to": label.strip(), "n_cells": int(mask.sum())})
    adata.obs["cell_type"] = pd.Categorical(current)
    SESSION.relabels.append({"changes": changes, "reason": reason.strip()})

    cleared = []
    if SESSION.composition is not None:
        SESSION.composition = None
        cleared.append("compare_composition")
    if SESSION.de_results:
        SESSION.de_results, SESSION.de_settings = {}, None
        cleared.append("run_pseudobulk_de")
    counts = adata.obs["cell_type"].value_counts()
    out: Summary = {
        "changes": changes,
        "cell_type_counts": {str(k): int(v) for k, v in counts.items()},
        "checkpoint": SESSION.checkpoint("after_relabel"),
    }
    if cleared:
        out["cleared"] = cleared
        out["message"] = f"Results from {cleared} used the old labels and were cleared; rerun them."
    return out


# --- Comparing conditions -------------------------------------------------------
#
# Both tools compare a condition (e.g. ctrl vs stim) across biological replicates (e.g.
# donors). The unit of replication is the sample, never the cell: cells from one donor are
# not independent, so testing cells directly overstates significance.

_DE_PADJ = 0.05  # significance threshold for DE counts, top-gene lists and figures


def _check_columns(adata, columns: list[str]) -> Summary | None:
    missing = [c for c in columns if c not in adata.obs]
    if missing:
        return {"error": "bad_column", "message": f"Not obs columns: {missing}. Columns: {list(adata.obs.columns)}"}
    return None


def _check_contrast(adata, contrast: list[str]) -> Summary | None:
    if len(contrast) != 3:
        return {"error": "bad_contrast", "message": "Contrast must be [factor, test, reference]."}
    factor, test, ref = contrast
    if err := _check_columns(adata, [factor]):
        return err
    levels = sorted(map(str, adata.obs[factor].unique()))
    for level in (test, ref):
        if level not in levels:
            return {"error": "bad_level", "message": f"'{level}' is not a level of '{factor}'. Levels: {levels}"}
    if test == ref:
        return {"error": "bad_contrast", "message": "Test and reference levels must differ."}
    return None


# Fewest samples per condition compare_composition will test; below this no test has power.
_MIN_COMPOSITION_SAMPLES = 3


def compare_composition(contrast: list[str], sample_key: str, celltype_key: str = "cell_type") -> Summary:
    """Compare cell-type proportions between two conditions across replicate samples.

    Proportions are computed per sample (sample_key x condition). If every sample appears in
    both conditions the test is paired (Wilcoxon signed-rank on per-sample proportions),
    otherwise unpaired (Mann-Whitney U); p-values are Benjamini-Hochberg adjusted across cell
    types. Guardrail: refuses with fewer than _MIN_COMPOSITION_SAMPLES samples per condition,
    where no test has power. Proportions are compositional (one type rising forces others down), which
    the summary notes. Non-mutating.
    """
    from scipy.stats import mannwhitneyu, wilcoxon
    from statsmodels.stats.multitest import multipletests

    adata = SESSION.require_adata()
    if err := _check_contrast(adata, contrast) or _check_columns(adata, [sample_key, celltype_key]):
        return err
    factor, test, ref = contrast
    if sample_key == factor:
        return {"error": "bad_sample_key", "message": "sample_key must identify replicates, not the condition."}

    obs = adata.obs[[factor, sample_key, celltype_key]].astype(str)
    obs = obs[obs[factor].isin([test, ref])]
    counts = obs.groupby([sample_key, factor, celltype_key]).size().unstack(fill_value=0)
    props = counts.div(counts.sum(axis=1), axis=0)
    per_level = props.index.get_level_values(factor).value_counts()
    if min(per_level.get(test, 0), per_level.get(ref, 0)) < _MIN_COMPOSITION_SAMPLES:
        return {
            "error": "too_few_samples",
            "message": f"Need at least {_MIN_COMPOSITION_SAMPLES} samples per condition; have {per_level.to_dict()}. "
            "Describe proportions without a test.",
        }

    by_level = {lvl: props.xs(lvl, level=factor) for lvl in (test, ref)}
    shared = by_level[test].index.intersection(by_level[ref].index)
    paired = len(shared) == len(by_level[test]) == len(by_level[ref])
    rows = []
    for ct in props.columns:
        t, r = by_level[test][ct], by_level[ref][ct]
        if paired:
            t, r = t.loc[shared], r.loc[shared]
            p = float(wilcoxon(t, r).pvalue) if (t - r).abs().sum() > 0 else 1.0
        else:
            p = float(mannwhitneyu(t, r).pvalue)
        eps = 1e-4  # keeps the log ratio finite when a type is absent from a sample
        rows.append({
            "cell_type": ct,
            f"mean_prop_{ref}": float(r.mean()),
            f"mean_prop_{test}": float(t.mean()),
            "log2_ratio": float(np.log2((t.mean() + eps) / (r.mean() + eps))),
            "n_samples_higher_in_test": int((t.values > r.values).sum()) if paired else None,
            "pvalue": p,
        })
    table = pd.DataFrame(rows).set_index("cell_type")
    table["padj"] = multipletests(table["pvalue"], method="fdr_bh")[1]
    table = table.sort_values("padj")

    SESSION.composition = {
        "table": table,
        "proportions": props,
        "settings": {"contrast": list(contrast), "sample_key": sample_key, "celltype_key": celltype_key,
                     "paired": paired, "test": "Wilcoxon signed-rank" if paired else "Mann-Whitney U"},
    }
    return {
        "contrast": contrast,
        "paired": paired,
        "test": SESSION.composition["settings"]["test"],
        "n_samples": {lvl: int(len(by_level[lvl])) for lvl in (test, ref)},
        "results": {
            ct: {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}
            for ct, row in table.to_dict("index").items()
        },
        "n_significant": int((table["padj"] < _DE_PADJ).sum()),
        "note": "Proportions are compositional: a large shift in one type moves the others. With few "
        "samples, the smallest attainable p-value is limited (paired Wilcoxon, n=8: p >= 0.0078).",
    }


def _pseudobulk(adata, counts, groups: list[str], min_cells: int):
    """Sum raw counts per group of cells. Returns (samples x genes counts, sample metadata
    with n_cells), dropping groups with fewer than min_cells cells."""
    keys = adata.obs[groups].astype(str).agg("|".join, axis=1)
    codes, uniques = pd.factorize(keys)
    indicator = sparse.csr_matrix(
        (np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniques), len(codes))
    )
    summed = indicator @ counts
    summed = summed.toarray() if sparse.issparse(summed) else np.asarray(summed)
    meta = pd.DataFrame([u.split("|") for u in uniques], columns=groups, index=uniques)
    meta["n_cells"] = np.bincount(codes)
    keep = meta["n_cells"].to_numpy() >= min_cells
    df = pd.DataFrame(summed[keep].round().astype(np.int64), index=uniques[keep], columns=adata.var_names)
    return df, meta[keep]


def run_pseudobulk_de(
    contrast: list[str],
    sample_key: str,
    covariates: list[str] | None = None,
    celltype_key: str = "cell_type",
    cell_types: list[str] | None = None,
    min_cells: int = 10,
) -> Summary:
    """Pseudobulk differential expression per cell type with PyDESeq2 (Wald test).

    For each cell type, raw counts are summed per sample (sample_key x condition). Samples
    with fewer than min_cells cells are dropped. DESeq2 is fitted with design
    ~ covariates + factor and the contrast [factor, test, reference] is tested with a Wald
    test. Adding the replicate column (e.g. donor) as a covariate gives a paired design.

    Guardrails: uses raw counts (layers['counts'], or X before normalize); refuses unknown
    columns or levels, covariates that vary within a pseudobulk sample, and confounded
    designs (a covariate that fully determines the condition). A cell type is skipped, with
    the reason, when fewer than 2 samples per condition remain or no residual degrees of
    freedom are left.
    """
    import logging

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    adata = SESSION.require_adata()
    covariates = list(covariates or [])
    if err := _check_contrast(adata, contrast) or _check_columns(adata, [sample_key, celltype_key, *covariates]):
        return err
    factor, test, ref = contrast
    if factor in covariates:
        return {"error": "bad_covariates", "message": f"'{factor}' is the contrast factor; don't list it as a covariate."}
    if sample_key == factor:
        return {"error": "bad_sample_key", "message": "sample_key must identify replicates, not the condition."}
    bad_names = [c for c in [factor, *covariates] if not str(c).isidentifier()]
    if bad_names:
        return {"error": "bad_column_name", "message": f"Design columns must be identifiers for the formula: {bad_names}"}

    if "counts" in adata.layers:
        counts = adata.layers["counts"]
    elif _looks_like_counts(adata.X):
        counts = adata.X
    else:
        return {"error": "no_raw_counts", "message": "Pseudobulk DE needs raw counts (layers['counts'])."}

    # A covariate must be constant within each pseudobulk sample, or summing mixes its levels.
    groups = [sample_key, factor]
    per_sample = adata.obs[groups + covariates].astype(str).drop_duplicates()
    varying = [c for c in covariates if c not in groups and per_sample.groupby(groups)[c].nunique().max() > 1]
    if varying:
        return {"error": "covariate_varies_within_sample",
                "message": f"{varying} vary within a {sample_key} x {factor} sample; pseudobulk can't adjust for them."}

    available = sorted(map(str, adata.obs[celltype_key].unique()))
    if cell_types is None:
        cell_types = available
    unknown = [c for c in cell_types if c not in available]
    if unknown:
        return {"error": "bad_cell_types", "message": f"Unknown cell types {unknown}. Available: {available}"}

    formula = "~" + " + ".join([*covariates, factor])
    in_contrast = adata.obs[factor].astype(str).isin([test, ref]).to_numpy()
    SESSION.de_results = {}
    results, skipped = {}, {}
    logging.getLogger("pydeseq2").setLevel(logging.ERROR)
    for ct in cell_types:
        mask = in_contrast & (adata.obs[celltype_key].astype(str) == ct).to_numpy()
        sub = adata[mask]
        df, meta = _pseudobulk(sub, counts[mask], sorted(set(groups + covariates)), min_cells)
        sizes = meta[factor].value_counts()
        if min(sizes.get(test, 0), sizes.get(ref, 0)) < 2:
            skipped[ct] = f"fewer than 2 samples per condition with >= {min_cells} cells ({sizes.to_dict()})"
            continue
        design = meta[[*covariates, factor]].astype(str)
        X = pd.get_dummies(design, drop_first=True).astype(float)
        X.insert(0, "intercept", 1.0)
        rank = int(np.linalg.matrix_rank(X.to_numpy()))
        if rank < X.shape[1]:
            return {"error": "confounded_design", "message": (
                f"In '{ct}', covariates {covariates} are confounded with '{factor}' (design rank {rank} < "
                f"{X.shape[1]}). Drop the confounded covariate; its effect can't be separated.")}
        if len(design) - rank < 1:
            skipped[ct] = "no residual degrees of freedom (too few samples for the design)"
            continue
        # Genes with at least 10 counts in at least as many samples as the smaller condition.
        min_group = int(min(sizes[test], sizes[ref]))
        df = df.loc[:, (df >= 10).sum(axis=0) >= min_group]

        dds = DeseqDataSet(counts=df, metadata=design, design=formula, quiet=True)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=[factor, test, ref], quiet=True)
        stats.summary()
        res = stats.results_df.copy()
        res.attrs["n_samples"] = {str(k): int(v) for k, v in sizes.items()}
        SESSION.de_results[ct] = res

        sig = res[res["padj"] < _DE_PADJ]
        results[ct] = {
            "n_samples": sizes.to_dict(),
            "n_cells": int(meta["n_cells"].sum()),
            "genes_tested": int(res["padj"].notna().sum()),
            "n_significant": int(len(sig)),
            "n_up": int((sig["log2FoldChange"] > 0).sum()),
            "n_down": int((sig["log2FoldChange"] < 0).sum()),
            "top_up": _top_genes(res, "up", 5),
            "top_down": _top_genes(res, "down", 5),
        }

    SESSION.de_settings = {"contrast": list(contrast), "sample_key": sample_key, "covariates": covariates,
                           "celltype_key": celltype_key, "min_cells": min_cells, "design": formula,
                           "skipped": skipped}
    out: Summary = {"design": formula, "contrast": contrast, "results": results, "skipped": skipped}
    if not results:
        out["warning"] = "No cell type could be tested."
    return out


def _top_genes(res: pd.DataFrame, direction: str, n: int) -> list[dict]:
    sig = res[res["padj"] < _DE_PADJ]
    sig = sig[sig["log2FoldChange"] > 0] if direction == "up" else sig[sig["log2FoldChange"] < 0] if direction == "down" else sig
    ranked = sig.assign(_abs=sig["log2FoldChange"].abs()).sort_values(["padj", "_abs"], ascending=[True, False])
    return [{"gene": str(g), "log2FC": round(float(r["log2FoldChange"]), 3), "padj": float(f"{r['padj']:.3g}")}
            for g, r in ranked.head(n).iterrows()]


def get_top_genes(cell_type: str, n: int = 20, direction: str = "both") -> Summary:
    """Top significant DE genes for one cell type, ranked by padj then |log2FC|. Non-mutating."""
    if cell_type not in SESSION.de_results:
        return {"error": "no_de_results", "message": f"No DE results for '{cell_type}'. "
                f"Tested: {sorted(SESSION.de_results)}"}
    return {"cell_type": cell_type, "direction": direction,
            "genes": _top_genes(SESSION.de_results[cell_type], direction, n)}


def query_genes(genes: list[str], cell_types: list[str] | None = None) -> Summary:
    """log2FC and padj for specific genes in each tested cell type (or the ones given), so the
    agent can check a gene before writing about it. Non-mutating."""
    if not SESSION.de_results:
        return {"error": "no_de_results", "message": "Run run_pseudobulk_de first."}
    cts = cell_types or sorted(SESSION.de_results)
    out: dict[str, dict] = {}
    for gene in genes:
        out[gene] = {}
        for ct in cts:
            res = SESSION.de_results.get(ct)
            if res is None or gene not in res.index or pd.isna(res.loc[gene, "padj"]):
                out[gene][ct] = "not tested"
                continue
            r = res.loc[gene]
            out[gene][ct] = {"log2FC": round(float(r["log2FoldChange"]), 3), "padj": float(f"{r['padj']:.3g}"),
                             "significant": bool(r["padj"] < _DE_PADJ)}
    return {"genes": out}


def summarize_findings() -> Summary:
    """Consolidate the final analysis state into one factual summary for the write-up.

    Non-mutating: reads the current adata (cell counts, cluster->cell-type mapping, top
    markers, cell-type counts) so the agent can compose the report narrative from a single
    joined view rather than re-reading each earlier tool result. Also returns the facts and
    tables the report cites through placeholders.
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
    # What the report may cite: placeholder names with their current values. The agent reads
    # the values to interpret them but writes only the placeholders.
    log = report.read_log(SESSION.paths.tool_log)
    out["facts"] = report.facts(log, adata)
    out["tables_available"] = sorted(report.tables(log, adata))
    out["tables_required"] = [t for t in report.REQUIRED_TABLES if t in out["tables_available"]]
    return out


# Rejections before generate_report writes the report anyway, with the failed checks shown.
_MAX_REPORT_REJECTIONS = 2


def generate_report(report_markdown: str) -> Summary:
    """Check and render the narrative, draw figures, and write the report and annotated .h5ad.

    Every number in the narrative must be a placeholder (see summarize_findings), so values
    come from code rather than the model. Any problem rejects the whole report with a list
    of all issues and nothing is written; after _MAX_REPORT_REJECTIONS it is written with
    the unresolved problems shown at the top. report.build_report then adds the run
    summary, decisions table, captioned figures and Methods from the tool log.
    """
    adata = SESSION.require_adata()
    paths = SESSION.paths
    log = report.read_log(paths.tool_log)
    narrative, problems = report.render(report_markdown, log, adata)
    if problems and SESSION.report_attempts < _MAX_REPORT_REJECTIONS:
        SESSION.report_attempts += 1
        return {
            "error": "report_rejected",
            "message": "Nothing was written. Fix ALL of these problems, then call generate_report "
            "again with the complete report.",
            "problems": problems,
            "available_facts": sorted(report.facts(log, adata)),
            "available_tables": sorted(report.tables(log, adata)),
        }

    text, figures = report.build_report(adata, narrative, _detect_batch_key(adata), warnings=problems)
    adata.write_h5ad(paths.annotated)
    paths.report.write_text(text)

    out: Summary = {
        "report_path": str(paths.report),
        "figures": [str(paths.figures / f.name) for f in figures],
        "annotated_h5ad": str(paths.annotated),
        "report_chars": len(report_markdown),
    }
    if problems:
        out["warning"] = f"Written after {SESSION.report_attempts} rejections with these problems shown: {problems}"
    return out
