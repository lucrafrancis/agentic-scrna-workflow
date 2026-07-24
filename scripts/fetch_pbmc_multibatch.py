"""Fetch a two-batch PBMC dataset into data/pbmc_multibatch.h5ad.

Uses scvi-tools' built-in CITE-seq PBMCs, which concatenates two independent 10x runs
(PBMC5k + PBMC10k) with a `batch` key. We keep only the RNA counts and subsample each
batch so scVI training stays fast for a demo, while preserving the real technical batch
structure — the case where the agent should choose scVI over PCA.

    uv run python scripts/fetch_pbmc_multibatch.py
"""

from __future__ import annotations

import numpy as np
import scvi

from agent import config

OUT = config.DATA_DIR / "pbmc_multibatch.h5ad"
PER_BATCH = 2000  # cells to keep per batch (fast scVI without losing batch structure)


def main() -> None:
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    adata = scvi.data.pbmcs_10x_cite_seq(save_path=str(config.DATA_DIR))

    # Subsample each batch reproducibly.
    rng = np.random.default_rng(config.SEED)
    keep = []
    for b in adata.obs["batch"].unique():
        cells = np.flatnonzero(adata.obs["batch"].values == b)
        keep.append(rng.choice(cells, size=min(PER_BATCH, cells.size), replace=False))
    adata = adata[np.sort(np.concatenate(keep))].copy()

    # Keep it a clean RNA object (drop the CITE-seq protein matrix).
    adata.obsm.pop("protein_expression", None)
    # NB: barcodes collide across the two batches. We deliberately leave that for the
    # workflow to handle (Session.load makes them unique and inspect_dataset reports it),
    # so the pipeline is robust to any messy input rather than relying on clean fetches.

    adata.write_h5ad(OUT)
    counts = {str(k): int(v) for k, v in adata.obs["batch"].value_counts().items()}
    print(f"Wrote {OUT}  ({adata.n_obs} cells x {adata.n_vars} genes); batches: {counts}")


if __name__ == "__main__":
    main()
