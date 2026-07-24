"""Fetch the pbmc3k example dataset into data/pbmc3k.h5ad.

pbmc3k is the scanpy-standard demo dataset (~2700 PBMCs, single 10x batch). We save the
*raw counts* object so the agent starts from a realistic pre-analysis state and makes its
own QC / normalization / dimensionality-reduction choices.

    uv run python scripts/fetch_pbmc3k.py
"""

from __future__ import annotations

import scanpy as sc

from agent import config

OUT = config.DATA_DIR / "pbmc3k.h5ad"


def main() -> None:
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    # Raw-counts version (not the pre-processed pbmc3k_processed) — we want the agent to
    # start before QC/normalization.
    adata = sc.datasets.pbmc3k()
    adata.write_h5ad(OUT)
    print(f"Wrote {OUT}  ({adata.n_obs} cells x {adata.n_vars} genes)")


if __name__ == "__main__":
    main()
