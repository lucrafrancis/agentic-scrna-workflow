"""Fetch the Kang et al. 2018 IFN-β PBMC dataset into data/kang.h5ad.

PBMCs from 8 donors, each split into a control and an IFN-β stimulated sample (Kang et al.,
Nat Biotechnol 2018). Donors are replicates and every donor appears in both conditions, so
this is the case for paired pseudobulk differential expression and composition testing.
The file is the copy pertpy distributes (raw counts, already filtered by the authors, with
mitochondrial genes removed upstream).

We keep only the raw counts and the two design columns, renamed to `condition` (ctrl/stim)
and `donor`, so the agent starts from a pre-analysis state. The authors' embeddings,
clusters and normalised values are dropped. Their cell-type labels are written to a
separate file (data/kang_reference_labels.csv) for benchmarking, and are never shown to the
agent.

    uv run python scripts/fetch_kang.py
"""

from __future__ import annotations

import urllib.request

import anndata as ad
import pandas as pd

from agent import config

URL = "https://scverse-exampledata.s3.eu-west-1.amazonaws.com/pertpy/kang_2018.h5ad"
RAW = config.DATA_DIR / "kang_2018_raw.h5ad"
OUT = config.DATA_DIR / "kang.h5ad"
LABELS = config.DATA_DIR / "kang_reference_labels.csv"


def main() -> None:
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)
    if not RAW.exists():
        print(f"Downloading {URL}")
        urllib.request.urlretrieve(URL, RAW)
    src = ad.read_h5ad(RAW)

    obs = pd.DataFrame(
        {
            "condition": pd.Categorical(src.obs["label"].astype(str)),
            "donor": pd.Categorical(src.obs["replicate"].astype(str).str.replace("patient_", "donor_")),
        },
        index=src.obs_names,
    )
    adata = ad.AnnData(X=src.X.copy(), obs=obs, var=pd.DataFrame(index=src.var_names))
    adata.write_h5ad(OUT)

    src.obs[["cell_type"]].rename(columns={"cell_type": "reference_cell_type"}).to_csv(
        LABELS, index_label="barcode"
    )

    design = adata.obs.groupby(["donor", "condition"], observed=True).size().unstack()
    print(f"Wrote {OUT}  ({adata.n_obs} cells x {adata.n_vars} genes)")
    print(f"Wrote {LABELS}  (published cell types, for benchmarking only)")
    print(design)


if __name__ == "__main__":
    main()
