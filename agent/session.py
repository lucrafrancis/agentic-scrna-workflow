"""Run-scoped state shared across tools (Option A in CLAUDE.md).

Tools operate on one AnnData for the duration of a run. The LLM only ever sees JSON
summaries, never the matrix, so tools reach the working object through this module-level
SESSION rather than receiving `adata` as an argument. The session also carries derived
run state (e.g. which representation clustering should use) and owns checkpointing.
"""

from __future__ import annotations

import anndata as ad

from agent import config


class Session:
    def __init__(self) -> None:
        self.adata: ad.AnnData | None = None
        # Set by the dimensionality-reduction tool; read by clustering. e.g. "X_pca" / "X_scVI".
        self.representation: str | None = None
        self._step = 0

    def load(self, path) -> ad.AnnData:
        self.adata = ad.read_h5ad(path)
        return self.adata

    def require_adata(self) -> ad.AnnData:
        if self.adata is None:
            raise RuntimeError("No dataset loaded into the session.")
        return self.adata

    def checkpoint(self, label: str) -> str:
        """Persist the current AnnData to a numbered checkpoint. Returns the path."""
        adata = self.require_adata()
        config.CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)
        self._step += 1
        path = config.CHECKPOINT_DIR / f"{self._step:02d}_{label}.h5ad"
        adata.write_h5ad(path)
        return str(path)


# Single run-scoped instance the tools import.
SESSION = Session()
