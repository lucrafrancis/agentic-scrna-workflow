"""Run-scoped state shared across tools.

Tools operate on one AnnData for the duration of a run. The LLM only ever sees JSON
summaries, never the matrix, so tools reach the working object through this module-level
SESSION rather than receiving `adata` as an argument. The session also carries derived
run state (e.g. which representation clustering should use) and owns checkpointing.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad

from agent.config import RunPaths


class Session:
    def __init__(self) -> None:
        self.adata: ad.AnnData | None = None
        # Set by load(): the dataset name (input filename stem) and its output locations.
        self.name: str | None = None
        self.paths: RunPaths | None = None
        # Set by the dimensionality-reduction tool; read by clustering. e.g. "X_pca" / "X_scVI".
        self.representation: str | None = None
        # Set by check_gene_identifiers; read by compute_qc.
        self.gene_format: str | None = None  # "symbol" | "ensembl" | "other"
        self.mito_prefix: str | None = None  # e.g. "MT-" (human) / "mt-" (mouse); None if unknown
        # Set by normalize; read by annotate_celltypes (CellTypist expects 1e4 + log1p).
        self.normalize_target_sum: float | None = None
        # Structural hygiene applied at load (reported by inspect_dataset), not agent decisions.
        self.n_duplicate_barcodes = 0
        self._step = 0

    def load(self, path) -> ad.AnnData:
        self.name = Path(path).stem
        self.paths = RunPaths(self.name)
        adata = ad.read_h5ad(path)
        # Deterministic hygiene (no judgement): cell barcodes must be unique for downstream
        # keying. Record how many collided so inspect_dataset can surface it in the run.
        self.n_duplicate_barcodes = int(adata.obs_names.duplicated().sum())
        if self.n_duplicate_barcodes:
            adata.obs_names_make_unique()
        self.adata = adata
        return self.adata

    def begin_run(self) -> None:
        """Clear this dataset's output dir of the previous run's provenance artifacts.

        The tool log is append-only and checkpoint numbering restarts at 1 each run, so
        re-running a dataset would otherwise interleave two runs in one tool_calls.jsonl and
        leave run 1's checkpoints (e.g. 06_after_scvi.h5ad) sitting next to run 2's as if
        they were one sequence. The audit trail is a deliverable — it must describe exactly
        one run. Figures, report and annotated .h5ad are simply overwritten in place.
        """
        self.paths.dir.mkdir(parents=True, exist_ok=True)
        self.paths.tool_log.unlink(missing_ok=True)
        for stale in self.paths.checkpoints.glob("*.h5ad"):
            stale.unlink()
        self._step = 0

    def require_adata(self) -> ad.AnnData:
        if self.adata is None:
            raise RuntimeError("No dataset loaded into the session.")
        return self.adata

    def checkpoint(self, label: str) -> str:
        """Persist the current AnnData to a numbered checkpoint. Returns the path."""
        adata = self.require_adata()
        self.paths.checkpoints.mkdir(parents=True, exist_ok=True)
        self._step += 1
        path = self.paths.checkpoints / f"{self._step:02d}_{label}.h5ad"
        adata.write_h5ad(path)
        return str(path)


# Single run-scoped instance the tools import.
SESSION = Session()
