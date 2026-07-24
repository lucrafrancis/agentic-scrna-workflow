"""Session.load: dataset name, per-run paths, and barcode hygiene."""

from __future__ import annotations

from agent import tools
from agent.session import SESSION


def test_name_and_paths_derive_from_filename(synthetic_h5ad):
    SESSION.load(synthetic_h5ad)
    assert SESSION.name == "synthetic"
    assert SESSION.paths.report.name == "report.md"
    assert "synthetic" in str(SESSION.paths.dir)
    assert SESSION.paths.checkpoints == SESSION.paths.dir / "checkpoints"


def test_duplicate_barcodes_are_fixed_and_reported(dup_barcode_h5ad):
    SESSION.load(dup_barcode_h5ad)
    assert SESSION.adata.obs_names.is_unique
    assert SESSION.n_duplicate_barcodes > 0
    assert tools.inspect_dataset()["duplicate_barcodes_fixed"] == SESSION.n_duplicate_barcodes


def test_clean_barcodes_report_zero_fixed(synthetic_h5ad):
    SESSION.load(synthetic_h5ad)
    assert SESSION.n_duplicate_barcodes == 0
    assert tools.inspect_dataset()["duplicate_barcodes_fixed"] == 0


def test_batch_key_detected(batched_h5ad):
    SESSION.load(batched_h5ad)
    summary = tools.inspect_dataset()
    assert summary["candidate_batch_key"] == "batch"
    assert summary["n_batches"] == 2
