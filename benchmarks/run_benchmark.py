"""Run the agent repeatedly on one prompt and keep each run's outputs for scoring.

    uv run python benchmarks/run_benchmark.py prompts/pbmc3k.txt --runs 5

Each run calls run.py as a fresh process (so no state carries over), then copies what the
scoring needs from outputs/<dataset>/ into benchmarks/runs/<dataset>/<run_id>/:

    tool_calls.jsonl   the decision log
    report.md          the agent's report
    prompt.txt         what the agent was asked
    usage.jsonl        tokens and estimated cost
    labels.csv         per cell: barcode, Leiden cluster, final label, CellTypist label
    run_meta.json      model, git commit, start time, wall time, exit code

Every run spends API tokens.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import anndata as ad

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from agent import config  # noqa: E402
from run import dataset_path  # noqa: E402

RUNS_DIR = ROOT / "benchmarks" / "runs"
COPIED = ("tool_calls.jsonl", "report.md", "prompt.txt", "usage.jsonl")


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True).stdout.strip()


def _save_labels(annotated: Path, out: Path) -> None:
    obs = ad.read_h5ad(annotated, backed="r").obs
    cols = [c for c in ("leiden", "cell_type", "cell_type_celltypist") if c in obs]
    obs[cols].to_csv(out, index_label="barcode")


def run_once(prompt_file: Path, dataset: str) -> Path:
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    dest = RUNS_DIR / dataset / run_id
    dest.mkdir(parents=True)
    commit, dirty = _git("rev-parse", "--short", "HEAD"), bool(_git("status", "--porcelain", "--", "agent", "run.py"))

    start = time.time()
    log = (dest / "trace.log").open("w")
    proc = subprocess.run(["uv", "run", "python", "run.py", str(prompt_file)], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    log.close()
    wall = time.time() - start

    # Only files this run wrote: a failed run must not pick up the previous run's outputs.
    src = config.OUTPUT_DIR / dataset

    def fresh(path: Path) -> bool:
        return path.exists() and path.stat().st_mtime >= start

    for name in COPIED:
        if fresh(src / name):
            shutil.copy(src / name, dest / name)
    if fresh(src / "annotated.h5ad"):
        _save_labels(src / "annotated.h5ad", dest / "labels.csv")
    meta = {
        "run_id": run_id,
        "dataset": dataset,
        "prompt": str(prompt_file),
        "model": config.MODEL,
        "git_commit": commit,
        "agent_code_modified": dirty,
        "started": datetime.fromtimestamp(start).isoformat(timespec="seconds"),
        "wall_time_s": round(wall, 1),
        "exit_code": proc.returncode,
        "report_written": (dest / "report.md").exists(),
    }
    (dest / "run_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return dest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("prompt", type=Path)
    parser.add_argument("--runs", type=int, default=5)
    args = parser.parse_args()

    dataset = dataset_path(args.prompt.read_text()).stem
    for i in range(args.runs):
        dest = run_once(args.prompt, dataset)
        usage = json.loads((dest / "usage.jsonl").read_text()) if (dest / "usage.jsonl").exists() else {}
        meta = json.loads((dest / "run_meta.json").read_text())
        print(f"[{i + 1}/{args.runs}] {dest.name}: {meta['wall_time_s']:.0f}s, "
              f"~${usage.get('estimated_cost_usd', float('nan')):.2f}, report={meta['report_written']}")


if __name__ == "__main__":
    main()
