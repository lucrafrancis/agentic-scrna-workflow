"""Entrypoint: uv run python run.py data/pbmc3k.h5ad

Loads the dataset, seeds everything, hands control to the agent loop, and leaves the
outputs (report, figures, annotated .h5ad, tool-call log) in outputs/.
"""

from __future__ import annotations

import sys
from pathlib import Path

from agent import config
from agent.loop import run_agent
from agent.session import SESSION


def main(h5ad_path: str) -> None:
    config.set_global_seed()

    # Load the dataset into the shared session; tools reach it via SESSION (see session.py).
    # The session derives the run name from the filename, so outputs land in outputs/<name>/.
    SESSION.load(h5ad_path)
    SESSION.paths.dir.mkdir(parents=True, exist_ok=True)

    prompt = (
        f"Analyze the single-cell dataset at {h5ad_path}. "
        "Inspect it first, then proceed through an appropriate analysis and produce a report."
    )
    run_agent(prompt)
    if SESSION.paths.report.exists():
        print(f"Done. Report: {SESSION.paths.report}")
    else:
        print("Done, but no report.md was generated (the agent did not call generate_report).")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python run.py <path-to.h5ad>")
    main(sys.argv[1])
