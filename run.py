"""Entrypoint: uv run python run.py data/pbmc3k.h5ad

Loads the dataset, seeds everything, hands control to the agent loop, and leaves the
outputs (report, figures, annotated .h5ad, tool-call log) in outputs/.
"""

from __future__ import annotations

import sys
from pathlib import Path

from agent import config
from agent.loop import run_agent


def main(h5ad_path: str) -> None:
    config.set_global_seed()
    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # TODO: load the AnnData and make it available to the tools (session state).
    #       Decide the sharing mechanism (module-level session object vs. explicit passing)
    #       when we implement the first real tool.

    prompt = (
        f"Analyze the single-cell dataset at {h5ad_path}. "
        "Inspect it first, then proceed through an appropriate analysis and produce a report."
    )
    run_agent(prompt)
    print(f"Done. See {config.REPORT_PATH}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python run.py <path-to.h5ad>")
    main(sys.argv[1])
