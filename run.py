"""Entrypoint: uv run python run.py prompts/pbmc3k.txt

The prompt file is the request in plain language: which dataset (a path to an .h5ad in the
text) and what to do with it. It is copied into the run's output folder, so every run records
exactly what it was asked. Loads the dataset, seeds everything, hands control to the agent
loop, and leaves the outputs (report, figures, annotated .h5ad, tool-call log) in outputs/.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

from agent import config
from agent.loop import run_agent
from agent.session import SESSION


def dataset_path(prompt: str) -> Path:
    """The .h5ad the prompt names. Exactly one is required."""
    paths = set(re.findall(r"[\w./-]+\.h5ad\b", prompt))
    if len(paths) != 1:
        sys.exit(f"The prompt must name exactly one .h5ad dataset; found {sorted(paths) or 'none'}.")
    return Path(paths.pop())


def main(prompt_file: str) -> None:
    prompt = Path(prompt_file).read_text().strip()
    h5ad_path = dataset_path(prompt)
    if not h5ad_path.exists():
        sys.exit(f"Dataset not found: {h5ad_path}. See scripts/ for the fetch scripts.")
    config.set_global_seed()

    # Load the dataset into the shared session; tools reach it via SESSION (see session.py).
    # The session derives the run name from the filename, so outputs land in outputs/<name>/.
    SESSION.load(h5ad_path)
    # Fresh run: the tool log and checkpoints must describe this run only (see begin_run).
    SESSION.begin_run()
    shutil.copyfile(prompt_file, SESSION.paths.dir / "prompt.txt")

    run_agent(prompt + "\n\nInspect the dataset first, then proceed through an appropriate analysis.")
    if SESSION.paths.report.exists():
        print(f"Done. Report: {SESSION.paths.report}")
    else:
        print(
            f"Finished without a report at {SESSION.paths.report}. See the trace above: the run "
            "either stopped early (stop_reason) or generate_report returned an error."
        )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python run.py <prompt.txt>   e.g. python run.py prompts/pbmc3k.txt")
    main(sys.argv[1])
