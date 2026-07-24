"""Central configuration. Everything tunable lives here so nothing else has to change.

The model name in particular is intentionally a single constant: swapping models is a
one-line edit and nothing else in the codebase depends on it.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# --- LLM ---------------------------------------------------------------------
# Switch to "claude-opus-5" if the agent's *reasoning* (tool choice) looks weak.
MODEL = "claude-sonnet-5"
# Headroom for the largest single generation in a run: generate_report's narrative, passed
# as one tool argument (~5k characters in the committed examples). Too low a ceiling
# truncates that tool_use block, and the run ends with stop_reason="max_tokens" instead of
# a report.
MAX_TOKENS = 16384
# Safety valve so a misbehaving agent can't loop forever.
MAX_TURNS = 40

# --- Reproducibility ---------------------------------------------------------
SEED = 0

# --- Paths -------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "outputs"
EXAMPLES_DIR = ROOT / "examples"


@dataclass(frozen=True)
class RunPaths:
    """Per-dataset output locations, namespaced by dataset name (the input filename stem).

    A run on data/pbmc3k.h5ad writes everything under outputs/pbmc3k/, so runs on different
    datasets never collide. The Session derives the name and holds an instance of this.
    """

    name: str

    @property
    def dir(self) -> Path:
        return OUTPUT_DIR / self.name

    @property
    def figures(self) -> Path:
        return self.dir / "figures"

    @property
    def checkpoints(self) -> Path:
        return self.dir / "checkpoints"

    @property
    def report(self) -> Path:
        return self.dir / "report.md"

    @property
    def annotated(self) -> Path:
        return self.dir / "annotated.h5ad"

    @property
    def tool_log(self) -> Path:
        return self.dir / "tool_calls.jsonl"


def set_global_seed(seed: int = SEED) -> None:
    """Seed every stochastic library we use, so each computation is reproducible.

    Called once at the start of a run. scvi-tools is seeded lazily inside its tool to avoid
    importing torch until it is actually needed. Note that this fixes the numerics, not the
    analysis path: the agent's tool choices are sampled from the model and can vary between
    runs on the same input.
    """
    import random

    import numpy as np
    import scanpy as sc

    random.seed(seed)
    np.random.seed(seed)
    sc.settings.seed = seed
