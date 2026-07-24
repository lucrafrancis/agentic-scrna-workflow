"""Central configuration. Everything tunable lives here so nothing else has to change.

The model name in particular is intentionally a single constant: swapping models is a
one-line edit and nothing else in the codebase depends on it.
"""

from __future__ import annotations

from pathlib import Path

# --- LLM ---------------------------------------------------------------------
# Switch to "claude-opus-4-8" if the agent's *reasoning* (tool choice) looks weak.
MODEL = "claude-sonnet-5"
MAX_TOKENS = 4096
# Safety valve so a misbehaving agent can't loop forever.
MAX_TURNS = 40

# --- Reproducibility ---------------------------------------------------------
SEED = 0

# --- Paths -------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
OUTPUT_DIR = ROOT / "outputs"
FIGURE_DIR = OUTPUT_DIR / "figures"
REPORT_PATH = OUTPUT_DIR / "report.md"
TOOL_LOG_PATH = OUTPUT_DIR / "tool_calls.jsonl"


def set_global_seed(seed: int = SEED) -> None:
    """Seed every stochastic library we use, so a run is reproducible.

    Called once at the start of a run. scvi-tools is seeded lazily inside its tool
    to avoid importing torch until it is actually needed.
    """
    import random

    import numpy as np
    import scanpy as sc

    random.seed(seed)
    np.random.seed(seed)
    sc.settings.seed = seed
