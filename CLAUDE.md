# Project: Agentic Single-Cell Analysis Workflow

## Goal

Build a lightweight proof-of-concept demonstrating an **agentic LLM workflow** for
single-cell RNA-seq analysis. The objective is **not** a deterministic pipeline — it is
to showcase an LLM acting as an autonomous analysis assistant that selects and executes
computational tools based on intermediate results, and reasons about what to do next.

The workflow starts from a preprocessed `AnnData (.h5ad)` object (not FASTQs), to focus on
the biological analysis where intelligent decision-making adds value. This is a portfolio
piece: the value is in *demonstrating* that an LLM can drive a real analysis, so favour
clarity and transparency of the agentic mechanism over feature completeness.

## Architecture

- **Manual agent loop on the raw `anthropic` Python SDK.** We deliberately do *not* use a
  framework (Agent SDK, LangChain). The loop is ~40 lines and lives in our code so the
  tool-calling mechanism is fully visible and explainable.
- **The core contract:** every tool returns a **structured, JSON-serializable summary**
  that the LLM reads to decide the next step. Side effects (writing `.h5ad`, saving
  figures) are secondary. A tool that returns `None` is a bug — the agent is then blind.
- Three separate pieces of code, do not conflate them:
  1. the **loop** (`agent/loop.py`) — sends messages + tool schemas to Claude, executes
     the tool calls Claude requests, feeds results back until Claude stops.
  2. the **tools** (`agent/tools.py`) — deterministic Python functions doing the analysis.
  3. the **runtime system prompt** (`agent/prompts.py`) — the instructions given to the
     *running* agent.

### CLAUDE.md vs. the runtime system prompt

Do not conflate them. **This file (`CLAUDE.md`)** steers Claude Code while *building* the
repo. The **runtime system prompt** in `agent/prompts.py` steers the LLM that *runs the
analysis* at execution time. Different audiences, different files.

## Repo layout

```
agentic-scrna-workflow/
├── agent/
│   ├── loop.py          # the manual agent loop
│   ├── tools.py         # tool functions (each returns a summary dict)
│   ├── schemas.py       # JSON tool definitions sent to Claude
│   ├── prompts.py       # the runtime system prompt
│   └── config.py        # model name, seeds, paths
├── data/                # input .h5ad (gitignored)
├── outputs/             # figures, annotated .h5ad, report.md
├── run.py               # entrypoint: python run.py data/pbmc3k.h5ad
└── pyproject.toml
```

## Tech stack

- **Env / packaging:** `uv` with `pyproject.toml` + `uv.lock` (reproducible). Run things
  with `uv run ...`.
- **Analysis:** `scanpy`, `anndata`, `scvi-tools` (scVI), `celltypist`, `scrublet`
  (doublets), `matplotlib`.
- **LLM:** `anthropic`.
- CPU-only is fine for pbmc3k; scVI's torch backend does not need a GPU here.

## Runtime model + reproducibility

- Model is a single constant in `agent/config.py`: `MODEL = "claude-sonnet-5"`. Switching
  models is a one-line change — nothing else depends on it. Bump to Opus only if the
  agent's *reasoning* looks weak.
- Seed everything (`numpy`, `scanpy`, `scvi`) from a single seed in `config.py`.
- **Log every tool call** (name + args + returned summary) so a run is reproducible and
  auditable. The log is part of the deliverable — it *is* the evidence the agent reasoned.

## Tool set

Each tool returns a structured summary dict. The agent decides order and arguments.

1. Inspect dataset (n cells/genes, organism, obs/var columns, batch key if any)
2. Compute QC metrics
3. Recommend QC filtering thresholds
4. Filter low-quality cells and genes
5. Detect doublets (scrublet)
6. Normalize & preprocess (log-normalize; HVGs; **stash raw counts first**)
7. Dimensionality reduction — **two sibling tools, the agent chooses:**
   - `run_pca` — standard, correct for a single clean batch
   - `run_scvi` — VAE latent with batch correction; reads raw counts
8. Cluster cells (neighbors + Leiden on the chosen representation)
9. Identify marker genes
10. Annotate cell types (CellTypist)
11. Summarize biological findings
12. Generate a Markdown analysis report

## Conventions / guardrails

Biological invariants the agent (and tools) must respect regardless of tool order:

- **Stash raw counts before normalizing:** `adata.layers["counts"] = adata.X.copy()`
  before any log-normalization. scVI reads `layers["counts"]`; PCA/clustering read the
  normalized matrix.
- **Never run PCA/clustering on raw counts**, and **never run scVI on log-normalized
  data.** Tools validate the state they expect and error loudly otherwise.
- **QC and filtering before doublet detection**, doublet detection before normalization.
- **Dimensionality-reduction decision rule** (encode in the system prompt): if `obs` has a
  batch key with >1 batch, prefer `run_scvi` for its correction; for a single clean batch
  (e.g. pbmc3k), `run_pca` is the correct, simpler choice. The agent picking PCA on pbmc3k
  is *good judgment*, not a missed feature.
- Tools should be defensive about the `adata` state they receive and return a clear error
  summary (not a stack trace) the agent can react to.

## Deliverables & definition of done

- Agentic workflow in Python (manual loop, tool-based).
- Auto-generated figures (QC plots, UMAP).
- Annotated `.h5ad` output.
- Markdown report + a tool-call log for the run.

**Done:** a single `uv run python run.py data/pbmc3k.h5ad` produces `outputs/report.md`
with QC, a UMAP, clusters, and CellTypist annotations — where the sequence of analysis
steps was chosen by the agent's own tool calls, not hard-coded.

## Why this project exists

Demonstrates experience building **LLM-powered agentic analysis workflows** for
computational biology, aligning with AI-driven research environments in industry. Keep the
agentic mechanism transparent and explainable — that is the thing being demonstrated.
