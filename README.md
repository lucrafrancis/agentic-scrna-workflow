# Agentic Single-Cell RNA-seq Analysis

This is a small project where an LLM (Claude) runs a single-cell RNA-seq analysis on its own.
It looks at the data, decides what to do next based on what it finds, calls the right tool,
and writes up the results at the end. There's no fixed script; the model chooses each step.

The point was to see whether an LLM could act like an analysis assistant that actually makes
sensible decisions, rather than just wrapping a chatbot around a pipeline.

![UMAP of the annotated PBMC3k result](examples/pbmc3k/umap.png)

*Result on the PBMC3k dataset: clusters annotated with their cell types. The full write-up the
agent produced is in [`examples/pbmc3k/report.md`](examples/pbmc3k/report.md).*

## What it decides

A few examples from the PBMC3k run:

- **PCA vs. scVI.** It noticed there was no batch information, so it picked PCA (scVI is only
  worth it when you need to correct for batches).
- **Doublet cutoff.** The score distribution wasn't cleanly split in two, so it used a
  `median + 3·MAD` cutoff, and even flagged that Scrublet's automatic threshold looked too loose.
- **Cell types.** The CellTypist labels matched the marker genes on every cluster.

## Running it

You'll need [`uv`](https://docs.astral.sh/uv/) and an Anthropic API key.

```bash
uv sync                                      # install dependencies

echo 'export ANTHROPIC_API_KEY=sk-ant-...' >> ~/.zshenv && source ~/.zshenv

uv run python scripts/fetch_pbmc3k.py        # download the example data
uv run python run.py data/pbmc3k.h5ad        # run the agent
```

As it runs you'll see its reasoning and the tools it calls. Results go to `outputs/`: the
report, the figures, the annotated `.h5ad`, and the decision log. A run costs a few cents.

## Layout

```
agent/
  loop.py       # the loop that talks to Claude and runs tools
  tools.py      # the analysis tools (each returns a summary)
  schemas.py    # tool descriptions Claude sees, and the name-to-function map
  prompts.py    # the instructions given to the agent
  session.py    # holds the dataset and saves checkpoints
  config.py     # model name, seed, paths
scripts/
  fetch_pbmc3k.py
run.py          # entry point
examples/       # saved example runs, one folder per dataset
CLAUDE.md       # notes on the design and why things are the way they are
```

It's a proof of concept: one dataset, and the interesting part is the decision-making rather
than the pipeline itself. More detail on the design is in [`CLAUDE.md`](CLAUDE.md).
