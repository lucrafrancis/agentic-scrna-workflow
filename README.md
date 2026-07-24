# Agentic Single-Cell RNA-seq Analysis

This is a small project where an LLM (Claude) runs a single-cell RNA-seq analysis. Each tool
hands back a summary of what it found; the model reads that, picks the next tool and the
parameters to call it with, and writes up the results at the end. The step order and the
thresholds come from the model, not from a script.

The point was to see how much of a real analysis an LLM can drive, and to keep the mechanism
visible: the loop is ~40 lines against the raw API, not a framework.

![UMAP of the annotated PBMC3k result](examples/pbmc3k/figures/umap.png)

*Result on the PBMC3k dataset: clusters annotated with their cell types. The full write-up the
agent produced is in [`examples/pbmc3k/report.md`](examples/pbmc3k/report.md).*

Each example folder also carries the decision log for that run —
[`examples/pbmc3k/tool_calls.jsonl`](examples/pbmc3k/tool_calls.jsonl) and
[`examples/pbmc_multibatch/tool_calls.jsonl`](examples/pbmc_multibatch/tool_calls.jsonl) — one
line per tool call, with the arguments the agent chose and the summary it read back.

## What it decides

A few examples from the two runs below:

- **PCA vs. scVI.** On PBMC3k it saw there was no batch information and picked PCA. On a second
  dataset with two batches it picked scVI instead, to correct for them.
- **Adjusting to the data.** On the two-batch set, the usual 5% mitochondrial cutoff would have
  thrown away 82% of the cells (their baseline mito was just higher), so it loosened the cutoff
  to 15% instead of blindly applying the default.
- **Doublet cutoff.** Nothing is hard-coded: the tool reports the score distribution and several
  candidate cutoffs, and the agent picks one and states why in the report.

Both runs happen to follow a fairly standard QC → cluster → annotate arc, which is what these
datasets call for. The variation is in the parameters and the choice of method at each step.

## A run with batches

On a dataset made of two separate 10x runs, the agent chose scVI and integrated them. Colouring
the result by batch shows the two runs mixed together within each cell type, which is what you
want to see when integration works:

![UMAP coloured by cell type, cluster, and batch](examples/pbmc_multibatch/figures/umap.png)

Full report: [`examples/pbmc_multibatch/report.md`](examples/pbmc_multibatch/report.md).

## Running it

You'll need [`uv`](https://docs.astral.sh/uv/) and an Anthropic API key.

```bash
uv sync                                      # install dependencies

echo 'export ANTHROPIC_API_KEY=sk-ant-...' >> ~/.zshenv && source ~/.zshenv

uv run python scripts/fetch_pbmc3k.py        # download the example data
uv run python run.py data/pbmc3k.h5ad        # run the agent
```

For the two-batch example instead:

```bash
uv run python scripts/fetch_pbmc_multibatch.py
uv run python run.py data/pbmc_multibatch.h5ad
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
tests/          # fast offline tests (no API, no downloads): uv run pytest
run.py          # entry point
examples/       # saved example runs, one folder per dataset
```