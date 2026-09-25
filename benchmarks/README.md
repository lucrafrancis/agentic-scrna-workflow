# Benchmarks

Repeated runs of the agent, to measure two things:

- **Accuracy:** do the agent's cell-type labels match published labels?
- **Consistency:** does the agent make the same choices when run again on the same data?

## Runs

5 runs each on PBMC3k and the two-batch PBMC set. Each run is saved to
`runs/<dataset>/<run_id>/`:

| File | Contents |
|---|---|
| `tool_calls.jsonl` | Decision log: every tool call, its arguments, and the summary returned |
| `report.md` | The agent's write-up |
| `labels.csv` | Per cell: barcode, cluster, agent label |
| `run_meta.json` | Model, temperature, git commit, timestamp, wall time, tokens in/out, cost, number of tool calls |

## Scoring

A script reads every run folder and writes `metrics.csv`, one row per run.

- **Decisions:** PCA vs. scVI, mitochondrial cutoff, doublet cutoff, clustering resolution,
  number of clusters, and the cell types found. These are taken from `tool_calls.jsonl`.
- **Accuracy (PBMC3k only):** the reference is the published PBMC3k labels (scanpy's processed
  PBMC3k, 2,638 cells), matched by barcode. Scores are the adjusted Rand index (ARI, a 0–1
  measure of how closely two groupings of the same cells agree) and per-cell-type agreement.
- **Label mapping:** agent and reference names differ (e.g. "Classical monocytes" vs.
  "CD14+ Monocytes"), so both are mapped to a fixed set of broad cell types. The mapping is
  committed before any benchmark run, so it cannot be tuned to the results.

The two-batch set has no reference labels, so it is scored on consistency only.
