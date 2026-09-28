# Benchmarks

Repeated runs of the agent on one prompt, to measure two things:

- **Accuracy:** do the agent's cell-type labels match published labels?
- **Consistency:** does the agent make the same choices when run again on the same data?

## Running

    uv run python benchmarks/run_benchmark.py prompts/pbmc3k.txt --runs 5   # spends API tokens
    uv run python benchmarks/score.py pbmc3k                                # offline

Each run is saved to `runs/<dataset>/<run_id>/`: the decision log (`tool_calls.jsonl`), the
report, the prompt, token use and estimated cost (`usage.jsonl`), per-cell labels
(`labels.csv`: cluster, final label, original CellTypist label), `run_meta.json` (model,
git commit, runtime), the run's console output (`trace.log`), and `replay.py`, which re-runs
the logged tool calls without the model.

## Scoring

`score.py` writes to `results/<dataset>/`:

- `metrics.csv`: one row per run. The decisions (embedding, mitochondrial cutoff, doublet
  threshold and % removed, number of clusters and cell types, clusters relabelled), then cells
  scored, % of reference cells removed, accuracy, ARI, report attempts, cost and runtime.
- `summary.md`: how many runs made each decision, and mean accuracy, ARI and cost.
- `confusion.png`: where each published cell type ended up, pooled over runs: an agent type,
  "Removed" (dropped by the agent's QC or doublet filter) or "Not scored" (ambiguous label).

**Reference labels** are the published labels (PBMC3k: scanpy's processed PBMC3k, 2,638 cells
from the Seurat tutorial; Kang: the authors' labels). They are one analyst's annotation, not
ground truth.

**Label mapping** (`label_mapping.csv`): agent and published names differ (e.g. "Classical
monocytes" vs "CD14+ Monocytes"), so both are mapped to shared broad types. A label mapped to
an empty broad type is ambiguous (e.g. MAIT cells, neither CD4 nor CD8 T) and its cells are
not scored. A label missing from the mapping stops scoring rather than being guessed.

**Accuracy** is the % of scored cells (kept by the agent, unambiguous label) whose broad type
matches the published one. **ARI** (adjusted Rand index) measures how closely the two
groupings of the same cells agree, from 0 (chance) to 1 (identical).

## Acceptance thresholds

Set before the benchmark runs. The same thresholds apply to PBMC3k and Kang, except cost:
Kang is about 9× larger and adds composition and DE steps.

| Metric | Acceptable |
|---|---|
| Accuracy | ≥ 90% |
| ARI | ≥ 0.80 |
| Embedding and annotation model | same in 5 of 5 runs |
| Mitochondrial cutoff, doublet threshold, number of clusters | same in ≥ 4 of 5 runs |
| Report attempts | ≤ 2 |
| Cost | < $0.50 per run (PBMC3k), < $1 per run (Kang) |
