"""Score benchmark runs against published labels.

    uv run python benchmarks/score.py pbmc3k

Reads every run in benchmarks/runs/<dataset>/, the published labels in
data/<dataset>_reference_labels.csv, and the label mapping in benchmarks/label_mapping.csv,
and writes to benchmarks/results/<dataset>/:

    metrics.csv    one row per run: the decisions, accuracy, ARI, cost
    summary.md     how often each decision was made, and mean accuracy and cost
    confusion.png  published vs agent broad type, pooled over runs

Labels are compared as broad types via the mapping file (columns: label, broad_type). A label
mapped to an empty broad type is ambiguous (e.g. MAIT cells, neither CD4 nor CD8) and its
cells are left out of scoring. A label missing from the mapping stops the script.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "benchmarks"

DECISIONS = ["embedding", "annotation_model", "mito_cutoff", "doublet_threshold", "pct_doublets_removed",
             "n_clusters", "n_cell_types", "n_relabelled"]


def _last(log: list[dict], tool: str) -> dict | None:
    for entry in reversed(log):
        if entry["tool"] == tool and "error" not in entry["summary"]:
            return entry
    return None


def decisions(log: list[dict], labels: pd.DataFrame) -> dict:
    filt, dbl, clus = _last(log, "filter_cells_and_genes"), _last(log, "filter_doublets"), _last(log, "cluster")
    ann = _last(log, "annotate_celltypes")
    return {
        "embedding": "scVI" if _last(log, "run_scvi") else "PCA",
        "annotation_model": (ann["args"].get("model") or "Immune_All_Low.pkl").removesuffix(".pkl") if ann else None,
        "mito_cutoff": filt["args"]["max_pct_mt"] if filt else None,
        "doublet_threshold": dbl["args"]["threshold"] if dbl else None,
        "pct_doublets_removed": dbl["summary"]["pct_removed"] if dbl else 0.0,
        "n_clusters": clus["summary"]["n_clusters"] if clus else None,
        "n_cell_types": labels["cell_type"].nunique(),
        "n_relabelled": sum(len(e["summary"]["changes"]) for e in log
                            if e["tool"] == "relabel_clusters" and "error" not in e["summary"]),
    }


def to_broad(labels: pd.Series, mapping: dict[str, str], source: str) -> pd.Series:
    unknown = sorted(set(labels) - set(mapping))
    if unknown:
        raise SystemExit(f"{source} labels missing from label_mapping.csv: {unknown}")
    return labels.map(mapping)


def score_run(run: Path, reference: pd.Series, mapping: dict[str, str]) -> tuple[dict, pd.DataFrame]:
    log = [json.loads(line) for line in (run / "tool_calls.jsonl").read_text().splitlines() if line.strip()]
    labels = pd.read_csv(run / "labels.csv", index_col="barcode")
    usage = json.loads((run / "usage.jsonl").read_text()) if (run / "usage.jsonl").exists() else {}
    meta = json.loads((run / "run_meta.json").read_text()) if (run / "run_meta.json").exists() else {}

    # Every reference cell ends up somewhere: an agent broad type, "Removed" (the agent's QC
    # or doublet filter dropped it) or "Not scored" (its agent label is ambiguous).
    ref_all = to_broad(reference, mapping, "Reference")
    ref_all = ref_all[ref_all != ""]
    kept = labels.index.intersection(ref_all.index)
    outcome = pd.Series("Removed", index=ref_all.index)
    agent_kept = to_broad(labels.loc[kept, "cell_type"].astype(str), mapping, f"Agent ({run.name})")
    outcome[kept] = agent_kept.replace("", "Not scored")
    scored = outcome[~outcome.isin(["Removed", "Not scored"])].index
    agent, ref = outcome[scored], ref_all[scored]

    row = {
        "run_id": run.name,
        **decisions(log, labels),
        "cells_scored": f"{len(agent)} of {len(ref_all)}",
        "pct_removed": round(100 * float((outcome == "Removed").mean()), 1),
        "accuracy": round(100 * float((agent == ref).mean()), 1),
        "ari": round(adjusted_rand_score(ref, agent), 3),
        "report_attempts": sum(e["tool"] == "generate_report" for e in log),
        "cost_usd": usage.get("estimated_cost_usd"),
        "runtime_s": meta.get("wall_time_s"),
    }
    return row, pd.DataFrame({"reference": ref_all, "agent": outcome})


def summary(metrics: pd.DataFrame) -> str:
    n = len(metrics)
    lines = [f"# Benchmark summary ({n} runs)", "", "| Decision | Most common | Runs |", "|---|---|---|"]
    for col in DECISIONS:
        value, count = Counter(metrics[col].astype(str)).most_common(1)[0]
        lines.append(f"| {col} | {value} | {count} of {n} |")
    cost = metrics["cost_usd"].dropna()
    lines += ["", "| Measure | Mean | Range |", "|---|---|---|"]
    for col, unit in (("accuracy", "%"), ("ari", ""), ("pct_removed", "%"), ("report_attempts", "")):
        lines.append(f"| {col} | {metrics[col].mean():.3g}{unit} | {metrics[col].min():.3g}–{metrics[col].max():.3g}{unit} |")
    if len(cost):
        lines.append(f"| cost_usd | ${cost.mean():.2f} | ${cost.min():.2f}–${cost.max():.2f} (total ${cost.sum():.2f}) |")
    return "\n".join(lines) + "\n"


def confusion_figure(pairs: pd.DataFrame, n_runs: int, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    table = pd.crosstab(pairs["reference"], pairs["agent"])
    types = sorted(set(table.index) | (set(table.columns) - {"Removed", "Not scored"}))
    cols = types + [c for c in ("Not scored", "Removed") if c in table.columns]
    table = table.reindex(index=types, columns=cols, fill_value=0)
    pct = 100 * table.div(table.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    fig, ax = plt.subplots(figsize=(0.8 * len(cols) + 2, 0.7 * len(types) + 1.5))
    ax.imshow(pct.to_numpy(), cmap="Blues", vmin=0, vmax=100)
    if len(cols) > len(types):  # separate the outcome columns from the agent types
        ax.axvline(len(types) - 0.5, color="#333333", linewidth=1)
    for i in range(len(types)):
        for j in range(len(cols)):
            if pct.iat[i, j] >= 0.5:
                ax.text(j, i, f"{pct.iat[i, j]:.0f}%", ha="center", va="center", fontsize=8,
                        color="white" if pct.iat[i, j] > 60 else "#333333")
    ax.set_xticks(range(len(cols)), cols, rotation=40, ha="right")
    ax.set_yticks(range(len(types)), types)
    ax.set_xlabel("Agent label (broad type), or what happened to the cell")
    ax.set_ylabel("Published label (broad type)")
    ax.set_title(f"Agent vs published labels, {n_runs} runs pooled\n(% of each published type)", fontsize=10)
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dataset")
    parser.add_argument("--runs-dir", type=Path, help="default: benchmarks/runs/<dataset>")
    parser.add_argument("--mapping", type=Path, default=BENCH / "label_mapping.csv")
    parser.add_argument("--out", type=Path, help="default: benchmarks/results/<dataset>")
    args = parser.parse_args()

    runs_dir = args.runs_dir or BENCH / "runs" / args.dataset
    out = args.out or BENCH / "results" / args.dataset
    reference = pd.read_csv(ROOT / "data" / f"{args.dataset}_reference_labels.csv", index_col="barcode")["reference_cell_type"].astype(str)
    mapping_df = pd.read_csv(args.mapping, keep_default_na=False)
    mapping = dict(zip(mapping_df["label"], mapping_df["broad_type"]))

    runs = sorted(p for p in runs_dir.iterdir() if (p / "labels.csv").exists())
    if not runs:
        raise SystemExit(f"No scored runs in {runs_dir}")
    rows, pairs = zip(*(score_run(r, reference, mapping) for r in runs))

    out.mkdir(parents=True, exist_ok=True)
    metrics = pd.DataFrame(rows)
    metrics.to_csv(out / "metrics.csv", index=False)
    (out / "summary.md").write_text(summary(metrics))
    confusion_figure(pd.concat(pairs), len(runs), out / "confusion.png")
    print(metrics.to_string(index=False))
    print(f"\nWrote {out}/metrics.csv, summary.md, confusion.png")


if __name__ == "__main__":
    main()
