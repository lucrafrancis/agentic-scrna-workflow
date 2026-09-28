"""Turn a run's tool log into a standalone replay script.

The agent's decisions (thresholds, embedding, resolution, relabels, contrasts, report text)
are the arguments of its tool calls in tool_calls.jsonl. Replaying the successful calls in
order, against the same code and seed, reproduces the analysis without the LLM.

Every successful call is replayed, read-only ones too: the report is built from the tool
log (e.g. the canonical-marker dotplot shows the genes the agent checked), so the replay
writes its own log the same way. Failed calls, including rejected report drafts, are left
out.

  uv run python -m agent.replay outputs/kang   — write outputs/kang/replay.py for a finished run
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from pprint import pformat

from agent import config


def _git_commit() -> str:
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=config.ROOT, capture_output=True,
                                text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "agent"], cwd=config.ROOT,
                               capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return f"{commit}-dirty" if dirty else commit


def _format_call(name: str, args: dict) -> str:
    if not args:
        return f'_run("{name}")'
    lines = [f'_run(\n    "{name}",']
    for key, value in args.items():
        formatted = pformat(value, width=88, compact=True, sort_dicts=False).replace("\n", "\n    ")
        lines.append(f"    {key}={formatted},")
    lines.append(")")
    return "\n".join(lines)


def write_replay_script(run_dir: Path, dataset: Path) -> Path | None:
    """Write run_dir/replay.py from run_dir/tool_calls.jsonl. None if nothing succeeded."""
    log_path = run_dir / "tool_calls.jsonl"
    if not log_path.is_file():
        return None
    entries = [json.loads(line) for line in log_path.read_text().splitlines() if line.strip()]
    calls = [e for e in entries if "error" not in e["summary"]]
    if not calls:
        return None

    name = run_dir.name
    out_path = run_dir / "replay.py"
    blocks = []
    for e in calls:
        block = _format_call(e["tool"], e["args"])
        # A report the agent got accepted only after the rejection limit was written with
        # warnings; start the replay at that limit so it is written the same way.
        if e["tool"] == "generate_report" and "warning" in e["summary"]:
            block = "SESSION.report_attempts = tools._MAX_REPORT_REJECTIONS\n" + block
        blocks.append(block)

    header = f'''"""Replay of the agent's tool calls for {name}, without the LLM.

Generated: {datetime.now().isoformat(timespec="seconds")}
Git commit: {_git_commit()}
Source log: {log_path.name} ({len(calls)} of {len(entries)} calls; failed calls omitted)

Re-runs the same tool functions with the same arguments and seed. Output goes to
outputs/{name}_replay/, with its own tool log and report. Package versions are pinned by
uv.lock at the commit above: check it out first if the code has changed since.

Caveats: annotate_celltypes downloads CellTypist models, and scVI training can differ
slightly between hardware (CPU vs GPU), which can change clusters downstream.

  uv run python {out_path.relative_to(config.ROOT) if out_path.is_relative_to(config.ROOT) else out_path}
"""

import os
import sys
from pathlib import Path

# The repo this script sits in (it may have been copied, e.g. into benchmarks/runs/), else
# the repo it was generated in.
ROOT = next((p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists()),
            Path({str(config.ROOT)!r}))
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)  # the dataset path below is relative to the repo root

from agent import config, tools
from agent.loop import _log_tool_call, _run_tool
from agent.session import SESSION


def _run(name, **args):
    print(f"-> {{name}}")
    summary = _run_tool(name, args)
    _log_tool_call(name, args, summary)
    if "error" in summary:
        sys.exit(f"   FAILED: {{summary['error']}}: {{summary.get('message', '')}}")
    return summary


config.set_global_seed()
SESSION.load({str(dataset)!r})
SESSION.paths = config.RunPaths({name + "_replay"!r})
SESSION.begin_run()
print(f"Replay output: {{SESSION.paths.dir}}")

'''
    out_path.write_text(header + "\n\n".join(blocks) + '\n\nprint(f"Done. Report: {SESSION.paths.report}")\n')
    return out_path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python -m agent.replay <run_dir>   e.g. python -m agent.replay outputs/kang")
    from run import dataset_path

    run_dir = Path(sys.argv[1]).resolve()
    script = write_replay_script(run_dir, dataset_path((run_dir / "prompt.txt").read_text()))
    print(f"Wrote {script}" if script else f"Nothing to replay in {run_dir / 'tool_calls.jsonl'}")
