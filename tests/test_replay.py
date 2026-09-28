"""The replay script re-runs a logged analysis without the LLM and gets the same result."""

from __future__ import annotations

import json
import runpy

from agent import config
from agent.replay import write_replay_script
from agent.session import SESSION


def test_replay_reproduces_the_run(logged_run, synthetic_h5ad, monkeypatch):
    logged_run("check_markers", genes=["GENE1", "GENE2"])
    logged_run("summarize_findings")
    original_dir = SESSION.paths.dir
    original_clusters = SESSION.adata.obs["leiden"].astype(str).tolist()
    # A failed call is left out of the replay.
    logged_run("cluster", resolution="not a number")

    script = write_replay_script(original_dir, synthetic_h5ad)
    text = script.read_text()
    assert text.count('_run(\n    "cluster"') == 1

    monkeypatch.chdir(config.ROOT)  # the script changes directory; restore it afterwards
    runpy.run_path(str(script), run_name="__main__")

    assert SESSION.paths.dir == config.OUTPUT_DIR / f"{original_dir.name}_replay"
    assert SESSION.adata.obs["leiden"].astype(str).tolist() == original_clusters
    replayed = [json.loads(line) for line in SESSION.paths.tool_log.read_text().splitlines()]
    original = [json.loads(line) for line in (original_dir / "tool_calls.jsonl").read_text().splitlines()]
    succeeded = [e for e in original if "error" not in e["summary"]]
    assert [(e["tool"], e["args"]) for e in replayed] == [(e["tool"], e["args"]) for e in succeeded]


def test_nothing_to_replay(tmp_path, synthetic_h5ad):
    assert write_replay_script(tmp_path, synthetic_h5ad) is None
    (tmp_path / "tool_calls.jsonl").write_text(json.dumps({"tool": "cluster", "args": {}, "summary": {"error": "x"}}) + "\n")
    assert write_replay_script(tmp_path, synthetic_h5ad) is None
