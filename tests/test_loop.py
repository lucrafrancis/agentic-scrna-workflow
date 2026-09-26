"""The agent loop, with a fake API client: no network, no tokens spent."""

from __future__ import annotations

import json
from types import SimpleNamespace

from agent import config, loop
from agent.session import SESSION


class _FakeClient:
    """Answers once, with no tool call, and records what it was sent."""

    def __init__(self):
        self.calls = []
        self.messages = self

    def create(self, **kwargs):
        self.calls.append(kwargs)
        return SimpleNamespace(
            content=[SimpleNamespace(type="text", text="Done.")],
            stop_reason="end_turn",
            usage=SimpleNamespace(input_tokens=1000, cache_creation_input_tokens=4000,
                                  cache_read_input_tokens=20000, output_tokens=500),
        )


def test_usage_is_logged_and_requests_are_cached(synthetic_h5ad, monkeypatch):
    fake = _FakeClient()
    monkeypatch.setattr(loop.anthropic, "Anthropic", lambda: fake)
    SESSION.load(synthetic_h5ad)
    SESSION.begin_run()

    loop.run_agent("Analyze it.", verbose=False)

    assert fake.calls[0]["cache_control"] == {"type": "ephemeral"}
    record = json.loads(SESSION.paths.usage_log.read_text())
    assert record["requests"] == 1
    assert record["cache_read_input_tokens"] == 20000
    price = config.PRICE_PER_MTOK[config.MODEL]
    expected = (1000 * price["input"] + 4000 * price["input"] * 1.25
                + 20000 * price["input"] * 0.1 + 500 * price["output"]) / 1e6
    assert abs(record["estimated_cost_usd"] - round(expected, 4)) < 1e-9
