"""The manual agent loop.

This is the heart of the project and intentionally framework-free: it is the ~40 lines
that turn "Claude wants to call a tool" into "the Python function ran and Claude saw the
result". Everything else (tools, schemas, prompt) plugs into this.

Flow:
  1. Send the conversation + tool schemas to Claude.
  2. If Claude requests tool calls, run each via the dispatch table, capture the summary.
  3. Append the tool results and loop.
  4. Stop when Claude responds without requesting a tool (or MAX_TURNS is hit).
"""

from __future__ import annotations

import json

import anthropic

from agent import config
from agent.prompts import SYSTEM_PROMPT
from agent.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS


def _log_tool_call(name: str, args: dict, summary: dict) -> None:
    """Append one tool call to the JSONL run log — the audit trail of the agent's reasoning."""
    config.TOOL_LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with config.TOOL_LOG_PATH.open("a") as f:
        f.write(json.dumps({"tool": name, "args": args, "summary": summary}) + "\n")


def _run_tool(name: str, args: dict) -> dict:
    """Dispatch to the registered tool, converting exceptions into an error summary the
    agent can react to rather than a crash."""
    try:
        return TOOL_FUNCTIONS[name](**args)
    except Exception as exc:  # deliberately broad: surface, don't crash the loop
        return {"error": type(exc).__name__, "message": str(exc)}


def run_agent(initial_user_message: str) -> list[dict]:
    """Drive the analysis to completion. Returns the full message transcript."""
    client = anthropic.Anthropic()
    messages: list[dict] = [{"role": "user", "content": initial_user_message}]

    for _turn in range(config.MAX_TURNS):
        response = client.messages.create(
            model=config.MODEL,
            max_tokens=config.MAX_TOKENS,
            system=SYSTEM_PROMPT,
            tools=TOOL_SCHEMAS,
            messages=messages,
        )
        messages.append({"role": "assistant", "content": response.content})

        if response.stop_reason != "tool_use":
            break  # Claude produced a final answer, no tool requested.

        tool_results = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            summary = _run_tool(block.name, block.input)
            _log_tool_call(block.name, block.input, summary)
            tool_results.append(
                {
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": json.dumps(summary),
                }
            )
        messages.append({"role": "user", "content": tool_results})

    return messages
