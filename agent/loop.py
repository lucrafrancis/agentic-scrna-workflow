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
from agent.session import SESSION


def _log_tool_call(name: str, args: dict, summary: dict) -> None:
    """Append one tool call to the JSONL run log — the audit trail of the agent's reasoning."""
    log_path = SESSION.paths.tool_log
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as f:
        f.write(json.dumps({"tool": name, "args": args, "summary": summary}) + "\n")


def _run_tool(name: str, args: dict) -> dict:
    """Dispatch to the registered tool, converting exceptions into an error summary the
    agent can react to rather than a crash."""
    try:
        return TOOL_FUNCTIONS[name](**args)
    except Exception as exc:  # deliberately broad: surface, don't crash the loop
        return {"error": type(exc).__name__, "message": str(exc)}


def _trace_text(text: str) -> None:
    print(f"\n\U0001f9e0 {text.strip()}")


def _trace_tool(name: str, args: dict, summary: dict) -> None:
    arg_str = ", ".join(f"{k}={v}" for k, v in args.items())
    print(f"\U0001f527 {name}({arg_str})")
    keys = "error" if "error" in summary else ", ".join(list(summary)[:6])
    print(f"   → {keys}")


def run_agent(initial_user_message: str, verbose: bool = True) -> list[dict]:
    """Drive the analysis to completion. Returns the full message transcript.

    With verbose=True, prints a live trace: the agent's reasoning before each tool call and
    the keys of each tool result. This is the readable record of the agent making decisions.
    """
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

        if verbose:
            for block in response.content:
                if block.type == "text" and block.text.strip():
                    _trace_text(block.text)

        if response.stop_reason != "tool_use":
            break  # Claude produced a final answer, no tool requested.

        tool_results = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            summary = _run_tool(block.name, block.input)
            _log_tool_call(block.name, block.input, summary)
            if verbose:
                _trace_tool(block.name, block.input, summary)
            tool_results.append(
                {
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": json.dumps(summary),
                }
            )
        messages.append({"role": "user", "content": tool_results})

    return messages
