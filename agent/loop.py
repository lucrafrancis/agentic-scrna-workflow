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


_USAGE_KEYS = ("input_tokens", "cache_creation_input_tokens", "cache_read_input_tokens", "output_tokens")


def _estimate_cost(totals: dict[str, int]) -> float | None:
    price = config.PRICE_PER_MTOK.get(config.MODEL)
    if price is None:
        return None
    return (
        totals["input_tokens"] * price["input"]
        + totals["cache_creation_input_tokens"] * price["input"] * config.CACHE_WRITE_MULTIPLIER
        + totals["cache_read_input_tokens"] * price["input"] * config.CACHE_READ_MULTIPLIER
        + totals["output_tokens"] * price["output"]
    ) / 1e6


def _log_usage(totals: dict[str, int]) -> None:
    """Print the run's token totals and estimated cost, and write them to usage.jsonl."""
    cost = _estimate_cost(totals)
    total_in = sum(totals[k] for k in _USAGE_KEYS[:3])
    cached = 100 * totals["cache_read_input_tokens"] / total_in if total_in else 0
    cost_str = f", ~${cost:.2f}" if cost is not None else ""
    print(f"\n\U0001f4ca {totals['requests']} requests, {total_in:,} input tokens ({cached:.0f}% from cache), "
          f"{totals['output_tokens']:,} output tokens{cost_str}")
    record = {"model": config.MODEL, **totals, "estimated_cost_usd": round(cost, 4) if cost is not None else None}
    SESSION.paths.usage_log.write_text(json.dumps(record) + "\n")


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
    totals = {"requests": 0, **{k: 0 for k in _USAGE_KEYS}}

    for _turn in range(config.MAX_TURNS):
        response = client.messages.create(
            model=config.MODEL,
            max_tokens=config.MAX_TOKENS,
            system=SYSTEM_PROMPT,
            tools=TOOL_SCHEMAS,
            messages=messages,
            # Cache the conversation so far: each turn resends it, and only the new tail is
            # billed at the full input price.
            cache_control={"type": "ephemeral"},
        )
        totals["requests"] += 1
        for key in _USAGE_KEYS:
            totals[key] += getattr(response.usage, key, None) or 0
        messages.append({"role": "assistant", "content": response.content})

        if verbose:
            for block in response.content:
                if block.type == "text" and block.text.strip():
                    _trace_text(block.text)

        if response.stop_reason != "tool_use":
            # "end_turn" is the agent finishing normally. Anything else (notably
            # "max_tokens", which truncates a tool_use block mid-generation) is the run
            # being cut short, and must not look like a clean finish.
            if response.stop_reason != "end_turn":
                print(f"\n⚠️  Run stopped early: stop_reason={response.stop_reason!r}.")
                if response.stop_reason == "max_tokens":
                    print(f"   The response hit MAX_TOKENS ({config.MAX_TOKENS}); raise it in agent/config.py.")
            break  # Claude produced a final answer, or was cut short.

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

    _log_usage(totals)
    return messages
