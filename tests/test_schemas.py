"""The tool schemas Claude sees must stay consistent with the dispatch table."""

from __future__ import annotations

from agent.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS


def test_every_schema_has_a_callable_function():
    for schema in TOOL_SCHEMAS:
        name = schema["name"]
        assert name in TOOL_FUNCTIONS, f"{name} exposed to Claude but not dispatchable"
        assert callable(TOOL_FUNCTIONS[name])


def test_schemas_are_well_formed():
    for schema in TOOL_SCHEMAS:
        assert {"name", "description", "input_schema"} <= schema.keys()
        assert schema["input_schema"]["type"] == "object"


def test_no_duplicate_tool_names():
    names = [s["name"] for s in TOOL_SCHEMAS]
    assert len(names) == len(set(names))
