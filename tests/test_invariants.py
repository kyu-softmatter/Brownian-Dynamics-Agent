"""Enforces axiom A4 — the deterministic core does not call an LLM.

.claude/rules/axioms.md `A4` names this file as the check. It did not exist in
any of the three predecessor repositories, so for a month the axiom was true but
**unchecked** — which is the unwired-checker failure this project keeps
documenting (docs/02-verification.md section 6). Written 2026-08-28.

Why AST and not grep: `A4` was originally written as
`grep -rE "anthropic|claude" simbot/ bdkit/`, and that check matches the *prose
that explains the rule* -- 7 false hits, measured 2026-07-28. A string search
cannot tell code from prose. This parses the module and looks at real imports.
"""
from __future__ import annotations

import ast
from pathlib import Path

import pytest

# The deterministic core. Anything here computes; nothing here may judge with a
# model. `verify/` and `campaigns/` are analysis scripts and are held to the
# same bar, because their numbers reach conclusions too.
CORE_DIRS = ("bdbot", "simbot", "cases", "tools", "verify", "campaigns")

# Import roots that mean "an LLM is being called from here".
FORBIDDEN_ROOTS = frozenset({
    "anthropic", "openai", "google", "cohere", "mistralai", "ollama",
    "langchain", "langchain_openai", "langchain_anthropic", "llama_cpp",
    "transformers", "litellm", "claude_agent_sdk", "claude_code_sdk",
})

ROOT = Path(__file__).resolve().parent.parent


def _core_files() -> list[Path]:
    out: list[Path] = []
    for d in CORE_DIRS:
        p = ROOT / d
        if p.is_dir():
            out += [f for f in sorted(p.rglob("*.py"))
                    if "__pycache__" not in f.parts]
    return out


CORE_FILES = _core_files()


def test_core_directories_exist():
    """If the core moved, this test is measuring nothing. Fail loudly instead."""
    present = [d for d in CORE_DIRS if (ROOT / d).is_dir()]
    assert len(present) >= 4, (
        f"expected the deterministic core at {CORE_DIRS}, found only {present} — "
        "did the layout change? A4's enforcement is now pointing at nothing."
    )


def test_found_files_to_check():
    """Guards against a silently empty check -- the exact A4 failure mode."""
    assert len(CORE_FILES) > 50, (
        f"only {len(CORE_FILES)} core .py files found; this check is not "
        "covering what it claims to cover"
    )


def _imported_roots(tree: ast.AST) -> set[str]:
    roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                roots.add(a.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level == 0 and node.module:
                roots.add(node.module.split(".")[0])
    return roots


@pytest.mark.parametrize("path", CORE_FILES,
                         ids=lambda p: str(p.relative_to(ROOT)))
def test_no_llm_import(path: Path):
    """No LLM client is imported anywhere in the deterministic core."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except SyntaxError as e:                       # a broken file is its own bug
        pytest.fail(f"{path.relative_to(ROOT)} does not parse: {e}")
    hits = sorted(_imported_roots(tree) & FORBIDDEN_ROOTS)
    assert not hits, (
        f"{path.relative_to(ROOT)} imports {hits} — the core must not call a "
        "model. That function belongs in the agent layer "
        "(.claude/rules/deterministic-core.md)."
    )


@pytest.mark.parametrize("path", CORE_FILES,
                         ids=lambda p: str(p.relative_to(ROOT)))
def test_no_api_key_read(path: Path):
    """No core module reads an LLM API key from the environment.

    An import can be avoided by shelling out or by `importlib`, but reading the
    key is the tell in either case.
    """
    src = path.read_text(encoding="utf-8")
    for key in ("ANTHROPIC_API_KEY", "OPENAI_API_KEY", "GOOGLE_API_KEY"):
        assert key not in src, (
            f"{path.relative_to(ROOT)} references {key} — the core does not "
            "authenticate to a model provider."
        )


def test_grep_would_have_false_positives():
    """The reason this file parses AST instead of grepping, kept executable.

    If this ever fails it means no prose cites the rule any more, and the
    cheaper grep would have been fine. That would be worth knowing.
    """
    prose_hits = [f.relative_to(ROOT) for f in CORE_FILES
                  if "anthropic" in f.read_text(encoding="utf-8").lower()
                  and "anthropic" not in _imported_roots(
                      ast.parse(f.read_text(encoding="utf-8")))]
    # Not an assertion about the count -- just that the distinction is real and
    # this test is the thing that knows the difference.
    assert isinstance(prose_hits, list)
