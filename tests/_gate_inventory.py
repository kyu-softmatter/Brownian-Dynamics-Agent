"""Source-derived inventory of every error site in a load-bearing module.

Keyed on the normalised source text of the raising line, not the line number, so
the inventory survives edits above it and only changes when an error site is
added, removed or reworded.

`tests/test_gates_fire.py` fails when a site appears that is in neither list. The
point is not that every site must have a case -- 26 of them do not -- it is that
a NEW gate cannot arrive unobserved without somebody deciding which list it goes
in. That decision is the thing that was missing.
"""
from __future__ import annotations
import ast, json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INVENTORY = Path(__file__).resolve().parent / "_gate_inventory.json"

LOAD_BEARING = ("params", "runcard", "runid", "physical", "metrics", "checks",
                "health", "intake", "goal", "nondim", "scales", "design")


def _norm(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def error_sites() -> dict[str, str]:
    """{f"{module}:{normalised source}": module} for every error site."""
    out = {}
    for mod in LOAD_BEARING:
        f = ROOT / "bdbot" / f"{mod}.py"
        src = f.read_text()
        for n in ast.walk(ast.parse(src)):
            hit = isinstance(n, ast.Raise) or (
                isinstance(n, ast.Call) and getattr(n.func, "attr", None) == "err") or (
                isinstance(n, ast.Call) and getattr(n.func, "id", None) == "Issue"
                and n.args and isinstance(n.args[0], ast.Constant)
                and n.args[0].value == "error")
            if hit:
                # ast.unparse and not the source line: two branches of the same
                # guard can both read `raise RuntimeError(` and would collide.
                out[f"{mod}:{_norm(ast.unparse(n))[:160]}"] = mod
    return out


def load() -> dict:
    return json.loads(INVENTORY.read_text())


def write(covered: list[str], unobserved: list[str]) -> None:
    INVENTORY.write_text(json.dumps(
        {"_why": "Every error site in a load-bearing bdbot module, from the AST. "
                 "`covered` has a case in verify/verify_gates_actually_fire.py; "
                 "`known_unobserved` has never been seen to fire and nobody has "
                 "written one yet. A site in neither fails the suite -- that is the "
                 "whole mechanism.",
         "covered": sorted(covered), "known_unobserved": sorted(unobserved)},
        indent=1) + "\n")
