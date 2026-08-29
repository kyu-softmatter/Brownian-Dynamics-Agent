# deterministic-core — the core does not call an LLM

`bdbot/`, `simbot/`, `cases/` and `tools/` carry no LLM dependency. The
dependency direction is one-way: the agent layer may call the core, the core may
never call the agent layer.

```
pytest tests/test_invariants.py     # this is the check
```

**Do not check this with `grep`.** `A4` was originally written as
`grep -rE "anthropic|claude" simbot/ bdkit/`, and that check **matches the prose
that explains the rule** — it fires once more every time a docstring cites this
file (7 false hits, measured 2026-07-28). A string search cannot tell code from
prose. `tests/test_invariants.py` parses the AST and looks at **real imports
only**.

> **Merge note (2026-08-28).** The predecessor design had a thin LLM layer at
> `agent/llm.py` that called the Anthropic API directly, and the rule named
> `bdkit/` as the core. Both are gone: the runtime is Claude Code itself, so
> there is no API-calling layer to isolate, and the core is now `bdbot/` +
> `simbot/`. **The rule got stronger, not weaker** — there is no sanctioned
> place for an LLM call in this repository at all. And it was found unenforced:
> `tests/test_invariants.py` did not exist in any of the three predecessors, so
> for a month `A4` was true but unchecked. That is the unwired-checker failure
> this project keeps documenting, and it is now wired.

**Why (the triggering incident):** an architecture decision, not an accident
(2026-07-27). There is one reason and it is sufficient: **when a result is wrong
you have to be able to tell whether the physics was wrong or the model reading
it was.** In simulation that is fatal — a plausible-but-wrong `g(r)` is
indistinguishable by eye, and with no grader in this domain there is no way to
learn you were wrong.

The separation has paid for itself at least once. When a trap's autocorrelation
`τ` came out 70 % above theory, there were exactly two candidates: the
simulation is wrong, or the estimator is wrong. Because the core is
deterministic, it was settled by **feeding in synthetic OU data whose answer was
known** — it was the estimator. With an LLM mixed into the core there would have
been three candidates, and the third would not even reproduce.

**How to apply:**
- If adding something to the core makes you want `import anthropic`, **the
  location is wrong.** That function belongs in the agent layer
- Anything the agent layer produces enters the core only after **a fixed schema
  plus a deterministic validator**. Never build a path where an LLM-supplied
  value reaches a config unchecked
- Test a new analysis estimator against **synthetic data whose answer is known**
  (`tests/test_s2_estimators.py`, `tests/test_s7_structure.py`). That is not
  physics verification — it asks "does the code answer correctly on an input
  whose answer it should know"
- `tests/test_invariants.py` enforces this rule. Do not disable it

**Anti-patterns explicitly forbidden:**
- **The convenience import** — letting an LLM into the core because "it is just
  one call here"
- **Bypassing the validator** — applying an LLM suggestion without the validator
  because "this one is obviously right"
- **Creating a third candidate** — adding an unreproducible cause to the
  debugging path alongside physics and estimator

See also: [axioms](axioms.md) ·
[docs/02-verification.md](../../docs/02-verification.md)
