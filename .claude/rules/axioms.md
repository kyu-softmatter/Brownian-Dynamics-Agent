# axioms — the four that cannot change without human approval

This is the L0 axis. Other rules can be retired when circumstances change; these
four need explicit human approval to change. All four are devices against
**silently wrong results**, and in a domain with a weak verification oracle that
is the most expensive failure there is.

**A1 · Every numerical claim is verified by three pieces of evidence of
different kinds.**
If two of them disagree, the number is fake. Three pieces of the *same* kind
count as one — three runs differing only in seed are one piece of evidence.

**A2 · Do not emit a number without an error bar.**
State block averaging alongside `τ_ac`. A bare value says nothing about how
confident you are, so it can neither be compared nor refuted.

**A3 · The cost gate must be cleared before running.**
If the estimate exceeds the budget, do not run. In a chat session, **ask**
rather than refuse. Either way, do not start without an estimate.

**A4 · The deterministic core does not call an LLM.**
`bdbot/`, `simbot/`, `cases/` and `tools/` carry no LLM dependency.
`pytest tests/test_invariants.py` must pass — **not `grep`.** A string search
matches the prose that *cites* this rule, so the test parses the AST and looks
at real imports only. Details in [deterministic-core](deterministic-core.md).

**Why (the triggering incident):** these were design decisions, not accidents
(2026-07-27). So unlike every other rule here they carry no date and no cost to
cite — **and that fact is not hidden.** Of the four, `A3` came from a predecessor
project where a multi-day job was started casually, and `A2` came from realizing,
while trying to compare a measured `D` against a published one, that the
comparison did not hold without an error bar. Neither left a path or a date
behind.

**How to apply:**
- Before putting a number in a conclusion, check that you are holding its error
  too (`A2`)
- Check that you counted three **different** kinds — literature, self-consistency,
  analytic limit (`A1`)
- Before a heavy run, look at the cost estimate first (`A3`)
- When adding code to the core, check that no LLM dependency crept in (`A4`)

**Anti-patterns explicitly forbidden:**
- **Inflating the evidence count** — counting runs that differ only in seed or
  box size as different kinds of evidence
- **Bolting an error bar on afterwards** — publishing the value first and
  computing the error later
- **Clearing the gate by raising the budget** — when the estimate exceeds the
  budget, raising the budget without looking at the reason. Raising it is
  allowed, but first say **why this system costs that much**

See also: [deterministic-core](deterministic-core.md) ·
[verify-against-literature](verify-against-literature.md)
