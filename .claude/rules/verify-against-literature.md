# verify-against-literature — literature is cited, not remembered

When comparing a simulation result against a published value, **open the
distillation in `knowledge/source/papers/` and cite it.** Do not write from
memory. If you need a number the distillation does not carry, distil it first.

The comparison has four layers and they must be **different in kind** (`A1`).
Three of the same kind count as one.

| | Evidence layer | Example |
|---|---|---|
| ① | analytic solution · limit | `MSD = 2·dim·D·t`, `D = kT/γ`, trap `⟨r²⟩ = dim·kT/k` |
| ② | self-consistency | input `D₀` ↔ measured `D_msd`, input `k` ↔ recovered `k` |
| ③ | literature benchmark | Carnahan–Starling `Z(φ)`, `φ_freeze = 0.494`, Zahn's melting `Γ` |
| ④ | convergence | does it stay put when `Δt`, `L`, `N` and the initial condition are shaken |

**If two pieces of evidence disagree, stop. Do not average them.**

**Why (the triggering incident):** this is design, not an accident (2026-07-27) —
and that is not hidden. But **why it is needed was experienced directly.** Three
"method numbers dressed as physics discrepancies" were caught in a row: a trap
`τ` 70 % high (truncation-choice bias), an MSD plateau 3 % high (a single sample
with only one or two time origins), and a `D` 9 % low (a structural downward bias
from the frame interval). In all three the simulation was fine and **the analysis
was wrong.**

Had only the literature been consulted, all three would have ended as "our system
differs from the literature." What separated them is that **the analytic solution
(①) and self-consistency (②) are independent of the literature (③).**

**How to apply:**
- When a published value goes into a conclusion, open
  `knowledge/source/papers/<slug>.md`. If it is not there, distil first
- When there is a disagreement, **suspect the analysis code first.** Re-test the
  estimator against synthetic data whose answer is known
- Record the comparison in `knowledge/wiki/findings/` — **agreement is a finding
  too.** The next run cites it
- Do not call a system "verified" when it is in no benchmark. No threshold, no
  verdict
- ⚠️ 38 of the 42 distillations here are the group's own published work, so the
  literature layer (③) is narrower than it looks. Weight it accordingly

**Anti-patterns explicitly forbidden:**
- **Citing from memory** — writing "Carnahan–Starling is roughly…". Open the
  distillation and take the number
- **Averaging a disagreement** — splitting the difference, or quietly dropping
  one side, when two pieces of evidence conflict
- **Filling the three with one kind** — counting three seed-only runs as three
  pieces of evidence
- **Citing a parameter marked `reproduced: no`** — that is a record of a fact, not
  a basis for one

See also: [axioms](axioms.md) · `knowledge/wiki/CLAUDE.md` ·
[docs/02-verification.md](../../docs/02-verification.md)
