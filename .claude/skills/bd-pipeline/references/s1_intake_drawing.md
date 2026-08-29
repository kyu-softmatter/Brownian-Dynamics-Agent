# S1 — the hand-drawing reading protocol

> **This document is the only content unique to the skill layer.** Every other
> stage is a core function call, but "what system is this drawing?" cannot be
> expressed as code.
>
> **If this is wrong, everything downstream is wrong.** It is the most expensive
> place to make a mistake.

## 0. First — what do you trust in this drawing?

> **Do not trust absolute sizes in a hand drawing.** People draw the box small
> and the particles large — otherwise the particles become invisible dots. Taking
> the drawn particle-to-box ratio as `φ` almost always overestimates it.

| Trust | Do not trust |
|---|---|
| **topology** — what is inside / beside / on top of what | absolute sizes |
| **ratios** — particle : box ≈ 1:20 | absolute arrow thickness or length |
| **counts** — three circles may mean three species | drawn particle count = simulation `N` |
| **symmetry** — is it left-right symmetric, is there an axis | precise angle values |
| **written numbers and units** — `R = 5 μm`, `T = 300 K` | any length read off by eye |
| **written equations** — `r = √(x²+y²)` states the dimensionality | |

**Mistaking the drawn particle count for `N` is the most common failure.** The
drawing is a sketch; `N` is **the number statistics requires**. In the first
example the drawing had one circle and `N = 1000` was used (non-interacting, so
one snapshot is 1000 independent samples).

## 1. Inventory — what is actually drawn

Sweep in this order, starting with what is easiest to miss:

1. **Text and numbers** — read these first. They are the most reliable and they
   constrain every other interpretation. **Always transcribe the unit.** Writing
   `10 pN/μm` down as `10` destroys the information at that instant.
2. **Written equations** — no `z` in `r = √(x²+y²)` is a 2D signal. **The
   equation states the dimensionality.**
3. **Particles** — count, size differences, species markers (colour, hatching,
   labels), any specially marked individual (a probe)
4. **Boundaries** — solid line (wall) vs dashed (periodic) vs absent. Slit,
   cylinder, sphere geometry. **The dimensionality**
5. **Arrows** — position, direction, length. **List the candidates for what it
   is: force / velocity / flow / time progression.** Do not settle on one
6. **Axes and labels** — two axes is a 2D signal, three is 3D
7. **Hand-drawn graphs** — an expected curve is a basis to compare against the S2
   prediction. Read the axis labels
8. **Captions and questions** — if the user wrote down what they want to know,
   that *is* the `question`

## 2. The three-way split — S1's central deliverable

| Grade | Definition | Example |
|---|---|---|
| `observation` | **read directly** off the drawing | "one circle", "an `R = 5 μm` leader line points at the circle" |
| `inference` | **derived** from the drawing plus physics knowledge | "no `z` in the equation and only x, y axes → 2D" |
| `assumption` | **absent from the drawing, so I supplied it** | "the medium is water, `η = 0.856 mPa·s` @ 300 K" |

Attach **`confidence: high/medium/low` plus a one-line basis** to every item.

> ★ **Only Opus may fill `inference` and `assumption`**
> (see [`.claude/README.md`](../../../README.md#authority-boundary)). Cheaper
> models produce `observation` and `derived` only.

Why the split matters: **the sensitivity analysis takes its target list from the
`assumption` entries.** Record an observation as an assumption and you test an
uncertainty that does not exist; record an assumption as an observation and you
miss one that does.

## 3. Ambiguity — state the candidates and build a discriminator

> Do not arbitrarily pick one. **Write down 2–3 candidates and predict how the
> result differs between them.** Then the user can decide immediately, and if they
> cannot, **run both** and hand back the criterion that decides it.

Format:

```markdown
### A1 — is this 2D or 3D?
- **A1-a (2D, adopted)** basis: `r = √(x²+y²)` has no `z`, two axes. confidence medium
- **A1-b (3D cross-section)** basis: a real optical trap confines in 3D. The drawing may be a section
- **Discriminator:** `⟨r²⟩(3D)/⟨r²⟩(2D) = 3/2`. Per-component `⟨x²⟩` **cannot** tell them apart
- **Cost:** running both branches doubles the runs. In this trap system, 2.3 s → 4.6 s (harmless)
```

The first example did exactly this, and running both dimensions gave `1.4955`,
**confirming `3/2` to 0.3 %** — an independent verification that running one
dimension could not have produced.

**An ambiguity with no discriminator goes to the user** (within the question
budget).

## 4. Gaps — required information that is absent

List what the simulation needs but the drawing lacks, and set a policy for each:

| Policy | When | Result |
|---|---|---|
| `ask_user` | it flips the conclusion, or not even the order of magnitude is known | spends question budget (three per round) |
| `fill_from_knowledge` | there is a basis in `knowledge/wiki/concepts/` | `provenance: from_knowledge` |
| `assume_and_flag` | the order of magnitude is known and sensitivity can check it | `provenance: assumed` + `confidence` |
| `sweep` | the value changes the regime and the candidates are few | thin runs at several values |

**Do not leave a blank as "I don't know."** Propose, proceed, and let sensitivity
check it. If an assumption turns out to be irrelevant, **report explicitly that
"this value does not change the conclusion"** — that is also a result.

⚠️ But `assume_and_flag` has a failure mode this project actually hit: an
`assumed` value written into the spec with a confident-looking `tier` becomes
indistinguishable from a measured one downstream. `T = 300 K` is recorded as
tier 1 across every case and is in fact a **choice** inherited from a sketch with
no temperature — worth −4 % to −14 % on every timescale. Mark the tier honestly.

## 5. The `question` — is it falsifiable?

| ✅ | ❌ |
|---|---|
| "How much slower does diffusion get?" | "What happens?" |
| "Quantify this trap's thermal fluctuation and compare against the analytic solution" | "Check whether it works" |
| "Does MIPS occur at `Pe = 45`?" | "Active matter simulation" |

**If the answer does not come out as a number or a true/false, rewrite it.**

## 6. Failure cases — suspect in this order

Recorded in `knowledge/wiki/findings/`, plus the common ones from the literature:

| Symptom | Cause |
|---|---|
| the result differs from the literature by orders of magnitude | an arrow was read as a force but was a velocity field |
| a 2D result disagrees with 3D literature | a 2D drawing was read as a 2D system but was a 3D section |
| the statistics are far too poor | the drawn particle count was used as `N` |
| the timescale is wrong by exactly 2× | `R` was read as a diameter when it was a radius (or the reverse) |
| `γ` is wrong by 2× | a diameter went into `6πηa` — **the most common mistake** |

★ The last two are **only revealed by measuring `τ_trap`.** `⟨x²⟩ = kT/k` is
analytically independent of `a` and `η`, so a static measurement cannot
discriminate. The first example had exactly this structure, and with no
experimental `f_c` available it **could not be closed by measurement and was
closed by parsimony** — and that fact was recorded.

## 7. Deliverable

The intake document — sections 2 through 5 above, verbatim. Format example:
[`01_intake.md`](../../../../runs_s1s8/2026-07-28_trap-2d-5um_2dfb9d/01_intake.md)

Put the source material in `intake/<case>/` and **record its sha256**. The
provenance chain is verified by hash.

## 8. Next — what you show the user

When the reading is done, show a **proposal table**. The form is not *asking* but
*being confirmed*:

```markdown
| Value | Basis | Confidence | Where to change it |
|---|---|---|---|
| `η = 0.856 mPa·s` | IAPWS interpolation at 300 K, concepts/water-298k.md | medium | `spec.yaml` `medium.eta_si` |
| `d = 2` | no z in the equation (ambiguity A1) | **medium** | `spec.yaml` `geometry.dim` |
| `N = 1000` | non-interacting → one snapshot = 1000 independent samples | high | `spec.yaml` `species[0].n_simulated` |
```

**Tell them where to change it.** Do not make a person hunt for the place to
edit.
