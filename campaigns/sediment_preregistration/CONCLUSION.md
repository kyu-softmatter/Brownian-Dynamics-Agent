# S8 — sediment-pmma-3d: is the Barker–Henderson mapping required?

**Status: proposed. `confirmed_by` is human-only** (CLAUDE.md), so nothing below
is a confirmed result.

Campaign `sedimentation-eos-S31`, 8 runs, 1.671 core-hours, sealed against
`prediction.yaml` **revision 4** (the source has since moved to revision 5 — see
its banner; the two differ in exactly two vestigial `role:` fields).

---

## 1 · The question, and the answer

> Does free-draining BD with a soft (WCA) core reproduce the hard-sphere equation
> of state — and does it do so only after the excluded volume is mapped to a
> Barker–Henderson effective diameter, or at the nominal diameter as well?

**The mapping is required.** Every discriminator says so, and they disagree with
the unmapped alternative by 4.5 to 18 σ.

| discriminator | measured (6 permitted runs) | vs CS(φ_eff) | vs CS(φ_nominal) |
|---|---|---|---|
| `Z_dev_wmean_vs_CS_eff` | **+1.55 ± 0.89 %** | 1.8 σ | **5.6 σ** |
| `l_g_fitted_dense` | **16.789 ± 0.169 d** | 1.8 σ | **4.5 σ** |
| `phi_wall_measured` | **0.11629 ± 0.00074** | **6.2 σ** | **18.4 σ** |

`Z` agrees with Carnahan–Starling at the effective diameter **inside the
pre-registered 2 % band**, and the same numbers over all 8 runs are +1.14 ±
0.70 % (1.6 σ / 7.6 σ), so the exclusion of two runs is not doing the work.

★ **The strongest single run is the one that started furthest away.** The
convergence arm began at the ideal-gas profile (φ_wall = 0.1589), travelled 89 %
of the way to CS(φ_eff) in 8 τ_sed — slightly past it — and returned
`Z_dev = +0.099 %`, the closest of all eight. A run that was never near either
candidate arrived at one of them.

★ **F7's bottom panel is the result, drawn.** The measured `arm B / arm A`
profile ratio is **flat at 1.0** from h = 2 d to ~45 d, while the ratio of the two
*candidate* profiles they started from runs from 0.93 at the wall to 1.15 at
h = 60 d. Both arms left their own starting profile and arrived at a common one.
Nothing in that panel depends on a threshold.

⚠ F7's **top** panel cannot show the discrimination and should not be read as if
it could: on a log axis spanning four decades a 9 % difference in φ(0) is
invisible. That is why the bottom panel exists.

★ **F7's middle panel explains the χ²/ν of 2.6.** Both arms scatter about 1 with
a slight *upward* trend confined to large h — the tail is still filling. The
residual drift is therefore in the dilute tail, and `Z_dev_wmean` is
**count-weighted**, so the tail contributes almost nothing to it. That is why the
convergence arm can be flagged `DRIFTING` (χ²/ν = 4.23) and still return the
closest `Z_dev` of all eight runs: the part of the profile that is still moving is
not the part the equation of state is read from. This was not foreseen, and it is
the reason §2's two verdicts can differ without either being wrong.

**Finite size: no effect.** 12 d against 24 d (4× the particles) agree to
0.7–1.6 σ on all three discriminators.

---

## 2 · ⚠ The pre-registered verdict is INCONCLUSIVE, and that stands

`decision_rule` step 5 made the verdict rest on **stationarity**: arm A
stationary and arm B drifting, or the mirror. Measured:

| | seeds with χ²/ν ≤ 2 | χ²/ν |
|---|---|---|
| arm A, CS(φ_eff) start | **1/3** | 2.61 ± 1.35 |
| arm B, CS(φ_nominal) start | **1/3** | 2.05 ± 0.68 |

The two arms differ by **0.6 σ**. Neither is stationary, so the sealed rule
returns `INCONCLUSIVE (neither stationary)`, and this document does not override
it. §1 is what the evidence supports; §2 is what was registered. They differ, and
the difference is the finding about the *method*:

> **The decision rule made stationarity the only route to a verdict, and treated
> "neither arm stationary" as a failure mode — when convergence from four
> different starting conditions to one state is stronger evidence than
> stationarity, and was available in the same data.**

Stationarity can only ever say *"we started here and nothing moved us"*, which is
one step from the circularity the two-arm design existed to remove. Convergence
says *"four different starting points arrive here"*. The rule preferred the
weaker of the two, and it was written before either had been measured.

⚠ The runs are also genuinely still relaxing, mildly: the mean half-window shift
is +1.5 % and `tools/postmortem.py` independently flags 2 of 8 as
`EQ_INSUFFICIENT`. 4 τ_sed is not full equilibration. That is a real limitation
and it is why §1's error bars are the seed-to-seed spread and not the within-run
block SEM.

---

## 3 · A residual the mapping does not explain

`phi_wall` sits **6.2 σ (+4.0 %) above** CS(φ_eff) — the only discriminator that
misses its own prediction significantly.

It is almost certainly **not** an EOS effect, and the reason is that it does not
propagate: `Z_dev_wmean` is measured over `h ≥ 2 d`, above the wall's reach
(`U_wall(2 d) = 4.5e-9 kT`), and reads only +1.55 %. The residual is confined to
the 0.5 d bin adjacent to the wall, which is exactly where the simulation has a
**structural layer** that the CS hydrostatic profile has no notion of — the
profile ODE is a local-density-functional statement and knows nothing about
layering against a hard boundary.

So the honest reading: the Barker–Henderson mapping is necessary and sufficient
for the EOS in the bulk of the column, and the wall-adjacent bin is outside what
either candidate profile models. **Testing that would need `g(z)` resolved below
0.5 d, which this campaign did not measure** — named here rather than asserted.

---

## 4 · What this campaign cannot claim

- **Not** that the paper's *"excellent agreement"* is reproduced. Newman &
  Yethiraj give no percentage, so there is nothing to reproduce *to*; the 2 %
  band is ours, inherited from [[wca-reproduces-carnahan-starling]].
- **Not** `K_L = −2.80`. Free-draining BD has no hydrodynamic interactions and
  cannot produce it — declared out of scope before the run (rule 7′).
- **Not** a general statement about hard spheres. One cell, one φ, one
  `eps_wca`, WCA rather than hard, monodisperse, no hydrodynamics. The finding
  transfers as evidence, not as a theorem.
- **Not** full equilibration. See §2.
- **Not** a finite-size exponent. Two cross-sections cannot fit one.

---

## 5 · The next experiment

Not more seeds in either arm — the discrimination is already 4.5–18 σ. In order:

1. **A second convergence seed, longer.** The ideal-gas arm is the strongest
   evidence and there is exactly one of it. 16 τ_sed × 2 seeds ≈ 1.2 h.
2. **`g(z)` at 0.1 d resolution in the first 2 d**, to settle §3's residual.
   Cheap — it is a re-analysis if the trajectory is stored, which it is not
   today (`Build.gsd_path` is declared in `bdbot/run.py` and never read).
3. **A dt/2 replicate.** Named as a deliberate omission in the sealed cost block.

---

## 6 · What the campaign cost in defects, and what caught each

Five revisions of the pre-registration, three stopped launches. Every defect was
found by execution or by reading one artefact against another — **none by
inspection of the thing itself**.

| # | defect | found by |
|---|---|---|
| 1 | `l_g_fitted` as one `implementation_check` against 13.375; the CS ODE gives 15.44 d there | integrating the ODE |
| 2 | two of three gates had a tolerance inside their own 1 σ — one 82 % likely per seed to declare IMPLEMENTATION FAILURE on working code | the design-power Monte Carlo revisions 1–2 never ran |
| 3 | the sealed plan required a half-window comparison the output could not produce | reading the plan against the code |
| 4 | those halves split 6 frames against 30; the statistic returned an ordinary 1.6356 | reading `frames_h1` out of a test run |
| 5 | `Z_dilute_tail` and the half-window χ² gated at equilibrium-only identities | a production-length pilot |
| 6 | the analyzer kept its own hardcoded gate list | running it after the demotion |
| 7 | two `role:` fields in the sealed document contradicted its own gate table | a new role-agreement test |
| 8 | `tools/postmortem.py` judged every predicted observable at a flat 5 % **ignoring `role` and the declared tolerance** — recording 4 of 8 runs as failed when the only real gate passed on all 8 | running the mandated post-mortem |

★ Defect 8 is the one worth keeping: it contradicted CLAUDE.md's own sentence —
*"Only `implementation_check` mismatches are FAIL; `hypothesis` mismatches are
reported as results"* — inside the tool the working practice makes mandatory, and
it wrote its verdict into `record.json`, which is the permanent record every
later query reads. It had been that way for every case before this one.

In hindsight, had revisions 1–3 been run as sealed:

| gate as sealed in r1–r3 | outcome on the 8 real runs |
|---|---|
| `l_g_fitted_tail` within 5 % | **8/8 FAIL** (range 11.6–19.4 d) |
| `Z_dilute_tail` within 5 % | 7/8 FAIL · at 10 %: 4/8 FAIL |
| `D_xy` within 10 % | 5/8 FAIL · at 20 %: **0/8 FAIL** |

The campaign would have reported IMPLEMENTATION FAILURE and produced no physics.

## See also

[`prediction.yaml`](prediction.yaml) · [`analysis_plan.yaml`](analysis_plan.yaml) ·
[`verify/verify_sediment_design_power.py`](../../verify/verify_sediment_design_power.py) ·
[`knowledge/wiki/systems/passive-sphere--sedimentation.md`](../../knowledge/wiki/systems/passive-sphere--sedimentation.md)
