# bdbot — the Brownian dynamics simulation pipeline

A system that reads a physical system out of a hand-drawn sketch, a note or the literature,
non-dimensionalizes it, runs it in HOOMD-blue, and feeds the successes and failures back as knowledge.
The full design is in [mater_plan.md](../../docs/history/2026-08_simulation_auto_master_plan.ko.md).

## Current status (2026-08-06)

| Phase | Status |
|---|---|
| 0 · a general environment + the HOOMD capability survey + 15 APIs verified | ✅ done |
| 8-min · knowledge capture (this file + 2 skills) | ✅ done |
| **1-A · `trap-2d-5um` end to end** | ✅ **done — 4 observables match their analytic solutions** |
| **1-B · `soft-r3-2d-A-sweep` end to end** | ✅ **done — 7 runs, 5 verifications + 2 convergence checks** |
| **1-C · abstraction (the `bdbot/` package)** | ✅ **done — metrics identical over 96 and 124 fields** |
| **the front end · the intake schema + checks + CLI + skill `bd-intake`** | ✅ **done — checks 27/27** |
| **L3 · `NondimSpec` (the L2↔L4 contract) + `specs/`** | ✅ **done — 3 cases migrated, adversarial checks 33/33** |
| **L3 · non-dimensionalizing the remaining 2 cases** | ✅ **done — 8 specs, adversarial checks 16/16** |
| **L4 · the numerical-health verdict** | ✅ **done — `bdbot/health.py`, adversarial checks 45/45, runs 82/82 HEALTHY** (but the step-resolution check applies from new runs onwards — see below) |
| **L5 · the assembly contract (`@RUN.builder`/`Build`) + the L6 run loop + L7 storage** | ✅ **done — `bdbot/run.py`, all 7 cases adopted** (the first 5 compared by re-running after migration) |
| **1-D · `chain-bend-2d-dlvo` (the alternative-hypothesis branch)** | ✅ **done — 145 runs, 3 protocols synthesized. The conclusion is protocol-independent** |
| **1-E · `network` (a 3D colloidal network — the 2-stage structure of rule 8)** | 🔶 **stage 1 (gelation) in progress — specs 3, 1 completed run** (the `sprout` topology, N=512). 2 compression-gelation runs incomplete · stage 2 (driving, G′(ω)) not started |
| **1-F · `chain-relax-2d-dlvo` (the static/undriven counterpart of chain-bend-2d-dlvo)** | 🔶 **the pilot is done (n=9, 1 seed) — awaiting ensemble extension** (see §1-F below) |
| KB · `record.json` + `kb/entries/` (SQLite later) | ✅ running — **353 entries** (227 with a run · 126 knowledge without one) |
| **literature · 2 books distilled (`docs/books/`)** | ✅ **done — entries 24, execution-verified 56/56** (see §literature below) |

**L5 is not a new layer but a correction of what was already there** (2026-08-05). The
`@RUN.builder("case name")` registry + the `Build` dataclass in `bdbot/run.py` are exactly the L5
design of mater_plan.md (`build(spec) -> ...`), and `execute()` doubles as L6 (the run loop) + L7 (storage).
`chain-bend-2d-oscill`, `trap-drag-2d-hex300` and **`trap-2d-5um`** (migrated and compared by re-running 2026-08-05)
already use this contract. `trap-2d-5um` had its own `HarmonicTrap` class replaced with `bdbot.traps.make_trap`
and was re-run, and 5 observables were compared against the old runs (`a49f2508556b`, `d724e8d507cc`) down to
display precision — they match completely (the same seed 20260803, the same physics, so the trajectory reproduces exactly).
Two bugs were caught during the migration (2 KB `tooling` entries): ① `MET.build(extra=...)` merges with `m.update(extra)`,
so `finalize()`'s `result` arrives not at `metrics["extra"]["result"]` but at the **top level**,
`metrics["result"]` — otherwise it is silently nan. ② if `sample()` returns the whole N-particle array
as it stands for every sample, `observables.npz` swells 448KB→148MB (330×) —
derived quantities have to be accumulated in a closure and the raw arrays filled only for the subset needed.

**`abp-rod-2d-run-flip` was also migrated and compared by re-running** (2026-08-05). This system has neither pair
interactions nor a trap (the active force is non-conservative), so there is no potential energy for
`pe_per_particle` to point at — instead it uses `⟨cos Δθ⟩` (the orientational correlation against the previous call),
which the original implementation already used as its equilibrium indicator: when normal, every call is a different
noise realization so it is never exactly constant, and if the rotational diffusion dies it locks to exactly 1.0 and
FROZEN is caught. The re-run result matches the old run's (N=1000, `74fb2d81066a`) 5 observables down to display
precision. All three cases had the same kind of hole during the migration —
**if the values build() needs (seed, L_star, Dr_star and so on) are missing from spec.params/numerics,
build(spec) cannot re-read the case YAML and cannot obtain them** (a KB
`tooling` entry). All three were added and the run_id changed (they are physics hash fields), but the existing runs
remain as they were.

**`soft-r3-2d-A-sweep` was also migrated and compared by re-running** (2026-08-05). A full-size re-run (A=100, N=400,
3.3 million steps) produced observables in `runs/soft-r3-2d-A-sweep__A100__30caa5c9e0` that are **bit-for-bit
identical** to the old run's (energy consistency 105.510722899358, the hexagonal NN distance
1.6170154574513436, ψ₆ 0.8854118251800537 — matching to 15 decimal places). This system has no analytic
solution, so it has post-hoc guards (`post_checks` — r_tab_min/min_sep margins and the like), but the
`checks` schema of `bdbot.run.execute()` still has only the two stages "design"/"post_run", so the post-hoc guards
`finalize()` produces went into `metrics["result"]["post_checks"]` (not top-level `checks[]`) — if a third
case needs this pattern, it gets promoted properly into `bdbot.run` then.

**The L5 migration finished all 5 cases that existed at the time** (2026-08-05). All three cases
(`trap-2d-5um`, `soft-r3-2d-A-sweep`, `abp-rod-2d-run-flip`) had re-run results matching the old runs down to
display precision (usually 15 decimal places), confirming that the refactor did not change the physics. The
`chain-bend-2d-dlvo` and `network` built afterwards were built on this contract from the start, so **all 7 cases
now** use `@RUN.builder` (the migration comparison applies only to the first 5). The details are in the
`run` entry of [`bdbot/__init__.py`](../../bdbot/__init__.py).

**When starting a new case** — the tools enforce the procedure:
```bash
$PY -m bdbot.cli status                      # where it is stuck
$PY -m bdbot.cli intake init  intake/<case>  # template → read the image and fill it in (skill bd-intake)
$PY -m bdbot.cli intake check intake/<case>  # FAIL / BLOCKED / READY
$PY -m bdbot.cli system check intake/<case>  # tier · derived_from · recompute the derived values
$PY -m bdbot.cli nondim spec  <case>         # L3 → specs/<run_id>.json  (does not run)
$PY -m bdbot.cli nondim show  <run_id>       # reproduce from the spec alone + verify the hash ← human check #3
$PY -m bdbot.cli health --gate specs/<run_id>.json   # L4 pre-run gate
$PY -m bdbot.cli health runs/<run_id>                # L4 numerical-health verdict
```

**The L3 output is one file, `specs/<run_id>.json`** ([`bdbot/nondim.py`](../../bdbot/nondim.py)).
It is the **only contract** between L2 and L4, and L4 reads only it, importing no case code.
Self-sufficiency is judged by `nondim show` — if the whole report can be drawn from the spec alone, it is self-sufficient.
The spec checks four things about itself (**a different layer** from the physics checks, `checks`):
the completeness of the ledger (4 required roles) · whether a dimensionless number really is the ratio of two ledger
entries · whether the reference is in the ledger · whether the inverse-transform anchor holds. The adversarial
checks are `scratch/verify_nondim_guards.py` (33/33).
The shared code is [`bdbot/`](../../bdbot/__init__.py) (20 modules). A case script writes **only the physics unique to it** —
take [`cases/soft_r3_2d.py`](../../cases/soft_r3_2d.py) as the model.
**`BLOCKED` is not a failure** — it means what is missing has been narrowed to one thing, and it must not be made
READY by inventing it (rule 3).

**What goes into `bdbot/` is judged solely by "has it appeared twice".** If it comes up again in a third case,
it goes in then. What is deliberately left out for now: the equilibrium criterion · the observables · the verification
strategy · the choice of the governing timescale · the initial placement · the sampling loop (skill `bd-physics` §6.3).

After a refactor, always confirm the results are unchanged with `$PY scratch/verify_1c_equivalence.py`.

**`run_id` hashes the physics fields only** ([`bdbot/runid.py`](../../bdbot/runid.py) `DOC_KEYS`). Fixing a comment,
a source or a derived value must not invalidate a run — there was a time when merely adding a `derived_from` field
changed 1-A's run_id. **The opposite direction is equally dangerous**: the 1-B spec had no physical system at all,
so changing `d` 5µm→0.5µm and `η` by 62× (a 16.1× difference in τ_B) left run_id **the same**, it was mistaken for
a completed run of a different system, and an old result was reported as the new system's. The hash has to cover
**everything that fixes the physics** and exclude documentation and derived values — both are needed (`scratch/verify_l3_spec_gaps.py`).

**All 8 cases are L0, L2 and L3 READY** (`bdbot.cli status`, measured 2026-08-06 — specs totalling
**278**, 6 cases with runs):

| Case | Script | Specs | Runs |
|---|---|---:|---:|
| `abp-rod-2d-run-flip` | end to end (O) | 3 | 4 |
| `chain-bend-2d-dlvo` | end to end (O) | 162 | 145 |
| `chain-bend-2d-oscill` | L3 | 15 | 0 |
| **`chain-relax-2d-dlvo`** | end to end (O) | 6 | 4 |
| **`network`** | **— (unregistered, below)** | 3 | 1 |
| `soft-r3-2d-A-sweep` | end to end (O) | 2 | 9 |
| `trap-2d-5um` | end to end (O) | 3 | 4 |
| `trap-drag-2d-hex300` | L3 | 84 | 0 |

Four have an end-to-end script (through L4 execution) (`trap-2d-5um` · `soft-r3-2d-A-sweep` ·
`abp-rod-2d-run-flip` · `chain-bend-2d-dlvo`), and two go **only as far as non-dimensionalization**
(shown as `L3` in `status`): `trap-drag-2d-hex300` (a commensurate hexagon, 17×18=306) ·
`chain-bend-2d-oscill` (an ω sweep, 7 points). ⚠️ **`network`'s script cell being `—` is not because there is no
script** — [`cases/network_3d.py`](../../cases/network_3d.py) exists and produced runs, but it is not registered in
the `CASE_SCRIPTS` registry of [`bdbot/cli.py`](../../bdbot/cli.py).
So `bdbot.cli run network` refuses with "there is no end-to-end script" (register it and it becomes `O`).

⚠️ **`status`'s run count counts only runs that have a `result.txt`** — `network` has 4 run directories
(including `__FIGS`) but only 1 is counted as complete. The 2 compression-gelation runs
(`st0.001-ag0.2` · `st0.01-ag0.2`) were stopped leaving only `traj_A.gsd` and are incomplete.
This convention is the same one as tooling incident ① (§1-D below).

`abp-rod`'s anisotropic translational friction was concluded to be **impossible in BD** (no HI — the isotropic
average γ̄ is used, [mater_plan §20 question 10](../../docs/history/2026-08_simulation_auto_master_plan.ko.md)).

**`trap-drag` L4 results — corrected by re-running 7 velocities × 9 seeds = 63 runs (2026-08-06)** ⭐️
`runs/trap-drag-2d-hex300__ENSEMBLE/vsweep_ensemble.{png,json}` ·
`scratch/trap_drag_vsweep.py`. **Every error is the scatter between realizations (the ensemble SEM)** —
the old values used a single run's block SEM per velocity.

| What | The old value (single run, block SEM) | **The new value (9 seeds, ensemble SEM)** |
|---|---|---|
| yield force `F(v→0)` | 35~45 kT/d, 11~20σ from 0 | **34~41 kT/d, 9.5~12.6σ** — the magnitude holds, the significance comes down |
| shear thinning `F/γv` | 16.9 → 1.87 | **16.8 → 1.84** (9.2×) — reproduced |
| **defect count** | **independent of v** (6.8~9.4) | ⛔ **there is a v dependence** (max−min = 5.2·SEM) |
| **recovery rate** | −0.5 → 67%, **monotonically increasing** | ⛔ **non-monotonic**, −13% → a maximum of 36% (γv=72.6) then decreasing |

⛔ **Two things were overturned**: ① the defect count is not independent of v — it **peaks at 9.43** around γv≈24
and falls to 6.8~7.0 at both ends (2.42, 1549), **a non-monotonic bump**. "There is no critical velocity"
holds, but "independent of v" does not. ② the recovery rate 67% is one realization's value, and the ensemble mean
has a maximum of 36%. **Both were misread because of a single run + an underestimated error.**

⚠️ **The factor by which the block SEM underestimates differs per observable and per velocity.** The widely cited
**2.35×** is the `ΔU/particle` value, and in the same run `F_drag` was **1.37×** — they must not be mixed.
Across the whole v-sweep `F_drag` is **1.09~2.28×** (median 1.63) and **grows with the velocity**
(1.1 at low speed → 2.28 at high). In a system that creates stochastic defects, citing an error without a
seed ensemble is overconfidence.

✅ **The observable definitions did not change in the L5 migration** (2026-08-06, an 8-seed comparison of the
v=0.5 group — `scratch/verify_trap_drag_generations.py`). The old and new runs' `result` scalars **55 + npz arrays
19 are all bit-identical**(relative difference <1e-12) and the step count is the same too (2,330,846).
The reason the spec hashes differed is only that `system.external.drag_velocity` changed from the single value 0.5 → a 7-point list — `params` and `numerics` have **0 differences**. That is, the correction above is not a code
change but **a single run + an underestimated error**. Only the wall clock got 1.65× faster (all 8 runs 1.64~1.68).
⚠️ For the other velocities the legacy was overwritten, so this comparison is **impossible** — the apparent 6×
speedup at v=1.5 and the like is unverified (there is no knowing under what load the legacy wall clock was measured).

**The yield force is an extrapolation, so it depends on the functional form** — rather than picking one, three are given:
ⓐ the lowest-speed measurement `40.7 ± 3.2` (12.6σ, no extrapolation — the most defensible) ⓑ log-linear
**diverges** at `v→0` and cannot define a yield force ⓒ Herschel–Bulkley `F_y = 34.2 ± 3.6` (9.5σ, n=0.90,
χ²/dof=2.97).

### The case where L4's core check had never once run (2026-08-05) ⭐️

`health.py` wrote down `step_health()` (the L4→L3 feedback) as **"the core of this module"**, and
**it did not execute in any of 81 runs.** The cause was a name mismatch —
`run.Guard` computed `dt·|F|max` and put it only inside `l4`, while `tools/health.py` was looking for
`numerics["step_rms_sigma"]`. When absent, it printed `"not measured"` and
**returned HEALTHY as it stood**.

Three lessons —
① **Silence is not success.** `82/82 HEALTHY` read like coverage. Now the number of unmeasured runs is
   printed separately, and it is stated explicitly that that HEALTHY means only "no divergence, no freezing, no
   collapse" and **not "dt is small enough"**.
② **Measure from the force.** `dt·|F|max/γ` is itself `dt/τ_fast`, so there is nothing to subtract.
   A position difference has to have the thermal noise `√(2·dim·dt)` subtracted in quadrature, and if the drift
   is 0.5% of the thermal noise, in a finite sample **the drift gets truncated to 0** (measured — `verify_health.py` ②b).
③ **Keep the worst value.** It used to store only the last guard sample. In measurement the maximum force was
   1062.9 against a last value of 244.2 kT/σ — a **4.4×** difference. A stability verdict is not the mean but the worst case.

**A second instance of the same pattern — the pre-run gate**: `health.gate()` judged on `verdict != "PASS"`
and so rejected `"PASS (3 warnings)"`. Measured, **80 of 83 specs were falsely rejected**,
and among them the number of genuine hard failures was **0**. bd-physics §4 defines statistics and finite size as
⚠ (not ❌) and `run.execute` was correctly looking at `startswith("FAIL")`.
**The reason the two were out of step and nobody knew is that `execute` does not call `gate()`** —
a checker that is not wired up does not reveal itself when it is wrong. What blocks now is only ① a hash mismatch
② `FAIL` ③ an L3 integrity **error**, and warnings and thin margins are **shown but not blocked** by `gate_notes()`
(they must be printed, because passing silently would make the gate meaningless).

**⛔ Step coverage for legacy runs cannot be filled in retroactively** (confirmed 2026-08-05). Both routes are closed —
① **re-running**: `run_id` is content-addressed, so as the code evolves the hash changes and **a new run under a
different id** is created while the old run is permanently unmeasured. Measured, **all** 79 unmeasured runs failed to
reproduce with the current code (`trap-drag`'s tag convention `tr0.117647-v4`→`v4`, and `trap-2d-5um` changed 3 times
as seed was added to the spec). ⚠️ The presence of `specs/<run_id>.json` is **not grounds for reproducibility** —
stale specs are left behind. ② **GSD replay** (recomputing only the forces from the trajectory) is **invalid under
time-dependent driving**: in `trap-drag` the replayed |F|max=17041 against the guard's true value of 1062.9, a
**16.03× overestimate**, and the frame median is already 1.58e4, so it is a systematic error (the trap anchor is fixed
at `t=0` — decomposed and confirmed as `HarmonicTrap` 1.64e4 vs pair forces 14~77). `sim.timestep` cannot be changed
after state creation.
→ So the coverage denominator is taken not as "what can be filled" but as **"runs executed after the measurement was
wired up"**, and the legacy is stated as **not retroactively fillable** (`scratch/probe_gsd_replay.py`, 2 KB `tooling` entries).

The first run in which the feedback actually turns: `dt/τ_fast` measured 1.76e-4 vs the L3 prediction 1e-2 → **0.02×**.
The ledger is complete (no fast scale missing), so it means **the design is 57× over-conservative**, and it is
reported as `dt margin (cost)` — not asserted (the maximum between guard samples may be missed, and the design's worst
approach distance may not have been reached). Cost has been a recurring problem in this project (a sweep 25 days→1.16 days).

**Two things L3 caught before running** (skill `bd-physics` §6.2b, §6.2c):
`trap-drag` passes all 7 hard checks and yet **the statistics do not come out** — with SNR=0.0985 the precision of
`⟨F_drag⟩` over one crossing is 5.5% (the target is 2%) and the lattice period only occurs 17 times.
For `chain-bend`, **the mode that sets dt is not what is being observed** (τ_fast/τ_chain = 2.2e-4).

### chain-bend's L3 gates (2026-08-05) — the case where the check itself was wrong ⭐️

**Running the two gates revealed that 5 of the spec's derivations were wrong.**
`scratch/verify_chain_bend_gates.py` (the gates) · `scratch/verify_angle_matrix.py` · `scratch/bench_chain_bend.py`

| What | The wrong value | The right value |
|---|---|---|
| the governing timescale | τ_chain = γ/κ_center | **τ_max = γ/λ_min** (9.18× longer) |
| De | ωτ_chain → "0.1~10" | ωτ_max → **0.99~91.8** |
| equilibration | 20 τ_chain (= 2.2 τ_max) | **20 τ_max** (K* shifts by 21%) |
| the SNR numerator | the driving amplitude `a` | **the response `|ŷ(ω)|`** (up to a 60× difference) |
| "the stiffness the trap feels" | 48EI/L³ (rigidly clamped) | **κ_drive = ×0.682** (both ends are trapped) |

Three lessons —
① **Setting dt from λ_max while not looking at λ_min** is the trap. It comes free from the same eigendecomposition,
and it is the governing scale. ② **The numerator of the SNR is not the driving but the observable's response**
(the same family as `trap-drag`'s tautology hole — the checks pass and the statistics do not come out).
③ **Deliberately putting the stiffness ratio near 1** (κ_end/k_t = 0.463) **guarantees that the boundary-condition assumption breaks.**

**What was resolved**: whether `τ_p/τ_fast = 0.60` (the fastest mode is not overdamped, ζ=0.65) contaminates the
observables was compared with `OverdampedViscous` vs `Langevin(kT=0)` across all 7 ω → **at most 0.159%**. It is fine to measure with BD.
★ A thermal comparison (Brownian vs Langevin, kT=1) **did not even have the power** to exclude less than 47% —
because `|ŷ|/ℓ_k < 1`. **When testing an integrator assumption, use a kT=0 deterministic difference** (noise 0,
the transient cancels as a common mode).

**★★ The cause of the 28% mismatch = a HOOMD bug (identified 2026-08-05)**: it is not the linear-response prediction
but **HOOMD that is wrong.** `md.angle.Harmonic`, while moving the torque into coordinates via `1/sin θ`,
**clamps `sin θ` at √2×10⁻³** (measured SMALL=1.414217e-3, standard deviation 1.4e-7).
Below that the force is scaled down by `sinθ/SMALL`, becoming `∝ κ(θ−π)²` — **quadratic rather than linear**.
`t0=π` has `sinθ=0` at equilibrium itself, so it bites **hardest in a stiff chain**.
**The energy is 0.0000% accurate throughout**, so an energy check does not catch it (traps 15·16·17 →
skill `bd-hoomd`). It was narrowed by a ladder — dt (the exact solution in the z region, 0.002%) · nonlinearity, the
x degree of freedom, bond stretching (an exact scipy minimization agreeing with the linear model to 0.32%) · the energy
(0.000%) were excluded in turn, leaving only the force.

⛔ **Conclusion: `chain-bend` cannot be run with `angle.Harmonic`.** **All 23 angles** of the response profile
are in the broken region (min|θ−π| = 7.3e-5 = 1/19 of SMALL), and even the largest angle is only 2.0× SMALL.
Trying to escape by raising the amplitude is blocked by `δ_max` (M<M_c).
A hard check was added, so **`nondim spec` refuses to write a spec** (0 specs is the normal state).
The workarounds: ① implement it directly with `force.Custom` (accurate but 26× slower) ② lower κ₀ to soften the
chain ([P2], a surfactant) ③ **this region is linear, so solve it analytically without MD** — the exact minimization
agrees to 0.32%, so with an OU process even G'(ω) is closed-form.

**What remains**: the lower bound of `omega_range` does not reach the quasi-static plateau (De(ω_min)=0.993, against
a criterion of 0.1). It is a tier 3 user-approved value in system.yaml, so **it was left unchanged and exposed by a check**.

**The cost was not the physics but the implementation**: `md.force.Custom` calls into Python every step and is
**26×** slower at n=25 (2,339 vs 61,264 steps/s). Switching to a ghost particle + `bond.Harmonic(r0=0)`
(= exactly a ½k r² harmonic trap, on the compiled path) gives 55,551 steps/s → the sweep goes from 25 days to
**1.16 days** (7 runs in parallel, 13.7 hours). ★ Before ruling "too expensive to run", measure which side the cost is on.

### 1-D · `chain-bend-2d-dlvo` — the alternative hypothesis built out of doubting the papers' beam-mechanics reading (2026-08-06) ⭐️

A branch the user started saying **"[P1][P2]'s beam-mechanics reading may be wrong; I want to check it as
colloidal beads and springs rather than a monolithic beam"**. The same 3-point-bending geometry as
`chain-bend-2d-oscill`, but with the bond changed to **the secondary minimum of a DLVO central pair potential**,
and **the angle (bending) potential deliberately left out** — that is the hypothesis.

★ This hypothesis is **exactly the same as the papers' own argument**. [P1] p.2: *"Such assumptions have
been made based on DLVO interactions between particles, which are centro-symmetric.
If particles did undergo free rotations, we would expect the aggregates to respond to
the bending moment by forming a **trianglelike structure**"* — the papers introduced JKR adhesive contact because
this prediction differed from what was measured. This case confirms **the prediction side** by running it.

**A structural fact derived before running (G1)**: with a straight chain + purely central forces + the bonds at their
natural length (U'(ℓ)=0), the coefficient of the O(y²) term in the transverse displacement is U'(ℓ), so
**the linear bending stiffness is exactly 0**.
So `K_prime`'s prior prediction was fixed at **0** (principle 9.2).

**The result — as predicted** (137 runs, ensembles of 6~9 seeds):

| Sweep | The cleanest condition | K' [kT/d²] |
|---|---|---|
| ω, 7 points (10~30000) | ω=10 (quasi-static, 4.6σ) | **41.8 ± 9.0** |
| chain length, 4 points (n=5~25) | n=25 | **12.7 ± 4.0** |
| amplitude, 5 points (a/d 0.07~1) | a=1d | **9.2 ± 2.6** |

Both sweeps **decrease monotonically and converge to 0**, and the more precise it gets (the error shrinking) the closer to 0.

**The JKR control group — the same geometry, the same seeds, the bending term alone ON/OFF** (n=9, ω=3000, a=632nm, 6 seeds):

| | DLVO-only | DLVO+JKR | Ratio |
|---|---|---|---|
| K' | 199.5 ± 69.9 | 43064 ± 2516 | **216×** |
| bow | 0.1135 ± 0.0048 d | 0.00639 ± 0.00011 d | 17.8× (**22.3σ**) |
| off-diagonal displacement correlation ⟨δy_i δy_j⟩ | **−0.015** | **+0.905** | — |

The response amplitude agrees with the linear-response prediction to **+1.3%**, which becomes independent evidence
that the implementation is right.
⚠️ JKR's K' is 4.5× the static-limit prediction (9628), but **that is not an error; it is because De=10.7 so it is
not quasi-static** (K''=120022 > K' corroborates this). A quasi-static quantitative comparison needs ω<28 rad/s,
which at 25 hours per seed is **not run**.

**Three observable lessons** —
① **The bow (the maximum deviation from the straight line joining the two ends) is the best discriminator.**
The trap has finite stiffness so the whole chain translates, and measuring by the y displacement lets that rigid-body
motion bury the bending (the full y range 3.4× vs the bow 19.6×). ② The `shape_localization` I made at first had
**0 discriminating power** (1.481 vs 1.470 = 0.1σ) — even though the difference is obvious to the eye. An indicator,
once made, has to be **run against the two extremes** to confirm it separates. ③ The tangent correlation ⟨t·t⟩ is
weak in this system — both are nearly straight along the x axis so cos θ≈1.

**Two driving modes** (`--kt-scale` / `--drive-mode`) — added after the user pointed out that in the default settings
**the driven bead follows the trap position by only 3.5~29%**:
`k_t`×100 → following 82% (almost no dt cost. ★ from ×300 up, the parenthesis of `K′=k_t(ŷ_c/ŷ−1)` shrinks to
0.09, so **raising the following further makes the measurement worse**) · `--drive-mode position`
→ force the bead's y directly, following **100%**, corresponding to strain control in rheology and the same as
[P1][P2]'s experimental protocol (move the middle and measure the force at both ends) → it produces `K_transfer`
(**a different quantity from the trap mode's K′** — direct comparison forbidden).

**★★ The 3-protocol synthesis — the conclusion is independent of the driving mode** (145 runs completed, `viz_three_cond.py`):

| Protocol | Following | DLVO stiffness from 0 | JKR/DLVO | Bow ratio |
|---|---|---|---|---|
| trap, default k_t | 27% | 2.28σ | 401× | 15.4× |
| trap, k_t×100 | 100% | **0.84σ** | 12,253× | 1.4× |
| position forcing | 100% | **0.66σ** | 2,525× | 1.9× |

**The better the following, the more DLVO converges to 0** — the default trap's "finite K′=105.7" was
not real elasticity but **an artefact created by the following failure**.
⚠️ Note that the bow's discriminating power collapses 15.4→1.4× — **when the deformation is forced, the shape
is not information**. The earlier paragraph's "the bow is the best discriminator (22.3σ)" holds **only for a weak trap
(free deformation)**.
The rule: free deformation, the shape; forced deformation, the force.

**★★ The system level (the whole system's energy budget) gives the strongest statement** (`system_moduli.py`):
measuring the total dissipation over one period, `K″_total = ∮F·dy/(π|ŷ|²)` — a bead alone (no chain) ωγ=18453 ·
**attaching the DLVO chain still gives 18380 = 0.996×** · JKR 75590 = 4.10× (the chain accounts for 75.8%).
⟹ **The DLVO chain is rheologically invisible** — looked at as a whole system, the dissipation is the same as
just shaking one bead in water. ★ A cross-check with no free parameters: the energy integral agrees with
(the lock-in `K″_chain` + ωγ) to **DLVO −0.11% · JKR −0.17%** (because `lockin.k_star` subtracts −iωγ by
definition — the initially apparent −43%/+32% discrepancy was exactly ωγ).

**Per-bead A·sin(ωt+φ) fitting** (`fit_mode_profile.py`) — JKR agrees with the linear-response prediction to
**+0.01~0.04%** in amplitude and **within 0.06°** in phase (no free parameters). But since
`force.Custom` uses the same `A_bend`, this is an **implementation_check** and not a physics verification.
**The stress propagation speed** (`propagation_speed.py`, the phase gradient `v=ωℓ/|dφ/ds|`, ω=3000): trap default
31,227 · k_t×100 **189,000** · position 30,144 µm/s. The default and position agree to within 3.6%
— it is not the driving mode but **the boundary condition (the stiffness of the traps at the two ends) that sets the speed**.
For DLVO the coherence is 0.46~0.59, so **v is not defined** (there is no stress to propagate). ⚠ Do not read the
amplitude decay length λ as a "penetration depth" — the static beam shape of 3-point bending dominates (at k_t×100
λ=2.6 while 1/k=42.5. For a diffusive form λ≈1/k).

**G′ and G″ cannot be extracted from a single chain** (`viz_moduli.py`). The only legitimate thing is the
`K*→EI→κ₀` conversion (JKR's 56.6 mN/m is a **round-trip confirmation** recovering the input 64 — and DLVO's −0.005 is
4~5 orders below).
GSER `G*=K*/(6πa)` is **invalid** — what the bead feels is not a medium but 2 neighbouring beads.
Producing a real G′(ω) requires building a **network** of chains and computing it directly from the stress tensor
(exactly as chain-bend-2d-oscill observation.yaml F1 forewarned).

**4 tooling incidents this session** (all KB `tooling`) — ① **`result.txt` is written by the case script**
(not by `bdbot.run`). Not writing it makes `status` count 0 runs, prevents skipping completed work, and makes
"clean up the incomplete" **delete completed runs** (6 were actually destroyed). ② **Editing code mid-batch**
wiped out 11 runs with `NameError` — xargs prints done even on a crash, so **the metrics.json count** has to be
counted (48/48 done and yet 37 metrics). ③ **τ_k missing from the dt candidates** — around `--kt-scale` 200 the
trap becomes the fastest mode and dt had been left as it was. ④ **Adding a spec field changed all 137 runs' ids** —
aggregation has to group by tag, not by id.

### 1-E · `network` — a 3D colloidal network (in progress, 2026-08-06)

The sketch ([intake/network/](../../intake/network)) has **not a single number in it** — only symbols (A, ω) and model
names (DLVO, JKR), which is why `missing_required` is unusually long (the result of not inventing). The properties are all
inherited from `chain-bend-2d-dlvo`. **It follows rule 8's 2-stage structure exactly**:

- **Stage 1 (now)** — as far as making the network. Scatter → aggregate via the DLVO secondary minimum → compress → φ=0.10 →
  after relaxation, measure the **structure** (z · independent loops · free ends · d_f · percolation · g(r) · the collapsed-bond fraction).
- **Stage 2 (not started)** — force one particle as `x(t)=A sin(ωt)` and sweep the driving **direction**.
  The driving, the stress propagation and `G′(ω)` are here.

There is one more reason for the split besides rule 8 — **HOOMD's bond/angle are a static topology and cannot be
declared during aggregation.** After gelation finishes, the contact list has to be extracted and the topology frozen for
stage 2 to happen (which is why the JKR bending is also stage 2).

**The 1 completed run is not compression but the `sprout` topology generator** (N=512, max_coord=4, 20.85 million steps,
1757s). z=4.92 · independent loops 748 · free ends 2.34% · d_f=2.51 · collapsed bonds 0. L4 is OK, but
**⟨U⟩/N drifts by 6.9 block SEMs over the relaxation stretch**, so the equilibration is insufficient. ⚠️ This topology
was **imposed** rather than created by DLVO dynamics, so the topology dependence only becomes visible by comparing
against the compression route (rule 7′) — and those 2 compression runs are exactly the ones left incomplete.

★ **The constraint on the compression was measured before running** (`scratch/verify_3d_boxresize.py`):
`update.BoxResize` affinely scales the coordinates (error 8.9e-16), so it also shrinks the bond length of pairs that
are already bonded. If the linear strain per trigger exceeds 0.703% the pair is pushed inside the barrier and
**collapses irreversibly into the primary minimum** (0.40% holds / 0.80% collapses). Hence 0.4% per stage · 178 stages.

### 1-F · `chain-relax-2d-dlvo` — the undriven static counterpart of chain-bend-2d-dlvo (a pilot, 2026-08-06)

A case that started from a user request ("implement a chain system as a situation with no inter-particle friction and
only attraction"). That physics is already what `chain-bend-2d-dlvo` has (DLVO secondary-minimum attraction only, no
explicit bending or friction), but that case started with trap + oscillation driving on top from the beginning —
**it skipped the "static stage" of rule 8 (the static system first)**. This case fills in that missing stage:
the trap, the ghost particle and the driving are all removed, leaving a minimal configuration of nothing but the
DLVO table potential (imported as it is from `chain_bend_dlvo_2d.py` — the same reuse precedent as network) + the
WCA core + Brownian.

Two experiments (`--init straight` / `--init kink`). The kink gives an exact turn angle to just the one middle bond
and **starts every bond exactly at its natural length** (so as not to mix in a stretching signal), perturbing pure
bending only.
With no trap the chain rotates and translates freely, so shape descriptors (the bow, the turn angle) must be computed
in **the chain's own body frame** (referenced to the axis joining the two ends) — chain-bend-2d-dlvo's `_bow(y)`
used lab coordinates as they were because the trap fixed the orientation, but that will not do here.

★★ **The golden test came out 4.6× outside the harmonic approximation — not a bug but the asymmetry of the well**.
Predicting the equipartition of the bond stretching (the same kind as a trap's `⟨x²⟩=kT/k`) naively as `kT/k_bond`
had the smoke run off by 6.67σ. The DLVO secondary minimum is an **asymmetric well**, steep on the inside (towards
the barrier) and far softer on the outside (h→∞, U→0⁻), so the harmonic approximation underestimates the true variance —
Boltzmann-integrating over the whole basin inside the barrier gives **4.57× the harmonic approximation** as the true
predicted value (the same trap as soft-r3's "the Einstein cage approximation goes off through anharmonicity",
bd-physics §6.2). Correcting the prediction to this integral made the error converge 6.67σ→2.15σ as the samples
grew (the smoke run's 50 → the full run's 20,000 τ_bond) — **no transplant bug**, just insufficient statistics (`bond_variance_boltzmann()`).

The n=9·1-seed pilot result (2 million steps, wall clock 1.8 seconds — very cheap, because with no trap and no driving
the only stiffness candidate is the bond stretching): the bond-stretching variance agrees with the corrected prediction
above (−13.7%, 2.15σ). The bending-angle thermal fluctuation is ⟨dθ²⟩=0.127 rad² (σ≈20° per bond, no prior prediction —
measurement, G1 says nothing about the magnitude). The kink (0.3rad) release experiment (300 thousand steps) had
bow_rms 0.334d→0.313d over 3000 τ_bond, **only 6%** — a qualitative signal that it is closer to
diffusive dissipation than to elastic recovery (a collapse to 0), but **with a single seed it is not decisive**.
The next step (awaiting user confirmation): settle the bow(t) relaxation curve statistically with an ensemble
(multiple seeds), and if necessary raise the signal with a larger kink angle.

### 1-D follow-up · "can friction or tension do it instead of JKR" — all 3 alternatives excluded (2026-08-06) ⭐️

A branch that came out of two user questions: **"is there no way to implement friction directly"** ·
**"I would think that even with no friction, a sufficiently strong DLVO would show similar motion"**. All three were
answered by running. The figure: `runs/_bending_model_compare/contact_mechanisms.png` (4 panels).

| Candidate | Implementation | Does it create a bending stiffness |
|---|---|---|
| `md.pair.friction` (3 built-in) | exists, on the compiled path | ⛔ **it cannot** — dissipative + exempt for rolling |
| rolling resistance + a tangential spring | `force.Custom(aniso=True)`, new | ✅ rolling only. The tangential is **exactly 0** |
| a strong DLVO (tension) | existing | △ **only as δ²** — not real bending; in this system 1~5% of the observed value |

**① `md.pair.friction` is powerless in this system** (`scratch/verify_pair_friction.py`, 23/23).
HOOMD 7.1 really does have `FrictionLJLinear`, `FrictionLJCoulomb` and `FrictionLJCoulombNewton`
(Hofmann et al. 2025, arXiv:2507.16388 — run for the first time in this project). Three structural
properties were measured: ⓐ the friction weight is `w(r) = −dU_WCA/dr`, so it is **exactly 0 outside the cutoff** —
the DLVO secondary minimum (`r*=1.00759`) is outside it, so `F_tan=|τ|=0.000e+00` (in the control at `r*=0.95` it is
6.19/18.58/1.00, so there is power) ⓑ **dissipative** — even with a tangential displacement, the force is 0 if `u=0`
(contributing to `K″` only, vanishing in the quasi-static limit) ⓒ **exempt for slip-free rolling** — at `ω=V/2R` all
three are exactly 0. **The bond deformation of 3-point bending is rolling rather than slipping, so this friction cannot stand the chain up.**
⚠️ The sensitivity to residual slip differs — Linear and CoulombNewton are proportional to `u`, but **Coulomb gives
full force even at `u/V=1e−4`** (18.578).

**② Rolling resistance is exactly equivalent to harmonic bending** (`scratch/rolling_contact.py` +
`verify_rolling_contact.py`, 24/24). Setting up
`U_r = k_rR²(2 − (a_i−a_j)·n̂)` with body-fixed contact-point markers gives, in the orientational relaxation limit,
**`κ_θ,eff = ½ k_r R²`** — the 3-point bending stiffness agrees with harmonic bending to **within 1e−5** (n=3~15).
The origin is **the internal particles being shared between both bonds and thereby frustrated** (an end particle relaxes fully).
Conversely, **with a tangential spring alone the bending stiffness is exactly 0** — the orientational unknowns (n)
outnumber the bonds (n−1), so there is always a solution that makes them all 0. At δ×100 the residual shrinks by 1/10⁴
(1.46 → 1.46e−4), settling that it is **cancellation round-off rather than a real stiffness**.
★ **Freezing the orientations makes it different physics** — it penalizes the bond's absolute rotation rather than the
curvature and is `(n²−2n+3)/3` times stiffer (measured n=5·9·15 → 6.0000·22.0002·66.0003).

**③ The dynamic comparison — in this system the two models are indistinguishable** (`compare_bending_models.py`,
a kT=0 deterministic ω sweep, strain control). The crossover frequency is `ω_c = 1/τ_rot = 8.35e6`, and
**the production run's `ω*=18453` is `ωτ_rot = 2.2e−3`**:

| ω* | ωτ_rot | harmonic K′ | rolling K′ | rolling/harmonic | harmonic/static |
|---|---|---|---|---|---|
| 3e5 | 0.036 | 923,451 | 934,173 | **1.012** | 1.011 |
| 1e6 | 0.120 | 1.020e6 | 1.139e6 | 1.116 | 1.116 |
| 3e6 | 0.359 | 1.772e6 | 2.467e6 | 1.393 | 1.939 |
| 1e7~1e8 | 1.2~12 | — | — | 1.20~1.31 | 5.7~9.0 |

The harmonic model **reproduces the static analytic solution (913,561) to 1.1%**, which verifies the MD setup itself.
The production ω is **16× lower than the lowest measurement point**, so the two models are effectively the same.
⚠️ The ratios in the `ωτ_rot > 1` region must not be read — `harmonic/static` is 5.7~9.0 there, meaning it is not the
bending but **the radial bond and the drag that dominate**, so the bending contribution is diluted. On halving dt both
models give +0.50% (the same direction, so the ratio converges far better).

**④ A strong DLVO creates an apparent stiffness through tension — but it is not bending**
(`scratch/dlvo_tension_stiffening.py`, 4/4). In 3-point bending the two ends are held, so pushing the middle by δ
lengthens the path by `2δ²/L` and the bonds **necessarily** stretch. The prediction with no free parameters,
**`K = 4 k_ext δ²/L²`** with `1/k_ext = (n−1)/k_bond + 2/k_t`, agrees with an exact minimization to a
**maximum deviation of 0.03%**, with a log-log slope of **1.9999**.
★ **The discriminator is not the magnitude but the δ scaling** — tension gives `K ∝ δ²` (0 as δ→0), while
real bending is δ-independent (1.008/1.004× over a 3000× change in δ). It is not a contradiction with G1 but **a consequence of G1**.
But in this system, because of two suppressors, it explains **only 1~5% of the observed K′**:
ⓐ the default trap's stretching compliance `2/k_t=3.83e−4` is **50×** the chain's `7.67e−6` → a 51× suppression
ⓑ the DLVO well's **maximum tensile force `F_max = 810.4 kT/d`** (h=14.6nm). At `k_t×100` and δ=0.43d the required
tension is **4.9×** F_max, and at δ=1d **26.6×** → bond rupture. This is consistent with 1-D having measured `K′≈0`
(0.84σ) under those conditions.
⚠️ Do not define DLVO's "escape from the well" by the `U=0` crossing — on the outside it only **asymptotes** to
`U→0⁻` and never crosses, so a `nan` comes out (experienced for real). The rupture verdict is made from the maximum tensile force.

⟹ **The 1-D conclusion is strengthened.** All three candidates for "mimicking bending with central forces alone" were
quantified and excluded, and the reason [P1] had to introduce JKR has been narrowed further —
what is needed is neither friction nor tension but **rolling resistance**, and that is precisely an adhesive contact area.
★ 2 new traps (`bd-hoomd` 19·20): the lever arm of `pair.friction` is `particles.diameter/2`
and independent of `sigma` · `Brownian`'s `velocity` is not 0 but **uncorrelated thermal noise**.

### Literature · 2 books distilled into `kb/` (2026-08-06) ⭐️
2 distillations in `docs/books/`, 24 KB `handbook` entries, and the verification `scratch/verify_book_claims.py`
(**56/56 PASS**). The verifications are marked in three kinds — `[BOOK]` reproducing a number the book itself reports
(is my reading right) · `[DERIV]` checking numerically that the book's own formulas are mutually consistent · `[OURS]` comparing against our values.

| Book | What it gives |
|---|---|
| [**[L] Leal, *Microstructural Rheology of Complex Fluids*** (Cambridge 2026)](docs/books/leal_microstructural_rheology.md) | microstructure → **the bulk constitutive equation**. It gives us not code but **the definition of an observable and that observable's entitlement to exist** |
| [**[W] Welty, *Momentum, Heat and Mass Transfer*** 5th ed.](docs/books/welty_transport.md) | the Newtonian baseline + **a citable property table** + non-dimensionalization methodology. The non-Newtonian classification is **not in** this book (the two are exactly complementary) |

**★ The status of the 1-D conclusion has risen.** [L] §2.2 proves that in creeping flow the hydrodynamic stress is
**exactly 0 the moment the flow is cut off**, so **the viscoelasticity comes entirely from the thermodynamic
component**. That is, "a degree of freedom with no restoring mechanism does not contribute to the stress" is
**a structure common to all complex fluids** — our G1 (central forces + natural length ⟹ bending stiffness 0) was a
derivation for our system, and it turned out to be an instance of a general theorem.
The dissipation ratio of 0.996 is no coincidence, and it is also why [P1] had to introduce JKR.

**★ "G′ cannot be extracted from a single chain" is also a theorem** — a bulk stress is by definition ⓐ a volume
average containing statistically many particles and ⓑ Batchelor's formula presupposes that the particles are
**force- and torque-free**. A trap is an external force so it breaks ⓑ, and a single chain breaks ⓐ.
**Our `K′` and `K″` are stiffnesses, not stresses, and renaming them `G` makes them wrong.** For a real `G′` you want
**periodic boundaries + a network of chains + shear imposed without traps + the Kramers stress**,
`T^(p) = n⟨F_s R⟩ − nkT I` — [L] Ch2+Ch11 gives every formula needed.
⚠️ Leave out `−nkT I` and it does not go to 0 at equilibrium, and that gets misread as **spurious elasticity**
(whether it is 0 at equilibrium is the golden test). To avoid **double-counting** the Brownian effective force
`F_Br = −kT∇ln P_N` with the instantaneous random force, it has to be compared against a `kT=0` run.

**★ 3 convention traps** (all execution-verified) — ① **the numerator of Pe/Wi is `|E| = √(2E:E)`**, not
`|∇u|` (vorticity cannot change an isotropic equilibrium). Using `√(E:E)` gives a **√2 factor**, a pure rotational flow
has `|∇u|`≠0 and yet **Pe=0**, and for a spheroid it is `|E|(r−1)/D_r`. ② **In SAOS the normal stress difference
oscillates at `2ω`** — our `lockin.k_star`, which locks in at `ω` only, **cannot see it in principle**. And `η′=G″/ω`,
`η″=G′/ω` are crossed over. ③ The orientational relaxation of a dilute axisymmetric particle is **`1/(6D_r)`**
(the pole of the loss term `De/(36+De²)` is at De=6, deviating from a single Maxwell mode with `λ=1/(6D_r)` by <1e-12).
**But the 6 is the 3D l=2 coefficient** — applying it to `abp-rod`'s 2D `⟨cos Δθ⟩` is wrong.

**★ Our `η=0.851 mPa·s @300 K` survived, but `T=300 K` is the weak point.** Log-linear interpolation of [W]'s
appendix I gives 0.8580 mPa·s → **+1.03%** in support (linear interpolation gives +2.91% — **the interpolation method
changes the result**). But **water's viscosity has a temperature sensitivity of 2.06 %/K**, and `T=300 K` is a value
inherited from 1-A while 1-A's sketch had no temperature — **it is written as tier 1 but it is a choice, not a measurement**.
If the truth is 298 K the η error is −4%, and at 293 K −14%, with τ_B following directly. The qualitative conclusion
does not change, but it has to be exposed when comparing quantitatively against the literature.

**★ 2 new case candidates** — ① **the `D_r` suppression of semi-dilute rods**, `D_r = β·D_r,dilute/(nL³)²`, β=1.3e3.
The suppression is **an excluded-volume (steric) effect and not HI**, and in the slender limit HI does not modify
Jeffery's rotational convection → **a rare case where HI-free BD is the right tool**, and being a quantitative
prediction with a single free parameter it becomes **the first rod case to have a real `hypothesis` (rule 7′)**
(`abp-rod` had 0 discoveries because all 5 of its predictions were `implementation_check`). ② **the flow distortion of
the hard-sphere pair distribution** — the **sign prediction** that `g` is high in the compressional quadrant and low in
the extensional one (Batchelor 1977), directly measurable with 2D BD.

**★ Boundary numbers** — dilute↔semi-dilute is `nL³=O(1)` in theory but **in experiment 10~50, 30 by convention** →
when putting it in a check it has to be **a warning (⚠)** rather than a hard failure (❌). A dilute suspension of spheres
is **exactly** Newtonian (if the structure is not changed by the flow it is Newtonian), and the `φ²` term reaches 10% of
the Einstein term at **φ≈4%**.

⚠️ [L] **does not treat microrheology or GSER at all** (0 hits in a full-text search) and **does not treat active matter
either** (the thermodynamic component is defined as "the mechanism that returns it to equilibrium"). `abp-rod` is
outside this book's frame, and it **neither conflicts with nor supports** our GSER-invalid conclusion (a different subject).

## The run environment

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
```
conda env `simulation_bot` · hoomd **7.1.0** (CPU, **no MPI**, no GPU) · macOS arm64
Defined in: [environment.yml](../../environment.yml)

## Absolute rules

1. **Dimensions come first.** Every system is fixed as an SI physical system first, then non-dimensionalized
   through a scale table. There is no path that starts from dimensionless values. If you must start from
   dimensionless values, the anchors (particle diameter, temperature, viscosity) must be stated explicitly.

2. **Do not hand-write a dimensionless spec.** Only what was derived from a physical system gets run.
   `specs/` and `runs/` are written **by tools only** (`bdbot.cli nondim spec`). A hand-edited spec is caught by
   `LoadedSpec.verify_hash()` — run_id is the hash of the content, so editing the content breaks it.

3. **Attach a provenance to every parameter.** State which it is: sketch, literature, handbook or estimate.
   If you do not know, say you do not know. **Never invent one.**

4. **Do not write HOOMD code from intuition.** There are several traps that are quietly wrong.
   Read skill `bd-hoomd` **before** writing the code.

5. **When reading a sketch, transcription comes first.** Copy down what is visible, then interpret.
   What could not be read and what is ambiguous must be stated. A value absent from the sketch stays `null`.

6. **Verify a physical claim before stating it.** This project has already been wrong several times
   (a missing minimum image → +1856% / "bind it rigid and it becomes anisotropic" → measured 1.000000 /
   the `pair.Table` grid convention → the force −1.65% / a cage-fluctuation formula missing √2·`a_NN`).
   Claims about HOOMD behaviour and about physics formulas are confirmed by execution, not by reasoning.

7. **Isolate independent elements and verify them one at a time.** If a system has several mutually independent elements,
   first run the **minimal configuration with only one element on**, compare against the analytic solution, and only then combine.
   Often only the isolated configuration has an analytic solution — sometimes isolation is the only
   way to make a ground truth. If the comparison fails, **split further** (ABP → Λ separately, v₀ separately).
   Details: mater_plan.md principle 9

7. **"Verified" does not mean it agrees with the prior hypothesis.** ⭐️
   The reason to compute these systems is that **they may differ from the prior hypothesis**. So calling a mismatch
   with the prediction a failure means calling discoveries failures. Attach a **role** to every comparison
   (`bdbot.metrics.ROLES`):

   | Role | Where the prediction comes from | What a mismatch means |
   |---|---|---|
   | `implementation_check` | **derived from the model I implemented** | **a bug** → fix it |
   | `hypothesis` | an assumption the simulation does **not** impose (continuum, dilute limit, effective medium, literature) | **a result** — report it |
   | `measurement` | no prediction | the simulation is the answer |

   When designing a case, **write down separately "the assumptions the simulation imposes" and "the assumptions the
   theory adds".** If the latter is empty the case is for validation and has no discovery —
   `abp-rod` actually was that (all 5 predictions `implementation_check`, 0 hypotheses).
   Details: skill `bd-physics` §7.5

   *(That there are two number 7s is not an error but a convention — the isolation one is cited as `rule 7` and this one
   as `rule 7'`. There are already 6 references, so they are not renumbered.)*

8. **Build the static system first, and add motion afterwards.** ⭐️
   **Before putting in a time-dependent term** like driving, oscillation or dragging, build the static configuration
   (the equilibrium arrangement, the bond structure, the energy minimum) and finish verifying with that alone.
   It is the time-axis version of rule 7 (isolation), and in this project there are **real cases where a static
   baseline caught a dynamic bug**:

   | What was done statically first | What it caught |
   |---|---|
   | an exact scipy minimization agreeing with the linear model to 0.32% | **`angle.Harmonic`'s force is up to 96% wrong** — the energy is 0.0000% accurate, so a dynamic run alone could not catch it |
   | `trap-2d-5um`'s `⟨x²⟩=kT/k` (static equilibrium) | the golden test of the trap implementation — before the PSD (dynamic) |
   | the eigendecomposition of `chain-bend`'s stiffness matrix | λ_max, which sets dt, and **the governing scale λ_min** differ by 9.18× |

   Why it has to be this order — ① a static problem **has an analytic solution or an exact minimization**, so a
   ground truth can be made. A dynamic response usually does not. ② the noise is 0 so systematic error is not buried
   (a thermal comparison did not even have the power — the same root as the lesson to use a `kT=0` deterministic difference).
   ③ the static structure is **the initial condition of the dynamic measurement**, so inverting the order mixes two errors inseparably.

   In practice: store the static stage as **its own run** (observables and verdict included) and have the dynamic stage
   read it and start — attach it to one script and the static stage cannot be re-run alone.
   `network` is in this form (stage 1 compression gelation → stage 2 driving).

9. **Build the system out of direct inter-particle interactions — do not substitute a macroscopic general formula.** ⭐️
   If an interaction is judged physically necessary, do not defer it to an approximate expression or a substitute
   relation but **put it straight into the simulation as a force that actually acts between particles — a pair
   potential, a bond, an angle.** The opposite direction is equally forbidden — to measure a macroscopic property
   (a modulus, a diffusion coefficient, a stress) of a system already built from inter-particle forces, do not take
   the shortcut of substituting into a **continuum, mean-field or literature general relation** (GSER, effective-medium
   theory, a continuum beam formula, and so on). Such relations have preconditions (a continuous medium, the dilute
   limit, rigid joints), and this project's systems break them often. **Get the system's macroscopic answer by running
   the inter-particle interactions** — do not derive it backwards from a formula.

   A real case — `chain-bend-2d-oscill`: GSER `G*=K*/(6πa)` was going to be used, and it was **invalid**.
   Because what the bead feels is not a continuous medium but **only 2 neighbouring beads** — for a real
   `G'(ω)` you have to build a **network** of chains and compute the stress tensor directly from the inter-particle
   forces (the conclusion already written in this file). `chain-bend-2d-dlvo` is the same
   principle applied from the opposite direction — rather than taking the bending stiffness `κ_θ=EI/ℓ` that
   [P1][P2] used (a general relation, namely continuum beam theory) as it stands, it built the chain out of
   **a single DLVO central pair force** actually acting between the beads and watched, by running, whether bending
   came out. `abp-rod`'s anisotropic translational friction also shows why this rule exists — BD has no
   inter-particle hydrodynamic interaction (HI), so only the isotropic average friction γ̄ can be used, and as a
   result the anisotropy **cannot be mimicked even with a general formula** (impossible in BD — mater_plan §20 question 10).

## Detailed knowledge

| When you need | What to read |
|---|---|
| using the shared code | [`bdbot/__init__.py`](../../bdbot/__init__.py) — 20 modules, what exists and what **does not** |
| reading a sketch or an image (**always first**) | skill **`bd-intake`** — transcription first + 8 rules against inventing |
| writing HOOMD code (**always first**) | skill **`bd-hoomd`** — 17 traps + verified snippets |
| defining a physical system, proposing parameters, non-dimensionalizing | skill **`bd-physics`** — how to write a scale table + 2 cases + a comparison table |
| what can be built | [docs/hoomd_capabilities.md](../../docs/hoomd_capabilities.md) — the measured capability matrix |
| **defining a stress, modulus or rheology observable** | [docs/books/leal_microstructural_rheology.md](../../knowledge/source/books/leal_microstructural_rheology.md) — **the [L] Leal 2026 distillation**. The two components of the stress · the Kramers form · the Pe/De conventions · rods |
| **when you need a source for a property value or a non-dimensionalization method** | [docs/books/welty_transport.md](../../knowledge/source/books/welty_transport.md) — **the [W] Welty 5th ed. distillation**. The water η(T) table · Stokes/Re · Buckingham |
| looking for past runs and lessons | `$PY tools/kb.py query --tags ... --origin ... --kind ...` |
| **leaving behind what was learned at the front end** | `$PY tools/kb.py add --origin intake\|paper\|tooling\|method\|handbook ...` |
| the full design and roadmap | [mater_plan.md](../../docs/history/2026-08_simulation_auto_master_plan.ko.md) |

## Working practice

- **Case-driven**: do not build the framework first. Drive one case end to end, and
  abstract only what has appeared twice (mater_plan.md §16 Phase 1).
- **Verify together**: scale tables and non-dimensionalization results get looked at with a human. Until knowledge
  accumulates, do not settle a physics judgment alone.
- **Keep what you learned even when there is no run.** Sketch readings, literature distillations and tooling lessons
  attach to no run, so there is no place for them in `record.json` → they go into `kb/entries/`:
  ```bash
  $PY tools/kb.py add --origin intake --kind pitfall --source "file#anchor" --claim "..."
  ```
  `origin`: `intake` (sketch reading) · `paper` · `tooling` · `method` ·
  `handbook` (**books and references** — the distillation goes in `docs/books/` and the entry's `source` is written as
  `distillation#section ← [abbrev] p.page` so it can be retraced)
- **Fix a role (`role`) for every observable.** The default is `measurement` (no verdict).
  Only an `implementation_check` mismatch is a FAIL; a `hypothesis` mismatch is **reported as a result**.
- **Show results as graphs and animations.** ⭐️ Do not report in words and tables alone —
  the default is **to produce both and attach them** (mater_plan.md §13.3).
  - **Graphs**: measurement vs prediction on the same axes. Overlay the check thresholds, the analytic solution and
    the literature values as lines. "What is wrong" has to be visible — [`scratch/viz_chain_bend.py`](../../verify/viz_chain_bend.py)'s
    6 panels are the model (the 5 spec errors correspond one to a panel).
  - **Animations**: how the system actually moves. Putting `kT=0` (the mode shape) and `kT>0` (how far it is buried in
    thermal fluctuation) **side by side** shows an SNR problem faster than a number.
    An animation may be made cheaply with a large `dt`, but **it must be stated that it is not the production measurement**.
  - ⚠️ **Write graph labels in English, not Korean.** matplotlib's default `DejaVu Sans` has no Hangul, so the labels
    all become `□`, while the fonts that do have Hangul (`AppleGothic`, `Apple SD Gothic Neo`,
    `NanumGothic`) are conversely missing symbols like `−` (U+2212) and `ŷ` (U+0177) (measured) — do not try to fix it
    by switching fonts; **write the axes, legends, titles and annotations in English from the start**.
    The physics and statistics terms in this field are often awkward or unfamiliar in Korean in the first place
    (for instance `de-correlation time`, `shear thinning`, `yield force`), and the English conveys the meaning more
    precisely and concisely — it is a rule to keep even without the font problem.
    After generating, confirm that the `missing from font` warnings are 0.
- **After a run finishes, always run the post-mortem.** Do not declare success or failure; decide it by measurement:
  ```bash
  $PY tools/postmortem.py runs/<run_id> --lesson "lesson::pitfall::coord=value"
  ```
  → a `record.json` (a tier 3 KB entry) is created. Queries are `$PY tools/kb.py list|query|lessons`.
  A case script must always emit `metrics.json` (the machine-readable result).
- **Leave verification scripts in `scratch/`.** They must be reproducible. All PASS:
  `verify_intake_guards`(27) · `verify_nondim_guards`(33) · `verify_l3_spec_gaps` ·
  `verify_l3_two_cases`(16) · `verify_1c_equivalence` · `verify_skill_snippets` ·
  `verify_pair_table` · `verify_bdbot` · `verify_health`(45) ·
  **`verify_angle_matrix`**(the matrix ↔ HOOMD `angle.Harmonic`, 0.55%) ·
  **`verify_chain_bend_gates`**(gate A the lock-in · gate B′ inertia) ·
  **`verify_book_claims`**(56 — claims extracted from the 2 books in `kb/`. `[BOOK]`/`[DERIV]`/`[OURS]` distinguished) ·
  **`verify_pair_friction`**(23 — the structural properties of the 3 `md.pair.friction` classes + traps 19·20) ·
  **`verify_rolling_contact`**(24 — rolling resistance ↔ harmonic bending equivalence, the analytic force == the energy gradient) ·
  **`dlvo_tension_stiffening`**(4 — the tension-hardening prediction `K ∝ δ²` to 0.03%)
  Benchmarks and visualizations are left alongside: `bench_chain_bend`(the step rate) · `viz_chain_bend`(graphs + video)
  The 1-D analysis scripts: `dlvo_ledger`(the DLVO curve in SI) · `hamaker_crosscheck`(the independent Hamaker check,
  −10.5%) · `analyze_correlations`(positional correlations) · **`fit_mode_profile`**(per-bead A and φ ↔ linear response,
  0.01%) · **`system_moduli`**(the whole system's energy budget, cross-checked to 0.17%) ·
  **`propagation_speed`**(the phase gradient → the propagation speed) · maintenance: `backfill_bow` ·
  `backfill_result_txt` · `resume_dlvo_runs.sh`
- **When you build a checker, deliberately break it and test it.** "Silently passing" and "not checking" are different —
  `verify_intake_guards.py` actually caught 1 crash bug.
- When you find a new trap, **add it to skill `bd-hoomd`** and leave a reproduction script.
