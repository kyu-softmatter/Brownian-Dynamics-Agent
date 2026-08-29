# 06 · Status and roadmap

Status as measured in the merged tree on **2026-08-28**, not as remembered. Where
a number here disagrees with a Korean design document in
[`history/`](history/), this page is the newer one.

---

## 1 · Where it actually is

| Layer | State |
|---|---|
| **L1 agent layer** | ✅ 6 skills · 9 subagents (model-tiered) · 4 rules · `settings.json` that refuses to edit a sealed document. Structure guarded by 64 tests |
| **L2 engine (`bdbot/`)** | ✅ 21 modules, L0→L7. All 8 cases on the `@RUN.builder` contract. The 5 pre-existing cases were re-run after migration and matched to 15 decimal places |
| **L2 pipeline half (`simbot/`)** | ✅ 19 modules, S2/S6/S7/S8 + sealing + `INCONCLUSIVE`. ⚠️ **one runner only** (`passive-sphere--harmonic-trap`) |
| **L3 knowledge** | 🔶 46 wiki pages · 44 distillations · 126 entries · 227 post-mortems — **but two unmerged schemas** |
| **L4 artifacts** | ✅ 278 specs · 261 run directories · 254 with `metrics.json` |
| **Tests** | ✅ **572 pass**, 2 skipped, 15 slow deselected, ~11 s |
| **Verification scripts** | ✅ 74 in [`verify/`](../verify/) |

Physics conclusions are in [04 Cases](04-cases.md). The headline one:
`chain-bend-2d-dlvo`, 145 runs across 3 driving protocols, concluding that a
DLVO-only colloidal chain has **no bending stiffness and is rheologically
invisible** — which is the prediction [P1] itself made before introducing JKR
adhesion, now executed.

---

## 2 · What is blocking progress

Mostly **facts and seams, not missing code.**

### The two-engine seam is the top item

A case run through `bdbot` gets a numerical-health verdict and no **sealed
prediction**. A case run through `cli.py` gets sealing and cannot use any of the
8 cases. So the strongest verification discipline in the repository — write the
answer down first, hash it, refuse to build the comparison table if the hash
breaks — **has never been applied to the physics that actually ran.**

That is not a small gap. Every conclusion in [04](04-cases.md) is defensible on
its numbers, and none of them is defensible against the charge of having been
*interpreted* after the fact. The fix is bounded: give `bdbot.run.execute()` a
pre-run hook that writes and seals a prediction file, and have `bdbot.metrics`
emit into `simbot.validate`'s comparison format.

### The knowledge base has two schemas

`wiki/` Markdown vs `entries/` JSON, queried by different tools, so a lesson
filed in one is invisible to a reader of the other ([03
§1](03-knowledge-base.md#1--two-schemas-not-yet-one)). Until unified, the honest
rule is *query both*, which nobody will remember to do.

### `result.txt` is written by the case script

87 real runs are invisible to `bdbot.cli status`, and the convention has already
caused completed runs to be deleted. `bdbot.run.execute()` should write it.

### `T = 300 K` is labelled tier 1 but is a choice

Inherited from a sketch that had no temperature. Water viscosity is 2.06 %/K
sensitive, so this propagates −5 % (at 298 K) to −14 % (at 293 K) into every
`τ_B`. Qualitative conclusions hold; **any quantitative comparison to literature
must expose it.** This one cannot be fixed by code — it needs the experiment's
actual temperature.

### The literature base is narrow

38 of 42 distillations are the group's own published work, so *"the literature
says…"* is weaker here than it looks ([NOTICE §1](../NOTICE.md)). And
bibliographic data is largely unchecked — confirm a citation before it reaches a
manuscript.

### `dt` gate logic exists in two places

`simbot.nondim.choose_dt` and `campaigns/chain_bend.py` reimplement the same
thresholds, and `choose_dt`'s displacement gate keys off `bool(spec.pair)` while
`SystemSpec` has no bond/angle field — so **a bond-only system silently turns the
gate off**. The measurement that needed it (`max|F*| = 1037.7`, `dt` cut 100×)
came from exactly such a system.

---

## 3 · Cases in flight

| Case | Next step |
|---|---|
| **`network`** | Stage 1 incomplete: the 2 compression-gelation runs stopped with only `traj_A.gsd`, and the one finished run is an *imposed* `sprout` topology whose relaxation drifts by 6.9 block-SEM. Finish compression, then contrast the two topologies (rule 7′). Stage 2 (driving, direction sweep, `G′(ω)`) not started |
| **`chain-relax-2d-dlvo`** | Pilot is n=9, **one seed**. Ensemble over seeds to fix the `bow(t)` relaxation curve; raise the kink angle if the signal is too small. The 6 % decay over 3000 τ_bond suggests diffusive loss rather than elastic recovery, but one seed cannot settle it |
| **`chain-bend-2d-dlvo`** | The JKR `K′` is 4.5× the static-limit prediction because De=10.7 — not quasi-static. A quantitative quasi-static comparison needs ω<28 rad/s, i.e. **25 h per seed**, and is unrun |
| **`chain-bend-2d-oscill`** | Past the `angle.Harmonic` bug via `force.Custom` linear bending, with 15 specs and 7 runs — but no end-to-end script registered, so it is still L3-only. The cheapest finish is ⓒ **solve it analytically**: the region is linear and exact minimization agrees to 0.32 %, so `G'(ω)` closes in OU form with no MD at all |
| **`trap-drag-2d-hex300`** | Complete but reports 0 runs in `status` for want of `result.txt` |

---

## 4 · Two new cases the books argued for

Both come out of the [L] Leal distillation, and both matter because they would be
the project's **first genuine `hypothesis` cases** in their families rather than
implementation checks.

**① `D_r` suppression in semi-dilute rods.** `D_r = β·D_r,dilute/(nL³)²`,
β=1.3e3. The suppression is an **excluded-volume effect, not HI**, and in the
slender limit HI does not modify Jeffery rotation — so this is **a rare case
where HI-free BD is the correct tool**, with a one-free-parameter quantitative
prediction. It would give the rod family a real hypothesis, which
`abp-rod-2d-run-flip` never had.

**② Flow distortion of the hard-sphere pair distribution.** A **sign** prediction
(Batchelor 1977): `g` elevated in the compressional quadrant, suppressed in the
extensional. Directly measurable in 2D BD, and a sign prediction is hard to
accidentally confirm.

**And the real `G′(ω)` route is now specified**, which the single-chain work
could not reach: periodic boundaries + a chain **network** + no traps + imposed
shear + the Kramers stress `T^(p) = n⟨F_s R⟩ − nkT I`. ⚠️ Omit `−nkT I` and
equilibrium stress is non-zero, which reads as **fake elasticity** — being zero
at equilibrium is the golden test. And `F_Br = −kT∇ln P_N` must not be
double-counted with the instantaneous random force; check against a `kT=0` run.

---

## 5 · Deliberately not done

| Not doing | Why |
|---|---|
| Hydrodynamic interactions (Oseen / RPY / Stokesian dynamics) | free-draining approximation. Where it breaks is documented, not hidden. It is also why anisotropic translational friction is **impossible** in this framework rather than merely unimplemented |
| GPU / MPI / cluster submission | the installed HOOMD is CPU-only with no MPI. `N ≲ 10⁴` |
| Chemistry, CFD coupling, full electromagnetics | out of scope |
| Monte Carlo / HPMC | would be a second engine to verify |
| Translating the Korean layers | skills, knowledge base and `history/` stay Korean. Translating 5,000+ lines of dense design history risks corrupting the corrections trail, which is the most valuable thing in it |
| Promoting anything to `bdbot/` on one appearance | the rule is only ever *"has it appeared twice?"* Currently held back on purpose: equilibrium criteria · observables · verification strategy · governing-timescale choice · initial placement · sampling loop |

---

## 6 · Future work — joining this agent to the microscope agent

The longer-term goal is not a better simulator. It is to close the loop between
**what to compute** and **what to measure**, which today are two separate agents
built by the same person for the same systems. The other half is
[**agentic-microscope**](https://github.com/kyu-softmatter/agentic-microscope):
the same architecture pointed at the instrument instead of the integrator. It
takes a research goal, decides what the microscope can actually record, and
returns executable settings with an evidence tier and a per-check margin —
refusing, with the one missing input named, when a gate's input was never
measured.

| | **Brownian-Dynamics Agent** (this repo) | **agentic-microscope** |
|---|---|---|
| Input | a sketch of a physical system | a research goal |
| Decides | what the system does, in silico | what the instrument can actually record |
| Refuses when | a number has no provenance | a gate's input was never measured |
| Produces | a dimensionless spec and a defended result | executable settings, an evidence tier, per-check margins |
| Its knowledge base | system cards, findings, benchmarks, post-mortems | instrument config, calibrations, tacit expertise, decisions |
| Its unit of doubt | `tier` and `derived_from` on every number, and a sealed prediction | `measured` vs `assumed`, and a falsifier on every prior |

The two architectures match because this one was built from that one's lessons:
hard gates that return `BLOCKED` naming the one missing input, a deterministic
core under a thin agent layer, and a knowledge base read before every decision
and written after every verdict. **Neither is finished, and coupling two moving
targets would be a mistake** — so this is future work, with a stated order of
preconditions.

→ the mirror of this section, from the instrument's side:
[agentic-microscope README](https://github.com/kyu-softmatter/agentic-microscope#future-work--joining-this-agent-to-the-simulation-agent)

### Why joining them is worth doing

**1 · Four of its gates currently ask a human for numbers this repo computes.**
Its decision order opens by taking *the physical quantity to measure and the
target precision* from a person, and its next step wants the system's
correlation time τ_c and correlation length ℓ_c — measured if measurable,
otherwise a theoretical estimate marked `evidence: assumed`. Those are outputs
here. They propagate through its committee: **G8** needs `D` or τ_c for the
motion-blur ceiling, **G5** needs ℓ_c and the task kind, **G11** needs a target
error, **G14** needs the trap stiffness κ. Fed from a spec instead of from a
person, four gates stop asking and start deriving — each number still carrying
its own provenance.

**2 · A measurement closes assumptions this repo cannot close by itself.**
The most damaging soft spot here is `T = 300 K`, labelled tier 1 but actually
inherited from a sketch that never stated a temperature (§2) — worth −4 % to
−14 % on every timescale, because water's viscosity is 2.06 %/K sensitive over
the 293–313 K segment that contains 300 K. A
thermometer reading ends that. The same holds for the particle size
distribution, the salt concentration and the surface potential: tier-1 *choices*
here, routine measurements there. Neither side has to lose track of what a
number is, because both already carry a tier and a falsifier on every stored
value — and that repo's `kb/literature/` exists precisely so a value nobody
there measured can let a gate **compute** while never setting
`evidence: measured`.

**3 · The central hypothesis needs both halves, and neither half can settle it
alone.** `chain-bend-2d-dlvo` found that a colloidal chain held together by DLVO
forces alone has no bending stiffness: bow **0.1135 ± 0.0048 d** without
adhesion against **0.00639 ± 0.00011 d** with JKR, **22.3 σ** apart
([04 · 1-D](04-cases.md)). At the bead diameter this case actually used,
d = 1.47 ± 0.01 µm from [P1] p.1, that is a transverse displacement of
**166.8 ± 7.1 nm** against **9.39 ± 0.16 nm** — and the two readings split apart:

| the experimental question | what it requires | against a 10 nm target precision |
|---|---|---|
| is it DLVO or JKR? | resolve the **157.5 nm difference** | **15.7 ×** margin — comfortable |
| is the JKR branch separable from zero? | resolve **9.39 nm** | **0.94 ×** — just under |

So *which model* is settled by physics and is easy; *whether the adhesive branch
is nonzero* is settled by photon count and frame count, at that repo's G11, and
sits right at the precision its own worked examples target. **That is the
question neither repository can answer alone** — and today it is answered by
consulting them separately and trusting that the two `d` mean the same thing in
the same units. (The ±0.01 µm on `d` itself contributes only ±1.14 nm, 0.68 %,
so the diameter is not the limiting uncertainty.)

**4 · Proposing the next experiment, not only checking the current one.**
The same result read the other way: bow separates DLVO from JKR at 22.3 σ under
a soft trap and at only **1.4 ×** under a stiff one, because forcing the
deformation destroys the shape signal (§ rule: *free deformation → shape;
imposed deformation → force*). **Discriminating power is a property of the
protocol, not of the effect.** So a sweep here, scored against that repo's
feasibility gates, ranks candidate experiments by predicted separation *per unit
of instrument time*, and the ones worth running are those whose predicted effect
clears the achievable precision by a stated margin. That pairing — predicted
separation against achievable precision — is a number both sides can compute and
neither can compute alone. It is also what turns a `BLOCKED` into a proposal
rather than a dead end: *the effect is below your localization precision; either
deepen the DLVO well or change objective.*

**5 · Its bias ledger is what makes the fifth evidence layer usable here.**
[02 §3](02-verification.md#3--result-verification--four-layers-of-evidence) lists
four kinds of evidence and deliberately leaves the fifth — comparison against
experiment — unadopted, because a mismatch there has too many candidate causes
(no HI, polydispersity, tracking error, or simply a different system). Its
lens 6 removes most of them: **G23** carries every bias that damages the
specific quantity being measured, **G24–G26** check that the calibrations behind
it exist. And the terms are already written down over there — a measured MSD
carries `−2D·t_exp/3` from motion blur and `+2ε²` from static localization
error, which **at short lags cancel into a plausible but wrong straight line**.
★ Those terms belong on this side of the comparison, **added to the prediction
rather than subtracted from the data** — the same discipline as
`implementation_check` versus `hypothesis` (rule 7′): correct the model for what
the instrument does, and let the residual mean something. In a domain with no
grader an independent measured oracle is the most valuable evidence there is,
but only when it arrives with its own bias ledger attached.

### What has to be true first

| Precondition | Where it stands |
|---|---|
| This repo **seals its predictions before running** | not yet — §2, and item 1 of §7 below. `simbot` has `SEALED.sha256`; the 8 `bdbot` cases do not go through it |
| The two `knowledge/` schemas here are **unified** | not yet — §2. Adding a third consumer to two divergent schemas would be a mistake |
| The instrument is **connected** on that side | not yet — the working PC and the microscope PC are separate. Its stages 5a–5d are built (2026-08-26) but exercised against a demo config only |
| Illumination power at the sample is **measured** | not yet — its top blocker, deferred by decision (2026-08-19). A power meter, not code |
| τ_c · ℓ_c have **somewhere to live** over there | not yet — its `kb/samples/` arrives with its Phase 4 |
| Computed values have a **provenance kind of their own** | they do not, on either side. That repo has `measured` and `assumed`; a simulated τ_c is neither a measurement of the sample nor a literature value. Giving it its own tier — with this repo's gate verdict as its falsifier — is the honest fix |
| A shared **quantity vocabulary** exists | it does not. Both sides already speak SI with a provenance and a tier, which is the hard half; a common serialization for *"particle diameter, measured, tier 1, ±3 %"* is the missing half |

**The order matters.** Sealing first, here, because an unsealed prediction handed
to an instrument produces an experiment designed around a post-hoc
rationalization. Then the vocabulary, because that is the actual interface and
nothing useful crosses until a number can cross with its provenance intact. The
wiring itself is small once those two exist.

**And one hazard to hold from the start:** a simulated number must never be
allowed to set `evidence: measured`. If it can, the loop closes on itself — this
repo supplies the threshold, the gate clears against it, and the experiment
confirms the simulation that designed it. That is the same failure this project
already has a name for in another form: a checker that was never wired up cannot
be wrong out loud (§2, and the `step_health` case in [05](05-pitfalls.md)). The
rule that keeps that repo's `kb/literature/` honest is the rule this interface
needs.

---

## 7 · What would most improve the science, in order

1. **Wire sealing into `bdbot`.** Everything else on this list is cheaper to fix
   and matters less. Until this lands, the project's central discipline is
   documented but not exercised on its own results.
2. **Finish `network` stage 1**, then stage 2. It is the only route to a real
   `G′(ω)`, and [04](04-cases.md) shows why a single chain cannot give one.
3. **Ensemble `chain-relax-2d-dlvo`.** One seed, 1.8 s per run — this is the
   cheapest unresolved question in the repository.
4. **Get the experiment's real temperature.** A one-number fix that removes a
   −14 % worst case from every timescale.
5. **Unify the two knowledge schemas**, before the count of misfiled lessons
   grows further.
