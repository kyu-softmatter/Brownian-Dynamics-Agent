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
sensitive, so this propagates −4 % (at 298 K) to −14 % (at 293 K) into every
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
built by the same person for the same systems:

| | [**Brownian-Dynamics Agent**](https://github.com/kyu-softmatter/Brownian-Dynamics-Agent) (this repo) | [**agentic-microscope**](https://github.com/kyu-softmatter/agentic-microscope) |
|---|---|---|
| Input | a sketch of a physical system | a research goal |
| Decides | what the system does, in silico | what the instrument can actually record |
| Refuses when | a parameter has no provenance | a measurement has no calibrated input |
| Produces | a dimensionless spec + a defended result | executable microscope settings + evidence tier |
| Its knowledge base | system cards, findings, benchmarks | instrument config, calibrations, tacit expertise |

They already share their architecture — a deterministic core under an agent
layer, hard gates that return `BLOCKED` naming the one missing input, and a
knowledge base that is read before every decision and written after every
verdict. That is not a coincidence; the second was built from the first's
lessons. **Neither is finished, and joining them before each stands alone would
couple two moving targets.** So this is future work, not a current task.

### Why joining them is worth doing

**1 · A simulation can tell the microscope what resolution the question needs.**
Right now the microscope agent's sample and detection lenses ask *"is this
measurable?"* against thresholds a human supplies. But whether 40 nm of
localization precision is enough is a property of **the physics being measured**,
not of the microscope. This repository computes exactly that: `chain-bend-2d-dlvo`
established that bow discriminates DLVO from JKR at 22.3σ under a soft trap and
at 1.4× under a stiff one — which is a *statement about what an experiment must
do to see the difference*. Fed forward, it becomes a required precision, a
required frame rate and a required trap stiffness, derived rather than guessed.

**2 · A measurement can close this repository's open assumptions.**
The single most damaging soft spot here is `T = 300 K`, labelled tier 1 but
actually inherited from a sketch with no temperature — worth −4 % to −14 % on
every timescale (§2). A thermometer reading closes it. The same is true of the
particle-size distribution, the salt concentration and the surface potential:
each is currently a tier-1 *choice*, and each is something the microscope side
either measures or has already measured. **Provenance is the shared currency** —
both repositories already carry `tier` and a falsifier on every stored value, so
a measured number can replace an assumed one without either side losing track of
which is which.

**3 · The two refusals compose into one honest answer.**
Today, asking "can I see this effect?" needs two separate consultations, and the
answers can silently contradict: the simulation says the effect is 0.1 d in
amplitude, the microscope says it can resolve 0.05 d, and nobody checks that the
two `d` mean the same thing in the same units. Joined, a single `BLOCKED` would
name the single missing input — *"the effect is below your localization
precision; either raise the DLVO well depth or change objective"* — with the
physics and the optics reconciled in one place.

**4 · Real trajectories become a fifth evidence layer.**
[02 §3](02-verification.md#3--result-verification--four-layers-of-evidence) lists
four kinds of evidence, and `BD_agent` reserved a fifth — comparison against
experiment — but deliberately did not adopt it, because a mismatch has too many
candidate causes (no HI, polydispersity, tracking error, or simply a different
system). **A microscope agent that reports its own calibration and evidence tier
removes most of those candidates**, which is precisely what makes the fifth layer
usable. In a domain with no grader, an independent measured oracle is the most
valuable evidence there is.

### What has to be true first

| Precondition | Where it stands |
|---|---|
| This repo's predictions are **sealed before running** | not yet — §2, and it is item 1 below |
| The microscope agent is **connected to the instrument** | not yet — it produces offline recommendations; the working PC and the microscope PC are separate |
| Illumination power at the sample is **measured** | not yet — its top blocker, and it needs a power meter, not code |
| The two `knowledge/` schemas here are **unified** | not yet — §2. Adding a third consumer to two divergent schemas would be a mistake |
| A shared **quantity vocabulary** exists | does not exist. Both sides speak SI with provenance and tiers, which is the hard half; a common serialization for "particle diameter, measured, tier 1, ±3 %" is the missing half |

**The order matters.** Sealing first, because an unsealed prediction fed to an
instrument produces an experiment designed around a post-hoc rationalization.
Then the shared vocabulary, because that is the actual interface. The connection
itself is small once those two exist.

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
