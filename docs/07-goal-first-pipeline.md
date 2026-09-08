# 07 · Goal-first pipeline — six additions to S1→S8

> **Status: 2026-09-02 — five of the six contracts are built.** Every measurement
> cited in §2 was taken against the tree at `4b4503a`.
>
> | § | Contract | Module | State |
> |---|---|---|---|
> | 3.1 | A · goal | `bdbot/goal.py` | **built** · 26 tests |
> | 3.2 | C · group registry | `bdbot/groups.py` | **built** · 8 groups, `St`/`Pe` collisions resolved |
> | 3.3 | D · abandoned-groups ledger | `bdbot/design.py` | **built** · coverage is a hard gate |
> | 3.4 | F · cost gate | `bdbot/cost.py` | **built** · budget required, `INFEASIBLE` must price the trade |
> | 3.5 | G · run card + seal | `bdbot/runcard.py` | **built** · seal verified **inside** `run.execute()` |
> | 3.6 | K · analysis plan | `bdbot/analysisplan.py` | **built** · estimators resolved at seal time |
> | — | J · monitor (FDT, MSD shape) | `bdbot/health.py` §4 | **built** · warn-only, never kills |
> | — | GSER estimator | `bdbot/microrheo.py` | **built** · recovers `w*/(3 pi)` from analytic input |
> | — | H · smoke | `bdbot/smoke.py` | **built** · 8 of 8 cases, one shared profile (3 unverified) |
> | — | K · plan runner | `bdbot/planrun.py` | **built** · figures come from the plan, not hand-coded |
> | — | **rule 10** · parameter manifest | `bdbot/params.py` | **built** · blocks the run card and `execute()` |
> | — | reference integrator | `bdbot/refbd.py` | **built** · force-free only; refuses forces by design |
>
> **1028 pytest passed** (1 pre-existing environment failure), 61/61 adversarial,
> `A4` intact. A→K has been driven end to end once, on `bead-water-3d` via
> `refbd`, and the sealed prediction **caught a real estimator bug** — the
> one-sided `alpha > 1` mask, which manufactured a 1.4 % storage modulus in a
> fluid with none.
>
> What is **not** built: 3 of the 8 smoke profiles have never executed
> (`verified: False` — they need HOOMD), and no *case script* calls these modules
> yet; the eight archived cases still run the old path.

**What this is not.** It is not a third pipeline. The design that motivated it
was described as eleven agents, A→K; five of those eleven already exist in the
core (four in `bdbot/`, one in `simbot/`), four are half-built, and two are
new. Building A→K as a separate
orchestrator would re-implement working code and turn the
[two-engines seam](00-merge-decisions.md#5--known-seams) into a three-engine
seam. So this document keeps S1→S8 and names **six additions** to it.

**What it is.** One ordering change and six contracts. The ordering change is
that the pipeline becomes **goal-first** rather than **sketch-first**: what the
experimentalist wants to know is fixed before the physical system is read, and
it is carried as an artifact to the end rather than as conversational memory.

---

## 1 · The crosswalk

The eleven-block design mapped onto the existing stages. `→` marks where a
proposed block lands.

| Block | Existing home | Module | Verdict |
|---|---|---|---|
| **A** goal | → **new, before S1** | — | **new artifact** (§3.1) |
| **B** extract physics | S1 + S3 | `intake.py` · `physical.py` · `provenance.py` | exists |
| **C** dimensionless | → S4 | `nondim.py` | **partial** — machinery yes, vocabulary no (§3.2) |
| **D** match design | → S4 | `nondim.py` · `scales.py` | **partial** — abandoned-*groups* ledger is new (§3.3) |
| **E** interactions | S3 | `interactions.py` · `traps.py` | exists |
| **F** cost | → **new gate, before S5** | `dt.py` · `cli.py cmd_calibrate` | **partial** (§3.4) |
| **G** run card + seal | → **new gate, before S5** | `simbot/io.py` — wired to `cli.py`, **not to `bdbot`** | **new as a gate** (§3.5) |
| **H** smoke | S5 | `--smoke` in **8 of 8** via `smoke.py` | **built** |
| **I** run | S5 | `run.py execute()` | exists |
| **J** monitor | S5 | `run.py StepGuard` · `health.Guard` · `health.judge_series` | mostly exists (§5) |
| **K** analysis | S6 + S7 + S8 | `simbot/viz.py` · `report.py` · `analysis/` | exists, plan is new (§3.6) |

Four exist, five partial, two new — **as measured at `4b4503a`**. The row states are now all `built`; the status table at the top of this file is the current state, and this table is kept as the starting point the design was argued from.

---

## 2 · The measured case for each addition

Recording the evidence first, because four of the six additions look like
tidiness until the numbers are in front of you.

**① ★ A blocker for the empty goal already exists, and two cases ran through it
anyway.** `stated_goals` is a `REQUIRED_TOP` field of
`bdbot.observation/0.1` (`intake.py:40`) and is **empty in 2 of 8 intakes** —
`intake/trap-2d-5um/observation.yaml:51` and
`intake/trap-drag-2d-hex300/observation.yaml:69`. The first is the case
[README](../README.md#how-it-is-being-built) names as the smallest complete
one. Both were run anyway: 4 run directories and 81 respectively.

And the gate was not missing. `.claude/skills/bd-intake/SKILL.md:116` says
*"§2.1 An empty `stated_goals` is also a blocker ⭐️"* and `:131` *"If
`stated_goals: []`, do not push it into `choice` — **ask the user.**"* That
section even tabulates `trap-2d-5um` as *"absent — settled by asking the
user"*, and records that `T_obs = 2000 tau_k` and the `tau_k/10` sample
interval were fixed only after the goal was settled conversationally.

So the problem is **not** that nothing reads the field. Three things do:
`interactions.py:225` folds it into the keyword corpus that drives
`recommend()`; `intake.py:277` counts it; `intake.py:166` errors if the key is
missing. The problem is that **the blocker lives in a skill as prose, the
resolution happened in a conversation, and neither reaches S8.** A `question`
field is set in **0 of 8**, and S8's instruction to *"answer S1's question"* is
satisfied by a model re-reading the intake and writing prose. The judgment that
settled `trap-2d-5um`'s goal exists only in a chat log — which is exactly the
failure [03 §5](03-knowledge-base.md#5--capturing-judgment-out-of-conversation)
calls *the real purpose of this project*.

**② There is no shared vocabulary of dimensionless groups.** Measured across
all 8 case scripts:

| Group | Cases using it |
|---|---:|
| `St` | 6 |
| `well_depth/kT` · `n_beads` · `k_bond_star` · `k*` | 3 each |
| `phi` · `L/d` · `Gamma` · `De_trap` · `A` · 7 others | 2 each |
| **56 further names** (`tau_k/tau_int`, `tau_bond/tau_chain_diff`, `sigma_bond/h_min`, …) | 1 each |

**73 distinct group names across 8 cases, 56 of them used exactly once.**

Péclet appears in **2 of 8** cases, Deborah in **2**, Weissenberg in **0**.
Reynolds exists only in `simbot` — the half that has produced no science —
where `simbot/spec.py:598-609` already computes and gates it
(`Check("stokes_reynolds", "pass" if Re < 1e-2 else "fail", …)`), at the same
`1e-2` threshold §3.2 proposes. Every case invents its groups by hand.

**③ ★ Reused group names are not consistently defined.** Two instances, and
they fail in opposite directions.

*One group, two names.* Péclet's two appearances are the **same physical
group**: `Group("Pe", …, "v d_eq/D̄ = tau_B/tau_v")` at `abp_rod_2d.py:155` and
`Group("Pe_drag", …, "v_x d/D_t")` at `trap_drag_2d.py:282`. Both are
`tau_B/tau_v`. A registry query for `Pe` finds one of the two.

*One name, two groups.* `St` is
`tau_p/tau_B` in `abp_rod_2d.py:169`, `chain_bend_2d.py:470`,
`chain_bend_dlvo_2d.py:397`, `soft_r3_2d.py:174` and `trap_drag_2d.py:312`, but
`tau_p/tau_bond` in `chain_relax_2d_dlvo.py:312` — 5 against 1. Both are
defensible Stokes numbers. They are not
the same quantity, they share a name, and `Group.mismatch()`
(`nondim.py:162`, called from `validate()` at `:268`) cannot catch it because
each is a true ratio of its own declared `num`/`den`.

So the two names with the most reuse — the only ones with enough uses to
support a cross-case comparison — are a collision and a split. Aggregating
either today gives a wrong answer, quietly.

**④ No fluctuation-dissipation check exists.** Zero hits for FDT across
`bdbot/`, `simbot/`, `cases/` and `verify/`.

**⑤ Sealing works, and is unreachable from the engine that runs.**
`simbot/io.py:341 write_seal()` is wired — `cli.py:214` calls it, and it has
produced **9 real `SEALED.sha256` files**. All 9 are under `runs_s1s8/`; **0
are under `runs/`.** The dependency is one-way `simbot → bdbot`
(`simbot/units.py:29`, `guards.py:29`, `nondim.py:138`, `estimators.py:129`),
so `bdbot` cannot call it. No `bdbot` run has ever carried a sealed
prediction — 254 of them.

**⑥ The cost gate is an axiom with no conversation.** `A3` requires an estimate
before running and says *"in a chat session, ask rather than refuse."* There is
an estimator (`cli.py cmd_run`, `cmd_calibrate`) and no structured place where
the user states a budget. The archived `trap-2d-5um` production run took
**3038 s** wall; nothing recorded whether anyone had 51 minutes.

---

## 3 · The six contracts

Each contract states the artifact, its gate, **what the stage may decide**,
**what it may not**, which existing function it calls, and how it fails. The
last two rows are the load-bearing ones: a stage that may decide anything is
not a stage.

### 3.1 · A — `goal.yaml`, and the spine that reaches S8

**Artifact.** `intake/<case>/goal.yaml`, schema `bdbot.goal/0.1`, written
**before** `observation.yaml`. It is the first artifact of the pipeline.

```yaml
schema: bdbot.goal/0.1
asked_by:        "experimentalist | self | paper"
question:        "does a DLVO-only bead chain have measurable bending stiffness?"
answering_quantity:
  symbol:        "K_prime"
  unit:          "kT/d^2"
  why:           "linear response coefficient of the O(y^2) term"
decisive_precision:
  value:         0.0
  tolerance:     null            # null => must be supplied at S2, not here
  basis:         "central forces at natural length give U'(l)=0 exactly"
  what_would_change_my_mind: "K' separable from 0 at >3 sigma"
analysis_implied:               # ← consumed by K (§3.6). NOT the plan itself.
  - "K'(omega) with error bars, zero line drawn"
  - "bow vs drive amplitude"
out_of_scope:
  - "absolute G'(omega) of the suspension"
provenance:
  source:        "conversation 2026-08-14"
  confirmed_by:  null            # human only
```

**Gate.** `question` non-empty · `answering_quantity.symbol` resolvable to a
`metrics.observable` name at S7 · `what_would_change_my_mind` non-empty. That
last field **is** the falsifiability gate S1 currently states in prose.

**May decide:** what the question is; which single quantity answers it; what
would refute it.
**May not decide:** any physical parameter value; the tolerance (that is S2's
job, and putting it here would let the goal author set their own pass mark);
`confirmed_by`.

**Calls:** nothing. This stage is pure judgment and has no arithmetic, which is
why it is the only one of the eleven blocks that is unambiguously an LLM.

**Fails how:** silently, if `analysis_implied` is allowed to become the analysis
plan. It is a *hint* to §3.6, written before anything is known about
feasibility. Keep them separate artifacts or the plan inherits the goal's
optimism.

> **Why goal-first rather than sketch-first.** The sketch is evidence about the
> system; the goal is a constraint on the whole pipeline. Reading the sketch
> first means the question gets fitted to whatever the sketch made easy —
> measured: `trap-2d-5um` and `trap-drag-2d-hex300` both have `stated_goals: []`
> and both were run anyway.

### 3.2 · C — a canonical group registry

**Artifact.** `bdbot/groups.py`, a registry of named groups with formula,
required inputs and an ambiguity guard. Cases keep the right to declare
bespoke groups; they lose the right to *name* one of these something else.

| Symbol | Formula | Requires | A mismatch means | Guard |
|---|---|---|---|---|
| `Re` | `rho_f v d / eta` | `rho_f, v, d, eta` | Stokes drag suspect above `1e-2` | must name which `v`; **adopt `simbot/spec.py:598`, do not rewrite** |
| `St` | `tau_p / tau_gov` | `m, gamma, tau_gov` | inertia not negligible above `1e-2` | ★ **must name `tau_gov`** (§2③) |
| `Pe_adv` | `v d / D_t` | `v, d, D_t` | advection dominates diffusion | ★ **absorbs `Pe` and `Pe_drag`** (§2③); distinct from `Pe_shear` |
| `Pe_shear` | `gammadot d^2 / (4 D_t)` | `gammadot, d, D_t` | shear dominates diffusion | — |
| `Wi` | `lambda gammadot` | `lambda, gammadot` | elastic stretch unrelaxed | **undefined for a Newtonian solvent** |
| `De` | `lambda / t_obs` (`lambda omega` oscillatory) | `lambda`, `omega` or `t_obs` | observation shorter than relaxation | — |
| `phi` | packing fraction | `N, d, box` | crowded | meaningless with no pair interaction |
| `k*` | `k d^2 / kT` | `k, d, kT` | trap dominates thermal motion | — |

**The overdamped/underdamped verdict is `St`, and it is deterministic.** It is
already computed as a *check* — `specs/trap-2d-5um__a5ef4f45d589.json:285`
carries `{"name": "inertia negligible  tau_p/tau_k", "value": 8.139e-04,
"limit": 0.01, "margin": 12.286, "hard": true}`. Note the check is named
`tau_p/tau_k` and `trap_2d_5um.py` declares **no `St` group at all**, so the
registry is introducing the name here, not reusing it. This is not an agent
decision and must not become one. What an LLM may do is *choose the governing
timescale* for the denominator; the comparison is `checks.verdict()`.

**Gate.** Every group named from the registry must recompute from the ledger
within `RATIO_RTOL` (already enforced by `NondimSpec.validate()`) **and** match
the registry's declared `num`/`den` roles. The second half is the new part and
is what catches §2③.

**May decide:** which groups are relevant to this system; which timescale
governs.
**May not decide:** a group's formula; a group's threshold; whether a group
passed.

**Calls:** `nondim.Group`, `Group.recompute`, `checks.Check`, `checks.verdict`.

**Fails how:** by growing. A registry that accepts every group anyone wants
becomes the per-case free-for-all it replaced. Admission rule: a group enters
the registry when **two** cases have used it — the repository's existing
promotion rule, unchanged.

### 3.3 · D — the abandoned-groups ledger

The design step most likely to be silently wrong, because matching all groups
at once is generally over-determined and the failure looks like success.

**Artifact.** Two new fields on `bdbot.nondim/0.1`:

```json
"groups_matched":   [{"symbol": "k*",  "experiment": 6.04e4, "simulation": 6.04e4, "rel_err": 0.0}],
"groups_abandoned": [{"symbol": "Re",  "experiment": 3.1e-7, "simulation": 2.0e-4,
                      "reason": "BD has no fluid inertia; both are << 1e-2 so Stokes drag holds in each",
                      "safe_because": "St and Re both clear their thresholds independently"}]
```

**Gate.** `groups_abandoned[].reason` and `.safe_because` non-empty — the
existing `declare_absent(role, reason)` behaviour (`scales.py:128`, raises on an
empty reason) extended from ledger *roles* to *groups*. Union of matched and
abandoned must cover every group the registry says this system supports; a
group in neither is an error, not a warning.

**May decide:** which groups to preserve; which to abandon; the argument for
why abandoning is safe.
**May not decide:** whether the abandonment *is* safe — that is a check with a
threshold, and a mismatch between the stated `safe_because` and a failing check
is a hard error.

**Calls:** `scales.ScaleLedger.declare_absent`, `groups.registry`.

**Fails how:** the agent matches the subset it can hit and reports "matched"
with the remainder unmentioned. The coverage requirement above is the only
thing that stops it, so it must be a hard gate.

### 3.4 · F — the cost gate, with a budget the user stated

**Artifact.** `runs/<run_id>/cost.json`, written before anything runs.

```json
{"schema": "bdbot.cost/0.1",
 "user_budget_s": 3600, "urgency": "same-day",
 "estimate_s": 3038, "estimate_basis": "measured 332 steps/s on this machine, 2026-09-02",
 "margin": 1.18,
 "verdict": "FEASIBLE",
 "if_infeasible_give_up": ["N 1000 -> 400 (statistical error x1.6)",
                           "T_obs 2000 tau_k -> 500 (tau fit degrades)"]}
```

**Gate.** `A3`: no run without an estimate. **New:** no run without
`user_budget_s`, and the source of that number is the user in conversation, not
a default. `INFEASIBLE` is a verdict, not a refusal — it must be accompanied by
`if_infeasible_give_up`, priced in the currency that matters (statistical error,
not step count).

**May decide:** nothing about physics.
**May not decide:** to shrink the system on its own authority. It proposes the
trade; the user takes it. Silently reducing `N` to fit a budget changes the
answer's error bar, which by `A2` changes whether the answer exists.

**Calls:** `dt.compare_criteria` (for the `dt` that actually binds),
`cli.py cmd_calibrate` (throughput), `run.py` step-rate history.

**Fails how:** by estimating from a stale throughput number. `cmd_calibrate`
already refuses to overwrite a global constant because the kernel differs per
system — keep that refusal and store the measurement per case.

### 3.5 · G — the run card, and sealing finally wired to `bdbot`

The gate the current pipeline does not have, and the cheapest place to close
[roadmap item 1](06-roadmap.md#7--what-would-most-improve-the-science-in-order).

**Artifact.** `runs/<run_id>/RUNCARD.md` — human-readable, one screen — plus
`prediction.yaml` and `SEALED.sha256` **in the same directory**, sealed at the
same moment.

The run card states, in this order: the question from §3.1 · the system in SI
with tiers · the groups matched and abandoned from §3.3 · `dt` and which
constraint bound it · the cost verdict from §3.4 · **the prediction, with each
item's role and tolerance** · **the analysis plan from §3.6** · what this design
cannot decide.

**Gate.** Human approval, recorded. Then `simbot.io` seals `prediction.yaml`
and `analysis_plan.yaml` together. `bdbot.run.execute()` gains one line beside
its existing `verify_hash()` call: if a seal exists and does not verify, raise.
`.claude/settings.json` already denies `Edit` on `SEALED.sha256` and
`02_prediction.md`; extend the deny list to the new paths.

**May decide:** how to present the card.
**May not decide:** to proceed without approval; to edit a sealed document; to
seal a prediction whose role is missing.

**Calls:** `simbot.io.write_seal(rundir, stages)` (`simbot/io.py:341`),
`simbot.io.verify_seal` (`:376`), `simbot.spec.load_prediction`
(`simbot/spec.py:832`), `bdbot.nondim.LoadedSpec.verify_hash`
(`bdbot/nondim.py:558`).

⚠️ **This is the one addition that requires `bdbot` to import `simbot`,
reversing the current dependency direction.** `simbot` already imports `bdbot`
in four modules, so a naive import creates a cycle. Either lift `write_seal` /
`verify_seal` into `bdbot/` (they depend on nothing but `hashlib` and paths), or
have `execute()` verify the seal by reading `SEALED.sha256` directly — it is
plain `sha256sum` format precisely so that no code is needed to check it. The
second is preferable and is what §3.5's failure mode (b) argues for anyway.

**Fails how:** two ways, both already seen in this repository. **(a)** The gate
is written and reachable only by hand. This is `health.gate()`'s exact
situation **today**: it is called from `tools/health.py:138` — the path behind
`bdbot.cli health --gate` — and `bdbot/run.py` contains **no reference to
`gate` at all**. So a human who runs the gate sees it, and a run that skips
that command is never gated. `gate()`'s own docstring
(`health.py:524`) records the consequence: *"the reason the two disagreed
unnoticed is that `execute` never called `gate()` — an unwired checker cannot
be wrong out loud."* A checker on an optional path is not enforcement. So the
seal check goes **inside `execute()`**, not beside it and not in a sibling
tool. **(b)** The
seal is verified by the code that wrote it, which verifies nothing; the CI job
uses plain `sha256sum` for exactly this reason and must keep doing so.

### 3.6 · K — the analysis plan, fixed before the data exists

**Artifact.** `runs/<run_id>/analysis_plan.yaml`, derived from §3.1, sealed at
§3.5.

```yaml
schema: bdbot.analysisplan/0.1
answers_goal: intake/chain-bend-2d-dlvo/goal.yaml
planned:
  - id: P1
    figure: "K'(omega), log x, zero line drawn, error bars from block SEM"
    estimator: bdbot.lockin.k_star          # bdbot/lockin.py:52
    decides: "is K' separable from 0"
    role: hypothesis
  - id: P2
    figure: "MSD vs t, analytic 2*dim*D*t overlaid"
    estimator: bdbot.stats.block_sem        # bdbot/stats.py
    decides: "implementation check on D"
    role: implementation_check
not_planned:
  - figure: "g(r)"
    reason: "no pair interaction between non-bonded beads; RDF is uninformative here"
post_hoc: []      # ← appended AFTER the run, and rendered as post-hoc
```

**Gate.** Every `planned[]` item must name what it `decides`, carry a role, and
name an **importable** estimator — resolved at seal time, not at plot time.
Every item in `goal.yaml:analysis_implied` must appear in `planned` **or** in
`not_planned` with a reason — the `FigureSet.skip()` discipline
(`simbot/viz.py:99`) raised to the plan level. Captions remain mandatory at
creation time (`simbot/viz.py:84` rejects an empty one), unchanged.

> The importability requirement is not pedantry: the first draft of this very
> document named `simbot.analysis.trap.lockin_k_star` as P1's estimator. No such
> function exists — the lock-in estimator is `bdbot/lockin.py:52`. A plan that
> names an unresolvable estimator is a plan nobody can execute, and it looks
> exactly like one that can.

**May decide:** the estimator, the axes, the overlays; and after the run, any
number of `post_hoc` figures.
**May not decide:** to move an item from `post_hoc` into `planned`; to drop a
`planned` item because it looked bad.

**Calls:** `simbot.viz.FigureSet`, `simbot.report.render`, `bdbot.metrics.judge`.

> **Why the plan is sealed and not chosen afterwards.** "Decide the best plot
> for this system from the goal" is principled. "Decide the best plot after
> seeing the trajectory" is a garden of forking paths, and with 145 runs in one
> case and 3,856 in one campaign there are enough branches for anything to look
> significant. Sealing the plan costs one file. `post_hoc` is not forbidden —
> it is labelled, which is the whole difference.

---

## 4 · Which of the eleven blocks needs a model

| Block | Needs judgment | Because |
|---|---|---|
| **A** goal | ★ **yes** | no arithmetic exists here at all |
| **B** extract physics | ★ **yes** + retrieval | reading a source, and filling gaps from the KB |
| **C** dimensionless | **no** | formulas and thresholds are fixed; only the choice of `tau_gov` is judgment |
| **D** match design | ★ **partly** | *what* to abandon is judgment; *whether* it is safe is a check |
| **E** interactions | ★ **partly** | the form is judgment (`interactions.recommend()` already scores it); the values are formulas |
| **F** cost | **no** + conversation | the estimate is a model; the budget comes from the user |
| **G** run card | **no** | assembly and a human gate |
| **H** smoke | **no** + triage | running it is mechanical; reading a failure is not |
| **I** run | **no** | — |
| **J** monitor | **no** + triage | the checks are deterministic; the kill decision is narrow (§5) |
| **K** analysis | ★ **yes** | choosing the estimator and reading the result |

Four unambiguous, three partial, four not. Wrapping the four "not" blocks in
subagents costs latency and money and, per
[deterministic-core](../.claude/rules/deterministic-core.md), adds a candidate
cause the next time a result is wrong. They ship as **tools the orchestrator
calls**, with a model narrating the choice.

**On concurrency:** A→B→C→D→E→F→G→H→I→K is a chain, not a wide DAG. Real
fan-out exists only inside **E** (trap, excluded volume and bonds are mutually
independent) and inside **K** (independent estimators over one finished
trajectory).

⚠️ **J is *not* concurrent with I, and must not be made so.** Both monitors run
**synchronously inside the integration loop**: `StepGuard.check` is called from
`execute()`'s `_loop` (`run.py:305-360`), and `health.Guard.as_action`
(`health.py:265`) attaches as a `hoomd.write.CustomWriter` on a `Periodic`
trigger — a blocking callback that calls `state.get_snapshot()`. They are
interleaved, not parallel, and `run.py:335-341` deliberately forces the guard
period equal to the sample period. That is the correct design: a monitor that
observes asynchronously cannot stop a step, and stopping the step is the whole
point.

Wall time is dominated by I regardless — 3038 s against seconds for the entire
specify-and-gate path — so parallelising the agent layer is the wrong target.
Do not design for it.

---

## 5 · The monitor: what it may and may not kill

Kill authority is narrow on purpose.

Two classes called `Guard` own different halves of this, and `run.py:160-173`
documents the split deliberately. Getting the owner wrong means adding a check
to the class that cannot see the quantity.

| Class | Signal | Action | Owner |
|---|---|---|---|
| unambiguous | non-finite **position** | **KILL** | `health.Guard` (`health.py:245`) |
| unambiguous | non-finite **PE/N or force array** | **KILL** | `run.StepGuard` (`run.py:188-205`) |
| unambiguous | `\|r\|max > 50 L` (box escape) | **KILL** | `health.Guard` (`health.py:249`) |
| unambiguous | `\|PE\| > 1e3 \|PE_0\|` (relative) | **KILL** | `health.Guard` (`health.py:257`) |
| unambiguous | `\|PE\|/N > 1e8 kT` (absolute) | **KILL** | `run.StepGuard` (`run.py:194`) |
| unambiguous | `dt*·\|F\|max > 0.1 sigma` | **KILL** | `run.StepGuard` (`run.py:208`) |
| statistical | FDT residual | **WARN + checkpoint** | ✗ new |
| statistical | MSD exponent vs *expected* shape | **WARN + checkpoint** | partly, `judge_series` |
| statistical | PE drift beyond `3 sigma` block SEM | **WARN** | `run.judge` (`run.py:62,263`) |

`health.Guard` sees positions; `run.StepGuard` sees the force arrays and
`PE/N`. An FDT check needs **both** displacements and forces, so it belongs
with neither as written — that is a real design question §3 does not settle,
and the honest answer is that it wants a third observer or a widening of one of
these two.

Two warnings about the two checks that were specifically requested.

**⚠️ FDT cannot tell you anything about the physics.** HOOMD's `Brownian`
integrator imposes `D = kT/gamma` by construction, so a fluctuation-dissipation
check is an `implementation_check` in the sense of rule 7′ — a mismatch is a
bug, never a result. It is still worth having: it is the class of check that
would have caught the missing minimum image (+1856 %). But it must be labelled
`implementation_check`, or it will be read as evidence about the system.

**⚠️ "Is the MSD diffusive" false-positives, and this repository has already
paid for it.** In a trap the MSD *must* plateau; at short times it is steeper
than diffusive; on log-spaced lags an index-based exponent exceeds any fixed
bound. `health.judge_series` grew its `cumulative` flag and its real `t` axis
because MSD was flagged `NUM_DIVERGE` for legitimately diffusing — *"1010×
larger in the second half."* So J cannot run a generic shape test. **The
expected MSD shape must be handed down from D/E as part of the design**, which
makes it another cross-stage contract, and the check compares against that.

**Statistical checks never kill.** They have real variance early in a run, and
the precedent is the gate that refused 80 of 83 specs with zero real hard
failures among them. A jittery kill on a 51-minute job is expensive; a
checkpoint plus a warning is not.

---

## 6 · What this deliberately does not do

- **It does not rename S1→S8** — but the usual reason given for that is wrong,
  so here is the real one. The 227 `runs/*/record.json` post-mortems contain
  **zero** `S1`–`S8` tokens; their schema has no stage field, and stage numbers
  live in `runs_s1s8/` (24 markdown files) instead. So renumbering would *not*
  invalidate the post-mortems. What it would invalidate is the agent layer:
  6 skills, 5 `bd-pipeline` reference files named `s1_…` through `s8_…`,
  9 subagent descriptions, and `docs/01`–`docs/06`. That is a large edit with no
  scientific return, which is the actual argument. The goal-first ordering is
  achieved by adding A *before* S1.
- **It does not merge the two engines.** `bdbot` still runs, `simbot` still
  seals; §3.5 wires the seal into `bdbot.run.execute()` and leaves the runners
  separate. Merging the pipelines is a human decision, not a patch
  ([00 §5](00-merge-decisions.md#5--known-seams)).
- **It does not unify the two knowledge schemas.** Still open, still
  [roadmap item 5](06-roadmap.md#7--what-would-most-improve-the-science-in-order).
- **It does not add rods or polymers.** §3.2's registry is written so that
  anisotropic friction and bond/angle groups can enter under the
  appeared-twice rule, and nothing more.
- **It adds no LLM dependency to the core.** `A4` unchanged;
  `tests/test_invariants.py` must keep passing.

---

## 7 · Order of implementation

Ordered by *evidence gained per unit of work*, not by pipeline position.

0. ~~**Promote `bd-intake` §2.1 from prose to a validator.**~~ **DONE
   2026-09-02.** `Observation.open_goal` + `Observation.blockers` in
   `bdbot/intake.py`, wired into `cli.py` (4 sites), 13 pytest cases in
   `tests/test_s1_intake_goal.py`, 18 adversarial checks added to
   `verify/verify_intake_guards.py`. `bdbot.cli status` now reports
   `trap-2d-5um` and `trap-drag-2d-hex300` as **BLOCKED · stated_goals**.

   ⚠️ **Implementing it proved this item's own prescription wrong twice**, which
   is worth keeping rather than editing away:

   - It said *"make `intake.validate()` return an `error` Issue."* Wrong level.
     `errors` renders as `VERDICT: FAIL -- schema error`, but the key **is**
     present; only the content is missing. That is `BLOCKED` semantics, so it
     had to be a property alongside `open_missing`, not an Issue.
   - More seriously: `ready_for_system` had **zero call sites**. All four real
     consumers re-derived readiness from `open_missing` alone. Adding the
     blocker to the property, exactly as this item prescribed, would have
     changed **nothing observable** — it would have been a seventh unwired
     checker, written by the document complaining about unwired checkers.
     Fixed by making `Observation.blockers` the single source and defining
     `ready_for_system` as `not blockers`.

   Both are filed in `knowledge/entries/` with `origin: tooling`.
1. **§3.1 `goal.yaml` + the S8 consumption gate.** Smallest artifact, and it is
   the spine the other five hang from. Backfill the 8 existing intakes; the two
   with `stated_goals: []` are the test of whether the schema is honest.
2. **§3.5 sealing wired into `execute()`.** One line inside `execute()` plus a
   deny-list entry. Closes the repository's top acknowledged hole and needs no
   new physics.
3. **§3.2 the group registry**, seeded with the groups that already appear
   twice, and `St` split into `St(tau_gov)` with the denominator named. Fixes
   §2③.
4. **§3.6 the analysis plan.** Depends on 1.
5. **§3.3 the abandoned-groups ledger.** Depends on 3.
6. **§3.4 the cost gate.** Cheapest to state, least likely to be wrong, and the
   only one whose absence has never yet produced a wrong result.

**Do one end to end before abstracting.** `trap-2d-5um` is the right vehicle —
it is the smallest complete case, it has an empty `stated_goals`, and its five
observables are all `implementation_check` against closed-form answers, so a
mistake in the new plumbing shows up as a number rather than as a plausible
figure.

---

## 8 · Registered risks

| Risk | Why it is real here | Mitigation |
|---|---|---|
| **The goal becomes decoration** | ★ it already has: an explicit blocker exists in `bd-intake` §2.1 and 2 of 8 cases ran past it | S8 cannot render `REPORT.md` without resolving `answering_quantity` to a measured observable — a *code* gate, since the prose one demonstrably did not hold |
| **A prose gate is mistaken for a gate** | the empty-goal blocker (`bd-intake` §2.1) and the original `A4` grep — both written, neither enforced | any rule that matters becomes a validator or a test; §7 item 0 is the first repayment |
| **The registry becomes the free-for-all it replaced** | 73 bespoke group names today, 56 used once | admission only at the second independent use |
| **A gate lands on an optional path** | ★ `health.gate()` **today**: called from `tools/health.py:138`, absent from `run.py` entirely | every gate in §3 lands *inside* the function that would otherwise proceed, never in a sibling tool, and arrives with a test that deliberately breaks it |
| **The analysis plan is written after a peek** | 145 runs in one case, 3,856 in one campaign | plan sealed with the prediction at §3.5; `post_hoc` is a separate, rendered-as-post-hoc list |
| **`decisive_precision` gets set to whatever the design can reach** | this is the design-power question, and the repository has answered it well exactly once (22.3 sigma / 0.94×) | the field is written at A, *before* F and G know what is affordable, and is immutable afterwards |
| **Six additions is still a rewrite** | 4 of 6 touch `nondim.py`, the single L2↔L4 contract | additions are new *fields* on `bdbot.nondim/0.1`; `physics_only()` must keep them out of the hash, or every one of the 278 specs re-ids |

★ That last row is the one to check first. `run_id` hashes
`{system, params, numerics}` with `physics_only()` applied
(`runid.py:41-48`, stripping the 19 `DOC_KEYS` at `:33-38`). `groups_matched`
and `groups_abandoned` (§3.3) are **documentation of a design decision, not
physics inputs**, so they belong in `DOC_KEYS`. Get that wrong and all 278
specs change identity at once — which has already happened once, to every
`chain-bend-2d-dlvo` run simultaneously
([04 §aggregate-by-tag](04-cases.md)).

---

## 9 · Provenance of this document

Drafted 2026-09-02 in conversation, then fact-checked against the tree at
`4b4503a`. **Six claims in the first draft were wrong** and are recorded here
rather than silently corrected, because the same failure modes are what §3
gates against:

| First draft said | Actually | Effect on the argument |
|---|---|---|
| `--smoke` in all 8 case scripts | **5 of 8** (`chain_bend_2d`, `chain_bend_dlvo_2d`, `network_3d` lack it) | block H is partial, not existing |
| `stated_goals` is "read by nothing" | **3 consumers**, and a prose blocker in `bd-intake` §2.1 | argument got **stronger** — a gate existed and was walked past |
| Péclet in 1 of 8 cases | **2 of 8**, as `Pe` and `Pe_drag` — one group, two names | argument got **stronger** — §2③ is a collision *and* a split |
| 227 post-mortems cite stage numbers | **zero** do; stage numbers are in `runs_s1s8/` | §6's stated reason not to renumber was false; replaced |
| J is concurrent with I | both monitors are **synchronous, blocking** | reversed into a design constraint |
| `simbot.io.seal`, `simbot.analysis.trap.lockin_k_star` | `write_seal`, `bdbot.lockin.k_star` | two of three named functions did not exist |
| `health.gate()` is never called | called from `tools/health.py:138`; absent from `run.py` | §3.5(a) narrowed from "unwired" to "on an optional path" |

⚠️ That last row was introduced **by the fact-check itself**, which reported
zero call sites for `health.gate()` and missed `tools/health.py:138`. It was
caught by a third pass. One reviewer is not a verification layer — which is
`A1`, and the reason it asks for three pieces of evidence of different kinds.

Three of those six were overstatements in the direction of making the proposal
look more necessary, and two of the three turned out to understate the problem
instead. That asymmetry is worth noticing: an unverified draft is not randomly
wrong, it is wrong in the direction its author was arguing.

---

**See also:** [01 Architecture](01-architecture.md) ·
[02 Verification](02-verification.md) ·
[00 Merge decisions §5](00-merge-decisions.md#5--known-seams) ·
[06 Roadmap §7](06-roadmap.md#7--what-would-most-improve-the-science-in-order) ·
[axioms](../.claude/rules/axioms.md)
