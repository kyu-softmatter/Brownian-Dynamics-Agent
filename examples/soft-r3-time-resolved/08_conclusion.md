# S8 — conclusion: at `A = 0.1, 1, 10`, what the structure becomes and **when**

run `2026-07-29_soft-r3-time-resolved` · card
[`soft-repulsive-2d--equilibrium-structure`](../../knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md)
· proposed verdict: `PASS_WITH_CAVEATS` · `confirmed_by: null`

---

## Answering S1's question directly

> **"What is the final arrangement at each `A`, and how do `g(r)`, the Voronoi
> structure and the defect structure change over time"**

### ① The final arrangement — none of the three `A` is a crystal

| `A` | `Γ = π^{3/2}A` | global `ψ₆` | defect fraction | coordination kinds | `g(r)` first peak | reading |
|---|---|---|---|---|---|---|
| 0.1 | 0.56 | `0.0477 ± 0.0019` | `0.6445 ± 0.0028` | 6 | `r=0.83, g=1.04` | **almost no structure** |
| 1 | 5.57 | `0.0597 ± 0.0025` | `0.5489 ± 0.0031` | 5 | `r=0.90, g=1.39` | weak liquid structure |
| 10 | 55.68 | `0.2479 ± 0.0080` | `0.2950 ± 0.0068` | **3** | `r=1.03, g=2.85` | **dislocation-rich hexatic-like** |

`ψ₆` falls short of `0.7` in all three cases → **there is no crystal.**
The `ψ₆` difference between `A = 0.1` and `A = 1` is below the finite-size floor
(`1/√N = 0.1`), so it can only be stated as **"no orientational order in either"**.

What separates `A = 10` is not `ψ₆` but **the character of the defects**: the number
of coordination kinds drops `6 → 5 → 3` and the 5-7 imbalance goes
`0.130 → 0.093 → 0.010`, effectively 0. The defects at `A = 10` are not scattered
liquid defects but **5-7 dislocation pairs**.

### ② The change over time — **it is over within `0.03–0.10 τ_d`** ★

This is the run's new result. The single-exponential relaxation of the defect
fraction from the initial placement to the steady state:

| `A` | defects `t=0` | → steady | amplitude | **`τ`** [τ_d] | **`τ` [s]** |
|---|---|---|---|---|---|
| 0.1 | `0.5150` | `0.6475` | `−0.145` | `0.0311 ± 0.0028` | **`71 ± 6 s`** |
| 1 | `0.5150` | `0.5561` | `−0.050` | `0.0420 ± 0.0100` | **`96 ± 23 s`** |
| 10 | `0.5150` | `0.2908` | `+0.174` | `0.0982 ± 0.0063` | **`225 ± 15 s`** |

- **`τ` grows with `A`** — `A=10` vs `A=1` `4.7σ` · vs `A=0.1` `9.7σ`.
  `A=1` vs `A=0.1` is `1.1σ` (**indistinguishable**, even with 16 seeds).
- **The sign of the relaxation amplitude splits.** The `t=0` defect fraction is
  exactly the same at all three `A` (same seed → same placement), and yet `A ≤ 1`
  **increases** while `A = 10` **drops by 41 %.**
  ⇒ The initial placement sits **between** the `A=1` and `A=10` steady states.
- `g(r)` **does not change** from the first time window (`0.2–20 τ_d`) to the last.
  `A=0.1`'s `g(0.5 d)` goes `0.739 → 0.746`. The `min_sep = 0.8 d` shell the initial
  placement enforced had already filled in before the first frame.
- **The last `99.9 %` of the `80 τ_d` is steady state.** In real experimental terms,
  the structure is decided within the first `4 minutes` of `51 hours`.

### ③ Voronoi — the figure shows the character; the time series shows the change

`figs/03–05_voronoi_A*.png` show the coordination-coloured tiling as a time sequence
(red `z=5` · blue `z=7` · grey `z=6`). At `A=10` red and blue appear **adjacent, in
pairs**, and at `A=0.1` several colours are scattered — the fraction alone cannot
show that the same "defect fraction" number is different physics.

> ⚠ **A difference in defect count between panels must not be read as a change in
> time.** The per-frame defect count at `A=10` fluctuates as `30 ± 6` in steady
> state. The panels' `50 → 25 → 34` is mostly that fluctuation. The basis for a change
> in time is the ensemble-mean time series (`figs/01`, `figs/06`).

---

## Assumption dependence — separated into parametric form

| assumption | value | does the conclusion depend on it | basis |
|---|---|---|---|
| **coverage** (`d/σ`) | `8.7266 %` (`d/σ=3`) | **No. Sensitivity exactly 0** | no hard core + `n* = 1` ⇒ the dimensionless physics is a function of `A` alone |
| `σ` | `5 µm` | **only** through `τ_d` in seconds (`τ_d ∝ d²σ`) | `units.scales_soft2d` |
| `T`, `η` | `298.15 K`, `0.890 mPa·s` | only through `τ_d` in seconds | IAPWS, `concepts/water-298k.md` |
| initial condition `min_sep = 0.8 d` | rejection sampling | **Not for the steady state. Yes for the amplitude and `τ`** | the sign reversal in ② above is the evidence |
| `Δt` | the force/thermal displacement gates | unchecked (no `dt/2` convergence) | §open items |
| **`N = 100`** | from the drawing | **Yes at `A=10`** — `βU(r_cut) = 0.09 kT` | card §9. `N ≥ 252` required |

**The `N`-convergence check is this conclusion's largest open item.** `A = 10` still
carries a truncation error and is simultaneously closest to Zahn's transition point
(`Γ = 55.68`, `−7 %` from the boundary) — two weaknesses overlapping at the same
point.

---

## Confidence

| claim | confidence | why |
|---|---|---|
| none of the three `A` is a crystal | **high** | `ψ₆ < 0.25`, 6-fold modulation `< 0.03`. Agrees with the 40 prior runs to within `1.2σ` on independent seeds |
| the defects at `A=10` are 5-7 dislocation pairs | **high** | imbalance `0.010`, 3 coordination kinds. The aggregate histogram matches the prior run |
| `τ` grows with `A` | **medium-high** | `4.7σ` and `9.7σ`. But `A=1` vs `A=0.1` is unresolved and `A=10` is not `N`-converged |
| the sign of the relaxation amplitude splits | **high** | amplitude/noise `4–12×`, starting from the same initial placement |
| the 5 µm discs overlap at `A=0.1` | **high** | minimum separation `0.619 σ`. Correct behaviour, since the model has no hard core |
| agreement with Zahn's phase boundary | **low** | `reproduced: no` → `[source, not reproduced]`. Our reading used observables alone |

## What cannot go in the conclusion

- **"It reached equilibrium"** — there is no threshold. The drift was measured and
  reported, nothing more.
- **Initial-condition independence** — this run used `random` only (`hex` in a square
  box creates artificial defects).
- **`N` convergence** — not run.
- **`Δt` convergence** — not run.

## Next experiments (cheapest first)

1. **A coverage `3.7 %` control** (`d/σ = 4.6`) — removes the `A=0.1` overlap. The
   dimensionless result has to come out unchanged, and **that is the verification of
   the "coverage sensitivity 0" claim.** Cost ≈ this run.
2. **`τ_relax` with 64 seeds** — resolves the `1.1σ` between `A=1` and `A=0.1`. The
   transient pass is 12 runs in `4 s`, so 192 runs is `~1 min`. **The cheapest open
   item.**
3. **`N = 256` convergence** — brings `A=10`'s `βU(r_cut)` down from `0.09` to
   `0.02 kT`. `A=10` costs `65 s/run`, so `N=256` is `~170 s/run` — about `12 min`
   with 4 seeds.
4. `g₆(r)`'s exponent `η₆` — the literature criterion for a hexatic verdict. Zahn's
   reproduction condition §6-3.

## knowledge commits

| item | kind |
|---|---|
| [[coarse-sampling-hides-the-whole-transient]] | finding — `stride ≲ τ_relax/5` |
| [[fraction-threshold-flips-meaning-between-per-frame-and-aggregate]] | finding — a per-frame `0.5 %` threshold is toothless at `N=100` |
| [[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]] | finding — the `t(3)` correction is `1.95×` |
| [[provenance-must-have-one-definition-and-three-capture-points]] | finding — one definition · three capture points (seal, trajectory, analysis) |
| card §8.3 + benchmarks S9–S12 | systems card update |
| `systems/_index.md` | 2 missing cards registered (`9/20`) |

## Reproducibility — what can be said about this run

| check | result |
|---|---|
| bit-for-bit reproduction | ✅ re-ran `A1_s6` → `traj`, `energy`, `max_force`, `init_pos` byte-identical |
| seal integrity | ✅ 3 files under `shasum -a 256 -c` (verified independently of our code) |
| prediction derivation regenerated | ✅ `soft2d_time_series_predict.py`'s output is byte-identical to the sealed copy |
| analysis reproduction | ✅ `--analyze-only` → the same `metrics.json`, PASS 12 / FAIL 5 |
| **git commit** | ❌ **`git_dirty: true`** — `git_rev: b8a3a04` alone does not pin the code |
| the trajectory manifest's `freud` version | ❌ **absent** — the `env` record applies only to runs **after** this fix |

⚠ Those last two rows are this run's reproducibility limit. The 60 trajectories
reproduce bit-for-bit, but **the code that made them is not committed.**
