"""Proposing a cutoff radius (`r_cut`) — derived from a **stated error tolerance**,
not from convention.

## Why this is its own module

`r_cut = 2.5 σ` is a convention from the LJ-liquid literature and has no basis for
our systems. The knowledge base shows what this actually costs:

- `knowledge/wiki/benchmarks/choi2020-interfacial-rdf.md` — *"the cutoff distance is
  undetermined. The original paper does not mention it. **First suspect if
  reproduction fails.**"*
- `knowledge/wiki/systems/interfacial-colloid--equilibrium-structure.md` §7 —
  `r_cut/d ≥ 21–69` is required while minimum image forces `r_cut ≤ L/2 = 65 d`.
  **Systems whose box cannot contain the interaction range are real.**

So this module always returns three things together:
  ① the `r_cut` value  ② **which constraint set it**  ③ the energy and force left
  at the cutoff

## Convention

| potential | `r_cut` | character |
|---|---|---|
| WCA | `2^{1/6} σ_LJ` | **a definition, not a choice.** The LJ minimum |
| truncated LJ · Morse · Yukawa · DLVO | numerical solution from the tolerance | a choice — the basis has to be recorded |

Absolute constraint: `r_cut ≤ L/2` (minimum image)
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Callable

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

WCA_RCUT_OVER_SIGMA_LJ = 2.0 ** (1.0 / 6.0)   # 1.1224620483...

# =============================================================================
# Tolerance presets — **the observable decides the tolerance**
# =============================================================================
# The literature convention `r_cut = 2.5 sigma` is a value tuned to structural
# calculations in LJ liquids. Measuring what is left at that point at beta*eps=1
# gives beta U = 1.63e-2 and beta F sigma = 3.90e-2. So the conventional value is a
# choice to say "beta U up to ~2e-2 is acceptable", and that is not sufficient for
# every observable.
#
# So the tolerances are keyed by observable. Choosing a value **records the basis
# along with it.**
TOLERANCE_PRESETS: dict[str, tuple[float, float]] = {
    # name          : (beta_u_tol, beta_f_sigma_tol)
    "convention":    (2.0e-2, 4.0e-2),   # equivalent to the 2.5 sigma convention.
                                         # For structural calculations
    "structure":     (1.0e-2, 2.0e-2),   # g(r), coordination number, S(q)
    "thermodynamics": (1.0e-3, 1.0e-2),  # B2, phase boundaries, osmotic pressure —
                                         # the default
    "precision":     (1.0e-4, 1.0e-3),   # near a critical point, precise free energy
}
DEFAULT_TOLERANCE_PRESET = "thermodynamics"
DEFAULT_BETA_U_TOL, DEFAULT_BETA_F_SIGMA_TOL = TOLERANCE_PRESETS[DEFAULT_TOLERANCE_PRESET]

# The cost-comparison baseline. In a cell list the neighbour count goes as r_cut^dim.
CONVENTIONAL_LJ_RCUT = 2.5


@dataclass
class CutoffProposal:
    """An `r_cut` proposal — the value + the dominant constraint + the residual
    error. The three are never separated."""

    potential: str
    r_cut_over_sigma: float
    criterion: str                     # which constraint set the value
    beta_U_at_cut: float               # potential left at the cutoff [kT]
    beta_F_sigma_at_cut: float         # force left at the cutoff [kT/sigma]
    min_box_over_sigma: float          # the minimum L minimum image requires
    exact: bool                        # is the value fixed by the potential's
                                       # definition
    tolerance_preset: str = ""         # which preset was used
    rationale: str = ""
    warnings: list[str] = field(default_factory=list)
    hoomd_note: str = ""

    def cost_vs_convention(self, dim: int = 3) -> float:
        """Neighbour-pair multiple relative to the conventional `2.5 sigma`
        (= the pair-force cost multiple).

        In a cell list the neighbour count goes as `r_cut^dim`. This value feeds
        straight into `run_policy.yaml`'s wall-time budget, so it is always
        reported alongside.
        """
        return (self.r_cut_over_sigma / CONVENTIONAL_LJ_RCUT) ** dim

    def box_ok(self, box_over_sigma: float) -> tuple[bool, str]:
        """The minimum-image check. `r_cut <= L/2` is required."""
        if box_over_sigma >= self.min_box_over_sigma:
            return True, ""
        return False, (
            f"minimum image violation: r_cut/sigma = {self.r_cut_over_sigma:.4g} "
            f"requires L/sigma >= {self.min_box_over_sigma:.4g} but it is "
            f"{box_over_sigma:.4g}. "
            f"Enlarge the box, or reduce r_cut and record the price."
        )


# =============================================================================
# WCA — not a choice
# =============================================================================
def wca_cutoff() -> CutoffProposal:
    """WCA's `r_cut` is **fixed by definition** at the LJ minimum `2^{1/6} σ_LJ`.

    The force is exactly 0 there, so **the force is continuous** -- there is no
    truncation artefact. The potential is off by `-ε`, so `mode="shift"` has to
    lift it for the interaction to be purely repulsive.

    ⚠ The lab code (Quah `graybox_abp_mpc`) has **no** `mode="shift"` --
      the force is right but a `-ε` offset is left in the energy. That does not
      affect the BD dynamics, but logging the potential energy or computing an
      equation of state gives wrong answers.
      Basis: knowledge/source/papers/2024-quah-graybox-abp-mpc-repo.md
      §"WCA as an LJ truncation"
    """
    return CutoffProposal(
        potential="WCA",
        r_cut_over_sigma=WCA_RCUT_OVER_SIGMA_LJ,
        criterion="potential_definition",
        beta_U_at_cut=0.0,
        beta_F_sigma_at_cut=0.0,
        min_box_over_sigma=2.0 * WCA_RCUT_OVER_SIGMA_LJ,
        exact=True,
        rationale=(
            "WCA = LJ cut at its minimum and lifted. r_cut = 2^{1/6} sigma_LJ is a "
            "definition, and F = 0 there so the force is continuous. There is no "
            "room for an error tolerance to enter."
        ),
        hoomd_note=(
            'md.pair.LJ(nlist, default_r_cut=2**(1/6), mode="shift"); '
            'params = dict(epsilon=eps, sigma=sigma_LJ).  '
            'To use the particle diameter d, set sigma_LJ = d / 2**(1/6) and '
            'r_cut = d.'
        ),
    )


def barker_henderson_diameter(beta_epsilon: float, sigma_lj: float = 1.0) -> float:
    """WCA's effective hard-sphere diameter
    d_BH = int_0^{r_min} [1 - exp(-beta u(r))] dr.

    This correction is required before putting `phi` into a hard-sphere theory such
    as Carnahan-Starling.
    Basis: knowledge/wiki/findings/wca-reproduces-carnahan-starling.md
      — *"do not put phi straight into CS. At beta*eps=1, d_eff = 1.017 sigma"*
    """
    if beta_epsilon <= 0:
        raise ValueError("beta_epsilon must be > 0")
    r_min = WCA_RCUT_OVER_SIGMA_LJ * sigma_lj

    def integrand(r: float) -> float:
        s6 = (sigma_lj / r) ** 6
        bu = beta_epsilon * (4.0 * (s6 * s6 - s6) + 1.0)
        return 1.0 - math.exp(-bu) if bu < 700 else 1.0

    # As r -> 0 the integrand -> 1 and is smooth. Leaving the lower bound at 0
    # overflows exp, so the analytically-1 region is split off.
    r_lo = 0.3 * sigma_lj
    val, _ = quad(integrand, r_lo, r_min, limit=200)
    return r_lo + val


# =============================================================================
# The general routine that solves r_cut from an error tolerance
# =============================================================================
def cutoff_from_tolerance(
    *,
    potential: str,
    beta_u: Callable[[float], float],
    beta_f_sigma: Callable[[float], float],
    r_start: float,
    r_max: float = 100.0,
    preset: str = DEFAULT_TOLERANCE_PRESET,
    beta_u_tol: float | None = None,
    beta_f_sigma_tol: float | None = None,
    hoomd_note: str = "",
) -> CutoffProposal:
    """The smallest `r_cut` (in sigma) satisfying `|βU(r_cut)| ≤ beta_u_tol`
    **and** `|βFσ(r_cut)| ≤ beta_f_sigma_tol`.

    Records which of the two constraints dominated -- without that, why `r_cut`
    took that value cannot be reconstructed later.

    Args:
        preset: a key of `TOLERANCE_PRESETS`. The observable decides the tolerance.
        beta_u_tol / beta_f_sigma_tol: stated explicitly, these override the preset.
    """
    if preset not in TOLERANCE_PRESETS:
        raise KeyError(f"unknown preset {preset!r}; known: {sorted(TOLERANCE_PRESETS)}")
    p_u, p_f = TOLERANCE_PRESETS[preset]
    beta_u_tol = p_u if beta_u_tol is None else beta_u_tol
    beta_f_sigma_tol = p_f if beta_f_sigma_tol is None else beta_f_sigma_tol
    def solve(f: Callable[[float], float], tol: float) -> float | None:
        """Find the **last crossing** such that `|f| <= tol` for all `r > solution`.

        ⚠ Do not find the first crossing. **The force of an attractive potential is
          non-monotonic** -- LJ's force is exactly 0 at its minimum
          `2^{1/6} sigma` and grows again beyond it. Finding the first crossing
          concludes "already negligible at the start" and returns `r_start`. That
          is how it was actually implemented on 2026-07-28, and a test caught it.
          (Morse's force is non-monotonic for the same reason.)
        """
        g = lambda r: abs(f(r)) - tol
        grid = np.geomspace(r_start, r_max, 3000)
        over = np.nonzero(np.array([g(r) for r in grid]) > 0.0)[0]
        if over.size == 0:
            return r_start                      # negligible over the whole range
        last = int(over[-1])
        if last == grid.size - 1:
            return None                         # still outside tolerance at r_max
        return float(brentq(g, grid[last], grid[last + 1], xtol=1e-13, rtol=1e-14))

    r_u = solve(beta_u, beta_u_tol)
    r_f = solve(beta_f_sigma, beta_f_sigma_tol)

    warnings: list[str] = []
    if r_u is None or r_f is None:
        warnings.append(
            f"the tolerance is still not satisfied at r/sigma = {r_max:g} — "
            f"this is a long-range interaction. Ewald, or an explicit truncation "
            f"compromise, is required."
        )
        r_cut = r_max
        criterion = "r_max_reached"
    elif r_u >= r_f:
        r_cut, criterion = r_u, "potential_tolerance"
    else:
        r_cut, criterion = r_f, "force_tolerance"

    return CutoffProposal(
        potential=potential,
        r_cut_over_sigma=r_cut,
        criterion=criterion,
        beta_U_at_cut=abs(beta_u(r_cut)),
        beta_F_sigma_at_cut=abs(beta_f_sigma(r_cut)),
        min_box_over_sigma=2.0 * r_cut,
        exact=False,
        rationale=(
            f"of |beta U| <= {beta_u_tol:g} (solution r/sigma = "
            f"{'n/a' if r_u is None else format(r_u, '.4g')}) and "
            f"|beta F sigma| <= {beta_f_sigma_tol:g} (solution "
            f"{'n/a' if r_f is None else format(r_f, '.4g')}), the "
            f"{'potential' if criterion == 'potential_tolerance' else 'force'} "
            f"constraint dominates."
        ),
        warnings=warnings,
        hoomd_note=hoomd_note,
    )


# =============================================================================
# Per-potential wrappers
# =============================================================================
def lj_cutoff(beta_epsilon: float, **kw) -> CutoffProposal:
    """Truncated LJ with attraction. U = 4eps[(s/r)^12 - (s/r)^6], distance in
    sigma_LJ.

    Reported against the conventional `2.5 sigma` -- whether to use the convention
    is decided from the basis.
    """
    bu = lambda r: beta_epsilon * 4.0 * (r ** -12 - r ** -6)
    bf = lambda r: beta_epsilon * 24.0 * (2.0 * r ** -13 - r ** -7)
    p = cutoff_from_tolerance(
        potential="LJ (truncated)", beta_u=bu, beta_f_sigma=bf,
        r_start=WCA_RCUT_OVER_SIGMA_LJ, r_max=50.0,
        hoomd_note='md.pair.LJ(..., mode="shift") or md.pair.ForceShiftedLJ '
                   '(the latter to remove the force discontinuity too)',
        **kw)
    p.rationale += (
        f"  [at the conventional r_cut = 2.5 sigma, beta U = {abs(bu(2.5)):.3e}, "
        f"beta F sigma = {abs(bf(2.5)):.3e}]"
    )
    return p


def yukawa_cutoff(beta_epsilon: float, kappa_sigma: float, **kw) -> CutoffProposal:
    """Screened Coulomb. U = eps * exp(-kappa r)/r (the HOOMD md.pair.Yukawa
    convention).

    Distance in sigma, `kappa_sigma` = kappa*sigma.
    ★ `r_cut` is governed by the screening length, so a small `kappa_sigma` (weak
      screening) makes the range explode -- that is the problem the interfacial
      colloid card ran into.
    """
    ks = kappa_sigma
    bu = lambda r: beta_epsilon * math.exp(-ks * r) / r
    bf = lambda r: beta_epsilon * math.exp(-ks * r) * (1.0 / r ** 2 + ks / r)
    return cutoff_from_tolerance(
        potential=f"Yukawa (kappa*sigma={ks:g})", beta_u=bu, beta_f_sigma=bf,
        r_start=1.0, r_max=200.0,
        hoomd_note="md.pair.Yukawa(nlist, default_r_cut=r_cut); "
                   'params = dict(epsilon=eps, kappa=kappa).  mode="shift" advised',
        **kw)


def morse_cutoff(beta_D0: float, alpha_sigma: float, r0_over_sigma: float = 1.0,
                 **kw) -> CutoffProposal:
    """Morse attraction. U = D0[exp(-2a(r-r0)) - 2exp(-a(r-r0))] (HOOMD convention).

    The default attraction for aggregating colloids (domain D).
    """
    a, r0 = alpha_sigma, r0_over_sigma
    def bu(r: float) -> float:
        e = math.exp(-a * (r - r0))
        return beta_D0 * (e * e - 2.0 * e)
    def bf(r: float) -> float:
        e = math.exp(-a * (r - r0))
        return beta_D0 * 2.0 * a * (e * e - e)
    return cutoff_from_tolerance(
        potential=f"Morse (alpha*sigma={a:g}, r0/sigma={r0:g})",
        beta_u=bu, beta_f_sigma=bf, r_start=r0, r_max=100.0,
        hoomd_note="md.pair.Morse(nlist, default_r_cut=r_cut); "
                   'params = dict(D0=D0, alpha=alpha, r0=r0).  mode="shift" advised',
        **kw)


# =============================================================================
# ★ The final proposal — look at the convention first, but extend when affordable
# =============================================================================
def neighbors_per_particle(r_cut_over_sigma: float, phi: float, dim: int = 3) -> float:
    """Neighbours falling inside `r_cut`, per particle.

    3D:  n = 6 phi/(pi sigma^3)  =>  N_nb = (4pi/3) n r_cut^3 = 8 phi (r_cut/sigma)^3
    2D:  n = 4 phi/(pi sigma^2)  =>  N_nb = pi n r_cut^2      = 4 phi (r_cut/sigma)^2

    ★ **This is the entire cost of `r_cut`.** In a dilute system (small `phi`) a
      large `r_cut` still catches almost no neighbours, so it is effectively free.
      For a single particle it is 0.
    """
    if not 0.0 <= phi:
        raise ValueError("phi must be >= 0")
    if dim == 3:
        return 8.0 * phi * r_cut_over_sigma ** 3
    if dim == 2:
        return 4.0 * phi * r_cut_over_sigma ** 2
    raise ValueError(f"dim must be 2 or 3, got {dim}")


# The condition Lambda (the throughput constant) was measured under. The cost
# baseline.
#   knowledge/wiki/findings/local-cpu-parallelism.md — WCA, phi=0.30, 3D
LAMBDA_BASELINE_RCUT = WCA_RCUT_OVER_SIGMA_LJ
LAMBDA_BASELINE_PHI = 0.30
DEFAULT_NEIGHBOR_BUDGET = 150.0     # cap on neighbours per particle. About 4x the
                                    # convention (phi=0.3, 2.5s)


def pair_cost_vs_lambda_baseline(r_cut_over_sigma: float, phi: float,
                                 dim: int = 3) -> float:
    """Pair-force cost multiple relative to the
    `estimators.THROUGHPUT_PARTICLE_STEPS_PER_S` baseline.

    ⚠ **The throughput constant was measured on WCA (`r_cut = 1.122 sigma`,
      `phi = 0.30`).** That condition has only 3.4 neighbours -- a very cheap
      regime. Enlarge `r_cut` and the wall-time estimate is optimistic by exactly
      this factor. This multiple has to be applied to `run_policy.yaml`'s budget
      check.
    """
    base = neighbors_per_particle(LAMBDA_BASELINE_RCUT, LAMBDA_BASELINE_PHI, dim)
    here = neighbors_per_particle(r_cut_over_sigma, phi, dim)
    return max(here, 1e-12) / max(base, 1e-12)


def propose_cutoff(
    tolerance_proposal: CutoffProposal,
    *,
    phi: float,
    n_particles: int,
    dim: int = 3,
    box_over_sigma: float | None = None,
    conventional: float = CONVENTIONAL_LJ_RCUT,
    neighbor_budget: float = DEFAULT_NEIGHBOR_BUDGET,
) -> CutoffProposal:
    """★ The final `r_cut` decision. Implements the user's rule (2026-07-28):

    > *"Consider the conventional value first, but for a dilute system (a single
    >   particle or few particles) or a long-range interaction it seems fine to
    >   consider going longer."*

    The implementation:

    1. **If it is fixed by definition, use it** (WCA). Discussion over.
    2. If the tolerance-based `r_cut` is **smaller** than the convention, use it --
       a free gain.
    3. If larger, ask **whether the neighbour count is affordable**:
       `N_nb = 8 phi (r_cut/sigma)^3`.
       - affordable (`<= neighbor_budget`) → **use the tolerance-based value.** The
         accuracy gain is close to free. Dilute systems, single particles and
         long-range interactions all land here.
       - not affordable → reduce to the largest `r_cut` the budget allows and
         **record the residual error as a known limit** (never pass over it
         quietly).
    4. `r_cut <= L/2` is absolute, so it is applied last. A violation leaves a
       warning -- there is no automatic way to resolve that case (interfacial
       colloid card §7).
    """
    p = tolerance_proposal
    if p.exact:
        p.rationale += "  [fixed by definition — no convention or cost argument]"
        return p

    r_tol = p.r_cut_over_sigma
    nb_tol = neighbors_per_particle(r_tol, phi, dim)
    nb_conv = neighbors_per_particle(conventional, phi, dim)
    notes: list[str] = list(p.warnings)

    single = n_particles <= 1
    if single:
        chosen, why = r_tol, "single_particle_free"
        notes.append("there is 1 particle — no pair interaction, so r_cut is "
                     "irrelevant to the dynamics.")
    elif r_tol <= conventional:
        chosen, why = r_tol, "tolerance_below_convention"
    elif nb_tol <= neighbor_budget:
        chosen, why = r_tol, "dilute_or_longrange_affordable"
        notes.append(
            f"larger than the conventional {conventional:g} sigma, but with only "
            f"{nb_tol:.1f} neighbours (budget {neighbor_budget:g}) it is "
            f"affordable — take the accuracy. phi={phi:g}."
        )
    else:
        # invert for the largest r_cut the budget allows
        if dim == 3:
            chosen = (neighbor_budget / (8.0 * phi)) ** (1.0 / 3.0)
        else:
            chosen = (neighbor_budget / (4.0 * phi)) ** 0.5
        chosen = max(chosen, conventional)
        why = "neighbor_budget_limited"
        notes.append(
            f"the tolerance requires r_cut = {r_tol:.3g} sigma but that is "
            f"{nb_tol:.0f} neighbours, over the budget {neighbor_budget:g}. "
            f"Reduced to {chosen:.3g} sigma — "
            f"★ residual error: beta U = {abs(p.beta_U_at_cut):.2e} was the target "
            f"and the actual value is larger. This is a **known limit** and the "
            f"cutoff sensitivity has to be checked in S7."
        )

    out = CutoffProposal(
        potential=p.potential,
        r_cut_over_sigma=chosen,
        criterion=why,
        beta_U_at_cut=p.beta_U_at_cut if chosen >= r_tol else float("nan"),
        beta_F_sigma_at_cut=p.beta_F_sigma_at_cut if chosen >= r_tol else float("nan"),
        min_box_over_sigma=2.0 * chosen,
        exact=False,
        tolerance_preset=p.tolerance_preset,
        rationale=(
            f"{p.rationale}  || final decision: {why}. "
            f"{neighbors_per_particle(chosen, phi, dim):.1f} neighbours "
            f"(at the conventional {conventional:g} sigma it is {nb_conv:.1f}). "
            f"Pair-force cost relative to the Lambda baseline "
            f"{pair_cost_vs_lambda_baseline(chosen, phi, dim):.2f}x."
        ),
        warnings=notes,
        hoomd_note=p.hoomd_note,
    )

    if box_over_sigma is not None:
        ok, msg = out.box_ok(box_over_sigma)
        if not ok:
            out.warnings.append(msg)
    return out


# =============================================================================
# The neighbour-list buffer — performance, not accuracy
# =============================================================================
def neighbor_list_buffer(r_cut_over_sigma: float,
                         max_step_displacement_over_sigma: float,
                         rebuild_every: int = 20) -> float:
    """Cell-list buffer proposal.
    `r_buff >= rebuild_every * (max displacement per step) * 2`.

    The factor 2 is because two particles can approach each other head-on.
    ⚠ The buffer is a **performance** question, not an accuracy one -- except that
      too small a buffer misses neighbours and the force is quietly wrong. HOOMD
      rebuilds automatically, but the buffer sets the performance.
    """
    if max_step_displacement_over_sigma <= 0:
        raise ValueError("the displacement has to be > 0")
    need = 2.0 * rebuild_every * max_step_displacement_over_sigma
    # clamped to a practical range: too large and the neighbour count grows as
    # (1+buff/r_cut)^3
    return float(min(max(need, 0.1 * r_cut_over_sigma), 0.5 * r_cut_over_sigma))
