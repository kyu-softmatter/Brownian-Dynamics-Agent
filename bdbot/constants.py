"""Physical constants and medium property tables -- **the single source of truth.**

Pure floats, `math` only. No pint, no numpy, no yaml -- so both halves of the
core can import it without paying for the other half's dependencies
(`import bdbot.constants` costs ~1 ms; `import bdbot` used to cost 190 ms, which
is why `bdbot/__init__` is lazy).

    bdbot.materials   the pint-carrying API      -> imports this
    simbot.units      the `*_si` float API       -> imports this

**Why this module exists (measured 2026-08-29).** The water viscosity table
existed twice, and the two copies had already diverged:

    bdbot.materials.WATER_VISCOSITY[300]   0.851   mPa*s   (3 discrete points)
    simbot.units.water_viscosity_si(300)   0.85566 mPa*s   (interpolated)
    relative difference                    0.545 %

0.545 % in `eta` is 0.545 % in `gamma`, in `D`, and in **every timescale derived
from them.** Neither copy was wrong-by-accident; they came from different
provenances, which is exactly the case that a silent duplicate hides. Both are
kept here, labelled, with the gap computed rather than reconciled by fiat
(`.claude/rules/verify-against-literature.md`: *do not average a disagreement*).

⚠ **Do not spend effort resolving the 0.545 % — resolve `T` first.**
`T = 300 K` is carried as tier 1 across the cases although **nothing measured it**:
it was inherited from the 1-A sketch, which stated no temperature at all. And
`eta(T)` is steep. Measured on `WATER_ETA_SI` below, which has **direct rows** at
both temperatures of interest, so no interpolation is involved:

    what 0.851e-3 is worth if the true T was...   (0.851 - eta)/eta
      298.15 K (25 C)   eta = 0.8900e-3  [row]         -4.38 %
      293.15 K (20 C)   eta = 1.0016e-3  [row]        -15.04 %

    ...and against this table's own eta(300 K) = 8.5566e-4   (eta(T) - base)/base
      298.15 K                                             +4.01 %
      293.15 K                                            +17.06 %

Either way the unmeasured `T` dominates the 0.545 % anchor gap by roughly an order
of magnitude. (`docs/06-roadmap.md` section 2.)

⚠ **A sensitivity number needs FOUR labels, not two: table, interpolation,
baseline direction, and temperature convention.** Three sessions produced three
different answers to one question on 2026-08-29, and every answer was internally
consistent -- which is why none of them was visible from its own side. Measured, at
25 C, one cell per combination:

    table        interp       eta(T) uPa*s   (0.851-eta)/eta   T used
    IAPWS 5 K    (row)             890.000        -4.38 %      298.15  <- use this
    Welty 20 K   log-linear        893.157        -4.72 %      298.15
    Welty 20 K   linear            906.738        -6.15 %      298.15
    Welty 20 K   log-linear        895.918        -5.01 %      298.00
    Welty 20 K   linear            909.250        -6.41 %      298.00

Each error is small and they compound in the same direction: **-4.38 % to -6.41 %
from four labelling choices**, a spread larger than the anchor gap being argued
about. Two of the four cost real time:
  * `25 C = 298.00 K` is wrong; it is **298.15 K**. Worth 0.3 pp here, and it is
    what made a `+4.01 %` in an earlier draft of this file read `+4.40 %`.
  * **Do not interpolate a 20 K table when a 5 K one is on disk.** `WATER_ETA_SI`
    has direct rows at 293.15 and 298.15 K, so at those two temperatures the
    interpolation axis disappears entirely. That is the reason to prefer it -- not
    that it is "more accurate" in the abstract.

Sources
    k_B                 SI 2019 definition (exact)
    eta(T), rho(T)       IAPWS, via knowledge/wiki/concepts/water-298k.md
    eta(300 K) = 0.851   knowledge/source/books/welty_transport.md section 1.1 --
                        an independent 20 K-spaced handbook table supports it to
                        1.03 % under **log** interpolation and 2.91 % under
                        linear, and that distillation concludes the denser source
                        is the better one. ⚠ Always cite the interpolation method
                        with the value: on that table the method is worth 1.9
                        percentage points, far more than the gap being argued about.
"""
from __future__ import annotations

import math

# ── fundamental ─────────────────────────────────────────────────────────────
K_B = 1.380649e-23           # J/K -- SI 2019 definition, exact
KELVIN_0C = 273.15           # K -- exact by definition of the Celsius scale


def celsius(t_c: float) -> float:
    """deg C -> K. **Use this instead of typing a Kelvin literal.**

    ⚠ The reason this function exists: `25 C` was written as `298.00 K` by two
      independent sessions on 2026-08-29, in the same exchange in which one of them
      diagnosed that exact slip in the other's work. It is worth 0.3 percentage
      points in `eta` -- small enough to survive review, and it also **moves 25 C
      off a direct table row**, which silently reintroduces an interpolation axis.
      A literal cannot be wrong in one place and right in another; a function can
      only be wrong everywhere at once.
    """
    return KELVIN_0C + t_c


T_20C = celsius(20.0)        # 293.15 K -- a direct row of WATER_ETA_SI
T_25C = celsius(25.0)        # 298.15 K -- a direct row of WATER_ETA_SI

# ── water, IAPWS (knowledge/wiki/concepts/water-298k.md) ────────────────────
#  5 K spacing over 293-308 K. Do not widen this without re-deriving the
#  interpolation error below -- the whole reason the table is dense is that
#  |d(ln eta)/dT| is >2 %/K here, so a 1 K slip in T is a >2 % slip in every
#  timescale. Per segment, measured: 293.15-298.15 **2.363**, 298.15-303.15
#  **2.202**, 303.15-308.15 **2.062** %/K.
#  ⚠ **Not a table-wide constant, and do not extrapolate a segment value.** The
#    Welty 20 K table, which reaches further, gives 2.058 %/K over 293-313 but
#    **2.957 %/K** over 273-293 -- carrying the 300 K figure down to 273 K would
#    understate the sensitivity by 1.4x. The 300 K segment differs between the two
#    tables by 1.07x (2.202 vs 2.058), which is the accuracy this question has.
WATER_ETA_SI = {                 # T [K] -> eta [Pa*s]
    293.15: 1.0016e-3,
    298.15: 0.8900e-3,
    303.15: 0.7972e-3,
    308.15: 0.7191e-3,
}
WATER_RHO_SI = {                 # T [K] -> rho [kg/m^3]
    293.15: 998.21,
    298.15: 997.05,
    303.15: 995.65,
    308.15: 994.03,
}
WATER_TABLE_RANGE_K = (293.15, 308.15)

# ── the separately-sourced anchor that disagrees with the table ─────────────
#  Read the module docstring before using either number. This is a **record of a
#  disagreement**, not a second table to pick from casually: the 8 `bdbot` cases
#  were run with 0.851 mPa*s (it is written into their spec files and into their
#  `run_id` hashes), and the `simbot` S2 prediction documents sealed the
#  interpolated 8.5566e-4. Changing either one now breaks something real -- a
#  content hash on one side, a `SEALED.sha256` on the other.
WATER_ETA_SOURCED_SI = {
    300.00: 0.851e-3,            # welty_transport.md section 1.1; the 8 bdbot cases use this
}

# A third, independent anchor: the Welty Appendix I table (20 K spacing,
# 273 -> 1794, 293 -> 993, 313 -> 658 uPa*s) interpolated to 300 K. It lands
# **between** the two anchors above and much closer to the IAPWS one.
#  Recomputed here (2026-08-29): log-linear 8.5980e-4, linear 8.7575e-4 Pa*s;
#  +1.03 % / +2.91 % against 0.851e-3, agreeing with the distillation.
#  ★ The two values disagreed when this was first written -- the distillation
#    printed `0.8580` alongside its own `+1.03 %`, and +1.03 % of 0.851 is 0.8598,
#    so the absolute value was the typo. Fixed upstream in `20da7a9`. Recording it
#    because of *why it survived*: the book-claims verifier computed `mu_log`
#    correctly, printed it, and asserted only `|mu_log - 0.851|/0.851 < 0.015` --
#    which passes at 1.03 % and at 0.82 % alike. Nothing tied the computed value to
#    the printed one, so the check had **no power to detect the typo at all.**
#    Same shape as this project's other unwired checkers (`A4` true but unchecked
#    for a month; L4 `step_health` printing `82/82 HEALTHY` while never running).
WATER_ETA_HANDBOOK_20K_SI = {
    "log_linear": 8.5980e-4,
    "linear": 8.7575e-4,
    "source": "welty_transport.md section 1.1, Appendix I table interpolated to 300 K",
}

#  The Appendix I rows themselves, so that nothing else has to copy them.
#   ★ Added 2026-08-29 after a peer session found that its own fix for a
#     table-divergence bug had introduced a **fourth copy of the water table**
#     inside the verifier, and after this file's own test had copied these three
#     Welty rows into a test body. A de-duplication pass that spawns copies while
#     running is the failure it is supposed to prevent.
#   ⚠ This is the **coarse** table. Prefer `WATER_ETA_SI`, which has direct rows at
#     `T_20C` and `T_25C` -- there, the interpolation axis does not exist.
WATER_ETA_HANDBOOK_20K_ROWS_SI = {      # T [K] -> eta [Pa*s]; welty_transport.md sec 1
    273.0: 1794e-6, 293.0: 993e-6, 313.0: 658e-6,
    333.0: 472e-6, 353.0: 352e-6, 373.0: 278e-6,
}

# Interpolation convention. `welty_transport.md` section 1.1 measures that on a
# **20 K**-spaced table linear vs log interpolation differ by 1.9 percentage
# points, and recommends log. On *this* 5 K table the same comparison is
# 0.140 % (measured 2026-08-29) -- inside the spread between the two
# provenances above, so the convention is not what limits accuracy here.
# `linear` is kept because it is what the sealed S2 documents were computed with.
INTERP_CONVENTION = "linear"
LIN_VS_LOG_AT_300K = 1.400e-3    # relative, measured on WATER_ETA_SI


# ── sphere geometry: ONE expression ─────────────────────────────────────────
#  ★ Merged 2026-08-29. Two algebraically identical expressions existed:
#      bdbot.materials.sphere_mass   rho * (pi/6) * d**3
#      simbot.spec.Species.mass_si   rho * (4/3)*pi * (d/2)**3
#    Equal on paper -- `(4/3)pi/8 == pi/6` to the last bit -- but **not equal in
#    floating point**, because the operation order differs. Measured: exactly
#    **1 ULP** apart (1.9626037930826046e-15 vs 1.962603793082605e-15).
#    That is the whole disagreement, and it is the kind that never shows up as a
#    wrong answer -- only as two numbers that will not compare equal, forever.
#  ⚠ Safe to unify because mass reaches **no hash**: it lands at
#    `.system.derived_scales.tau_p` in a spec, `derived_scales` is in
#    `bdbot.runid.DOC_KEYS`, and `physics_only()` strips it before the run_id is
#    taken. Verified on all 278 specs (0 mismatches) before and after.
def sphere_mass_si(rho_si: float, d_si: float) -> float:
    """Mass of a sphere from **density and DIAMETER** [kg].

    ⚠ `d` is a diameter. `simbot` stores a radius, so its caller passes `2*a` --
      the same factor-2 surface as `sphere_drag` vs `stokes_drag_si`.
    """
    return rho_si * (math.pi / 6.0) * d_si**3


def interp_table(table: dict[float, float], T_si: float) -> tuple[float, bool]:
    """Linear interpolation in a `{T: value}` table.

    Returns `(value, extrapolated)`. **The flag is the point**: outside the table
    the caller has to lower the provenance tier from `derived` to `assumed`
    (CLAUDE.md rule 3), and it cannot do that if the function quietly extrapolates.
    """
    ts = sorted(table)
    if T_si in table:
        return table[T_si], False
    if T_si < ts[0] or T_si > ts[-1]:
        lo, hi = (ts[0], ts[1]) if T_si < ts[0] else (ts[-2], ts[-1])
        f = (T_si - lo) / (hi - lo)
        return table[lo] + f * (table[hi] - table[lo]), True
    for lo, hi in zip(ts, ts[1:]):
        if lo <= T_si <= hi:
            f = (T_si - lo) / (hi - lo)
            return table[lo] + f * (table[hi] - table[lo]), False
    raise AssertionError("unreachable")


def interp_table_log(table: dict[float, float], T_si: float) -> tuple[float, bool]:
    """Log-linear interpolation -- the convention `welty_transport.md` recommends.

    Not the default (see `INTERP_CONVENTION`). Exposed so the sensitivity is
    computable instead of asserted.
    """
    ts = sorted(table)
    if T_si in table:
        return table[T_si], False
    extrap = T_si < ts[0] or T_si > ts[-1]
    if extrap:
        lo, hi = (ts[0], ts[1]) if T_si < ts[0] else (ts[-2], ts[-1])
    else:
        lo, hi = next((a, b) for a, b in zip(ts, ts[1:]) if a <= T_si <= b)
    f = (T_si - lo) / (hi - lo)
    return math.exp(math.log(table[lo]) + f * (math.log(table[hi]) - math.log(table[lo]))), extrap


def water_eta_row(T_si: float) -> float:
    """The table's **exact row** at `T_si`, or a diagnosis of why there isn't one.

    ⚠ Do not reach into `WATER_ETA_SI[298.15]` directly. A missing row then raises
      `KeyError: 298.15`, which says nothing about *why* -- and a peer session
      measured that this exact pattern made an adversarial fork die instead of
      report, discarding 63 other results and burying the diagnostic that had been
      recorded one line earlier. Silently passing < crashing < **reporting**.
    """
    if T_si in WATER_ETA_SI:
        return WATER_ETA_SI[T_si]
    nearest = min(WATER_ETA_SI, key=lambda t: abs(t - T_si))
    raise KeyError(
        f"WATER_ETA_SI has no direct row at {T_si} K; nearest is {nearest} K "
        f"(off by {abs(nearest - T_si)} K). If you meant a conventional "
        f"temperature, build it with celsius() -- 25 C is {T_25C}, not 298.0. "
        f"If you meant an interpolated value, call water_viscosity_si().")


def water_viscosity_si(T_si: float) -> tuple[float, bool]:
    """Water dynamic viscosity [Pa*s]. Returns `(value, extrapolated)`.

    Outside 293-308 K this extrapolates and says so -- lower the tier if it does.
    """
    return interp_table(WATER_ETA_SI, T_si)


def water_density_si(T_si: float) -> tuple[float, bool]:
    """Water density [kg/m^3]. Returns `(value, extrapolated)`."""
    return interp_table(WATER_RHO_SI, T_si)


def water_viscosity_sensitivity_per_K(T_si: float) -> float:
    """`|d(ln eta)/dT|` [1/K] from the table -- how much a slip in `T` costs.

    Evaluated on the table interval whose midpoint is nearest `T_si`: 2.363 %/K
    near 295 K, 2.202 %/K near 300 K, 2.062 %/K near 305 K.
    This is the number that makes "room temperature" an
    unusable specification: 293 vs 300 K is a 14 % difference in every timescale
    (`welty_transport.md` section 1.2).
    """
    ts = sorted(WATER_ETA_SI)
    lo, hi = min(zip(ts, ts[1:]), key=lambda p: abs((p[0] + p[1]) / 2 - T_si))
    return abs(math.log(WATER_ETA_SI[hi] / WATER_ETA_SI[lo])) / (hi - lo)


def water_viscosity_provenance_gap(T_si: float) -> dict | None:
    """The gap between the IAPWS table and a separately-sourced anchor at `T_si`.

    `None` when no independent anchor exists at that temperature. Report this
    rather than resolving it: at 300 K it is 0.545 %, and both sides are in use
    by runs that already exist.
    """
    anchor = WATER_ETA_SOURCED_SI.get(round(T_si, 6))
    if anchor is None:
        return None
    table, extrap = water_viscosity_si(T_si)
    out = {
        "T_si": float(T_si),
        "eta_table_si": table,
        "eta_sourced_si": anchor,
        "rel_gap": abs(table - anchor) / anchor,
        "table_extrapolated": extrap,
        # The number that puts the gap in proportion. See the module docstring:
        # arguing about 0.545 % on top of an unmeasured T is false precision.
        "T_uncertainty_cost_per_K": water_viscosity_sensitivity_per_K(T_si),
        "note": ("IAPWS 5 K table (linear interp, what simbot S2 documents sealed) "
                 "vs welty_transport.md section 1.1 (what the 8 bdbot cases ran). "
                 "Do not average -- pick one and record which. And check T first: "
                 "1 K of T is worth ~2.2 %, four times this gap."),
    }
    if round(T_si, 6) == 300.00:
        out["eta_handbook_20K_si"] = WATER_ETA_HANDBOOK_20K_SI
    return out


__all__ = ["K_B", "KELVIN_0C", "celsius", "T_20C", "T_25C",
           "WATER_ETA_SI", "WATER_RHO_SI", "WATER_TABLE_RANGE_K", "water_eta_row",
           "WATER_ETA_HANDBOOK_20K_ROWS_SI",
           "WATER_ETA_SOURCED_SI", "WATER_ETA_HANDBOOK_20K_SI",
           "INTERP_CONVENTION", "LIN_VS_LOG_AT_300K",
           "interp_table", "interp_table_log", "water_viscosity_si",
           "water_density_si", "water_viscosity_sensitivity_per_K",
           "water_viscosity_provenance_gap", "sphere_mass_si"]
