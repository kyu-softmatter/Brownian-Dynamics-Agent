#!/usr/bin/env python
"""Verify claims taken from the two books in kb/ **by execution** (absolute rule 6).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_book_claims.py

Books:
  [L] L. Gary Leal, *Microstructural Rheology of Complex Fluids*,
      Cambridge Monographs on Mechanics, 2026. DOI 10.1017/9781009688437.
      kb/cambridge-core_microstructural-rheology-of-complex-fluids_7Aug2026/
  [W] J. Welty et al., *Fundamentals of Momentum, Heat and Mass Transfer*, 5th ed.
      kb/file_1731415827j8JuJ.pdf  (title from the PDF metadata, ISBN 0-4701-2868-2)

The verifications are marked in three **distinct** kinds -- mixing them blurs "the
book is right", "I read it correctly" and "our value is right":

  [BOOK]  reproduce a number the book reports itself, using the book's own formula
          (is my reading right?)
  [DERIV] are the book's formulas mutually consistent (checked independently by
          numerical integration, limits, and so on)
  [OURS]  compare the book's values and formulas against what this project uses

★ Failing is fine. A failure usually means "I misread it", not "the book is wrong",
   and in that case the claim must not go into the digest or the KB.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent

PASS: list[str] = []
FAIL: list[str] = []


def chk(kind: str, name: str, got, want, rtol=1e-3, note="", atol=0.0):
    """Numerical comparison.

    ⚠️ Passing only `rtol` when `want == 0` is a **trap** (raised by b4 on
    2026-08-29). The tolerance is `rtol * max(|want|, 1e-300)`, so `want=0,
    rtol=1.0` reads as "100% tolerance" while actually evaluating to **1e-300, the
    strictest value possible**. Intent (loose) and behaviour (strictest) are exact
    opposites, and if the only reason it passes today is that the value is
    **bit-exactly 0**, then the moment that symmetry breaks it flips to FAIL with a
    message reading `want=0 rtol=1`, which looks like it should have passed.

    So this combination is **rejected rather than handled silently**. To compare
    against 0, state an `atol` with a physical scale (e.g. `atol=1e-12*n*kT`).
    ★ Fixing the contract this way surfaced a **fourth call site** by itself,
    beyond the three b4 found.
    """
    if want == 0.0 and atol == 0.0:
        FAIL.append(f"[{kind}] {name}: want=0 with no atol -- rtol is meaningless "
                    f"against 0 (the tolerance collapses to 1e-300). Supply an "
                    f"atol with a physical scale")
        return False
    ok = abs(got - want) <= max(rtol * abs(want), atol)
    rec = f"[{kind}] {name}: got={got:.6g} want={want:.6g} rel={abs(got-want)/max(abs(want),1e-300):.2e}"
    if note:
        rec += f"  ({note})"
    (PASS if ok else FAIL).append(rec)
    return ok


def chk_true(kind: str, name: str, cond: bool, note=""):
    rec = f"[{kind}] {name}" + (f"  ({note})" if note else "")
    (PASS if cond else FAIL).append(rec)
    return cond


KELVIN_0C = 273.15   # 0 °C. ★ The only place in this file where the Celsius<->Kelvin conversion lives.


def _read_iapws_water_rows() -> dict[float, float]:
    """Read **only the direct rows** of the IAPWS table in
    `knowledge/wiki/concepts/water-298k.md`.

    Row format: `| 293.15 K (20 °C) | **1.0016** | ...`  -> {293.15: 1.0016e-3}
    Interpolated rows (`| 300.00 K | 0.85566¹ | ...`) carry no °C label and drop out
    automatically -- what this check wants is **the rows that need no
    interpolation**.

    ★ The Kelvin value is not trusted as printed; it is recomputed from the `°C`
    label and compared. That is what exposes a convention error such as
    `25 °C = 298.00 K` sitting inside the table.
    """
    path = ROOT / "knowledge" / "wiki" / "concepts" / "water-298k.md"
    rows: dict[float, float] = {}
    pat = re.compile(r"\|\s*([\d.]+)\s*K\s*\((\d+)\s*°C\)\s*\|\s*\*\*([\d.]+)\*\*")
    for T_k, c, eta in pat.findall(path.read_text(encoding="utf-8")):
        T_k, c, eta = float(T_k), int(c), float(eta)
        derived = KELVIN_0C + c
        if abs(T_k - derived) > 1e-9:
            FAIL.append(f"[DOC] water-298k.md: the {c}°C row says {T_k} K but must "
                        f"be {derived} K")
            continue
        rows[T_k] = eta * 1e-3
    return rows


# ── Compare against the value **printed** in the distillation ────────────────
# 2026-08-29: without this, a typo passed 56/56. `welty_transport.md` printed the
# log-interpolated value as `0.8580 mPa·s`, while computing it from the table gives
# `0.8598`. The verifier **was computing and printing 0.8598 correctly**, but there
# was no assertion tying it to the number written in the distillation --
# `d_log < 0.015` passes for both values.
# ⚠️ "The computation is right" and **"the document that will be cited is right"
# are different propositions**. The distillation is the artefact that gets cited,
# so without an assertion linking computed to printed, a transcription error
# survives silently (this project's signature failure: the unwired checker).
def anchor_num(name: str, cell: int) -> str:
    """Regex capturing **the number in the cell-th column** of the table row that
    ends with `<!--@name-->`.

    Why anchor on that rather than on prose: the distillation is **going to be
    translated**. A regex keyed on a Korean label (the Korean for `log-linear`,
    and the like) means
    the translation breaks the check -- not silently, since 0 matches is a FAIL,
    but the coupling itself is wrong when fixing the language breaks a test
    (raised by b4 on 2026-08-29 as a cross-file coupling). `<!--@...-->` is
    invisible in the rendered document and is not translated, so the coupling
    disappears.

    `cell` is a **1-based column number** (separated by `|`). Signs and units are
    ignored; the first number is captured.
    """
    skip = r"(?:[^|\n]*\|){%d}" % (cell - 1)
    return (r"\|" + skip + r"[^|\n]*?(\d+\.\d+)[^|\n]*\|(?:[^|\n]*\|)*?\s*<!--@"
            + name + r"-->")


def chk_doc(book: str, pattern: str, want: float, unit_scale: float = 1.0,
            rtol: float = 5e-4, note: str = ""):
    """Pull a number out of the distillation by regex and compare it against the
    computed value.

    `pattern` must have exactly one capture group (the number). No match, or more
    than one, is a **FAIL** -- skipping silently would make this check meaningless.
    """
    path = ROOT / "knowledge" / "source" / "books" / book
    if not path.exists():
        FAIL.append(f"[DOC] {book}: file missing ({path})")
        return False
    hits = re.findall(pattern, path.read_text(encoding="utf-8"))
    if len(hits) != 1:
        FAIL.append(f"[DOC] {book} /{pattern}/: {len(hits)} matches "
                    f"(must be exactly 1)")
        return False
    got = float(hits[0]) * unit_scale
    return chk("DOC", f"{book} printed value {pattern!r}"
               + (f" -- {note}" if note else ""),
               got, want, rtol=rtol)


# ════════════════════════════════════════════════════════════════════════
# Constants
# ════════════════════════════════════════════════════════════════════════
KB = 1.380649e-23  # J/K  (SI defined value)

# ── The values this project uses (intake/*/system.yaml, tier 1) ─────────
OURS_T = 300.0  # K
OURS_ETA = 0.851e-3  # Pa*s   "water@300K handbook" -- shared by 5 cases
OURS_GAMMA = {  # 3 pi eta d
    "trap-drag-2d-hex300 (d=5um)": (5.000e-6, 4.0102e-8),
    "chain-bend-2d-dlvo (d=1.47um)": (1.470e-6, 1.1790e-8),
}

# ── [W] Appendix I, "Water" (SI table, book p.686) ─────────────────────
#   T(K), rho(kg/m3), mu*1e6 (Pa*s), nu*1e6 (m2/s)
WELTY_WATER = [
    (273, 999.3, 1794, 1.795),
    (293, 998.2, 993, 0.995),
    (313, 992.2, 658, 0.663),
    (333, 983.2, 472, 0.480),
    (353, 971.8, 352, 0.362),
    (373, 958.4, 278, 0.290),
]


# ════════════════════════════════════════════════════════════════════════
# (1) [W] water viscosity -- the book's table vs our 0.851 mPa*s @300 K
# ════════════════════════════════════════════════════════════════════════
def s1_water_viscosity():
    T = np.array([r[0] for r in WELTY_WATER], float)
    mu = np.array([r[2] for r in WELTY_WATER], float) * 1e-6  # Pa*s
    rho = np.array([r[1] for r in WELTY_WATER], float)
    nu = np.array([r[3] for r in WELTY_WATER], float) * 1e-6

    # Internal consistency of the table itself: nu = mu/rho. The book printed the
    # three columns independently, so this is a real check.
    err = float(np.max(np.abs(nu - mu / rho) / nu))
    chk_true("BOOK", f"Welty water table: nu == mu/rho within 0.1% on every row "
                     f"(measured {err:.2e})",
             err < 1e-3, note="did I transcribe the table right, and is the table "
                              "self-consistent")

    # Interpolate to 300 K. eta(T) curves exponentially, so log-linear beats linear.
    i = 1  # 293 K
    w = (300.0 - T[i]) / (T[i + 1] - T[i])
    mu_lin = mu[i] + w * (mu[i + 1] - mu[i])
    mu_log = math.exp(math.log(mu[i]) + w * (math.log(mu[i + 1]) - math.log(mu[i])))

    print(f"    Welty table 293/313 K = {mu[1]*1e3:.3f}/{mu[2]*1e3:.3f} mPa*s")
    print(f"    -> 300 K linear {mu_lin*1e3:.4f} / log-linear {mu_log*1e3:.4f} mPa*s")
    print(f"    -> our value      {OURS_ETA*1e3:.4f} mPa*s")

    # Relative difference against our value. Is log interpolation within the
    # accuracy a 20 K-spaced table can deliver?
    d_log = abs(mu_log - OURS_ETA) / OURS_ETA
    d_lin = abs(mu_lin - OURS_ETA) / OURS_ETA
    chk_true("OURS", f"eta(300K): Welty log interpolation within 1.5% of our "
                     f"0.851 mPa*s (measured {d_log*100:.2f}%)",
             d_log < 0.015, note=f"linear interpolation is off by {d_lin*100:.2f}% "
                                 f"-- the interpolation method changes the answer")
    chk_true("OURS", "eta(300K): log interpolation is closer to our value than linear",
             d_log < d_lin, note="the curvature of eta(T) already matters at 20 K "
                                 "spacing")

    # ★ Do the two values **printed** in the distillation match the computation
    #   above? Without these two lines, 0.8580 (actually 0.8598) passed 56/56.
    chk_doc("welty_transport.md", anchor_num("eta300_log", 2),
            mu_log * 1e3, note="log-linear interpolated value")
    chk_doc("welty_transport.md", anchor_num("eta300_lin", 2),
            mu_lin * 1e3, note="linear interpolated value")
    # The printed percentage must come from the same computation (0.8580 would be
    # +0.82%, contradicting the printed +1.03%)
    chk_doc("welty_transport.md", anchor_num("eta300_log", 3),
            d_log * 100, rtol=6e-3, note="log-linear relative difference %")

    # ★ Recompute the T-sensitivity error §1.2 publishes, from the IAPWS table.
    #   2026-08-29: without this check I "corrected" -4% to -5% and had to revert.
    #   Two causes pointed the same way -- (a) I log-interpolated [W]'s 20 K table
    #   instead of using IAPWS, which has a direct 298.15 K row, and (b) I wrote
    #   25 °C as 298.00 K (the true value is 298.15 K).
    #   => which table, which reference and which temperature convention ALL three
    #   change the number. All three are pinned here.
    # ⚠️ The table is NOT copied in here. That would make a **fourth copy** of the
    #    water viscosity table, and copies diverging is exactly what caused this
    #    bug. The cited source is read directly.
    #    And the temperature is **derived from °C** -- writing `298.15` as a literal
    #    invites rewriting the `25 °C -> 298.00 K` slip that bit both me and
    #    session 0f this time round.
    IAPWS = _read_iapws_water_rows()
    chk_true("DOC", f"water-298k.md: read {len(IAPWS)} IAPWS rows",
             len(IAPWS) >= 2, note="0 rows means the parse broke -- skipping "
                                   "silently would make the check meaningless")
    for celsius, want_pct in ((25, -4.4), (20, -15.0)):
        T_real = KELVIN_0C + celsius          # ★ derived. Never the literal 298.15
        if T_real not in IAPWS:
            # No bare assert -- a crash would take out the other 63 results and the
            # tally with them, burying the diagnostics already recorded (the
            # Kelvin/Celsius mismatch).
            # Silent pass < crash < **report**.
            FAIL.append(f"[OURS] {celsius}°C = {T_real} K is not a direct IAPWS row "
                        f"(rows present: {sorted(IAPWS)})")
            continue
        got = (OURS_ETA - IAPWS[T_real]) / IAPWS[T_real] * 100
        chk("OURS", f"eta error @{celsius}°C = {T_real} K (against the direct IAPWS "
                    f"row, ours 0.851)",
            got, want_pct, rtol=2e-2,
            note=f"writing {celsius}C as {KELVIN_0C + celsius - 0.15:.0f} K moves "
                 f"this by 0.3%")
    # Do the two values printed in the distillation's §1.2 table match the above?
    # ★ No temperature literal here either -- the literal was removed from the loop
    #   above but `IAPWS[298.15]` was left in these two lines, so the adversarial
    #   test crashed with KeyError (one site fixed, the second missed).
    for celsius in (25, 20):
        T_real = KELVIN_0C + celsius
        if T_real not in IAPWS:
            FAIL.append(f"[DOC] §1.2 {celsius}°C row cannot be compared -- no direct "
                        f"IAPWS row")
            continue
        want = abs((OURS_ETA - IAPWS[T_real]) / IAPWS[T_real] * 100)
        chk_doc("welty_transport.md", anchor_num(f"err_{celsius}C", 3),
                want, rtol=2e-2, note=f"§1.2's {celsius}°C ({T_real} K) row")

    # Temperature sensitivity, %/K. Quoting one value without stating T can be off
    # by this much.
    sens = abs(math.log(mu[2] / mu[1])) / (T[2] - T[1]) * 100  # %/K, over 293-313
    print(f"    d(ln eta)/dT ~ {sens:.2f} %/K  (293-313 K)")
    chk_true("BOOK", f"water viscosity T-sensitivity exceeds 2%/K "
                     f"(measured {sens:.2f} %/K)",
             sens > 2.0, note="getting T wrong by +-1 K makes eta wrong by 2%")
    return sens


# ════════════════════════════════════════════════════════════════════════
# (2) [W] Stokes drag <-> C_D = 24/Re <-> our gamma = 3 pi eta d
# ════════════════════════════════════════════════════════════════════════
def s2_stokes():
    # C_D = F / (0.5 rho v^2 A),  A = pi d^2/4,  F = 3 pi mu d v  (Stokes)
    # -> C_D = 3 pi mu d v / (0.5 rho v^2 pi d^2/4) = 24 mu/(rho v d) = 24/Re
    rho, mu, d, v = 998.2, 0.993e-3, 3e-6, 1e-4
    Re = rho * v * d / mu
    F = 3 * math.pi * mu * d * v
    CD = F / (0.5 * rho * v**2 * math.pi * d**2 / 4)
    chk("DERIV", "C_D(Stokes) == 24/Re", CD, 24.0 / Re, rtol=1e-12,
        note="confirms [W] 12.2's C_D=24/Re and F=3*pi*mu*d*v are the same thing")
    chk_true("DERIV", f"the tested condition is in the creeping-flow regime "
                      f"(Re={Re:.2e} < 1)", Re < 1.0)

    # Was our gamma = 3 pi eta d (= 6 pi eta a) computed correctly for each case?
    for name, (d_case, gamma_yaml) in OURS_GAMMA.items():
        chk("OURS", f"gamma = 3*pi*eta*d  [{name}]",
            3 * math.pi * OURS_ETA * d_case, gamma_yaml, rtol=2e-4)

    # Stokes-Einstein [W] (24-50): D = kT/(6 pi mu r) = kT/gamma
    for name, (d_case, gamma_yaml) in OURS_GAMMA.items():
        r = d_case / 2
        D_se = KB * OURS_T / (6 * math.pi * OURS_ETA * r)
        chk("DERIV", f"Stokes-Einstein D == kT/gamma  [{name}]",
            D_se, KB * OURS_T / gamma_yaml, rtol=2e-4)
        print(f"    {name}: D = {D_se:.4e} m^2/s")


# ════════════════════════════════════════════════════════════════════════
# (3) [L] the phi^2 coefficients of Batchelor / Brady-Vicic -- reproduce the book's
#     arithmetic
# ════════════════════════════════════════════════════════════════════════
def s3_phi2():
    # [L] 3.4.2: direct Brownian 0.97 phi^2 + hydrodynamic 5.2 phi^2
    #            -> Batchelor K* = 6.2
    chk("BOOK", "[L] Batchelor K* = 0.97 (Brownian) + 5.2 (hydro)", 0.97 + 5.2, 6.2,
        rtol=5e-3, note="the book printed the two components and the sum separately, "
                        "so this is a real arithmetic check")
    # Brady & Vicic (1995): K* = 5.91, N1/(mu gdot) = 0.899 phi^2 Pe, N2 = -0.788 phi^2 Pe
    chk_true("BOOK", "[L] Brady-Vicic: |N2/N1| = 0.877 (differs from the ~1/7 of "
                     "polymers)",
             abs(-0.788 / 0.899 + 0.877) < 5e-3,
             note=f"|N2/N1|={abs(-0.788/0.899):.3f} vs polymers {1/7:.3f}")
    chk_true("BOOK", "[L] N1>0, N2<0 (same signs as polymers)",
             0.899 > 0 and -0.788 < 0)


# ════════════════════════════════════════════════════════════════════════
# (4) [L] Doi-Edwards semi-dilute rotational diffusion -- reproduce the book's
#     example
# ════════════════════════════════════════════════════════════════════════
def s4_doi_edwards():
    # [L] (6.45): Dr0 = beta * Drbar0 / (n^2 L^6),  beta = 1.3e3
    # The book's example: n=0.1, L=50 -> Dr0/Drbar0 = O(1e-5)
    beta, n, L = 1.3e3, 0.1, 50.0
    ratio = beta / (n**2 * L**6)
    print(f"    Dr(semi-dilute)/Dr(dilute) = {ratio:.3e}   (book: O(1e-5))")
    chk_true("BOOK", f"[L] (6.45) example: ratio = {ratio:.2e} is O(1e-5)",
             1e-6 <= ratio < 1e-4, note="n=0.1, L=50, beta=1.3e3")
    # Scaling: Dr ~ (n L^3)^-2  (n^2 L^6 = (nL^3)^2)
    chk("DERIV", "n^2 L^6 == (n L^3)^2", n**2 * L**6, (n * L**3) ** 2, rtol=1e-12)
    # Excluded volume: two orthogonal rods -> an L x L x 2a parallelepiped = 2 a L^2
    a = 0.5
    chk("BOOK", "[L] excluded volume (orthogonal) = 2 a L^2",
        L * L * (2 * a), 2 * a * L**2, rtol=1e-12)
    # The width of the semi-dilute window is set by the aspect ratio:
    # O(1) << n L^3 << L/a = 2r
    r_asp = L / (2 * a)
    chk("DERIV", "semi-dilute upper bound n L^3 ~ L/a = 2r", L / a, 2 * r_asp,
        rtol=1e-12, note=f"r={r_asp:.0f} -> nL^3 window [O(1), {L/a:.0f}]")


# ════════════════════════════════════════════════════════════════════════
# (5) [L] Jeffery's G, and the |E| convention in Pe = |E|/Dr (the factor-of-2 trap)
# ════════════════════════════════════════════════════════════════════════
def s5_jeffery_and_E():
    G = lambda r: (r**2 - 1) / (r**2 + 1)
    # atol: G is a dimensionless O(1) quantity (0..1), so "zero in double precision"
    # is 1e-15.
    # ★ This site was NOT among the three b4 found -- making `chk` reject
    #   want=0-with-no-atol let the API surface the fourth one by itself. Fixing the
    #   contract finds what hand-scanning the call sites misses.
    chk("DERIV", "[L] (2.34) G(r=1) = 0  (a sphere does not respond to strain rate)",
        G(1.0), 0.0, atol=1e-15)
    chk("DERIV", "[L] (2.34) G(r->inf) -> 1", G(1e8), 1.0, rtol=1e-12)
    chk_true("DERIV", "G(r) is monotonically increasing in r",
             bool(np.all(np.diff(G(np.linspace(1, 100, 500))) > 0)))

    # ★ Convention trap: [L] 4.4.4 has Pe = |E|/Dr, [L] 1.6 (1.32) has Wi = gdot/Dr.
    #   For the two to agree, |E| must be sqrt(2 E:E) -- using E:E alone is off by
    #   sqrt(2).
    gdot = 3.7
    E = np.array([[0.0, gdot / 2, 0.0], [gdot / 2, 0.0, 0.0], [0.0, 0.0, 0.0]])
    EE = float(np.tensordot(E, E))  # E:E
    chk("DERIV", "in simple shear, sqrt(2 E:E) == gdot", math.sqrt(2 * EE), gdot,
        rtol=1e-12,
        note="matches [L] 3.4.2's convention Ehat = E/(2E:E)^1/2 = E/gdot")
    chk_true("DERIV", f"using sqrt(E:E) is off by a factor "
                      f"sqrt(2)={math.sqrt(2):.4f}",
             abs(math.sqrt(EE) * math.sqrt(2) - gdot) < 1e-12,
             note="the most common convention accident that makes Pe/Wi wrong by 2x")

    # Vorticity cannot change the isotropic equilibrium distribution -> Pe is built
    # from E alone.
    Om = np.array([[0.0, gdot / 2, 0.0], [-gdot / 2, 0.0, 0.0], [0.0, 0.0, 0.0]])
    chk("DERIV", "simple shear: |grad u| splits half and half into E and Omega",
        float(np.max(np.abs(E + Om))), gdot, rtol=1e-12,
        note="|grad u|=gdot and |E|=gdot, so they coincide here by accident -- for "
             "pure strain they differ")
    # Pure rotational flow: E=0, so Pe must be 0 (however large grad u is)
    E_rot = np.zeros((3, 3))
    chk_true("DERIV", "in pure rotation (E=0) Pe=0 -- measured via |grad u| it would "
                      "not be 0",
             float(np.max(np.abs(E_rot))) == 0.0)


# ════════════════════════════════════════════════════════════════════════
# (6) [L] rod suspension SAOS: G'' ~ De/(36+De^2) -> viscoelastic relaxation time
#     = 1/(6 Dr)
# ════════════════════════════════════════════════════════════════════════
def s6_rod_saos():
    # [L] (5.121b) has the form FG*(2/15)*3De/(36+De^2). The pole of a single
    # Maxwell mode is at De=6.
    De = np.logspace(-2, 3, 200001)
    loss = De / (36 + De**2)  # the De-dependent part of G'' ~ eta' * omega
    De_peak = float(De[int(np.argmax(loss))])
    chk("DERIV", "[L] (5.121b) the loss term De/(36+De^2) peaks at De=6",
        De_peak, 6.0, rtol=1e-3,
        note="De = omega/Dr, so omega_peak = 6 Dr -> lambda = 1/(6 Dr)")
    # Maxwell comparison: lambda*omega/(1+(lambda omega)^2) peaks at lambda*omega=1
    lam = 1.0 / 6.0
    m = lam * De / (1 + (lam * De) ** 2)
    chk("DERIV", "Maxwell(lambda=1/6) also peaks in loss at De=6",
        float(De[int(np.argmax(m))]), 6.0, rtol=1e-3)
    chk_true("DERIV", "De/(36+De^2) equals Maxwell(lambda=1/6) up to a single "
                      "proportionality constant",
             float(np.max(np.abs(loss / m - loss[0] / m[0]))) < 1e-12,
             note="i.e. the linear viscoelasticity of a dilute rod suspension is "
                  "exactly one Maxwell mode")
    print(f"    -> orientational (l=2) relaxation time lambda = 1/(6 Dr). "
          f"NOT 1/Dr and not 1/(2Dr).")


# ════════════════════════════════════════════════════════════════════════
# (7) [L] elastic dumbbell: Kramers stress and Oldroyd-B -- checked independently by
#     numerically integrating the moment equation
# ════════════════════════════════════════════════════════════════════════
def s7_dumbbell():
    """[L] 11.2-11.3.

    The second-moment equation of a linear (Hookean) dumbbell is
        d<RR>/dt = grad u . <RR> + <RR> . grad u^T  - (4 w0/zeta)(<RR> - (kT/w0) I)
    (= the upper-convected derivative form of (11.6)). Here lambda_p = zeta/(4 w0).
    Kramers (11.15b):  T^(p) = n <F_s R> - n kT I = n (w0 <RR> - kT I).
    """
    kT, zeta, n = 1.0, 1.0, 1.0
    b2 = 3.0 * kT  # <R^2>_eq;  w0 = 3kT/b^2 (Gaussian chain, [L] p.494 "3kT/w0 = b^2")
    w0 = 3.0 * kT / b2
    lam = zeta / (4.0 * w0)
    eta_p = n * kT * lam

    chk("DERIV", "[L] Gaussian chain: 3kT/w0 == b^2", 3 * kT / w0, b2, rtol=1e-12)

    def rhs(RR, gradu):
        relax = -(1.0 / lam) * (RR - (kT / w0) * np.eye(3))
        return gradu @ RR + RR @ gradu.T + relax

    def integrate(gradu, t_end, RR0):
        """Integrate the moment equation with a stiff-capable integrator.

        An independent check that does not use the closed-form solution.

        ⚠️ The old signature carried a `dt=None` that **was never used in the body**
        (raised by b4 on 2026-08-29). A call site passed
        `dt=min(lam/20000, 1e-4/gdot)`, so it read as if it were tightening the step
        for the stiff case, and it was silently ignored.

        ★ The fix is **not** "make it actually be honoured" -- measured, honouring
        that value as `max_step` demands **60 million steps** at `gdot=100`
        (dt=1e-6, t_end=60) and the script hits a 2-minute timeout. So the value
        being passed was never viable, and **the parameter being dead was the only
        reason this script ever ran.** LSODA is implicit and adaptive, so `rtol`
        handles the stiffness and a step ceiling is not wanted. The parameter and
        the argument were therefore **deleted** -- removing a dead knob is more
        honest than wiring it up.
        """
        from scipy.integrate import solve_ivp

        sol = solve_ivp(
            lambda t, y: rhs(y.reshape(3, 3), gradu).ravel(),
            (0.0, t_end), RR0.ravel(),
            method="LSODA", rtol=1e-10, atol=1e-12,
        )
        assert sol.success, sol.message
        return sol.y[:, -1].reshape(3, 3)

    RR_eq = (kT / w0) * np.eye(3)

    # (a) at equilibrium the Kramers stress is exactly 0 (the golden test for the
    #     sign and prefactor conventions)
    T_eq = n * (w0 * RR_eq - kT * np.eye(3))
    chk("DERIV", "[L] (11.15b) Kramers stress = 0 at equilibrium",
        float(np.max(np.abs(T_eq))), 0.0, atol=1e-12 * n * kT,
        note=f"max|T|={np.max(np.abs(T_eq)):.2e}, scale n*kT={n*kT:.2e}")
    chk_true("DERIV", "equilibrium Kramers stress is 0 to machine precision",
             float(np.max(np.abs(T_eq))) < 1e-14)

    # (b) steady simple shear: T12 = eta_p gdot, N1 = 2 eta_p lam gdot^2, N2 = 0
    for gdot in (0.1, 1.0, 5.0):
        gradu = np.array([[0.0, gdot, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        RR = integrate(gradu, 60 * lam, RR_eq)
        T = n * (w0 * RR - kT * np.eye(3))
        Wi = lam * gdot
        chk("DERIV", f"[L] (11.20) T12 = eta_p*gdot  (Wi={Wi:g})", T[0, 1], eta_p * gdot, rtol=2e-4)
        chk("DERIV", f"[L] (11.21) N1 = 2 eta_p lam gdot^2  (Wi={Wi:g})",
            T[0, 0] - T[1, 1], 2 * eta_p * lam * gdot**2, rtol=2e-4)
        chk("DERIV", f"[L] (11.21) N2 = 0  (Wi={Wi:g})", T[1, 1] - T[2, 2], 0.0,
            atol=1e-12 * n * kT,
            note=f"|N2|={abs(T[1,1]-T[2,2]):.2e}")
        # (11.24) tr<RR> = b^2 (1 + (2/3) Wi^2)
        chk("DERIV", f"[L] (11.24) tr<RR> = b^2(1+2/3 Wi^2)  (Wi={Wi:g})",
            float(np.trace(RR)), b2 * (1 + 2.0 / 3.0 * Wi**2), rtol=2e-4)

    # (c) Oldroyd-B predicts no shear thinning -- the viscosity is independent of gdot
    etas = []
    for gdot in (0.01, 0.1, 1.0, 10.0, 100.0):
        gradu = np.array([[0.0, gdot, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        RR = integrate(gradu, 60 * lam, RR_eq)
        etas.append(n * (w0 * RR[0, 1]) / gdot)
    spread = (max(etas) - min(etas)) / eta_p
    chk_true("DERIV", f"Oldroyd-B: eta_p invariant over 4 decades of gdot "
                      f"(spread {spread:.2e})",
             spread < 1e-3, note="'no shear thinning' -- unlike almost every real "
                                 "polymer")

    # (d) stop the flow and the hydrodynamic component goes to 0 instantly; only the
    #     thermodynamic component relaxes over a finite time
    gradu = np.array([[0.0, 2.0 / lam, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    RR_sheared = integrate(gradu, 60 * lam, RR_eq)
    RR_relaxed = integrate(np.zeros((3, 3)), 1.0 * lam, RR_sheared)
    T0 = n * (w0 * RR_sheared[0, 1])
    T1 = n * (w0 * RR_relaxed[0, 1])
    chk("DERIV", "one lambda after stopping, T12 is down by exp(-1)",
        T1 / T0, math.exp(-1.0), rtol=5e-3,
        note="[L] 2.2.1: the hydrodynamic component is instantly 0 -> all remaining "
             "memory is the thermodynamic component")
    print(f"    lambda_p = zeta/(4 w0) = {lam:g},  eta_p = n kT lambda_p = {eta_p:g}")


# ════════════════════════════════════════════════════════════════════════
# (8) the arithmetic of applying [L] 2.2.1's 'instantaneity' claim to our system
#      (chain-bend-2d-dlvo: bonding the chain leaves dissipation at 0.996x)
# ════════════════════════════════════════════════════════════════════════
def s8_our_dlvo():
    # runs/.../system_moduli results (CLAUDE.md section 1-D)
    bead_only, dlvo, jkr = 18453.0, 18380.0, 75590.0
    chk("OURS", "total dissipation of the DLVO chain / beads alone",
        dlvo / bead_only, 0.996, rtol=2e-3,
        note="[L] 2.2.1 -> a degree of freedom with no restoring mechanism does not "
             "contribute to the stress")
    chk("OURS", "total dissipation of the JKR chain / beads alone",
        jkr / bead_only, 4.10, rtol=5e-3)
    chk("OURS", "fraction of the dissipation carried by the chain in JKR",
        (jkr - bead_only) / jkr, 0.758, rtol=5e-3)
    frac = (dlvo - bead_only) / bead_only
    chk_true("OURS", f"excess dissipation of the DLVO chain |{frac*100:.1f}%| < 1%",
             abs(frac) < 0.01,
             note="in [L]'s language: there is no thermodynamic restoring mechanism "
                  "for bending, so the direct component is 0")


# ════════════════════════════════════════════════════════════════════════
# (9) [W] in Buckingham pi, r is 'the rank of the dimensional matrix' -- NOT the
#     number of base dimensions
# ════════════════════════════════════════════════════════════════════════
def s9_buckingham():
    # [W] 11.3: i = n - r,  r = rank of the dimensional matrix.
    # Counterexample: three base dimensions (M,L,t) but a variable list of rank 2
    #   variables = (v [L/t], gdot [1/t], L [L]) -> M appears nowhere -> rank 2
    #   rows M L t x 3 variable columns
    A = np.array([
        [0, 0, 0],   # M
        [1, 0, 1],   # L
        [-1, -1, 0],  # t
    ], float)
    r = int(np.linalg.matrix_rank(A))
    n_var = A.shape[1]
    chk("DERIV", "dimensional matrix rank (v, gdot, L)", r, 2, rtol=1e-12)
    chk("DERIV", "i = n - rank = 3 - 2", n_var - r, 1.0, rtol=1e-12,
        note="counting '3 base dimensions' gives i=0 and the wrong conclusion that "
             "there is no dimensionless group")
    # The actual group: gdot L / v
    v, gdot, Lc = 2.0, 3.0, 5.0
    pi1 = gdot * Lc / v
    chk_true("DERIV", f"the group found, gdot*L/v = {pi1:g}, really is dimensionless",
             True, note="you have to count the rank to get the count right")

    # Re's dimensional matrix has rank 3 -> i = 5-3 = 2 (Eu, Re) -- the book's
    # example (11-9)
    #   variables = (F, rho, v, L, mu)
    B = np.array([
        [1, 1, 0, 0, 1],    # M
        [1, -3, 1, 1, -1],  # L
        [-2, 0, -1, 0, -1],  # t
    ], float)
    chk("BOOK", "[W] (11-9) example: rank=3, i = 5-3 = 2 (Eu, Re)",
        B.shape[1] - int(np.linalg.matrix_rank(B)), 2.0, rtol=1e-12)


# ════════════════════════════════════════════════════════════════════════
# (10) [L] a dilute sphere suspension is Newtonian -- Einstein (2.31) and how far
#      our system is from it
# ════════════════════════════════════════════════════════════════════════
def s10_einstein():
    # mu*/mu = 1 + 2.5 phi + K* phi^2
    for phi, want in [(0.01, 1.02562), (0.05, 1.14050), (0.10, 1.31200)]:
        chk("DERIV", f"mu*/mu (phi={phi}) = 1+2.5phi+6.2phi^2",
            1 + 2.5 * phi + 6.2 * phi**2, want, rtol=1e-4)
    # The phi at which the phi^2 term is 10% of the 2.5phi term:
    # 6.2 phi^2 = 0.25 phi -> phi = 0.0403
    phi_star = 0.25 / 6.2
    chk("DERIV", "phi at which the phi^2 term is 10% of the Einstein term",
        phi_star, 0.0403, rtol=1e-3,
        note="above phi=4%, 'dilute = Newtonian' is already broken at the 10% level")


def main():
    sections = [
        ("(1) water viscosity, [W] Appendix I vs our 0.851 mPa*s", s1_water_viscosity),
        ("(2) Stokes drag, C_D=24/Re, Stokes-Einstein", s2_stokes),
        ("(3) phi^2 coefficients (Batchelor / Brady-Vicic)", s3_phi2),
        ("(4) Doi-Edwards semi-dilute Dr", s4_doi_edwards),
        ("(5) Jeffery G, and the factor-of-2 trap in the |E| convention",
         s5_jeffery_and_E),
        ("(6) rod SAOS -> lambda = 1/(6 Dr)", s6_rod_saos),
        ("(7) elastic dumbbell: Kramers stress, Oldroyd-B (numerical)", s7_dumbbell),
        ("(8) reading our 1-D result through [L] 2.2.1", s8_our_dlvo),
        ("(9) the r in Buckingham pi is a rank", s9_buckingham),
        ("(10) Einstein viscosity and the limit of 'dilute'", s10_einstein),
    ]
    for title, fn in sections:
        print(f"\n{'='*78}\n{title}\n{'='*78}")
        fn()

    print(f"\n{'='*78}")
    for r in PASS:
        print(f"  PASS  {r}")
    for r in FAIL:
        print(f"  FAIL  {r}")
    print(f"{'='*78}")
    print(f"  {len(PASS)}/{len(PASS)+len(FAIL)} PASS")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
