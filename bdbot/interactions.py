"""A catalogue of colloidal interactions plus a recommender -- for when the sketch
has no `U_ij`.

**Why it is needed**: 2 of the 5 sketches left the pair potential blank
(`chain-bend` had a squiggle, `trap-drag` had nothing at all). Rule 3 in CLAUDE.md
forbids inventing one, but stopping at "it is absent" leaves the person with no
idea what to choose.

-> **Recommend, and ask.** Present the standard candidates with their basis, the
   parameters they need and the dimensionless coupling they produce, plus the YAML
   to paste into `resolution`. Anything else, the person writes directly. The
   decision stays with the human; the tool only **shows the options and their
   consequences.**

The standard skeleton for colloids in a single Newtonian medium is **DLVO**:
    U = U_electrostatic (screened Coulomb) + U_van der Waals   (+ excluded volume)
Depletion attraction, field-induced dipoles and contact tangential forces attach to
that, depending on the system.

WARNING: the functional forms are standard, but **the actual numbers (Hamaker
   constant, effective charge, lambda_D) must come from a handbook or the
   literature.** The catalogue supplies not values but "what you have to know"
   (rule 3).
"""
from __future__ import annotations

from dataclasses import dataclass, field

SCHEMA = "bdbot.interactions/0.1"


@dataclass
class Interaction:
    key: str
    name: str
    form: str
    needs: tuple                 # the dimensional parameters that must be specified
    coupling: str                # the dimensionless coupling this interaction produces
    use_when: str
    avoid_when: str
    hoomd: str
    verified: str = ""           # verified by execution in this project (with the basis)
    not_verified: str = ""       # what is NOT confirmed, and hard constraints found by measurement
    tags: frozenset = frozenset()
    note: str = ""


# ══════════════════════════════════════════════════════════════════════
# the catalogue -- standard interactions for colloids in a single medium
# ══════════════════════════════════════════════════════════════════════
CATALOG: dict[str, Interaction] = {
    "pair.none": Interaction(
        key="pair.none", name="no interaction (independent particles)",
        form="U_ij = 0",
        needs=(),
        coupling="—",
        use_when="the observable is a single-particle property (MSD, MSAD, the distribution "
                 "inside a trap) and the density is low. Purely to gain statistics "
                 "from an ensemble of independent particles.",
        avoid_when="the goal is structure (g(r), Voronoi) or collective behaviour.",
        hoomd="just leave the pair potential out of the forces list.",
        verified="trap-2d-5um -- 4 observables matched their analytic solutions",
        tags=frozenset({"single_particle", "dilute"})),

    "pair.wca": Interaction(
        key="pair.wca", name="WCA (pure excluded volume)",
        form="U = 4ε[(σ/r)¹² − (σ/r)⁶] + ε,  r < 2^(1/6)σ",
        needs=("sigma (= particle diameter)", "epsilon (kT by convention)"),
        coupling="none -- effectively hard-core. Structure is set by phi alone.",
        use_when="the charge is well screened so the repulsion is short and all you need is "
                 "to prevent overlap. Almost always used as the core of another "
                 "interaction too.",
        avoid_when="a long-range repulsion or attraction dominates the structure (this alone "
                   "is insufficient).",
        hoomd="md.pair.LJ(r_cut=2^(1/6)sigma, mode='shift')  <- no dedicated class (trap 4)",
        verified="soft-r3 -- used as the core, two-particle force comparison 0.000%",
        tags=frozenset({"excluded_volume", "screened", "structure"})),

    "pair.yukawa": Interaction(
        key="pair.yukawa", name="Yukawa / screened Coulomb (the DLVO repulsive branch)",
        form="U(r)/kT = ε_Y · exp(−κ(r−σ)) / (r/σ),    κ = 1/λ_D",
        needs=("effective surface charge or zeta potential",
               "ionic strength -> Debye length lambda_D",
               "medium permittivity epsilon_r"),
        coupling="eps_Y = U(sigma)/kT (contact coupling) . kappa*sigma = sigma/lambda_D (screening vs particle)",
        use_when="charged colloids plus an electrolyte. **The most standard choice for an "
                 "aqueous single-medium system**, and lambda_D lets the repulsion "
                 "range be tuned experimentally.",
        avoid_when="at very high salt, lambda_D << sigma and it is indistinguishable from "
                   "WCA. At very low salt, lambda_D exceeds the box and minimum image "
                   "breaks.",
        hoomd="md.pair.Yukawa(nlist, default_r_cut=...)  — params: kappa, epsilon",
        verified="",   # not used in this project yet
        tags=frozenset({"charged", "screened", "repulsive", "structure", "dlvo"}),
        note="lambda_D enters the ledger as a new length scale, and a separation check "
             "`lambda_D <= L/4` is added."),

    "pair.soft_power": Interaction(
        key="pair.soft_power", name="soft power-law repulsion r^-n",
        form="U(r)/kT = A (sigma/r)^n     (n=3 is the standard model for a 2D dipolar system)",
        needs=("amplitude A (dimensionless, or length^n)", "exponent n", "cutoff r_c"),
        coupling="Gamma = U(a_mean)/kT = A (sigma/a_mean)^n   * not A alone but A combined with density",
        use_when="field-induced dipoles (paramagnetic or dielectric colloids), or a soft "
                 "corona. n=3 in 2D is the classic model for interfacial "
                 "paramagnetic colloids.",
        avoid_when="for n <= dim the energy integral diverges and a cutoff cannot fix it. "
                   "Small n gives a long tail, and minimum image binds first.",
        hoomd="md.pair.Table (endpoint=False! trap 10) + a WCA core (trap 11)",
        verified="all of soft-r3 -- two-particle 0.000%, energy consistency +0.00 to 0.67%, "
                 "hexagonal crystallization confirmed between Gamma 3 and 30",
        tags=frozenset({"repulsive", "soft", "field_induced", "structure", "dipolar"}),
        note="in 2D, r^-3 tail energy falls off only as 1/r_c, so **an absolute kT cutoff "
             "criterion does not work.** Set the cutoff as a multiple of a_mean and "
             "check convergence."),

    "pair.dlvo": Interaction(
        key="pair.dlvo", name="full DLVO (screened repulsion + van der Waals attraction)",
        form="U = U_Yukawa(r) - (A_H/12)*(sigma/(r-sigma)) type + excluded volume",
        needs=("lambda_D and effective charge (as for Yukawa above)", "Hamaker constant A_H",
               "a contact cut-off distance (to prevent the vdW divergence)"),
        coupling="eps_Y, kappa*sigma, and A_H/kT  (the primary minimum depth)",
        use_when="looking at the onset of aggregation or gelation. When the secondary minimum "
                 "and the energy barrier are the core of the physics.",
        avoid_when="vdW diverges at contact so a core or cut-off is mandatory, and if the "
                   "primary minimum is deep it becomes effectively irreversible "
                   "adhesion and the pair-potential picture collapses "
                   "(-> see contact.adhesive_bending).",
        hoomd="composed with md.pair.Table (Yukawa + vdW + core in one table)",
        verified="",
        tags=frozenset({"charged", "attractive", "aggregation", "dlvo"}),
        note="A_H differs per material and medium combination (roughly 1e-21 to 1e-19 J). "
             "Take it from a handbook or the literature and attach a tier -- do not "
             "invent it."),

    "pair.ao_depletion": Interaction(
        key="pair.ao_depletion", name="Asakura-Oosawa depletion attraction",
        form="U_dep(r) < 0, effective range ~ 2R_g (the polymer size). Depth is "
             "proportional to the polymer osmotic pressure",
        needs=("polymer radius of gyration R_g", "polymer concentration (osmotic pressure)"),
        coupling="depth/kT  .  range ratio 2R_g/sigma  (how short-ranged it is)",
        use_when="adding non-adsorbing polymer to tune the attraction's **range and depth "
                 "independently.** The standard tool for colloid-polymer phase "
                 "separation and gels.",
        avoid_when="if the polymer adsorbs (steric stabilization) the sign flips and the model "
                   "does not apply.",
        hoomd="implemented with md.pair.Table",
        verified="",
        tags=frozenset({"attractive", "short_range", "aggregation", "tunable"})),

    "contact.adhesive_bending": Interaction(
        key="contact.adhesive_bending",
        name="adhesive contact + tangential (bending) stiffness  * NOT a pair potential",
        form="radial:     adhesive contact (JKR) -- effectively an inextensible bond\n"
             "        tangential: bond-angle stiffness kappa_theta = EI/l,  EI = pi*E*a_c^4/4\n"
             "        above a critical moment M_c: slipping / rolling (plastic)",
        needs=("particle modulus E, Poisson ratio nu",
               "adhesion energy W_SL (or contact radius a_c)",
               "critical bending moment M_c"),
        coupling="kappa_theta/kT (is it thermally stiff) . M/M_c (is it inside linear elasticity)",
        use_when="particles **touch and form aggregates, chains or gels.** A central pair "
                 "potential alone cannot reproduce a single bond sustaining a torque.",
        avoid_when="a dispersed system where particles do not touch. Use an ordinary pair "
                   "potential there.",
        hoomd="md.bond.Harmonic(stiff, r0=sigma) + md.angle.Harmonic(k=kappa_theta, t0=pi) "
              "(+ WCA). Plasticity above M_c is not built in -> custom, or restrict "
              "to M<M_c. STOP: **angle.Harmonic gets the force wrong on a nearly "
              "straight chain** -- see not_verified below",
        verified="the kappa_theta = EI/l mapping was verified against discrete chain bending "
                 "(verify/chain_bend_from_papers.py: -0.08% against the beam formula "
                 "at n=51, -0.35% at n=25). The discrete matrix also agrees statically "
                 "with HOOMD angle.Harmonic to 0.55% (verify/verify_angle_matrix.py) "
                 "-- but that comparison is **on energy**",
        not_verified=(
            "STOP: **a HOOMD hard constraint.** md.angle.Harmonic clamps sin theta at "
            "SMALL = sqrt(2)*10^-3 (measured 1.414217e-3). Below that the force is "
            "scaled by sinTheta/SMALL, becoming proportional to kappa*(theta-pi)^2 -- "
            "**quadratic, not linear** -- so the chain is softer and more nonlinear "
            "than it should be. **The energy is exact throughout**, so an energy check "
            "passes. With t0=pi the equilibrium itself is at sin theta=0, so **the "
            "stiffer the chain (larger kappa_theta/kT) the worse it gets** -- and this "
            "interaction lives in exactly that regime. When designing, confirm at the "
            "response amplitude that **every** bond angle has min|theta-pi| > SMALL "
            "(not just the maximum -- a chain's end angles are an order of magnitude "
            "smaller than the middle). chain-bend-2d-oscill had all 23 angles inside "
            "the broken region and execution was refused. Workarounds: implement "
            "force.Custom directly (26x slower), lower kappa_theta, or solve "
            "analytically if the regime is linear. Neither angle.Table nor "
            "CosineSquared is a substitute. Reproduce: "
            "verify/verify_angle_force_small_theta.py (skill bd-hoomd trap 15)"),
        tags=frozenset({"bonded", "aggregate", "contact", "gel", "tangential"}),
        note="the central result of Pantina & Furst (PRL 94 138301; Langmuir 24 1141). A "
             "single bond sustaining a torque is the measurement, and a central-force "
             "model like DLVO cannot produce it -- the paper says so explicitly. This "
             "project then executed the prediction side: see chain-bend-2d-dlvo in "
             "docs/04-cases.md."),
}


# ══════════════════════════════════════════════════════════════════════
# the recommender
# ══════════════════════════════════════════════════════════════════════
INTERACTION_SYMBOLS = frozenset({
    "U_ij", "pair_potential", "interaction", "interactions", "pair", "U", "potential",
})


def looks_like_interaction(symbol: str) -> bool:
    s = (symbol or "").strip()
    return s in INTERACTION_SYMBOLS or "potential" in s.lower() or s.lower().startswith("u_")


def infer_tags(obs) -> set:
    """Extract context tags from observation.yaml. Only what is written, never a guess.

    WARNING: the keyword lists below deliberately contain **both English and Korean**
    terms. The `intake/*/observation.yaml` files were written in Korean and are input
    records that are not rewritten, so removing the Korean keywords would silently
    stop matching every existing case. Add English terms; do not remove Korean ones.
    """
    raw = obs.raw if hasattr(obs, "raw") else dict(obs)
    text = " ".join([
        str(raw.get("system_guess", "")),
        str(raw.get("raw_transcription", "")),
        " ".join(str(e) for e in (raw.get("entities") or [])),
        " ".join(str(g) for g in (raw.get("stated_goals") or [])),
    ]).lower()
    tags: set[str] = set()
    kinds = {str((e or {}).get("kind", "")).lower() for e in (raw.get("entities") or [])
             if isinstance(e, dict)}

    if "bonded_interaction" in kinds or any(w in text for w in
                                           ("chain", "사슬", "aggregate", "응집", "비드", "bead")):
        tags |= {"bonded", "aggregate", "contact"}
    if any(w in text for w in ("rdf", "voronoi", "structure", "구조", "결정", "configuration")):
        tags.add("structure")
    if any(w in text for w in ("msd", "msad", "확산", "diffusion")):
        tags.add("single_particle")
    if any(w in text for w in ("charge", "전하", "zeta", "salt", "mgcl", "nacl", "전해질")):
        tags.add("charged")
    if any(w in text for w in ("magnetic", "자성", "dipole", "쌍극자", "field", "외부장")):
        tags.add("field_induced")
    # * The word 'rheology' alone must not recommend a contact or gel model.
    #   trap-drag is written up as "active microrheology" but it is a dragged probe,
    #   not a contact system (a false positive caught by running the recommender over
    #   the real cases).
    #   Only raise a tangential interaction as a candidate when a modulus measurement
    #   appears **together with evidence of contact or aggregation.**
    rheo = any(w in text for w in ("g'", "g''", "modulus", "탄성률", "rheolog", "레올로지"))
    if rheo and ({"bonded", "aggregate", "contact"} & tags):
        tags |= {"gel", "tangential"}
    elif rheo:
        tags.add("rheology")
    if any(w in text for w in ("hexagonal", "육방", "lattice", "격자")):
        tags |= {"structure", "repulsive"}
    if any(w in text for w in ("trap", "트랩", "tweezer", "집게")):
        tags.add("trapped")
    return tags


def recommend(obs, top: int = 3) -> list:
    """A list of (Interaction, score, basis). The score is the tag-match count plus a
    bonus for having been verified in this project.
    """
    tags = infer_tags(obs)
    scored = []
    for it in CATALOG.values():
        hit = tags & it.tags
        if not hit:
            continue
        score = len(hit) + (0.5 if it.verified else 0.0)
        why = f"context tags matched: {', '.join(sorted(hit))}"
        if it.verified:
            why += " . verified in this project"
        scored.append((it, score, why))
    scored.sort(key=lambda x: -x[1])
    return scored[:top], tags


def yaml_snippet(it: Interaction, symbol: str = "U_ij") -> str:
    """A draft to paste into `observation.yaml`'s resolution."""
    needs = "\n".join(f"    #   {n}" for n in it.needs) or "    #   (no extra parameters)"
    return (f'  - symbol: {symbol}\n'
            f'    resolution: |\n'
            f'      {it.key} — {it.name}\n'
            f'      {it.form.splitlines()[0]}\n'
            f'      coupling: {it.coupling}\n'
            f'      * the physical system is not settled until these are filled in '
            f'(each needs a source and a tier):\n'
            f'{needs}\n'
            f'    confirmed_by: user        # <- fill this in once confirmed\n')


def render_suggestion(obs, case_name: str = "") -> str:
    """Recommend, and ask. The decision is the human's."""
    recs, tags = recommend(obs)
    open_int = [m for m in obs.open_missing
                if looks_like_interaction(str(m.get("symbol", "")))]

    L: list[str] = []
    w = L.append
    w("=" * 80)
    w(f"intake suggest — {case_name or obs.path.parent.name}")
    w("=" * 80)
    w("the standard skeleton for colloids in a single medium is DLVO:")
    w("  U = U_electrostatic (screened Coulomb) + U_van der Waals  (+ excluded volume)")
    w("depletion attraction, field-induced dipoles and contact tangential forces "
      "attach to that, depending on the system.")
    w("")
    w(f"context tags read from this case: {', '.join(sorted(tags)) or '(none)'}")

    if not open_int:
        w("")
        w("no interaction is recorded as an unresolved gap.")
        if obs.raw.get("missing_required"):
            w("(either it is already resolved, or this system has no interaction at all)")
    else:
        w("")
        w(f"{len(open_int)} unspecified interaction(s): "
          + ", ".join(str(m.get("symbol")) for m in open_int))

    w("")
    w("─" * 80)
    w("RECOMMENDATIONS (by basis)")
    w("─" * 80)
    for i, (it, score, why) in enumerate(recs, 1):
        mark = "* top" if i == 1 else f"  #{i}"
        w("")
        w(f"{mark}  {it.key}  —  {it.name}")
        w(f"    form     {it.form.splitlines()[0]}")
        for extra in it.form.splitlines()[1:]:
            w(f"             {extra.strip()}")
        w(f"    coupling {it.coupling}")
        w(f"    use when   {it.use_when}")
        w(f"    avoid when {it.avoid_when}")
        w(f"    HOOMD    {it.hoomd}")
        w(f"    verified {it.verified or 'not used in this project yet (unverified)'}")
        if it.needs:
            w(f"    * values to fill in: {', '.join(it.needs)}")
        if it.note:
            for ln in _wrap(it.note, 72):
                w(f"    note     {ln}" if ln == _wrap(it.note, 72)[0] else f"             {ln}")
        w(f"    basis    {why}")

    w("")
    w("=" * 80)
    w("QUESTION -- how would you like to proceed?")
    w("=" * 80)
    if recs:
        top = recs[0][0]
        w(f"  (a) use the top candidate `{top.key}`")
        w("      -> give me the 'values to fill in' above and I will write system.yaml.")
        w("  (b) use a different candidate -- pick from the list above, or name one not in "
          "the catalogue")
        w("  (c) no interaction (`pair.none`) -- when only single-particle observables are wanted")
        w("")
        w("  a draft to paste in (the missing_required entry in observation.yaml):")
        w("")
        for ln in yaml_snippet(top, str(open_int[0].get("symbol")) if open_int else "U_ij").splitlines():
            w(f"  {ln}")
    else:
        w("  the context tags did not narrow the candidates. See the full catalogue:")
        w("    $PY -m bdbot.cli interactions list")
    w("")
    w("  WARNING: either way, **the actual numbers are not invented.** Supply handbook "
      "or literature")
    w("     values, or if I propose one it is recorded in system.yaml marked tier 3 "
      "(arbitrary assumption).")
    w("=" * 80)
    return "\n".join(L)


def render_catalog() -> str:
    L = ["=" * 80, "COLLOIDAL INTERACTION CATALOGUE (single medium)", "=" * 80,
         "DLVO = screened Coulomb + van der Waals is the standard skeleton.", ""]
    for it in CATALOG.values():
        v = "verified" if it.verified else "unverified"
        if it.not_verified:
            v += " . HAS CONSTRAINTS"
        L.append(f"  {it.key:<28} {it.name}")
        L.append(f"  {'':<28} {it.form.splitlines()[0]}")
        L.append(f"  {'':<28} coupling: {it.coupling}   [{v}]")
        if it.not_verified:                    # a hard constraint is never folded away
            for ln in _wrap(it.not_verified, 74):
                L.append(f"  {'':<28} {ln}")
        L.append("")
    L.append("details: $PY -m bdbot.cli intake suggest <case>")
    L.append("=" * 80)
    return "\n".join(L)


def _wrap(s: str, width: int) -> list:
    words, lines, cur = s.split(), [], ""
    for wd in words:
        if len(cur) + len(wd) + 1 > width:
            lines.append(cur)
            cur = wd
        else:
            cur = f"{cur} {wd}".strip()
    if cur:
        lines.append(cur)
    return lines or [""]


__all__ = ["SCHEMA", "Interaction", "CATALOG", "recommend", "infer_tags",
           "render_suggestion", "render_catalog", "yaml_snippet",
           "looks_like_interaction", "INTERACTION_SYMBOLS"]
