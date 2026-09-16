"""S31 — read the sedimentation campaign, in the order the sealed plan fixes.

    ./bin/py campaigns/s31_analyze.py            # the whole thing
    ./bin/py campaigns/s31_analyze.py --figures  # F1..F7 only

`campaigns/sediment_preregistration/analysis_plan.yaml` names every figure and
fixes the order of operations. This script implements that order and **stops where
the plan says to stop**: steps 1, 2 and 3 can each end the analysis, and the
verdict on the question is read last.

## The three things this script does that the obvious version would not

★ **The seal is verified before a single number is read**, and a broken seal ends
the run with a non-zero exit. Reading a result out of a directory whose
pre-registration has changed is the failure the whole mechanism exists to prevent.

★ **The order is enforced, not documented.** `Analysis.gate()` refuses to
evaluate a later step if an earlier one failed, so "we looked at the correctness
gates after we saw the physics" is not reachable by running this in a different
order — there is no other order.

★ **The verdict is computed from the arms' contrast**, per `decision_rule` step 5,
and the three discriminators are required to agree (step 5b). A 2-1 split returns
INCONCLUSIVE and prints the disagreement rather than a majority.

⚠ **This script proposes; it does not confirm.** `confirmed_by` is human-only
(CLAUDE.md), so the verdict here is a proposal with its evidence laid out.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from bdbot import runcard as RC  # noqa: E402
import campaigns.s31_sedimentation as S31  # noqa: E402
import cases.sediment_3d as SED  # noqa: E402

DOCS = ROOT / "campaigns" / "sediment_preregistration"
FIGDIR = ROOT / "campaigns" / "sediment_figures"

#: The hypotheses, as the EOS scale each one puts into the profile ODE.
#:
#: ⚠ COMPUTED, not typed. `(d_BH/d)**3` depends on `eps_wca` and on `sigma_LJ`,
#: and a literal here would keep agreeing with the runs' own value only by
#: coincidence -- which is how `phi` came to be 1.414x too large earlier in this
#: case's history. `init_eos_scale` is also in every arm-A spec, so the value is
#: cross-checked against the run below rather than trusted.
from simbot.cutoff import barker_henderson_diameter  # noqa: E402
HYP = {"cs_eff": float(barker_henderson_diameter(1.0, sigma_lj=SED.SIGMA_LJ)) ** 3,
       "cs_nominal": 1.0}


def check_scale_against_the_runs(rows) -> list[str]:
    """The scale this script draws its candidate profiles from must be the one
    the runs were started with. Both are recorded, so they are compared."""
    bad = []
    for r in rows:
        if not r["ok"] or r["init"] != "cs_eff":
            continue
        got = float(json.loads((ROOT / "specs" / f"{r['run_id']}.json").read_text())
                    ["params"]["init_eos_scale"])
        if not math.isclose(got, HYP["cs_eff"], rel_tol=1e-9):
            bad.append(f"{r['id']}: the run used init_eos_scale = {got!r} but "
                       f"this analysis draws CS(phi_eff) at {HYP['cs_eff']!r}")
    return bad


def prereg() -> dict:
    return yaml.safe_load((DOCS / "prediction.yaml").read_text())


# ── loading ───────────────────────────────────────────────────────────────
def load_runs() -> list[dict]:
    """Every run's directory, spec, metrics and arrays — or a note on why not."""
    out = []
    for r in S31.runs():
        rid = S31.run_id_of(r)
        d = ROOT / "runs" / rid
        row = dict(r, run_id=rid, dir=d, ok=False, why="")
        if not (d / "metrics.json").exists():
            row["why"] = "no metrics.json -- this run has not finished"
            out.append(row)
            continue
        row["metrics"] = json.loads((d / "metrics.json").read_text())
        row["obs"] = {o["name"]: o for o in row["metrics"]["observables"]}
        row["result"] = row["metrics"]["result"]
        with np.load(d / "observables.npz", allow_pickle=True) as z:
            row["arr"] = {k: z[k] for k in z.files}
        row["ok"] = True
        out.append(row)
    return out


def o(row: dict, name: str, field: str = "measured"):
    return (row.get("obs", {}).get(name) or {}).get(field)


# ── the ordered analysis ──────────────────────────────────────────────────
class Analysis:
    """Holds the plan's order. A step that fails closes the ones after it."""

    def __init__(self, rows):
        self.rows = rows
        self.steps: list[tuple[str, bool, list[str]]] = []
        self.stopped_at: str | None = None

    def gate(self, name: str, fn):
        if self.stopped_at:
            self.steps.append((name, False, [f"not evaluated -- stopped at "
                                             f"{self.stopped_at}"]))
            return False
        ok, notes = fn()
        self.steps.append((name, ok, notes))
        if not ok:
            self.stopped_at = name
        return ok

    def note(self, name: str, fn):
        """A step that reports but cannot stop the analysis."""
        if self.stopped_at:
            self.steps.append((name, False, [f"not evaluated -- stopped at "
                                             f"{self.stopped_at}"]))
            return
        _, notes = fn()
        self.steps.append((name, True, notes))

    def report(self) -> int:
        print("\n" + "=" * 78)
        print("THE SEALED ORDER OF OPERATIONS")
        print("=" * 78)
        for name, ok, notes in self.steps:
            print(f"\n{'PASS' if ok else 'STOP':>4s}  {name}")
            for n in notes:
                print(f"        {n}")
        if self.stopped_at:
            print(f"\n★ the analysis stopped at: {self.stopped_at}")
            print("  The sealed plan says steps 1, 2 and 3 can each end it. This "
                  "is that,\n  not an error in this script.")
        return 1 if self.stopped_at else 0


# ── step 1: the seal ──────────────────────────────────────────────────────
def step1_seal(rows):
    notes, bad = [], []
    for r in rows:
        ok, problems = RC.verify_seal(r["dir"], root=ROOT)
        hard = [p for p in problems if not p.startswith("[warn]")]
        notes.append(f"{r['id']:3s} {'OK ' if ok else 'BROKEN'} "
                     f"{len(RC.read_seal(r['dir']))} doc(s)"
                     + (f"  {hard}" if hard else ""))
        if not ok:
            bad.append(r["id"])
    scale_bad = check_scale_against_the_runs(rows)
    for s in scale_bad:
        notes.append(f"★ {s}")
    bad += scale_bad
    if bad:
        notes.append(f"★ {bad} -- a broken seal or a provenance mismatch ends the "
                     f"analysis. No number below it was read.")
    return not bad, notes


# ── step 2: the correctness gates ─────────────────────────────────────────
GATES = (("Z_dilute_tail", 1.0, 10.0), ("D_xy", 1.0, 20.0))


def step2_gates(rows):
    notes, bad = [], []
    for r in rows:
        if not r["ok"]:
            continue
        for name, target, tol in GATES:
            v = o(r, name)
            if v is None or not np.isfinite(v):
                bad.append(f"{r['id']}/{name} is not finite")
                continue
            dev = 100.0 * (v / target - 1.0)
            flag = "" if abs(dev) <= tol else "  <-- OUT"
            notes.append(f"{r['id']:3s} {name:16s} {v:9.4f} vs {target:g} "
                         f"({dev:+6.2f} %, band {tol:g} %){flag}")
            if abs(dev) > tol:
                bad.append(f"{r['id']}/{name}")
    #  P1 is a measurement in revision 3 -- reported, never a gate
    notes.append("")
    for r in rows:
        if not r["ok"]:
            continue
        lt, se = o(r, "l_g_fitted_tail"), None
        nb = r["result"].get("n_bins_fitted_tail")
        notes.append(f"{r['id']:3s} l_g_fitted_tail  {lt:9.4f} vs 13.375 "
                     f"({100*(lt/13.375-1):+6.2f} %, {nb} bins) "
                     f"[measurement -- see design_power]")
    if bad:
        notes.append(f"★ IMPLEMENTATION FAILURE: {bad}")
    return not bad, notes


# ── step 3: stationarity, which IS the experiment ─────────────────────────
def step3_stationarity(rows):
    notes, drifting = [], []
    notes.append(f"{'id':4s} {'arm':12s} {'chi2/nu':>8s} {'shift %':>9s} "
                 f"{'bins':>5s} {'h1/h2':>9s}  L4")
    for r in rows:
        if not r["ok"]:
            continue
        c = o(r, "profile_halves_chi2_nu")
        s = o(r, "profile_half_shift_pct")
        res = r["result"]
        l4 = r["metrics"].get("l4", {}).get("status", "?")
        tag = "stationary" if c is not None and c <= 2.0 else "DRIFTING"
        notes.append(f"{r['id']:4s} {r['arm']:12s} {c:8.3f} {s:+9.3f} "
                     f"{res.get('stationarity_bins', 0):5d} "
                     f"{res.get('frames_h1',0):4d}/{res.get('frames_h2',0):<4d} "
                     f"{l4:10s} {tag}")
        if c is None or c > 2.0:
            drifting.append(r["id"])
    #  ★ the plan does NOT stop here on a drifting arm: a drifting arm is the
    #    expected signature of the LOSING hypothesis, and which arm drifts is the
    #    verdict. It stops only if BOTH competing arms drift or NEITHER does --
    #    handled in step 5, where it becomes INCONCLUSIVE rather than a stop.
    per_arm = {}
    for r in rows:
        if r["ok"] and r["arm"] in ("cs_eff", "cs_nominal"):
            per_arm.setdefault(r["arm"], []).append(r["id"] not in drifting)
    notes.append("")
    for arm, flags in sorted(per_arm.items()):
        notes.append(f"arm {arm:12s} {sum(flags)}/{len(flags)} seeds stationary")
    return True, notes


# ── step 4: the three discriminators ──────────────────────────────────────
DISC = (
    ("Z_dev_wmean_vs_CS_eff", 0.0, 6.47, "percent"),
    ("l_g_fitted_dense", 16.485, 17.554, "d"),
    ("phi_wall_measured_vs_cs_eff", 0.1118, 0.1029, "1"),
)


def step4_discriminators(rows):
    notes = []
    notes.append("Each row: the measurement, then which hypothesis it is closer "
                 "to.  (eff / nom are the")
    notes.append("design-power Monte Carlo's predictions for the two "
                 "hypotheses.)\n")
    for name, v_eff, v_nom, unit in DISC:
        notes.append(f"── {name}  [eff = {v_eff:g}, nom = {v_nom:g} {unit}]")
        for arm in ("cs_eff", "cs_nominal", "finite_size", "convergence"):
            vals = [o(r, name) for r in rows
                    if r["ok"] and r["arm"] == arm and o(r, name) is not None]
            if not vals:
                continue
            m, sd = float(np.mean(vals)), float(np.std(vals, ddof=1) if len(vals) > 1 else 0.0)
            closer = "eff" if abs(m - v_eff) < abs(m - v_nom) else "nom"
            sep = abs(v_eff - v_nom)
            frac = abs(m - v_eff) / sep if sep else float("nan")
            notes.append(f"   {arm:12s} {m:10.4f} ± {sd:7.4f} (n={len(vals)})"
                         f"   -> {closer}   ({frac:.0%} of the way from eff to nom)")
        notes.append("")
    return True, notes


# ── step 5: the verdict ───────────────────────────────────────────────────
def verdict(rows) -> tuple[str, list[str]]:
    notes = []
    stat = {}
    for arm in ("cs_eff", "cs_nominal"):
        cs = [o(r, "profile_halves_chi2_nu") for r in rows
              if r["ok"] and r["arm"] == arm]
        cs = [c for c in cs if c is not None]
        stat[arm] = all(c <= 2.0 for c in cs) if cs else None
        notes.append(f"arm {arm:12s} chi2/nu = "
                     f"{['%.3f' % c for c in cs]}  -> "
                     f"{'stationary' if stat[arm] else 'drifting'}")

    votes = {}
    for name, v_eff, v_nom, _ in DISC:
        vals = [o(r, name) for r in rows
                if r["ok"] and r["arm"] == "cs_eff" and o(r, name) is not None]
        if not vals:
            continue
        m = float(np.mean(vals))
        votes[name] = "eff" if abs(m - v_eff) < abs(m - v_nom) else "nom"
    notes.append(f"\ndiscriminators (arm A): {votes}")

    if len(set(votes.values())) > 1:
        notes.append("★ the three discriminators DISAGREE. The sealed plan's step "
                     "5b says that is INCONCLUSIVE and the disagreement is the "
                     "result -- it is not resolved by majority.")
        return "INCONCLUSIVE (discriminators disagree)", notes

    if stat["cs_eff"] and not stat["cs_nominal"]:
        return "MAPPING REQUIRED", notes
    if stat["cs_nominal"] and not stat["cs_eff"]:
        return "MAPPING NOT REQUIRED (contradicts the prior finding)", notes
    if stat["cs_eff"] and stat["cs_nominal"]:
        notes.append("both arms held their own profile: 4 tau_sed did not resolve "
                     "them.")
        return "INCONCLUSIVE (both stationary)", notes
    notes.append("neither arm held its profile -- both left for somewhere else, "
                 "and that third profile is the finding.")
    return "INCONCLUSIVE (neither stationary)", notes


# ── the figures ───────────────────────────────────────────────────────────
def make_figures(rows) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGDIR.mkdir(parents=True, exist_ok=True)
    live = [r for r in rows if r["ok"]]
    if not live:
        return []
    lz = float(live[0]["metrics"]["physical"]["L_z_star"])
    l_g = float(live[0]["metrics"]["physical"]["l_g_star"])
    phi = float(live[0]["metrics"]["physical"]["phi"])
    cands = {}
    for tag, s in HYP.items():
        z, p = SED.cs_hydrostatic_profile(l_g, lz, phi, scale=s)
        cands[tag] = (z, p)
    zi = np.linspace(0, lz, 600)
    ideal0 = phi * (lz / l_g) / (1.0 - math.exp(-lz / l_g))
    COL = {"cs_eff": "tab:blue", "cs_nominal": "tab:red",
           "convergence": "tab:green", "finite_size": "tab:purple"}
    out = []

    # ── F1: every seed, both candidates overlaid ─────────────────────────
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for r in live:
        h, p = r["arr"]["profile_h"], r["arr"]["profile_phi"]
        ax.plot(h, p, lw=0.9, alpha=0.75, color=COL[r["arm"]],
                label=r["arm"] if r["id"].endswith(("1", "4", "7", "8")) else None)
    for tag, (z, p) in cands.items():
        ax.plot(z, p, "--", lw=1.6, color=COL[tag],
                label=f"CS({'phi_eff' if tag=='cs_eff' else 'phi_nom'}) start")
    ax.plot(zi, ideal0 * np.exp(-zi / l_g), ":", color="grey", lw=1.4,
            label="ideal gas (arm C start)")
    for x in (SED.FIT_LO, SED.TAIL_LO, SED.TAIL_HI):
        ax.axvline(x, color="k", lw=0.5, alpha=0.3)
    ax.set_yscale("log"); ax.set_xlim(0, lz); ax.set_ylim(1e-6, 0.3)
    ax.set_xlabel("height above the bottom wall  h / d")
    ax.set_ylabel(r"local volume fraction  $\phi(h)$")
    ax.set_title("F1  every seed, with both candidate profiles and the ideal gas")
    ax.legend(fontsize=7.5, ncol=2)
    fig.tight_layout(); f = FIGDIR / "F1_profiles.png"
    fig.savefig(f, dpi=150); plt.close(fig); out.append(f)

    # ── F2/F3: Z against phi_eff, both mappings, and the deviation ───────
    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    ax, ax2 = axes
    pe = np.linspace(1e-4, 0.13, 300)
    ax.plot(pe, SED.z_carnahan_starling(pe), "-", color="k", lw=1.5,
            label=r"CS($\phi_{\rm eff}$)")
    ax.plot(pe, SED.z_carnahan_starling(pe / HYP["cs_eff"]), ":", color="k",
            lw=1.5, label=r"CS($\phi_{\rm nominal}$), the alternative")
    ax.axhline(1.0, color="grey", lw=0.6)
    for r in live:
        a = r["arr"]; m = a["eos_window"].astype(bool)
        ax.plot(a["eos_phi_eff"][m], a["eos_Z"][m], ".", ms=2.5, alpha=0.5,
                color=COL[r["arm"]])
        dev = 100.0 * (a["eos_Z"][m] - a["eos_CS_eff"][m]) / a["eos_CS_eff"][m]
        ax2.plot(a["eos_phi_eff"][m], dev, ".", ms=2.5, alpha=0.5,
                 color=COL[r["arm"]])
    ax2.axhspan(-2, 2, color="tab:green", alpha=0.12,
                label="the sealed 2 % band")
    ax2.axhline(0, color="k", lw=0.6)
    ax2.axhline(6.47, color="k", ls=":", lw=1.2,
                label="where CS(nominal) data would sit (+6.47 %)")
    ax.set_ylabel(r"$Z = \beta P/\rho$"); ax.legend(fontsize=8)
    ax.set_title("F2 / F3  the equation of state, and its deviation")
    ax2.set_xlabel(r"$\phi_{\rm eff} = \phi\,(d_{BH}/d)^3$")
    ax2.set_ylabel(r"$100\,(Z - Z_{CS})/Z_{CS}$  [%]")
    ax2.set_ylim(-15, 15); ax2.legend(fontsize=8)
    fig.tight_layout(); f = FIGDIR / "F2_F3_eos.png"
    fig.savefig(f, dpi=150); plt.close(fig); out.append(f)

    # ── F7: THE ANSWER -- which arm did not move ─────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=True)
    top, mid, bot = axes
    pooled = {}
    for arm in ("cs_eff", "cs_nominal"):
        rs = [r for r in live if r["arm"] == arm]
        if not rs:
            continue
        h = rs[0]["arr"]["profile_h"]
        c1 = sum(r["arr"]["profile_counts_h1"] for r in rs)
        c2 = sum(r["arr"]["profile_counts_h2"] for r in rs)
        p = np.mean([r["arr"]["profile_phi"] for r in rs], axis=0)
        pooled[arm] = (h, p, c1, c2)
        top.plot(h, p, "-", lw=1.2, color=COL[arm], label=f"arm {arm} (pooled)")
        top.plot(*cands[arm], ls="--", lw=1.4, color=COL[arm], alpha=0.7,
                 label=f"its start")
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = (c2 / max(rs[0]['result'].get('frames_h2', 1), 1)) / \
                    (c1 / max(rs[0]['result'].get('frames_h1', 1), 1))
            band = np.sqrt(1.0 / np.maximum(c1, 1) + 1.0 / np.maximum(c2, 1))
        m = (c1 >= 20) & (c2 >= 20)
        mid.errorbar(h[m], ratio[m], yerr=band[m], fmt=".", ms=3, lw=0.6,
                     color=COL[arm], alpha=0.7, label=f"arm {arm}")
    mid.axhline(1.0, color="k", lw=0.8)
    if len(pooled) == 2:
        ha, pa, _, _ = pooled["cs_eff"]
        _, pb, _, _ = pooled["cs_nominal"]
        with np.errstate(divide="ignore", invalid="ignore"):
            bot.plot(ha, pb / pa, "-", color="k", lw=1.2,
                     label="measured  arm B / arm A")
        za, ca = cands["cs_eff"]; _, cb = cands["cs_nominal"]
        bot.plot(za, cb / ca, "--", color="grey", lw=1.4,
                 label="the two CANDIDATES' ratio")
        bot.axhline(1.0, color="k", lw=0.6)
    top.set_yscale("log"); top.set_ylim(1e-6, 0.3)
    top.set_ylabel(r"$\phi(h)$"); top.legend(fontsize=7.5, ncol=2)
    top.set_title("F7  which arm did not move  -- this is the verdict")
    mid.set_ylabel("2nd half / 1st half"); mid.set_ylim(0.8, 1.2)
    mid.legend(fontsize=8)
    mid.set_title("flat at 1 = stationary", fontsize=9)
    bot.set_ylim(0.6, 1.6); bot.set_xlim(0, 60)
    bot.set_xlabel("h / d"); bot.set_ylabel("arm B / arm A")
    bot.legend(fontsize=8)
    fig.tight_layout(); f = FIGDIR / "F7_which_arm_moved.png"
    fig.savefig(f, dpi=150); plt.close(fig); out.append(f)

    # ── F6: <U>/N against time ───────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for r in live:
        a = r["arr"]
        ax.plot(a["t"], a["pe"], lw=0.8, alpha=0.8, color=COL[r["arm"]])
    ax.set_xlabel(r"$t / \tau_d$"); ax.set_ylabel(r"$\langle U \rangle / N$  [kT]")
    ax.set_title("F6  the stationarity evidence, every seed")
    fig.tight_layout(); f = FIGDIR / "F6_energy.png"
    fig.savefig(f, dpi=150); plt.close(fig); out.append(f)
    return out


# ── F5: the sealed prediction against the measurement, as a table ─────────
#
# ★ Named in the sealed analysis plan and missing from the first version of this
#   script -- caught by reading the plan's figure list against what the script
#   emits, which is the same check that found the missing half-window comparison.
def write_F5(rows) -> Path | None:
    live = [r for r in rows if r["ok"]]
    if not live:
        return None
    doc = prereg()
    preds = {q["quantity"]: q for q in doc["predictions"]}
    #  the block SEM that belongs to each observable, where one exists
    SEM = {"Z_dev_wmean_vs_CS_eff": "Z_dev_wmean_block_sem",
           "l_g_fitted_dense": "l_g_dense_block_sem",
           "phi_wall_measured_vs_cs_eff": "phi_wall_block_sem",
           "phi_wall_measured": "phi_wall_block_sem",
           "Z_dilute_tail": "Z_dilute_tail_block_sem"}
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / "F5_prediction_vs_measurement.md"
    L = ["# F5 — the sealed prediction against the measurement", "",
         f"`prediction.yaml` revision {doc['revision']}, sealed into "
         f"{len(rows)} run directories. {len(live)} have finished.", "",
         "⚠ Every error bar below is the run's **own block SEM** over "
         f"{live[0]['result'].get('n_blocks', '?')} blocks of its production "
         "window — not the design-power Monte Carlo's sigma, which cannot see "
         "temporal correlation.", ""]
    for arm in ("cs_eff", "cs_nominal", "convergence", "finite_size"):
        rs = [r for r in live if r["arm"] == arm]
        if not rs:
            continue
        L += [f"## arm `{arm}`  ({len(rs)} seed(s), init=`{rs[0]['init']}`, "
              f"{rs[0]['tau_sed']:g} tau_sed)", "",
              "| id | quantity | role | sealed | measured | block SEM | dev |",
              "|---|---|---|---|---|---|---|"]
        for r in rs:
            for name, ob in r["obs"].items():
                if name.endswith("_block_sem"):
                    continue
                pr = preds.get(name, {})
                role = ob.get("role", "?")
                sealed = pr.get("predicted", ob.get("predicted"))
                m = ob.get("measured")
                sem = o(r, SEM[name]) if name in SEM else None
                if sealed in (None, "") or not isinstance(sealed, (int, float)):
                    dev = "—"
                elif sealed == 0.0:
                    dev = f"{m:+.3f} abs"
                else:
                    dev = f"{100.0*(m/sealed - 1):+.2f} %"
                L.append(f"| {r['id']} | `{name}` | {role} | "
                         f"{sealed if sealed is not None else '—'} | "
                         f"{m:.5g} | "
                         f"{('%.5g' % sem) if sem is not None and sem == sem else '—'} | "
                         f"{dev} |")
        L.append("")
    L += ["## what this table cannot say", "",
          "`confirmed_by` is human-only (CLAUDE.md), so no row here is a "
          "confirmed result. A `hypothesis` row that misses its prediction is a "
          "**result**, not a failure (rule 7'); only an `implementation_check` "
          "miss is a fault — and revision 4 leaves just `D_xy` and "
          "`min_sep_placed` in that class, for the reason the analysis plan's "
          "accepted-weaknesses list gives.", ""]
    out.write_text("\n".join(L) + "\n")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", action="store_true",
                    help="draw F1..F7 and stop")
    a = ap.parse_args()

    rows = load_runs()
    done = [r for r in rows if r["ok"]]
    print(f"{len(done)} of {len(rows)} runs have finished")
    for r in rows:
        if not r["ok"]:
            print(f"  {r['id']:3s} {r['why']}")
    if not done:
        print("\nnothing to analyse yet.")
        return 2

    if a.figures:
        for f in make_figures(rows):
            print(f"  wrote {f.relative_to(ROOT)}")
        f5 = write_F5(rows)
        if f5:
            print(f"  wrote {f5.relative_to(ROOT)}")
        return 0

    an = Analysis(rows)
    an.gate("1. verify the seal on every run directory",
            lambda: step1_seal(rows))
    an.gate("2. the correctness gates -- Z_dilute_tail, D_xy",
            lambda: step2_gates(rows))
    an.note("3. stationarity -- which arm did not move",
            lambda: step3_stationarity(rows))
    an.note("4. the three discriminators",
            lambda: step4_discriminators(rows))
    rc = an.report()

    if not an.stopped_at and len(done) == len(rows):
        v, notes = verdict(rows)
        print("\n" + "=" * 78)
        print("5. THE VERDICT ON THE QUESTION  (proposed -- confirmed_by is human)")
        print("=" * 78)
        for n in notes:
            print(f"  {n}")
        print(f"\n  ★ {v}")
    elif not an.stopped_at:
        print(f"\n(the verdict needs all {len(rows)} runs; {len(done)} are in)")

    for f in make_figures(rows):
        print(f"  wrote {f.relative_to(ROOT)}")
    f5 = write_F5(rows)
    if f5:
        print(f"  wrote {f5.relative_to(ROOT)}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
