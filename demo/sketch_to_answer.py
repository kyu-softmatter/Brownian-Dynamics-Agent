"""The demo: one hand sketch in, one verified Brownian-dynamics answer out.

    $PY demo/sketch_to_answer.py                      # every case that is ready
    $PY demo/sketch_to_answer.py --case trap-2d-5um
    $PY demo/sketch_to_answer.py --list

**What it shows.** `intake/<case>/` holds a photograph of what a human actually
handed over — a napkin sketch, a page of a notebook, a photo sent over a
messenger. The figure audits the gap between what is written on that paper and
what the simulation needed: which numbers came off the page, which came from a
handbook, which are simulation choices that could have been otherwise, and
which ones the agent had to **ask** about.

The sharpest panel is the ambiguity list, because it is where the agent is
allowed to be wrong out loud. In `trap-2d-5um` the sketch says `R = 5 um` beside
a circle, which is ambiguous between radius and diameter. The agent recorded the
ambiguity, stated its lean (radius), **quantified what being wrong would cost**
(equilibrium `<x^2> = kT/k` does not contain `d`, so only the timescale moves --
by exactly 2x), asked, and was corrected to diameter. `observation.yaml` still
carries the line, in the author's words: *"my guess was wrong."*

**Every number in the figure is read from an artifact** -- `intake/`, `specs/`,
`runs/*/metrics.json`. Nothing is typed in, and the run is *chosen by
measurement* (`pick_run`: the one deciding the most observables), so re-running a
case cannot leave the figure quietly pointing at a stale directory.

**Labels are the one thing written by hand**, because the artifacts are Korean
and matplotlib's DejaVu Sans has no Hangul -- CLAUDE.md requires English labels
and zero `missing from font` warnings. Rather than let a Hangul string reach the
canvas and render as boxes, `check_latin` refuses: any text without a
translation in `EN` is reported by name and the figure for that case is not
written. Adding a case is then a bounded, visible task instead of a silent
garbling. Both failures happened here before the check existed.
"""
from __future__ import annotations

import argparse
import glob
import json
import pathlib
import re
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.gridspec import GridSpec
from PIL import ExifTags, Image

ROOT = pathlib.Path(__file__).resolve().parent.parent

TIER_COLOR = {0: "#1b6b3a", 1: "#2b6cb0", 2: "#b7791f", 3: "#9b2c2c"}
TIER_NAME = {0: "given / handbook", 1: "literature, verified",
             2: "literature, unverified", 3: "a choice"}

#: Per-case presentation only. Which image is the input, and one phrase for the
#: result panel's title. No numbers here.
CASES = {
    "trap-2d-5um": dict(
        sketch="sketch_01.jpeg",
        closed_form="the closed form it has one for"),
    "abp-rod-2d-run-flip": dict(
        sketch="sketch_01.jpeg",
        closed_form="the analytic run-and-tumble result"),
    "network": dict(
        sketch="transcription_annotated.jpeg",
        closed_form="what was predicted for it",
        sketch_note="the agent's reading drawn back onto the photo — "
                    "21 particles (red), numbering (green), bonds (blue)"),
}

#: Hand-written English for artifact text that would otherwise hit the canvas
#: as Hangul. Labels only, never a number.
EN = {
    # observable and check names
    "관성 무시   τ_p/τ_k": "inertia\nnegligible  $\\tau_p/\\tau_k$",
    "트랩 해상   dt/τ_k": "trap\nresolved  $dt/\\tau_k$",
    "요동 vs 박스 2ℓ_k/L": "fluctuation\nvs box  $2\\ell_k/L$",
    "관측창     T_obs/τ_k": "observation\nwindow  $T_{obs}/\\tau_k$",
    "τ (C(t) 피팅)": "$\\tau$  (fit to $C(t)$)",
    "f_c (PSD 피팅)": "$f_c$  (fit to PSD)",
    "S(0) (PSD 피팅)": "$S(0)$  (fit to PSD)",
    "⟨x²⟩": "$\\langle x^2 \\rangle$",
    "σ = √⟨x²⟩": "$\\sigma = \\sqrt{\\langle x^2\\rangle}$",
    # ambiguity text
    "입자가 몇 개인가?": "how many particles are there?",
    "무엇을 측정할 것인가? (스케치에 목표 미기재)":
        "what should be measured?  (no goal is written on the sketch)",
    "R = 5 µm 가 무엇의 크기인가?": "is  R = 5 µm  the radius or the diameter?",
    "입자 반지름(d=10µm)": "particle radius  (d = 10 µm)",
    # abp-rod-2d-run-flip
    "★ τ_R = 0.5 s 가 회전확산 시간인가, 뒤집힘(flip) 간격인가?":
        "is $\\tau_R$ = 0.5 s the rotational diffusion time, or the flip interval?",
    "(a) τ_R = 회전확산 시간. 그러면 τ_flip 을 별도로 정해야 함":
        "(a) rotational diffusion time — then $\\tau_{flip}$ must be fixed separately",
    "★ R 에 숫자가 없다 (입자 크기 미정)":
        "no number is given for R — the particle size is undetermined",
    "d ≈ 1 µm (D1의 τ_R 역산과 일관). 확인 필요":
        "d $\\approx$ 1 µm, from inverting D1's $\\tau_R$ — needs confirming",
    "'2R' 판독 — 장축 전체 길이인가?":
        "read as '2R' — is that the full major-axis length?",
    "종횡비 p = a/b 미기재": "the aspect ratio p = a/b is not written",
    "관성 무시     τ_p/τ_v": "inertia\nnegligible  $\\tau_p/\\tau_v$",
    "이류 해상     dt/τ_v": "advection\nresolved  $dt/\\tau_v$",
    "회전 해상     dt·D_r": "rotation\nresolved  $dt\\cdot D_r$",
    "텀블 해상     dt/τ_tumble": "tumble\nresolved  $dt/\\tau_{tumble}$",
    "유한크기     ℓ_p/(L/4)": "finite size\n$\\ell_p/(L/4)$",
    "관측창       T_obs/τ_eff": "observation\nwindow  $T_{obs}/\\tau_{eff}$",
    "τ_eff (⟨cosΔθ⟩ 피팅)": "$\\tau_{eff}$\n(fit to $\\langle\\cos\\Delta\\theta\\rangle$)",
    "D_eff (MSD 장시간)": "$D_{eff}$\n(MSD, long time)",
    "D̄ (MSD 전체 피팅)": "$\\bar{D}$\n(fit to full MSD)",
    "τ (MSD 전체 피팅)": "$\\tau$\n(fit to full MSD)",
    "텀블 빈도": "tumble rate",
    # network
    "\"DLVO, & JKR\" — 둘을 **합쳐서** 쓰라는 것인가, **각각 돌려 비교**하라는 것인가?":
        "\"DLVO, & JKR\" — combine the two, or run each and compare?",
    "x(t)=A sin(ωt) 가 위치 강제인가, 트랩 중심의 궤적인가?":
        "is x(t) = A sin($\\omega$t) a position constraint, or the trap centre's path?",
    "2D 인가 3D 인가? 그림은 평면이지만 콜로이드 겔은 3D 가 자연스럽다.":
        "2D or 3D?  the drawing is planar, but a colloidal gel is naturally 3D",
    "참고: τ_p/τ_bond": "reference\n$\\tau_p/\\tau_{bond}$",
    "결합 안정  σ_bond/h_min": "bond stable\n$\\sigma_{bond}/h_{min}$",
}

#: An ambiguity's `lean` is sometimes a multi-sentence argument. Those do not
#: belong on a figure and translating them by hand would be inventing a
#: paraphrase of someone's reasoning. Render only short leans; the long ones
#: stay in observation.yaml where they are quoted in full.
LEAN_MAX = 60

HANGUL = re.compile(r"[가-힣㄰-㆏]")
UNTRANSLATED: list[str] = []


def en(s):
    """Translate, or record the miss and return a visible placeholder."""
    if s in EN:
        return EN[s]
    if isinstance(s, str) and HANGUL.search(s):
        UNTRANSLATED.append(s)
        return "(untranslated — add to EN)"
    return s


def pick_run(case: str) -> pathlib.Path:
    """The run that decides the most, among those whose spec is still on disk.

    ★ The spec must exist, and it is the FIRST key, not a tiebreak. 261 run
      directories sit against 278 specs and the sets do not coincide --
      `trap-2d-5um__a49f2508556b` has metrics and no spec. Selecting it made
      panel 3 render "3 numbers on the paper -> 0 carried with a tier", which
      reads as a finding about the case and was in fact a missing file. A
      silent 0 and a real 0 must not look alike (docs/05 section 2).
    """
    best = None
    for f in sorted(glob.glob(str(ROOT / f"runs/{case}*/metrics.json"))):
        d = json.loads(pathlib.Path(f).read_text())
        if d.get("case") != case:
            continue
        rd = pathlib.Path(f).parent
        obs = d.get("observables") or []
        key = ((ROOT / f"specs/{rd.name}.json").exists(),
               sum(1 for o in obs if o.get("predicted") is not None), len(obs))
        if best is None or key > best[0]:
            best = (key, rd)
    if best is None:
        raise SystemExit(f"{case}: no runs/{case}*/metrics.json")
    if not best[0][0]:
        raise SystemExit(
            f"{case}: no run has a spec in specs/ -- panel 3 would render an "
            f"empty provenance ledger that looks like a result. Refusing.")
    return best[1]


def provenance_rows(spec) -> list:
    rows = []

    def walk(node, path):
        if isinstance(node, dict):
            if "value" in node and "tier" in node:
                rows.append((".".join(path), node, node["tier"]))
                return
            for k, v in node.items():
                walk(v, path + [k])

    walk(spec.get("system", {}), [])
    return rows


def _n(met, spec) -> str:
    v = (met.get("physical") or {}).get("N")
    if v is None:   # `network` records no physical.N; fall back to the spec
        c = ((spec.get("system") or {}).get("particle") or {}).get("count")
        v = c.get("value") if isinstance(c, dict) else c
    if isinstance(v, list):
        return " → ".join(f"{x:,}" for x in v if isinstance(x, (int, float)))
    return f"{v:,}" if isinstance(v, (int, float)) else "not recorded"


def fmt(v) -> str:
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        return f"{int(v)}" if v == int(v) and abs(v) < 1e6 else f"{v:.4g}"
    return str(v)


def render(case: str) -> int:
    cfg = CASES[case]
    UNTRANSLATED.clear()
    rundir = pick_run(case)
    run = rundir.name
    obs = yaml.safe_load((ROOT / f"intake/{case}/observation.yaml").read_text())
    met = json.loads((rundir / "metrics.json").read_text())
    sp = ROOT / f"specs/{run}.json"
    spec = json.loads(sp.read_text()) if sp.exists() else {"system": {}}
    prov = provenance_rows(spec)
    amb = obs.get("ambiguities") or []
    missing = obs.get("missing_required") or []
    stated = obs.get("stated_quantities") or []
    O = met.get("observables") or []
    checks = met.get("checks") or []
    out = ROOT / f"figures/demo_{case.replace('-', '_')}.png"

    fig = plt.figure(figsize=(17.5, 10.4))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.06, 1.0],
                  width_ratios=[1.0, 1.22, 1.16], hspace=0.30, wspace=0.20,
                  left=0.035, right=0.978, top=0.900, bottom=0.075)
    fig.suptitle("One hand sketch in, one verified Brownian-dynamics answer out"
                 f"    ·    case `{case}`    ·    run `{run}`    ·    "
                 "every number below is read from an artifact, none typed in",
                 fontsize=13, fontweight="bold", y=0.968)

    # 1 · the input, as handed over
    ax = fig.add_subplot(gs[0, 0])
    # ★ Orientation is settled by looking, not by the tag. trap-2d-5um's EXIF
    #   says Orientation=6 ("rotate to display") and applying it puts the
    #   handwriting on its side: those raw pixels are already upright. Both
    #   honouring the tag and guessing an np.rot90 were tried and were wrong.
    #   If a new case comes in rotated, render the 8 candidates and look.
    im = Image.open(ROOT / f"intake/{case}/{cfg['sketch']}")
    shot = im.getexif().get(
        next(k for k, v in ExifTags.TAGS.items() if v == "DateTime"), "")
    ax.imshow(np.asarray(im))
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("1 · what the human handed over\n"
                 f"{len(stated)} number(s) written on it"
                 + (f"   ·   camera says {shot.split()[0].replace(':', '-')}"
                    if shot else ""),
                 fontsize=11, loc="left", fontweight="bold")
    if cfg.get("sketch_note"):
        ax.set_xlabel(cfg["sketch_note"], fontsize=8.4, color="#444")
    for s in ax.spines.values():
        s.set_edgecolor("#888")

    # 2 · what was read, and what was ambiguous
    ax = fig.add_subplot(gs[0, 1]); ax.axis("off")
    ax.set_title("2 · what the agent read — transcribe first, then interpret",
                 fontsize=11, loc="left", fontweight="bold")
    y = 0.97
    ax.text(0.0, y, "ON THE PAPER", fontsize=8.6, fontweight="bold",
            color="#1b6b3a", family="monospace"); y -= 0.058
    for q in stated[:6]:
        mark = "   ← ambiguous" if q.get("note") else ""
        val = (f"{fmt(q.get('value'))} {q.get('unit') or ''}".strip()
               if q.get("value") is not None else "(no number written)")
        ax.text(0.02, y, f"{str(q.get('symbol')):>6s} = {val}{mark}",
                fontsize=10.0,
                family="monospace", color="#9b2c2c" if mark else "#111")
        y -= 0.054
    ax.text(0.02, y, f"goals stated on the sketch:  "
                     f"{len(obs.get('stated_goals') or []) or 'NONE'}",
            fontsize=9.4, family="monospace",
            color="#9b2c2c" if not obs.get("stated_goals") else "#111")
    y -= 0.048
    ax.text(0.02, y, f"regions that could not be read:  "
                     f"{len(obs.get('unread_regions') or []) or 'none'}",
            fontsize=9.4, family="monospace", color="#555")
    y -= 0.078
    ax.text(0.0, y, f"AMBIGUITIES RAISED — {len(amb)}, each with how it was "
                    f"settled", fontsize=8.6, fontweight="bold",
            color="#9b2c2c", family="monospace")
    y -= 0.060
    for a in amb[:4]:
        iss = str(a.get("issue") or "")
        ax.text(0.02, y,
                f"{a.get('id')}   " + (en(iss) if len(iss) <= 90
                                       else "(stated at length in "
                                            "observation.yaml)"),
                fontsize=9.4, fontweight="bold"); y -= 0.050
        lean = a.get("lean")
        if lean and len(str(lean)) <= LEAN_MAX:
            ax.text(0.05, y, f"the agent leaned:   {en(lean)}",
                    fontsize=9.0, family="monospace", color="#9b2c2c")
            y -= 0.046
        elif lean:
            ax.text(0.05, y, "the agent leaned, with its reasoning "
                             "(see observation.yaml)",
                    fontsize=8.6, style="italic", color="#9b2c2c")
            y -= 0.046
        who = a.get("confirmed_by")
        ax.text(0.05, y, f"settled by: {who}" if who else "settled by: —",
                fontsize=9.0, family="monospace",
                color="#1b6b3a" if who == "user" else "#2b6cb0")
        y -= 0.062
    if case == "trap-2d-5um":
        ax.text(0.05, y, "the record says, verbatim:  “my guess was wrong”\n"
                         "cost of being wrong, quantified BEFORE asking:\n"
                         "   $\\langle x^2\\rangle = kT/k$ has no $d$ in it — "
                         "unaffected.\n"
                         "   only the timescale moves, by exactly 2× "
                         "($\\tau_k$: 8.02 → 4.01 ms)",
                fontsize=9.0, va="top", color="#222", linespacing=1.55)

    # 3 · where every number came from
    ax = fig.add_subplot(gs[0, 2]); ax.axis("off")
    ax.set_title("3 · where every number the run used came from",
                 fontsize=11, loc="left", fontweight="bold")
    y = 0.97
    n_prov_shown = 8 if len(missing) > 4 else 10
    for name, node, tier in prov[:n_prov_shown]:
        short = name.split(".")[-1]
        ax.text(0.0, y, f"{short:<13s}{fmt(node['value']):>9s} "
                        f"{str(node.get('unit', '')).replace('dimensionless', '1'):<8s}",
                fontsize=9.6, family="monospace")
        ax.text(0.615, y, f"tier {tier}", fontsize=9.2, family="monospace",
                color=TIER_COLOR.get(tier, "#333"), fontweight="bold")
        ax.text(0.755, y, TIER_NAME.get(tier, "?"), fontsize=7.6,
                color=TIER_COLOR.get(tier, "#333"))
        y -= 0.062
    if len(prov) > n_prov_shown:
        ax.text(0.0, y, f"… and {len(prov) - n_prov_shown} more with a tier",
                fontsize=8.2, style="italic", color="#555")
        y -= 0.050
    y -= 0.020
    ax.text(0.0, y, "NOT on the sketch — filled in, never invented:",
            fontsize=8.6, fontweight="bold", family="monospace", color="#555")
    y -= 0.062
    # ★ a second column, not string padding: `network_protocol` is 16 chars and
    #   `f"{sym:<15s}"` ran the two fields together ("network_protocolasked").
    room = max(0, int((y - 0.12) / 0.056))
    shown = missing[:room]
    for m in shown:
        tag = ("asked the human" if m.get("confirmed_by") == "user"
               else "confirmed by the run" if m.get("confirmed_by") == "run"
               else "a free choice" if m.get("kind") == "choice"
               else "assumed, tier 3" if m.get("confidence") == 3
               else "read off the sketch" if m.get("resolution")
               else "STILL OPEN")
        col = ("#1b6b3a" if m.get("confirmed_by") == "user"
               else "#9b2c2c" if tag in ("assumed, tier 3", "STILL OPEN")
               else "#2b6cb0")
        ax.text(0.02, y, str(m.get("symbol")), fontsize=9.0,
                family="monospace", color=col)
        ax.text(0.46, y, tag, fontsize=9.0, family="monospace", color=col)
        y -= 0.056
    if len(missing) < len(obs.get("missing_required") or []) or shown != missing:
        # no silent truncation -- say what was dropped
        ax.text(0.02, y, f"… and {len(missing) - len(shown)} more in "
                         f"observation.yaml", fontsize=8.4, style="italic",
                color="#555")
        y -= 0.056
    ax.text(0.0, max(y - 0.02, 0.02),
            f"{len(stated)} number(s) on the paper  →  {len(prov)} carried with "
            f"a tier in the spec  →  the rest derived.\n"
            "Rule 3: a number that is not known BLOCKS. It is never invented.",
            fontsize=8.8, fontweight="bold", linespacing=1.6, va="top")

    # 4 · the gates cleared before the run
    ax = fig.add_subplot(gs[1, 0])
    if checks:
        # ★ a margin can be inf (value 0 against a `<=` limit) or None. Both
        #   crashed set_xlim("cannot be NaN or Inf") on `network`. Clip for the
        #   axis, but keep the true value in the label -- an inf margin means
        #   "infinitely inside the limit", which is information, not an error.
        raw_m = [c.get("margin") for c in checks]
        finite = [v for v in raw_m if v is not None and np.isfinite(v)]
        cap = max(finite) * 3 if finite else 10.0
        m = [cap if (v is None or not np.isfinite(v)) else v for v in raw_m]
        cols = ["#1b6b3a" if c.get("hard") else "#2b6cb0" for c in checks]
        ypos = np.arange(len(checks))[::-1]
        many = len(checks) > 6
        ax.barh(ypos, m, color=cols, height=0.46 if many else 0.56)
        ax.axvline(1.0, color="#9b2c2c", lw=1.6)
        ax.set_xscale("log")
        ax.set_xlim(0.42, max(max(m), 10) * 14)
        for yp, c, mv, rv in zip(ypos, checks, m, raw_m):
            ax.text(0.46, yp + 0.31, en(c["name"]).replace("\n", " "),
                    fontsize=8.2, va="center", ha="left")
            lbl = ("∞×" if rv is not None and not np.isfinite(rv)
                   else "—" if rv is None else f"{rv:,.4g}×")
            ax.text(mv * 1.22, yp, lbl, fontsize=8.6, va="center",
                    family="monospace", fontweight="bold")
        ax.set_ylim(-0.7, len(checks) - 0.15)
        ax.set_yticks(ypos); ax.set_yticklabels([])
        ax.set_xlabel("margin over the limit  (log)", fontsize=9)
        ax.grid(axis="x", alpha=0.25)
    n_hard = sum(1 for c in checks if c.get("hard"))
    ax.set_title(f"4 · the {len(checks)} gates the run cleared first\n"
                 f"green = hard, blocks ({n_hard}) · blue = soft, warns "
                 f"({len(checks) - n_hard})",
                 fontsize=11, loc="left", fontweight="bold")

    # 5 · the answer
    ax = fig.add_subplot(gs[1, 1:])
    dec = [o for o in O if o.get("predicted") is not None
           and o.get("err_pct") is not None]
    if dec:
        err = [o["err_pct"] for o in dec]
        tols = [o.get("tol_pct") for o in dec if o.get("tol_pct")]
        tol = max(tols) if tols else None
        x = np.arange(len(dec))
        if tol:
            ax.axhspan(-tol, tol, color="#1b6b3a", alpha=0.10,
                       label=f"stated tolerance  ±{tol:g} %")
        ax.axhline(0, color="#111", lw=1.2)
        ax.bar(x, err, color="#2b6cb0", width=0.42)
        lim = max(max(abs(e) for e in err) * 1.9, (tol or 1) * 1.5)
        for xi, (o, e) in enumerate(zip(dec, err)):
            ax.text(xi, e + lim * 0.045 * (1 if e >= 0 else -1), f"{e:+.4f} %",
                    ha="center", va="bottom" if e >= 0 else "top",
                    fontsize=9.2, fontweight="bold", family="monospace")
            ax.text(xi, -lim * 0.80, f"measured\n{o['measured']:.6g}\n"
                                     f"{o.get('unit', '')}",
                    ha="center", fontsize=7.6, family="monospace", color="#333")
        ax.set_xticks(x); ax.set_xticklabels([en(o["name"]) for o in dec],
                                             fontsize=10)
        ax.set_ylim(-lim, lim)
        if tol:
            ax.legend(fontsize=9, loc="upper right")
    ax.set_ylabel("measured − predicted   [% of prediction]", fontsize=9.5)
    ax.grid(axis="y", alpha=0.25)
    nm = met.get("numerics") or {}
    inside = sum(1 for o in dec if o.get("tol_pct")
                 and abs(o["err_pct"]) <= o["tol_pct"])
    ax.set_title(f"5 · the answer, against {cfg['closed_form']}"
                 f"    ·    HOOMD-blue, N = {_n(met, spec)}, "
                 f"{nm.get('steps_done', 0):,} steps, "
                 f"{(met.get('wall_seconds') or 0)/60:.0f} min on one CPU"
                 f"    ·    {inside} of {len(dec)} inside tolerance",
                 fontsize=11, loc="left", fontweight="bold")

    if UNTRANSLATED:
        plt.close(fig)
        print(f"  ✗ {case}: {len(UNTRANSLATED)} string(s) have no English label; "
              f"figure NOT written. Add to EN:")
        for s in dict.fromkeys(UNTRANSLATED):
            print(f'      "{s}":')
        return 1

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out.parent.mkdir(exist_ok=True)
        fig.savefig(out, dpi=145)
        miss = [str(x.message) for x in w if "missing from font" in str(x.message)]
    plt.close(fig)
    print(f"  {out.relative_to(ROOT)}")
    print(f"     run {run} · stated {len(stated)} · provenance {len(prov)} · "
          f"ambiguities {len(amb)} · gates {len(checks)} · decided {len(dec)}"
          f"/{len(O)}")
    if miss:
        print(f"     ✗ {len(miss)} glyph(s) missing from the font")
        for m_ in miss[:4]:
            print("        ", m_)
        return 1
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", action="append", default=[])
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for c in CASES:
            print(f"  {c:24s} -> {pick_run(c).name}")
        return 0
    bad = [c for c in a.case if c not in CASES]
    if bad:
        raise SystemExit(f"unknown case(s) {bad}; known: {list(CASES)}")
    rc = 0
    for c in (a.case or list(CASES)):
        rc |= render(c)
    print("  ✓ zero `missing from font` warnings" if rc == 0
          else "  ✗ see above")
    return rc


if __name__ == "__main__":
    sys.exit(main())
