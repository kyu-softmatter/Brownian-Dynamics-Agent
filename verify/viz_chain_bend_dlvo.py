"""Check the first `chain-bend-2d-dlvo` run -- graphs plus animation.

**Do not report a result in prose alone** (working practice). This is the visual
check of whether the run supports hypothesis G1: a straight chain with purely
central-force DLVO, equilibrated at its natural length, has exactly zero linear
bending stiffness.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/viz_chain_bend_dlvo.py [run_id]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "scratch" / "_viz"
OUT.mkdir(parents=True, exist_ok=True)

for _f in ("Arial Unicode MS", "Apple SD Gothic Neo", "AppleGothic", "NanumGothic"):
    if _f in {f.name for f in matplotlib.font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = _f
        break
plt.rcParams["axes.unicode_minus"] = False

C_MEAS, C_ZERO, C_GOOD, C_BAD = "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e"

RUN_ID = sys.argv[1] if len(sys.argv) > 1 else "chain-bend-2d-dlvo__n9-w1000-a100__cbe816dbec24"
RUN = ROOT / "runs" / RUN_ID


def load():
    m = json.loads((RUN / "metrics.json").read_text())
    npz = np.load(RUN / "observables.npz", allow_pickle=True)
    return m, npz


def make_figure(m, npz):
    n = int(m["physical"]["n_beads"])
    kp = next(o for o in m["observables"] if o["name"] == "K_prime")
    shape_mean = np.array(m["result"]["shape_profile_mean"])
    t, yb, yg = npz["t"], npz["y_bead"], npz["y_ghost"]
    shape_y = npz["shape_y"]

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    fig.suptitle(f"chain-bend-2d-dlvo -- first run (n={n}, DLVO central force only, "
                 f"no bending term)  "
                 f"run_id={RUN_ID}", fontsize=11, y=1.03)

    # (1) mean shape profile -- smooth curvature (parabola-like) vs a localized
    #     kink (triangle-like)
    a = ax[0]
    idx = np.arange(n)
    a.plot(idx, shape_mean * 1e3, "o-", color=C_MEAS,
           label="measured <y_i> (time-averaged)")
    # If the curvature were smooth and beam-like a quadratic would fit well --
    # overlaid as a reference
    coef = np.polyfit(idx, shape_mean, 2)
    a.plot(idx, np.polyval(coef, idx) * 1e3, "--", color=C_GOOD, alpha=0.7,
           label="least-squares quadratic (a smooth shape would sit near this)")
    a.set_xlabel("bead index i")
    a.set_ylabel("<y_i>  [d, ×1e-3]")
    a.set_title(f"shape profile (shape_localization={m['result']['shape_localization']:.3f}"
                f"  .  1=smooth, >>1=localized kink)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # (2) drive vs response time series (the last few cycles)
    a = ax[1]
    T_period = 2 * np.pi / m["physical"]["omega_star"] / (t[1] - t[0]) * (t[1] - t[0])
    mask = t > t.max() - 6 * (2 * np.pi / (m["physical"]["omega_star"] / m["numerics"].get("dt_star", 1)))
    # Play it safe: just the last 400 samples
    sl = slice(max(0, len(t) - 400), len(t))
    a.plot(t[sl] * 1e3, yg[sl], color="#888", lw=1, label="drive (ghost, y_ghost)")
    a.plot(t[sl] * 1e3, yb[sl], color=C_MEAS, lw=1.3,
           label="centre-bead response (y_bead)")
    a.set_xlabel("t  [ms]")
    a.set_ylabel("y  [d]")
    a.set_title("drive vs centre-bead response (last 400 samples)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # (3) the G1 check: measured K' vs the prediction of 0, in units of sigma
    a = ax[2]
    meas, sig = kp["measured"], kp["sigma"]
    z = kp["err_sigma"]
    a.errorbar([0], [meas], yerr=[sig], fmt="o", color=C_MEAS, capsize=6, ms=9,
               label=f"measured K' = {meas:.1f} ± {sig:.1f} kT/d^2")
    a.axhline(0, color=C_ZERO, ls="--",
              label="prediction (G1: linear bending stiffness = 0)")
    ok = abs(z) < (kp.get("tol_sigma") or 3.0)
    a.set_title(f"G1 verdict: {'✓ consistent' if ok else '✗ inconsistent'} "
                f"({z:+.2f} sigma, criterion ±3 sigma)",
               color=C_GOOD if ok else C_BAD)
    a.set_xlim(-1, 1)
    a.set_xticks([])
    a.set_ylabel("K′  [kT/d²]")
    a.legend(fontsize=8, loc="upper right")
    a.grid(alpha=0.3)

    fig.tight_layout()
    p = OUT / "chain_bend_dlvo_results.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


def make_animation(npz):
    """How the chain shape evolves in time -- to see whether a triangular buckle
    actually appears.

    ★ Taken from the real production data (the first 300 of 2000 samples, ~3 cycles)
    -- this is this run's actual trajectory, not a separate cheap run.
    """
    shape_y = npz["shape_y"][:300]
    t = npz["t"][:300]
    n = shape_y.shape[1]
    idx = np.arange(n)

    fig, ax = plt.subplots(figsize=(5, 4))
    ymax = np.abs(shape_y).max() * 1.15 + 1e-6
    line, = ax.plot([], [], "o-", color=C_MEAS, lw=2, ms=6)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(-ymax, ymax)
    ax.set_xlabel("bead index i")
    ax.set_ylabel("y  [d]")
    title = ax.set_title("")
    ax.grid(alpha=0.3)

    def update(fr):
        line.set_data(idx, shape_y[fr])
        title.set_text(f"chain-bend-2d-dlvo  t={t[fr]*1e3:.3f} ms  "
                       f"(production trajectory)")
        return line, title

    anim = FuncAnimation(fig, update, frames=len(t), interval=60, blit=False)
    p = OUT / "chain_bend_dlvo_motion.gif"
    anim.save(p, writer=PillowWriter(fps=16))
    plt.close(fig)
    return p


def main():
    m, npz = load()
    png = make_figure(m, npz)
    gif = make_animation(npz)
    print(f"graph: {png.relative_to(ROOT)}")
    print(f"animation: {gif.relative_to(ROOT)}  (this run's actual production "
          f"trajectory, first ~3 cycles)")
    missing = {w for w in ("missing from font",) if False}
    return 0


if __name__ == "__main__":
    sys.exit(main())
