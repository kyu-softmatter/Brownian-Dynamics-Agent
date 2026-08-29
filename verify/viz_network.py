#!/usr/bin/env python
"""Turn the stage-1 (compression gelation) `network` results into figures --
**show** the network that was built.

CLAUDE.md: "show results as graphs and animations. Do not report in prose and
tables alone."
⚠️ All labels are in **English** (matplotlib's default font has no Hangul, and the
   Hangul fonts have no minus sign or y-hat).

What it makes:
  network_structure.png  6 panels -- 3D network, per-axis projections, time series,
                         g(r), collapse/minimum distance, topology
  network_compare.png    compression-rate comparison (A8: rate independence) --
                         only when both runs exist
  network_compress.mp4   animation of the compression (traj_A.gsd)

Run:
  $PY scratch/viz_network.py                       # auto-discovers runs/network__*
  $PY scratch/viz_network.py --run <run_id>
  $PY scratch/viz_network.py --no-anim
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402
import numpy as np                                                # noqa: E402
from matplotlib.collections import LineCollection                 # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection           # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "runs" / "network__FIGS"
ELL_STAR = 1.0075925903599383        # DLVO secondary-minimum bond distance (from the ledger)
H_MIN_STAR = 0.0075925903599383


def load(run_dir: Path) -> dict:
    z = np.load(run_dir / "observables.npz")
    m = json.loads((run_dir / "metrics.json").read_text())
    d = {k: z[k] for k in z.files}
    d["_metrics"] = m
    d["_id"] = run_dir.name
    d["_stage_tau"] = float(m["numerics"].get("stage_tau", float("nan"))) \
        if isinstance(m.get("numerics"), dict) else float("nan")
    return d


def phase_edges(d: dict):
    """Find the aggregation / compression / post-relaxation boundaries from phi(t).

    phi only changes during compression.
    """
    phi = d["series_phi"]
    t = d["series_t"]
    moving = np.abs(np.diff(phi)) > 1e-12
    if not moving.any():
        return t[0], t[-1]
    i0 = int(np.argmax(moving))
    i1 = int(len(moving) - np.argmax(moving[::-1]))
    return float(t[i0]), float(t[min(i1, len(t) - 1)])


# ══════════════════════════════════════════════════════════════════════
def panel_network_3d(ax, pos, pairs, L, title):
    deg = np.zeros(len(pos), dtype=int)
    for i, j in pairs:
        deg[i] += 1
        deg[j] += 1
    segs = []
    for i, j in pairs:
        a, b = pos[i], pos[j]
        if np.abs(b - a).max() > L / 2:          # bonds crossing the PBC are not drawn
            continue
        segs.append([a, b])
    if segs:
        ax.add_collection3d(Line3DCollection(segs, colors="0.45", linewidths=0.7,
                                             alpha=0.7))
    s = ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=deg, cmap="viridis",
                   s=14, vmin=0, vmax=max(6, deg.max()), depthshade=True,
                   edgecolors="none")
    ax.set_xlim(-L / 2, L / 2); ax.set_ylim(-L / 2, L / 2); ax.set_zlim(-L / 2, L / 2)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x / d", fontsize=8); ax.set_ylabel("y / d", fontsize=8)
    ax.set_zlabel("z / d", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(title, fontsize=9)
    return s


def panel_slab(ax, pos, pairs, L, thick=3.0):
    """z-slab projection -- shows the connectivity better than a 3D scatter."""
    sel = np.abs(pos[:, 2]) < thick / 2
    idx = np.where(sel)[0]
    keep = set(idx.tolist())
    segs = [[pos[i, :2], pos[j, :2]] for i, j in pairs
            if i in keep and j in keep and np.abs(pos[j] - pos[i]).max() < L / 2]
    if segs:
        ax.add_collection(LineCollection(segs, colors="0.35", linewidths=1.0, alpha=0.8))
    ax.scatter(pos[idx, 0], pos[idx, 1], s=26, c="#2a6ebb", edgecolors="w", linewidths=0.4)
    ax.set_xlim(-L / 2, L / 2); ax.set_ylim(-L / 2, L / 2); ax.set_aspect("equal")
    ax.set_xlabel("x / d"); ax.set_ylabel("y / d")
    ax.set_title(f"z-slab |z| < {thick/2:.1f} d   ({len(idx)} of {len(pos)} beads)",
                 fontsize=9)


def panel_series(ax, d, keys_labels, ylabel, logy=False):
    t = d["series_t"]
    tc0, tc1 = phase_edges(d)
    ax.axvspan(t[0], tc0, color="#eef3fa", zorder=0)
    ax.axvspan(tc0, tc1, color="#fdf0e6", zorder=0)
    ax.axvspan(tc1, t[-1], color="#eef7ee", zorder=0)
    for k, lab, c in keys_labels:
        ax.plot(t, d[f"series_{k}"], lw=1.6, label=lab, color=c)
    ax.set_xlabel(r"time / $\tau_B$")
    ax.set_ylabel(ylabel)
    if logy:
        ax.set_yscale("log")
    ax.legend(fontsize=7, framealpha=0.9)
    ax.grid(alpha=0.25)
    return tc0, tc1


def figure_structure(d: dict, out: Path):
    pos, pairs = d["final_positions"], d["final_pairs"].astype(int)
    L = float(d["_metrics"]["result"]["L_final"])
    n = len(pos)
    obs = {o["name"]: o for o in d["_metrics"]["observables"]}
    g = lambda k: obs.get(k, {}).get("measured", float("nan"))

    fig = plt.figure(figsize=(16.5, 10.0))
    gs = fig.add_gridspec(2, 3, hspace=0.30, wspace=0.28)

    ax = fig.add_subplot(gs[0, 0], projection="3d")
    s = panel_network_3d(ax, pos, pairs, L,
                         f"Final network  N={n}, $\\phi$={g('phi_final'):.3f}, "
                         f"L={L:.1f} d\n"
                         f"z={g('coordination_number'):.2f}, "
                         f"loops={int(g('independent_loops'))}, "
                         f"percolation={g('percolation'):.2f}")
    cb = fig.colorbar(s, ax=ax, shrink=0.62, pad=0.10)
    cb.set_label("coordination number", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    panel_slab(fig.add_subplot(gs[0, 1]), pos, pairs, L)

    ax = fig.add_subplot(gs[0, 2])
    tc0, tc1 = panel_series(ax, d, [("z", "coordination z", "#1f4e9c"),
                                    ("dangling", "dangling fraction", "#c0392b"),
                                    ("largest_cluster", "largest cluster / N", "#1e8449")],
                            "value")
    ax2 = ax.twinx()
    ax2.plot(d["series_t"], d["series_phi"], ls="--", lw=1.2, color="0.4")
    ax2.set_ylabel(r"$\phi$  (dashed)", fontsize=8)
    ax2.tick_params(labelsize=7)
    ax.set_title("Topology vs time  (z is still rising -> not yet converged)", fontsize=9)

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(d["rdf_fine_r"], d["rdf_fine_g"], lw=1.4, color="#1f4e9c",
            label="g(r), fine bins")
    ax.axvline(ELL_STAR, color="#c0392b", ls="--", lw=1.3,
               label=f"DLVO 2nd min  $\\ell$={ELL_STAR:.5f} d")
    ax.axvline(1.0, color="0.5", ls=":", lw=1.1, label="contact r = d  (h = 0)")
    pk = g("rdf_first_peak")
    ax.axvline(pk, color="#1e8449", ls="-.", lw=1.1,
               label=f"measured peak {pk:.5f} d")
    ax.set_xlim(0.97, 1.06)
    ax.set_xlabel("r / d"); ax.set_ylabel("g(r)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25)
    err = 100 * (pk / ELL_STAR - 1) if pk == pk else float("nan")
    ax.set_title(f"Bond length check (implementation_check)   error {err:+.3f}%",
                 fontsize=9)

    ax = fig.add_subplot(gs[1, 1])
    panel_series(ax, d, [("crushed", "crushed bond fraction", "#c0392b")], "fraction")
    ax3 = ax.twinx()
    ax3.plot(d["series_t"], d["series_min_sep"], lw=1.4, color="#1f4e9c")
    ax3.axhline(1.0, color="0.5", ls=":", lw=1.1)
    ax3.axhline(ELL_STAR, color="#c0392b", ls="--", lw=1.0)
    ax3.set_ylabel("min separation / d  (blue)", fontsize=8)
    ax3.tick_params(labelsize=7)
    ax.set_title("Crush / overlap check — no pair may sit inside the DLVO barrier\n"
                 rf"$\epsilon_{{max}}$={d['_metrics']['result']['eps_max']*100:.3f}% "
                 rf"< $\epsilon_{{crit}}$={d['_metrics']['result']['eps_crit']*100:.3f}%",
                 fontsize=9)

    ax = fig.add_subplot(gs[1, 2])
    panel_series(ax, d, [("loops", "independent loops", "#7d3c98"),
                         ("n_components", "connected components", "#d68910")],
                 "count", logy=True)
    ax4 = ax.twinx()
    ax4.plot(d["series_t"], d["series_d_f"], lw=1.4, color="#117864")
    ax4.axhline(1.8, color="0.5", ls=":", lw=1.0)
    ax4.set_ylabel(r"fractal dim $d_f$  (teal, dotted=DLCA 1.8)", fontsize=8)
    ax4.tick_params(labelsize=7)
    ax.set_title(f"Loops / components / $d_f$    final $d_f$={g('fractal_dimension'):.2f}",
                 fontsize=9)

    gen = ("sprout generator + DLVO relaxation" if "sprout" in d["_id"]
           else "compression gelation")
    fig.suptitle(f"network — stage 1, {gen}   {d['_id']}\n"
                 f"sketch drew a 21-bead TREE (0 loops) — schematic; "
                 f"the real topology emerges from the DLVO relaxation",
                 fontsize=11)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=135, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.relative_to(ROOT)}")


def figure_compare(runs: list, out: Path):
    """A8 -- compression-rate independence.

    Same seed, same aggregation; only the compression rate differs.
    """
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.2))
    keys = [("z", "coordination z"), ("loops", "independent loops"),
            ("largest_cluster", "largest cluster / N"), ("crushed", "crushed fraction")]
    for ax, (k, lab) in zip(axes, keys):
        for d in runs:
            st = d["_metrics"]["numerics"].get("stage_tau", float("nan"))
            ax.plot(d["series_phi"], d[f"series_{k}"], lw=1.6,
                    label=rf"$T_{{stage}}$ = {st:g} $\tau_B$")
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    # Tabulate the final values -- if it is rate-independent the two must agree
    txt = []
    for d in runs:
        o = {x["name"]: x.get("measured") for x in d["_metrics"]["observables"]}
        st = d["_metrics"]["numerics"].get("stage_tau", float("nan"))
        txt.append(f"T_stage={st:g}: z={o['coordination_number']:.3f}  "
                   f"loops={int(o['independent_loops'])}  "
                   f"perc={o['percolation']:.2f}  d_f={o['fractal_dimension']:.2f}")
    fig.suptitle("A8 — is the structure independent of compression rate?   "
                 "same seed, same aggregation, compression speed differs 10x\n"
                 + "     |     ".join(txt), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    fig.savefig(out, dpi=135, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out.relative_to(ROOT)}")


def animate(run_dir: Path, out: Path, stride: int = 1):
    """Animation of the compression -- z-slab projection. The box visibly shrinks."""
    try:
        import gsd.hoomd
        from matplotlib.animation import FFMpegWriter, PillowWriter
    except Exception as e:                                        # noqa: BLE001
        print(f"  skipping the animation: {e}")
        return
    p = run_dir / "traj_A.gsd"
    if not p.exists():
        print("  skipping the animation: no traj_A.gsd")
        return
    with gsd.hoomd.open(str(p), mode="r") as tr:
        frames = [(np.array(f.particles.position), float(f.configuration.box[0]))
                  for f in tr[::stride]]
    if not frames:
        return
    L0 = frames[0][1]
    fig, ax = plt.subplots(figsize=(6.2, 6.4))
    writer = (FFMpegWriter(fps=15) if matplotlib.animation.writers.is_available("ffmpeg")
              else PillowWriter(fps=12))
    if isinstance(writer, PillowWriter):
        out = out.with_suffix(".gif")
    with writer.saving(fig, str(out), dpi=110):
        for pos, L in frames:
            ax.clear()
            sel = np.abs(pos[:, 2]) < 1.5
            ax.scatter(pos[sel, 0], pos[sel, 1], s=30, c="#2a6ebb",
                       edgecolors="w", linewidths=0.4)
            ax.add_patch(plt.Rectangle((-L / 2, -L / 2), L, L, fill=False,
                                       ec="#c0392b", lw=1.6))
            ax.set_xlim(-L0 / 2 * 1.05, L0 / 2 * 1.05)
            ax.set_ylim(-L0 / 2 * 1.05, L0 / 2 * 1.05)
            ax.set_aspect("equal")
            ax.set_xlabel("x / d"); ax.set_ylabel("y / d")
            phi = len(pos) * math.pi / 6.0 / L ** 3
            ax.set_title(f"compression gelation   L = {L:.2f} d   "
                         rf"$\phi$ = {phi:.4f}", fontsize=10)
            writer.grab_frame()
    plt.close(fig)
    print(f"  → {out.relative_to(ROOT)}  ({len(frames)} frames)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=None, help="run_id (default: every completed network run)")
    ap.add_argument("--no-anim", action="store_true")
    ap.add_argument("--stride", type=int, default=2)
    args = ap.parse_args()

    dirs = ([ROOT / "runs" / args.run] if args.run else
            sorted(p for p in (ROOT / "runs").glob("network__N*")
                   if (p / "metrics.json").exists()))
    done = []
    for p in dirs:
        try:
            d = load(p)
            if "final_positions" not in d:
                print(f"  skipped (finalize incomplete): {p.name}")
                continue
            done.append(d)
        except Exception as e:                                    # noqa: BLE001
            print(f"  skipped {p.name}: {e}")
    if not done:
        print("no completed network runs.")
        return 1

    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(done)} runs")
    for d in done:
        st = d["_metrics"]["numerics"].get("stage_tau", "?")
        figure_structure(d, OUT / f"structure_st{st}.png")
    if len(done) >= 2:
        figure_compare(done, OUT / "compare_rate.png")
    if not args.no_anim:
        import matplotlib.animation                                # noqa: F401
        for d in done:
            st = d["_metrics"]["numerics"].get("stage_tau", "?")
            animate(ROOT / "runs" / d["_id"], OUT / f"compress_st{st}.mp4",
                    stride=args.stride)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
