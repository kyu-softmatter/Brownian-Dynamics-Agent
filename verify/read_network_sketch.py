#!/usr/bin/env python
"""Reproduce the transcription of the intake/network sketch (for L0 audit).

What it does:
  (1) correct EXIF orientation + rotate 90 deg CCW (the original was shot sideways)
  (2) **hard-code** the positions and radii of the 21 beads counted in the image and
      overlay them on the original
  (3) reconstruct the links by a contact test and print the topology (loop count,
      branch points, free ends)

★ The bead coordinates are NOT automatic detection -- they were **read by eye by a
  person** (the session that wrote this script). Automatic Hough detection was
  unavailable because this environment has neither opencv nor skimage. So this
  script's job is not "I counted them" but **"make the count auditable by someone
  else"** -- putting out/network_sketch_annotated.png next to the original makes any
  omission or duplicate visible.

  Transcription rule (skill bd-intake §1): if the particle count is small (<~50),
  count it exactly; if large, do not count, propose a suitably large N instead. Here
  there are 21, so it **was counted exactly**. That does not mean 21 is the
  simulation's N -- see observation.yaml A6.

Run:
  $PY scratch/read_network_sketch.py
"""
from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "intake/network/KakaoTalk_Photo_2026-08-06-14-36-42.jpeg"
OUT = ROOT / "intake/network/transcription_annotated.png"

# Rotation/crop convention: exif_transpose -> ROTATE_90(CCW) ->
# crop(260,600,1920,2260)
# Beads (x, y, r) read in that coordinate system (1660x1660). r is eyeballed -- it is
# used only for the contact test.
CROP = (260, 600, 1920, 2260)
BEADS = [
    (700, 260, 150), (795, 510, 155),                                   # upper branch
    (1070, 665, 95), (1215, 570, 80), (1335, 520, 85),                  # upper-right branch
    (940, 730, 105), (1015, 975, 105), (775, 865, 75), (650, 975, 80),  # centre
    (505, 1075, 65), (425, 1035, 45), (370, 1015, 35),                  # left branch
    (305, 975, 55), (185, 950, 72),
    (395, 1205, 62), (320, 1330, 58),                                   # lower-left branch
    (1105, 1175, 100), (1225, 1150, 65), (1345, 1055, 62), (1465, 1000, 68),  # lower-right branch
    (1160, 1420, 100),                                                  # lower branch
]
# The links (pairs of touching circles) are **also read by eye**. 1-based index pairs,
# ordered as BEADS above.
# ⚠️ Automatic reconstruction (centre distance < r_i+r_j) **does not work** with these
# coordinates -- see below.
LINKS_READ = [
    (1, 2), (2, 6), (6, 3), (3, 4), (4, 5),          # upper + upper-right branch
    (6, 7), (6, 8),                                   # 6 = 4-way branch point
    (8, 9), (9, 10), (10, 11), (11, 12), (12, 13), (13, 14),   # left branch
    (10, 15), (15, 16),                               # 10 = 3-way branch point
    (7, 17), (17, 18), (18, 19), (19, 20), (17, 21),  # 17 = 3-way branch point
]


def contact_ratio(beads, i, j):
    """Centre distance / (r_i + r_j). Below 1, the circles overlap."""
    xi, yi, ri = beads[i]
    xj, yj, rj = beads[j]
    return ((xi - xj) ** 2 + (yi - yj) ** 2) ** 0.5 / (ri + rj)


def main() -> int:
    if not SRC.exists():
        print(f"original not found: {SRC}", file=sys.stderr)
        return 1

    img = ImageOps.exif_transpose(Image.open(SRC)).transpose(Image.ROTATE_90)
    dia = img.crop(CROP).convert("RGB")

    lk = [(a - 1, b - 1) for a, b in LINKS_READ]
    deg = [0] * len(BEADS)
    for i, j in lk:
        deg[i] += 1
        deg[j] += 1

    # Connected components (one BFS, no union-find -- fine for 21 nodes)
    adj = {i: [] for i in range(len(BEADS))}
    for i, j in lk:
        adj[i].append(j)
        adj[j].append(i)
    seen, comps = set(), 0
    for s in range(len(BEADS)):
        if s in seen:
            continue
        comps += 1
        stack = [s]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            stack.extend(v for v in adj[u] if v not in seen)

    n, e = len(BEADS), len(lk)
    loops = e - n + comps          # cycle rank (Betti-1)

    print(f"beads              {n}")
    print(f"links              {e}   (read from the image)")
    print(f"components         {comps}")
    print(f"independent loops  {loops}   <- 0 means a tree "
          f"(no closed path to carry stress)")
    print(f"free ends (deg=1)  {sum(1 for d in deg if d == 1)}")
    print(f"3-way (deg=3)      {sum(1 for d in deg if d == 3)}")
    print(f"4-way (deg=4)      {sum(1 for d in deg if d == 4)}")
    print(f"mean coordination z {2 * e / n:.3f}")

    # ── Why automatic reconstruction is unusable (instead of tuning TOL, it is
    #    **measured and reported**) ──────
    lset = {frozenset(p) for p in lk}
    lr = sorted(contact_ratio(BEADS, i, j) for i, j in lk)
    nr = sorted(contact_ratio(BEADS, i, j)
                for i in range(n) for j in range(i + 1, n)
                if frozenset((i, j)) not in lset)
    print("\ncentre distance/(r_i+r_j) -- quality diagnostic for the eyeballed radii")
    print(f"  linked pairs     ({len(lr):2d})  min {lr[0]:.2f} . "
          f"median {lr[len(lr)//2]:.2f} . max {lr[-1]:.2f}")
    print(f"  non-linked pairs ({len(nr):3d})  min {nr[0]:.2f}")
    print(f"  -> max linked ({lr[-1]:.2f}) "
          f"{'>' if lr[-1] > nr[0] else '<'} min non-linked ({nr[0]:.2f})"
          f": distance alone is "
          f"{'NOT separable' if lr[-1] > nr[0] else 'separable'}")
    print("  ★ BEADS' r was chosen **small, to fit the index labels on top**, so it")
    print("    is not the semi-axis of the drawn ellipse. Automating the contact test")
    print("    therefore misses links, and raising TOL to fit would be fitting the")
    print("    answer we wanted -- so links are read by eye (CLAUDE.md rule 6).")

    dr = ImageDraw.Draw(dia)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 52)
    except OSError:
        font = None
    for i, j in lk:                                    # links first (drawn under the
                                                       # beads)
        dr.line([BEADS[i][:2], BEADS[j][:2]], fill=(0, 120, 255), width=7)
    for idx, (x, y, r) in enumerate(BEADS):
        dr.ellipse([x - r, y - r, x + r, y + r], outline=(255, 0, 0), width=5)
        dr.text((x - 18, y - 26), str(idx + 1), fill=(0, 150, 0), font=font)
    dia.save(OUT)
    print(f"\nannotated figure -> {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
