"""Verify that the code in the skill documents actually runs.

A skill containing a broken snippet is worse than no skill at all. Run this every
time the documents are edited.

  PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
  $PY scratch/verify_skill_snippets.py
"""
import ast
import math
import re
import sys
from pathlib import Path

import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

ROOT = Path(__file__).resolve().parent.parent
SKILL = ROOT / ".claude/skills/bd-hoomd/SKILL.md"

failures = []

# ── 1. syntax-check every python block ──────────────────────────────────
blocks = re.findall(r"```python\n(.*?)```", SKILL.read_text(), re.S)
print(f"{len(blocks)} python blocks")
for i, b in enumerate(blocks, 1):
    try:
        ast.parse(b)
    except SyntaxError as e:
        failures.append(f"block {i} syntax error: {e}")
        print(f"  ✗ block {i}: {e}")
print(f"  syntax: "
      f"{len(blocks) - len([f for f in failures if 'syntax error' in f])}"
      f"/{len(blocks)} OK")


# ── 2. does the harmonic-trap snippet run, and is the physics right? ────
# The class written in the document is exec'd as-is. Not copy-pasting the code here
# is the whole point -- if the document is wrong, this test must break.
print("\nharmonic-trap snippet (running the document's code verbatim)")
try:
    trap_src = next(b for b in blocks if "class HarmonicTrap" in b)
except StopIteration:
    print("  ✗ could not find the HarmonicTrap block -- has the document changed?")
    sys.exit(1)

ns = {"md": md, "np": np, "hoomd": hoomd}
exec(trap_src, ns)
HarmonicTrap = ns["HarmonicTrap"]

# k=2 is the condition most vulnerable to trap 1 (minimum image). A deliberately
# weak trap is used for the check.
N, L, kT, gamma, k = 256, 40.0, 1.0, 1.0, 2.0
n = int(math.ceil(math.sqrt(N)))
a = L / n
f = gsd.hoomd.Frame()
f.particles.N = N
f.particles.position = [[(i % n + .5) * a - L / 2, (i // n + .5) * a - L / 2, 0.]
                        for i in range(N)]
f.particles.typeid = [0] * N
f.particles.types = ["A"]
f.configuration.box = [L, L, 0, 0, 0, 0]
f.configuration.dimensions = 2

sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=21)
sim.create_state_from_snapshot(f)
anchors = np.array(f.particles.position)
trap = HarmonicTrap(k, anchors, L, dimensions=2)
bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=kT, default_gamma=gamma)
tau = gamma / k
dt = tau / 2000
integ = md.Integrator(dt=dt, methods=[bd], forces=[trap])
integ.integrate_rotational_dof = False
sim.operations.integrator = integ

sim.run(int(20 * tau / dt))
samples = []
for _ in range(200):
    sim.run(int(0.5 * tau / dt))
    samples.append((trap.displacements(sim.state)[:, :2] ** 2).mean(axis=0))

mean = float(np.array(samples).mean())
pred = kT / k
err = 100 * (mean - pred) / pred
n_nan = int(np.isnan(np.array(trap.forces)).sum())
max_fz = float(np.abs(np.array(trap.forces)[:, 2]).max())

print(f"  <x^2> = {mean:.5f}   predicted kT/k = {pred:.5f}   error {err:+.2f}%")
print(f"  force NaN = {n_nan} (trap 7)   max z component = {max_fz:.1e} "
      f"(0, being 2D)")

if abs(err) >= 5:
    failures.append(f"trap physics error {err:+.2f}% (limit 5%)")
if n_nan:
    failures.append(f"{n_nan} NaN in the force array -- trap 7 has recurred")
if max_fz > 1e-12:
    failures.append(f"2D but the z force is not 0: {max_fz:.1e}")

# ── result ──────────────────────────────────────────────────────────────
print()
if failures:
    print("✗ FAIL")
    for f_ in failures:
        print(f"   - {f_}")
    sys.exit(1)
print("✓ PASS -- the skill documents' code runs verbatim and is physically correct")
