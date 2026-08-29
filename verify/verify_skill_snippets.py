"""스킬 문서의 코드가 실제로 동작하는지 검증한다.

깨진 스니펫이 든 스킬은 없느니만 못하다. 문서를 고칠 때마다 이걸 돌린다.

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

# ── 1. 모든 python 블록 문법 검사 ────────────────────────────────────────
blocks = re.findall(r"```python\n(.*?)```", SKILL.read_text(), re.S)
print(f"python 블록 {len(blocks)}개")
for i, b in enumerate(blocks, 1):
    try:
        ast.parse(b)
    except SyntaxError as e:
        failures.append(f"블록 {i} 문법 오류: {e}")
        print(f"  ✗ 블록 {i}: {e}")
print(f"  문법: {len(blocks) - len([f for f in failures if '문법' in f])}/{len(blocks)} OK")


# ── 2. 조화 트랩 스니펫 실동작 + 물리 정확도 ─────────────────────────────
# 문서에 적힌 클래스를 그대로 exec 해서 쓴다. 여기 코드를 복붙하지 않는 것이 요점 —
# 문서가 틀리면 이 테스트가 깨져야 한다.
print("\n조화 트랩 스니펫 (문서 코드 그대로 실행)")
try:
    trap_src = next(b for b in blocks if "class HarmonicTrap" in b)
except StopIteration:
    print("  ✗ HarmonicTrap 블록을 찾지 못함 — 문서가 바뀌었나?")
    sys.exit(1)

ns = {"md": md, "np": np, "hoomd": hoomd}
exec(trap_src, ns)
HarmonicTrap = ns["HarmonicTrap"]

# k=2 는 함정 1(최소 이미지)에 가장 취약한 조건. 일부러 약한 트랩으로 검증한다.
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

print(f"  <x²> = {mean:.5f}   예측 kT/k = {pred:.5f}   오차 {err:+.2f}%")
print(f"  force NaN = {n_nan} (함정 7)   z성분 최대 = {max_fz:.1e} (2D라 0)")

if abs(err) >= 5:
    failures.append(f"트랩 물리 오차 {err:+.2f}% (한계 5%)")
if n_nan:
    failures.append(f"force 배열에 NaN {n_nan}개 — 함정 7 재발")
if max_fz > 1e-12:
    failures.append(f"2D인데 z 힘이 0이 아님: {max_fz:.1e}")

# ── 결과 ────────────────────────────────────────────────────────────────
print()
if failures:
    print("✗ FAIL")
    for f_ in failures:
        print(f"   - {f_}")
    sys.exit(1)
print("✓ PASS — 스킬 문서의 코드가 그대로 동작하고 물리적으로 정확함")
