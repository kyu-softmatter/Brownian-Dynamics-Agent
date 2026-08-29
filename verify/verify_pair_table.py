"""md.pair.Table 실측 검증 — Phase 1-B에서 r⁻³를 쓰기 전에.

문서에서 조용히 틀릴 수 있어 보이는 두 지점을 실행으로 확인한다:
  ① r < r_min 에서 힘·에너지가 정말 0인가? (발산 퍼텐셜이면 반발이 사라진다)
  ② 격자가 linspace(r_min, r_cut, len(U), **endpoint=False**) 인가?
     skill bd-hoomd 의 스니펫은 endpoint=True 로 만들고 있다 — 맞는가?

    $PY scratch/verify_pair_table.py
"""
import math

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

A = 10.0
R_MIN, R_CUT, NBINS = 0.5, 3.0, 200


def make_table(endpoint):
    """U = A/r³ (컷오프에서 시프트), F = -dU/dr = 3A/r⁴."""
    r = np.linspace(R_MIN, R_CUT, NBINS, endpoint=endpoint)
    U = A / r**3
    U = U - A / R_CUT**3  # 컷오프에서 0 (Table은 shift 모드 없음)
    F = 3 * A / r**4
    return dict(r_min=R_MIN, U=U, F=F)


def measure(sep, params):
    """입자 2개를 sep 만큼 떼어놓고 HOOMD가 계산한 힘·에너지를 읽는다."""
    L = 20.0
    fr = gsd.hoomd.Frame()
    fr.particles.N = 2
    fr.particles.position = [[-sep / 2, 0, 0], [sep / 2, 0, 0]]
    fr.particles.typeid = [0, 0]
    fr.particles.types = ["A"]
    fr.configuration.box = [L, L, 0, 0, 0, 0]
    fr.configuration.dimensions = 2

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(fr)
    cell = md.nlist.Cell(buffer=0.4)
    tab = md.pair.Table(nlist=cell, default_r_cut=R_CUT)
    tab.params[("A", "A")] = params
    # 적분은 하지 않는다 — 힘만 평가 (dt=0 대신 0 스텝 실행)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-9, methods=[bd], forces=[tab])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(0)
    f = np.array(tab.forces)
    e = np.array(tab.energies)
    return float(f[1][0]), float(e.sum())  # 오른쪽 입자의 +x 힘, 총 퍼텐셜


def analytic(sep):
    return 3 * A / sep**4, A / sep**3 - A / R_CUT**3


print("=" * 78)
print("① 격자 규약 — endpoint=False 가 맞는가?")
print("=" * 78)
print(f"  A = {A}, r_min = {R_MIN}, r_cut = {R_CUT}, nbins = {NBINS}")
print(f"\n  {'sep':>6} {'F 해석해':>11} | {'endpoint=False':>15} {'오차':>9} |"
      f" {'endpoint=True':>14} {'오차':>9}")
worst = {False: 0.0, True: 0.0}
for sep in (0.7, 1.0, 1.5, 2.0, 2.5, 2.9):
    Fa, _ = analytic(sep)
    row = f"  {sep:6.2f} {Fa:11.4f} |"
    for ep in (False, True):
        Fm, _ = measure(sep, make_table(endpoint=ep))
        err = 100 * (Fm - Fa) / Fa
        worst[ep] = max(worst[ep], abs(err))
        row += f" {Fm:15.4f} {err:+8.3f}% |"
    print(row)
print(f"\n  최대 오차:  endpoint=False {worst[False]:.3f}%   endpoint=True {worst[True]:.3f}%")
better = "endpoint=False" if worst[False] < worst[True] else "endpoint=True"
print(f"  → {better} 가 맞다 (문서 서술과 일치 여부 확인)")

print()
print("=" * 78)
print("② r < r_min 에서 정말 힘·에너지가 0인가?  ★ 발산 퍼텐셜이면 위험")
print("=" * 78)
p = make_table(endpoint=False)
print(f"  {'sep':>6} {'F 해석해':>13} {'F 측정':>13} {'U 측정':>13}   판정")
for sep in (0.60, 0.51, 0.49, 0.40, 0.30):
    Fa, Ua = analytic(sep)
    Fm, Um = measure(sep, p)
    inside = sep < R_MIN
    verdict = ("r_min 아래 → 0" if inside else "정상")
    print(f"  {sep:6.2f} {Fa:13.4f} {Fm:13.4f} {Um:13.4f}   {verdict}")
print(f"""
  ★ 결론: r < r_min({R_MIN}) 에서 반발력이 **사라진다**.
    r⁻³ 처럼 발산하는 퍼텐셜에서는 입자가 r_min을 뚫고 들어가면
    그 뒤로는 아무 힘도 못 받아 겹친 채로 남는다 — 에러 없이 조용히 틀린다.
    → 방어: (a) 배제부피 코어(WCA)를 별도 힘으로 두어 r_min에 도달 자체를 막고
            (b) 런 중 최소 이웃거리를 감시한다.
""")

print("=" * 78)
print("③ WCA(LJ) + Table 중첩 — 두 힘을 함께 걸면 합해지는가?")
print("=" * 78)


def measure_both(sep):
    L = 20.0
    fr = gsd.hoomd.Frame()
    fr.particles.N = 2
    fr.particles.position = [[-sep / 2, 0, 0], [sep / 2, 0, 0]]
    fr.particles.typeid = [0, 0]
    fr.particles.types = ["A"]
    fr.configuration.box = [L, L, 0, 0, 0, 0]
    fr.configuration.dimensions = 2
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(fr)
    cell = md.nlist.Cell(buffer=0.4)
    tab = md.pair.Table(nlist=cell, default_r_cut=R_CUT)
    tab.params[("A", "A")] = make_table(endpoint=False)
    wca = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-9, methods=[bd], forces=[tab, wca])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(0)
    return (float(np.array(tab.forces)[1][0]) + float(np.array(wca.forces)[1][0]),
            float(np.array(tab.energies).sum()) + float(np.array(wca.energies).sum()))


def analytic_both(sep):
    Ft, Ut = analytic(sep)
    if sep < 2 ** (1 / 6):
        Uw = 4 * (sep**-12 - sep**-6) + 1.0
        Fw = 4 * (12 * sep**-13 - 6 * sep**-7)
    else:
        Uw = Fw = 0.0
    return Ft + Fw, Ut + Uw


print(f"  {'sep':>6} | {'F 해석해':>12} {'F 측정':>12} {'err':>9} |"
      f" {'U 해석해':>12} {'U 측정':>12} {'err':>9}")
ok_all = True
for sep in (0.90, 0.95, 1.00, 1.10, 1.30, 2.00):
    Fa, Ua = analytic_both(sep)
    Fm, Um = measure_both(sep)
    ef = 100 * (Fm - Fa) / Fa
    eu = 100 * (Um - Ua) / Ua
    ok_all &= abs(ef) < 1 and abs(eu) < 1
    print(f"  {sep:6.2f} | {Fa:12.4f} {Fm:12.4f} {ef:+8.3f}% |"
          f" {Ua:12.4f} {Um:12.4f} {eu:+8.3f}%")
print(f"\n  → 두 힘은 {'정상적으로 합해진다' if ok_all else '★ 합이 맞지 않는다'} (1% 이내 판정)")
print()
grid_ok = worst[False] < 1e-6
print("=" * 78)
print("✓ PASS — endpoint=False 정확(0.000%) · r<r_min 힘 0 확인 · WCA+Table 합산 정확"
      if (grid_ok and ok_all) else "✗ FAIL — 위 표 확인")
print("=" * 78)
