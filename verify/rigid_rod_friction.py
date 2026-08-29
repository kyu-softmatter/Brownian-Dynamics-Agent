"""강체 막대의 병진 마찰이 실제로 이방성인가? — 결정적 실측.

주장: "md.constrain.Rigid로 구를 여러 개 묶어 막대를 만들면 마찰이 자연스럽게 이방성이 된다"
검증: 방향을 고정한 막대에 일정한 힘을 축방향(∥)과 수직방향(⊥)으로 각각 걸고
      종단속도를 잰다. OverdampedViscous(열잡음 없음)에서 v = F/γ 이므로
      γ = F/v 를 직접 얻는다.

      γ∥ ≠ γ⊥  → 이방성 성립 (주장 옳음)
      γ∥ = γ⊥  → 등방 (주장 틀림, 대안 필요)

이론값(prolate spheroid, Perrin): 종횡비 p에서 γ⊥/γ∥ ≈ 1.4~2.0 (p=2~10)
막대 근사(slender body): γ⊥/γ∥ → 2 (p → ∞)
"""
import math
import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

N_BEADS = 5          # 막대를 이루는 구의 개수
BEAD_SEP = 1.0       # 구 간격 (body 좌표)
GAMMA_CENTER = 1.0   # 중심 입자에 지정할 스칼라 gamma
FORCE = 1.0
L_BOX = 200.0
N_STEPS = 2000
DT = 1e-3


def build(force_dir):
    """막대 강체 1개. force_dir 방향으로 일정한 힘. 회전 고정."""
    f = gsd.hoomd.Frame()
    f.particles.N = 1
    f.particles.types = ["R", "A"]          # R=중심(막대), A=구성 구
    f.particles.typeid = [0]
    f.particles.position = [[0.0, 0.0, 0.0]]
    f.particles.orientation = [[1.0, 0.0, 0.0, 0.0]]   # 막대축 = body x축 = lab x축
    f.particles.moment_inertia = [[0.0, 100.0, 100.0]]
    f.particles.body = [0]                  # 중심 입자는 자기 tag
    f.particles.mass = [float(N_BEADS)]
    f.configuration.box = [L_BOX, L_BOX, L_BOX, 0, 0, 0]
    f.configuration.dimensions = 3

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(f)

    # 강체 정의: body x축을 따라 N_BEADS개 구를 배치
    offs = (np.arange(N_BEADS) - (N_BEADS - 1) / 2) * BEAD_SEP
    rigid = md.constrain.Rigid()
    rigid.body["R"] = {
        "constituent_types": ["A"] * N_BEADS,
        "positions": [[float(o), 0.0, 0.0] for o in offs],
        "orientations": [[1.0, 0.0, 0.0, 0.0]] * N_BEADS,
    }
    rigid.create_bodies(sim.state)

    const = md.force.Constant(filter=hoomd.filter.Rigid(("center",)))
    const.constant_force["R"] = tuple(float(x) * FORCE for x in force_dir)
    const.constant_torque["R"] = (0.0, 0.0, 0.0)
    const.constant_force["A"] = (0.0, 0.0, 0.0)
    const.constant_torque["A"] = (0.0, 0.0, 0.0)

    ov = md.methods.OverdampedViscous(filter=hoomd.filter.Rigid(("center",)),
                                      default_gamma=GAMMA_CENTER)
    integ = md.Integrator(dt=DT, methods=[ov], forces=[const],
                          rigid=rigid, integrate_rotational_dof=False)
    sim.operations.integrator = integ
    return sim


def terminal_gamma(force_dir):
    sim = build(force_dir)
    p0 = np.array(sim.state.get_snapshot().particles.position[0], dtype=float)
    sim.run(N_STEPS)
    p1 = np.array(sim.state.get_snapshot().particles.position[0], dtype=float)
    disp = p1 - p0
    t = N_STEPS * DT
    v = float(np.dot(disp, np.array(force_dir, dtype=float)) / t)
    return FORCE / v if abs(v) > 1e-14 else float("inf"), disp, v


print("=" * 84)
print("강체 막대 병진 마찰 이방성 검증")
print(f"막대: 구 {N_BEADS}개, 간격 {BEAD_SEP}, body x축 방향  |  중심 gamma={GAMMA_CENTER} (스칼라)")
print("=" * 84)

g_par, d_par, v_par = terminal_gamma((1, 0, 0))     # 축방향 ∥
g_perp, d_perp, v_perp = terminal_gamma((0, 1, 0))  # 수직방향 ⊥

print(f"  ∥ (축방향)   변위={np.array2string(d_par, precision=5)}  v={v_par:.6f}  γ∥ ={g_par:.6f}")
print(f"  ⊥ (수직방향) 변위={np.array2string(d_perp, precision=5)}  v={v_perp:.6f}  γ⊥ ={g_perp:.6f}")
print()
ratio = g_perp / g_par
print(f"  γ⊥/γ∥ = {ratio:.6f}")
print()
if abs(ratio - 1.0) < 0.02:
    print("  ✗ 등방입니다. 강체로 묶어도 병진 마찰은 중심 입자의 스칼라 gamma 하나로 결정됩니다.")
    print("    → '강체로 묶으면 이방성이 자연히 나온다'는 주장은 HOOMD에서 성립하지 않습니다.")
else:
    print(f"  ✓ 이방성입니다 (γ⊥/γ∥ = {ratio:.3f}). 이론 기대치는 1.4~2.0 (slender body → 2).")
print("=" * 84)

# 참고: 구성 입자에도 gamma를 주면 달라지는가?
print()
print("보조 확인: 구성 입자 타입 'A'에 gamma를 따로 줘도 병진에 반영되는가?")
sim = build((1, 0, 0))
ov = sim.operations.integrator.methods[0]
ov.gamma["A"] = 10.0          # 구성 입자에 큰 항력
p0 = np.array(sim.state.get_snapshot().particles.position[0], dtype=float)
sim.run(N_STEPS)
p1 = np.array(sim.state.get_snapshot().particles.position[0], dtype=float)
v2 = float((p1 - p0)[0] / (N_STEPS * DT))
print(f"  gamma['A']=10 설정 후  v∥={v2:.6f}  (앞의 v∥={v_par:.6f})")
print("  → 값이 같으면 구성 입자의 gamma는 병진에 전혀 관여하지 않는다는 뜻.")
