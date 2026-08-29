"""(1) gamma_r 축별 이방성은 실제로 동작하는가?
   (2) 자유배수(free-draining) 비드 사슬은 이방성 병진 마찰을 주는가?
"""
import numpy as np, gsd.hoomd, hoomd, hoomd.md as md

# ── (1) 회전 항력 gamma_r 축별 동작 확인 ────────────────────────────────
print("=" * 84)
print("(1) gamma_r 축별 이방성 — 일정 토크에 대한 각속도 (Langevin, kT=0 → 결정론)")
print("=" * 84)

def spin(gamma_r, torque_axis, steps=2000, dt=1e-3):
    f = gsd.hoomd.Frame()
    f.particles.N = 1
    f.particles.types = ["R"]
    f.particles.typeid = [0]
    f.particles.position = [[0, 0, 0]]
    f.particles.orientation = [[1, 0, 0, 0]]
    f.particles.moment_inertia = [[1.0, 1.0, 1.0]]
    f.particles.mass = [1.0]
    f.configuration.box = [50, 50, 50, 0, 0, 0]
    f.configuration.dimensions = 3
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=3)
    sim.create_state_from_snapshot(f)

    const = md.force.Constant(filter=hoomd.filter.All())
    const.constant_force["R"] = (0, 0, 0)
    const.constant_torque["R"] = tuple(float(x) for x in torque_axis)

    # kT=0 인 Brownian → 잡음 없음, dq/dt = τ/γ_r
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0,
                             default_gamma=1.0, default_gamma_r=tuple(gamma_r))
    integ = md.Integrator(dt=dt, methods=[bd], forces=[const],
                          integrate_rotational_dof=True)
    sim.operations.integrator = integ
    sim.run(steps)
    q = np.array(sim.state.get_snapshot().particles.orientation[0], dtype=float)
    # 회전각 = 2*arccos(w)
    return 2 * np.arccos(np.clip(abs(q[0]), -1, 1))

for gr, label in [((1.0, 1.0, 1.0), "등방 (1,1,1)"),
                  ((1.0, 1.0, 5.0), "z축 5배 (1,1,5)")]:
    ax = spin(gr, (1.0, 0, 0))
    az = spin(gr, (0, 0, 1.0))
    print(f"  gamma_r={gr:} [{label}]")
    print(f"     τ∥x → 회전각 {ax:.5f} rad     τ∥z → 회전각 {az:.5f} rad     비 z/x = {az/ax if ax else 0:.4f}")
print("  → z성분을 5배로 하면 z축 회전이 1/5로 느려져야 함 (비 ≈ 0.2)")

# ── (2) 자유배수 비드 사슬: 이방성이 나오는가? ──────────────────────────
print()
print("=" * 84)
print("(2) 자유배수 비드 막대 (강체 아님, 각 비드가 독립 항력) — 이방성 나오는가?")
print("=" * 84)

def bead_rod_gamma(force_dir, n_beads=5, k_bond=2000.0, steps=4000, dt=1e-4):
    f = gsd.hoomd.Frame()
    f.particles.N = n_beads
    f.particles.types = ["A"]
    f.particles.typeid = [0] * n_beads
    offs = (np.arange(n_beads) - (n_beads - 1) / 2) * 1.0
    f.particles.position = [[float(o), 0.0, 0.0] for o in offs]
    f.particles.mass = [1.0] * n_beads
    f.configuration.box = [200, 200, 200, 0, 0, 0]
    f.configuration.dimensions = 3
    f.bonds.N = n_beads - 1
    f.bonds.types = ["b"]
    f.bonds.typeid = [0] * (n_beads - 1)
    f.bonds.group = [[i, i + 1] for i in range(n_beads - 1)]
    f.angles.N = n_beads - 2
    f.angles.types = ["a"]
    f.angles.typeid = [0] * (n_beads - 2)
    f.angles.group = [[i, i + 1, i + 2] for i in range(n_beads - 2)]

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=5)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic(); bond.params["b"] = dict(k=k_bond, r0=1.0)
    ang = md.angle.Harmonic(); ang.params["a"] = dict(k=k_bond, t0=np.pi)
    const = md.force.Constant(filter=hoomd.filter.All())
    # 막대 전체에 총 힘 F가 걸리도록 비드당 F/n
    const.constant_force["A"] = tuple(float(x) / n_beads for x in force_dir)
    const.constant_torque["A"] = (0, 0, 0)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.All(), default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=dt, methods=[ov],
                                              forces=[bond, ang, const])
    com0 = np.array(sim.state.get_snapshot().particles.position).mean(axis=0)
    sim.run(steps)
    com1 = np.array(sim.state.get_snapshot().particles.position).mean(axis=0)
    v = float(np.dot(com1 - com0, np.array(force_dir, float)) / (steps * dt))
    return 1.0 / v if abs(v) > 1e-14 else float("inf")   # F_total=1 이므로 γ=1/v

gp = bead_rod_gamma((1.0, 0, 0))
gt = bead_rod_gamma((0, 1.0, 0))
print(f"  γ∥ = {gp:.5f}   γ⊥ = {gt:.5f}   γ⊥/γ∥ = {gt/gp:.6f}")
print(f"  (비드 5개 × γ_bead=1 → 총 γ = 5 예상)")
print()
print("  → 자유배수에서는 각 비드가 독립적으로 Stokes 항력을 받으므로")
print("    총 항력 = N·γ_bead 로 방향에 무관. 이방성은 '유체역학적 상호작용(HI)'의 효과이고,")
print("    HI가 없는 BD에서는 기하만으로 이방성이 생기지 않는다.")
