"""탐색 — chain-relax-2d-dlvo 사슬의 중앙 비드를 y(t)=A sin(ωt) 로 직접 강제.

chain-relax-2d-dlvo(마찰 없음, DLVO 인력만, 명시적 굽힘 없음, **양끝도 트랩 없이 자유**)
에 국소 구동 하나만 얹어 "이 구동이 사슬을 따라 얼마나 전파되는가"를 본다.
chain-bend-2d-dlvo 의 `--drive-mode position` 과 같은 강제 방식(트랩이 아니라 비드
위치를 직접 씀, `_move_ghost_action` 재사용)이지만 **양끝을 트랩으로 잡지 않는다** —
그래서 트랩 컴플라이언스가 변위를 흡수해주지 않는다. G1(중심력+자연장 ⟹ 횡방향 결합에너지가
O(y⁴))에 의해, 작은 진폭에서는 이웃 결합의 신장이 O(y²)로만 늘어난다 — 즉 **진폭이
작으면 구동이 거의 전파되지 않을 것**이라는 예측이 나온다(사전 예측, `--amp` 로 확인).

★★ 진폭 안전 한계 — 트랩 컴플라이언스가 없어 100% 가 결합 신장으로 간다:
    induced stretch ≈ amp² / (2·ell)   (기하, O(y²))
    F_max 여유(h_min→h_min+3.46nm) 안에 있어야 파단하지 않는다
    amp=50nm → stretch 0.84nm (안전) / amp=632nm(chain-bend-2d-dlvo 관례값) → 134.8nm (파단)
→ 기본 진폭을 **50nm**로 낮춰 잡았다 (사용자가 큰 진폭을 원하면 --amp 로 올리되 파단
가능성을 감수해야 한다 — 그 자체도 결과다).

★ 시각화/탐색 전용 스크립트다 — L3 스펙·health 게이트를 거치지 않는다(생산 런 아님).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/drive_chain_relax_center.py --omega 3000 --amp 50 --cycles 20
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

from chain_bend_dlvo_2d import (  # noqa: E402
    SIGMA_CORE_STAR, build_table_arrays, _move_ghost_action, find_well, dlvo_reduced_params,
    F_h_star,
)
from chain_relax_2d_dlvo import kink_positions, bend_angles, bow_metrics  # noqa: E402
from bdbot import materials as M, sim as SIM  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
from bdbot.units import Q  # noqa: E402
import yaml  # noqa: E402


def load_physics():
    raw = yaml.safe_load((ROOT / "intake/chain-bend-2d-dlvo/system.yaml").read_text())
    P = load_node
    sys_ = {"d": P(raw["particle"]["diameter"]), "T": P(raw["medium"]["temperature"]),
           "eta": P(raw["medium"]["viscosity"]), "rho_p": P(raw["particle"]["density"]),
           "eps_r": P(raw["medium"]["relative_permittivity"]), "psi0": P(raw["particle"]["surface_potential"]),
           "ionic_strength": P(raw["medium"]["ionic_strength"]),
           "A_H": P(raw["interactions"][0]["hamaker_constant"])}
    return sys_


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=9)
    ap.add_argument("--omega", type=float, default=3000.0, help="rad/s")
    ap.add_argument("--amp", type=float, default=50.0, help="nm")
    ap.add_argument("--cycles", type=float, default=20.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--eq-tau-bond", type=float, default=200.0)
    ap.add_argument("--kT", type=float, default=1.0, help="0.0 = 결정론(OverdampedViscous), "
                    "1.0 = 실제 열적(Brownian) — kT=0 vs kT=1 나란히 비교하려면 두 번 실행")
    ap.add_argument("--vdw-amp", type=float, default=None,
                    help="vdW 진폭 재조정(무차원, = A_H/(12kT)). 기본은 문헌값(0.2113). "
                         "0.2000 이면 U(h*=0.1)=-1.00kT (장벽·2차극소 유지 확인됨, "
                         "barrier=427.6kT well=-10.98kT). edl_amp/kappa_star(이온강도)는 그대로 둔다")
    ap.add_argument("--r-cut", type=float, default=None,
                    help="표 컷오프 r*=r/d 재조정(기본 1.06). h=1d 까지 닿으려면 2.0")
    ap.add_argument("--r-min", type=float, default=None,
                    help="표 하한 r*(기본 1.000001, 함정11: r<r_min 이면 F=0). ★ 장벽이 없는 "
                         "vdw_amp(예: 2.0)와 같이 쓸 때 필수 — 그 근처 힘이 h→0 로 발산해서 "
                         "r_min 을 물리 접촉(1.0)에서 살짝 띄워야 dt 가 유한해진다(예: 1.01)")
    args = ap.parse_args()

    sys_ = load_physics()
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    gamma, tau_B = b["gamma"], b["tau_B"]
    p = dlvo_reduced_params(sys_)
    if args.vdw_amp is not None:
        p = dict(p, vdw_amp=args.vdw_amp)
    w = find_well(p)
    print(f"DLVO: vdw_amp={p['vdw_amp']:.4f}  barrier_h*={w['barrier_h']:.5f} "
         f"({w['barrier_U']:.1f}kT)  h_min*={w['h_min']:.5f}  well={w['U_min']:.2f}kT  "
         f"k_bond*={w['k_bond_star']:.4e}")

    no_barrier_mode = args.r_min is not None
    if no_barrier_mode:
        # ★★ 장벽이 없으면(2차극소도 없음) find_well() 의 k_bond_star 는 무의미하다 —
        #   테이블 하한(r_min)에서의 최대 힘으로 직접 dt 를 역산한다(안전여유 1%,
        #   bd-hoomd Guard 의 step_disp_max 관례와 같은 정신). r_min 을 물리 접촉(r*=1)
        #   보다 살짝 바깥에 둬서 진짜 h→0 발산을 안 만나게 한다(그 근처는 WCA+표 하한의
        #   "죽은 구간"이 대신 막아준다 — 표는 r<r_min 에서 F=0, WCA 는 r<1 에서만 반발).
        r_min_star = args.r_min
        F_at_rmin = abs(float(F_h_star(r_min_star - 1.0, p)))
        dt = 0.01 / F_at_rmin
        ell_star = r_min_star
        print(f"장벽 없음 모드 — r_min*={r_min_star:.4f} 에서 |F|={F_at_rmin:.2f} kT/d, "
             f"dt(1%% 여유)={dt:.4e} τ_B")
    else:
        h_min_star = w["h_min"]
        ell_star = 1.0 + h_min_star
        k_bond = Q(w["k_bond_star"], "dimensionless") * b["kT"] / d ** 2
        tau_bond = (gamma / k_bond).to("s")
        dt_star = 1e-2                                  # dt/tau_bond, 이 계의 게이트와 동일
        dt = dt_star * float((tau_bond / tau_B).to("dimensionless").magnitude)  # dt*=dt/tau_B

    omega_star = float((Q(args.omega, "1/s") * tau_B).to("dimensionless").magnitude)
    amp_star = float((Q(args.amp, "nm") / d).to("dimensionless").magnitude)

    tau_period_star = 2 * math.pi / omega_star     # 주기(τ_B 단위)
    period_steps = int(round(tau_period_star / dt))
    n_cycles = args.cycles
    n_prod = int(round(n_cycles * tau_period_star / dt))
    frames_per_cycle = 24
    n_frames = int(round(n_cycles * frames_per_cycle))
    capture_every = max(1, n_prod // n_frames)
    # ZOH(bd-hoomd 함정17) 비율 UPDATE_EVERY/period_steps 을 1% 이하로 — dt 가 계마다
    # 크게 달라져서(장벽 없음 모드는 100배 큼) 고정값 50은 period_steps 가 작을 때
    # 사인파를 계단으로 뭉갠다(실측: period=332 에서 50이면 한 주기에 6~7번만 갱신).
    UPDATE_EVERY = max(1, min(50, period_steps // 100))

    n, mid = args.n, args.n // 2
    print(f"n={n} mid={mid}  omega={args.omega:.0f} rad/s (omega*={omega_star:.4e})  "
         f"amp={args.amp:.0f} nm (amp*={amp_star:.4e})")
    print(f"dt*={dt:.4e}  period={period_steps:,} steps  n_prod={n_prod:,} steps "
         f"({n_cycles:g} cycles)  capture every {capture_every} steps -> {n_prod//capture_every} frames")

    import hoomd
    import hoomd.md as md

    r_cut_star = args.r_cut if args.r_cut is not None else 1.0 + 0.06
    r_min_used = args.r_min if args.r_min is not None else 1.0 + 1e-6
    box_star = max(4.0 * (n - 1) * ell_star, 4.0 * r_cut_star)
    pos0 = kink_positions(n, ell_star, 0.0)             # 직선, 자연장(또는 r_min)에서 시작
    sim = SIM.make_sim(SIM.frame_2d(pos0, box_star), seed=args.seed)

    cell = md.nlist.Cell(buffer=0.2)
    reduced = {k: p[k] for k in ("kappa_star", "edl_amp", "vdw_amp", "a_star")}
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min_used, r_cut_star)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut_star)
    tab.params[("A", "A")] = dict(r_min=r_min_used, U=U_arr, F=F_arr)
    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * 2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)

    # ── ① 짧은 평형화 — 전 입자 자유(구동 아직 없음). kT=0 이면 뜻이 없어 건너뛴다 ──
    n_eq = 0
    if args.kT > 0:
        bd_all = md.methods.Brownian(filter=hoomd.filter.All(), kT=args.kT, default_gamma=1.0)
        integ = md.Integrator(dt=dt, methods=[bd_all], forces=[tab, wca])
        integ.integrate_rotational_dof = False
        sim.operations.integrator = integ
        if no_barrier_mode:
            # tau_bond 개념이 없다(2차극소가 없어 k_bond_star 도 없음) — r_min 경계의
            # "죽은 구간"(표 F=0 · WCA F=0) 안에서 자리잡는 데 필요한 스텝수를 그냥
            # dt 배수로 잡는다(고정, ★탐색적 어림 — 엄밀한 시간척도 아님)
            # ★ 짧게 — 초기 배치(정확히 r_min 경계)의 국소 지터만 재우는 목적이다.
            #   길게 두면(예: 2e4) 굽힘강성이 없는 사슬이 그 사이에 이미 랜덤 코일로
            #   재배열해버려서(실측: bead 위치가 ±400nm 로 흩어짐 — 개별 결합의 죽은
            #   구간 0.01d 폭보다 훨씬 큼, 다관절 누적 회전 탓) 구동 실험의 "직선에서
            #   시작" 전제가 깨진다.
            n_eq = int(args.eq_tau_bond * 2)
        else:
            n_eq = int(round(args.eq_tau_bond
                            * (tau_bond / (dt * tau_B)).to("dimensionless").magnitude))
        sim.run(n_eq)
    print(f"평형화 {n_eq:,} 스텝 완료 (⟨U⟩/N 안정화용, 구동 없음)"
         + (" — kT=0 이라 건너뜀" if args.kT == 0 else ""))

    # ── ② 중앙 비드를 적분에서 빼고 위치를 직접 강제 ─────────────────────
    bd_filter = hoomd.filter.SetDifference(hoomd.filter.All(), hoomd.filter.Tags([mid]))
    if args.kT > 0:
        bd = md.methods.Brownian(filter=bd_filter, kT=args.kT, default_gamma=1.0)
    else:
        bd = md.methods.OverdampedViscous(filter=bd_filter, default_gamma=1.0)
    integ2 = md.Integrator(dt=dt, methods=[bd], forces=[tab, wca])
    integ2.integrate_rotational_dof = False
    sim.operations.integrator = integ2
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=_move_ghost_action(mid, amp_star, omega_star, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))

    # ── ③ 구동 구간 — 프레임 캡처 ─────────────────────────────────────
    frames = [np.array(sim.state.get_snapshot().particles.position, dtype=float)[:, :2].copy()]
    done = 0
    while done < n_prod:
        chunk = min(capture_every, n_prod - done)
        sim.run(chunk)
        done += chunk
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, dtype=float)[:, :2]
        img = np.array(snap.particles.image, dtype=float)[:, :2]
        frames.append(pos + img * box_star)
    frames = np.array(frames)
    print(f"구동 완료 — 프레임 {len(frames)}개 캡처")

    t_star = np.arange(len(frames)) * capture_every * dt          # τ_B 단위
    t_over_period = t_star / tau_period_star
    y_all = frames[:, :, 1] - frames[:1, :, 1].mean()              # 대략 중심 보정(정보용)

    out_dir = ROOT / "runs" / "_scratch_drive_chain_relax"
    out_dir.mkdir(exist_ok=True)
    tag = f"n{n}-w{args.omega:g}-a{args.amp:g}-kT{args.kT:g}"
    if no_barrier_mode:
        tag += f"-vdw{p['vdw_amp']:g}-rc{r_cut_star:g}-rm{r_min_used:g}"
    np.savez_compressed(out_dir / f"{tag}.npz", frames=frames, t_over_period=t_over_period,
                       omega=args.omega, amp_nm=args.amp, mid=mid, ell_star=ell_star,
                       box_star=box_star)
    make_plots(frames, t_over_period, mid, n, args, out_dir, tag)


def make_plots(frames, t_over_period, mid, n, args, out_dir, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.animation as anim
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    d_nm = 1470.0
    y_nm = frames[:, :, 1] * d_nm

    # ① y_i(t) — 비드마다, 중앙에서 가장자리로 갈수록 색을 옅게
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    cmap = plt.cm.viridis
    for i in range(n):
        dist = abs(i - mid)
        ax1.plot(t_over_period, y_nm[:, i], lw=1.3 if dist == 0 else 0.9,
                 color=cmap(1 - dist / max(1, mid)), alpha=1.0 if dist == 0 else 0.85,
                 label=f"bead {i}" + ("  (driven)" if i == mid else f"  (|Δi|={dist})"))
    ax1.set(xlabel="t / period", ylabel="y [nm]",
           title=f"chain-relax-2d-dlvo, driven center bead — n={n}, "
                 f"ω={args.omega:g} rad/s, A={args.amp:g} nm, kT={args.kT:g}")
    ax1.legend(fontsize=7, ncol=2)
    ax1.grid(alpha=.3)
    fig1.tight_layout()
    fig1.savefig(out_dir / f"{tag}_y_of_t.png", dpi=140)
    plt.close(fig1)

    # ② 진폭 전파 — 각 비드의 y 진폭(정상상태 구간 RMS*sqrt(2)) vs |i-mid|
    steady = t_over_period > t_over_period.max() * 0.4
    amp_meas = y_nm[steady].std(axis=0) * math.sqrt(2)
    fig2, ax2 = plt.subplots(figsize=(6, 4.5))
    dist = np.abs(np.arange(n) - mid)
    ax2.semilogy(dist, amp_meas, "o-")
    ax2.axhline(args.amp, color="gray", ls=":", lw=1, label=f"driven amplitude ({args.amp:g} nm)")
    ax2.set(xlabel="|bead index - driven bead|", ylabel="oscillation amplitude [nm] (log)",
           title="Amplitude vs distance from drive\n(no bending stiffness -> expect fast decay)")
    ax2.legend(fontsize=8); ax2.grid(alpha=.3, which="both")
    fig2.tight_layout()
    fig2.savefig(out_dir / f"{tag}_propagation.png", dpi=140)
    plt.close(fig2)
    print("전파 진폭 [nm] vs |Δi|:", dict(zip(dist.tolist(), np.round(amp_meas, 3).tolist())))

    # ③ 애니메이션 — 몸좌표계(양끝 대신, 구동 비드를 원점으로) 라벨과 함께
    xmax = float(np.abs(frames[:, :, 0] - frames[:, mid:mid+1, 0]).max()) * 1.15
    ymax = max(args.amp / d_nm * 1.4, float(np.abs(y_nm).max()) / d_nm * 1.15)
    fig3, ax3 = plt.subplots(figsize=(7, 3.6))
    ax3.set_xlim(-xmax, xmax); ax3.set_ylim(-ymax, ymax)
    ax3.set_aspect("equal"); ax3.grid(alpha=.3)
    ax3.set_xlabel("x - x_driven [d]"); ax3.set_ylabel("y [d]")
    ax3.set_title(f"n={n}, ω={args.omega:g} rad/s, A={args.amp:g} nm, kT={args.kT:g}")
    line, = ax3.plot([], [], "o-", color="tab:purple", ms=6, lw=1.5)
    driven_pt, = ax3.plot([], [], "o", color="tab:red", ms=9)
    txt = ax3.text(0.02, 0.92, "", transform=ax3.transAxes, fontsize=9)

    def update(i):
        xr = frames[i, :, 0] - frames[i, mid, 0]
        yr = frames[i, :, 1]
        line.set_data(xr, yr)
        driven_pt.set_data([xr[mid]], [yr[mid]])
        txt.set_text(f"t = {t_over_period[i]:.2f} periods")
        return line, driven_pt, txt

    ani = anim.FuncAnimation(fig3, update, frames=len(frames), interval=60, blit=False)
    ani.save(out_dir / f"{tag}_anim.gif", writer="pillow", fps=15)
    plt.close(fig3)
    print("저장:", out_dir / f"{tag}_y_of_t.png", "/", out_dir / f"{tag}_propagation.png",
         "/", out_dir / f"{tag}_anim.gif")


if __name__ == "__main__":
    main()
