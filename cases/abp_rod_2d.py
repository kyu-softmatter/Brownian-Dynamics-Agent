"""Phase 1-D — abp-rod-2d-run-flip 관통 (run-and-tumble 타원체).

물리계(SI) → 스케일 표 → 분리 검사 → 무차원화 → 실행 → 역변환 → 해석해 대조.

공통 부분은 `bdbot/` (Phase 1-C). 여기 남은 것은 이 계 고유의 것뿐입니다:
타원체 마찰(Perrin) · run-and-tumble 커스텀 updater · MSD·MSAD 관측량 · 해석해 2종.

★ 이 케이스에는 **해석해 예측이 두 개** 있습니다 (1-B에는 없었습니다):
    방향 상관   ⟨cos Δθ(t)⟩ = exp(−t/τ_eff),   1/τ_eff = D_r + 1/τ_tumble
    유효 확산   D_eff = D̄ + v²τ_eff/2          (2D)
  둘 다 텀블과 열적 회전확산이 **독립 과정**이라는 가정에서 나옵니다 — 그게 이 런의 검증 대상입니다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/abp_rod_2d.py              # 전체 (~14만 스텝 × N=1000)
    $PY cases/abp_rod_2d.py --smoke      # 짧게
    $PY cases/abp_rod_2d.py --report     # 리포트만
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bdbot import Q, checks as C, metrics as MET, physical as PH, report as R  # noqa: E402
from bdbot import nondim as ND, run as RUN, runid as RID, scales as SC  # noqa: E402
from bdbot import sim as SIM  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CASE = ROOT / "intake/abp-rod-2d-run-flip"
f_ = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)

# 표본 간격 사이 방향각이 이 이상 뛰면 텀블로 본다 (열적 회전확산만으로는 한 표본
# 간격 안에 이렇게 못 튄다 — Dr*sample_dt ≪ 1). 저장 시점에 한 번 계산해 npz에 넣는다 —
# 그래야 make_plots든 나중에 노트북에서 열어보든 같은 판정을 본다 (추정을 여러 곳에 안 흩뜬다).
TUMBLE_JUMP_RAD = 1.0


def detect_tumbles(theta, threshold=TUMBLE_JUMP_RAD):
    """theta: (n_t, n_particles). 반환: (n_t-1, n_particles) bool — 표본 간 텀블 추정."""
    dth = np.abs((np.diff(theta, axis=0) + math.pi) % (2 * math.pi) - math.pi)
    return dth > threshold


# ════════════════════════════════════════════════════════════════════════
# ① 물리계 — bdbot.physical 이 읽고 검사한다 (tier·derived_from·형상)
# ════════════════════════════════════════════════════════════════════════
def node(s, *path):
    cur = s.raw
    for k in path:
        cur = cur[k]
    u = cur.get("unit") or "dimensionless"
    return Q(cur["value"], u)


# ════════════════════════════════════════════════════════════════════════
# ② 스케일 표 — ★ 케이스 고유: 타원체라 길이가 3개(장·단·등가), 시간에 텀블이 들어온다
# ════════════════════════════════════════════════════════════════════════
def build_ledger(s, dt=None, T_obs=None):
    """원장. ★ `dt`·`T_obs` 도 원장에 올립니다 (bdbot.scales.MANDATORY_ROLES).

    이 계는 dt 를 `numerics.dt` 로 **직접 받습니다** — 1-A(편향 역산)·1-B(τ_int 역산)와
    달리 사람이 정한 값입니다. 원장에 올려야 시간척도 정렬에서 위치가 보입니다.
    """
    d = node(s, "particle", "diameter").to("m")                 # 등가부피 구 지름 (기준)
    semi = node(s, "particle", "semi_axes").to("m")
    kT = node(s, "derived_scales", "kT").to("J") if "kT" in s.raw["derived_scales"] \
        else (Q(1.380649e-23, "J/K") * node(s, "medium", "temperature")).to("J")
    gbar = node(s, "friction", "gamma_bar_2d").to("kg/s")
    gr = node(s, "friction", "gamma_rot_z").to("kg*m^2/s")
    v = node(s, "active", "speed").to("m/s")
    tau_tumble = node(s, "active", "tumble_interval").to("s")
    L = node(s, "geometry", "box_length").to("m")

    D_bar = (kT / gbar).to("m^2/s")
    tau_B = (d**2 / D_bar).to("s")
    D_r = (kT / gr).to("1/s")
    tau_r = (1 / D_r).to("s")
    # ★ 독립 과정 가정: 감쇠율이 더해진다
    tau_eff = (1 / (D_r + 1 / tau_tumble)).to("s")
    l_p = (v * tau_eff).to("m")
    tau_v = (d / v).to("s")
    rho = node(s, "particle", "density").to("kg/m^3")
    m = (rho * (4 * math.pi / 3) * semi[0] * semi[1] * semi[2]).to("kg")
    tau_p = (m / gbar).to("s")

    lg = SC.ScaleLedger()
    lg.add_length("2a_minor", (2 * semi[1]).to("m"), "단축 길이")
    lg.add_length("d_eq", d, "등가부피 구 지름 (기준)")
    lg.add_length("l_p", l_p, "v·τ_eff 지속길이", star=True)
    lg.add_length("2a_major", (2 * semi[0]).to("m"), "장축 길이")
    lg.add_length("L", L, "box", role="box")
    lg.add_time("tau_p", tau_p, "m/γ̄ 관성", role="inertia")
    if dt is not None:
        lg.add_time("dt", dt, "적분 스텝 (numerics.dt, 사람이 정함)", role="dt")
    lg.add_time("tau_v", tau_v, "d_eq/v 이류", star=True)
    lg.add_time("tau_eff", tau_eff, "방향 상관", star=True)
    lg.add_time("tau_tumble", tau_tumble, "텀블 간격")
    lg.add_time("tau_r", tau_r, "1/D_r 열적 회전확산")
    lg.add_time("tau_B", tau_B, "d_eq²/D̄ 확산 (기준)")
    if T_obs is not None:
        lg.add_time("T_obs", T_obs, "observation window", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("fa_d", (gbar * v * d).to("J"), "f_a·d_eq 자기추진 일")
    lg.derived = dict(d=d, kT=kT, gbar=gbar, gr=gr, v=v, D_bar=D_bar, tau_B=tau_B,
                      D_r=D_r, tau_r=tau_r, tau_eff=tau_eff, tau_tumble=tau_tumble,
                      tau_v=tau_v, tau_p=tau_p, l_p=l_p, L=L, semi=semi)
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " 기준 길이는 등가부피 구 지름 d_eq — 타원체는 장축/단축이 "
        "따로 있으나 무차원화 기준은 하나여야 한다. 병진 마찰은 BD 제약상 등방 평균 γ̄ "
        "(조화평균 — 장시간 MSD가 D̄로 정해지도록).",
        length_symbol="d_eq")                     # ★ 원장 기호와 일치시켜야 한다
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ③ 무차원수 + ④ 분리 검사
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(s, lg, N):
    """무차원수와 분리 검사. dt·T_obs 는 원장에서 읽습니다."""
    D = lg.derived
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    Pe = f_(D["v"] * D["d"] / D["D_bar"])
    # Pe 는 τ_B/τ_v 와 **같은 값**입니다 (v d/D̄ = (d²/D̄)/(d/v)) — 원장으로 교차검증됩니다.
    groups = [
        ND.Group("Pe", Pe, ("times", "tau_B"), ("times", "tau_v"),
                 "v d_eq/D̄ = τ_B/τ_v", "advection vs diffusion"),
        ND.Group("D_r*", f_(D["D_r"] * D["tau_B"]), ("times", "tau_B"), ("times", "tau_r"),
                 "D_r τ_B", "회전 vs 병진"),
        ND.Group("l_p/d_eq", f_(D["l_p"] / D["d"]), ("lengths", "l_p"), ("lengths", "d_eq"),
                 "", "지속길이"),
        ND.Group("p", f_(D["semi"][0] / D["semi"][1]), ("lengths", "2a_major"),
                 ("lengths", "2a_minor"), "a_maj/a_min", "종횡비"),
        ND.Group("tau_tumble/tau_r", f_(D["tau_tumble"] / D["tau_r"]),
                 ("times", "tau_tumble"), ("times", "tau_r"), "", "★ 텀블 vs 회전확산"),
        ND.Group("L/d_eq", f_(D["L"] / D["d"]), ("lengths", "L"), ("lengths", "d_eq"),
                 "", "box"),
        ND.Group("T_obs/tau_eff", f_(T_obs / D["tau_eff"]), ("times", "T_obs"),
                 ("times", "tau_eff"), "", "observation window"),
        ND.Group("St", f_(D["tau_p"] / D["tau_B"]), ("times", "tau_p"), ("times", "tau_B"),
                 "tau_p/tau_B", "inertia vs diffusion"),
    ]
    ck = [
        C.Check("model", "관성 무시     τ_p/τ_v", f_(D["tau_p"] / D["tau_v"]), C.GATE, "<=",
                "τ_dyn = 관심 최속 척도 = τ_v (이류). BD 타당성, dt와 무관"),
        C.Check("integration", "advection resolved   dt/tau_v", f_(dt / D["tau_v"]), C.GATE, "<=",
                f"한 스텝에 d_eq의 {100*f_(dt/D['tau_v']):.1f}% 이동"),
        C.Check("integration", "회전 해상     dt·D_r", f_(dt * D["D_r"]), C.GATE, "<=",
                "열적 회전확산 해상"),
        C.Check("integration", "텀블 해상     dt/τ_tumble", f_(dt / D["tau_tumble"]), C.GATE, "<=",
                "포아송 근사 (p = dt/τ_tumble ≪ 1)"),
        C.Check("geometry", "유한크기     ℓ_p/(L/4)", f_(D["l_p"] / (D["L"] / 4)), 1.0, "<=",
                "지속길이가 박스의 1/4 이내 (액티브 인공효과)", hard=False),
        C.Check("statistics", "관측창       T_obs/τ_eff", f_(T_obs / D["tau_eff"]), 100.0, ">=",
                f"방향 상관시간 기준. 독립 입자 {N}개가 통계 배수", hard=False),
    ]
    return groups, ck


def report_blocks(s, lg, n_eq, n_prod, N, sample_every):
    D = lg.derived
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    inp = [
        R.kv("semi_axes", "(1.0,0.25,0.25) µm", 1, "사용자 확정 2026-08-04", val_w=20),
        R.kv("medium", "water @300K", 1, "사용자 확정(매질)+1-A 승계(온도)", val_w=20),
        R.kv("v", f"{D['v'].to('um/s'):~.3gP}", 0, "sketch 'v ≤ 5 µm/s' 상한", val_w=20),
        R.kv("tau_tumble", f"{D['tau_tumble']:~.3gP}", 1, "사용자 확정", val_w=20),
        R.kv("N", str(N), 3, "★제안 독립 앙상블", val_w=20),
        R.kv("tumble", "균등 무작위 (2D)", 3, "★제안 (스케치는 'flip'=180°)", val_w=20),
    ]
    der = [
        f"  γ̄(2D)  = {D['gbar']:~.4eP}   D̄ = {D['D_bar'].to('um^2/s'):~.4fP}  (Perrin, 구 극한 검증)",
        f"  γ_r,z   = {D['gr']:~.4eP}   D_r = {D['D_r']:~.4fP}   τ_r = {D['tau_r']:~.4fP}",
        f"  ★ 실제 마찰은 이방(ζ⊥/ζ∥=1.287)이나 BD는 등방만 — 단시간 MSD 이방성 손실",
        f"  τ_eff   = 1/(D_r + 1/τ_tumble) = {D['tau_eff']:~.4fP}   ← 독립 과정 가정",
    ]
    plan = [
        f"  dt      = {dt.to('ms'):~.4gP}  = {f_(dt/D['tau_B']):.3e} τ_B",
        f"  T_obs   = {T_obs:~.4gP} = {f_(T_obs/D['tau_eff']):.0f} τ_eff",
        f"  steps   = eq {n_eq:,} + prod {n_prod:,}  (표본 {n_prod//sample_every:,}개)  × N={N}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# ⑤ 해석해 — ★ 케이스 고유. 2D run-and-tumble + 열적 회전확산
# ════════════════════════════════════════════════════════════════════════
def analytic(lg):
    D = lg.derived
    tau_eff, v, D_bar = D["tau_eff"], D["v"], D["D_bar"]
    return {
        "tau_eff": tau_eff,                              # ⟨cos Δθ⟩ 감쇠시간
        "D_eff": (D_bar + v**2 * tau_eff / 2).to("um^2/s"),   # 2D 장시간
        "D_bar": D_bar.to("um^2/s"),                     # 단시간
        "ratio": f_((D_bar + v**2 * tau_eff / 2) / D_bar),
    }


def msd_analytic(t, D_bar, v, tau):
    """2D: MSD = 4D̄t + 2v²τ²(t/τ − 1 + e^{−t/τ})"""
    x = t / tau
    return 4 * D_bar * t + 2 * v**2 * tau**2 * (x - 1 + np.exp(-x))


def analytic_star(Pe, Dr_star, tau_tumble_star) -> dict:
    """`analytic(lg)` 와 같은 식을 무차원(σ=d_eq, τ=τ_B)으로. D̄*=1 은 정의상 참이다 —
    D̄ = kT/γ̄ 가 곧 기준 확산(τ_B := d²/D̄)이므로 별도 유도가 아니라 항등식이다.
    """
    tau_eff_star = 1.0 / (Dr_star + 1.0 / tau_tumble_star)
    D_bar_star = 1.0
    D_eff_star = D_bar_star + Pe ** 2 * tau_eff_star / 2
    return {"tau_eff_star": tau_eff_star, "D_bar_star": D_bar_star,
            "D_eff_star": D_eff_star, "ratio": D_eff_star / D_bar_star}


# ════════════════════════════════════════════════════════════════════════
# ⑥ L4 — 스펙만 읽고 계를 세운다 (bdbot.run 이 돌리고 판정한다)
#    ★ 케이스 고유: 액티브 힘 + 열적 회전확산 updater + run-and-tumble updater
# ════════════════════════════════════════════════════════════════════════
@RUN.builder("abp-rod-2d-run-flip")
def build(spec, outdir=None) -> RUN.Build:
    """스펙 → HOOMD 계.

    ★ 이 계는 페어 상호작용도 트랩도 없어서 `pe_per_particle` 이 가리킬 퍼텐셜
      에너지가 없습니다(활성힘은 비보존력). 대신 원래(1-C 이전) 구현이 이미 평형
      지표로 썼던 `⟨cos Δθ⟩`(표본 간 방향 상관)을 씁니다 — 회전확산이 죽으면
      1.0 에서 안 변하고(정확히 상수 → FROZEN), 정상이면 매 호출마다 다른 노이즈
      실현이라 절대 정확히 상수가 되지 않습니다.
    """
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    N, Pe = int(P["N"]), float(P["Pe"])
    Dr_star = float(P["Dr_star"])
    tau_tumble_star, L_star = float(P["tau_tumble_star"]), float(P["L_star"])
    dt_star, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_eq, n_prod = int(Nm["n_eq"]), int(Nm["n_prod"])
    sample_every = int(Nm["sample_every"])

    rng = np.random.default_rng(seed)
    n_side = int(math.ceil(math.sqrt(N)))
    a = L_star / n_side
    pos = np.array([[(i % n_side + .5) * a - L_star / 2,
                     (i // n_side + .5) * a - L_star / 2, 0.0] for i in range(N)])
    fr = SIM.frame_2d(pos, L_star, orientation=True)
    th0 = rng.uniform(0, 2 * math.pi, N)
    fr.particles.orientation = np.c_[np.cos(th0 / 2), np.zeros(N), np.zeros(N),
                                     np.sin(th0 / 2)]
    sim = SIM.make_sim(fr, seed=seed)

    active = md.force.Active(filter=hoomd.filter.All())
    active.active_force["A"] = (Pe, 0.0, 0.0)      # 무차원: f_a* = Pe (입자 로컬 프레임)
    active.active_torque["A"] = (0.0, 0.0, 0.0)
    SIM.attach_brownian(sim, dt_star, [active])   # γ=1, kT=1 (thermal 규약)

    # 열적 회전확산 — 적분기가 아니라 updater로 (bd-hoomd 함정 3)
    sim.operations.updaters.append(
        active.create_diffusion_updater(trigger=hoomd.trigger.Periodic(1),
                                        rotational_diffusion=Dr_star))

    # ★ run-and-tumble: 포아송 과정으로 방향을 **균등 무작위 재배향**
    #   bd-hoomd 의 run-and-flip 스니펫과 다르다 — 그건 180° 반전이고 이건 완전 무작위.
    class RunAndTumble(hoomd.custom.Action):
        def __init__(self, rate, dt, seed=7):
            self.p = rate * dt                      # 스텝당 텀블 확률 (≪1 이어야 함)
            self.rng = np.random.default_rng(seed)
            self.n_tumbles = 0
            self.n_steps = 0

        def act(self, timestep):
            self.n_steps += 1
            with self._state.cpu_local_snapshot as snap:
                q = np.array(snap.particles.orientation, copy=True)
                hit = self.rng.random(len(q)) < self.p
                k = int(hit.sum())
                if k:
                    self.n_tumbles += k
                    th = self.rng.uniform(0, 2 * math.pi, k)
                    q[hit] = np.c_[np.cos(th / 2), np.zeros(k), np.zeros(k), np.sin(th / 2)]
                    snap.particles.orientation[:] = q

    tumble = RunAndTumble(1.0 / tau_tumble_star, dt_star, seed=seed + 1)
    sim.operations.updaters.append(
        hoomd.update.CustomUpdater(action=tumble, trigger=hoomd.trigger.Periodic(1)))
    gsd = (Path(outdir) / "traj_A.gsd") if outdir else None
    SIM.add_trajectory_writer(sim, gsd, max(1, n_prod // 200))

    # ── 표본 누적기 — 원래(1-C 이전) 구현과 같은 메모리 사용량을 유지한다 ──
    #   sample()이 매 표본 전체 N-입자의 위치·방향을 그대로 cols로 돌려주면 그게
    #   그대로 (n_samp,N,...)로 쌓여 observables.npz에 저장된다. 원래도 MSD·상관함수는
    #   전체 N으로 계산했지만 **저장은 대표 입자 8개뿐**이었다(1급 산출물 궤적) — 그
    #   비대칭을 유지하려고 전체 배열은 클로저에, npz에 갈 것만 finalize에서 슬라이스한다.
    n_samp = n_prod // sample_every if sample_every else 0
    xy = np.empty((max(n_samp, 1), N, 2), dtype=np.float64)     # 언랩 위치 (MSD용)
    th_arr = np.empty((max(n_samp, 1), N), dtype=np.float64)    # 방향각
    i_sample = [0]
    prev_theta = [None]

    def pe_pp():
        """⟨cos Δθ⟩ (직전 호출 대비) — 이 계의 유일한 '항상 변해야 하는' 신호."""
        sn = sim.state.get_snapshot()
        q = np.array(sn.particles.orientation)
        theta = 2 * np.arctan2(q[:, 3], q[:, 0])
        if prev_theta[0] is None:
            prev_theta[0] = theta
            return 1.0
        val = float(np.cos(theta - prev_theta[0]).mean())
        prev_theta[0] = theta
        return val

    def sample(timestep, phase):
        sn = sim.state.get_snapshot()
        p = np.array(sn.particles.position)[:, :2]
        img = np.array(sn.particles.image)[:, :2]
        i = i_sample[0]
        if i < n_samp:
            xy[i] = p + img * L_star                        # 언랩 — 주기 wrap 제거
            q = np.array(sn.particles.orientation)
            th_arr[i] = 2 * np.arctan2(q[:, 3], q[:, 0])     # z축 회전각
        i_sample[0] = i + 1
        return {}

    def finalize(cols):
        from scipy import optimize

        n = i_sample[0]
        XY, TH = xy[:n], th_arr[:n]
        dt_sample = dt_star * sample_every
        v_star = Pe

        n_t = len(XY)
        lags = np.unique(np.round(np.logspace(0, math.log10(n_t // 2), 40)).astype(int))
        msd, cth, msad = [], [], []
        for L in lags:
            d = XY[L:] - XY[:-L]
            msd.append(float((d ** 2).sum(axis=2).mean()))
            dth = TH[L:] - TH[:-L]
            cth.append(float(np.cos(dth).mean()))
            # MSAD: 텀블이 최대 π 점프를 만들어 언랩이 모호하다 → [-π,π]로 접어 계산하고
            #       "접힌 MSAD"임을 명시한다. 방향 상관 C(t)가 모호하지 않은 관측량이다.
            msad.append(float((((dth + math.pi) % (2 * math.pi) - math.pi) ** 2).mean()))
        t = lags * dt_sample
        msd, cth, msad = np.array(msd), np.array(cth), np.array(msad)

        # ★ 변위 상관 C(t) = ⟨Δr(t₀)·Δr(t₀+t)⟩ / ⟨Δr²⟩
        #   BD는 과감쇠라 **속도에 물리적 의미가 없다** (bd-hoomd 함정 5) → 속도 자기상관 금지.
        #   대신 고정 간격 δ 의 변위끼리 상관을 본다 — 이산 시간의 정당한 대응물이고
        #   run-and-tumble 이면 exp(−t/τ_eff) 로 감쇠해야 한다 (또 하나의 독립 검증).
        dr = XY[1:] - XY[:-1]                       # (n_t-1, N, 2) 연속 변위
        dlags = np.unique(np.round(np.logspace(
            0, math.log10(max(2, len(dr) // 3)), 28)).astype(int))
        norm = float((dr * dr).sum(axis=2).mean())
        cdr = np.array([float((dr[L:] * dr[:-L]).sum(axis=2).mean()) / norm
                        for L in dlags])
        disp_corr_t = dlags * dt_sample

        # τ_eff: ⟨cos Δθ⟩ 를 지수로 피팅 (표본평균 빼지 않음 — 참 평균이 0이 아님)
        m = cth > 0.05
        if m.sum() >= 4:
            popt, _ = optimize.curve_fit(lambda tt, A, tau: A * np.exp(-tt / tau),
                                         t[m], cth[m],
                                         p0=[1.0, t[m][len(t[m]) // 2]], maxfev=20000)
            tau_eff_fit, cos_A = float(abs(popt[1])), float(popt[0])
        else:
            tau_eff_fit = cos_A = float("nan")

        # ★ MSD 는 **해석식 전체를 2모수 피팅**한다 (D̄, τ 자유) — bd-physics 교훈 참조
        #   (단시간·장시간 분리 피팅은 활성 오염으로 D̄가 +37% 틀렸다).
        def msd_model(tt, D_bar, tau):
            x = tt / tau
            return 4 * D_bar * tt + 2 * v_star ** 2 * tau ** 2 * (x - 1 + np.exp(-x))

        tau0 = tau_eff_fit if np.isfinite(tau_eff_fit) else t[len(t) // 3]
        try:
            popt, _ = optimize.curve_fit(msd_model, t, msd, p0=[1.0, tau0],
                                         sigma=msd, absolute_sigma=False, maxfev=40000)
            D_bar_fit, tau_msd_fit = float(popt[0]), float(abs(popt[1]))
            msd_resid_pct = float(100 * np.sqrt(np.mean(
                ((msd - msd_model(t, *popt)) / msd) ** 2)))
        except Exception as e:
            D_bar_fit = tau_msd_fit = msd_resid_pct = float("nan")
            print(f"    (MSD 피팅 실패: {e})")
        tail = t > 5 * tau0
        D_eff_fit = (float(np.polyfit(t[tail], msd[tail], 1)[0] / 4)
                    if tail.sum() >= 3 else float("nan"))
        head = t < 0.1 * tau0
        D_short_naive = (float(np.polyfit(t[head], msd[head], 1)[0] / 4)
                        if head.sum() >= 3 else float("nan"))
        if head.sum():
            th_ = float(t[head][-1]); xh = th_ / tau0
            head_contam_pct = float(
                100 * 2 * v_star ** 2 * tau0 ** 2 * (xh - 1 + math.exp(-xh))
                / (4 * D_bar_fit * th_))
        else:
            head_contam_pct = float("nan")

        rate_measured = tumble.n_tumbles / (tumble.n_steps * dt_star * N)
        ana = analytic_star(Pe, Dr_star, tau_tumble_star)

        # ★ 역할을 정직하게 붙인다 (bdbot.metrics.ROLES). 다섯 개 전부 내가 구현한
        #   모델(활성힘+회전확산 updater+포아송 텀블+BD)에서 해석적으로 따라 나온다 —
        #   implementation_check. 일치 = 코드 검증. 이 케이스에는 가설 검증이 없다.
        ROLE, TOL = "implementation_check", 5.0
        rows = [
            ("τ_eff (⟨cosΔθ⟩ 피팅)", tau_eff_fit, ana["tau_eff_star"], 0, 1, "s", ROLE, TOL,
             "1/τ_eff = D_r + 1/τ_tumble — 두 과정을 독립으로 **구현했으니** 따라 나온다"),
            ("D_eff (MSD 장시간)", D_eff_fit, ana["D_eff_star"], 2, -1, "um^2/s", ROLE, TOL,
             "D_eff = D̄ + v²τ_eff/2 — 같은 모델의 결과"),
            ("D̄ (MSD 전체 피팅)", D_bar_fit, ana["D_bar_star"], 2, -1, "um^2/s", ROLE, TOL,
             "kT/γ̄ — 적분기에 넣은 값이 그대로 나오는지 (가장 순환적)"),
            ("τ (MSD 전체 피팅)", tau_msd_fit, ana["tau_eff_star"], 0, 1, "s", ROLE, TOL,
             "⟨cosΔθ⟩ 와 MSD 두 경로가 같은 τ 를 주는지 — 내부 일관성"),
            ("텀블 빈도", rate_measured, 1.0 / tau_tumble_star, 0, -1, "1/s", ROLE, 3.0,
             "updater가 지정한 rate 를 내는지"),
        ]
        obs = []
        for name, meas_star, pred_star, L, T, u, role, tol, note in rows:
            meas = float(spec.physical(meas_star, L=L, T=T).to(u).magnitude)
            pred = float(spec.physical(pred_star, L=L, T=T).to(u).magnitude)
            # scope=composite + implementation_check 이므로 derivation 필수 (원칙 9.1) —
            # note 가 곧 유도 근거고, 그게 전부 "내가 넣은 것이 그대로 나오는가" 라는
            # 사실이 원칙 8이 지적한 순환성이다 (가설 검증 0건).
            obs.append(MET.observable(name, meas, pred, u, "analytic_from_model",
                                      role=role, tol_pct=tol, note=note,
                                      scope="composite", derivation=note))

        n_tr = min(8, N)
        traj_theta = TH[:, :n_tr]
        tumble_mask = detect_tumbles(traj_theta)
        cos_series = np.cos(TH[1:] - TH[:-1]).mean(axis=1)
        return {"observables": obs,
                "extra": {"tumble_rate_star": rate_measured, "n_tumbles": tumble.n_tumbles,
                          "D_eff_over_Dbar": D_eff_fit / D_bar_fit,
                          "ratio_model": ana["ratio"], "msd_resid_pct": msd_resid_pct,
                          "D_short_naive_star": D_short_naive,
                          "head_contam_pct": head_contam_pct,
                          "n_particles_traced": n_tr, "n_particles_total": N,
                          "sample_dt_star": dt_sample,
                          "tumble_jump_threshold_rad": TUMBLE_JUMP_RAD},
                "arrays": {"t": t, "msd": msd, "cos_theta": cth, "msad_folded": msad,
                          "cos_theta_series": cos_series,
                          "disp_corr_t": disp_corr_t, "disp_corr": cdr,
                          "traj_xy": XY[:, :n_tr, :], "traj_theta": traj_theta,
                          "tumble_mask": tumble_mask,
                          "final_xy": XY[-1], "final_theta": TH[-1]}}

    return RUN.Build(
        sim=sim, forces=[active], n_particles=N,
        sample=sample, pe_per_particle=pe_pp,
        n_eq=n_eq, n_prod=n_prod, sample_every=sample_every,
        tags=["2D", "active", "run_and_tumble", "ellipsoid", "newtonian",
             "no_pair_interaction", "MSD"],
        physical={"N": N, "Pe": Pe, "Dr_star": Dr_star, "tau_tumble_star": tau_tumble_star},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true",
                    help="L3 스펙을 specs/<run_id>.json 으로 쓰고 종료")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--N", type=int, default=None)
    ap.add_argument("--tobs", type=float, default=None, help="T_obs/τ_eff")
    args = ap.parse_args()

    s = PH.load(CASE)
    if s.errors:
        print(PH.render_check(s))
        return 1
    lg0 = build_ledger(s)                      # dt·T_obs 를 정하려면 τ_eff 가 먼저 필요
    D = lg0.derived
    N = args.N or int(node(s, "particle", "count").magnitude)
    dt = node(s, "numerics", "dt").to("s")
    T_obs = Q(s.raw["numerics"]["production_s"], "s")
    if args.tobs:
        T_obs = args.tobs * D["tau_eff"]
    if args.smoke:
        N, T_obs = min(N, 200), 30 * D["tau_eff"]
    lg = build_ledger(s, dt=dt, T_obs=T_obs)   # 원장 완성 (필수 역할 전부 채움)
    D = lg.derived

    n_prod = int(round(f_(T_obs / dt)))
    n_eq = max(1, int(round(f_(5 * D["tau_eff"] / dt))))     # 방향 상관 5배로 평형화
    sample_every = max(1, n_prod // 3000)
    n_prod = (n_prod // sample_every) * sample_every

    groups, ck = analyze_scales(s, lg, N)
    # ★ L4(build())는 스펙만 읽습니다 — Dr_star·tau_tumble_star·L_star 도 여기 넣어야
    #   케이스 YAML을 다시 읽지 않고 계를 세울 수 있습니다 (L2↔L4 계약은 스펙 하나).
    spec = ND.NondimSpec(
        case=s.label, system=s.raw, reference=lg.ref, ledger=lg, groups=groups, checks=ck,
        params={"N": N, "Pe": groups[0].value,
                "Dr_star": f_(D["D_r"] * D["tau_B"]),
                "tau_tumble_star": f_(D["tau_tumble"] / D["tau_B"]),
                "L_star": f_(D["L"] / D["d"])},
        numerics={"dt_star": f_(dt / D["tau_B"]), "dt_over_tau_B": f_(dt / D["tau_B"]),
                  "n_eq": n_eq, "n_prod": n_prod, "sample_every": sample_every,
                  "seed": 20260804},
        tag="smoke" if args.smoke else None, nhex=12)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(s, lg, n_eq, n_prod, N, sample_every)
    report, verdict = R.render(
        title=f"DimensionlessReport — {s.label}   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=ck,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)
    if spec.errors:
        print(f"\n❌ L3 무결성 오류 {len(spec.errors)}건 — 무차원화가 성립하지 않습니다.")
        return 1
    if verdict == "FAIL":
        print("\n❌ 하드 분리 검사 실패 — 실행하지 않습니다.")
        return 1
    p = spec.write(ROOT / "specs" / f"{run_id}.json")
    if args.spec or args.report:
        if args.spec:
            print(f"\nL3 스펙: {p.relative_to(ROOT)}")
        return 0

    # ── L4 — 디스크의 스펙을 되읽어서 실행한다 (bdbot.run, 해시 검증이 그때 걸린다) ──
    outdir = ROOT / "runs" / run_id
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    print(RUN.render_verdict(v))
    if v["status"] == "skipped":
        return 0
    if v["status"] != RUN.OK:
        return 1
    # ★ report.txt는 execute() 이후에 쓴다 — 먼저 쓰면 prepare_outdir가 지운다.
    (outdir / "report.txt").write_text(report)

    # ── 해석해 대조 — finalize()가 이미 계산해 metrics.json에 넣어둔 것을 표로 찍는다 ──
    m = json.loads((outdir / "metrics.json").read_text())
    obs_out = m["observables"]
    res_extra = m.get("result", {})
    ana = analytic(lg)

    verdict_o, bad_impl, dev_hypo, meas = MET.judge(obs_out)
    lines = ["", "=" * R.W, "결과 — 물리 단위로 역변환 후 대조", "=" * R.W,
             f"{'관측량':<24}{'측정':>15}{'예측':>15}{'오차':>9}  역할        판정"]
    for o in obs_out:
        mark = "✓" if o not in bad_impl and o not in dev_hypo else "✗"
        r_short = {"implementation_check": "구현검사", "hypothesis": "가설",
                   "measurement": "측정"}[o["role"]]
        lines.append(f"{o['name']:<24}{o['measured']:>15.6g}{o['predicted']:>15.6g}"
                     f"{o['err_pct']:>+8.2f}%  {r_short:<10}  {mark}")
    d_um, tau_B_s = f_(D["d"] / Q(1, "um")), f_(D["tau_B"] / Q(1, "s"))
    d_short = res_extra["D_short_naive_star"]
    lines += ["",
              f"  D_eff/D̄  측정 {res_extra['D_eff_over_Dbar']:.3f}"
              f"  vs 모델 {res_extra['ratio_model']:.3f}   (활성이 확산을 몇 배 키웠나)",
              f"  MSD 전체 피팅 잔차 RMS {res_extra['msd_resid_pct']:.2f}%",
              (f"  ⚠ naive 단시간 기울기 D̄={d_short*d_um**2/tau_B_s:.4f}"
               f" µm²/s — 그 창의 활성 오염 {res_extra['head_contam_pct']:.0f}%. 쓰면 안 된다"
               if math.isfinite(d_short)
               else "  (표본 간격이 0.1τ_eff 보다 커서 naive 단시간 창은 아예 없음)"),
              "",
              "★ 이 런이 확인한 것과 확인하지 못한 것",
              "  확인함  : 구현이 의도한 모델과 일치한다 (활성힘·회전확산·텀블·BD 적분).",
              "            다섯 예측이 전부 그 모델에서 유도되므로, 일치는 **코드 검증**이다.",
              "  확인 못함: 이 계의 실제 물리. 예측이 시뮬레이션의 가정 위에 서 있어서",
              "            '다를 수 있는' 것을 시험하지 않았다 — 가설 검증이 0건이다.",
              "  가설이 생기려면 시뮬레이션이 부과하지 않은 가정이 있어야 한다. 예:",
              "    · 유한 밀도로 올려 상호작용을 넣는다 → 희박 D_eff 공식이 깨지는지",
              "    · 가둠(벽·채널)을 넣는다 → 축적·유영이 나오는지",
              "    · 병진 이방성 → **BD로는 불가**(HI 부재). MPCD 등 다른 방법이 필요",
              "=" * R.W,
              f"VERDICT: {verdict_o}",
              "=" * R.W]
    all_ok = not bad_impl
    result = "\n".join(lines)
    print(result)
    (outdir / "result.txt").write_text(report + "\n" + result)

    make_plots(lg, ana, outdir)
    print("\n".join(RID.list_artifacts(outdir, ROOT)))
    return 0 if all_ok else 1


def make_plots(lg, ana, outdir):
    """① 단일 입자 궤적  ② MSD  ③ 상관함수 2종  ④ 접힌 MSAD"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

    D = lg.derived
    tau_B = f_(D["tau_B"] / Q(1, "s"))
    d_um = f_(D["d"] / Q(1, "um"))
    te = f_(ana["tau_eff"] / Q(1, "s"))
    Db = f_(ana["D_bar"] / Q(1, "um^2/s"))
    v_um = f_(D["v"] / Q(1, "um/s"))
    z = np.load(outdir / "observables.npz")
    obs = z
    res_extra = json.loads((outdir / "metrics.json").read_text())["result"]
    msd_resid_pct = res_extra["msd_resid_pct"]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))

    # ── ① 단일 입자 궤적 ─────────────────────────────────────────────
    tr = z["traj_xy"] * d_um            # (n_t, n_tr, 2) [µm], 언랩
    tumble_mask = z["tumble_mask"]      # 저장 시점에 확정된 텀블 판정 (detect_tumbles) — 재추정 안 함
    dts = float(res_extra["sample_dt_star"]) * tau_B
    cmap = plt.cm.turbo(np.linspace(0.08, 0.92, tr.shape[1]))
    for k in range(tr.shape[1]):
        xy = tr[:, k, :] - tr[0, k, :]
        ax[0, 0].plot(xy[:, 0], xy[:, 1], "-", lw=0.9, color=cmap[k], alpha=.85)
        jump = np.where(tumble_mask[:, k])[0]
        ax[0, 0].plot(xy[jump + 1, 0], xy[jump + 1, 1], ".", ms=3.5,
                      color=cmap[k], alpha=.7)
    ax[0, 0].plot(0, 0, "k+", ms=9, mew=1.6)
    ax[0, 0].set(xlabel="x [µm]", ylabel="y [µm]", aspect="equal",
                 title=f"① 단일 입자 궤적 {tr.shape[1]}개 (점 = 텀블 추정 지점)")
    ax[0, 0].grid(alpha=.25)
    ax[0, 0].text(.02, .98, f"ℓ_p = {f_(D['l_p']/Q(1,'um')):.2f} µm\n"
                            f"표본 간격 {dts*1e3:.0f} ms",
                  transform=ax[0, 0].transAxes, va="top", fontsize=8,
                  bbox=dict(fc="w", alpha=.7, lw=0))

    # ── ② MSD ────────────────────────────────────────────────────────
    t_s = obs["t"] * tau_B
    ax[0, 1].loglog(t_s, obs["msd"] * d_um**2, "o", ms=4.5, label="측정")
    tt = np.logspace(math.log10(t_s[0]), math.log10(t_s[-1]), 300)
    ax[0, 1].loglog(tt, msd_analytic(tt, Db, v_um, te), "-", lw=2,
                    label="해석해 (2D run-and-tumble)")
    ax[0, 1].loglog(tt, 4 * Db * tt, "--", lw=1.2, c="gray",
                    label=r"단시간 $4\bar{D}t$")
    ax[0, 1].loglog(tt, 4 * f_(ana["D_eff"] / Q(1, "um^2/s")) * tt, ":", lw=1.4,
                    c="tab:red", label=r"장시간 $4D_{eff}t$")
    ax[0, 1].axvline(te, ls=":", c="k", alpha=.6)
    ax[0, 1].text(te, ax[0, 1].get_ylim()[0] * 2, r" $\tau_{eff}$", fontsize=9)
    ax[0, 1].set(xlabel="t [s]", ylabel=r"MSD [µm²]",
                 title=f"② MSD — 잔차 RMS {msd_resid_pct:.2f}%")
    ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=.3, which="both")

    # ── ③ 상관함수 2종 ───────────────────────────────────────────────
    ax[1, 0].semilogy(t_s, np.clip(obs["cos_theta"], 1e-3, None), "o", ms=4.5,
                      label=r"방향 $\langle\cos\Delta\theta\rangle$")
    dct = z["disp_corr_t"] * tau_B
    ax[1, 0].semilogy(dct, np.clip(z["disp_corr"], 1e-3, None), "s", ms=4.5,
                      mfc="none", label=r"변위 $\langle\Delta r\cdot\Delta r\rangle$ (정규화)")
    ax[1, 0].semilogy(tt, np.exp(-tt / te), "-", lw=2,
                      label=f"exp(−t/τ_eff), τ={te:.3f} s")
    ax[1, 0].set(xlabel="t [s]", ylabel="상관 (정규화)", ylim=(1e-3, 2),
                 xlim=(0, min(6 * te, t_s[-1])),
                 title="③ 상관함수 — 두 경로가 같은 τ_eff 를 줘야 한다")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=.3)
    ax[1, 0].text(.98, .95, "★ BD는 과감쇠라 속도에 의미가 없다 (함정 5)\n"
                            "   → 속도 자기상관 대신 **변위 상관**",
                  transform=ax[1, 0].transAxes, ha="right", va="top", fontsize=7.5,
                  bbox=dict(fc="w", alpha=.75, lw=0))

    # ── ④ 접힌 MSAD ──────────────────────────────────────────────────
    ax[1, 1].semilogx(t_s, obs["msad_folded"], "o", ms=4.5, label="측정 (접힘)")
    ax[1, 1].axhline(math.pi**2 / 3, ls="--", c="gray",
                     label=r"균등분포 극한 $\pi^2/3$")
    ax[1, 1].axvline(te, ls=":", c="k", alpha=.6)
    ax[1, 1].set(xlabel="t [s]", ylabel=r"접힌 MSAD [rad²]",
                 title="④ 접힌 MSAD — 텀블 점프로 언랩이 모호")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=.3, which="both")

    for a in ax.ravel():
        a.title.set_fontsize(10)
    fig.suptitle(f"abp-rod-2d — run-and-tumble 타원체 (2µm × 500nm, 물 300K, "
                 f"Pe={f_(D['v']*D['d']/D['D_bar']):.2f})", fontsize=12.5)
    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
