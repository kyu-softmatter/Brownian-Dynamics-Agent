"""L3 무차원화 리포트 — 새로 확정된 3케이스 (실행하지 않음).

사용자 지시 (2026-08-04): "제안값으로 채우고 리포트 먼저 보여줘"

`bdbot` 의 공통 부분(ScaleLedger · Check · report.render)을 그대로 쓴다 — 1-C에서 뽑은
추상화가 **세 케이스 더**에 통하는지 보는 시험대이기도 하다.
케이스마다 다른 것(어떤 스케일이 들어가는지 · 어떤 검사가 필요한지)은 여기서 채운다.

    $PY scratch/nondim_3cases.py [abp|trap-drag|chain]
"""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import Q, checks as C, physical as P, report as R, scales as SC  # noqa: E402

f = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)


def node(s, *path):
    """system.yaml 의 (value, unit) 노드 → pint Quantity."""
    cur = s.raw
    for k in path:
        cur = cur[k]
    u = cur.get("unit", "dimensionless")
    return Q(cur["value"], u if u else "dimensionless")


def ds(s, key):
    return node(s, "derived_scales", key)


# ══════════════════════════════════════════════════════════════════════
def abp_rod():
    s = P.load(ROOT / "intake/abp-rod-2d-run-flip")
    d = node(s, "particle", "diameter").to("m")
    kT = Q(1.380649e-23, "J/K") * node(s, "medium", "temperature")
    gbar = node(s, "friction", "gamma_bar_2d").to("kg/s")
    v = node(s, "active", "speed").to("m/s")
    tau_tumble = node(s, "active", "tumble_interval").to("s")
    tau_r, tau_B, tau_v = ds(s, "tau_r").to("s"), ds(s, "tau_B").to("s"), ds(s, "tau_v").to("s")
    tau_eff, l_p = ds(s, "tau_eff").to("s"), ds(s, "l_p").to("m")
    tau_p, D_r = ds(s, "tau_p").to("s"), ds(s, "D_r").to("1/s")
    L = node(s, "geometry", "box_length").to("m")
    dt = node(s, "numerics", "dt").to("s")
    T_obs = Q(s.raw["numerics"]["production_s"], "s")
    N = int(node(s, "particle", "count").magnitude)
    semi = node(s, "particle", "semi_axes").to("m")

    lg = SC.ScaleLedger()
    lg.lengths = {
        "d_eq     등가부피 구 지름 (기준)": d,
        "2a_minor 단축 길이": 2 * semi[1],
        "2a_major 장축 길이": 2 * semi[0],
        "l_p      v·τ_eff 지속길이 ★": l_p,
        "L        박스": L,
    }
    lg.times = {
        "tau_p    m/γ̄ 관성 이완": tau_p,
        "dt       적분 스텝": dt,
        "tau_v    d_eq/v 이류 ★": tau_v,
        "tau_eff  방향 상관 (텀블+회전확산)": tau_eff,
        "tau_tumble 텀블 간격": tau_tumble,
        "tau_r    1/D_r 열적 회전확산": tau_r,
        "tau_B    d_eq²/D̄ 확산 (기준)": tau_B,
        "T_obs    관측창": T_obs,
    }
    lg.energies = {
        "kT       열에너지 (기준)": kT.to("J"),
        "f_a·d_eq 자기추진 일": (gbar * v * d).to("J"),
    }
    lg.derived = {"gamma": gbar, "d": d, "kT": kT.to("J")}
    lg.ref = SC.thermal_reference(
        d, kT.to("J"), tau_B,
        SC.THERMAL_RATIONALE + " 기준 길이는 등가부피 구 지름 d_eq — 타원체라 장축/단축이 "
        "따로 있으나 무차원화 기준은 하나여야 한다. 병진 마찰은 BD 제약상 등방 평균 γ̄.")
    lg.rationale = lg.ref["rationale"]

    groups = {
        "Pe     = v d_eq/D̄      이류 vs 확산": f(v * d * gbar / kT.to("J")),
        "D_r*   = D_r τ_B       회전 vs 병진": f(D_r * tau_B),
        "l_p/d_eq               지속길이": f(l_p / d),
        "p      = a_major/a_minor 종횡비": f(semi[0] / semi[1]),
        "tau_tumble/tau_r       ★ 텀블 vs 회전확산": f(tau_tumble / tau_r),
        "zeta_perp/zeta_par     실제 마찰 이방성 (BD는 못 냄)": 1.287,
        "L/d_eq                 박스": f(L / d),
        "St     = tau_p/tau_B   관성 vs 확산": f(tau_p / tau_B),
    }
    ck = [
        C.Check("model", "관성 무시     τ_p/τ_v", f(tau_p / tau_v), C.GATE, "<=",
                "τ_dyn = 관심 최속 척도 = τ_v (이류). BD 타당성, dt와 무관"),
        C.Check("integration", "이류 해상     dt/τ_v", f(dt / tau_v), C.GATE, "<=",
                f"한 스텝에 d_eq의 {100*f(dt/tau_v):.1f}% 이동"),
        C.Check("integration", "회전 해상     dt·D_r", f(dt * D_r), C.GATE, "<=",
                "방향 동역학 해상"),
        C.Check("integration", "텀블 해상     dt/τ_tumble", f(dt / tau_tumble), C.GATE, "<=",
                "포아송 텀블 근사가 성립하려면 (bd-hoomd run-and-flip 스니펫)"),
        C.Check("geometry", "유한크기     ℓ_p/(L/4)", f(l_p / (L / 4)), 1.0, "<=",
                "지속길이가 박스의 1/4 이내여야 액티브 인공효과가 없다", hard=False),
        C.Check("statistics", "관측창       T_obs/τ_eff", f(T_obs / tau_eff), 100.0, ">=",
                f"방향 상관시간 기준. 독립 입자 {N}개로 통계 배수 확보", hard=False),
    ]
    inp = [
        R.kv("semi_axes", "(1.0, 0.25, 0.25) µm", 1, "사용자 확정 2026-08-04", val_w=22),
        R.kv("medium", "water @300K", 1, "사용자 확정 (매질) + 1-A 승계 (온도)", val_w=22),
        R.kv("v", f"{v.to('um/s'):~.3gP}", 0, "sketch 'v ≤ 5 µm/s' 상한", val_w=22),
        R.kv("tau_tumble", f"{tau_tumble:~.3gP}", 1, "사용자 확정 '매 0.5초마다 텀블'", val_w=22),
        R.kv("N", str(N), 3, "★제안 독립 앙상블", val_w=22),
        R.kv("tumble_angle", "균등 무작위 (2D)", 3, "★제안 — 스케치는 'flip'(180°)", val_w=22),
    ]
    der = [
        f"  γ̄(2D)   = {gbar:~.4eP}   (Perrin, 구 극한 검증 완료)",
        f"  D̄       = {(kT.to('J')/gbar).to('um^2/s'):~.4fP}",
        f"  γ_r,z    = {node(s,'friction','gamma_rot_z'):~.4eP}   D_r = {D_r:~.4fP}",
        f"  ★ 실제 마찰은 이방(ζ⊥/ζ∥=1.287)이나 BD는 등방만 가능 → 단시간 MSD 이방성 손실",
    ]
    plan = [
        f"  dt      = {dt.to('ms'):~.4gP}  = {f(dt/tau_B):.3e} τ_B",
        f"  T_obs   = {T_obs:~.4gP} = {f(T_obs/tau_eff):.0f} τ_eff",
        f"  steps   = {f(T_obs/dt):,.0f}   × N={N}",
    ]
    return "abp-rod (run-and-tumble 타원체)", s, lg, groups, ck, inp, der, plan


# ══════════════════════════════════════════════════════════════════════
def trap_drag():
    s = P.load(ROOT / "intake/trap-drag-2d-hex300")
    d = node(s, "particle", "diameter").to("m")
    kT = ds(s, "kT").to("J")
    gam = ds(s, "gamma").to("kg/s")
    k_t = node(s, "external", "stiffness").to("N/m")
    v = node(s, "external", "drag_velocity").to("m/s")
    tau_k, tau_int = ds(s, "tau_k").to("s"), ds(s, "tau_int").to("s")
    tau_v, tau_B, tau_p = ds(s, "tau_v").to("s"), ds(s, "tau_B").to("s"), ds(s, "tau_p").to("s")
    l_k, dr_ss = ds(s, "l_k").to("m"), ds(s, "dr_ss").to("m")
    a_mean, a_nn = node(s, "geometry", "a_mean").to("m"), node(s, "geometry", "a_nn").to("m")
    L = node(s, "geometry", "box_length").to("m")
    r_c = node(s, "interactions", 0, "cutoff").to("m")
    dt = node(s, "numerics", "dt").to("s")
    T_obs = Q(s.raw["numerics"]["production_s"], "s")
    N = int(node(s, "particle", "count").magnitude)
    A = f(node(s, "interactions", 0, "amplitude_A"))
    phi = f(node(s, "geometry", "area_fraction"))
    Gamma = A / f(a_mean / d) ** 3

    lg = SC.ScaleLedger()
    lg.lengths = {
        "dr_ss    γv/k 끌림 지연 ★": dr_ss,
        "l_k      √(kT/k) 트랩 요동": l_k,
        "d        입자 지름 (기준)": d,
        "a_mean   평균 간격": a_mean,
        "a_NN     육방 최근접": a_nn,
        "r_c      컷오프": r_c,
        "L        박스": L,
    }
    lg.times = {
        "tau_p    m/γ 관성": tau_p,
        "dt       적분 스텝": dt,
        "tau_k    γ/k_t 트랩 이완 ★": tau_k,
        "tau_int  γ/U''(r_min) 페어": tau_int,
        "tau_v    d/v_x 끌기": tau_v,
        "tau_B    d²/D_t 확산 (기준)": tau_B,
        "T_obs    관측창 (박스 횡단)": T_obs,
    }
    lg.energies = {
        "kT       열에너지 (기준)": kT,
        "Γ·kT     평균간격 페어 결합": (Gamma * kT).to("J"),
        "k_t d²   트랩 강성": (k_t * d**2).to("J"),
    }
    lg.derived = {"gamma": gam, "d": d, "kT": kT}
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " 이 계에는 **두 개의 강성**(트랩 k_t, 페어 U'')이 있고 "
        "트랩이 218배 빠르다 → dt는 트랩이 정한다.")
    lg.rationale = lg.ref["rationale"]

    groups = {
        "Gamma  = A(d/a_mean)³  페어 결합 ★": Gamma,
        "k*     = k_t d²/kT     트랩 vs 열요동": f(k_t * d**2 / kT),
        "phi                    밀집도": phi,
        "Pe     = v d/D_t       이류 vs 확산": f(v * d * gam / kT),
        "SNR    = Δr_ss/ℓ_k     ★ 신호 vs 잡음": f(dr_ss / l_k),
        "r_c/a_mean             이웃 껍질 수": f(r_c / a_mean),
        "tau_k/tau_int          ★ 두 강성의 비": f(tau_k / tau_int),
        "L/d                    박스": f(L / d),
        "St     = tau_p/tau_B": f(tau_p / tau_B),
    }
    ck = [
        C.Check("model", "관성 무시     τ_p/τ_k", f(tau_p / tau_k), C.GATE, "<=",
                "τ_dyn = 최속 관심척도 = τ_k (트랩)"),
        C.Check("integration", "트랩 해상     dt/τ_k", f(dt / tau_k), C.GATE, "<=",
                f"편향 ≈ {C.bias_from_dt(dt, tau_k):.3f}% (선형계 기준)"),
        C.Check("integration", "페어 해상     dt/τ_int", f(dt / tau_int), C.GATE, "<=",
                "페어는 트랩보다 218배 느려 여유가 크다"),
        C.Check("geometry", "컷오프       r_c/(L/2)", f(r_c / (L / 2)), 1.0, "<=",
                "최소 이미지 (bd-hoomd 함정 6). 1-B에서 이게 진짜 게이트였다"),
        C.Check("statistics", "관측창       T_obs/τ_k", f(T_obs / tau_k), 100.0, ">=",
                "트랩 상관시간 기준", hard=False),
        C.Check("statistics", "측정가능성   SNR", f(dr_ss / l_k), 1.0, ">=",
                "★ 끌림 신호가 열요동보다 커야 1회 표본으로 보인다. "
                "SNR<1 이면 평균화로만 볼 수 있다", hard=False),
    ]
    inp = [
        R.kv("d", f"{d.to('um'):~.3gP}", 1, "sketch 'R=5µm' + 1-A 지름 관례", val_w=22),
        R.kv("k_t", f"{k_t.to('pN/um'):~.3gP}", 0, "sketch", val_w=22),
        R.kv("v_x", f"{v.to('um/s'):~.3gP}", 0, "sketch", val_w=22),
        R.kv("N", str(N), 0, "sketch 'N ~ 300'", val_w=22),
        R.kv("A", f"{A:g}", 3, "★제안 — 1-B에서 육방 결정 확인된 값", val_w=22),
        R.kv("phi", f"{phi:g}", 3, "★제안 — 1-B와 동일", val_w=22),
    ]
    der = [
        f"  γ = {gam:~.4eP}   D_t = {(kT/gam).to('um^2/s'):~.4fP}",
        f"  Γ = {Gamma:.2f} → 1-B 실측으로 육방 결정 영역 (ψ₆=0.885, 6배위 98.7%)",
        f"  Δr_ss = γv/k_t = {dr_ss.to('nm'):~.3fP}  vs  ℓ_k = {l_k.to('nm'):~.2fP}",
    ]
    plan = [
        f"  dt      = {dt.to('us'):~.4gP}  = {f(dt/tau_B):.3e} τ_B   (트랩이 정함)",
        f"  T_obs   = {T_obs:~.4gP} (박스 횡단 L/v_x)",
        f"  steps   = {f(T_obs/dt):,.0f}   × N={N}",
    ]
    return "trap-drag (육방 격자 + 이동 트랩)", s, lg, groups, ck, inp, der, plan


# ══════════════════════════════════════════════════════════════════════
def chain_bend():
    s = P.load(ROOT / "intake/chain-bend-2d-oscill")
    d = node(s, "particle", "diameter").to("m")
    kT = ds(s, "kT").to("J")
    gam = ds(s, "gamma").to("kg/s")
    k_t = node(s, "external", "stiffness").to("N/m")
    amp = node(s, "external", "amplitude").to("m")
    tau_k, tau_chain = ds(s, "tau_k").to("s"), ds(s, "tau_chain").to("s")
    tau_fast, tau_B, tau_p = (ds(s, "tau_fast").to("s"), ds(s, "tau_B").to("s"),
                              ds(s, "tau_p").to("s"))
    l_k, dmax = ds(s, "l_k").to("m"), ds(s, "delta_max").to("m")
    k_end = ds(s, "kappa_end").to("N/m")
    L = node(s, "geometry", "contour_length").to("m")
    dt = node(s, "numerics", "dt").to("s")
    T_obs = Q(s.raw["numerics"]["production_s"], "s")
    n = int(s.raw["geometry"]["n_beads"])
    kth = node(s, "interactions", 0, "angle_stiffness").to("J")
    Mc = node(s, "interactions", 0, "critical_moment").to("N*m")
    om_lo, om_hi = s.raw["external"]["omega_range"]["value"]
    tau_w_lo = Q(1.0 / om_lo, "s")

    lg = SC.ScaleLedger()
    lg.lengths = {
        "l_k      √(kT/k_t) 트랩 요동": l_k,
        "a        구동 진폭 ★": amp,
        "delta_max M_c 선형 한계 ★": dmax,
        "d        비드 지름 (기준)": d,
        "L        사슬 윤곽길이": L,
    }
    lg.times = {
        "tau_p    m/γ 관성": tau_p,
        "dt       적분 스텝": dt,
        "tau_fast γ/λ_max 최속 굽힘 모드 ★": tau_fast,
        "tau_k    γ/k_t 트랩": tau_k,
        "tau_chain γ/κ_center 집단 굽힘 ★": tau_chain,
        "tau_w    1/ω_min 최저 구동": tau_w_lo,
        "tau_B    d²/D_t 확산 (기준)": tau_B,
        "T_obs    관측창": T_obs,
    }
    lg.energies = {
        "kT       열에너지 (기준)": kT,
        "kappa_th 결합각 강성 ★": kth,
    }
    lg.derived = {"gamma": gam, "d": d, "kT": kT}
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " ★ 이 계는 시간척도가 4자릿수에 걸쳐 있다 — 최속 굽힘 모드"
        "(0.28 µs)가 dt를 정하는데 관심 모드는 τ_chain(1.27 ms)이다.")
    lg.rationale = lg.ref["rationale"]

    groups = {
        "kappa_th/kT            ★ 사슬이 열적으로 뻣뻣": f(kth / kT),
        "a/l_k                  진폭 vs 열요동 = SNR": f(amp / l_k),
        "a/delta_max            ★ 선형 탄성 여유": f(amp / dmax),
        "kappa_end/k_t          사슬 vs 트랩 강성": f(k_end / k_t),
        "De_max = w_max·tau_k   최고 Deborah": om_hi * f(tau_k / Q(1, "s")),
        "tau_fast/tau_chain     ★ 척도 분리 폭": f(tau_fast / tau_chain),
        "n_beads": float(n),
        "St     = tau_p/tau_B": f(tau_p / tau_B),
    }
    ck = [
        C.Check("model", "관성 무시     τ_p/τ_chain", f(tau_p / tau_chain), C.GATE, "<=",
                "τ_dyn = **관심** 최속 척도 = τ_chain (집단 굽힘). "
                "G'(ω)를 재는 대역이 여기다"),
        C.Check("model", "참고: τ_p/τ_fast", f(tau_p / tau_fast), C.GATE, "<=",
                "★ 최속 굽힘 모드는 실제로 **과감쇠가 아니다**(비 0.6). BD는 그 모드를 "
                "과감쇠로 다룬다 — 관측 대역(τ_chain)에서 4570배 떨어져 있어 G'(ω)에 "
                "영향은 없을 것으로 보이나 **확인하지 않았다**", hard=False),
        C.Check("integration", "최속 모드 해상 dt/τ_fast", f(dt / tau_fast), C.GATE, "<=",
                "강성 행렬 최대 고유값 기준. 이걸 놓치면 발산한다"),
        C.Check("integration", "구동 해상     dt/τ_w", f(dt / tau_w_lo), C.GATE, "<=",
                "최저 구동 주기 해상"),
        C.Check("geometry", "선형 탄성     a/δ_max", f(amp / dmax), 1.0, "<=",
                "★ M < M_c. 넘으면 결합이 미끄러져 조화 angle 퍼텐셜이 무효 "
                "(논문 [P2] 결론)"),
        C.Check("statistics", "SNR          a/ℓ_k", f(amp / l_k), 3.0, ">=",
                "진폭이 트랩 안 열요동보다 충분히 커야 위상 추출이 된다", hard=False),
        C.Check("statistics", "관측창       T_obs·ω_min/2π", f(T_obs / (2 * math.pi * tau_w_lo)),
                100.0, ">=", "최저 주파수에서 주기 수", hard=False),
    ]
    inp = [
        R.kv("d", f"{d.to('um'):~.3gP}", 2, "[P1] '2a = 1.47 µm' — 스케치 R=5µm와 3.4배 불일치",
             val_w=22),
        R.kv("kappa_0", "64 mN/m", 2, "[P1] 10 mM MgCl₂ 실험값", val_w=22),
        R.kv("M_c", f"{Mc.to('pN*um'):~.3gP}", 2, "[P1] 고염 플래토", val_w=22),
        R.kv("k_t", f"{k_t.to('pN/um'):~.3gP}", 0, "sketch (논문은 ~40)", val_w=22),
        R.kv("n_beads", str(n), 3, "★제안 [P1] Fig.4 — n=11은 진폭 창이 닫힘", val_w=22),
        R.kv("amplitude", f"{amp.to('nm'):~.3gP}", 3, "★제안 — 창 20nm≪a<429nm", val_w=22),
    ]
    der = [
        f"  EI = {node(s,'interactions',0,'flexural_rigidity'):~.4eP}"
        f"   κ_θ = EI/ℓ = {kth:~.4eP} = {f(kth/kT):.3e} kT",
        f"  κ_end = {k_end.to('pN/um'):~.3fP} (논문 정의, 끝 힘)"
        f"   κ_center = {(2*k_end).to('pN/um'):~.3fP} (구동 트랩이 느낌)",
        f"  ★ 힘 정의를 섞으면 논문값과 정확히 2배 어긋난다 (검증으로 확인)",
        f"  δ_max = M_c L²/(12EI) = {dmax.to('nm'):~.0fP}   ℓ_k = {l_k.to('nm'):~.2fP}",
    ]
    plan = [
        f"  dt      = {dt.to('ns'):~.4gP}  = {f(dt/tau_B):.3e} τ_B   (최속 굽힘 모드가 정함)",
        f"  ω 스윕  = {om_lo:.0f} ~ {om_hi:.0f} rad/s  (De = ωτ_k 0.1~10)",
        f"  T_obs   = {T_obs:~.4gP} (최저 ω에서 100주기)",
        f"  steps   = {f(T_obs/dt):,.0f}   × n={n}   ★ 비싸다 — 아래 VERDICT 참조",
    ]
    return "chain-bend (3점 굽힘 마이크로레올로지)", s, lg, groups, ck, inp, der, plan


# ══════════════════════════════════════════════════════════════════════
BUILDERS = {"abp": abp_rod, "trap-drag": trap_drag, "chain": chain_bend}
want = sys.argv[1:] or list(BUILDERS)
for key in want:
    title, s, lg, groups, ck, inp, der, plan = BUILDERS[key]()
    txt, verdict = R.render(
        title=f"DimensionlessReport — {s.label}   [{title}]",
        ref=lg.ref, ledger=lg, groups=groups, checks=ck,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(txt)
    print()
