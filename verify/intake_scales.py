"""intake/ 잔여 4케이스 — 스케일 표 + 무차원화 초안 (2026-08-03).

CLAUDE.md 규칙 1(차원이 먼저) · 3(모든 값에 출처) 을 따른다.
  READ   = 스케치에 적혀 있는 값 (tier 0)
  ANCHOR = 스케치에 없어서 채운 값 (tier 표시 — 무차원화가 결정되지 않으므로 필수)
  OPEN   = 채울 수 없는 값 (지어내지 않음. 사람 확인 대상)

용법: $PY scratch/intake_scales.py
"""

import math

import pint

u = pint.UnitRegistry()
Q = u.Quantity
PI = math.pi

kB = Q(1.380649e-23, "J/K")

# ══════════════════════════════════════════════════════════════════════════
# 공통 앵커 — 4개 스케치 전부 매질·온도를 적지 않았다.
#   1-A(trap-2d-5um)에서 사용자가 "물, 300 K"를 확인했으므로 그 관례를 잇는다.
#   같은 노트에 연달아 그린 스케치이므로 같은 실험계로 읽는 것이 타당하나,
#   스케치 자체에 근거가 없으므로 tier 1(관례 승계)로 표시한다.
# ══════════════════════════════════════════════════════════════════════════
T = Q(300, "K")  # ANCHOR tier1 — 1-A 확인값 승계
ETA = Q(0.851, "mPa*s").to("Pa*s")  # ANCHOR tier0(핸드북 물@300K) / 매질 가정 tier1
RHO_P = Q(2000, "kg/m^3")  # ANCHOR tier3 — 실리카 가정. τ_p(모델 검사)에만 쓰임
KT = (kB * T).to("J")


def stokes(d):
    """구 기준 물성 유도 (bd-physics §2). d = 지름."""
    d = d.to("m")
    gamma = (3 * PI * ETA * d).to("kg/s")
    D_t = (KT / gamma).to("um^2/s")
    tau_B = (d**2 / D_t).to("s")
    m = (RHO_P * (PI / 6) * d**3).to("kg")
    tau_p = (m / gamma).to("s")
    # 회전 (구에만 성립 — 비구형에 쓰지 말 것)
    D_r = (KT / (PI * ETA * d**3)).to("1/s")
    tau_r = (1 / D_r).to("s")
    return dict(d=d, gamma=gamma, D_t=D_t, tau_B=tau_B, m=m, tau_p=tau_p, D_r=D_r, tau_r=tau_r)


def hr(title, ch="="):
    print(f"\n{ch * 78}\n{title}\n{ch * 78}")


def chk(name, value, limit, kind="≤", note=""):
    """분리 검사 한 줄. 여유(margin) 함께 보고 (bd-physics §4)."""
    ok = value <= limit if kind == "≤" else value < limit
    margin = limit / value if value > 0 else float("inf")
    mark = "✓" if ok else "✗"
    print(f"    {mark} {name:26s} = {value:9.3e}  {kind} {limit:7.1e}   여유 {margin:8.1f}×  {note}")
    return ok


print(f"공통 앵커:  T = {T:~P}   η = {ETA.to('mPa*s'):~.3fP}(물@300K)   kT = {KT:~.4eP}")
print(f"           ρ_p = {RHO_P:~P} (실리카 가정, τ_p 검사 전용)")
print("※ 4개 스케치 모두 매질·온도 미기재 → 1-A에서 확인된 '물 300 K' 관례를 승계함")


# ══════════════════════════════════════════════════════════════════════════
# CASE 1-B ─ soft-r3-2d-A-sweep
#   READ : U_ij/kT = A/r_ij³ ,  A = 0.1, 1, 10, 100 ,  N = 100 ,  Lx = Ly, 2D
#   OPEN : 밀도(φ 또는 L 또는 a_mean) — 없으면 A만으로 물리가 결정되지 않는다
#   OPEN : 배제부피 코어 유무 (r⁻³만 있으면 A≲1에서 입자가 서로 관통)
#   ANCHOR: d — A의 무차원 해석을 고정하기 위해 필요
# ══════════════════════════════════════════════════════════════════════════
hr("CASE 1-B  soft-r3-2d-A-sweep   ⬅ mater_plan §16 Phase 1-B 대상")

d_soft = Q(1.0, "um")  # ANCHOR tier3 (아래 근거)
s = stokes(d_soft)
print(f"""
READ    U_ij/kT = A/r_ij³   A ∈ {{0.1, 1, 10, 100}}   N = 100   Lx = Ly (2D)
        목표: final configuration / rdf / voronoi / structure analysis
ANCHOR  d = {d_soft:~P}  (tier 3)
        ★ A는 무차원인가? U/kT = A/r³ 이 차원적으로 성립하려면
            (i) r을 지름 단위로 읽어 A 무차원  →  U/kT = A (d/r)³
            (ii) r을 µm로 읽어 [A] = µm³      →  U/kT = A[µm³]/r[µm]³
          d = 1 µm 앵커를 쓰면 (i)과 (ii)가 **수치적으로 완전히 일치**한다.
          → 앵커를 1 µm로 잡으면 이 모호성이 해소된다 (그래서 1 µm을 골랐다)
OPEN    밀도(φ / L / a_mean) 미기재 — A 단독으로는 구조가 결정되지 않는다 (아래 증명)
OPEN    배제부피 코어(WCA 등) 유무 미기재
""")

print(f"DERIVED (d = {d_soft:~P}, 물@300K)")
print(f"    γ   = {s['gamma']:~.4eP}      D_t = {s['D_t']:~.4fP}")
print(f"    τ_B = d²/D_t = {s['tau_B']:~.4fP}   ★ 이 계의 지배 시간척도 (1-A와 대조: 거기선 τ_k)")
print(f"    τ_p = m/γ    = {s['tau_p'].to('us'):~.3fP}  (모델 검사 전용)")

# ── 이 계에만 있는 새 길이척도: U(r) = c·kT 가 되는 거리 ────────────────
print(f"\n새 길이척도 — r_u(A,c) : U(r)/kT = c 인 거리 = A^(1/3) c^(-1/3) d")
print(f"    {'A':>6} | {'r(U=kT)':>9} {'r(U=5kT)':>9} {'r(U=0.01kT)':>12}   해석")
for A in (0.1, 1.0, 10.0, 100.0):
    r1 = A ** (1 / 3)
    r5 = (A / 5) ** (1 / 3)
    rc = (A / 0.01) ** (1 / 3)
    tag = "코어 없으면 관통" if r5 < 1.0 else "r⁻³만으로 접촉 차단됨"
    print(f"    {A:6.1f} | {r1:8.3f}d {r5:8.3f}d {rc:11.2f}d   {tag}")
print("    → A ≲ 5 는 r⁻³ 반발만으로 겹침을 막지 못한다 (r(U=5kT) < d).")
print("      배제부피 코어가 있는지 확인 필요 — 없으면 '점입자 모델'로 명시해야 한다.")

# ── 진짜 제어 파라미터: 평균간격에서의 결합세기 ────────────────────────
print(f"\n무차원수 — A 단독이 아니라 밀도와의 조합이 구조를 정한다")
print(f"    Γ ≡ U(a_mean)/kT = A (d/a_mean)³      a_mean = ρ^(-1/2) = L/√N  (2D)")
print(f"    φ = (π/4)(d/a_mean)²                  L = a_mean √N")
print(f"\n    {'a/d':>5} {'φ':>7} | " + " ".join(f"{'Γ(A=' + str(A) + ')':>12}" for A in (0.1, 1, 10, 100)))
for a_over_d in (1.2, 1.5, 2.0, 3.0, 4.3):
    phi = (PI / 4) / a_over_d**2
    row = " ".join(f"{A / a_over_d**3:12.3f}" for A in (0.1, 1, 10, 100))
    print(f"    {a_over_d:5.1f} {phi:7.3f} | {row}")
print("    → 같은 A라도 밀도에 따라 Γ가 10배 넘게 바뀐다. 밀도 없이는 스윕이 정의되지 않는다.")

# ── 컷오프 vs 최소이미지: N=100의 한계 (함정 7 영역) ────────────────────
print(f"\n기하 검사 — r⁻³는 꼬리가 길다. 절단오차 u_c 를 지키면서 r_c < L/2 가 가능한가?")
print(f"    r_c = (A/u_c)^(1/3) d   ,   L/2 = (√N/2) a_mean")
print(f"    두 조건을 동시에 만족하는 결합세기 상한:  Γ_max = N^(3/2) u_c / 8   (A와 무관!)")
N_soft = 100
PHI_MAX = 0.9069  # 2D 육방 최밀충전 — 이보다 촘촘할 수 없다
A_MIN = math.sqrt(PI / (4 * PHI_MAX))  # a_mean/d 하한 = 0.934
for u_c in (0.01, 0.1):
    G_max = N_soft**1.5 * u_c / 8
    print(f"\n    u_c = {u_c:<5} (절단 시 남는 퍼텐셜)   N = {N_soft} →  Γ_max(컷오프) = {G_max:.3f}")
    for A in (0.1, 1.0, 10.0, 100.0):
        rc = (A / u_c) ** (1 / 3)  # in units of d
        a_need = max(2 * rc / math.sqrt(N_soft), A_MIN)  # 컷오프 제약 ∨ 최밀충전 제약
        phi_max = (PI / 4) / a_need**2
        G_at = A / a_need**3
        bind = "컷오프가" if 2 * rc / math.sqrt(N_soft) > A_MIN else "최밀충전이"
        print(
            f"      A={A:6.1f}  r_c = {rc:6.2f}d  → a_mean ≥ {a_need:5.2f}d "
            f"(φ ≤ {phi_max:.3f})  Γ ≤ {G_at:7.3f}   ← {bind} 제약"
        )
    for G_t in (10, 60):
        N_need = (8 * G_t / u_c) ** (2 / 3)
        print(f"      Γ = {G_t:3d} 를 달성하려면 N ≥ {N_need:8.0f}")
print("\n    ★ N = 100(스케치)은 Γ ≲ 1 밖에 못 낸다. '강한 결합에서의 최종 구조'를 보려면")
print("      N을 키우거나(수백~수천), 절단오차를 키우거나, r_c = L/2 + 꼬리보정을 써야 한다.")
print("      → 함정 7(최소이미지 위반, 과거 +1856% 오차)이 정확히 이 지점이다.")

# ── dt: 곡률(=국소 강성)에서 역산 ───────────────────────────────────────
print(f"\n적분 해상 — τ_int 를 국소 곡률에서 정의한다 (트랩의 τ_k = γ/k 와 같은 구조)")
print(f"    U'' (r) = 12 A kT d³/r⁵   →   τ_int(r) = γ/U'' = τ_B (r/d)⁵ / (12A)")
print(f"    dt ≤ 10⁻² τ_int(r_min) ,  r_min = 가장 가까이 접근하는 거리")
print(f"\n    {'A':>6} | {'코어없음: r_min(U=5kT)':>22} {'dt/τ_B':>10} {'dt [ms]':>9} |"
      f" {'코어 d: r_min=d':>16} {'dt/τ_B':>10} {'dt [ms]':>9}")
tau_B_ms = s["tau_B"].to("ms").magnitude
for A in (0.1, 1.0, 10.0, 100.0):
    rmin_nc = (A / 5) ** (1 / 3)
    dt_nc = 1e-2 * rmin_nc**5 / (12 * A)
    dt_c = 1e-2 * 1.0 / (12 * A)
    print(
        f"    {A:6.1f} | {rmin_nc:20.3f}d {dt_nc:10.2e} {dt_nc * tau_B_ms:9.4f} |"
        f" {1.0:15.1f}d {dt_c:10.2e} {dt_c * tau_B_ms:9.4f}"
    )
print("    ★ 코어가 없으면 **작은 A가 수치적으로 더 어렵다** (깊이 관통 → 강성 급증).")
print("      코어가 있으면 반대로 큰 A가 어렵다. 이 분기는 물리 선택에 달려 있다.")
print(f"    참고: 열적 한 스텝 √(2 D_t dt) ≤ 0.01 r  →  dt/τ_B ≤ 5e-5 (r/d)²")

# ── 비용 ────────────────────────────────────────────────────────────────
print(f"\n비용 (T_obs = 100 τ_B = {100 * s['tau_B'].to('s').magnitude:.0f} s 가정)")
for A, dt_star in ((1.0, 1e-2 * 1.0 / 12), (100.0, 1e-2 / 1200)):
    steps = 100 / dt_star
    print(f"    A = {A:5.1f} (코어 d)  dt* = {dt_star:.2e}  →  {steps:.2e} steps")


# ══════════════════════════════════════════════════════════════════════════
# CASE ─ trap-drag-2d-hex300
#   READ : N ≈ 300, R = 5 µm, k_t = 10 pN/µm, v_x = 0.5 µm/s, 육방 초기배치
#          U_trap = ½ k_t (Δr)², Δr = |r_trap − r_i|, r_trap(t) = r_trap(0) + v t
#   OPEN : 페어 퍼텐셜 미기재 — '육방 평형'을 만들 수 있는 힘이 스케치에 없다
#   OPEN : 격자 간격 / 박스 / 관측 목표
# ══════════════════════════════════════════════════════════════════════════
hr("CASE  trap-drag-2d-hex300")

d_dr = Q(5.0, "um")  # READ 'R = 5µm' + 1-A에서 확인된 '지름' 관례
k_t = Q(10, "pN/um").to("N/m")
v_x = Q(0.5, "um/s").to("m/s")
t = stokes(d_dr)
tau_k = (t["gamma"] / k_t).to("s")
l_k = ((KT / k_t) ** 0.5).to("nm")
dr_ss = (t["gamma"] * v_x / k_t).to("nm")  # 등속 끌기 시 정상상태 지연
snr = float((dr_ss / l_k).to("dimensionless"))
Pe_drag = float((v_x * d_dr / t["D_t"]).to("dimensionless"))
tau_v = (d_dr / v_x).to("s")

print(f"""
READ    N ≈ 300   R = 5 µm   k_t = 10 pN/µm   v_x = 0.5 µm/s   육방 초기배치
        U_trap = ½ k_t (Δr)² ,  r_trap(t) = r_trap(0) + v_x t   (트랩 1개만 이동)
ANCHOR  d = 5 µm — 'R=5µm'을 지름으로 읽음 (1-A에서 사용자 확인된 관례, tier 1)
OPEN    ★ 페어 퍼텐셜이 없다. 육방 격자를 유지할 반발력이 스케치에 없음
        (soft-r3의 A/r³와 같은 계인가? 아니면 트랩 배열로 격자를 고정하는가?)
OPEN    격자 간격 a / 박스 L / 측정 목표(끌림힘? 격자 변형? 결함 생성?)
""")

print("DERIVED / 스케일 원장 (작은 것부터)")
print(f"    τ_p = {t['tau_p'].to('us'):~.3fP}   (모델 검사)")
print(f"    dt  ≤ 10⁻² τ_k = {(0.01 * tau_k).to('us'):~.2fP}")
print(f"    τ_k = γ/k_t = {tau_k.to('ms'):~.3fP}   ★ 가장 빠른 물리 시간 → dt를 정한다")
print(f"    τ_B = {t['tau_B']:~.1fP}   (트랩이 붙잡아 실현되지 않음 — 1-A와 동일)")
print(f"    τ_v = d/v_x = {tau_v:~.1fP}   ★ 끌기 시간 → T_obs를 정한다")
print(f"    길이:  ℓ_k = √(kT/k) = {l_k:~.2fP}  <  d = {d_dr:~P}")
print(f"           Δr_ss = γv/k_t = {dr_ss:~.3fP}   ← 끌림이 만드는 결정론적 지연")

print("\n무차원수")
print(f"    k*  = k d²/kT      = {float((k_t * d_dr**2 / KT).to('')):9.3e}   트랩 vs 열요동 (매우 뻣뻣)")
print(f"    Pe  = v d/D_t      = {Pe_drag:9.3f}     이류 vs 확산")
print(f"    τ_k/τ_v            = {float((tau_k / tau_v).to('')):9.3e}   트랩이완 vs 끌기 (극단 분리)")
print(f"    ★ SNR = Δr_ss/ℓ_k  = {snr:9.4f}     끌림 신호 vs 트랩 안 열요동")
print(f"      = γv/√(k kT)  →  SNR ∝ v/√k")

v_snr1 = (((k_t * KT) ** 0.5) / t["gamma"]).to("um/s")
print(f"\n    ★★ 이 조건에서 결정론적 지연(={dr_ss:~.2fP})이 열요동(={l_k:~.2fP})의 1/10 이다.")
print(f"       SNR = 1 이 되는 속도:  v = √(k kT)/γ = {v_snr1:~.2fP}   (스케치 값의 10배)")
n_ind = (1 / (0.01 * snr)) ** 2
T_need = (2 * tau_k * n_ind).to("s")
print(f"       Δr_ss 를 1% 정확도로 재려면 독립표본 ≈ {n_ind:.2e} 개")
print(f"         → T_obs ≈ 2 τ_k × 표본 = {T_need:~.0fP} = {T_need.to('hour'):~.2fP} (물리시간)")
print(f"         → dt = 10⁻²τ_k 이면 {float((T_need / (0.01 * tau_k)).to('')):.2e} steps")
print("       (끌리는 입자 1개 기준. 격자 300개의 응답을 함께 쓰면 줄어든다)")

print("\n비용 — 박스를 한 번 횡단하는 데 필요한 시간 (격자간격 a 가정별)")
for a_over_d in (1.5, 2.0, 3.0):
    a = a_over_d * d_dr
    L = (a * math.sqrt(300 * math.sqrt(3) / 2)).to("um")  # 육방: 입자당 면적 (√3/2)a²
    t_cross = (L / v_x).to("s")
    steps = float((t_cross / (0.01 * tau_k)).to(""))
    print(
        f"    a = {a_over_d:.1f}d = {a:~.1fP}  →  L = {L:~.0fP}   횡단 {t_cross:~.0fP}"
        f"   =  {steps:.2e} steps   (ℓ_k/L = {float((l_k / L).to('')):.1e})"
    )


# ══════════════════════════════════════════════════════════════════════════
# CASE ─ chain-bend-2d-oscill
#   READ : 비드 사슬, optical tweezers로 트랩, k_t = 10 pN/µm, R = 5 µm,
#          y = a sin(ωt) (y방향 진동), 목표 G' & G'' vs ω, "Eric Furst 논문 참고"
#   OPEN : ★ U_ij (사슬 결합 퍼텐셜)이 스케치에서 **빈칸**
#   OPEN : a, ω(범위), 비드 개수, 어느 비드가 트랩되는가
# ══════════════════════════════════════════════════════════════════════════
hr("CASE  chain-bend-2d-oscill")

c = stokes(Q(5.0, "um"))
tau_k_c = (c["gamma"] / k_t).to("s")
l_k_c = ((KT / k_t) ** 0.5).to("nm")
omega_1 = (1 / tau_k_c).to("1/s")

print(f"""
READ    비드 사슬(원 ~6개 + "…"), 일부에 × 표시 = "trapped by optical tweezers"
        U = ½ k_t (Δr)²   k_t = 10 pN/µm   R = 5 µm
        y = a sin(ωt)  "oscillation in y-dir"      목표: G' & G'' = ?  (vs ω 그래프)
        "particles are connected by U_ij = ～～"   ← **빈칸으로 남겨져 있음**
        "Eric Furst, 논문 참고"                    ← 문헌 지정 (KB 시드 후보)
OPEN    ★ U_ij 미기재 → 사슬의 탄성·이완 스펙트럼을 계산할 수 없다 (이 케이스의 핵심 물리)
OPEN    진폭 a · 주파수 ω(범위) · 비드 개수 N · 트랩되는 비드 (×는 1·4·끝 비드로 보임)
OPEN    굽힘 강성(폴더명 'chain-bend')은 스케치에 없음 — 각도 퍼텐셜 미기재
""")

print("계산 가능한 부분 (트랩 + 구동만)")
print(f"    τ_k = γ/k_t = {tau_k_c.to('ms'):~.3fP}     ℓ_k = {l_k_c:~.2fP}")
print(f"    새 시간척도 τ_ω = 1/ω  →  무차원수 De = ω τ_k (Deborah)")
f_1 = (omega_1 / (2 * PI)).to("Hz").magnitude
print(f"    ★ De = 1 인 주파수:  ω = 1/τ_k = {omega_1.magnitude:.1f} rad/s  (f = {f_1:.1f} Hz)")
print(f"      G'/G'' 교차는 이 근방에서 나온다 → ω 스윕은 {omega_1.magnitude / 100:.1f} ~ "
      f"{omega_1.magnitude * 10:.0f} rad/s 정도를 덮어야 한다")
print(f"\n    진폭 창 (두 무차원수가 양쪽에서 조인다)")
print(f"      a/ℓ_k ≫ 1  (열요동 위로 신호를 내려면)   ℓ_k = {l_k_c:~.2fP}")
print(f"      a/d   ≪ 1  (선형응답을 유지하려면)        d   = 5 µm")
for a_nm in (20, 100, 200, 500):
    a_q = Q(a_nm, "nm")
    print(
        f"      a = {a_q:~4.0fP}  →  a/ℓ_k = {float((a_q / l_k_c).to('')):6.2f} (SNR)"
        f"   a/d = {float((a_q / Q(5.0, 'um')).to('')):7.4f}"
    )
print("      → a ≈ 100~500 nm 이 두 조건을 동시에 만족하는 창")

print("\n★ 물리 해석에서 짚어야 할 점 (규칙 6: 검증하고 말한다 — 이건 모델의 성질)")
print("    BD + 뉴턴 용매(물)에서는 **용매의 G' = 0, G'' = ηω** 이다 (탄성 없음).")
print("    따라서 여기서 나오는 G'·G''는 '매질'이 아니라 **사슬(결합)의 탄성**이다.")
print("    매질의 점탄성을 재려는 의도라면 BD로는 안 되고 medium.* 모듈이 필요하다")
print("    (1-A observation.yaml의 followup '비뉴턴 매질' 요청과 같은 항목).")

print("\n비용 (dt = 10⁻²τ_k, 최저 ω에서 100주기 관측)")
for de in (0.1, 1.0):
    omega = de / tau_k_c
    period = (2 * PI / omega).to("s")
    T_obs = 100 * period
    print(
        f"    De = {de:4.1f}  ω = {omega.to('1/s'):~.1fP}  주기 {period.to('ms'):~.1fP}"
        f"  T_obs = {T_obs:~.2fP}  →  {float((T_obs / (0.01 * tau_k_c)).to('')):.2e} steps"
    )


# ══════════════════════════════════════════════════════════════════════════
# CASE ─ abp-rod-2d-run-flip
#   READ : ellipsoid "active" particle, run-and-flip motion
#          τ_R : rotation time 0.5 s ,  v : speed  v ≤ 5 µm/s
#          measure MSD, MSAD          그림에 치수 표시 "2R"·"R" (숫자 없음)
#   OPEN : R 값 · 종횡비 · N · 박스 · (아래) 뒤집힘 빈도
# ══════════════════════════════════════════════════════════════════════════
hr("CASE  abp-rod-2d-run-flip")

tau_R = Q(0.5, "s")
v_a = Q(5.0, "um/s")

print(f"""
READ    "ellipsoid, 'active' particle    run-and-flip motion"
        τ_R : rotation time = 0.5 s        v : speed  v ≤ 5 µm/s
        "measure MSD. MSAD."
        그림: 타원 + 치수선 "2R"(판독 다소 불확실) + 아래 화살표와 "R"
OPEN    ★ R에 숫자가 없다 (다른 3개 스케치는 모두 'R = 5 µm'라고 적혀 있음)
OPEN    종횡비 p = a/b · 입자 수 N · 박스 L
""")

# ── τ_R을 회전확산으로 읽으면 크기가 역산된다 (구 근사) ────────────────
d_impl = ((KT * tau_R / (PI * ETA)) ** (1 / 3)).to("um")
print("★ τ_R = 0.5 s 를 '회전확산 시간'으로 읽으면 크기가 역산된다 (구 근사)")
print(f"    τ_r = 1/D_r = πηd³/kT   →   d = (kT τ_R/πη)^(1/3) = {d_impl:~.3fP}")
print("    → 물@300K에서 τ_R=0.5s 는 지름 ≈0.92 µm 구와 정확히 일치한다.")
print("      (타원체는 회전마찰이 더 커서 등가지름은 이보다 약간 작아진다)")
print("    ⇒ 'R = 5 µm' 관례를 여기에 그대로 적용하면 τ_R이 스케치 값과 안 맞는다 (아래 표)")

print(f"\n크기 후보별 무차원수 (v = {v_a:~P}, ℓ_p = v·τ_R = {(v_a * tau_R).to('um'):~.2fP} 로 고정됨)")
print(f"    {'d':>8} {'D_t':>10} {'τ_B':>9} {'τ_v=d/v':>9} {'Stokes τ_r':>11} "
      f"{'Pe=vd/D_t':>10} {'ℓ_p/d':>7} {'D_r*':>6}")
for d_um in (0.5, d_impl.magnitude, 1.0, 5.0, 10.0):
    a = stokes(Q(d_um, "um"))
    Pe = float((v_a * a["d"] / a["D_t"]).to(""))
    l_p = (v_a * tau_R).to("um")
    print(
        f"    {d_um:6.2f}µm {a['D_t'].magnitude:8.4f}µm²/s {a['tau_B'].to('s').magnitude:8.2f}s"
        f" {(a['d'] / v_a).to('s').magnitude:8.3f}s {a['tau_r'].to('s').magnitude:10.3f}s"
        f" {Pe:10.2f} {float((l_p / a['d']).to('')):7.2f} {float((a['D_r'] * a['tau_B']).to('')):6.2f}"
    )
print("    → d를 5 µm로 잡으면 ℓ_p/d = 0.5 : 몸길이 절반도 못 가고 방향이 바뀐다(활성 거의 안 보임)")
print("      d ≈ 1 µm면 ℓ_p/d ≈ 2.5, Pe ≈ 9 — MSD에 활성 구간이 뚜렷이 나온다")
print("      ★ 스케치의 τ_R과 자기일관적인 것은 d ≈ 1 µm 쪽이다")

print("\n★ 관측량과 파라미터의 대응 — 여기서 빠진 파라미터가 하나 더 있다")
print("    타원체의 '몸 방향' n 과 '추진 방향' ±n 은 별개다.")
print("      · MSAD  ← 몸의 회전확산 (τ_R)                       ... 스케치에 있음")
print("      · MSD   ← 추진속도 v + **뒤집힘 빈도 τ_flip**        ... τ_flip이 없음")
print("    두 해석 중 어느 쪽이든 파라미터 하나가 빈다:")
print("      (a) τ_R = 회전확산 시간  → τ_flip 미기재  [크기 역산과 일관 → 이쪽으로 기울어짐]")
print("      (b) τ_R = 뒤집힘 간격    → 회전확산 미기재. 그러면 MSAD가 자명해진다")
print("    run-and-flip의 속도 상관: ⟨v(0)·v(t)⟩ = v² exp(−2t/τ_flip)  (포아송 ±반전)")
print("      → 유효 지속시간 = τ_flip/2 이고 이것이 MSD의 크로스오버를 정한다")
print("    ※ mater_plan §20 질문 10의 모호성을 이 형태로 좁혔다 (해소는 사람 확인)")

print("\n하드 제약 (bd-hoomd 실측): BD에는 HI가 없어 **병진 마찰이 등방**이다")
print("    → MSAD ✓ / 장시간 MSD ✓ / **단시간 MSD의 이방성 ✗**")
print("    스케치는 'MSD, MSAD'라고만 적었고 이방성을 명시하지 않았다")
print("    → mater_plan §20 옵션 A(등방 평균 마찰 + 정확한 γ_r)로 진행 가능. 확인 필요")

hr("요약 — 무차원화가 '완결'된 케이스는 없다. 막힌 지점", "─")
print("""
  soft-r3        밀도(φ/L) 없음 + 코어 유무 불명 → A 스윕이 정의되지 않음
                 게다가 N=100은 Γ ≲ 1 밖에 못 냄 (최소이미지 제약)
  trap-drag      페어 퍼텐셜 없음(육방 격자를 만들 힘이 없음) + 격자간격 없음
                 v=0.5µm/s에서 신호가 열요동의 1/10 (측정 설계 문제)
  chain-oscill   U_ij가 빈칸 + a, ω 없음 → 사슬 이완 스펙트럼 계산 불가
  abp-rod        R 값 없음 + τ_flip(또는 회전확산) 중 하나가 없음
""")
