"""Phase 1-B `soft-r3-2d-A-sweep` 설계 계산 — 사용자 확정값 반영 (2026-08-03).

확정된 선택 (사용자, 2026-08-03):
    밀도      φ = 0.35  (a_mean = 1.5 d)
    코어      WCA를 d에 추가 (ε = kT)
    꼬리      N을 늘려 r_c를 제대로 확보
    앵커      d = 5 µm  (다른 스케치의 R=5µm 관례와 통일)

아직 미확정: A가 무차원인가 µm³ 차원인가 (§0에서 두 해석을 나란히 계산해 판정 근거를 냄)

용법: $PY scratch/soft_r3_design.py
"""

import math

import numpy as np
import pint

u = pint.UnitRegistry()
Q = u.Quantity
PI = math.pi

kB = Q(1.380649e-23, "J/K")

# ── 확정 앵커 ──────────────────────────────────────────────────────────
d = Q(5.0, "um")  # 사용자 확정 (다른 스케치 R=5µm 관례)
T = Q(300, "K")  # 1-A 확인값 승계
ETA = Q(0.851, "mPa*s").to("Pa*s")  # 물@300K
RHO_P = Q(2000, "kg/m^3")  # 실리카 가정, τ_p 전용
PHI = 0.35  # 사용자 확정
A_LIST = (0.1, 1.0, 10.0, 100.0)  # 스케치
EPS_WCA = 1.0  # WCA 깊이 [kT] — 관례. 코어 위치만 정하는 역할

KT = (kB * T).to("J")
gamma = (3 * PI * ETA * d).to("kg/s")
D_t = (KT / gamma).to("um^2/s")
tau_B = (d**2 / D_t).to("s")
m = (RHO_P * (PI / 6) * d**3).to("kg")
tau_p = (m / gamma).to("s")

a_mean = d * math.sqrt(PI / (4 * PHI))  # 2D: φ = (π/4)(d/a)²
a_star = float(a_mean / d)


def hr(t, ch="="):
    print(f"\n{ch * 78}\n{t}\n{ch * 78}")


# ══════════════════════════════════════════════════════════════════════
hr("0.  A의 차원 해석 — d = 5 µm 에서는 이 선택이 물리를 125배 바꾼다")
print(f"""
  해석 (i)  A 무차원, r을 지름 단위로:   U/kT = A (d/r)³
  해석 (ii) [A] = µm³, r을 µm로:        U/kT = A[µm³] / r[µm]³

  평균간격에서의 결합세기 Γ ≡ U(a_mean)/kT   (구조를 정하는 실제 제어 파라미터)
    a_mean = {a_mean:~.2fP} = {a_star:.3f} d
""")
print(f"    {'A':>7} | {'Γ  해석(i)':>12} | {'Γ  해석(ii)':>13} | 해석(i) 체제")
for A in A_LIST:
    G_i = A / a_star**3
    G_ii = A / a_mean.to("um").magnitude ** 3
    regime = "약상관 유체" if G_i < 1 else ("상관 유체" if G_i < 10 else "강결합(결정 후보)")
    print(f"    {A:7.1f} | {G_i:12.4f} | {G_ii:13.5f} | {regime}")
print(f"""
  ★ 해석 (ii)면 A=100 에서도 Γ = 0.24 — 스윕 전 구간이 kT 아래다.
    그러면 rdf·voronoi에 WCA 배제부피 말고는 아무 구조도 안 나온다.
    스케치가 'final configuration? / rdf / voronoi / structure analysis'를
    목표로 적었으므로, **해석 (i)이 목표와 일관된 유일한 읽기**다.
    → 이하 (i)로 진행. 확인 필요 (observation.yaml A1).
""")

# ══════════════════════════════════════════════════════════════════════
hr("1.  물리계 (SI) — 확정 앵커에서 유도")
print(f"    d     = {d:~P}          (사용자 확정)")
print(f"    T     = {T:~P}          η = {ETA.to('mPa*s'):~.3fP}   kT = {KT:~.4eP}")
print(f"    γ     = 3πηd = {gamma:~.4eP}")
print(f"    D_t   = kT/γ = {D_t:~.4fP}")
print(f"    τ_B   = d²/D_t = {tau_B:~.2fP}   ★ 기준 시간 (이 계의 지배 척도)")
print(f"    τ_p   = m/γ = {tau_p.to('us'):~.3fP}   (모델 검사 전용)")
print(f"    φ     = {PHI}  →  a_mean = {a_mean:~.3fP} = {a_star:.3f} d")


# ══════════════════════════════════════════════════════════════════════
# 퍼텐셜 (무차원: 길이는 d, 에너지는 kT)
#   U*(r*) = U_WCA(r*) + A / r*³      r* = r/d
# ══════════════════════════════════════════════════════════════════════
R_WCA = 2 ** (1 / 6)


def U_star(rs, A):
    w = np.where(rs < R_WCA, 4 * EPS_WCA * (rs**-12.0 - rs**-6.0) + EPS_WCA, 0.0)
    return w + A / rs**3


def U2_star(rs, A):
    """U''(r) [kT/d²] — 국소 강성."""
    w = np.where(rs < R_WCA, 4 * EPS_WCA * (156 * rs**-14.0 - 42 * rs**-8.0), 0.0)
    return w + 12 * A / rs**5


def solve_r_pair(A, u_max=5.0):
    """U*(r) = u_max 가 되는 거리 (희박/쌍 접근 기준). 이분법."""
    lo, hi = 0.5, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if U_star(np.array([mid]), A)[0] > u_max:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


hr("2.  최근접 접근거리 r_min — dt는 여기서 나온다")
LINDEMANN = 0.15  # u_rms/a 가 이보다 크면 케이지(결정) 그림이 무효 → 유체
print(f"""  두 기준을 계산해 **더 작은 쪽**(= 더 뻣뻣한 쪽)을 쓴다. 보수적 선택.
    (a) 쌍 기준   : U(r) = 5 kT 인 거리        — 희박·약결합에서 지배
    (b) 진동 기준 : a_mean − 3 u_rms           — 조밀·강결합에서만 유효
        2D 육방 케이지 강성 k_cage = 3[U''(a) + U'(a)/a],  u_rms = √(kT/k_cage)
        ⚠ u_rms/a > {LINDEMANN} 이면 케이지가 녹은 것이므로 (b)를 쓰지 않는다
          (그대로 쓰면 r_cage가 음수까지 나온다 — 첫 계산에서 실제로 그랬다)
""")
print(f"    {'A':>7} {'Γ':>8} | {'(a) r_pair':>11} | {'u_rms/a':>8} {'(b) r_cage':>11} |"
      f" {'r_min':>8} {'U″(r_min)':>10} {'τ_int/τ_B':>11} {'dt/τ_B':>10}  기준")
design = {}
for A in A_LIST:
    G = A / a_star**3
    r_pair = solve_r_pair(A)
    # 케이지 강성 (a_mean=1.5d 에서 WCA는 0이므로 r⁻³ 성분만 기여)
    Upp_a = float(U2_star(np.array([a_star]), A)[0])
    Up_a = -3 * A / a_star**4  # U'(a) [kT/d]
    k_cage = 3 * (Upp_a + Up_a / a_star)
    u_rms = math.sqrt(1.0 / k_cage) if k_cage > 0 else float("inf")
    ur = u_rms / a_star
    crystalline = ur < LINDEMANN
    r_cage = a_star - 3 * u_rms
    if crystalline and r_cage < r_pair:
        r_min, which, state = r_cage, "진동", "결정"
    else:
        r_min, which, state = r_pair, "쌍", ("결정" if crystalline else "유체")
    Upp = float(U2_star(np.array([r_min]), A)[0])
    tau_int_star = 1.0 / Upp  # τ_int/τ_B = (γ/U'')/(γd²/kT) = 1/U''*
    dt_star = 1e-2 * tau_int_star
    design[A] = dict(G=G, r_min=r_min, dt_star=dt_star, which=which, u_rms=ur, state=state)
    cage_txt = f"{r_cage:10.3f}d" if crystalline else "   (무효)  "
    print(
        f"    {A:7.1f} {G:8.3f} | {r_pair:10.3f}d | {ur:8.4f} {cage_txt} |"
        f" {r_min:7.3f}d {Upp:10.1f} {tau_int_star:11.3e} {dt_star:10.2e}  {which}({state})"
    )
print(f"""
    ★ A=0.1·1 은 **WCA 코어가 dt를 정한다** (r⁻³는 거의 무의미, Γ≤0.3).
      A=100 은 케이지에 갇혀 접근 자체가 제한되어 오히려 dt가 느슨하다.
      → 이 스윕에서 가장 비싼 것은 물리적으로 가장 심심한 A=0.1 이다.
    ★ dt는 A마다 다르게 쓴다 (bd-physics §1.1: 단위계는 고정, dt는 계마다).
    ⚠ 비선형계라 편향 닫힌형태가 없다 → dt 절반 수렴 확인 필요 (mater_plan §20 B)
    ⚠ u_rms/a 는 Lindemann 지표이기도 하다: A=100 에서 {design[100.0]['u_rms']:.3f}
      → 결정 후보. 실제 판정은 실행 후 voronoi/ψ₆ 로 (선언하지 않는다)
""")


# ══════════════════════════════════════════════════════════════════════
hr("3.  컷오프와 N — r⁻³ 꼬리의 정직한 처리")
print("""  2D에서 r⁻³ 꼬리 에너지는 ∫r⁻³·r dr ~ 1/r_c 로 **느리게** 수렴한다.
    r_c 밖 에너지 / 최근접이웃 에너지 = 2π a /(3 r_c)
  → 절대 에너지·압력은 꼬리보정이 필요하고, **구조(rdf·voronoi)는 영향이 작다**
    (균질한 먼 장에서 오는 힘은 대칭으로 상쇄. 단 ξ < r_c 여야 함)
  → 컷오프는 '절대 kT 기준'이 아니라 **a_mean의 배수**로 잡는 것이 맞다.
""")
print(f"    {'r_c':>10} {'U(r_c)/kT @A=100':>18} {'/Γ':>8} {'꼬리/이웃':>10} | 최소이미지 N 하한")
for k in (3, 5, 7, 10):
    rc = k * a_star
    u_rc = 100.0 / rc**3
    G100 = 100.0 / a_star**3
    tail = 2 * PI / (3 * k)
    N_min = (2 * k) ** 2  # r_c < L/2 = (√N/2)a  →  √N > 2k
    print(f"    {k:2d} a_mean {rc:5.2f}d {u_rc:14.4f} {u_rc / G100:8.4f} {tail:10.3f} | N > {N_min}")

N = 400
L = a_mean * math.sqrt(N)
rc_star = 5 * a_star
print(f"""
  권장:  N = {N}   →  L = {L:~.1fP} = {float(L / d):.1f} d = {math.sqrt(N):.0f} a_mean
         r_c = 5 a_mean = {rc_star:.2f} d = {(rc_star * d).to('um'):~.1fP}
         최소이미지 여유 = (L/2)/r_c = {math.sqrt(N) / 2 / 5:.2f}×
         ★ N=400이면 r_c = 5·7·9 a_mean 을 전부 최소이미지 안에서 시험할 수 있다
           → r_c 수렴 확인을 같은 N으로 할 수 있음 (A=100에서 rdf 비교)
         이웃 수 ≈ π r_c²/a_mean² = {PI * 25:.0f} 개/입자
""")


# ══════════════════════════════════════════════════════════════════════
hr("4.  분리 검사 (모델 / 적분 / 기하 / 통계) — 1-A와 같은 분류")
T_OBS_STAR = 100.0  # τ_B 단위 — 통계 검사 T_obs ≥ 10² τ_B 를 만족시키는 값
for A in A_LIST:
    dd = design[A]
    print(f"\n  A = {A}   (Γ = {dd['G']:.3f},  dt* = {dd['dt_star']:.2e})")
    rows = [
        ("model", "관성 무시  τ_p/τ_int",
         float(tau_p / tau_B) / (dd["dt_star"] / 1e-2), 1e-2),
        ("integration", "상호작용 해상  dt/τ_int", 1e-2, 1e-2),
        ("geometry", "컷오프  r_c/(L/2)", rc_star / (math.sqrt(N) * a_star / 2), 1.0),
        ("statistics", "관측창  τ_B/T_obs", 1.0 / T_OBS_STAR, 1e-2),
    ]
    for kind, name, val, lim in rows:
        ok = val <= lim
        mark = "✓" if ok else ("⚠" if val <= 1.5 * lim else "✗")
        print(f"    {mark} [{kind}] {name:26s} = {val:9.3e}  ≤ {lim:7.1e}"
              f"   여유 {lim / val:6.1f}×")
    steps = T_OBS_STAR / dd["dt_star"]
    print(f"      비용: T_obs = {T_OBS_STAR:.0f} τ_B = {(T_OBS_STAR * tau_B).to('s'):~.0fP}(물리)"
          f"  →  {steps:.2e} steps")

print(f"""
  네 검사 모두 통과. 다만 통과가 정확도를 보장하지는 않는다 (bd-physics §4):
  ⚠ 통계 검사는 τ_B 기준으로 T_obs ≥ 10²τ_B 를 본 것이다. 강결합(A=100)에서
    결함 어닐링은 τ_B보다 느릴 수 있어 이 검사만으로는 부족할 수 있다.
    → 평형을 선언하지 말고 사후분석(tools/postmortem.py) EQ 진단으로 측정한다.
  ⚠ 적분 검사의 여유가 정확히 1.0× 인 것은 dt를 한계값으로 잡았기 때문이다.
    0.5% 편향에 해당한다 (선형계 기준). 비선형계이므로 dt 절반 수렴 확인 필요.
""")

total = sum(T_OBS_STAR / design[A]["dt_star"] for A in A_LIST)
print(f"  4개 A 합계 ≈ {total:.2e} steps (독립 런 4개 → 병렬 실행 가능, §원칙 6)")
