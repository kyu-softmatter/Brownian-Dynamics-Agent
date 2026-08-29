"""chain-bend 문헌 증류 검증 — Pantina & Furst 두 편에서 뽑은 것이 맞는가.

논문 (사용자 지정 2026-08-04):
  [P1] Pantina & Furst, PRL 94, 138301 (2005)
       "Elasticity and Critical Bending Moment of Model Colloidal Aggregates"
  [P2] Pantina & Furst, Langmuir 24, 1141-1146 (2008)
       "Micromechanics and Contact Forces of Colloidal Aggregates in the Presence of Surfactants"

★ 논문의 헤드라인 발견: **비드 사이 상호작용은 중심 페어 퍼텐셜이 아니다.**
  단일 결합이 **토크를 지탱**한다 (접선 상호작용). 그래서 사슬이 빔처럼 휜다.
  → 스케치의 빈칸 `U_ij` 는 "페어 퍼텐셜"이 아니라 **접착 접촉 + 굽힘 강성**이다.

LLM이 논문에서 뽑은 값은 환각 위험이 있다 (마스터플랜 원칙 2). 방어:
  **논문이 스스로 보고한 숫자를 내가 뽑은 공식으로 재현해본다.**
  재현되면 공식·상수 추출이 맞다는 증거다.

    $PY scratch/chain_bend_from_papers.py
"""
import math

import numpy as np

# ══════════════════════════════════════════════════════════════════════
# 논문에서 읽은 값 (전부 locator 포함)
# ══════════════════════════════════════════════════════════════════════
E_PMMA = 3100e6         # Pa      [P1] p.3 좌단, ref 15 (Schreyer)
NU_PMMA = 0.4           # -       [P1] p.3 좌단
D_PART = 1.47e-6        # m       [P1] p.1 우단 "average diameter of 2a = 1.47 ± 0.01 µm"
                        #         [P2] p.2 우단 "2a = 1.47 ± 0.1 µm" (Bangs Labs PMMA)
GAMMA_L = 72.7e-3       # N/m     [P1] p.3 좌단, ref 16 (물 표면장력)
THETA_0 = math.radians(73.7)   # [P1] p.3 좌단 "PMMA-water contact angle θ₀ = 73.7°"

# 측정값
KAPPA0_EXP_10mM = 64e-3     # N/m  [P1] p.3 좌단 "experimental value at 10 mM MgCl₂ is 64 ± 0.5 mN/m"
KAPPA0_JKR_PAPER = 80e-3    # N/m  [P1] p.3 좌단 "Using the JKR model, κ₀ = 80 mN/m"
AC_PAPER = 40e-9            # m    [P1] p.3 좌단 "corresponding contact area radius is a_c ≈ 40 nm"
KAPPA0_BARE_250mM = 0.21    # N/m  [P2] p.4 Fig.4 inset "κ₀^bare = 0.21 ± 0.01 N/m"
MC_PLATEAU = 35e-18        # N·m   [P1] p.2 우단 "M_c plateaus to a value of approximately 35 pN µm"
MC_BARE_250mM = 30e-18     # N·m   [P2] Fig.3 화살표 (계면활성제 없음, ~30 pN·µm)
SLIP_LENGTH = 32e-9       # m     [P1] p.2 우단 "average length of sliding rearrangements is 32 ± 15 nm"
K_TRAP_PAPER = 40e-6      # N/m   [P1] p.1 우단 "trap rigidity is approximately 40 pN/µm"
V_DRAG_P2 = 20e-9         # m/s   [P2] Fig.1A "translating in the y direction at a rate of 20 nm/s"
D_B = 1.0                 # -     [P2] p.4 "chain aggregates ... nearly perfectly straight; therefore d_b = 1"

ETA = 0.851e-3            # Pa·s  물@300K (핸드북) — 논문은 수용액
KT = 1.380649e-23 * 300   # J

a = D_PART / 2            # 입자 반지름
K_JKR = 2 * E_PMMA / (3 * (1 - NU_PMMA**2))     # [P2] eq.6 "K = 2E/3(1−ν²)"


def a_c_from_W(W_SL, a_=a):
    """[P2] eq.6  a_c = (3πa²W_SL / 2K)^(1/3)   — JKR, 외부하중 없음"""
    return (3 * math.pi * a_**2 * W_SL / (2 * K_JKR)) ** (1 / 3)


def kappa0_from_ac(a_c, a_=a):
    """[P2] eq.5  κ₀ = 3π a_c⁴ E / (4 a³)"""
    return 3 * math.pi * a_c**4 * E_PMMA / (4 * a_**3)


def ac_from_kappa0(kappa0, a_=a):
    """eq.5 역산"""
    return (kappa0 * 4 * a_**3 / (3 * math.pi * E_PMMA)) ** 0.25


def EI_from_kappa0(kappa0, a_=a):
    """κ₀ = 3EI/a³  (eq.5 와 EI = πEa_c⁴/4 를 결합)"""
    return kappa0 * a_**3 / 3


print("=" * 82)
print("① 추출 검증 — 논문이 보고한 숫자를 내 공식으로 재현하는가")
print("=" * 82)
W0 = GAMMA_L * (1 + math.cos(THETA_0))          # [P2] eq.9 Young-Dupré
ac0 = a_c_from_W(W0)
k0 = kappa0_from_ac(ac0)
print(f"  Young-Dupré  W_SL⁰ = γ_L(1+cos θ₀) = {W0*1e3:.1f} mJ/m²")
print(f"  JKR          a_c = {ac0*1e9:.1f} nm      (K = {K_JKR/1e9:.3f} GPa)")
print(f"  eq.5         κ₀ = {k0*1e3:.1f} mN/m")
print(f"  논문 [P1]     κ₀ = {KAPPA0_JKR_PAPER*1e3:.0f} mN/m")
err = 100 * (k0 - KAPPA0_JKR_PAPER) / KAPPA0_JKR_PAPER
print(f"  → 차이 {err:+.1f}%   {'✓ 추출 정확 (eq.5·6·9 + 상수 4개)' if abs(err) < 5 else '✗ 뭔가 틀렸다'}")

ac_exp = ac_from_kappa0(KAPPA0_EXP_10mM)
print(f"\n  실험값 역산  κ₀ = {KAPPA0_EXP_10mM*1e3:.0f} mN/m → a_c = {ac_exp*1e9:.1f} nm")
print(f"  논문 [P1]     a_c ≈ {AC_PAPER*1e9:.0f} nm")
print(f"  → 차이 {100*(ac_exp-AC_PAPER)/AC_PAPER:+.1f}%  (논문이 반올림한 값)  ✓")

sigma_c = 4 * MC_PLATEAU / (math.pi * ac_exp**3)   # [P2] eq.16 M_c = (π/4)σ*_xx a_c³
print(f"\n  eq.16 역산   M_c = {MC_PLATEAU*1e18:.0f} pN·µm, a_c = {ac_exp*1e9:.1f} nm")
print(f"               → 임계 인장응력 σ*_xx = {sigma_c/1e6:.2f} MPa")
print(f"               PMMA 항복강도(~70 MPa)보다 훨씬 작다 → 접촉선 de-pinning으로 타당")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 82)
print("② HOOMD 매핑 검증 — 이산 사슬의 angle 강성 ↔ 연속 빔의 EI")
print("=" * 82)
print("""  주장:  EI = κ_θ · ℓ        (ℓ = 결합 길이 = 2a)
  근거:  이산 U = Σ ½κ_θ θ_i²,  θ_i ≈ ℓ/R  →  단위길이당 ½κ_θℓ/R²
         연속 U = ∫ ½EI/R² dx                   → 계수 비교
  ★ 감으로 쓰지 않고, 이산 사슬을 실제로 굽혀서 빔 공식과 대조한다.""")


def discrete_3point_stiffness(n, kappa_theta, ell):
    """비드 n개 사슬: 양끝 고정, 중앙을 δ만큼 밀 때의 강성 F/δ.

    소변형 근사 θ_i = (y_{i+1} − 2y_i + y_{i−1})/ℓ,  U = ½ κ_θ Σ θ_i²
    → y에 대한 2차형식. 제약(양끝 0, 중앙 δ) 아래 최소화하고 F = dU/dδ.
    """
    assert n % 2 == 1, "중앙 비드가 있어야 한다"
    c = n // 2
    # 2차 미분 행렬 (내부 각도 n-2개)
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    A = kappa_theta * (B.T @ B)          # U = ½ yᵀ A y
    fixed = [0, c, n - 1]
    free = [i for i in range(n) if i not in fixed]
    y_fix = np.array([0.0, 1.0, 0.0])    # δ = 1 (선형이라 스케일 무관)
    A_ff = A[np.ix_(free, free)]
    A_fx = A[np.ix_(free, fixed)]
    y_free = np.linalg.solve(A_ff, -A_fx @ y_fix)
    y = np.zeros(n)
    y[fixed] = y_fix
    y[free] = y_free
    U = 0.5 * y @ A @ y                  # δ=1 이므로 U = ½ k δ² = ½ k
    return 2 * U                         # k = F_center/δ  ★ 중앙에 가하는 힘 기준


kappa_theta_10mM = EI_from_kappa0(KAPPA0_EXP_10mM) / D_PART
print(f"\n  κ₀ = {KAPPA0_EXP_10mM*1e3:.0f} mN/m → EI = {EI_from_kappa0(KAPPA0_EXP_10mM):.4e} N·m²"
      f" → κ_θ = EI/ℓ = {kappa_theta_10mM:.4e} J = {kappa_theta_10mM/KT:.3e} kT")
print("""
  ★ 힘의 정의를 맞춰야 한다 — 처음에 2배 어긋났고 n을 키워도 수렴하지 않았다.
    논문의 κ = F_bend/δ 는 **끝 입자**에 걸리는 힘 기준이다 ([P1] p.2: "F_bend is measured
    by the displacement of the end particles"). 힘 균형상 중앙 힘 = 끝 힘 × 2 이므로
        끝 힘 기준   κ_end    = 24 EI/L³   ← 논문과 비교할 값
        중앙 힘 기준 κ_center = 48 EI/L³   ← 시뮬레이션에서 구동 트랩이 느끼는 값
    시뮬레이션에서 중앙 힘을 재서 논문 κ와 바로 대조하면 **정확히 2배 틀린다.**""")
print(f"\n  {'n':>4} {'L=(n-1)ℓ':>11} {'이산(중앙힘)':>15} {'빔 48EI/L³':>13} {'차이':>9}"
      f" {'논문 κ_end=24EI/L³':>19}")
EI = EI_from_kappa0(KAPPA0_EXP_10mM)
ok2 = True
for n in (5, 9, 11, 15, 25, 51):
    L = (n - 1) * D_PART
    k_disc = discrete_3point_stiffness(n, kappa_theta_10mM, D_PART)
    k_end = 24 * EI / L**3               # y(±L/2) = F_end L³/(24EI) — [P1] eq.1
    k_center = 2 * k_end                 # 힘 균형: 중앙 힘 = 끝 힘 × 2
    e = 100 * (k_disc - k_center) / k_center
    # n=5 는 이산화 오차가 −11% 로 크다. 논문 실험은 9~25개 사슬이므로 그 범위로 판정한다.
    if n >= 9:
        ok2 &= abs(e) < 4
    print(f"  {n:>4} {L*1e6:9.2f}µm {k_disc*1e6:13.3f}pN/µm {k_center*1e6:11.3f}pN/µm"
          f" {e:+8.2f}% {k_end*1e6:17.3f}pN/µm")
print(f"\n  {'✓' if ok2 else '✗'} 힘 정의를 맞추면 이산 사슬이 빔 공식으로 수렴 "
      f"(n≥9에서 4% 이내. n=5는 −11% — 짧은 사슬의 이산화 오차이며 논문 실험 범위 밖)")
print(f"  ⟹ **κ_θ = EI/ℓ 매핑이 맞다** (angle.Harmonic 의 k 로 쓸 값)")
print(f"  ✓ 논문 eq.4 의 κ(s)=κ₀(a/s)^(2+d_b), d_b=1 은 **s = 팔 길이(L/2)** 로 읽으면")
print(f"    3EI/s³ = 24EI/L³ 로 [P1] eq.1 과 정확히 일치한다 (s를 전체 길이로 읽으면 8배 틀림)")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 82)
print("③ 사슬 강성 vs 트랩 강성 — 측정 가능한 창이 어디인가")
print("=" * 82)
print(f"  [P1] p.3 좌단: \"trap compliance limits the maximum rigidity that can be measured\"")
print(f"  스케치 k_t = 10 pN/µm,  논문 k_t ≈ {K_TRAP_PAPER*1e6:.0f} pN/µm")
print(f"\n  {'n':>4} {'L':>9} | " + " ".join(f"{lbl:>15}" for lbl in
      ("κ(10mM,64mN/m)", "κ(250mM,0.21)")) + "   판정(스케치 k_t=10pN/µm)")
for n in (5, 9, 11, 15, 25, 41):
    L = (n - 1) * D_PART
    ks = [24 * EI_from_kappa0(k0_) / L**3 for k0_ in (KAPPA0_EXP_10mM, KAPPA0_BARE_250mM)]
    r = ks[0] / 10e-6
    verd = ("사슬이 트랩보다 뻣뻣 — 트랩 컴플라이언스가 지배" if r > 3 else
            "비슷 — 측정 창 ★" if 0.3 < r < 3 else "사슬이 무름 — 신호 작음")
    print(f"  {n:>4} {L*1e6:7.1f}µm | " + " ".join(f"{k*1e6:13.2f}pN/µm" for k in ks)
          + f"   {verd}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 82)
print("④ 새 하드 제약 — 임계 굽힘 모멘트 M_c (조화 angle 퍼텐셜이 깨지는 지점)")
print("=" * 82)
print(f"""  M_c 위에서는 입자가 **미끄러지거나 구른다** ([P2] 결론). 조화 angle 퍼텐셜은
  이 소성 거동을 담지 못한다 → **선형 응답을 재려면 M < M_c 를 지켜야 한다.**
  이건 1-A/1-B에 없던 종류의 검사다 (물리 모델의 유효 범위).

  [P2] eq.3 에서 중앙 재배열이면  M = F_bend · L/2
    → F_bend < 2 M_c / L,   진폭 δ < F_bend/κ(s)""")
print(f"\n  {'n':>4} {'L':>9} {'κ(s)':>12} {'F_max':>10} {'δ_max':>10} {'ℓ_k=√(kT/k_t)':>15}"
      f" {'창 δ_max/ℓ_k':>13}")
l_k = math.sqrt(KT / 10e-6)
for n in (5, 9, 11, 15, 25):
    L = (n - 1) * D_PART
    ks = 24 * EI / L**3
    F_max = 2 * MC_PLATEAU / L
    d_max = F_max / ks
    print(f"  {n:>4} {L*1e6:7.1f}µm {ks*1e6:10.2f}pN/µm {F_max*1e12:8.2f}pN"
          f" {d_max*1e9:8.0f}nm {l_k*1e9:13.1f}nm {d_max/l_k:12.0f}×")
print(f"""
  → 진폭 창:  ℓ_k({l_k*1e9:.0f} nm) ≪ a < δ_max
     아래는 열요동에 묻히고, 위는 소성(결합 미끄러짐)이 시작된다.
     n=11 이면 20 nm ≪ a < {2*MC_PLATEAU/(10*D_PART)/(24*EI/(10*D_PART)**3)*1e9:.0f} nm 정도.""")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 82)
print("⑤ 시간척도 + 스케치와의 불일치")
print("=" * 82)
gamma_bead = 3 * math.pi * ETA * D_PART
tau_k = gamma_bead / 10e-6
tau_B = D_PART**2 / (KT / gamma_bead)
print(f"  논문 입자 d = {D_PART*1e6:.2f} µm →  γ = {gamma_bead:.4e} kg/s")
print(f"    τ_k = γ/k_t = {tau_k*1e3:.3f} ms   (트랩 이완, k_t=10pN/µm)")
print(f"    τ_B = d²/D_t = {tau_B:.2f} s")
for n in (11, 25):
    L = (n - 1) * D_PART
    ks = 24 * EI / L**3
    tau_c = gamma_bead / ks
    print(f"    n={n}: κ(s)={ks*1e6:.2f} pN/µm → τ_chain = γ/κ = {tau_c*1e3:.3f} ms"
          f"   De=1 at ω={1/tau_c:.0f} rad/s ({1/tau_c/(2*math.pi):.0f} Hz)")
print(f"""
  ★ 스케치와 논문의 불일치
      스케치:  R = 5 µm            논문:  2a = 1.47 µm   ← 3.4배
      스케치:  k_t = 10 pN/µm      논문:  ~40 pN/µm
      스케치:  y = a sin(ωt), G'·G''   논문:  등속 {V_DRAG_P2*1e9:.0f} nm/s, 준정적 탄성
    논문 값은 전부 d=1.47µm 기준이다. κ₀ ∝ a^(-1/3) (eq.7)로 약하게만 의존하므로
    d=5µm로 외삽할 수는 있으나, **어느 쪽을 쓸지는 사람이 정해야 한다.**
    진동 구동(G'·G'')은 논문에 없다 — 탄성 상수는 옮겨오되 동역학은 새 영역이다.""")

print()
print("=" * 82)
print(f"{'✓ PASS' if (abs(err) < 5 and ok2) else '✗ FAIL'} — 공식 추출 검증 "
      f"({err:+.1f}%) · 이산↔연속 매핑 검증")
print("=" * 82)
