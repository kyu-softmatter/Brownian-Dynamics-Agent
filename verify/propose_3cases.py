"""막힌 3케이스의 제안값 산출 — system.yaml 에 넣을 숫자를 한 곳에서 계산한다.

사용자 지시 (2026-08-04): "제안값으로 채우고 리포트 먼저 보여줘"
→ 스케치·논문·사용자 확정에서 온 값은 tier 0/1/2, 내가 고른 값은 **tier 3** 로 표시한다.

    $PY scratch/propose_3cases.py
"""
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scipy import integrate  # noqa: E402

ETA = 0.851e-3
KT = 1.380649e-23 * 300
PI = math.pi
HEX_NN = math.sqrt(2 / math.sqrt(3))


def hr(t):
    print(f"\n{'=' * 84}\n{t}\n{'=' * 84}")


# ══════════════════════════════════════════════════════════════════════
hr("A.  abp-rod-2d-run-tumble")
# ── 타원체 마찰 (scratch/perrin_friction.py 와 같은 계산, 구 극한 검증 완료) ──
def chi(axes):
    a = np.asarray(axes, float)
    dl = lambda l: math.sqrt((a[0]**2 + l) * (a[1]**2 + l) * (a[2]**2 + l))
    sc = float(a.max())**2
    sub = lambda f: integrate.quad(lambda u: f(sc * (1 / u - 1)) * sc / u**2,
                                   1e-14, 1.0, limit=400)[0]
    return sub(lambda l: 1 / dl(l)), [sub(lambda l, i=i: 1 / ((a[i]**2 + l) * dl(l)))
                                      for i in range(3)]


a1, a2 = 1.0e-6, 0.25e-6                       # 사용자 확정: 장축 2µm, 단축 500nm
chi0, chis = chi((a1, a2, a2))
zt = [16 * PI * ETA / (chi0 + x**2 * c) for x, c in zip((a1, a2, a2), chis)]
zr_perp = (16 * PI * ETA * (a1**2 + a2**2)) / (3 * (a1**2 * chis[0] + a2**2 * chis[1]))
d_eq = 2 * (a1 * a2 * a2) ** (1 / 3)
gamma_bar = 2.0 / (1 / zt[0] + 1 / zt[1])      # 2D 면내 조화평균
D_bar = KT / gamma_bar
D_r = KT / zr_perp
tau_r = 1 / D_r
TAU_TUMBLE = 0.5                                # 사용자 확정
V = 5e-6                                        # 스케치 상한
tau_eff = 1 / (1 / tau_r + 1 / TAU_TUMBLE)      # 균등 무작위 재배향 (제안, tier 3)
l_p = V * tau_eff
tau_B_abp = d_eq**2 / D_bar
tau_v = d_eq / V
RHO_P = 1050.0                                  # 제안 tier 3 (폴리스티렌류)
m_ell = RHO_P * (4 * PI / 3) * a1 * a2 * a2
tau_p_abp = m_ell / gamma_bar

print(f"  확정   반축 (1.0, 0.25, 0.25) µm · p=4 · 물@300K · τ_tumble=0.5s · v≤5µm/s")
print(f"  유도   d_eq = {d_eq*1e6:.4f} µm   γ̄(2D) = {gamma_bar:.4e} kg/s"
      f"   D̄ = {D_bar*1e12:.4f} µm²/s")
print(f"         γ_r,z = {zr_perp:.4e} kg·m²/s   D_r = {D_r:.4f} 1/s   τ_r = {tau_r:.4f} s")
print(f"         τ_B(d_eq) = {tau_B_abp:.4f} s   τ_v = d_eq/v = {tau_v:.4f} s"
      f"   τ_p = {tau_p_abp*1e6:.4f} µs")
print(f"  제안   텀블 = 균등 무작위 재배향(2D)  → τ_eff = {tau_eff:.4f} s")
print(f"         ℓ_p = v·τ_eff = {l_p*1e6:.4f} µm = {l_p/d_eq:.3f} d_eq = {l_p/(2*a1):.3f} 몸길이")
print(f"         Pe = v d_eq/D̄ = {V*d_eq/D_bar:.3f}    D_r* = D_r τ_B = {D_r*tau_B_abp:.3f}")
N_ABP, L_ABP_D = 1000, 32.0                     # 제안 tier 3
L_abp = L_ABP_D * d_eq
dt_abp = 1e-2 * min(tau_v, tau_r, TAU_TUMBLE)
print(f"         N = {N_ABP} (독립 앙상블, 1-A와 같은 방식) · L = {L_ABP_D:.0f} d_eq"
      f" = {L_abp*1e6:.2f} µm")
print(f"  dt     가장 빠른 관련 척도 = τ_v = {tau_v:.4f} s  →  dt = {dt_abp*1e3:.4f} ms"
      f" = {dt_abp/tau_B_abp:.3e} τ_B")
T_ABP = 200 * max(tau_r, TAU_TUMBLE)
print(f"  T_obs  {T_ABP:.0f} s = {T_ABP/max(tau_r,TAU_TUMBLE):.0f}×max(τ_r,τ_tumble)"
      f"  →  {T_ABP/dt_abp:.2e} steps")
print(f"  검사   ℓ_p/(L/4) = {l_p/(L_abp/4):.4f}  ·  τ_p/τ_v = {tau_p_abp/tau_v:.3e}"
      f"  ·  dt·D_r = {dt_abp*D_r:.3e}  ·  dt/τ_tumble = {dt_abp/TAU_TUMBLE:.3e}")

# ══════════════════════════════════════════════════════════════════════
hr("B.  trap-drag-2d-hex300")
d_td = 5.0e-6
gam_td = 3 * PI * ETA * d_td
Dt_td = KT / gam_td
tauB_td = d_td**2 / Dt_td
k_t = 1e-5                                      # 10 pN/µm (스케치)
tau_k_td = gam_td / k_t
l_k = math.sqrt(KT / k_t)
v_x = 0.5e-6
A_TD, PHI_TD, N_TD = 100.0, 0.35, 300           # A·φ 제안 tier 3 (1-B에서 육방 결정 확인된 조합)
a_mean_td = d_td * math.sqrt(PI / (4 * PHI_TD))
a_nn_td = HEX_NN * a_mean_td
L_td = a_mean_td * math.sqrt(N_TD)
r_c_td = 5 * a_mean_td
Gamma_td = A_TD / (a_mean_td / d_td) ** 3
# 페어 국소 강성 (1-B와 같은 정의): U''(r) = 12A kT d³/r⁵ + WCA
r_min_td = a_nn_td / d_td - 3 * math.sqrt(2 / (3 * (12 * A_TD / (a_nn_td/d_td)**5
                                                   - 3 * A_TD / (a_nn_td/d_td)**5)))
Upp_pair = 12 * A_TD / r_min_td**5              # [kT/d²]
tau_int_td = tauB_td / Upp_pair
print(f"  확정   d=5µm(관례) · k_t=10pN/µm · v_x=0.5µm/s · N~300 (전부 스케치)")
print(f"         페어 = A/r³ + WCA (사용자 확정 2026-08-04)")
print(f"  제안   A = {A_TD:.0f} · φ = {PHI_TD} → Γ = {Gamma_td:.2f}"
      f"   (1-B에서 이 Γ가 육방 결정 ψ₆=0.885 로 확인됨)")
print(f"         a_mean = {a_mean_td*1e6:.3f} µm   a_NN = {a_nn_td*1e6:.3f} µm"
      f"   L = {L_td*1e6:.2f} µm = {L_td/d_td:.1f} d")
print(f"         r_c = 5 a_mean = {r_c_td*1e6:.2f} µm    최소이미지 여유 "
      f"{(L_td/2)/r_c_td:.2f}×")
print(f"  유도   γ = {gam_td:.4e} kg/s  D_t = {Dt_td*1e12:.4f} µm²/s  τ_B = {tauB_td:.2f} s")
print(f"         τ_k = γ/k_t = {tau_k_td*1e3:.3f} ms   ℓ_k = {l_k*1e9:.2f} nm")
print(f"         τ_int(페어, r_min={r_min_td:.3f}d) = {tau_int_td*1e3:.3f} ms")
print(f"         τ_v = d/v_x = {d_td/v_x:.1f} s   Δr_ss = γv/k = {gam_td*v_x/k_t*1e9:.3f} nm")
dt_td = 1e-2 * min(tau_k_td, tau_int_td)
binder = "트랩 τ_k" if tau_k_td < tau_int_td else "페어 τ_int"
print(f"  dt     ★ 두 강성이 경쟁한다 — 더 빠른 쪽이 {binder}")
print(f"         dt = 10⁻²×{binder.split()[1]} = {dt_td*1e6:.2f} µs = {dt_td/tauB_td:.3e} τ_B")
T_TD = L_td / v_x
print(f"  T_obs  박스 횡단 = L/v = {T_TD:.0f} s  →  {T_TD/dt_td:.2e} steps")
print(f"  검사   SNR = Δr_ss/ℓ_k = {gam_td*v_x/k_t/l_k:.4f}  ← 신호가 잡음의 1/10 (기존 경고)")

# ══════════════════════════════════════════════════════════════════════
hr("C.  chain-bend-2d-oscill")
d_cb = 1.47e-6                                  # 논문 [P1][P2] — 스케치 5µm과 불일치
gam_cb = 3 * PI * ETA * d_cb
Dt_cb = KT / gam_cb
tauB_cb = d_cb**2 / Dt_cb
tau_k_cb = gam_cb / k_t
KAPPA0 = 64e-3                                  # [P1] 10 mM MgCl₂, tier 2
a_cb = d_cb / 2
EI = KAPPA0 * a_cb**3 / 3
kappa_theta = EI / d_cb
# ★ n=11 로 잡았더니 δ_max = 74 nm 로 진폭 창이 거의 닫혔다 (ℓ_k=20nm 바로 위).
#   δ_max = M_c L²/(12 EI) ∝ n² 이므로 사슬을 늘리면 창이 열린다.
#   [P1] Fig.4 는 25입자 사슬을 쓴다 → n=25 로 제안.
N_BEADS = 25                                    # [P1] Fig.4 (25-particle aggregate)
L_cb = (N_BEADS - 1) * d_cb
k_end = 24 * EI / L_cb**3
k_center = 2 * k_end
MC = 35e-18                                     # [P1] tier 2
# 결합(반경) 강성: 열적 신축을 지름의 0.1% 이하로
K_BOND = KT / (1e-3 * d_cb) ** 2
print(f"  논문   d = {d_cb*1e6:.2f} µm · E=3100MPa · ν=0.4 · κ₀={KAPPA0*1e3:.0f} mN/m"
      f" · M_c={MC*1e18:.0f} pN·µm  (tier 2)")
print(f"  유도   EI = {EI:.4e} N·m²   κ_θ = EI/ℓ = {kappa_theta:.4e} J"
      f" = {kappa_theta/KT:.3e} kT")
print(f"         n={N_BEADS} → L = {L_cb*1e6:.2f} µm  κ_end = {k_end*1e6:.2f} pN/µm"
      f"  κ_center = {k_center*1e6:.2f} pN/µm")
print(f"         γ = {gam_cb:.4e} kg/s  τ_B = {tauB_cb:.3f} s  τ_k = {tau_k_cb*1e3:.3f} ms")
print(f"  제안   k_bond = {K_BOND*1e6:.0f} pN/µm (열적 신축 ≤ 0.1% d) · tier 3")

# ★ dt: 강성 행렬의 최대 고유값에서 (가장 빠른 모드가 dt를 정한다)
n = N_BEADS
B = np.zeros((n - 2, n))
for i in range(n - 2):
    B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
B /= d_cb
A_bend = kappa_theta * (B.T @ B)                 # 횡방향 굽힘
G = np.zeros((n - 1, n))
for i in range(n - 1):
    G[i, i], G[i, i + 1] = -1.0, 1.0
G /= 1.0
A_bond = K_BOND * (G.T @ G)                      # 종방향 신축
lam_bend = float(np.linalg.eigvalsh(A_bend).max())
lam_bond = float(np.linalg.eigvalsh(A_bond).max())
lam = max(lam_bend, lam_bond)
tau_fast = gam_cb / lam
dt_cb = 1e-2 * tau_fast
print(f"\n  ★ dt 는 **강성 행렬의 최대 고유값**에서 뽑는다 (가장 빠른 모드가 정한다)")
print(f"    λ_max(굽힘) = {lam_bend:.4e} N/m  → τ = {gam_cb/lam_bend*1e6:.4f} µs")
print(f"    λ_max(신축) = {lam_bond:.4e} N/m  → τ = {gam_cb/lam_bond*1e6:.4f} µs")
print(f"    지배: {'굽힘' if lam_bend > lam_bond else '신축'}  →  τ_fast = {tau_fast*1e6:.4f} µs")
print(f"    dt = 10⁻² τ_fast = {dt_cb*1e9:.2f} ns = {dt_cb/tauB_cb:.3e} τ_B")
print(f"    참고: 집단 굽힘 모드(κ_center) 는 {gam_cb/k_center*1e6:.1f} µs — "
      f"최속 모드보다 {(gam_cb/k_center)/tau_fast:.0f}배 느리다")

# 진폭 창
F_max = 2 * MC / L_cb
d_max = F_max / k_end
print(f"\n  진폭 창  ℓ_k = {l_k*1e9:.1f} nm ≪ a < δ_max = {d_max*1e9:.0f} nm"
      f"   (M<M_c: F<{F_max*1e12:.2f} pN)")
print(f"  {'n':>4} {'L':>9} {'κ_end':>11} {'δ_max':>9}   진폭 창")
for nb in (11, 15, 25, 41):
    Lb = (nb - 1) * d_cb
    ke = 24 * EI / Lb**3
    dm = MC * Lb**2 / (12 * EI)
    print(f"  {nb:>4} {Lb*1e6:7.1f}µm {ke*1e6:9.2f}pN/µm {dm*1e9:7.0f}nm   "
          f"{'닫힘 (ℓ_k 바로 위)' if dm < 4*l_k else f'20nm ≪ a < {dm*1e9:.0f}nm'}")
A_AMP = 200e-9
print(f"  제안     a = {A_AMP*1e9:.0f} nm  → a/ℓ_k = {A_AMP/l_k:.1f} (SNR)"
      f"  a/δ_max = {A_AMP/d_max:.3f} (선형 여유 {d_max/A_AMP:.1f}×) · tier 3")

# 주파수 창 + 비용
print(f"\n  주파수  De = ωτ 의 τ 는 무엇인가 — 트랩(τ_k) vs 사슬 집단모드")
tau_chain = gam_cb / k_center
for lbl, tau_ in (("τ_k(트랩)", tau_k_cb), ("τ_chain(집단)", tau_chain)):
    print(f"    {lbl:<16} = {tau_*1e3:8.4f} ms  →  De=1 at ω = {1/tau_:9.1f} rad/s"
          f" ({1/tau_/(2*PI):8.1f} Hz)")
OMEGA_LO, OMEGA_HI = 0.1 / tau_k_cb, 10 / tau_chain
print(f"  제안     ω 스윕 {OMEGA_LO:.0f} ~ {OMEGA_HI:.0f} rad/s (De 0.1~10 을 덮음) · tier 3")
n_cyc = 100
T_cb = n_cyc * 2 * PI / OMEGA_LO
print(f"\n  ★ 비용  최저 ω에서 {n_cyc}주기 = {T_cb:.2f} s  →  {T_cb/dt_cb:.2e} steps")
print(f"    1-B 실측 처리량(N=400, 이웃 79개)이 ~6600 steps/s 였다. 여기는 N={N_BEADS}로")
print(f"    훨씬 작아 빠르겠지만, 스텝 수 자체가 {T_cb/dt_cb:.1e} 개다.")
print(f"    → **논문 강성을 그대로 쓰면 직접 BD가 비싸다.** 이유: κ_θ = {kappa_theta/KT:.1e} kT 로")
print(f"      사슬이 열적으로 완전히 뻣뻣해서, 측정과 무관한 최속 모드가 dt를 정한다.")
print(f"    선택지: (a) 그대로 감수  (b) κ₀를 낮춰(계면활성제 조건, [P2] Fig.4) 스윕")
print(f"            (c) 높은 ω만 재고 낮은 ω는 준정적 극한으로 대체")

hr("요약 — system.yaml 에 넣을 tier")
print("""  tier 0  스케치에 적힌 값        k_t · v_x · N~300 · τ_R=0.5s · v≤5µm/s · A={0.1,1,10,100}
  tier 1  사용자 확정 / 관례 승계    타원체 형상 · 매질(물 300K) · d=5µm 관례 · 페어=A/r³+WCA
  tier 2  논문 추출 (미검증)         d=1.47µm · E · ν · κ₀=64mN/m · M_c=35pN·µm
  tier 3  ★ 내 제안                텀블 각도분포 · N_abp=1000 · L · A=100·φ=0.35 (trap-drag)
                                    n_beads=11 · k_bond · 진폭 200nm · ω 스윕 범위 · ρ_p""")
