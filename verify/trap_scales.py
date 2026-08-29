"""trap-2d-5um 스케일 표 초안 — pint로 단위 붙여 계산."""
import pint

u = pint.UnitRegistry()
Q = u.Quantity

kB = Q(1.380649e-23, "J/K")

# ── 스케치에서 읽은 값 (차원 있음) ────────────────────────────────────
k_t = Q(10, "pN/um").to("N/m")
T = Q(300, "K")
kT = (kB * T).to("J")

# ── 스케치에 없어서 채워야 하는 값 ────────────────────────────────────
eta = Q(0.851, "mPa*s").to("Pa*s")     # 물 @300K (핸드북). 매질 미기재 — 가정
rho_p = Q(2000, "kg/m^3")              # 실리카 가정. τ_p 계산에만 쓰임


def table(d, label):
    d = d.to("m")
    gamma = (3 * 3.141592653589793 * eta * d).to("kg/s")
    D_t = (kT / gamma).to("um^2/s")
    tau_B = (d**2 / D_t).to("s")
    tau_k = (gamma / k_t).to("s")
    l_k = ((kT / k_t) ** 0.5).to("nm")
    x2 = (kT / k_t).to("um^2")
    m = (rho_p * (3.141592653589793 / 6) * d**3).to("kg")
    tau_p = (m / gamma).to("s")
    k_star = float((k_t * d**2 / kT).to("dimensionless"))

    print(f"\n{'='*74}\n{label}   d = {d.to('um'):~.3fP}\n{'='*74}")
    print(f"  γ      = {gamma:~.4eP}")
    print(f"  D_t    = {D_t:~.4fP}")
    print(f"  m      = {m:~.3eP}   (ρ_p = {rho_p:~P} 가정)")
    print()
    print(f"  시간척도 (작은 것부터)")
    print(f"    τ_p = m/γ    = {tau_p.to('us'):~.3fP}   관성 이완")
    print(f"    τ_k = γ/k    = {tau_k.to('ms'):~.3fP}   ★ 트랩 이완 — 지배 시간척도")
    print(f"    τ_B = d²/D_t = {tau_B:~.1fP}  확산 (트랩이 있어 실현 안 됨)")
    print()
    print(f"  길이척도")
    print(f"    ℓ_k = √(kT/k) = {l_k:~.2fP}   트랩 안 요동 폭 (d와 무관!)")
    print(f"    d              = {d.to('um'):~.1fP}")
    print(f"    ℓ_k/d          = {float((l_k/d).to('dimensionless')):.2e}")
    print()
    print(f"  무차원수")
    print(f"    k* = k d²/kT   = {k_star:.3e}   트랩 vs 열요동")
    print(f"    τ_p/τ_k        = {float((tau_p/tau_k).to('dimensionless')):.3e}   관성 vs 트랩")
    print(f"    τ_k/τ_B        = {float((tau_k/tau_B).to('dimensionless')):.3e}")
    print()
    print(f"  해석해 (골든 테스트)")
    print(f"    자유도당 ⟨x²⟩ = kT/k = {x2:~.5eP}  → rms {l_k:~.2fP}")
    print(f"    2D ⟨r²⟩       = 2kT/k = {(2*x2):~.5eP}")
    print(f"    이완시간 τ    = γ/k = {tau_k.to('ms'):~.3fP}")

    # dt 창
    dt_max = 0.01 * tau_k
    print()
    print(f"  dt 선택")
    print(f"    트랩 해상 dt ≤ 0.01 τ_k = {dt_max.to('us'):~.2fP}")
    print(f"    권장 dt = τ_k/2000      = {(tau_k/2000).to('us'):~.3fP}   (골든 검증에서 쓴 값)")
    print(f"    τ_p와 비교: τ_p = {tau_p.to('us'):~.3fP}")
    ratio = float((tau_p / (tau_k / 2000)).to("dimensionless"))
    print(f"    → τ_p/dt = {ratio:.2f}")
    if ratio > 0.01:
        print(f"       ⚠️ bd-physics의 '관성 무시: τ_p/dt ≤ 1e-2' 검사에 걸림 (아래 논의)")
    return dict(gamma=gamma, tau_k=tau_k, tau_p=tau_p, l_k=l_k, k_star=k_star)


print(f"스케치에서 읽은 값:  k_t = {k_t:~.3eP}   T = {T:~P}   kT = {kT:~.4eP}")
print(f"채운 값:            η = {eta:~P} (물@300K, 가정)   ρ_p = {rho_p:~P} (실리카, 가정)")

A = table(Q(10, "um"), "해석 A — R=5µm 이 입자 반지름  → d = 10 µm")
B = table(Q(5, "um"), "해석 B — 입자 지름이 5 µm       → d = 5 µm")

print(f"\n{'='*74}")
print("두 해석의 차이")
print(f"{'='*74}")
print(f"  ⟨x²⟩ (해석해)  : 동일 — kT/k 는 d에 무관")
print(f"  τ_k (이완시간) : {A['tau_k'].to('ms'):~.3fP}  vs  {B['tau_k'].to('ms'):~.3fP}   (정확히 2배)")
print(f"  k*             : {A['k_star']:.3e}  vs  {B['k_star']:.3e}")
print("  → 평형 결과는 같고 시간척도만 2배 차이. 골든 검증에는 영향 없음.")
