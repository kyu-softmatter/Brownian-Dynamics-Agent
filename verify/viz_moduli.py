import json, glob, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

kT = 1.380649e-23*300; d = 1.47e-6; a = d/2; eta = 0.851e-3
KTD2 = kT/d**2
ELL = 1.0076*d

def grab(pat, key, minprod=None):
    v=[]
    for dd in sorted(glob.glob(pat)):
        try:
            m=json.load(open(f"{dd}/metrics.json")); s=json.load(open(f"{dd}/spec.json"))
        except FileNotFoundError: continue
        if minprod and s["numerics"]["n_prod"]<minprod: continue
        o={x["name"]:x["measured"] for x in m["observables"]}
        if o.get(key) is not None: v.append(o[key])
    return np.array(v)
def ms(x): return (x.mean(), x.std(ddof=1)/np.sqrt(len(x))) if len(x)>1 else (x.mean() if len(x) else np.nan, np.nan)

def kappa0(K_star, n=9):
    """K*[kT/d²] → κ₀[N/m].  κ=48EI/L³ (강체고정 3점굽힘) · κ₀=3EI/a³"""
    L = (n-1)*ELL
    EI = K_star*KTD2*L**3/48
    return 3*EI/a**3

# ── DLVO ω 스윕 (a=632nm, k_t 기본) ──
ws, kp, kpe = [], [], []
for w in (10,100,300,1000,3000,10000,30000):
    v = grab(f"runs/chain-bend-2d-dlvo__n9-w{w}-a632__*", "K_prime")
    if not len(v): continue
    m,e = ms(v); ws.append(w); kp.append(m); kpe.append(e)
ws, kp, kpe = np.array(ws,float), np.array(kp), np.array(kpe)

# ── kt100, a=1d, ω=3000 ──
Pk = "runs/chain-bend-2d-dlvo__n9-w3000-a1470"
Dk, Dke = ms(grab(f"{Pk}-kt100__*","K_prime"))
Jk, Jke = ms(grab(f"{Pk}-jkr-kt100__*","K_prime"))
Dk2,_   = ms(grab(f"{Pk}-kt100__*","K_doubleprime"))
Jk2,_   = ms(grab(f"{Pk}-jkr-kt100__*","K_doubleprime"))

fig, ax = plt.subplots(1, 2, figsize=(14, 5.4))

# ── ① κ₀(ω) — 정당한 환산 ──
A = ax[0]
A.errorbar(ws, np.abs(kappa0(kp))*1e3, yerr=kappa0(kpe)*1e3, fmt="o-",
           color="#1f77b4", capsize=4, label="DLVO-only, ω 스윕 (k_t 기본, a=632nm)")
A.errorbar([3000], [abs(kappa0(Dk))*1e3], yerr=[kappa0(Dke)*1e3], fmt="v",
           color="#2ca02c", ms=13, capsize=5, label="DLVO-only (k_t×100, a=1d) ← 추종 100%")
A.errorbar([3000], [kappa0(Jk)*1e3], yerr=[kappa0(Jke)*1e3], fmt="s",
           color="#d62728", ms=13, capsize=5, label="DLVO+JKR (k_t×100, a=1d)")
A.axhline(64, color="k", ls="--", lw=1.8, label="[P1] 실측 κ₀ = 64 mN/m (10 mM MgCl₂)")
A.set_xscale("log"); A.set_yscale("log")
A.set_xlabel("ω [rad/s]"); A.set_ylabel("|κ₀| [mN/m]")
A.set_ylim(1e-4, 3e2)
A.set_title("① 정당한 환산 —  K′ → EI → κ₀ (논문의 접선 스프링상수)\n"
            "JKR 56.6 mN/m 은 **입력 64 를 되찾은 것**(왕복확인, −12%는 유한트랩+De=10.7)\n"
            "★ DLVO 는 4~5 자릿수 아래 — 이게 결과다", fontsize=10)
A.legend(fontsize=7.5, loc="lower right"); A.grid(alpha=0.3, which="both")

# ── ② GSER 명목값 — 정당화되지 않음 ──
B = ax[1]
lab = ["DLVO-only\n(k_t×100)", "DLVO+JKR\n(k_t×100)"]
x = np.arange(2); w_ = 0.35
gp = [Dk*KTD2/(6*np.pi*a), Jk*KTD2/(6*np.pi*a)]
gpp= [Dk2*KTD2/(6*np.pi*a), Jk2*KTD2/(6*np.pi*a)]
B.bar(x-w_/2, np.abs(gp), w_, color="#1f77b4", label="|G′| 명목")
B.bar(x+w_/2, np.abs(gpp), w_, color="#ff7f0e", label="|G″| 명목")
B.axhline(eta*3000, color="k", ls=":", lw=1.8, label=f"순수 물 G″=ηω={eta*3000:.2f} Pa (G′=0)")
B.set_yscale("log"); B.set_xticks(x); B.set_xticklabels(lab, fontsize=9)
B.set_ylabel("|G| [Pa]  (로그)")
for xi, v in zip(x-w_/2, gp): B.text(xi, abs(v)*1.5, f"{v:.3f}", ha="center", fontsize=8)
for xi, v in zip(x+w_/2, gpp): B.text(xi, abs(v)*1.5, f"{v:.3f}", ha="center", fontsize=8)
B.set_title("② GSER 를 억지로 적용한 명목값  G*=K*/(6πa)\n"
            "★★ 정당화되지 않는다 — GSER 는 프로브가 **연속 매질에 잠겨** 있을 때 성립.\n"
            "여기서 비드가 느끼는 것은 매질이 아니라 **이웃 비드 2개**다", fontsize=10)
B.legend(fontsize=8); B.grid(alpha=0.3, axis="y", which="both")

fig.suptitle("G′·G″ 를 추정할 수 있는가 — 무엇이 정당하고 무엇이 아닌가", fontsize=12)
fig.tight_layout()
out="/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/moduli.png"
fig.savefig(out, dpi=140); print("saved", out)
