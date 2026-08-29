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
    """K*[kT/d^2] -> kappa_0[N/m].

    kappa = 48EI/L^3 (rigidly clamped three-point bending), kappa_0 = 3EI/a^3
    """
    L = (n-1)*ELL
    EI = K_star*KTD2*L**3/48
    return 3*EI/a**3

# ── DLVO omega sweep (a=632nm, default k_t) ──
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

# ── (1) kappa_0(omega) -- the legitimate conversion ──
A = ax[0]
A.errorbar(ws, np.abs(kappa0(kp))*1e3, yerr=kappa0(kpe)*1e3, fmt="o-",
           color="#1f77b4", capsize=4, label="DLVO-only, omega sweep (default k_t, a=632nm)")
A.errorbar([3000], [abs(kappa0(Dk))*1e3], yerr=[kappa0(Dke)*1e3], fmt="v",
           color="#2ca02c", ms=13, capsize=5, label="DLVO-only (k_t x100, a=1d) <- 100% tracking")
A.errorbar([3000], [kappa0(Jk)*1e3], yerr=[kappa0(Jke)*1e3], fmt="s",
           color="#d62728", ms=13, capsize=5, label="DLVO+JKR (k_t×100, a=1d)")
A.axhline(64, color="k", ls="--", lw=1.8, label="[P1] measured kappa_0 = 64 mN/m (10 mM MgCl2)")
A.set_xscale("log"); A.set_yscale("log")
A.set_xlabel("ω [rad/s]"); A.set_ylabel("|κ₀| [mN/m]")
A.set_ylim(1e-4, 3e2)
A.set_title("(1) the legitimate conversion --  K' -> EI -> kappa_0 (the papers' "
            "tangential spring constant)\n"
            "JKR's 56.6 mN/m is **recovering the input 64** (a round trip; the -12% "
            "is the finite trap plus De=10.7)\n"
            "★ DLVO is 4-5 orders of magnitude below -- that IS the result",
            fontsize=10)
A.legend(fontsize=7.5, loc="lower right"); A.grid(alpha=0.3, which="both")

# ── (2) nominal GSER values -- NOT justified ──
B = ax[1]
lab = ["DLVO-only\n(k_t×100)", "DLVO+JKR\n(k_t×100)"]
x = np.arange(2); w_ = 0.35
gp = [Dk*KTD2/(6*np.pi*a), Jk*KTD2/(6*np.pi*a)]
gpp= [Dk2*KTD2/(6*np.pi*a), Jk2*KTD2/(6*np.pi*a)]
B.bar(x-w_/2, np.abs(gp), w_, color="#1f77b4", label="|G'| nominal")
B.bar(x+w_/2, np.abs(gpp), w_, color="#ff7f0e", label="|G''| nominal")
B.axhline(eta*3000, color="k", ls=":", lw=1.8, label=f"pure water G''=eta*omega={eta*3000:.2f} Pa (G'=0)")
B.set_yscale("log"); B.set_xticks(x); B.set_xticklabels(lab, fontsize=9)
B.set_ylabel("|G| [Pa]  (log)")
for xi, v in zip(x-w_/2, gp): B.text(xi, abs(v)*1.5, f"{v:.3f}", ha="center", fontsize=8)
for xi, v in zip(x+w_/2, gpp): B.text(xi, abs(v)*1.5, f"{v:.3f}", ha="center", fontsize=8)
B.set_title("(2) nominal values from forcing GSER on:  G* = K*/(6*pi*a)\n"
            "★★ NOT justified -- GSER holds when the probe is **immersed in a "
            "continuous medium**.\n"
            "What the bead feels here is not a medium but **its two neighbouring "
            "beads**", fontsize=10)
B.legend(fontsize=8); B.grid(alpha=0.3, axis="y", which="both")

fig.suptitle("Can G' and G'' be estimated -- what is legitimate and what is not",
             fontsize=12)
fig.tight_layout()
out="/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/moduli.png"
fig.savefig(out, dpi=140); print("saved", out)
