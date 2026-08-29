import json, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

def grab(pat, key):
    v = []
    for d in sorted(glob.glob(pat)):
        try: m = json.load(open(f"{d}/metrics.json"))
        except FileNotFoundError: continue
        o = {x["name"]: x["measured"] for x in m["observables"]}
        if o.get(key) is not None: v.append(o[key])
    return np.array(v)
def ms(a):
    if len(a) > 1: return a.mean(), a.std(ddof=1)/np.sqrt(len(a)), len(a)
    return (a.mean() if len(a) else np.nan), np.nan, len(a)

P = "runs/chain-bend-2d-dlvo__n9-w3000-a1470"
CONDS = [("k_t default\n(trap)",   f"{P}__*",        f"{P}-jkr__*"),
         ("k_t ×100\n(trap)",   f"{P}-kt100__*",  f"{P}-jkr-kt100__*")]

fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))

# ── (1) tracking ratio ──
a0 = ax[0]
x = np.arange(len(CONDS)); w = 0.35
for off, br, col, lab in [(-w/2, 1, "#1f77b4", "DLVO-only"), (w/2, 2, "#d62728", "DLVO+JKR")]:
    m = [ms(grab(c[br], "y_response"))[0] for c in CONDS]
    e = [ms(grab(c[br], "y_response"))[1] for c in CONDS]
    a0.bar(x+off, m, w, yerr=e, capsize=4, color=col, label=lab)
    for xi, v in zip(x+off, m):
        if np.isfinite(v): a0.text(xi, v+0.03, f"{100*v:.0f}%", ha="center", fontsize=9)
a0.axhline(1.0, color="gray", ls=":", lw=1.2, label="perfect tracking")
a0.set_xticks(x); a0.set_xticklabels([c[0] for c in CONDS], fontsize=9)
a0.set_ylabel("|y_hat| / a  (tracking ratio)"); a0.set_ylim(0, 1.18)
a0.set_title("(1) making the trap 100x stiffer recovered the tracking ratio\n"
             "design prediction 82% vs measured 81.7% (JKR)", fontsize=10)
a0.legend(fontsize=8); a0.grid(alpha=0.3, axis="y")

# ── (2) K' (log; DLVO is near 0 so it is annotated separately) ──
a1 = ax[1]
for off, br, col, lab in [(-w/2, 1, "#1f77b4", "DLVO-only"), (w/2, 2, "#d62728", "DLVO+JKR")]:
    m = [ms(grab(c[br], "K_prime"))[0] for c in CONDS]
    e = [ms(grab(c[br], "K_prime"))[1] for c in CONDS]
    mm = [abs(v) if np.isfinite(v) and v != 0 else np.nan for v in m]
    a1.bar(x+off, mm, w, yerr=e, capsize=4, color=col, label=lab)
    for xi, v, ee in zip(x+off, m, e):
        if np.isfinite(v):
            z = abs(v/ee) if (np.isfinite(ee) and ee) else np.nan
            txt = f"{v:.4g}" + (f"\n({z:.1f}σ from 0)" if np.isfinite(z) else "")
            a1.text(xi, abs(v)*1.6, txt, ha="center", fontsize=7.5)
a1.set_yscale("log"); a1.set_xticks(x); a1.set_xticklabels([c[0] for c in CONDS], fontsize=9)
a1.set_ylabel("|K'| [kT/d^2]  (log)")
a1.set_ylim(1, 1e7)
a1.axhline(1.391e6, color="k", ls=":", lw=1.4, label="JKR bending stiffness kappa_theta*=1.4e6")
a1.set_title("★ (2) as the trap stiffened, DLVO K' collapsed to 0\n"
             "-9.2±10.9 = **0.84 sigma from 0** (indistinguishable). JKR is a clear "
             "1.13e5", fontsize=10)
a1.legend(fontsize=8); a1.grid(alpha=0.3, axis="y")

# ── (3) bow -- its discriminating power inverts with the condition ──
a2 = ax[2]
bows = []
for c in CONDS:
    d = grab(c[1], "bow_rms"); 
    if not len(d): d = grab(c[1], "bow_rms_gsd")
    j = grab(c[2], "bow_rms")
    if not len(j): j = grab(c[2], "bow_rms_gsd")
    bows.append((ms(d), ms(j)))
for off, idx, col, lab in [(-w/2, 0, "#1f77b4", "DLVO-only"), (w/2, 1, "#d62728", "DLVO+JKR")]:
    m = [b[idx][0] for b in bows]; e = [b[idx][1] for b in bows]
    a2.bar(x+off, m, w, yerr=e, capsize=4, color=col, label=lab)
    for xi, v in zip(x+off, m):
        if np.isfinite(v): a2.text(xi, v*1.15, f"{v:.4f}", ha="center", fontsize=8)
a2.set_yscale("log"); a2.set_xticks(x); a2.set_xticklabels([c[0] for c in CONDS], fontsize=9)
a2.set_ylabel("bow RMS [d]  (log)")
r = [b[0][0]/b[1][0] for b in bows]
a2.set_title(f"(3) the discriminating power of bow inverts with the condition\n"
             f"weak trap {r[0]:.1f}x -> stiff trap {r[1]:.2f}x\n"
             "once the deformation is imposed, shape carries no information",
             fontsize=10)
a2.legend(fontsize=8); a2.grid(alpha=0.3, axis="y")

fig.suptitle("DLVO vs JKR against trap stiffness (n=9, omega=3000 rad/s, "
             "a=1470nm=1d, 6 seeds; only 1 for default-k_t JKR)", fontsize=11)
fig.tight_layout()
out = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/kt100_results.png"
fig.savefig(out, dpi=140); print("saved", out)
