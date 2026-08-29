"""`chain-bend-2d-oscill` 결과 확인 — 그래프 + 애니메이션.

**결과를 말로만 보고하지 않는다** (작업 관행). 두 산출물을 만든다:

  scratch/_viz/chain_bend_results.png   6패널 — 관문 결과와 스펙 오류를 눈으로 확인
  scratch/_viz/chain_bend_motion.gif    사슬 운동 — kT=0(모드 형태) vs kT=1(왜 SNR 이 문제인지)

애니메이션은 **kT=0 결정론 + 큰 dt** 로 만든다 (dt·λ_max = 0.22, 안정 한계의 11%).
운동의 모양을 보이는 것이 목적이고 **생산 측정이 아니다** — 생산 dt 는 4.53e-10 이다.
열요동 진폭 ℓ_k = √(kT/k_t) 는 정적 평형량이라 이 dt 에서도 옳게 나온다 (dt·k_t/γ = 5e-5).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/viz_chain_bend.py
"""
from __future__ import annotations

import glob
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "verify" / "_viz"
GATES = ROOT / "verify" / "_gates"
sys.path.insert(0, str(ROOT))

C_MEAS, C_THEORY, C_BAD, C_GOOD, C_GREY = "#1f77b4", "#d62728", "#ff7f0e", "#2ca02c", "#888888"

# 한글 폰트. DejaVu Sans 에는 한글 글리프가 없어서 그냥 두면 라벨이 전부 □ 로 나온다.
# ★ 한글 폰트 대부분에 U+2212(−) 와 ŷ(U+0177) 글리프가 **없다** — 실측: AppleGothic ·
# Apple SD Gothic Neo · NanumGothic 셋 다 없음. Arial Unicode MS 만 한글+기호를 다 갖는다.
for _f in ("Arial Unicode MS", "Apple SD Gothic Neo", "AppleGothic", "NanumGothic"):
    if _f in {f.name for f in matplotlib.font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = _f
        break
plt.rcParams["axes.unicode_minus"] = False


def load_specs():
    """스펙이 있으면 읽고, 없으면 system.yaml 에서 직접 만든다.

    ★ L3 하드 검사가 실패하면 `nondim spec` 은 스펙을 **쓰지 않는다** (정상 동작).
    결과 확인 그림이 그것 때문에 못 그려지면 안 되므로 케이스 모듈로 되돌아간다.
    """
    sp = [json.loads(Path(q).read_text())
          for q in glob.glob(str(ROOT / "specs" / "chain-bend-2d-oscill__*.json"))]
    if sp:
        return sorted(sp, key=lambda s: s["params"]["omega_star"])
    import argparse as _ap
    from cases import chain_bend_2d as CB
    sys_ = CB.load_system(ROOT / "intake" / "chain-bend-2d-oscill" / "system.yaml")
    args = _ap.Namespace(dt_scale=1.0, cycles=CB.N_CYCLES, samples=2000)
    lo, hi = sys_["omega_range"]
    out = []
    for om in np.geomspace(lo, hi, CB.N_SWEEP):
        _, spec, _, _, _, _ = CB.build_spec(sys_, float(om), args)
        out.append(json.loads(json.dumps(spec.to_doc(), default=str))
                   if hasattr(spec, "to_doc") else
                   dict(params=spec.params, numerics=spec.numerics))
    return sorted(out, key=lambda s: s["params"]["omega_star"])


def chain_matrices(p):
    n = int(p["n_beads"]); trapped = sorted(int(t) for t in p["trapped"])
    kth = float(p["kappa_theta_star"]); kt = float(p["k_t_star"])
    ell = float(p["L_chain_star"]) / (n - 1)
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    A = kth * (B.T @ B)
    T = np.zeros((n, n))
    for e in trapped:
        T[e, e] += kt
    return n, trapped, kth, kt, ell, A, T


# ════════════════════════════════════════════════════════════════════════
# 그래프
# ════════════════════════════════════════════════════════════════════════
def make_figure():
    specs = load_specs()
    p = specs[0]["params"]
    n, trapped, kth, kt, ell, A, T = chain_matrices(p)
    amp = float(p["amp_star"]); mid = trapped[len(trapped) // 2]
    lk = math.sqrt(1.0 / kt)
    em = np.eye(n)[mid]
    ev = np.linalg.eigvalsh(A + T)
    kappa_center = 48 * (kth * ell) / float(p["L_chain_star"]) ** 3
    Kstat = kt * (1.0 / np.linalg.solve(A + T, kt * em)[mid] - 1.0)

    det, deq = {}, {}
    for f in glob.glob(str(GATES / "det_*.json")):
        r = json.load(open(f))
        if r["method"] == "ov":
            det[round(r["de"], 3)] = r
    for f in glob.glob(str(GATES / "deq_*.json")):
        r = json.load(open(f))
        deq[round(r["de"], 3)] = r

    om_all = np.array([s["params"]["omega_star"] for s in specs])
    de_all = np.array([s["params"]["De"] for s in specs])
    snr_all = np.array([s["params"]["snr_response"] for s in specs])

    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.2))
    fig.suptitle("chain-bend-2d-oscill — L3 관문 결과 (2026-08-05)  ·  "
                 "스펙 유도 5곳의 오류를 실행으로 잡아낸 기록", fontsize=13, y=0.98)

    # ① K*(ω) — 실측 vs 선형응답
    a = ax[0, 0]
    om_f = np.geomspace(om_all.min() * 0.7, om_all.max() * 1.4, 200)
    KA = np.array([kt * amp / np.linalg.solve(1j * w * np.eye(n) + A + T, kt * amp * em)[mid]
                   - kt - 1j * w for w in om_f])
    a.loglog(om_f, KA.real, "-", c=C_THEORY, lw=1.6, label="K′ 선형응답 (정확 최소화와 0.32%)")
    a.loglog(om_f, KA.imag, "--", c=C_THEORY, lw=1.6, label="K″ 선형응답")
    if deq:
        k = sorted(deq)
        a.loglog([deq[d]["omega"] for d in k], [deq[d]["K_re"] for d in k], "o",
                 c=C_MEAS, ms=8, label="K′ HOOMD (평형화 20τ_max)")
        a.loglog([deq[d]["omega"] for d in k], [deq[d]["K_im"] for d in k], "s",
                 c=C_MEAS, ms=7, mfc="none", label="K″ HOOMD")
    a.axhline(kappa_center, c=C_GREY, ls=":", lw=1.2)
    a.text(om_f[-1], kappa_center * 1.15, "48EI/L³ (강체 고정)", fontsize=8, c=C_GREY, ha="right")
    a.axhline(Kstat, c=C_GOOD, ls=":", lw=1.4)
    a.text(om_f[-1], Kstat * 1.15, f"κ_drive={Kstat:.0f} (트랩 경계)", fontsize=8,
           c=C_GOOD, ha="right")
    a.set_ylim(Kstat * 0.55, None)
    a.set_xlabel("ω*  [1/τ_B]"); a.set_ylabel("K*  [kT/d²]")
    a.set_title("① K*(ω) — 실측이 28% 낮다. 원인은 HOOMD angle.Harmonic 의\nsinθ 클램프 (힘만 틀림, 에너지는 정확)", fontsize=10)
    a.legend(fontsize=7.5, loc="upper left"); a.grid(alpha=0.25, which="both")

    # ② SNR — 스펙이 검사한 것 vs 실제
    a = ax[0, 1]
    a.loglog(de_all, snr_all, "o-", c=C_MEAS, ms=7, label="실제  |ŷ(ω)|/ℓ_k")
    a.axhline(amp / lk, c=C_BAD, ls="-", lw=2,
              label=f"스펙이 검사한 a/ℓ_k = {amp/lk:.2f} (ω 무관)")
    a.axhline(3, c=C_GOOD, ls="--", lw=1.3, label="검사 기준 3")
    a.axhline(1, c=C_GREY, ls=":", lw=1.3, label="열요동과 같아지는 선")
    for d, s in zip(de_all, snr_all):
        if s < 1:
            a.plot(d, s, "x", c=C_BAD, ms=11, mew=2.2)
    a.set_xlabel("De = ω τ_max"); a.set_ylabel("SNR")
    a.set_title("② SNR 검사가 분자를 틀렸다 — 최대 60배 과대평가\n"
                "(× = 응답이 열요동보다 작은 점, 7점 중 4점)", fontsize=10)
    a.legend(fontsize=7.5, loc="lower left"); a.grid(alpha=0.25, which="both")

    # ③ 이완 스펙트럼 + 스윕 범위
    a = ax[0, 2]
    a.semilogy(range(1, n + 1), ev, "o", c=C_MEAS, ms=5)
    a.axhline(ev[0], c=C_GOOD, ls="--", lw=1.4)
    a.text(1.4, ev[0] * 1.5, f"λ_min={ev[0]:.0f} → τ_max  ★지배 척도", fontsize=8, c=C_GOOD)
    a.axhline(kappa_center, c=C_BAD, ls="--", lw=1.4)
    a.text(1.4, kappa_center * 1.5, f"κ_center={kappa_center:.0f} → τ_chain (스펙이 쓴 것)",
           fontsize=8, c=C_BAD)
    a.axhline(ev[-1], c=C_GREY, ls=":", lw=1.2)
    a.text(1.4, ev[-1] * 0.35, f"λ_max={ev[-1]:.2e} → dt 를 정한다", fontsize=8, c=C_GREY)
    a.fill_between([1, n], om_all.min(), om_all.max(), color=C_MEAS, alpha=0.13)
    a.text(n * 0.97, om_all.max() * 1.6, "ω 스윕 범위", fontsize=8.5, c=C_MEAS, ha="right")
    a.set_xlabel("모드 번호"); a.set_ylabel("고유값 λ  [kT/d²]")
    a.set_title(f"③ 스윕이 준정적 영역에 못 들어간다\nτ_max/τ_chain = "
                f"{kappa_center/ev[0]:.2f}배", fontsize=10)
    a.grid(alpha=0.25, which="both")

    # ④ 관문 A — 공칭 진폭이 부호까지 틀린다
    a = ax[1, 0]
    de_A = [0.11, 0.23, 0.49, 1.04, 2.21, 4.70, 10.0]
    Kn = [4805.4, 4840.1, 4823.4, 4777.6, 4076.3, 2106.5, -6559.1]
    Km = [4807.8, 4846.7, 4863.7, 4913.2, 4765.0, 4712.8, 5863.2]
    Ks = [57.8, 62.6, 61.0, 77.9, 108.3, 493.4, 1715.2]
    ks_true = 4830.66
    a.semilogx(de_A, Kn, "o-", c=C_BAD, ms=7, label="공칭 진폭 a 사용 → 붕괴")
    a.errorbar(de_A, Km, yerr=Ks, fmt="s-", c=C_MEAS, ms=6, capsize=3,
               label="측정 위상자 ŷ_c 사용 → 평평")
    a.axhline(ks_true, c=C_GOOD, ls="--", lw=1.5, label=f"해석해 k_s = {ks_true:.0f}")
    a.axhline(0, c="k", lw=0.8)
    a.annotate("부호까지 틀림\n(오차 236%)", xy=(10.0, -6559), xytext=(1.5, -5200),
               fontsize=8.5, c=C_BAD, arrowprops=dict(arrowstyle="->", color=C_BAD, lw=1.2))
    a.set_xlabel("De (관문 A 단독 비드)"); a.set_ylabel("K′  [kT/d²]")
    a.set_title("④ 구동의 영차유지(ZOH) — 추정량에 공칭 진폭을\n쓰면 조용히 틀린다", fontsize=10)
    a.legend(fontsize=7.5, loc="lower left"); a.grid(alpha=0.25)

    # ⑤ 평형화 — σ 가 1000배 줄었다
    a = ax[1, 1]
    if det and deq:
        k = sorted(deq)
        x = np.arange(len(k)); w = 0.36
        a.bar(x - w/2, [det[d]["K_sem"] for d in k], w, color=C_BAD,
              label="평형화 5τ_chain = 0.54τ_max")
        a.bar(x + w/2, [max(deq[d]["K_sem"], 1e-2) for d in k], w, color=C_GOOD,
              label="평형화 20τ_max")
        a.set_yscale("log")
        a.set_xticks(x); a.set_xticklabels([f"{d:g}" for d in k])
        for i, d in enumerate(k):
            r = det[d]["K_sem"] / max(deq[d]["K_sem"], 1e-9)
            a.text(i - w/2, det[d]["K_sem"] * 1.7, f"×{r:,.0f}", ha="center", fontsize=8)
        a.set_ylim(None, max(det[d]["K_sem"] for d in k) * 12)
    a.set_xlabel("De (예전 정의)"); a.set_ylabel("블록 산포 σ(K*)  [kT/d²]")
    a.set_title("⑤ 평형화가 2.2τ_max 로 부족했다\nτ_max 로 고치니 산포가 ~1000배 감소", fontsize=10)
    a.legend(fontsize=7.5); a.grid(alpha=0.25, axis="y")

    # ⑥ 관문 B′ — 관성의 영향
    a = ax[1, 2]
    rows = []
    for f in glob.glob(str(GATES / "det_*.json")):
        rows.append(json.load(open(f)))
    pair = {}
    for r in rows:
        pair.setdefault(round(r["de"], 3), {})[r["method"]] = r
    dd = sorted([d for d, v in pair.items() if len(v) >= 2])
    if dd:
        dre = [100 * abs(pair[d]["ov"]["K_re"] - pair[d]["lang0"]["K_re"])
               / abs(pair[d]["ov"]["K_re"]) for d in dd]
        dim = [100 * abs(pair[d]["ov"]["K_im"] - pair[d]["lang0"]["K_im"])
               / abs(pair[d]["ov"]["K_im"]) for d in dd]
        a.loglog(dd, dre, "o-", c=C_MEAS, ms=7, label="K′ 차이")
        a.loglog(dd, dim, "s--", c=C_GOOD, ms=6, label="K″ 차이")
    a.axhline(47, c=C_BAD, ls="-", lw=2, label="열적 비교의 검정력 한계 (47%)")
    a.axhline(1, c=C_GREY, ls=":", lw=1.2, label="1% 기준")
    a.set_xlabel("De (예전 정의)"); a.set_ylabel("|과감쇠 − 관성| / K  [%]")
    a.set_title("⑥ τ_p/τ_fast=0.60 은 무해하다 — 최대 0.159%\n"
                "kT=0 결정론 차분이 열적 비교보다 300배 예민", fontsize=10)
    a.legend(fontsize=7.5, loc="upper left"); a.grid(alpha=0.25, which="both")

    fig.tight_layout(rect=[0, 0, 1, 0.955])
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "chain_bend_results.png"
    fig.savefig(out, dpi=125)
    plt.close(fig)
    return out


# ════════════════════════════════════════════════════════════════════════
# 애니메이션 — kT=0 (모드 형태) vs kT=1 (SNR 문제를 눈으로)
# ════════════════════════════════════════════════════════════════════════
def simulate_frames(kT, de_target, n_cycles=3, n_frames=90):
    """유령 트랩 사슬을 큰 dt 로 돌려 프레임을 모은다. 애니메이션 전용."""
    import gsd.hoomd, hoomd, hoomd.md as md

    specs = load_specs()
    sp = min(specs, key=lambda s: abs(s["params"]["De"] - de_target))
    p = sp["params"]
    n, trapped, kth, kt, ell, A, T = chain_matrices(p)
    amp = float(p["amp_star"]); omega = float(p["omega_star"])
    lam_max = float(np.linalg.eigvalsh(A + T)[-1])
    lam_min = float(np.linalg.eigvalsh(A + T)[0])
    dt = 0.22 / lam_max                      # 안정 한계 2/λ 의 11%

    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in trapped:
        pos.append(list(pos[g])); typeid.append(1)
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.configuration.box = [4 * float(p["L_chain_star"])] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(trapped)]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(trapped)
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=7)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=float(p["k_bond_star"]), r0=ell)
    bond.params["trap"] = dict(k=kt, r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=kth, t0=math.pi)
    filt = hoomd.filter.Type(["A"])
    meth = (md.methods.OverdampedViscous(filter=filt, default_gamma=1.0) if kT == 0
            else md.methods.Brownian(filter=filt, kT=kT, default_gamma=1.0))
    integ = md.Integrator(dt=dt, methods=[meth], forces=[bond, angle])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    mid = trapped[len(trapped) // 2]
    ghost_mid = n + trapped.index(mid)

    class Move(hoomd.custom.Action):
        def act(self, timestep):
            y = amp * math.sin(omega * timestep * dt)
            with self._state.cpu_local_snapshot as s:
                tg = np.array(s.particles.tag, copy=True)
                loc = np.flatnonzero(tg == ghost_mid)
                if len(loc):
                    s.particles.position[loc[0], 1] = y

    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=Move(), trigger=hoomd.trigger.Periodic(10)))

    n_eq = int(round(10.0 / (lam_min * dt)))          # 10 τ_max
    sim.run(n_eq)
    period = 2 * math.pi / omega
    total = int(round(n_cycles * period / dt))
    every = max(1, total // n_frames)
    frames, drive = [], []
    for _ in range(n_frames):
        sim.run(every)
        snap = sim.state.get_snapshot()               # 전역 스냅샷 = tag 순서
        q = np.array(snap.particles.position)
        frames.append(q[:n, :2].copy())
        drive.append(q[ghost_mid, 1])
    return dict(frames=np.array(frames), drive=np.array(drive), n=n, trapped=trapped,
                amp=amp, ell=ell, lk=math.sqrt(1.0 / kt),
                de=float(p["De"]), dt=dt, n_eq=n_eq)


def make_animation():
    hot = simulate_frames(1.0, 9.5)
    cold = simulate_frames(0.0, 9.5)
    nf = len(cold["frames"])
    ylim = max(np.abs(hot["frames"][:, :, 1]).max(), cold["amp"]) * 1.35

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 6.4), sharex=True)
    fig.suptitle(f"chain-bend 사슬 운동  ·  De = {cold['de']:.1f}  ·  "
                 f"구동 진폭 a = {cold['amp']:.3f} d  ·  ℓ_k = {cold['lk']:.4f} d",
                 fontsize=11.5)
    arts = []
    for a, dat, ttl, col in ((axes[0], cold, "kT = 0 (결정론) — 3점 굽힘 모드 형태", C_GOOD),
                             (axes[1], hot, "kT = 1 (열적) — 응답이 열요동에 묻힌다", C_MEAS)):
        a.set_xlim(-cold["ell"] * cold["n"] / 2 * 1.06, cold["ell"] * cold["n"] / 2 * 1.06)
        a.set_ylim(-ylim, ylim)
        a.axhline(0, c=C_GREY, lw=0.7, ls=":")
        a.axhspan(-dat["lk"], dat["lk"], color=C_GREY, alpha=0.18)
        ln, = a.plot([], [], "-", c=col, lw=1.4, zorder=2)
        pt, = a.plot([], [], "o", c=col, ms=6.5, zorder=3)
        tr, = a.plot([], [], "v", c=C_BAD, ms=11, zorder=4)
        a.set_ylabel("y  [d]"); a.set_title(ttl, fontsize=9.5, loc="left")
        a.grid(alpha=0.2)
        arts.append((ln, pt, tr, dat))
    axes[0].text(0.995, 0.05, "회색 띠 = ±ℓ_k (열요동)", transform=axes[0].transAxes,
                 ha="right", fontsize=8, c=C_GREY)
    axes[1].set_xlabel("x  [d]   (▼ = 구동 트랩 중심)")

    def upd(i):
        out = []
        for ln, pt, tr, dat in arts:
            q = dat["frames"][i]
            ln.set_data(q[:, 0], q[:, 1]); pt.set_data(q[:, 0], q[:, 1])
            trp = dat["trapped"]
            tr.set_data([q[trp[len(trp) // 2], 0]], [dat["drive"][i]])
            out += [ln, pt, tr]
        return out

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    anim = FuncAnimation(fig, upd, frames=nf, blit=True, interval=55)
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "chain_bend_motion.gif"
    anim.save(out, writer=PillowWriter(fps=18))
    plt.close(fig)
    return out, cold


def main() -> int:
    png = make_figure()
    print(f"그래프  → {png.relative_to(ROOT)}")
    if "--fig-only" in sys.argv:           # 라벨을 고칠 때 시뮬레이션을 다시 돌리지 않는다
        return 0
    gif, meta = make_animation()
    print(f"애니메이션 → {gif.relative_to(ROOT)}   "
          f"(De={meta['de']:.1f}, dt={meta['dt']:.2e}, 평형화 {meta['n_eq']:,} 스텝)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
