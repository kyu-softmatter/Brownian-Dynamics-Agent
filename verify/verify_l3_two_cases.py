"""L3 적대적 검사 — `trap-drag-2d-hex300` · `chain-bend-2d-oscill` (2026-08-04).

`verify_nondim_guards.py` 가 `NondimSpec` **자체**를 시험한다면, 여기서는 새로 만든
**두 케이스의 원장·무차원수·해시**를 일부러 망가뜨려 시험합니다. "조용히 통과"와
"검사를 안 함"은 다릅니다 (CLAUDE.md 작업 관행).

특히 이 두 케이스에만 있는 것:
  · chain-bend 는 주기경계가 없어 `box` 역할을 `declare_absent` 로 비웁니다 —
    비우는 것을 **잊으면** 필수역할 누락으로 잡혀야 합니다.
  · chain-bend 의 `λ_max` 는 굽힘·신축 두 블록의 **큰 쪽**입니다. 더하면(분리된
    자유도를 섞으면) dt 가 18% 과소평가되어 비용만 늘어납니다.
  · trap-drag 의 물리계가 run_id 를 덮는가 — 1-B에서 실제로 뚫렸던 구멍입니다.

    $PY scratch/verify_l3_two_cases.py
"""
import sys, json, copy, argparse, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))
import chain_bend_2d as CB, trap_drag_2d as TD
from bdbot import nondim as ND

A = argparse.Namespace(dt_scale=1.0, cycles=CB.N_CYCLES, samples=2000, spec=False, report=True)
ok = lambda b: "✓" if b else "✗ 잡지 못함"
res = []

# ── ① chain-bend: box 를 declare_absent 하지 않으면 필수역할 누락으로 잡히는가
s = CB.load_system(ROOT/"intake/chain-bend-2d-oscill/system.yaml")
lg, spec, g, c, *_ = CB.build_spec(s, 85.0, A)
lg.absent.pop("box")
errs = [i for i in spec.validate() if i.level == "error"]
res.append(("chain: box 미선언 → 필수역할 누락", any("ledger.box" in i.where for i in errs)))
lg.declare_absent("box", "복원")

# ── ② chain-bend: dt_star 를 원장과 어긋나게 하면 잡는가 (HOOMD가 딴 스텝으로 돈다)
spec.numerics["dt_star"] *= 1.0000001
res.append(("chain: dt_star 를 1e-7 어긋냄", any("numerics.dt_star" in i.where
            for i in spec.validate() if i.level == "error")))
spec.numerics["dt_star"] /= 1.0000001
res.append(("chain: 되돌리면 통과", not spec.errors))

# ── ③ chain-bend: 무차원수 값만 살짝 틀면 원장 재계산이 잡는가
De = next(x for x in spec.groups if x.name == "De"); orig = De.value
De.value = orig * 1.000001
res.append(("chain: De 를 1e-6 어긋냄", any("groups.De" in i.where
            for i in spec.validate() if i.level == "error")))
De.value = orig

# ── ④ trap-drag: 물리계를 바꾸면 run_id 가 바뀌는가 (1-B에서 실제로 뚫렸던 구멍)
t = TD.load_system(ROOT/"intake/trap-drag-2d-hex300/system.yaml")
Ta = argparse.Namespace(dt_scale=1.0, traverse=1.0, samples=2000)
lg2 = TD.build_ledger(t, dt_scale=1.0, n_traverse=1.0)
g2, c2, Gam, ex = TD.analyze_scales(t, lg2)
mk = lambda raw: ND.NondimSpec(case=t["label"], system=raw, reference=lg2.ref, ledger=lg2,
        groups=g2, checks=c2, params={"A": t["A"]},
        numerics={"dt_star": lg2.ratio("times","dt","tau_B"), "n_prod": 1}, nhex=12).run_id()
base = mk(t["_raw"])
alt = copy.deepcopy(t["_raw"]); alt["particle"]["diameter"]["value"] = 0.5   # 5µm → 0.5µm
res.append(("trap-drag: d 를 10배 바꾸면 run_id 변화", base != mk(alt)))
doc = copy.deepcopy(t["_raw"]); doc["description"] = "주석만 고침"
doc["particle"]["diameter"]["source"] = "출처 문구만 고침"
res.append(("trap-drag: 주석·출처만 고치면 run_id 불변", base == mk(doc)))

# ── ⑤ 저장된 스펙을 손으로 고치면 해시 검증이 잡는가
# ★ run_id 를 하드코딩하지 않는다 — 물리계가 바뀌면(N 300→306 정합 육방) 바뀌는 값이다.
p = max(ROOT.glob("specs/trap-drag-2d-hex300__*.json"), key=lambda f: f.stat().st_mtime)
raw = json.loads(p.read_text()); raw["params"]["A"] = 42.0
tmp = Path(tempfile.mkdtemp()) / "tampered.json"
tmp.write_text(json.dumps(raw))
res.append(("trap-drag: 스펙 손수정 → 해시 불일치", not ND.load(tmp).verify_hash()[0]))
res.append(("trap-drag: 원본은 해시 일치", ND.load(p).verify_hash()[0]))

# ── ⑥ trap-drag: 정합 육방 가드 — 격자가 주기박스와 안 맞으면 멈추는가
#     ⚠️ 처음에 φ 를 격자에서 재계산해 대조했는데 **항등적으로 통과**했다 (L_x·L_y 를 φ 에서
#        유도하므로). 돌지 않는 검사였다. 의미 있는 대조는 **L2가 적어둔 박스**와의 대조다.
def broken(mut):
    tt = copy.deepcopy(t); mut(tt)
    try:
        TD.build_ledger(tt); return False
    except ValueError:
        return True

from bdbot.units import Q as _Q
res.append(("trap-drag: n_y 홀수 → 거부", broken(lambda x: x.update(n_y=17, N=17*17))))
res.append(("trap-drag: N ≠ n_x·n_y → 거부", broken(lambda x: x.update(N=300))))
res.append(("trap-drag: φ 를 바꿔 박스와 불일치 → 거부", broken(lambda x: x.update(phi=0.30))))
res.append(("trap-drag: YAML box_length_x 만 손댐 → 거부",
            broken(lambda x: setattr(x["box_x"], "value", _Q(129.73, "um")))))
_lg = TD.build_ledger(t); _D = _lg.derived
res.append(("trap-drag: L_x/a_NN 가 정수 n_x",
            abs(_D["Lx_star"]/_D["a_nn_star"] - _D["n_x"]) < 1e-12))
res.append(("trap-drag: L_y/(√3/2 a_NN) 가 짝수 n_y",
            abs(_D["Ly_star"]/(3**0.5/2*_D["a_nn_star"]) - _D["n_y"]) < 1e-12 and _D["n_y"] % 2 == 0))

# ── ⑦ chain-bend: λ_max 는 두 블록의 큰 쪽인가 (더하면 dt 가 과도하게 작아진다)
D = lg.derived
res.append(("chain: λ_max = max(굽힘, 신축)", abs(D["lam_max"] - max(D["lam_bend"], D["lam_bond"])) < 1e-18))
res.append(("chain: 굽힘이 신축보다 빠르다", D["lam_bend"] > D["lam_bond"]))

print("="*72); print("적대적 검사 — 새 두 케이스의 L3 가드"); print("="*72)
for name, good in res: print(f"  {ok(good):<10} {name}")
n = sum(1 for _, g_ in res if g_)
print("="*72); print(f"{'✓ PASS' if n==len(res) else '✗ FAIL'}  {n}/{len(res)}"); print("="*72)
sys.exit(0 if n == len(res) else 1)
