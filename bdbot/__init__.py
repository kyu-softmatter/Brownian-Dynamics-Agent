"""bdbot — Brownian dynamics 파이프라인의 공통 부분 (모듈 20개).

**Phase 1-C 산출물 + 이후 증분.** 여기 있는 것은 전부 최소 두 케이스에서
**실제로 두 번 나온 것**입니다. 판정 근거는 skill `bd-physics` §6.3 대조표입니다.
(`pairpot`·`traps`·`lockin`·`run`은 1-C 이후에 두 번째/세 번째 케이스가 나오면서
나중에 올라왔습니다 — 아래 각 모듈 docstring에 몇 번째로 올라왔는지 적혀 있습니다.)

    cli          ⭐️ 앞단 진입점 — status · intake · system · nondim · run
    intake       L0 Observation 스키마 + 검사 (스케치 5장 실사용에서 도출)
    interactions 콜로이드 상호작용 카탈로그 + 추천기 (U_ij 가 빈칸일 때)
    physical     L2 PhysicalSystem 로더 + tier·derived_from·유도값 재계산 검증
    units        단일 pint 레지스트리 (섞이면 pint가 거부한다)
    provenance   Provenanced — 모든 숫자에 출처 + tier
    materials    γ=3πηd, D_t=kT/γ, τ_B=d²/D_t, m, τ_p  (구 + 뉴턴 유체)
    pairpot      소프트 반발 페어의 수치 — U, U'', 최근접 접근거리 r_min (dt를 정하는 물리)
    scales       ScaleLedger — 길이·시간·에너지 원장 + 기준 + 근거 + 필수역할 완전성
    nondim       ⭐️ L3 NondimSpec — L2와 L4 사이의 유일한 계약. 원장·무차원수·검사·역변환
    checks       Check(모델/적분/기하/통계) + 하드/소프트 판정 + dt = 10⁻²·γ/(국소 강성)
    traps        조화 트랩 하나로 고정·등속·진동 구동을 표현 (세 번 나와서 올라옴)
    sim          2D 프레임 · Brownian 적분기 · GSD 라이터 · 시드 처리 · 최소이미지
    run          ⭐️ L5 조립 계약(`@builder`/`Build`) + L6 실행 루프 + L7 산출물 저장.
                 케이스는 `build(spec) -> Build` 만 주면 나머지(평형화/생산 루프,
                 가드 호출, metrics.json 저장)는 공통이다. **5개 케이스 전부** 이
                 계약을 쓴다 (2026-08-05 완료). `trap-2d-5um`·`soft-r3-2d-A-sweep`·
                 `abp-rod-2d-run-flip` 세 케이스는 이관하며 재실행해 옛 결과와
                 대조했고, 관측량이 표시 정밀도(소수 15자리)까지 완전히 일치했다 —
                 리팩터가 물리를 바꾸지 않았다는 뜻이다.
    health       L4 수치 건전성 판정 — Guard·judge·step_health (L3 되먹임)
    lockin       진동 구동계의 복소강성 K*(ω) 추정 (두 번 나와서 올라옴)
    report       DimensionlessReport 렌더러 (케이스별 블록은 케이스가 준다)
    runid        콘텐츠 주소 지정 run_id + 재실행 방지
    metrics      metrics.json 스키마 (postmortem의 유일한 입력)
    stats        블록 평균 · 자기상관 보정 · 불편 자기상관

`cli` 는 마스터플랜 §15.10의 방어책입니다 — Claude Code 세션 밖에서도 같은 명령으로
같은 결과가 나와야 하고, Phase 8-rest의 훅이 가로챌 대상도 이 명령들입니다.
`bdbot.cli` 는 무거운 의존성(hoomd·freud)을 임포트하지 않습니다 — 앞단만 쓸 때 빠릅니다.
`run`도 같은 이유로 hoomd를 함수 안에서만(`execute` 내부) 지연 임포트한다.

**의도적으로 넣지 않은 것** (한 번만 나왔거나 계마다 다름):
    평형 판정 지표 · 관측량 · 검증 전략 · 지배 시간척도 선택 · 초기 배치 · 표본 수집 루프
    → 케이스 스크립트에 남깁니다. 세 번째 케이스에서 또 나오면 그때 올립니다.

절대 규칙은 [CLAUDE.md](../CLAUDE.md), 물리 절차는 skill `bd-physics`,
HOOMD 함정은 skill `bd-hoomd` 를 보세요.
"""
from . import (checks, intake, interactions, materials, metrics, nondim, pairpot, physical,
               provenance, report, runid, scales, sim, stats)
from .checks import GATE, Check, bias_from_dt, dt_from_bias, dt_from_gate, relaxation_time, verdict
from .nondim import Group, NondimSpec, Reference
from .provenance import Provenanced, load_node
from .scales import Scale, ScaleLedger, thermal_reference
from .units import Q, kB, u

__version__ = "0.1.0"

__all__ = [
    "u", "Q", "kB",
    "Provenanced", "load_node",
    "ScaleLedger", "Scale", "thermal_reference",
    "NondimSpec", "Reference", "Group",
    "Check", "GATE", "verdict", "relaxation_time", "dt_from_gate", "dt_from_bias",
    "bias_from_dt",
    "units", "provenance", "materials", "scales", "nondim", "checks", "report", "runid",
    "metrics", "stats", "sim", "intake", "physical", "interactions", "pairpot",
]
