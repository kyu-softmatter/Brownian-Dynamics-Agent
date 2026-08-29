"""락인(lock-in) 응답함수 추정 — 진동 구동계의 복소 강성 `K*(ω)`.

**두 번 나와서 올렸습니다**: `scratch/verify_chain_bend_gates.py`(관문 A)가 이 추정량을
해석해에 대조해 검증했고, `cases/chain_bend_2d.py` 의 L4 생산 런이 같은 것을 씁니다.
함수 본문은 검증된 것을 **그대로** 옮겼습니다 — 여기서 손대면 관문 A의 검증이 무효가 됩니다.

## ★★ 공칭 진폭을 추정량에 쓰면 안 됩니다 (관문 A가 FAIL하면서 발견)

구동 유령을 `U` 스텝마다 옮기면 구동이 **영차 유지(zero-order hold)** 가 되어
기본파가 `sinc(ωΔt/2)` 배로 줄고 위상이 `ωΔt/2` 늦습니다 (`Δt = U·dt`).
실측(De=10, 관문 A): `|ŷ_c|/a = 0.98999` · 위상 `−0.2522 rad`
                     ZOH 예측 `0.99040` · `−0.2404 rad`

**유령 위치를 같이 재서 측정된 위상자 `ŷ_c` 를 쓰면 ZOH 감쇠가 분자·분모에서 상쇠됩니다** —
비드는 공칭 사인이 아니라 유령이 **실제로 있는 곳**에 반응하기 때문입니다.
공칭을 쓴 `K′` 은 De=10 에서 −6559 (부호까지 틀림, 오차 236%), 측정 위상자로는 5863 (21%).
**에러 없이 조용히 틀립니다.**

→ 그래서 `k_star()` 는 `drive_hat` 을 **인자로 요구**합니다. 기본값을 두지 않은 것이
  의도입니다: 공칭을 쓰려면 호출부가 명시적으로 그렇게 써야 합니다.
"""
from __future__ import annotations

import math

import numpy as np


def lockin_blocks(t, s, omega: float, *, harmonic: int = 1, n_blocks: int = 10):
    """블록마다 `ŝ = s_in + i·s_qu`.

    규약: `s(t) = Im[ŝ e^{iωt}]`, 구동 `y_c = a sin(ωt)` → `ŷ_c = a`.
    블록으로 쪼개는 이유는 SEM 을 그 산포에서 얻기 위함입니다 (`agg`).
    """
    t = np.asarray(t, dtype=float)
    s = np.asarray(s, dtype=float)
    ph = harmonic * omega * t
    out = []
    for bt, bs in zip(np.array_split(ph, n_blocks), np.array_split(s, n_blocks)):
        out.append(complex(2.0 * np.mean(bs * np.sin(bt)), 2.0 * np.mean(bs * np.cos(bt))))
    return np.array(out)


def k_star(y_hat: complex, drive_hat: complex, k_t: float, omega: float,
           gamma: float = 1.0, mass: float = 0.0) -> complex:
    """시료의 복소 강성. 궤적 두 개(비드·유령)만 있으면 됩니다 — 힘 로깅이 불필요합니다.

    `m ÿ + γ ẏ = −k_t(y − y_c) + F_sample` 를 ω 성분으로 쓰면
      `(iωγ − mω²) ŷ = −k_t(ŷ − ŷ_c) + F̂_sample`
    → `K* ≡ −F̂_sample/ŷ = k_t·ŷ_c/ŷ − k_t − iωγ + mω²`

    ★ `drive_hat` 은 **측정된** 유령 위상자여야 합니다 (모듈 도크스트링 참조).
      BD 는 `mass=0`. 잡음은 ω 성분에 코히런트하게 기여하지 않으므로 기댓값에서 정확합니다.
    """
    return k_t * drive_hat / y_hat - k_t - 1j * omega * gamma + mass * omega ** 2


def agg(vals) -> tuple[complex, float]:
    """블록 값들 → `(전체 추정, SEM)`. SEM 은 실·허 중 큰 쪽을 보수적으로 씁니다."""
    vals = np.asarray(vals)
    n = len(vals)
    sem = max(vals.real.std(ddof=1), vals.imag.std(ddof=1)) / math.sqrt(n)
    return complex(vals.mean()), float(sem)


def zoh_factor(omega: float, dt_update: float) -> complex:
    """영차 유지 보정 인자 `sinc(ωΔt/2)·e^{−iωΔt/2}` — **진단용**입니다.

    측정 위상자를 쓰면 보정이 필요 없습니다. 이 함수는 "구동을 정량적으로 이해했는가"를
    확인할 때(관문 A ①) 측정된 `ŷ_c/a` 를 이 예측과 대조하는 데 씁니다.
    """
    x = 0.5 * omega * dt_update
    s = 1.0 if x == 0 else math.sin(x) / x
    return s * complex(math.cos(x), -math.sin(x))


__all__ = ["lockin_blocks", "k_star", "agg", "zoh_factor"]
