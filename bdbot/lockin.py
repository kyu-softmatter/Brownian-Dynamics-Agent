"""Lock-in response-function estimation -- the complex stiffness `K*(omega)` of
an oscillatory-driven system.

**Promoted after appearing twice**: `verify/verify_chain_bend_gates.py` (gate A)
verified this estimator against an analytic solution, and the L4 production run
in `cases/chain_bend_2d.py` uses the same thing. The function bodies were carried
over **unchanged** from the verified version -- touching them here invalidates
gate A's verification.

## ** Do not put the nominal amplitude into the estimator (found when gate A FAILed)

Moving the driving ghost every `U` steps makes the drive a **zero-order hold**:
the fundamental shrinks by `sinc(omega*dt/2)` and the phase lags by
`omega*dt/2` (`dt = U*dt_step`).
Measured (De=10, gate A): `|y_hat_c|/a = 0.98999` · phase `-0.2522 rad`
                          ZOH prediction `0.99040` · `-0.2404 rad`

**Measure the ghost position as well and use the measured phasor `y_hat_c`: the
ZOH attenuation then cancels between numerator and denominator** -- because the
bead responds to where the ghost **actually is**, not to the nominal sine.
Using the nominal value, `K'` at De=10 came out -6559 (wrong even in sign, 236%
error); with the measured phasor, 5863 (21%).
**It fails silently, with no error.**

-> So `k_star()` **requires** `drive_hat` as an argument. Having no default is
   deliberate: using the nominal value means the caller has to say so explicitly.
"""
from __future__ import annotations

import math

import numpy as np


def lockin_blocks(t, s, omega: float, *, harmonic: int = 1, n_blocks: int = 10):
    """`s_hat = s_in + i*s_qu`, per block.

    Convention: `s(t) = Im[s_hat e^{i*omega*t}]`, and for a drive
    `y_c = a sin(omega*t)` this gives `y_hat_c = a`.
    The reason for splitting into blocks is to obtain the SEM from their spread
    (`agg`).
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
    """The sample's complex stiffness. Only two trajectories are needed (bead and
    ghost) -- no force logging.

    Writing `m*y'' + gamma*y' = -k_t(y - y_c) + F_sample` in its omega component
    gives
      `(iωγ − mω²) ŷ = −k_t(ŷ − ŷ_c) + F̂_sample`
    → `K* ≡ −F̂_sample/ŷ = k_t·ŷ_c/ŷ − k_t − iωγ + mω²`

    * `drive_hat` must be the **measured** ghost phasor (see the module
      docstring). BD has `mass=0`. Noise does not contribute coherently to the
      omega component, so this is exact in expectation.
    """
    return k_t * drive_hat / y_hat - k_t - 1j * omega * gamma + mass * omega ** 2


def agg(vals) -> tuple[complex, float]:
    """Block values -> `(overall estimate, SEM)`.

    The SEM conservatively takes the larger of the real and imaginary parts.
    """
    vals = np.asarray(vals)
    n = len(vals)
    sem = max(vals.real.std(ddof=1), vals.imag.std(ddof=1)) / math.sqrt(n)
    return complex(vals.mean()), float(sem)


def zoh_factor(omega: float, dt_update: float) -> complex:
    """The zero-order-hold correction factor `sinc(w*dt/2)*e^{-i*w*dt/2}` --
    **for diagnostics.**

    No correction is needed when the measured phasor is used. This function exists
    to check "have we understood the drive quantitatively" (gate A step 1) by
    comparing the measured `y_hat_c/a` against this prediction.
    """
    x = 0.5 * omega * dt_update
    s = 1.0 if x == 0 else math.sin(x) / x
    return s * complex(math.cos(x), -math.sin(x))


__all__ = ["lockin_blocks", "k_star", "agg", "zoh_factor"]
