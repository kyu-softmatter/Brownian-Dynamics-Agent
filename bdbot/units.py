"""단일 pint 레지스트리.

★ 반드시 이 모듈의 `u`/`Q`를 쓰세요. 레지스트리가 다르면 pint가 두 Quantity를
  섞지 못합니다 (1-A·1-B가 각자 `pint.UnitRegistry()`를 만들고 있었습니다).
"""
import pint

u = pint.UnitRegistry()
Q = u.Quantity

kB = Q(1.380649e-23, "J/K")      # CODATA 2019 (정의값)

__all__ = ["u", "Q", "kB"]
