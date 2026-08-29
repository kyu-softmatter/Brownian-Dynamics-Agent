"""One pint registry.

* Always use this module's `u`/`Q`. With different registries, pint refuses to
  combine two Quantities -- the first two cases had each created their own
  `pint.UnitRegistry()`.
"""
import pint

u = pint.UnitRegistry()
Q = u.Quantity

kB = Q(1.380649e-23, "J/K")      # CODATA 2019 (a defined value)

__all__ = ["u", "Q", "kB"]
