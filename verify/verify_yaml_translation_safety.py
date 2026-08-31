#!/usr/bin/env python
"""Prove that translating a YAML file changed ONLY string values.

The Python token checker cannot help here, so this is the YAML analogue: parse
both files, walk them in lockstep, and require that

  - the key sets match at every level, in the same order
  - every non-string leaf (int, float, bool, None) is byte-identical
  - only `str` leaves are allowed to differ

A key rename, a reordered list, a changed number or a changed type is a
structural change and is reported. That is exactly the class of mistake a
translation must not make -- and unlike prose, it is silently load-bearing here:
these files feed spec hashes and gate declarations.

    $PY yamlsafe.py before.yaml after.yaml
"""
import sys

import yaml


def walk(a, b, path=""):
    bad = []
    if type(a) is not type(b):
        # int/float are interchangeable only when equal in value
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            if a != b:
                bad.append(f"{path}: number {a!r} -> {b!r}")
            return bad
        bad.append(f"{path}: type {type(a).__name__} -> {type(b).__name__}")
        return bad
    if isinstance(a, dict):
        if list(a) != list(b):
            only_a = [k for k in a if k not in b]
            only_b = [k for k in b if k not in a]
            if only_a or only_b:
                bad.append(f"{path}: keys changed -a={only_a} +b={only_b}")
            else:
                bad.append(f"{path}: key ORDER changed")
            return bad
        for k in a:
            bad += walk(a[k], b[k], f"{path}.{k}" if path else str(k))
    elif isinstance(a, list):
        if len(a) != len(b):
            bad.append(f"{path}: list length {len(a)} -> {len(b)}")
            return bad
        for i, (x, y) in enumerate(zip(a, b)):
            bad += walk(x, y, f"{path}[{i}]")
    elif isinstance(a, str):
        pass                                     # the only thing allowed to change
    elif a != b:
        bad.append(f"{path}: {a!r} -> {b!r}")
    return bad


def main(p1: str, p2: str) -> int:
    a = yaml.safe_load(open(p1, encoding="utf-8"))
    b = yaml.safe_load(open(p2, encoding="utf-8"))
    bad = walk(a, b)
    ko = sum(1 for line in open(p2, encoding="utf-8")
             if any("가" <= c <= "힣" for c in line))
    if bad:
        print(f"CHANGED {p2} — {len(bad)} structural difference(s):")
        for x in bad[:20]:
            print("  " + x)
        return 1
    print(f"OK   {p2}: structure identical, only string values differ, "
          f"Hangul lines left {ko}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1], sys.argv[2]))
