"""HOOMD truncates `Simulation(seed=...)` to 16 bits — a repo-wide gate.

★ Why this file exists. Measured 2026-09-16, while timing the sedimentation
campaign:

    *Warning*: Provided seed 20260916 is larger than 65535. Truncating to 10292.

The consequence is specific and silent. `run_id` is the content hash of the spec
and the spec stores the UNTRUNCATED seed, so two runs whose seeds differ by
exactly 65536 get **two different run_ids and one thermal-noise stream**. A
seed-to-seed standard deviation computed over such a pair is not an error bar —
it is two copies of one number, and `A2` would be satisfied by a value that
carries no information.

This is the same defect class the repository keeps recording: nothing fails, no
verdict changes, and the symptom is a number that is too *good*. HOOMD does emit
a warning, but `cases/sediment_3d.py` runs at `notice_level=0` and 296 archived
specs went by without anyone reading it.

⚠ **No archived run is affected** — that was audited, not assumed, and
`test_no_archived_replicate_group_collides` is that audit kept running. The point
of the gate is the next campaign, not this one.
"""
from __future__ import annotations

import collections
import json
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

#: HOOMD's own limit. Asserted against the library below rather than trusted.
SEED_BITS = 0xFFFF


def _specs():
    for p in sorted((ROOT / "specs").glob("*.json")):
        yield p, json.loads(p.read_text())


def test_hoomd_really_does_truncate_the_seed():
    """★ The load-bearing measurement. If a HOOMD upgrade widens the seed, this
    test fails and the whole file can go — which is the correct outcome, and is
    why the behaviour is asserted instead of being written in a comment."""
    hoomd = pytest.importorskip("hoomd")
    dev = hoomd.device.CPU(notice_level=0)
    for asked in (1, 65535, 65536, 20260916, 20260916 + 65536):
        sim = hoomd.Simulation(device=dev, seed=asked)
        assert sim.seed == asked & SEED_BITS, (asked, sim.seed)
    #  and the collision itself, stated as an executable fact
    a = hoomd.Simulation(device=dev, seed=20260916).seed
    b = hoomd.Simulation(device=dev, seed=20260916 + 65536).seed
    assert a == b == 10292


def test_no_archived_replicate_group_collides():
    """Group every archived spec by everything EXCEPT the seed; within a group,
    two DISTINCT seeds that truncate to the same value are the defect.

    ⚠ Two specs in DIFFERENT groups sharing a truncated seed is fine and is not
    flagged — a parameter sweep that reuses one noise stream across its points is
    correlated sampling, which is deliberate. Only replicates are wrong, and this
    is the distinction the first version of this audit got wrong: it reported 6
    'collisions' in `soft-r3-2d-A-sweep` that were sweep points, not replicates.
    """
    groups: dict[tuple, set[int]] = {}
    for _, d in _specs():
        num = dict(d.get("numerics", {}))
        seed = num.pop("seed", None)
        if not isinstance(seed, int):
            continue
        key = (d["case"], json.dumps(d.get("params", {}), sort_keys=True),
               json.dumps(num, sort_keys=True))
        groups.setdefault(key, set()).add(seed)

    assert groups, "no archived spec carries an integer seed -- the audit is vacuous"
    collisions = []
    for key, seeds in groups.items():
        byt = collections.defaultdict(list)
        for s in seeds:
            byt[s & SEED_BITS].append(s)
        for trunc, ss in byt.items():
            if len(ss) > 1:
                collisions.append(f"{key[0]}: {sorted(ss)} all become {trunc}")
    assert not collisions, "\n".join(collisions)


def test_the_audit_would_catch_a_collision_if_one_existed():
    """The audit above passes, so on its own it is indistinguishable from an
    audit that checks nothing. This runs it against a synthetic group that DOES
    collide."""
    seeds = {20260916, 20260916 + 65536, 5}
    byt = collections.defaultdict(list)
    for s in seeds:
        byt[s & SEED_BITS].append(s)
    bad = {t: ss for t, ss in byt.items() if len(ss) > 1}
    assert bad == {10292: [20260916, 20326452]} or \
           list(bad) == [10292] and sorted(bad[10292]) == [20260916, 20326452], bad


def test_at_least_one_archived_case_uses_seeds_in_the_dangerous_range():
    """Otherwise the audit is a gate over an empty set. 8-digit date seeds are
    this repository's habit, and they are all above 65535."""
    big = [s for _, d in _specs()
           if isinstance(s := d.get("numerics", {}).get("seed"), int) and s > SEED_BITS]
    assert len(big) >= 50, f"only {len(big)} archived specs use a seed > 65535"


# ── the sedimentation campaign's sealed seed list ──────────────────────────

def test_the_sealed_sediment_seeds_survive_truncation():
    doc = yaml.safe_load(
        (ROOT / "campaigns/sediment_preregistration/prediction.yaml").read_text())
    entries = doc["runs"]["entries"]
    per_arm = collections.defaultdict(list)
    for e in entries:
        per_arm[(e["arm"], e["lxy"], e.get("init", doc["runs"]["common"]["init"]))
                ].append(int(e["seed"]))
    for key, seeds in per_arm.items():
        trunc = [s & SEED_BITS for s in seeds]
        assert len(set(trunc)) == len(set(seeds)), (
            f"{key}: seeds {seeds} collapse to {sorted(set(trunc))} inside HOOMD")
    #  and the two COMPETING arms must each have 3 independent streams, with no
    #  seed shared between them -- the verdict is a comparison between the arms,
    #  so a shared RNG stream would correlate them in a way the seed-to-seed sd
    #  could not see
    a = per_arm[("cs_eff", 12, "cs_eff")]
    b = per_arm[("cs_nominal", 12, "cs")]
    assert len(a) == 3 and len({s & SEED_BITS for s in a}) == 3, a
    assert len(b) == 3 and len({s & SEED_BITS for s in b}) == 3, b
    assert not (set(a) & set(b)), f"arms A and B share seeds {set(a) & set(b)}"
    assert not ({s & SEED_BITS for s in a} & {s & SEED_BITS for s in b}), (
        "arms A and B share a seed INSIDE hoomd even though their spec seeds differ")
    #  no production seed may be 1 -- that is the disclosed pilots' seed
    assert all(s != 1 for seeds in per_arm.values() for s in seeds)


def test_the_campaign_driver_refuses_a_colliding_seed_list():
    """★ Mutation test on the gate that actually blocks `--prepare`."""
    import campaigns.s31_sedimentation as S31
    assert S31.SEED_BITS == SEED_BITS
    assert S31.seed_collisions() == [], S31.seed_collisions()

    real_runs = S31.runs
    try:
        S31.runs = lambda: [
            {"id": "P1", "arm": "primary", "lxy": 12.0, "seed": 20260916,
             "N": 584, "tau_sed": 4.0, "init": "cs"},
            {"id": "P2", "arm": "primary", "lxy": 12.0, "seed": 20260916 + 65536,
             "N": 584, "tau_sed": 4.0, "init": "cs"},
        ]
        bad = S31.seed_collisions()
        assert len(bad) == 1 and "10292" in bad[0], bad
        #  ⚠ and the same two seeds in DIFFERENT arms must NOT be flagged: the
        #     arms are never pooled into one sd, so sharing a stream there is
        #     only correlated sampling
        S31.runs = lambda: [
            {"id": "P1", "arm": "primary", "lxy": 12.0, "seed": 20260916,
             "N": 584, "tau_sed": 4.0, "init": "cs"},
            {"id": "C2", "arm": "convergence", "lxy": 12.0,
             "seed": 20260916 + 65536, "N": 584, "tau_sed": 8.0, "init": "ideal"},
        ]
        assert S31.seed_collisions() == []
    finally:
        S31.runs = real_runs
    assert S31.seed_collisions() == []
