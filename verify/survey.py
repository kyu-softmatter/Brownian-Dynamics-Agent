"""HOOMD capability survey — what physics modules can we actually build?

Introspects the installed HOOMD to enumerate every class we might map a
PhysicsModule onto. Output feeds the module registry design (master plan §5.6, §11).
"""
import inspect
import hoomd
import hoomd.md as md


def classes(mod, skip_private=True):
    out = []
    for name, obj in sorted(vars(mod).items()):
        if skip_private and name.startswith("_"):
            continue
        if inspect.isclass(obj):
            out.append(name)
    return out


def submodules(mod):
    out = []
    for name, obj in sorted(vars(mod).items()):
        if name.startswith("_"):
            continue
        if inspect.ismodule(obj):
            out.append(name)
    return out


def section(title):
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


print(f"hoomd {hoomd.version.version}   gpu={hoomd.version.gpu_enabled}   mpi={hoomd.version.mpi_enabled}")
print(f"compile flags: {hoomd.version.compile_flags}")

section("hoomd top-level")
print("  submodules:", ", ".join(submodules(hoomd)))
print("  classes   :", ", ".join(classes(hoomd)))

section("hoomd.md submodules")
print("  ", ", ".join(submodules(md)))

for sub in ["pair", "bond", "angle", "dihedral", "improper", "special_pair",
            "external", "force", "methods", "update", "constrain", "many_body",
            "nlist", "minimize", "compute", "tune", "long_range", "mesh"]:
    try:
        m = getattr(md, sub)
    except AttributeError:
        print(f"\n[md.{sub}] -- ABSENT --")
        continue
    section(f"md.{sub}")
    subs = submodules(m)
    if subs:
        print("  submodules:", ", ".join(subs))
    print("  classes:", ", ".join(classes(m)) or "(none)")
    # one level deeper for pair.aniso, methods.thermostats, external.field/wall
    for s in subs:
        inner = getattr(m, s)
        print(f"    md.{sub}.{s}: {', '.join(classes(inner)) or '(none)'}")

section("hoomd.variant / trigger / filter / write / custom / util")
for name in ["variant", "trigger", "filter", "write", "custom", "util", "tune",
             "logging", "state", "device", "box", "data", "error", "wall", "mesh",
             "communicator", "operation", "simulation", "update", "hpmc", "mpcd"]:
    try:
        m = getattr(hoomd, name)
    except AttributeError:
        continue
    if inspect.ismodule(m):
        subs = submodules(m)
        line = ", ".join(classes(m)) or "(none)"
        print(f"\n  hoomd.{name}: {line}")
        for s in subs:
            inner = getattr(m, s)
            print(f"    .{s}: {', '.join(classes(inner)) or '(none)'}")

section("KEY SIGNATURES (모듈 설계에 직접 필요한 것들)")
targets = [
    ("md.methods.Brownian", md.methods.Brownian),
    ("md.methods.Langevin", getattr(md.methods, "Langevin", None)),
    ("md.methods.OverdampedViscous", getattr(md.methods, "OverdampedViscous", None)),
    ("md.force.Active", md.force.Active),
    ("md.force.Custom", getattr(md.force, "Custom", None)),
    ("md.update.ActiveRotationalDiffusion", getattr(md.update, "ActiveRotationalDiffusion", None)),
    ("md.pair.LJ", md.pair.LJ),
    ("md.pair.Table", getattr(md.pair, "Table", None)),
    ("md.bond.Harmonic", getattr(md.bond, "Harmonic", None)),
    ("md.bond.FENEWCA", getattr(md.bond, "FENEWCA", None)),
    ("md.angle.Harmonic", getattr(md.angle, "Harmonic", None)),
    ("md.constrain.Rigid", getattr(md.constrain, "Rigid", None)),
    ("md.minimize.FIRE", getattr(md.minimize, "FIRE", None)),
    ("hoomd.write.GSD", hoomd.write.GSD),
    ("hoomd.write.Burst", getattr(hoomd.write, "Burst", None)),
    ("hoomd.write.HDF5Log", getattr(hoomd.write, "HDF5Log", None)),
    ("hoomd.variant.Ramp", getattr(hoomd.variant, "Ramp", None)),
    ("hoomd.variant.Cycle", getattr(hoomd.variant, "Cycle", None)),
    ("hoomd.custom.Action", getattr(hoomd.custom, "Action", None)),
    ("hoomd.Simulation", hoomd.Simulation),
]
for label, cls in targets:
    if cls is None:
        print(f"\n  {label}\n      -- ABSENT --")
        continue
    try:
        sig = inspect.signature(cls.__init__)
        params = str(sig).replace("(self, ", "(").replace("(self)", "()")
    except (TypeError, ValueError):
        params = "(introspection failed)"
    print(f"\n  {label}{params}")

section("ANISOTROPIC / SHAPE 관련 (타원체·막대 지원 여부)")
for path in ["md.pair.aniso", "md.constrain"]:
    obj = md
    for part in path.split(".")[1:]:
        obj = getattr(obj, part, None)
        if obj is None:
            break
    print(f"  {path}: {', '.join(classes(obj)) if obj else '-- ABSENT --'}")

section("VARIANT 종류 (시간 의존 구동)")
print("  hoomd.variant:", ", ".join(classes(hoomd.variant)))
for s in submodules(hoomd.variant):
    print(f"    .{s}:", ", ".join(classes(getattr(hoomd.variant, s))))
