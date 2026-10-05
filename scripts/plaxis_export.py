"""Build the reconstructed geometry in PLAXIS 3D Input through its Python
API (plxscripting): elastic plate and beam materials, a planar surface with
a plate for every polygon of the exchange file, a line with a beam for
every bar piece.

Usage (PLAXIS Input open, remote scripting server enabled in
Expert > Configure remote scripting server):

    python plaxis_export.py model.plaxis.json --port 10000 --password PASS
    python plaxis_export.py model.plaxis.json --dry-run      # print commands

Run it with the Python distribution shipped with PLAXIS (it has
plxscripting), or any Python 3 with `pip install plxscripting`.

Property names of material data sets changed between PLAXIS versions
(Bentley, "Material property changes for Python scripting"): Plates 3D
use Identification/d/Isotropic/StructNu12/Gamma from V22.02 on, D3d in
V22.00, MaterialName/d/IsIsotropic/Nu12/w in V21; beams Identification/
CrossSectionType "User-defined"/A/I2/I3/E/Gamma (V21: MaterialName/BeamType/
Izz/Iyy/w). Each material is tried with these sets in turn, and the set that
worked is reported. Materials are assigned with `setmaterial`, as PLAXIS
logs it. A beam's local axis 2 (PLAXIS: section height direction) is set to
the LIRA Z1 axis for rectangular sections. On an error every object created
by this run is deleted again, so a retry does not duplicate them.

Loads (kN, m): when the exchange file has them, every point, line and surface
load is created (`pointload`, `lineload`, `surfload`: Fx..Mz, qx/qy/qz_start/
_end, sigx/sigy/sigz) and, in Staged construction, a base phase activates the
whole structure and one phase per LIRA load case (named after it) activates
that case's loads. A load PLAXIS refuses is counted and reported, not fatal.
`--no-phases` creates the loads without phases.
"""
import argparse
import json
import sys
import time

# Plate material property sets, newest first.
PLATE_SETS = [
    dict(name="Identification", d="d", nu="StructNu12", iso=("Isotropic", True), gamma="Gamma",
         extra=("MaterialType", "Elastic")),                                  # V22.02 .. 2024.3
    dict(name="Identification", d="D3d", nu="StructNu12", iso=("Isotropic", True), gamma="Gamma",
         extra=("MaterialType", "Elastic")),                                  # V22.00
    dict(name="MaterialName", d="d", nu="Nu12", iso=("IsIsotropic", True), gamma="w",
         extra=None),                                                         # V21
]
BEAM_SETS = [
    dict(name="Identification", gamma="Gamma", i2="I2", i3="I3",
         extra=("MaterialType", "Elastic", "CrossSectionType", "User-defined")),  # V22.00 ..
    dict(name="MaterialName", gamma="w", i2="Izz", i3="Iyy",
         extra=("BeamType", "User-defined")),                                     # V21
]


def plate_properties(m, names):
    props = [names["name"], m["name"]]
    if names["extra"]:
        props += list(names["extra"])
    props += [names["d"], m["d"], names["iso"][0], names["iso"][1], "E1", m["e"], names["nu"], m["nu"]]
    if m["gamma"] > 0:
        props += [names["gamma"], m["gamma"]]   # unit weight (kN/m3) in every version
    return props


def beam_properties(m, names):
    props = [names["name"], m["name"]]
    if names["extra"]:
        props += list(names["extra"])
    props += ["E", m["e"], "A", m["a"], names["i2"], m["i2"], names["i3"], m["i3"]]
    if m["gamma"] > 0:
        props += [names["gamma"], m["gamma"]]
    return props


def last(result):
    """The object a command created: commands creating points with a line
    or surface return a list with the main object last."""
    return result[-1] if isinstance(result, (list, tuple)) else result


def make_material(g_i, create, sets, properties, m, log):
    """Create a material with the first property set the running PLAXIS
    version accepts; returns it and the index of that set."""
    for index, names in enumerate(sets):
        material = create()
        try:
            material.setproperties(*properties(m, names))
            return material, index
        except Exception as error:  # a property name unknown to this version
            log.append(f"{m['name']} (set {index}): {type(error).__name__}: {error}")
            try:
                g_i.delete(material)
            except Exception:
                pass
    raise RuntimeError(f"no known property set accepted for material {m['name']}: {log[-1]}")


def assign(g_i, obj, material):
    """`setmaterial` as PLAXIS logs it; the Material property otherwise."""
    try:
        g_i.setmaterial(obj, material)
    except Exception:
        obj.Material = material


def orient(line, axis2):
    """Local axis 2 of a line (Axis function Manual); False if refused."""
    try:
        line.AxisFunction = "Manual"
        line.Axis2x, line.Axis2y, line.Axis2z = axis2
        return True
    except Exception:
        return False


def build(g_i, data, shift=(0.0, 0.0, 0.0), progress=print, phases=True):
    """Create materials, plates and beams; returns a report. On an error the
    objects created so far are deleted and the error is raised again."""
    def point(p):
        return tuple(round(p[k] - shift[k], 9) for k in range(3))

    created = []

    def keep(obj):
        created.append(obj)
        return obj

    try:
        return _build(g_i, data, point, keep, progress, phases)
    except Exception:
        for obj in reversed(created):
            try:
                g_i.delete(obj)
            except Exception:
                pass
        raise


def _build(g_i, data, point, keep, progress, phases):
    g_i.gotostructures()
    log = []
    plate_mats, beam_mats, used = {}, {}, set()
    for m in data["plate_materials"]:
        plate_mats[m["name"]], index = make_material(g_i, g_i.platemat, PLATE_SETS, plate_properties, m, log)
        keep(plate_mats[m["name"]])
        used.add(f"plate set {index}")
    rectangular = set()
    for m in data["beam_materials"]:
        beam_mats[m["name"]], index = make_material(g_i, g_i.beammat, BEAM_SETS, beam_properties, m, log)
        keep(beam_mats[m["name"]])
        used.add(f"beam set {index}")
        if abs(m["width"] - m["height"]) > 1e-9:
            rectangular.add(m["name"])
    plates = beams = oriented = 0
    orientation_refused = False
    structure = []      # plates and beams, activated in the base phase
    polygons = {}       # outline of a plate polygon -> its surface
    started = time.time()
    total = sum(len(p["polygons"]) for p in data["plates"]) + len(data["beams"])
    done = 0
    for p in data["plates"]:
        for polygon in p["polygons"]:
            points = [point(x) for x in polygon]
            surface = keep(last(g_i.surface(*points)))
            polygons[outline(points)] = surface
            plate = last(g_i.plate(surface))
            structure.append(plate)
            if p["material"]:
                assign(g_i, plate, plate_mats[p["material"]])
            plates += 1
            done += 1
            if done % 200 == 0:
                progress(f"{done}/{total} objects, {time.time() - started:.0f} s")
    for b in data["beams"]:
        line = keep(last(g_i.line(point(b["start"]), point(b["end"]))))
        beam = last(g_i.beam(line))
        structure.append(beam)
        if b["material"]:
            assign(g_i, beam, beam_mats[b["material"]])
            if b["material"] in rectangular and not orientation_refused:
                if orient(line, b["axis2"]):
                    oriented += 1
                else:
                    orientation_refused = True
        beams += 1
        done += 1
        if done % 200 == 0:
            progress(f"{done}/{total} objects, {time.time() - started:.0f} s")
    rect_beams = sum(1 for b in data["beams"] if b["material"] in rectangular)
    warnings = list(data.get("warnings", []))
    if orientation_refused:
        warnings.append(f"PLAXIS refused the line local axis (AxisFunction/Axis2): {rect_beams - oriented} "
                        "rectangular beams keep the automatic axes - check their orientation")
    load_report = create_loads(g_i, data, point, keep, progress, phases, structure, polygons, warnings)
    return dict(plates=plates, beams=beams, plate_materials=len(plate_mats), loads=load_report,
                beam_materials=len(beam_mats), property_sets=sorted(used),
                rejected_property_sets=log, oriented_beams=oriented, rectangular_beams=rect_beams,
                seconds=round(time.time() - started, 1),
                missing_materials=data.get("missing_materials", []), warnings=warnings)


def outline(points):
    """A polygon outline independent of where it starts and of its direction."""
    keyed = [tuple(round(c, 6) for c in p) for p in points]
    return frozenset(keyed)


def create_loads(g_i, data, point, keep, progress, phases, structure, polygons, warnings):
    """Loads of the exchange file, and the phases that activate them. Every
    failure is counted and reported (`warnings`), none stops the export."""
    loads = data.get("loads") or []
    if not loads:
        return None
    made = {"point": 0, "line": 0, "surface": 0}
    refused = {}
    first_error = {}
    by_case = {}

    def refuse(kind, error):
        refused[kind] = refused.get(kind, 0) + 1
        first_error.setdefault(kind, f"{type(error).__name__}: {error}")

    def configure(obj, values):
        for name, value in values:
            setattr(obj, name, value)

    started = time.time()
    for n, load in enumerate(loads):
        kind = load["kind"]
        try:
            if kind == "point":
                obj = last(g_i.pointload(point(load["at"])))
                f, m = load["force"], load["moment"]
                configure(obj, [("Fx", f[0]), ("Fy", f[1]), ("Fz", f[2]), ("Mx", m[0]), ("My", m[1]), ("Mz", m[2])])
            elif kind == "line":
                obj = last(g_i.lineload(point(load["start"]), point(load["end"])))
                a, b = load["q_start"], load["q_end"]
                values = [("Distribution_x", "Linear"), ("Distribution_y", "Linear"), ("Distribution_z", "Linear")]
                configure(obj, [v for v in values if a != b])
                configure(obj, [("qx_start", a[0]), ("qy_start", a[1]), ("qz_start", a[2])])
                if a != b:
                    configure(obj, [("qx_end", b[0]), ("qy_end", b[1]), ("qz_end", b[2])])
            elif kind == "surface":
                for polygon in load["polygons"]:
                    points = [point(x) for x in polygon]
                    existing = polygons.get(outline(points))
                    try:
                        obj = last(g_i.surfload(*([existing] if existing is not None else points)))
                    except Exception:
                        if existing is not None:
                            raise
                        obj = last(g_i.surfload(keep(last(g_i.surface(*points)))))
                    s = load["sigma"]
                    configure(obj, [("sigx", s[0]), ("sigy", s[1]), ("sigz", s[2])])
                    keep(obj)
                    by_case.setdefault(load["case"], []).append(obj)
                    made["surface"] += 1
                continue
            else:
                continue
            keep(obj)
            by_case.setdefault(load["case"], []).append(obj)
            made[kind] += 1
        except Exception as error:   # a command or property this PLAXIS version does not know
            refuse(kind, error)
            if sum(refused.values()) >= 20 and not sum(made.values()):
                break                # nothing works: do not try every load
        if n % 200 == 199:
            progress(f"loads {n + 1}/{len(loads)}, {time.time() - started:.0f} s")
    for kind, count in refused.items():
        warnings.append(f"PLAXIS refused {count} {kind} loads (first error: {first_error[kind]})")
    report = dict(created=made, refused=refused, phases=[])
    if phases and sum(made.values()):
        try:
            report["phases"] = create_phases(g_i, data, by_case, structure, progress)
        except Exception as error:
            warnings.append(f"phases were not created ({type(error).__name__}: {error}); "
                            "the loads exist but are not activated in any phase")
    return report


def create_phases(g_i, data, by_case, structure, progress):
    """Staged construction: a base phase with the structure, a phase per load
    case derived from it with that case's loads."""
    g_i.gotostages()
    base = last(g_i.phase(g_i.InitialPhase))
    try:
        base.Identification = "Structure"
    except Exception:
        pass
    for obj in structure:
        g_i.activate(obj, base)
    names = {c[0]: c[1] for c in data.get("load_cases", [])}
    created = []
    for case in sorted(by_case):
        phase = last(g_i.phase(base))
        label = f"{case} {names.get(case, '')}".strip()
        try:
            phase.Identification = label
        except Exception:
            pass
        for obj in by_case[case]:
            g_i.activate(obj, phase)
        created.append(label)
        progress(f"phase {label}: {len(by_case[case])} loads")
    return created


class Recorder:
    """A stand-in for PLAXIS Input (`--dry-run`, tests): records every
    command as PLAXIS logs it and returns named objects."""

    def __init__(self, reject=(), fail_on=None, reject_commands=()):
        self.reject_commands = set(reject_commands)
        self.commands = []
        self.counts = {}
        self.reject = set(reject)
        self.fail_on = fail_on

    def _object(self, kind):
        self.counts[kind] = self.counts.get(kind, 0) + 1
        return _Object(self, f"{kind}_{self.counts[kind]}")

    def _log(self, name, args):
        def text(a):
            if isinstance(a, tuple):
                return "(" + " ".join(repr(x) for x in a) + ")"
            if isinstance(a, str):
                return f'"{a}"'
            return str(a)
        self.commands.append(" ".join([name] + [text(a) for a in args]))

    def __getattr__(self, name):
        kinds = {"platemat": "PlateMat", "beammat": "BeamMat", "plate": "Plate", "beam": "Beam",
                 "pointload": "PointLoad", "lineload": "LineLoad", "surfload": "SurfaceLoad",
                 "phase": "Phase"}

        if name == "InitialPhase":
            return self._object("InitialPhase")

        def command(*args):
            if name == self.fail_on or name in self.reject_commands:
                raise RuntimeError(f"{name} failed")
            self._log(name, args)
            if name == "surface":
                return self._object("Polygon")

            if name == "line":
                return [self._object("Point"), self._object("Point"), self._object("Line")]
            if name in kinds:
                return self._object(kinds[name])
            return None
        return command


class _Object:
    def __init__(self, recorder, name):
        object.__setattr__(self, "_recorder", recorder)
        object.__setattr__(self, "_name", name)

    def __repr__(self):
        return self._name

    def setproperties(self, *args):
        names = args[0::2]
        bad = [n for n in names if n in self._recorder.reject]
        if bad:
            raise ValueError(f"unknown property {bad[0]}")
        self._recorder._log(f"setproperties {self._name}", args)

    def __setattr__(self, key, value):
        if key in self._recorder.reject:
            raise ValueError(f"unknown property {key}")
        self._recorder.commands.append(f"set {self._name}.{key} {value!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("exchange", help="exchange file written by the editor (*.plaxis.json)")
    parser.add_argument("--host", default="localhost")
    parser.add_argument("--port", type=int, default=10000)
    parser.add_argument("--password", default="")
    parser.add_argument("--new", action="store_true", help="start a new PLAXIS project first")
    parser.add_argument("--shift-to-origin", action="store_true",
                        help="move the model so that its lower corner is at (0, 0, top)")
    parser.add_argument("--no-phases", action="store_true",
                        help="create the loads but no Staged construction phases")
    parser.add_argument("--dry-run", action="store_true", help="print the commands instead")
    args = parser.parse_args()
    with open(args.exchange, encoding="utf-8") as f:
        data = json.load(f)
    if data.get("format") != "topo-plaxis-1":
        sys.exit(f"{args.exchange}: not a topo-plaxis-1 exchange file")
    shift = (0.0, 0.0, 0.0)
    if args.shift_to_origin:
        points = [x for p in data["plates"] for poly in p["polygons"] for x in poly]
        points += [b[k] for b in data["beams"] for k in ("start", "end")]
        if points:
            shift = (min(p[0] for p in points), min(p[1] for p in points), 0.0)
    if args.dry_run:
        g_i = Recorder()
        report = build(g_i, data, shift, phases=not args.no_phases)
        print("\n".join(g_i.commands))
    else:
        try:
            from plxscripting.easy import new_server
        except ImportError:
            sys.exit("plxscripting is not available: run with the PLAXIS Python distribution "
                     "or `pip install plxscripting`")
        s_i, g_i = new_server(args.host, args.port, password=args.password)
        if args.new:
            s_i.new()
        report = build(g_i, data, shift, progress=lambda text: print(text, flush=True),
                       phases=not args.no_phases)
    report["shift"] = shift
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
