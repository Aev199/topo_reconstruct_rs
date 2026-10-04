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
(2022: Identification/D3d/StructNu12, older: MaterialName/d/nu12); each
material is tried with the known sets in turn and the names that worked are
reported.
"""
import argparse
import json
import sys
import time

# Plate material property sets, newest first: (identification, thickness,
# Poisson ratio, isotropy, unit weight).
PLATE_SETS = [
    dict(name="Identification", d="D3d", nu="StructNu12", iso=("Isotropic", True), gamma="Gamma",
         extra=("MaterialType", "Elastic")),
    dict(name="Identification", d="d", nu="nu12", iso=("IsIsotropic", True), gamma="Gamma",
         extra=("MaterialType", "Elastic")),
    dict(name="MaterialName", d="d", nu="nu12", iso=("IsIsotropic", True), gamma="w", extra=None),
]
BEAM_SETS = [
    dict(name="Identification", gamma="Gamma", i2="I2", i3="I3",
         extra=("MaterialType", "Elastic", "CrossSectionType", "User-defined")),
    dict(name="Identification", gamma="Gamma", i2="I2", i3="I3", extra=("MaterialType", "Elastic")),
    dict(name="MaterialName", gamma="w", i2="I2", i3="I3", extra=None),
]


def plate_properties(m, names):
    props = [names["name"], m["name"]]
    if names["extra"]:
        props += list(names["extra"])
    props += [names["d"], m["d"], names["iso"][0], names["iso"][1], "E1", m["e"], names["nu"], m["nu"]]
    if m["gamma"] > 0:
        # The plate's unit weight: the material's, times the thickness for
        # versions taking a weight per area (w).
        props += [names["gamma"], m["gamma"] if names["gamma"] == "Gamma" else m["gamma"] * m["d"]]
    return props


def beam_properties(m, names):
    props = [names["name"], m["name"]]
    if names["extra"]:
        props += list(names["extra"])
    props += ["E", m["e"], "A", m["a"], names["i2"], m["i2"], names["i3"], m["i3"]]
    if m["gamma"] > 0:
        props += [names["gamma"], m["gamma"] if names["gamma"] == "Gamma" else m["gamma"] * m["a"]]
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


def build(g_i, data, shift=(0.0, 0.0, 0.0), progress=print):
    """Create materials, plates and beams; returns a report."""
    def point(p):
        return tuple(round(p[k] - shift[k], 9) for k in range(3))

    g_i.gotostructures()
    log = []
    plate_mats, beam_mats, used = {}, {}, set()
    for m in data["plate_materials"]:
        plate_mats[m["name"]], index = make_material(g_i, g_i.platemat, PLATE_SETS, plate_properties, m, log)
        used.add(f"plate set {index}")
    for m in data["beam_materials"]:
        beam_mats[m["name"]], index = make_material(g_i, g_i.beammat, BEAM_SETS, beam_properties, m, log)
        used.add(f"beam set {index}")
    plates = beams = 0
    started = time.time()
    total = sum(len(p["polygons"]) for p in data["plates"]) + len(data["beams"])
    done = 0
    for p in data["plates"]:
        for polygon in p["polygons"]:
            surface = last(g_i.surface(*[point(x) for x in polygon]))
            plate = last(g_i.plate(surface))
            if p["material"]:
                plate.Material = plate_mats[p["material"]]
            plates += 1
            done += 1
            if done % 200 == 0:
                progress(f"{done}/{total} objects, {time.time() - started:.0f} s")
    for b in data["beams"]:
        line = last(g_i.line(point(b["start"]), point(b["end"])))
        beam = last(g_i.beam(line))
        if b["material"]:
            beam.Material = beam_mats[b["material"]]
        beams += 1
        done += 1
        if done % 200 == 0:
            progress(f"{done}/{total} objects, {time.time() - started:.0f} s")
    return dict(plates=plates, beams=beams, plate_materials=len(plate_mats),
                beam_materials=len(beam_mats), property_sets=sorted(used),
                rejected_property_sets=log, seconds=round(time.time() - started, 1),
                missing_materials=data.get("missing_materials", []))


class Recorder:
    """A stand-in for PLAXIS Input (`--dry-run`, tests): records every
    command as PLAXIS logs it and returns named objects."""

    def __init__(self, reject=()):
        self.commands = []
        self.counts = {}
        self.reject = set(reject)

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
        kinds = {"platemat": "PlateMat", "beammat": "BeamMat", "plate": "Plate", "beam": "Beam"}

        def command(*args):
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
        report = build(g_i, data, shift)
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
        report = build(g_i, data, shift, progress=lambda text: print(text, flush=True))
    report["shift"] = shift
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
