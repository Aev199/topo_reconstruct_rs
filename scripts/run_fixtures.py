"""Tier C gate: run the reconstruction on every private fixture and audit it.

Usage: python3 scripts/run_fixtures.py FIXTURE_DIR BINARY OUT_DIR [--jobs N] [--only NAME...]

For each `*.txt` in FIXTURE_DIR (never committed) the binary writes the
mesh-preview report (frame cached per fixture in --cache or OUT_DIR), then the
assembly checker, the strict global audit and the PLAXIS profile run on
it. A summary of the metrics that matter is written to OUT_DIR/summary.json
and printed as a table. Reports stay in OUT_DIR (never commit them).
"""
import argparse
import concurrent.futures
import json
import math
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def angles(mesh):
    vertices = mesh["vertices"]
    counts = [0, 0, 0]
    smallest = 180.
    for t in mesh["triangles"]:
        p = [vertices[i] for i in t["vertices"]]
        tri_min = 180.
        for k in range(3):
            a, b, c = p[k], p[(k + 1) % 3], p[(k + 2) % 3]
            u = [b[i] - a[i] for i in range(3)]
            v = [c[i] - a[i] for i in range(3)]
            nu, nv = math.sqrt(sum(x * x for x in u)), math.sqrt(sum(x * x for x in v))
            if nu == 0 or nv == 0:
                tri_min = 0.
                continue
            cos = max(-1., min(1., sum(x * y for x, y in zip(u, v)) / (nu * nv)))
            tri_min = min(tri_min, math.degrees(math.acos(cos)))
        smallest = min(smallest, tri_min)
        for i, limit in enumerate((1., 5., 20.)):
            if tri_min < limit:
                counts[i] += 1
    return counts, smallest


def run(model, binary, out, cache, extra, legacy=False):
    name = model.stem
    report = out / f"{name}.json"
    started = time.time()
    if legacy:
        cmd = [str(binary), str(model), "--v2-mesh-preview-json", str(report),
               "--v2-frame-cache", str(cache / f"{name}.frame.json"), *extra]
    else:
        cmd = [str(binary), str(model), "-o", str(report), "--mesh",
               "--frame-cache", str(cache / f"{name}.frame.json"), *extra]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    result = {"model": name, "seconds": round(time.time() - started, 1)}
    if proc.returncode != 0:
        result["error"] = proc.stderr.strip()[-500:]
        return result
    data = json.loads(report.read_text())
    topology = data["topology"]
    model_data = topology["preview"]
    result["surfaces"] = len(model_data["surfaces"])
    result["axes"] = len(topology.get("axis_assembly", {}).get("axes", []))
    mesh = data.get("mesh")
    if mesh:
        result["mesh_valid"] = mesh["topology_valid"]
        result["mesher_ready"] = mesh["external_mesher_ready"]
        result["triangles"] = len(mesh["triangles"])
        counts, smallest = angles(mesh)
        result["below_1_5_20"] = "/".join(map(str, counts))
        result["min_angle"] = round(smallest, 2)
    else:
        result["mesh_error"] = data.get("mesh_error")
    checks = {
        "assembly": [sys.executable, str(HERE / "check_v2_assembly.py"), str(report)],
        "strict": [sys.executable, str(HERE / "check_v2_global_geometry.py"), str(report),
                   "--output", str(out / f"{name}.audit.json"), "--strict"],
        "plaxis": [sys.executable, str(HERE / "check_plaxis_profile.py"), str(report),
                   "--output", str(out / f"{name}.plaxis.json")],
    }
    for key, check in checks.items():
        proc = subprocess.run(check, capture_output=True, text=True)
        result[key] = proc.returncode == 0
        if proc.returncode != 0 and key != "plaxis":
            result[f"{key}_error"] = (proc.stderr or proc.stdout).strip()[-300:]
    audit_path = out / f"{name}.audit.json"
    if audit_path.exists():
        audit = json.loads(audit_path.read_text())
        kinds = {}
        for issue in audit.get("issues", []) + audit.get("point_and_bar_issues", []):
            kinds[issue.get("kind", issue.get("type", "?"))] = kinds.get(issue.get("kind", issue.get("type", "?")), 0) + 1
        result["audit_issues"] = kinds
    plaxis_path = out / f"{name}.plaxis.json"
    if plaxis_path.exists():
        plaxis = json.loads(plaxis_path.read_text())
        result["plaxis_items"] = {k: v for k, v in plaxis.get("counts", {}).items() if v}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixtures", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--cache", type=Path, help="frame cache directory (default OUT)")
    parser.add_argument("--only", nargs="*")
    parser.add_argument("--legacy-cli", action="store_true",
                        help="binary built before 2026-10-04 (--v2-mesh-preview-json)")
    parser.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                        help="further arguments for the binary")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    cache = args.cache or args.out
    cache.mkdir(parents=True, exist_ok=True)
    models = sorted(args.fixtures.glob("*.txt"), key=lambda p: p.stat().st_size, reverse=True)
    if args.only:
        models = [m for m in models if any(o in m.stem for o in args.only)]
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        results = list(pool.map(lambda m: run(m, args.binary, args.out, cache, args.extra, args.legacy_cli), models))
    (args.out / "summary.json").write_text(json.dumps(results, ensure_ascii=False, indent=2))
    for r in results:
        print(json.dumps(r, ensure_ascii=False))


if __name__ == "__main__":
    main()
