# MIDAS GTS NX adapter

Status: solver-mesh planning implemented; direct FPN serialization intentionally
not enabled yet.

## Boundary

The MIDAS adapter consumes only:

`topo-reconstruct-solver-mesh-v1`

It must not depend on OpenCASCADE entity tags, Gmsh physical groups, v2 assembly
internals or private-model IDs.

The solver mesh contains:

- global zero-based vertices;
- triangular shell elements;
- two-node bar elements;
- surface regions with stiffness and source-element provenance;
- bar regions with stiffness, source axis/span and source-element provenance.

The Gmsh backend can write this package directly with:

```sh
python scripts/gmsh_occ_trial.py input.json \
  --solver-mesh-output model.solver.json
```

The backend refuses this standalone output if its solver-mesh coverage audit is
not clean.

## Deterministic MIDAS import plan

`scripts/midas_fpn.py` validates the solver package and produces
`topo-reconstruct-midas-gts-plan-v1`.

The plan freezes all choices that do not depend on FPN serialization details:

- node ID = solver vertex index + 1;
- shell element IDs start at 1;
- bar element IDs continue after the last shell element;
- all element IDs are globally unique;
- shell connectivity becomes planned `TRIA` connectivity;
- bar connectivity becomes planned `LINE` connectivity;
- mesh sets are separated by dimensionality and stiffness:
  `SHELL_K<stiffness>` and `BAR_K<stiffness>`;
- each mesh set contains deterministic element and node membership;
- every planned element retains its solver region;
- every region retains source stiffness and source-element provenance.

Example:

```sh
python scripts/midas_fpn.py model.solver.json \
  --plan-output model.midas-plan.json
```

This plan is suitable for testing, diffing and later FPN serialization without
rerunning reconstruction or Gmsh.

## What is verified about GTS NX FPN

MIDAS documentation describes FPN as a text neutral exchange format used by GTS
NX and other MIDAS products, and GTS NX supports importing and exporting mesh
through FPN.

A public parser of an actual GTS NX FPN confirms these record families:

- `NODE`;
- `LINE`;
- `TRIA`;
- `RECT`;
- `MSET`;
- `MSETE`;
- `MSETN`.

It also demonstrates:

- NODE id in field 1 and X/Y/Z in fields 2/3/4;
- LINE node ids in fields 3/4;
- TRIA node ids in fields 3/4/5;
- RECT node ids in fields 3/4/5/6;
- MSET name in field 2;
- MSETE/MSETN count in field 2;
- mesh-set member IDs written on following lines, commonly in blocks of eight.

References used during implementation:

- MIDAS GTS NX user documentation / neutral-format import and export;
- the public GTS NX FPN-to-FLAC3D parser by 木子水星;
- older FX+/Pre-Neutral documentation only as historical context.

## Why the writer is not enabled yet

The exact modern GTS NX field semantics between element ID and node IDs are not
sufficiently documented publicly for us to call a generated file reliable.
Older FX+/Pre-Neutral syntax is closely related but differs from modern GTS NX
at least in mesh-set membership records and may differ in element/property
fields.

Guessing here would create the worst kind of failure: a syntactically accepted
file with incorrectly assigned properties or element interpretation.

Therefore `scripts/midas_fpn.py` currently creates the complete import plan but
does not claim to create an import-ready FPN.

## One-file calibration path

The adapter can inspect a tiny FPN exported by the target GTS NX version:

```sh
python scripts/midas_fpn.py \
  --inspect-fpn tiny-reference.fpn \
  --inspection-output tiny-reference.profile.json
```

The inspector detects text encoding, section headings, record families, field
counts and a few record examples.

The ideal calibration fixture is deliberately tiny:

- one triangle shell with three nodes;
- one two-node line element;
- different mesh-set names for shell and line;
- known property/stiffness IDs if GTS NX requires them in the element records.

Once that fixture is available, the remaining work is serialization mapping,
not geometry reconstruction.

## Acceptance for the future FPN writer

Do not mark the writer verified until all of the following pass:

1. generated file imports into the target GTS NX version without repair;
2. node and element counts match the solver plan;
3. shell and bar mesh sets contain exactly the planned element IDs;
4. no shell/bar group names collide;
5. element connectivity is unchanged;
6. stiffness/property assignment is either imported correctly or deliberately
   left unassigned and reported as such;
7. export-back-to-FPN can be compared structurally with the generated input;
8. the same writer passes a small regression fixture before the private model.

Until then, `format_profile.writer_verified=false` is intentional.
