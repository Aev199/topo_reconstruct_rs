# AGENTS.md — topo_reconstruct_rs

## Mission

Reconstruct robust structural geometry from imperfect FE meshes for downstream
geotechnical workflows. The target is valid, connected, explainable geometry
suitable for a quality MIDAS mesh and PLAXIS geometry/remeshing. Exact
one-to-one reproduction of the source FE discretization is secondary.

Do not optimize for one private example. Rules must be geometry- and
topology-based, invariant to IDs, input ordering, translation/rotation, and
reasonable changes of scale.

## Canonical project context

Read these before substantial geometry work:

1. `docs/DEV_STATE.md` — compact current status and next development batch.
2. `docs/GEOTECHNICAL_GEOMETRY.md` — current acceptance contract and verified
   full-model findings.
3. `docs/RECONSTRUCTION_V2.md` — detailed reconstruction requirements,
   architecture, provenance and invariants.
4. Relevant source/tests only after the above. Do not reconstruct project state
   from chat history when repository state is sufficient.

When documents disagree, prefer the newer verified statement in
`GEOTECHNICAL_GEOMETRY.md`, then update `DEV_STATE.md` after the work.

## Non-negotiable geometry rules

- Valid topology is more important than matching the original FE tessellation.
- Preserve structural/material regions, meaningful openings, structural
  junctions, loads/properties/provenance when available.
- Never hard-code source element/node IDs or coordinates from a private fixture.
- Shared geometric junctions must have shared topological identity; coordinate
  coincidence alone is insufficient.
- Panels must be planar by construction. Real folds become multiple panels.
- Do not create degenerate faces, self-intersections, duplicate faces/edges,
  parasitic short edges, or accidental near-coincident boundaries.
- Rods/axes may and should be split at directly connected structures. Global
  straightness over unrelated spans is not required.
- Adaptive repair tolerance is allowed and may exceed 50 mm when justified by
  local scale/context. Numerical precision, repair tolerance and target mesh
  size are separate concepts.
- Do not delete or merge a feature merely because it is small. Any
  simplification needs explicit geometric/topological evidence and provenance.
- A good source-deviation metric does not override a failed topology/geometry
  audit.
- MIDAS/PLAXIS suitability is the product goal; a prettier reconstruction is not
  a goal by itself.

## Development loop

Work in coherent batches instead of chat-sized micro-edits.

For each non-trivial defect:

1. Define the defect and an objective acceptance condition.
2. Reproduce it with the smallest useful synthetic or extracted regression case
   when possible.
3. Add/adjust the regression test before or with the fix.
4. Implement the smallest general rule that solves the class of defects.
5. Run targeted tests first, then the normal Rust suite.
6. Run independent Python/audit checks when the changed subsystem has one.
7. Run the private full model only after a coherent batch is green locally.
8. Compare before/after metrics and inspect only meaningful residuals.
9. Commit a finished batch; avoid a chain of speculative micro-commits.
10. Update `docs/DEV_STATE.md` when priorities, verified metrics, or blockers
    materially change.

Do not weaken an acceptance threshold simply to make the current fixture pass.
If a hypothesis fails, keep diagnostics and change the hypothesis.

## Test tiers

### Tier A — fast development

Use targeted Rust tests while editing, followed by:

```sh
cargo test --offline
```

Synthetic regression tests should cover invariance where relevant: reordered
input, transformed coordinates, reversed normals/orientations, and scale
variation.

### Tier B — subsystem audit

Run the independent checker(s) relevant to the changed subsystem, for example:

```sh
python3 scripts/check_v2_assembly.py <report.json>
python3 -m unittest discover -s scripts -p test_v2_global_geometry.py
python3 scripts/check_v2_global_geometry.py <report.json> --output <audit.json> --strict
```

Do not install optional audit dependencies or expand CI merely to run a check
that can be performed in the existing development environment.

### Tier C — private full-model gate

Use the private full `скала1` fixture only when the batch is already green on
Tier A/B. The fixture and full diagnostic reports must not be committed.

A full-model run is an integration/acceptance test, not the inner development
loop.

## GitHub Actions / free-account policy

This repository is public. Standard GitHub-hosted runners for public
repositories do not consume the private-repository minutes allowance, but CI
still costs developer time and artifact/cache storage. Keep Actions lean.

- Do not use Actions as the normal edit-test loop.
- Do not add scheduled workflows, large matrices, macOS jobs, larger runners, or
  repeated full-model jobs without a concrete need.
- Never put the private full-model fixture in Actions.
- Documentation/state-only commits should use `[skip ci]`.
- Iterative commits that were already validated locally may use `[skip ci]`
  when no cross-platform build is needed.
- Prefer rerunning only a failed job rather than an entire successful matrix.
- Upload binaries/artifacts only when they are useful to the user; keep
  retention short.
- Avoid duplicate Windows/Linux release builds merely to prove geometry logic;
  geometry/unit tests should normally run on one platform, with cross-platform
  builds reserved for integration/release checks.
- Before changing CI, check the current workflow and expected storage/runtime
  impact.

## Mobile-first collaboration

The user operates primarily from the ChatGPT mobile app. Do not make them act as
a terminal operator unless unavoidable.

Use connected GitHub access to inspect repository state directly. Return compact
status after a batch:

- what changed;
- tests/audits run and their result;
- before/after geometry metrics that matter;
- remaining blocker;
- commit SHA.

Do not dump routine logs or require repeated "continue" prompts. Continue
through implementation/test/review inside the same coherent batch when possible.

## Safety for repository writes

- Never force-push unless the user explicitly requests it.
- Preserve unrelated work.
- Do not rewrite history just to make it cleaner.
- Keep commits reviewable and scoped to one coherent batch.
- Use `[skip ci]` for documentation-only infrastructure commits unless CI is
  intentionally being tested.
