# Changelog

All notable changes to Stratonix. Versions follow semantic versioning; the description schema
has its own version (`stratonix_schema`).

## [0.4.0] - 2026-11-13 - Fault framework and kinematics (Prompt 4)

### Added
- Fault network as its own object: planar faults by point/strike/dip or trace/dip; tips by
  horizon-relative rules (horizon + offset, fraction through a unit, depth, free); along-strike
  length with rectangular or elliptical tip line; throw profiles (sine, flat-topped, elliptical,
  linear) down dip and along strike.
- Kinematics plugins (dip_slip, vertical_offset) and hanging-wall field plugins (rigid, decay,
  linear_between_faults); ages, abutting, and the frame rule (fixed well or point).
- Representations: fault surfaces clipped to the tip line, fault sticks, hanging-wall and
  footwall cut-off lines refined by bisection, throw tables.
- Checks C3 (no offset outside tips), C4 (throw equals profile), C6 (top/base unoffset),
  C7 (faults in domain, mesh quality), C8 (fault intersections resolved).
- Benchmarks: constant-throw fault, blind fault, two faults with different extents, intersecting
  faults with abutting, tip inside a unit; dip-slip versus vertical offset comparison.
- Prototype regression harness (skipped until the reference tables are supplied).

### Changed (found while testing)
- The throw profile is read at the restored position of the material, at the fault point of the
  same elevation; restoration solves a bracketed 1D equation per fault. Before, a rigid hanging
  wall moved material above a blind fault's upper tip, and the throw at a horizon depended on
  which block the frame held fixed.
- Tips of a younger fault are resolved on the stratigraphy already displaced by older faults.
- Synthetic example: F1 made truly blind (lower tip inside Unit_D); closure with sealing faults.

## [0.3.0] - 2026-11-13 - Structure engine (Prompt 3)

Delivered in the same commit as 0.4.0.

### Added
- Parametric structure builder: stratigraphy from a reference well or point (tops or
  thicknesses); fold primitives (anticline, dome, syncline, basin) with axis azimuth, plunge,
  half-length, limb half-widths, radii, rim synclines; relief per horizon by anchors; flat or
  conformable base; imported surfaces (IRAP classic ASCII, ZMAP+, regular XYZ) and isochore
  stacking; control wells honoured by construction or by a reported correction field.
- Grid importers; horizon grids per block and faulted with gaps; well intersections; closure
  analysis (crest, spill point, closure area); contours.
- Checks C1 (thickness), C2 (crossing), C5 (wells), C9 (closure against targets, reported, never
  tuned); checks report (text and CSV); build fails when a check fails.
- Benchmarks: flat layers, symmetric dome with known closure, asymmetric anticline with known limb
  dips, second control well, imported surface with isochore, closure targets.

## [0.2.0] - 2026-11-13 - Core (Prompt 2)

### Added
- CMake project with one target per architectural layer and a configure-time layer check.
- Core: geometry, regular grids, triangulated meshes, units, CRS record, SHA-256,
  locale-independent number formatting, deterministic `parallel_for`, build information.
- Model descriptions: YAML/JSON reader (libyaml), JSON Schema v1 (embedded), validator,
  migration hook, typed description with evidence class on every parameter, semantic rules
  (E1, U1, S1-S4, F1-F8, W1-W2, R1-R2) with messages naming the field and the rule.
- Model interface (`Model`), flat-layers reference model, interface conformance suite (S1-S9).
- Plugin registry with compile-time self-registration; plugin C ABI header for later.
- Command line: validate, build, plugins, suite, verify, schema, version (section, export,
  package, bench reserved for later stages).
- Run manifest (inputs with SHA-256, resolved description, version, commit, compiler, platform,
  timestamp).
- C API (`stratonix_c.h`) with a C test.
- CI workflow (Windows MSVC + vcpkg, Linux GCC/Clang/sanitisers), CPack zip and NSIS installer.

### Known limitations
- Windows builds are defined in CI but were not run in the build environment (no network, no MSVC).
- The vcpkg baseline is to be pinned on the first networked CI run.
