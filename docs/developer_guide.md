---
title: Stratonix developer guide
doc_id: STX-DG
version: 1.0
date: 2026-11-13
---

# Stratonix developer guide

# 1. Architecture

Stratonix is one C++20 library (`libstratonix`, namespace `stratonix`) split into CMake targets by
layer, a command-line program and a C API. A layer may link only layers below it;
`cmake/LayerCheck.cmake` fails the configure step otherwise, and the `layer_check` test repeats the
check.

| Level | Target | Holds |
|---|---|---|
| 0 | `stratonix_core` | geometry, grids, meshes, units, CRS record, wells, root finders, threads, SHA-256, logging, number formatting |
| 1 | `stratonix_desc` | YAML/JSON reader, embedded JSON Schema and validator, typed description, semantic rules, migration |
| 1 | `stratonix_registry` | plugin registry (kinds, entries, factories) and the plugin C ABI header |
| 2 | `stratonix_model` | the `Model` interface, the flat-layers reference model, the interface conformance suite |
| 3 | `stratonix_engines` | structure builder, structure field, fault network and kinematics, parametric model, wells |
| 4 | `stratonix_repr` | grids with gaps, per-block grids, fault surfaces and sticks, cut-offs, throws, well tops, sections, voxels, closure, contours |
| 5 | `stratonix_checks` | checks C1-C9, run manifest |
| 5 | `stratonix_io` | file formats, exporters, importers |
| 5 | `stratonix_figures` | display list, SVG, Cairo raster, PNG/TIFF encoders, PPTX, software 3D renderer, the figure set |
| 6 | `stratonix_package` | build pipeline, exports, figures, benchmarks, delivery package |
| 7 | `stratonix_capi`, `stratonix_c`, `stratonix` | C API (static and shared), command line |

Choices that keep later stages open (specification section 2): everything reads a model only
through the `Model` interface; the fault network is its own object; every parameter is a typed
`Param` with value, unit, evidence and source; representations are derived, never stored as the
model; plugins are registered by kind; figures are drawn on an abstract display list; the C API
uses opaque handles and status codes.

# 2. The model interface and its test suite

`include/stratonix/model/model.hpp`. A model answers, for any batch of points, the unit, the fault
block (bit f set = hanging wall of fault f) and one scalar per horizon (positive above the
horizon, in metres of elevation along the vertical in the model's frame); `scalar_in_block`
evaluates a horizon as if the point were in a given block; `evaluate_fault` gives signed distance,
extent coordinate and displacement for a fault; `restore` maps present-day points to restored
ones. Models are immutable after building and every call is `const` and thread-safe.

`model_suite.hpp` runs nine conformance checks on any model (S1 consistent info, S2 field sizes,
S3 unit from scalars, S4 no crossing, S5 per-block scalars, S6 block bits against fault distance,
S7 thread and batch invariance, S8 fixed block unmoved, S9 scalar monotone down a vertical).
`stratonix suite <description>` runs it on a description's model; every new engine must pass it.

# 3. Build and dependencies

- C++20, CMake 3.24 or later, Ninja or Visual Studio 2022. Presets: `linux-gcc-release`,
  `linux-gcc-debug`, `linux-clang-release`, `linux-asan` (address and undefined-behaviour
  sanitisers), `windows-msvc-release`, `windows-msvc-debug`.
- Dependencies: libyaml (MIT), zlib (zlib licence), Cairo (LGPL 2.1, dynamically linked; optional:
  without it only SVG figures are written). On Windows they come from vcpkg (`vcpkg.json`); the
  vcpkg baseline must be pinned on the first networked build (`vcpkg x-update-baseline
  --add-initial-baseline`). No GPL code. See `THIRD_PARTY.md`.
- Tests: `ctest` runs the unit tests (own minimal harness, `tests/stx_test.hpp`), the model suite,
  parity tests of optimised kernels, round trips of every exporter and importer, golden SHA-256
  of the synthetic Petrel set and package (so any platform writing different bytes fails), the
  C API test and example, command-line tests, and independent Python checks of figure and export
  files (skipped when NumPy, Pillow, tifffile or python-pptx are missing).
- Continuous integration: `.github/workflows/ci.yml` builds and tests Linux (GCC, Clang, sanitisers)
  and Windows (MSVC, release and debug), packages the portable zips, the Windows installer and the
  Linux archive, smoke-tests the Windows zip on a stripped PATH, and on a version tag publishes a
  release with `SHA256SUMS.txt`.
- API reference: `doxygen docs/Doxyfile` where Doxygen is installed; otherwise
  `python3 tools/gen_api_reference.py` writes `docs/api_reference.md` from the same header comments.

# 4. Extending Stratonix

Every extension is a plugin: derive from the kind's interface in the layer that owns it, register
it with `STRATONIX_REGISTER(id, kind, name, version, description, factory)` in its `.cpp`, and add
`STRATONIX_PLUGIN_ANCHOR(id);` to that layer's registration function so the linker keeps it.
`stratonix plugins` must then list it. The registry is filled once by `register_all_plugins()`.

## 4.1 A structure builder

Derive from `StructureBuilder` (`engines/builder.hpp`): `accepts(description)` says whether it can
build the description; `build(description, inputs, log)` returns a `ModelPtr`. Register with kind
`StructureBuilder`. The model it returns must pass the model suite (add a test that runs
`run_model_suite` on it) and the checks C1-C9 that apply. `build_model` picks the first builder that
accepts a description, or the one named with `--builder`.

## 4.2 A kinematics model or hanging-wall field

Derive from `Kinematics` (`engines/faults.hpp`): `slip_vector(down_dip, dip)` returns the slip
direction scaled so its vertical component is -1 (the throw then moves material by `throw` metres
vertically). Derive from `HangingWallField` for a new decay of displacement with distance:
`scale(sigma, next, influence)` in [0, 1]. Register with kinds `Kinematics` or `HangingWall`; the
description's `slip_mode` and `hanging_wall.type` select them by name. Add a test in
`test_engines.cpp` with a case that has an analytical answer, and run the synthetic example's
checks (C3, C4) with it.

## 4.3 An importer

Derive from `io::Importer`: `data_kind()` ("grid", "points", "lines", "meshes", "tops") and
`read(path)` returning `io::ImportedData` with coordinates as in the file and `z_is` set. Register
with kind `Importer`. Add a round trip against the matching exporter in `test_interop.cpp`.

## 4.4 An exporter

Derive from `io::Exporter`: `id()`, `petrel()`, and `write(context, folder)` returning the files
written, each with format, Petrel import option and a one-line description (they go into the
export README). Take every representation from the `ExportContext` (it caches grids, fault
surfaces, cut-offs and tops). Rules: z as elevation unless the format says otherwise; `'\n'` line
endings; numbers through `fmt_fixed` (locale independent); names through `petrel_name`. Add a
writer and a reader, a round trip against the model, and an entry in the independent Python check
(`tests/tools/check_petrel_set.py` or `check_neutral_set.py`). If the exporter is part of the
Petrel set, add it to `petrel_set_ids()` and refresh `tests/data/petrel_synthetic_sha256.txt`
deliberately (the test writes the new list to the test folder).

## 4.5 An optimised kernel

Keep the reference implementation; add the optimised one behind the same function, selected by
`optimised_kernels()`. Parallelise only through `parallel_for` (fixed chunks, results independent
of the thread count). Add a parity test in `test_parity.cpp` on the synthetic models and randomised
points, including points on fault planes, at tips and outside the domain, with a stated
tolerance, and a benchmark row in `bench.cpp`. See `docs/performance_and_threading.md`.

## 4.6 A check

Derive from `Check` (`checks/checks.hpp`), return one `CheckResult` per subject with value,
tolerance (and where it came from) and pass/fail; register with kind `Check`. Checks report; they
never change the model.

# 5. Conventions

- No site, field or project names in code, examples or documentation.
- Every number in a description has a unit, evidence class and source; defaults that Stratonix
  supplies are labelled "Stratonix default" wherever they are reported.
- Output folders must be new or empty; nothing overwrites an issued output.
- Timestamps appear only in manifests.
