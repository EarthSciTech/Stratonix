---
title: Stratonix developer guide
doc_id: STX-DG
version: 4.1
date: 2026-11-28
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
| 3 | `stratonix_engines` | structure builder, structure field (series, salt, compaction), fault network and kinematics (curved faults, trishear, fault-parallel flow), event timeline and states (`history.hpp`), parametric model, implicit model (dense and partition-of-unity solvers, series, isochore horizons), hybrid conditioning (`conditioning.hpp`), wells |
| 4 | `stratonix_repr` | grids with gaps, per-block grids, fault surfaces and sticks, cut-offs, throws, well tops, sections, voxels, closure, contours, thickness and termination lines, flexural unfolding, fault rock and curvature (`faultrock.hpp`), corner-point and tetrahedral grids, properties, transmissibility and upscaling (`grid.hpp`) |
| 5 | `stratonix_checks` | checks C1-C25, balance rows, fault-rock, stress and grid inputs from a description (`make_grid_products`), run manifest |
| 5 | `stratonix_io` | file formats, exporters (GRDECL, grid VTK, RESQML through FESAPI when built with it), importers |
| 5 | `stratonix_figures` | display list, SVG, Cairo raster, PNG/TIFF encoders, PPTX, software 3D renderer, the figure set |
| 6 | `stratonix_package` | build pipeline, exports, figures, benchmarks, delivery package, ensembles (`uncertainty.hpp`) |
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

**Capabilities (2026.2).** Layers above the engines must not include engine headers to reach
engine-specific queries. A model offers them through `capability(name)`, returning a pointer to an
interface declared in the model layer, or null:

| Name | Interface | Offered by | Used by |
|---|---|---|---|
| `structure` | `Structure` (engines) | parametric | checks C12-C17, exporters (porosity), history |
| `faults` | `FaultNetwork` (engines) | parametric | checks C11, C12 |
| `stratigraphy` | `StratigraphyQuery` (model) | parametric | termination lines, isopachs (repr) |
| `column` | `ColumnQuery` (model) | parametric | fault rock: restored stratigraphic positions (repr) |
| `data_misfit` | (engines) | implicit, conditioned | reports, C21 |
| `conditioning` | `ConditioningResult` (engines) | conditioned parametric | reports, C22 |
| `posterior` | `PosteriorResult` (engines) | conditioned parametric, `method: bayesian` | posterior tables and figures, C28 |
| `porosity` | `PorosityQuery` (model) | parametric with compaction | grid porosity (repr) |

A new engine that cannot answer a capability returns null; the callers then report the dependent
attribute as missing, never a default.

**Restoration-based evaluation.** Every mechanism added in v2 (curved faults, shear and
fault-parallel-flow hanging walls, trishear, series relations, salt, compaction) has an exact
inverse. Evaluation restores a point through the fault network (`FaultNetwork::restore`). It then
evaluates the stratigraphy in the restored frame (`Structure::z_all`). New mechanisms must keep
both directions consistent; `test_engines`, the model suite (S8, S9) and the analytical tests
check it.

# 3. Build and dependencies

- C++20, CMake 3.24 or later, Ninja or Visual Studio 2022. Presets: `linux-gcc-release`,
  `linux-gcc-debug`, `linux-clang-release`, `linux-asan` (address and undefined-behaviour
  sanitisers), `windows-msvc-release`, `windows-msvc-debug`.
- Dependencies: libyaml (MIT), zlib (zlib licence), Cairo (LGPL 2.1, dynamically linked; optional:
  without it only SVG figures are written). FESAPI (Apache 2.0, optional, `STRATONIX_WITH_FESAPI`
  and `FESAPI_ROOT`) for the RESQML exporter; the CI job `resqml` builds it from source with HDF5,
  Minizip and Boost, and reads the result back with resqpy. On Windows they come from vcpkg (`vcpkg.json`); the
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

The implicit builder (`engines/implicit.cpp`, plugin `implicit-rbf`) is an example of a builder
that reuses the parametric engine: it interpolates a scalar field from data, turns it into
horizon surfaces in the restored frame, and hands them to `ParametricModel` with the fault network;
`ImplicitModel` adds provenance and the capability `data_misfit`.

## 4.2 A kinematics model or hanging-wall field

Derive from `Kinematics` (`engines/faults.hpp`): `slip_vector(down_dip, dip)` returns the slip
direction scaled so its vertical component is -1 (the throw then moves material by `throw` metres
vertically). Derive from `HangingWallField` for a new decay of displacement with distance:
`scale(sigma, next, influence)` in [0, 1]. Register with kinds `Kinematics` or `HangingWall`; the
description's `slip_mode` and `hanging_wall.type` select them by name. Add a test in
`test_engines.cpp` with a case that has an analytical answer, and run the synthetic example's
checks (C3, C4) with it.

### Layer-rigid restoration (1.2)

With `layer_rigid` hanging walls the displacement of material depends on its stratigraphic level
s (horizon index, fractional inside a unit, from the unfaulted structure), not on its position.
`FaultNetwork::solve_layer` finds the one unknown s of a present point by bracketing: s must equal
the level of the point restored through every fault with each fault's throw read at level s. Each
fault keeps a table of where every level meets its plane in its own frame (exact at the horizons,
sampled every 10 m down dip between them), resolved oldest first. With vertical slips the column
does not move, so its horizons are read once per point. The build refuses a fault whose throw
changes faster than the layers thicken at the fault (layers would overlap).

### Section sheets (1.2)

`figures/sheet.cpp` draws a sheet from `SheetSpec` (built from the description by the pipeline).
Per column of the section it finds the fault crossings, splits the column into block intervals and
finds each horizon by root finding in that block, so offsets are sharp; unit fills are polygons
per block and unit, horizon lines are clipped to their block. Key numbers are computed from the
same geometry and the model interface (well intersections, restoration, fault evaluation); the
unfaulted levels come from the engine through a callback. Since 1.3 every key number carries a
basis and an `on_sheet` flag; the pipeline writes all of them to the summary table
(`quantity,value,basis`) and the sheet draws the flagged ones right-aligned in the bottom-right
corner. A new key-number kind is a branch in `section_sheet` (value, text and basis placeholders),
a name in the description rules (`description.cpp`, rule O1) and the schema enum, a default
rounding in `pipeline.cpp`, and a check in `test_sheet.cpp` against an independent calculation.

### Events and states (2026.2)

`engines/history.hpp`. `event_ages(d)` lists the ages of a description's events, oldest first.
`state_description(d, full, t)` cuts a description at age t:

- it removes the later horizons, units, series, compaction laws and events;
- it scales the events under way;
- it fixes the reference depths from the full model.

`build_state` builds the result with the normal builder, so a state is a full model. A new
mechanism that can be timed needs three things:

- an event type (schema enum, `EventDesc`, rule F26);
- a branch in `state_description` that scales it;
- a test in `test_v2_history` that restates an untimed example as events and compares the
  horizons.

### Fault rock and stress (2026.2)

`repr/faultrock.hpp` takes plain inputs (`FaultRockInput`, `StressInput`), never description
types. `fault_rock_input` and `stress_input` (checks layer) convert a description; checks,
figures and exports all use them, so they see the same values. Every attribute is NaN with a
reason in `missing` when an input is absent. A new attribute is:

- a field in `FaultRockPoint`;
- its computation in `fault_rock_at`;
- a figure call in `standard_figures`;
- a column in the CSV exporter and a VTK array;
- a summary row in `fault_rock_rows`;
- an analytical case in `test_v2_faultrock`.

### Grids and properties (2026.3)

`repr/grid.hpp` takes plain inputs (`GridOptions`, `PropertyInputs`, `FaultRockInput`).
`make_grid_products` (checks layer) builds the grid through the `corner_point` extractor plugin,
the properties, the transmissibility and the tetrahedral grid once per run; checks C23-C25,
`grid/cells.csv`, the grid figures and the grid exporters share them (`CheckContext::grid_products`).
A new property is:

- a branch of `grid_properties` (or a transform), with unit, evidence class and source, and a
  missing count with the reason;
- a column in `grid/cells.csv` and a keyword in GRDECL (automatic for properties in the set);
- an analytical case in `test_v3_grids`.

Other grid extractors (kind `Extractor`, interfaces `CornerPointExtractor`,
`TetrahedralExtractor`) and upscalers (kind `Upscaler`, `UpscalerPlugin`) register like any plugin.

### Ensembles (2026.3)

`package/uncertainty.hpp`. An ensemble edits the description document by path
(`unc::set_parameter`, scenario `set`/`remove`), validates it again and builds it through the
usual pipeline. Members run in parallel; results are written in member order. Rules for code
touched by ensembles:

- random numbers only from `unc::Rng` streams seeded by `member_seed(seed, member)`;
- no use of `std::` distributions (their output differs between standard libraries);
- outputs written in member order, numbers through `fmt_fixed`/`fmt_shortest`;
- `test_v3_uncertainty` U9 compares files from 1 and 4 threads.

### Real data (2026.4)

Readers live in the io layer and read files as delivered: `io/welldata.hpp` (LAS 2.0/3.0, DLIS),
`io/seismic.hpp` (SEG-Y, picks, sticks), `io/maps.hpp` (shapefile, GeoJSON, GeoPackage, GeoTIFF,
ESRI ASCII), `io/velocity.hpp` (the description's velocity model). Every reader returns its
values with a `report` of what it met (nulls, units, encodings, lines not read); a reader never
assumes a value the file does not state. Inputs go through `io::read_input`, which an ensemble
caches (key: format, path, size, modification time). The engines (`engines/surveys.hpp`,
`engines/ages.hpp`, `engines/timedepth.hpp`) take plain values. New readers are cross-checked in
CI against an independent reader used as a test tool only (lasio, dlisio, segyio, GeostatsPy).

### Geostatistics (2026.4)

`engines/geostat.hpp` works on points in geological-time coordinates; `package/geostat_props.cpp`
maps a unit's cells to (u, v, w) and back. Rules for code that produces simulated values:

- the stream is SplitMix64 with Marsaglia's polar method; `exp` and `log` that enter simulated
  values are the deterministic versions in `geostat.cpp` (IEEE + - * / and sqrt only);
- neighbours are sorted by (scaled distance, node index); no unordered container is iterated;
- `test_v4_scale` S7 compares a realisation's checksum with the value from the Linux build, so CI
  on Windows is the cross-platform test.

### Python (2026.4)

`python/stratonix/__init__.py` binds the C API with `ctypes`; a function added to the C API gets
its `argtypes`/`restype` there and a test in `python/tests`. The bindings never call the C++
library directly, so the C API's stability is theirs.

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

RESQML exports are validated in CI against the Energistics RESQML 2.0.1 schemas, which are licensed
and never committed or hosted. The owner packs them into the repository secrets
`STX_RESQML_XSD_1` to `_8` with `tools/resqml_xsd_secrets.py pack <package>` (documentation
removed, one tar.xz archive, base64 parts of at most 48 000 characters, a header with the part
count and the archive checksum); the `resqml` job rebuilds the folder with `unpack` and refuses
parts that do not rebuild the packed archive. `tests/tools/check_resqml_xsd_secrets.py` tests the
packing on a made-up package.

## 4.5 An optimised kernel

Keep the reference implementation; add the optimised one behind the same function, selected by
`optimised_kernels()`. Parallelise only through `parallel_for` (fixed chunks, results independent
of the thread count). Add a parity test in `test_parity.cpp` on the synthetic models and randomised
points, including points on fault planes, at tips and outside the domain, with a stated
tolerance, and a benchmark row in `bench.cpp`. See `docs/performance_and_threading.md`.

Shared representations (2026.4.1). A build's checks, figures and exports share one
`repr::FaultedSurfaceCache`: faulted horizon surfaces, fault surfaces and cut-off lines, each
computed once for its exact arguments (`CheckContext::surfaces()`, `FigureContext::surfaces`,
`ExportContext::share_surfaces`). Take fault surfaces and cut-offs from it rather than calling
`repr::fault_surface` again. Where only the distance to a fault (or its extent coordinate) is
needed, call `Model::evaluate_fault_geometry`: it gives the same values as `evaluate_fault`
without the displacement, which needs a restoration. Keep per-cell data compact: grids are sized
in bytes per cell (about 150 with GRDECL); write large tables as they are made (as GRDECL and
`grid/cells.csv` are) and build optional products on demand (the tetrahedral grid). Check a change
with `STRATONIX_TIMING=1`, which prints the time and memory after each check and grid step.

## 4.6 A check

Derive from `Check` (`checks/checks.hpp`), return one `CheckResult` per subject with value,
tolerance (and where it came from) and pass/fail; register with kind `Check`. Checks report; they
never change the model. A check that only applies to some descriptions (for example C17 to models with events,
C19 to models with `fault_rock`) returns no results otherwise, so the checks report and the golden
hashes of older descriptions do not change. Set `warn_only` for a check that reports a risk rather
than an error (C10, C18).

# 5. Conventions

- No site, field or project names in code, examples or documentation.
- Every number in a description has a unit, evidence class and source; defaults that Stratonix
  supplies are labelled "Stratonix default" wherever they are reported.
- Output folders must be new or empty; nothing overwrites an issued output.
- Timestamps appear only in manifests.
