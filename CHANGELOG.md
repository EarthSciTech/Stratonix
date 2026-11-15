# Changelog

All notable changes to Stratonix. Versions follow semantic versioning; the description schema
has its own version (`stratonix_schema`).

## [Unreleased] - v2 (2026.2) in progress

### Added
- Stage 2A specification (technical note STX-TN-05, kept with the project documents): curved
  faults with vertical- and inclined-shear hanging walls, drag, splays and relay ramps,
  displacement-length scaling, and reverse, oblique and strike-slip movement.
- Curved faults (Prompt 16): `geometry.type: listric` (exponential profile onto a detachment)
  and `ramp_flat` (profile polyline), with `hanging_wall: vertical_shear` or `inclined_shear`
  (shear angle and antithetic or synthetic sense) and the displacement given as `throw.heave`.
  The hanging wall follows the inclined-shear construction exactly, in both directions, so
  sections, maps, exports and checks work on curved faults unchanged. Per horizon the build
  reports heave, throw, rollover width and maximum rollover dip. Rules F10 (tips free), F11
  (fixed point in the footwall), F12 (shear lines cross the fault once) and F14 (curved faults
  are normal dip-slip, take heave, no down-dip profile) refuse what the construction cannot
  invert. New tolerance `checks.tolerances.rollover`. Worked example `listric_rollover.yaml`.
- Fault interaction (Prompt 17): `drag: {type: normal | reverse, fraction, width}` bends both
  walls along the slip, scaled by the local throw (rule F13 keeps horizons from crossing the
  fault); `abuts: {..., transfer: true}` makes a linked splay, whose slip follows the branch line
  so the block beyond it stays in contact with both faults (rule F15 refuses branch lines
  plunging less than 10 degrees); the branch line is reported (trend, plunge, point);
  `relay: {fault: F2}` gives two overlapping faults complementary throw ramps over the overlap.
  Relations are recorded in `FaultInfo.relations`.
- Checks C10 (displacement-length scaling, a warning; range `checks.fault_scaling`, default
  0.001 to 0.1) and C11 (fault-system consistency across branch lines and relays). The checks
  report is now in numeric order of the check ids.
- Figures: a displacement map on every fault surface with the horizons' cut-offs, and the throw
  along strike of each linked fault system. Worked example `fault_system_relay_splay.yaml`.

- Movement sense (Prompt 18): `throw.sense: normal | reverse` on dip-slip and vertical-offset
  faults (absent means normal, as before); `slip_mode: oblique` with the net slip `throw.slip`
  and a `rake` (Aki and Richards); `slip_mode: strike_slip` with `throw.slip` and
  `sense: left_lateral | right_lateral`. Rule F16 keeps faults with a strike component of slip on
  rigid hanging walls. Reverse faults report the repeated thickness per horizon; faults whose
  slip is not plain normal dip-slip report their slip vector. Worked examples
  `thrust_repeated_section.yaml`, `strike_slip_offset.yaml` and `oblique_fault.yaml`.
- Overlaps under reverse faults: the faulted surface keeps the upper sheet (now also in the
  optimised tracer) and counts overlap nodes; per-block grids carry each block's own sheet, so
  the hanging wall and footwall over a thrust are exported separately; sections draw each
  block's horizon. `FaultInfo.overlapping` marks faults that can override their footwall.

- Stage 2B specification (technical note STX-TN-06): trishear, fault-bend folds by
  fault-parallel flow, detachment folds, parallel and similar fold mechanics, en echelon arrays.
- Trishear (Prompt 20): `fold: {type: trishear, apical_angle, propagation_slip_ratio}` on a
  planar dip-slip fault with an upper tip in the model folds the layers ahead of the propagating
  tip by the symmetric linear velocity field of Zehnder and Allmendinger (2000), integrated in
  both directions by fourth-order Runge-Kutta. Reverse slip gives fault-propagation folds and
  normal slip forced folds. The build reports the slip, the initial and present tip depths and
  each folded horizon's maximum forelimb dip. Rule F17 sets the admissible range (upper tip in
  the model, constant slip behind it, rigid hanging wall, fixed point outside the zone's sweep).
  Worked example `trishear_fold.yaml`.
- Fault-parallel flow (Prompt 21): `hanging_wall: {type: fault_parallel_flow}` on `ramp_flat`
  and `listric` faults, with `throw.slip` and `throw.sense: normal | reverse`, moves the hanging
  wall along flow lines parallel to the fault that meet on the bisectors of its bends (Egan,
  Buddin and Kane 1997), exactly and in closed form: contractional and extensional fault-bend
  folds. Ramp-flat profiles may now start with a flat and need not start at the reference
  point. Rule F18 refuses offsets where the bisectors of neighbouring bends cross; the crossing
  offset is reported. Worked example `fault_bend_fold.yaml`.
- Detachment folds (Prompt 21): primitive `type: detachment_fold` with `half_width`,
  `shortening` and `detachment: <horizon>`; each horizon's relief follows from the area balance
  (shortening times depth to the detachment, divided by the half-width) and is reported.
- Fold mechanics (Prompt 21): `mechanics: similar` gives every horizon the anchor's relief;
  `mechanics: parallel` makes the horizons below the anchor its offsets at their regional
  thickness, measured normal to the folded surface (cusps kept). `array: {count, step_along,
  step_across}` repeats a primitive en echelon. Rule F19 sets their admissible range. Worked
  example `detachment_folds.yaml`.
- Check C12: area conserved by the restoration of every fault with fold kinematics or a curved
  hanging wall (Jacobian determinant in the dip plane), tolerance `checks.tolerances.area`
  (default 0.01).

- Stage 2C specification (technical note STX-TN-07): stratigraphic series with conformable,
  erosional and onlap relations, stratal terminations, thickness trends and pinch-outs, growth
  strata, clinoforms; checks C13 and C14; rules F20 to F23. New interfaces, additive:
  `Structure::missing_section`, `repr::unit_thickness_grid`, `repr::termination_lines`.

### Changed
- Check C4 skips trishear and fault-parallel-flow faults, whose throw is not a profile; C12
  checks them instead.
- Check C4 measures faults whose slip has a strike component against the footwall horizon where
  the moving cut-off came from.
- Tests: the LeakSanitizer suppressions (font library caches) now reach every test; tests
  defined after them, such as `cli_implicit_section`, failed in the sanitiser job on exit.
- Golden package hashes: the checks report and summary gain the C10 lines and every figure's
  key numbers the new check count; the figure set gains the fault-plane maps. Exports are
  unchanged.
- Analytical tests for stage 2A (`test_v2_curved`, `test_v2_interaction`, `test_v2_movement`)
  and stage 2B (`test_v2_folding`, against independent trishear and fault-parallel-flow
  integrations); stage 2C (`test_v2_strat`), expected to fail until its features are built.

## [2026.1] - 2026-11-14 - First public release

Version numbers now put the release year first: `<year>.<major>`, with a third number only
for a patch release (`2026.1.1`). This release is the content of 1.3.1 under its public number;
the 1.x numbers below are the development history.

### Changed
- Version 2026.1 (was 1.3.1). The C API keeps its own version (1.1), which changes only with
  its binary interface.
- Container images are tagged with the release number taken from the git tag (`2026.1`) and
  `latest`.

## [1.3.1] - 2026-11-14 - Builds on Windows (MSVC) and Ubuntu 22.04; container image

### Fixed
- Windows (MSVC): plugin anchors have C linkage and are declared at file scope; before, MSVC
  gave block-scope declarations a different linkage from their definitions (unresolved externals).
- Linux: static libyaml and zlib are linked only when the whole archive is position-independent;
  otherwise the shared libraries (Ubuntu 22.04's archives cannot go into the C API shared library).
- Clang 14: no structured binding captured in a lambda (benchmark runner).
- Windows checkouts keep LF line endings (`.gitattributes`), so the golden hash files and
  examples compare byte for byte.

### Added
- Container image of the command-line program and C API (`docker/Dockerfile`), published to
  GitHub Packages as `ghcr.io/earthscitech/stratonix` when a release is published.
- Build and test presets for `linux-clang-release` and `linux-gcc-debug`.

## [1.3.0] - 2026-11-14 - Summary tables and faults that end inside the model

### Added
- Section-sheet summary table with the columns `quantity,value,basis` (was `key,value`): every
  key number with its basis; `basis` takes placeholders for the location and block letters.
- Twelve key-number kinds: `text`, `completion`, `vertical_exaggeration`, `fold_relief`,
  `limb_dip`, `fault_dip`, `fault_tip_units`, `fault_cutoff_position`, `fault_cutoff_depths`,
  `fault_overlap`, `throw_ratio`, `fault_gap`. Key numbers can be kept off the sheet
  (`on_sheet: false`); the heading is `key_numbers_title` (default "Key numbers").
- Tests S4 (every new kind against an independent calculation) and L4 (a layer-rigid fault
  shorter than the model, throw tapering to zero at its ends).

### Changed
- The sheet draws its key numbers right-aligned in the bottom-right corner under "Key numbers"
  (no file name); notes, if any, go under the legend.
- Layer-rigid faults that end inside the model: the build note now says whether the block offset
  dies out (strike profile tapering to zero) or stops on a line that is not a fault (no taper).

## [1.2.0] - 2026-11-14 - Section sheets and layer-rigid faults (application model reference)

### Added
- Hanging wall `layer_rigid`: every layer moves as a rigid slab by the throw where it meets the
  fault, faults acting oldest first on the displaced geometry. Restoration solves for the
  material's level; each fault tabulates where every level meets its plane. Build refuses throws
  that change faster than the layers thicken, and notes layer-rigid faults that end inside the
  model along strike.
- Tip rules can read their level on the unfaulted stratigraphy (`stratigraphy: unfaulted`).
- Well completions; unit display labels (UTF-8); units `top: model_top` / `base: model_base`
  that name the outer units (role, colour, label).
- Section sheets (`outputs.section_sheets`): main and detail panels with sharp per-block fills,
  unit labels, leader / fault-level / completion annotations, legend by role, computed key numbers
  (eight kinds) with a CSV, notes; DejaVu Sans text with built-in metrics.
- Rules W4 (completion order), O1 (sheet references), O2 (sheet windows); F1 refuses a
  within-unit tip in an outer unit. Worked example 6 (`layered_faults_sheet.yaml`); test_sheet.

## [1.1.1] - 2026-11-14 - Figure fixes found on a 23-horizon model (Prompt 12)

### Fixed
- Legends with many entries (more than fit one column above the key-numbers box) are set two to a
  row in a smaller font, so the key-numbers box is no longer pushed past the footer.
- Legend groups keep stratigraphic order (order of first appearance) instead of alphabetical order;
  units without a role stay last.
- 3D views: the see-through opacity of shallower horizons now falls with their number so that at
  least 20 % of the light from the faults and the deepest horizon passes through all of them
  (opacity = min(0.4, 1 - 0.2^(1/n))). With 22 shallower horizons the faults were invisible.
  Golden hashes of the two 3D figures in the synthetic package are updated for this change.

## [1.1.0] - 2026-11-13 - Opening v2: data-conditioned implicit builder (Prompt 14)

### Added
- Structure builder `implicit-rbf`: one stratigraphic scalar field interpolated by Hermite radial
  basis functions (cubic kernel, linear drift) from surface points, orientations (tangent
  constraints: direction only) and control-well tops, with a vertical stretch for anisotropy;
  surfaces extracted as level sets in the restored frame; faults from the v1 fault network, data
  restored through it until no datum moves by more than 1 cm. Faults may be given by sticks
  (plane fitted to the stick points). The model passes the model-interface suite and every v1
  command (build, section, export, package, suite) works on it unchanged.
- Misfit at every data point (`data_misfit.csv`; `data/data_misfit.csv` in exports) with the
  evidence class of its input; provenance of every input in export and package READMEs.
- Dense solvers in core: reference Gaussian elimination and a blocked, parallel LU with parity
  tests (identical for any thread count); 3 x 3 symmetric eigen-solver.
- Importers `csv_points` and `csv_orientations` (GemPy layout); description section `implicit`
  (rules I1-I3), fault geometry type `sticks` (rule F9); units `m2` and `km2` for closure-area
  targets (an area target could not be written before: no area unit existed).
- `examples/implicit_from_samples.yaml` with data sampled from the synthetic parametric example;
  tests J1-J5 and command-line tests on the implicit model. Method note STX-TN-04.

### Fixed
- Imported surfaces and isochores hold their edge value outside their grid
  (`Grid2D::sample_clamped`); before, fault-tip resolution near the domain edge met undefined
  values and stopped the build.

## [1.0.0] - 2026-11-13 - Release 1.0 (Prompt 11)

### Added
- C API 1.1 (additive to 1.0): horizon grids (faulted with gaps, or per block), export and package
  implemented; `stx_status_name`, `stx_model_fault_name`; every status code documented in
  `stratonix_c.h`. Shared library `stratonix_c` (SONAME 1). C example
  `capi/examples/stratonix_example.c` built and run as a test; `capi_test` covers build, evaluate,
  surfaces, export, package and error codes.
- Worked examples `examples/one_fault.yaml` and `examples/dome_blind_fault.yaml` (built as tests).
- Documentation: user guide, developer guide, roadmap v2-v5, API reference generated from the
  header comments (`tools/gen_api_reference.py`; `docs/Doxyfile` for Doxygen), `tools/build_docs.py`.
- Continuous integration: Windows packages install the vcpkg DLLs next to the program; a release
  job publishes the Windows installer and zip and the Linux archive with `SHA256SUMS.txt` on a tag.

### Release artefacts
- Linux: `stratonix-1.0.0-Linux-x86_64.tar.gz` and `.zip` with SHA-256, built here.
- Windows installer and portable zip: built by the CI release job (no Windows host here).

## [0.10.0] - 2026-11-13 - Delivery package (Prompt 10)

### Added
- `stratonix package <description> -o <folder> [--config package.yaml] [--zip]`: a self-contained
  delivery folder with README (model and version, DESIGN or interpreted status, coordinate system,
  vertical datum, z convention, units, grid spacing, horizons, units and faults, evidence classes in
  plain words, assumptions, known limitations, terms of use, contact), `parameters.csv` (every
  parameter with value, unit, evidence class, plain-words meaning and source), checks report,
  summary table, the model description, exports in the chosen formats, figures, a portable
  `manifest.json` with SHA-256 of every file, and optionally a deterministic zip.
- Inputs and their sharing status are listed first; an input with status `unknown` stops the
  command before anything is written. A failed check stops it too.
- Identifier scan: built-in patterns (register identifiers, work-breakdown codes, internal note
  identifiers, absolute local paths) plus configurable regular expressions are searched in every
  text file, every file name and the package folder name; internal labels are translated at the
  boundary first (configurable map). A package with a hit is renamed `<folder>_REJECTED_SCAN` with
  `SCAN_FAILED.txt` and no zip.
- Portable manifests (no compiler or platform), so the same package differs between machines only
  by the manifest timestamp; golden SHA-256 of every text file of the synthetic package
  (`tests/data/package_synthetic_sha256.txt`).
- Tests (`test_package`): the example packages and verifies; a contaminated description fails the
  scan and passes once translated; unknown sharing status is refused; two runs (1 and all threads)
  are identical apart from the manifest timestamp; configuration file parsing.
- `docs/package_config.md`.

### Fixed
- `stratonix verify manifest.json` with a relative path.

### Waiting for the model owner
- Terms-of-use text and the list of project codes to ban: packages are marked DRAFT until terms
  are supplied.

## [0.9.0] - 2026-11-13 - Interoperability (Prompt 9)

### Added
- Exporters: `vtk` (VTU surfaces with depth or displacement, wells as polylines, VTI voxel units
  and blocks, PVD collection, ParaView Python state with vertical exaggeration), `obj`, `stl`,
  `csv` (nodes with elevation, depth and block; XYZ; well tops), `geojson` (domain, fault polygons,
  wells; model coordinates, flagged as not WGS 84), `gempy` and `loopstructural` (interface points,
  orientations and fault data sampled from the model away from fault gaps).
- GOCAD TSurf now carries the coordinate-system name and a vertex property (`displacement` on
  faults, `depth` on horizons).
- Importers: `xyz_points`, `irap_lines`, `petrel_points` (grouped into sticks), `gocad_tsurf`,
  `well_tops` (grids: `irap_ascii`, `zmap`, `xyz` as before).
- Tests: `test_interop` round trips for every exporter and importer; GemPy and LoopStructural
  points lie on the model's horizons and faults; `tests/tools/check_neutral_set.py` checks the
  neutral files independently; `tests/tools/crosscheck_implicit.py` rebuilds the horizons from
  the LoopStructural set (LoopStructural when installed, otherwise a SciPy RBF substitute) and
  reports depth differences near and away from faults (see STX-TN-03).

### Not done here
- RESQML 2.0.1 (FESAPI and HDF5 cannot be built without network): deferred.
- ParaView, VTK and LoopStructural are not installable here: files checked with own parsers,
  cross-check run with the substitute interpolator.

## [0.8.0] - 2026-11-13 - Petrel export set (Prompt 8)

### Added
- Exporter plugins (`stratonix plugins` lists them): `irap_ascii` (faulted horizons with fault
  gaps; per-block surfaces for horizons a fault offsets), `zmap` (faulted horizons),
  `fault_polygons_irap` (per horizon, IRAP classic lines), `fault_sticks_irap`, `petrel_points`
  (stick vertices with Fault and Stick attributes), `gocad_tsurf` (fault planes clipped to the tip
  line; faulted horizons), `well_heads` (and paths of deviated wells), `well_tops`. The set name
  `petrel` selects all of them.
- Writers and readers for IRAP classic lines, GOCAD TSurf, Petrel points with attributes, well
  heads and well tops (IRAP and ZMAP+ grids already had both).
- `stratonix export <description> -o <folder> --formats petrel [--grid-spacing m]`; `build` writes
  the description's export list into `exports/` when every check passes. Every export folder has a
  `README.txt`: model, evidence basis, coordinate system, vertical datum, units, z convention,
  grid geometry, undefined values, each file with its Petrel import option, block codes.
- Tests: round trips against the model for every format (grids within 0.5e-4 m with identical gaps,
  sticks and TSurf vertices on the fault within 0.01 m, tops on their horizons); determinism under
  a decimal-comma locale; LF only, ASCII only; golden SHA-256 of the whole synthetic set
  (`tests/data/petrel_synthetic_sha256.txt`, version string normalised) so any platform that
  writes different bytes fails; `tests/tools/check_petrel_set.py` cross-checks the formats
  independently (IRAP against ZMAP+, TSurf against IRAP, sticks against planes, every grid node
  inside a fault polygon is a gap and every gap lies at a polygon, tops against grids).
- `docs/format_samples/petrel_synthetic_250m`: a small export to inspect by eye.
- Petrel import checklist: STX-TN-02 (project outputs).

### Not done here
- xtgeo is not installable here (no network); the independent check uses own Python parsers.
- Byte identity with Windows is enforced by the golden-hash test but has not been run on Windows.
- Import into Petrel not yet run (checklist ready; Prompt 13).

## [0.7.0] - 2026-11-13 - Figures (Prompt 7)

### Added
- Figures layer (`stratonix_figures`, level 5): a display list drawn to SVG by Stratonix itself
  (text kept as text) and rasterised with Cairo to PNG (300 dpi, pHYs) and TIFF (300 dpi, LZW,
  RGB, no alpha, own encoder); PPTX written by Stratonix (one 16:9 slide per figure with title,
  picture, caption and speaker notes) plus a deck of the whole set.
- Sections along any polyline (unit fills per fault block without seams, fault traces between
  tips, wells and tops within 300 m, optional detail panel, legend grouped by role, key-numbers box
  from the summary table, marked bends); maps per horizon (depth shading and contours, fault gaps
  drawn where the model has them, cut-off polygons, wells, closure outline, crest and spill point);
  3D views at two camera presets with vertical exaggeration (software renderer, see-through
  shallower horizons); fence diagram through the wells with tops and a depth scale.
- Every figure carries a footer with coordinate system, vertical datum and evidence basis
  (DESIGN models are also stamped "DESIGN MODEL - not an interpretation"); names follow
  `STX_<model>_v<version>_<kind>_<name>`.
- `stratonix build` writes `figures/` and `summary_table.csv`; `stratonix section` draws one
  section by name or along `--line x1,y1;x2,y2;...`.
- Tests: encoders (PNG pHYs, TIFF LZW round trip, PPTX parts), figure set of the example;
  `tests/tools/check_figures.py` checks every figure independently (Pillow, tifffile,
  python-pptx); `cli_build_example` requires the example to build with every check passing.

### Fixed
- Fault restoration: when the throw sits on the profile's plateau the root of the restoration
  equation is at the end of the bracket, and rounding could give it the wrong sign there, so the
  restoration fell back to the throw at the present position. Ends within 1e-9 of zero are now
  taken as the root. This made check C4 fail at five strike positions of F2 in the example.

### Added (validation)
- Build refuses a fault whose throw changes by 1 m or more per metre of slip (the wall would
  overlap itself and the model would not be unique); a note is written above 0.8.

### Not done here
- VTK is not available in the build environment; 3D views use Stratonix's own renderer.
- Windows rendering not run (no Windows host); SVG is platform independent, rasters depend on the
  fonts Cairo finds.

## [0.6.0] - 2026-11-13 - Performance: optimised kernels and parity (Prompt 6)

### Added
- Optimised kernels behind the reference ones, selected by `set_optimised_kernels()`
  (`--reference-kernels`, `bench --mode reference|optimised|both`): structure terms shared across
  horizons; per-block and faulted horizon grids by Illinois iteration seeded from a coarse grid;
  cut-off lines by Illinois on the bisection bracket; column-wise voxel grids.
- Parity tests K1-K5 (`test_parity`): point evaluation bitwise equal and identical for 1, 3 and
  8 threads; grids within 1e-6 m with identical gap masks; cut-off lines within 1e-4 m; voxel
  grids identical. Edge cases: points on fault planes, at tips, on horizons, outside the domain.
- `docs/performance_and_threading.md`: kernels, determinism and thread-safety rules.

### Measured (synthetic production model, 40 x 40 km, ten faults, 2 threads, GCC 13, Linux)
| Spacing | Reference full build + checks | Optimised | Speed-up | Peak memory |
|---|---|---|---|---|
| 100 m | 165.3 s | 7.3 s | 22.6x | 110 MB |
| 50 m | 465.6 s | 14.8 s | 31.5x | 144 MB |
| 25 m | not run (reference) | 32.9 s | - | 414 MB |

Before this release (0.5.0, optimised path of that time) the 50 m build took 32.6 s.
Windows timings: pending the CI runner (no Windows host in the build environment).

## [0.5.0] - 2026-11-13 - Evaluation and representations (Prompt 5)

### Added
- Horizon meshes from grids (no triangles across gaps); fault meshes checked for degenerate and
  non-manifold triangles; fault sticks.
- Tests: faulted and per-block surfaces lie on their horizons in the blocks they report; voxel
  grid agrees with per-block surfaces; well tops on vertical and deviated wells; sections are
  slices of the 3D model; contours; mesh quality.
- `stratonix bench`: full build and kernel timings with peak memory on the synthetic
  40 x 40 km, ten-fault production model (`--write-model` writes its description).
- Bracketed root finding (Illinois) shared by the engines; tip tables resolved in parallel.

### Changed
- Cut-off lines: bracketing segments must lie in one block and on the intended wall, so a cut-off
  is never interpolated across another fault.
- Check C4 measures throw in the fault's own frame (younger faults undone), where the fault is a
  plane, so faults cut by younger faults are measured correctly.

### Measured
- 40 x 40 km at 50 m, ten faults, 2 threads (GCC 13, Linux): full build and checks 32.6 s,
  peak memory 139 MB; faulted grids of six horizons (801 x 801) 8.3 s.

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
