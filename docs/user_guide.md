---
title: Stratonix user guide
doc_id: STX-UG
version: 1.3
date: 2026-11-14
---

# Stratonix user guide

Stratonix builds parametric structural models (horizons and faults in 2D sections and 3D), checks
them, draws them and exports them for other software, Petrel first. Version 1 builds **design**
models: every number in a model is stated in its description with a unit, an evidence class and a
source, and Stratonix refuses a description where one is missing.

# 1. Installing

## 1.1 Windows

Two forms are produced by the release build on Windows (MSVC):

- **Installer** `stratonix-1.0.0-Windows-AMD64.exe` (NSIS). Run it, accept the licence, choose a
  folder (default `C:\Program Files\Stratonix`) and choose to add Stratonix to the PATH if you
  want to run `stratonix` from any command prompt. Uninstall from Windows Settings, Apps.
- **Portable zip** `stratonix-1.0.0-Windows-AMD64.zip`. Unzip anywhere; run `bin\stratonix.exe`.
  Nothing is written outside that folder.

Check the download against `SHA256SUMS.txt` before installing:
`certutil -hashfile stratonix-1.0.0-Windows-AMD64.exe SHA256`.

Then open a command prompt and run `stratonix version`.

## 1.2 Linux

Unpack `stratonix-1.0.0-Linux-x86_64.tar.gz` anywhere (`tar xzf ...`) and run `bin/stratonix version`.
Check it first with `sha256sum -c SHA256SUMS.txt`. The program needs Cairo (`libcairo2`) for
raster figures; without Cairo it still builds, checks and exports, and writes SVG figures only.

## 1.3 From source

```
cmake --preset linux-gcc-release        # or windows-msvc-release
cmake --build --preset linux-gcc-release
ctest --preset linux-gcc-release
```

Dependencies (libyaml, zlib, Cairo) come from vcpkg on Windows and from the system on Linux; see
the developer guide.

# 2. Writing a model description

A description is a YAML file (JSON also works). `stratonix schema` prints the full JSON Schema.
Every numeric parameter is written in full or compact form:

```yaml
throw_max: {value: 60, unit: m, evidence: assumed, source: "design case 2"}
throw_max: "60 m | assumed | design case 2"
```

Evidence classes: `measured`, `interpreted`, `calibrated`, `upscaled`, `literature_derived`,
`assumed`, `user_instruction`. Lists are written `"[3000, 3000] m | ..."`; a flag is
`"true flag | ..."`. Depths in descriptions are positive down; exported files use elevation.

The sections of a description:

| Section | What it holds |
|---|---|
| `model` | name, version, `status: DESIGN` or `INTERPRETED`, title, contact |
| `crs` | coordinate system name (and EPSG code or WKT when there is one), vertical datum, length unit (m) |
| `domain` | origin, size, depth range and grid spacing of the model box |
| `inputs` | files the model reads (imported surfaces, isochores), each with its **sharing status** (internal, shareable, restricted or unknown) and source |
| `wells` | name, role (`control`, `check`, `planned`), head, optional datum, total depth or survey stations, tops |
| `stratigraphy` | reference (a control well, or a point with depth and thicknesses), horizons top to base, units between them; a horizon may come from an imported grid or an isochore |
| `structure` | regional surface (flat or planar), flat base if wanted, fold primitives (anticline, dome, syncline, basin) with relief per horizon, well-correction range |
| `faults` | per fault: age (order), planar geometry (point, depth, strike, dip, or trace and dip), extent (upper and lower tip rules, strike length, tip-line shape), throw (maximum, profiles down dip and along strike), slip mode, hanging-wall field, abutting |
| `frame` | which block stays fixed: a fixed well or point |
| `checks` | tolerances and closure targets |
| `outputs` | sections, vertical exaggeration, voxel spacing, stick spacing, export formats |

Run `stratonix validate model.yaml` while you write: it checks the schema and the semantic rules
(evidence on every parameter, units, horizons in order, tips in the right order, wells inside the
domain and so on) and names the field at fault.

## 2.1 Worked example 1: one fault

`examples/one_fault.yaml`: three horizons dipping 2 degrees east, cut by one normal fault that
reaches above the top horizon and below the base. Points to note:

- `regional: {type: planar, dip: ..., dip_azimuth: ...}` tilts the whole stack.
- The fault's tips are placed relative to the stratigraphy: `upper_tip: {rule: horizon, horizon:
  Top_A, offset: "-50 m | ..."}` puts the tip 50 m above Top_A (negative offset is upwards).
- Because the fault cuts the top and base horizons, those horizons carry `intended_offset: true`.
  Without it, check C6 fails ("top horizon unoffset") — the build stops and tells you which flag
  to set. This is deliberate: an offset top or base is usually a mistake.

## 2.2 Worked example 2: two faults with different extents

`examples/synthetic_fold_two_faults.yaml`: an asymmetric anticline with a blind fault F1 (tips
within Unit_B and Unit_D, elliptical tip line, 40 m) and a younger fault F2 reaching H1 (60 m,
rectangular tip line, displacement decaying over 1500 m in the hanging wall). Points to note:

- Ages order the faults; younger faults are restored first, and a younger fault's tips refer to
  the stratigraphy as displaced by older faults.
- `hanging_wall: {type: decay, influence_distance: ...}` lets the displacement die out away from
  the fault; `rigid` moves the whole block.
- `frame: {fixed: {well: W1}}` keeps the block of the control well fixed, so the well tops are
  honoured exactly.
- `checks.closure` reports crest, spill point, closure and area on H2 with faults sealing (an
  assumption, labelled `assumed`).

## 2.3 Worked example 3: a dome with a blind fault

`examples/dome_blind_fault.yaml`: a dome of 150 m relief on H1, decreasing to zero at the flat
base, with a blind fault on its east flank. Points to note:

- The throw profile down dip is `sine`. An `elliptical` dip profile is refused: its throw changes
  infinitely fast at the tips, so the hanging wall would overlap itself (the build names the fault
  and the remedy). Any profile is refused when the throw changes by 1 m or more per metre of slip.
- `faults_seal: "false flag | ..."`: with the fault open, the closure on H2 spills at the fault
  (44.3 m of closure, 5.72 km2) instead of at the dome's own rim.

## 2.4 Data-conditioned models (implicit builder, from 1.1)

`examples/implicit_from_samples.yaml` builds the same anticline from data instead of parameters:
surface points and orientations (CSV with X, Y, Z, surface and, for orientations, G_x, G_y, G_z),
the control-well tops, and fault sticks (`geometry: {type: sticks, input: ..., input_format:
irap_lines}`; tips and throw as in v1). The `implicit` section lists the data inputs and the
interpolation settings (nugget, vertical stretch, surface spacing, maximum restoration passes),
each with unit, evidence and source. The build writes `data_misfit.csv`: the misfit at every data
point with the evidence class of its input; exports and packages carry it and the provenance of
every input. Data are never adjusted to improve the fit. Method and accuracy: the note on the
implicit builder.

## 2.5 Presentation sheets, layer-rigid faults and other 1.2 additions

`examples/layered_faults_sheet.yaml` shows the features added in 1.2:

- **Layer-rigid hanging walls** (`hanging_wall: {type: layer_rigid}`). Every layer of the moving
  block moves as a rigid slab by the throw where that layer meets the fault, so a horizon keeps
  one offset across its whole block, also far from the fault. With `rigid` the throw is read where
  the point's own elevation meets the fault. Faults act in age order (oldest first) on the
  geometry the older ones displaced. Use it for all faults of a model or for none, and make the
  faults longer than the model along strike: a block offset that ends at a fault's strike end
  inside the model ends along a line that is not a fault (the build notes it).
- **Tips on the unfaulted stratigraphy** (`stratigraphy: unfaulted` in a tip rule). By default a
  younger fault's tip is placed on the stratigraphy as displaced by older faults; `unfaulted`
  places it on the stratigraphy before any faulting.
- **Completions** on wells (`completions: [{name, top, base}]`, depths positive down) are drawn
  on sections and sheets.
- **Unit labels** (`label`, UTF-8, for example "Våle Fm") are used in legends and on panels;
  `top: model_top` or `base: model_base` names the unit above the first or below the last horizon
  so it gets a role, a colour and a label.
- **Section sheets** (`outputs.section_sheets`): one page per section with a main panel at a stated
  vertical exaggeration, an optional detail panel (window and exaggeration), unit labels,
  annotations (`leader` to a horizon, `horizon_at_fault` with the computed depth, `completion`),
  a legend by role, the key numbers Stratonix computes from the model and your notes. Titles take
  placeholders `{units}`, `{surfaces}`, `{faults}`, `{faults_word}`, `{ve}`, and the detail title
  `{s0_km}`, `{s1_km}`; annotation texts take `{value:N}`, `{top:N}`, `{base:N}`. Key-number kinds:
  `horizon_at_well`, `horizon_highest`, `block_offset`, `fault_tips`, `fault_throw`,
  `well_fault_distance` (on the unfaulted unit), `block_highest_at_fault`, `relief_to_end`. The
  sheet is written in the four figure formats beside a CSV of its key numbers; text is set in
  DejaVu Sans (if the font is not installed the system substitutes another sans-serif font).

# 3. Running

| Command | What it does |
|---|---|
| `stratonix validate model.yaml` | schema and semantic rules only |
| `stratonix build model.yaml -o out/run1` | build, run every check, write summary, checks report, summary table, figures and the exports listed in the description, and `manifest.json` |
| `stratonix section model.yaml -o out/s1 --line "x1,y1;x2,y2;..."` | one section along any polyline (or `--name` of a section in the description) |
| `stratonix export model.yaml -o out/e1 --formats petrel,vtk` | exports only |
| `stratonix package model.yaml -o out/delivery --config package.yaml --zip` | delivery package (section 6) |
| `stratonix verify out/run1/manifest.json` | checks every file against its SHA-256 |
| `stratonix plugins` | lists structure builders, kinematics, hanging-wall fields, exporters, importers, checks |
| `stratonix bench --spacing 100,50` | timings of the kernels |

Options common to building commands: `--threads N` (results never depend on it),
`--reference-kernels` (plain reference implementations, slower; for checking), `--no-figures`.
Output folders must be new or empty: Stratonix never overwrites an issued output.

Exit codes: 0 success; 2 a check failed (nothing exported); 3 invalid description or input
sharing status unknown; 4 a check failed in packaging; 5 function not available; 6 package
identifier scan failed; 1 any other error. The message names the cause.

# 4. Reading the checks report

`checks_report.txt` lists every check with its value, tolerance (and where the tolerance came
from: your description or a Stratonix default) and PASS, FAIL or INFO. A single FAIL stops the
exports and figures.

| Check | Fails when |
|---|---|
| C1 thickness | a unit has less than the minimum thickness somewhere in a fault block |
| C2 crossing | two horizons cross |
| C3 offset outside tips | a horizon is offset beyond a fault's tip line |
| C4 throw | the throw measured at a horizon differs from the fault's profile |
| C5 wells | a control-well top is missed by more than the tolerance (check wells are reported only) |
| C6 top/base unoffset | the top or base horizon is offset and `intended_offset` is not set |
| C7 faults in domain | a fault plane leaves the domain or its mesh is degenerate |
| C8 fault intersections | two faults intersect without an abutting rule |
| C9 closure | INFO: crest, spill point, closure, area; FAIL only against targets you set |

Targets are reported, never used to tune the model.

# 5. Exporting

Choose formats with `--formats` or the description's `outputs.exports`. `petrel` selects the
Petrel set. Every export folder has a `README.txt` giving coordinate system, vertical datum,
units, z convention (elevation, negative below the datum), grid geometry, undefined values and,
for each file, the Petrel import option.

| Id | Files |
|---|---|
| `irap_ascii` | faulted horizons with fault gaps (undefined 9999900.0); per-block surfaces for horizons a fault offsets |
| `zmap` | faulted horizons (null -99999.0) |
| `fault_polygons_irap` | fault polygons per horizon (IRAP classic lines) |
| `fault_sticks_irap`, `petrel_points` | fault sticks as IRAP lines, and as points with Fault and Stick attributes |
| `gocad_tsurf` | fault planes clipped to the tip line, faulted horizons; vertex properties |
| `well_heads`, `well_tops` | Petrel well heads (and paths of deviated wells), well tops |
| `vtk` | VTU surfaces and wells, VTI voxel grid, PVD collection, ParaView Python state |
| `obj`, `stl`, `csv`, `geojson` | neutral surfaces, nodes with depth and block, outlines |
| `gempy`, `loopstructural` | input sets sampled from the model for those modellers |

Load into Petrel with the Petrel import checklist: it gives the dialog option for each file and
what to look at to confirm a correct load.

# 6. Packaging for delivery

`stratonix package` writes a folder for an external recipient: README (status, coordinate system,
datum, z convention, units, grid spacing, horizons, units and faults, evidence classes in plain
words, assumptions, limitations, terms of use, contact), `parameters.csv`, checks report, summary,
the description, exports, figures and `manifest.json`; optionally a zip. Before building it lists
every input's sharing status and stops on `unknown`. After writing it translates internal labels
and scans every text file and file name for banned patterns (register identifiers, internal note
numbers, local paths, and your own list); a hit rejects the package. Without approved terms of use
the package is marked DRAFT. Settings: `docs/package_config.md`.

# 7. Known limitations of the kinematics

- Faults are planar. Listric or curved faults are not in v1.
- Throw follows the stated profiles; there is no strain-compatible restoration, no fault-zone width
  and no drag.
- Dip-slip moves the hanging wall down the dip; vertical offset moves it straight down. They differ
  by up to a few metres on dipping horizons (3.5 m on the synthetic example).
- Hanging-wall fields are rigid, layer-rigid (1.2), linear decay to a distance, or linear decay to
  the next fault. Layer-rigid fields need every fault of the model to use them.
- A fault's throw must change by less than 1 m per metre of slip; elliptical dip profiles are
  refused for that reason.
- Reverse faults are evaluated, but where a horizon occurs twice the upper occurrence is kept on
  grids (reported as overlap nodes).
- At most 64 faults per model.
