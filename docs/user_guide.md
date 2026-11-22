---
title: Stratonix user guide
doc_id: STX-UG
version: 3.0
date: 2026-11-15
---

# Stratonix user guide

Stratonix builds parametric structural models (horizons and faults in 2D sections and 3D), checks
them, draws them and exports them for other software, Petrel first. It builds **design** models:
every number in a model is stated in its description with a unit, an evidence class and a source,
and Stratonix refuses a description where one is missing.

Release 2026.2 (version 2) adds:

- curved faults and fault systems;
- folding mechanisms;
- stratigraphic architecture;
- salt and compaction;
- deformation history and restoration;
- faults as rock with the stress field.

Sections 2.6 to 2.11 describe them. A description written for 2026.1 builds unchanged.

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
Units: m, km, ft, m2, km2, 1/m, 1/km, deg, rad, fraction, ratio, percent, ordinal, flag, Ma, Pa,
kPa, MPa, bar, psi, MPa/km, kPa/m, bar/m, psi/ft, mD, D; from 2026.4 also kg/m3, g/cm3, mGal,
uGal, m/s2, nT, T and SI (susceptibility).

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
| `bodies` | (2026.4) salt sheets, sills, saucer sills, dykes and sand injectites |
| `geophysics` | (2026.4) unit densities and susceptibilities, background values, inducing field, gravity and magnetic stations |

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
  geometry the older ones displaced. Use it for all faults of a model or for none. A fault may end
  inside the model along strike if its strike profile tapers the throw to zero at the ends (for
  example `profile_strike: {type: flat_topped, taper: 0.5}`, a cosine bell over the whole length):
  the block offset then dies out smoothly towards the end. Without a taper (`taper: 0`) the offset
  would stop along a line that is not a fault; the build notes either case.
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
  `{s0_km}`, `{s1_km}`; annotation texts take `{value:N}`, `{top:N}`, `{base:N}`.

  **Key numbers and the summary table (1.3).** Each entry of `key_numbers` gives a `kind`, a
  `label` (the quantity) and, optionally, a `basis` (where the number comes from). The sheet's
  summary table `<model>_v<version>_<sheet>_summary.csv` has the columns `quantity,value,basis`
  and lists every entry in order; the sheet draws the entries with `on_sheet` true (the default)
  in its bottom-right corner under the heading `key_numbers_title` (default "Key numbers").
  `decimals` sets the rounding. Kinds:

  | Kind | Value | Needs |
  |---|---|---|
  | `text` | the `value` as given | `value` |
  | `completion` | top-base of a completion (m) | `well`, `completion` |
  | `vertical_exaggeration` | the main panel's exaggeration | - |
  | `fold_relief` | unfaulted horizon: crest elevation minus its elevation at the section `end` | `horizon`, `end` |
  | `limb_dip` | unfaulted horizon: steepest dip between the crest and the section `end` (deg) | `horizon`, `end` |
  | `horizon_at_well` | depth where a well meets a horizon | `horizon`, `well` |
  | `horizon_highest` | shallowest depth of a horizon along the section, after faulting | `horizon` |
  | `block_offset` | how far the block between two faults is down-thrown relative to a well's block | `horizon`, `faults` (2), `well` |
  | `fault_tips` | upper / lower tip depths on the section | `fault` |
  | `fault_tip_units` | labels of the units holding the upper / lower tip | `fault` |
  | `fault_dip` | dip and dip direction, e.g. "60 deg towards north" | `fault` |
  | `fault_throw` | vertical throw at a horizon | `fault`, `horizon` |
  | `fault_cutoff_position` | where a horizon meets the fault on one `side`, as `coordinate` x, y or distance | `fault`, `horizon` |
  | `fault_cutoff_depths` | footwall / hanging-wall cut-off depths of a horizon | `fault`, `horizon` |
  | `fault_overlap` | 1 - throw / unit thickness at the fault | `fault`, `horizon`, `unit` |
  | `throw_ratio` | throw / thickness between horizons `from` and `to` at the fault | `fault`, `horizon`, `from`, `to` |
  | `fault_gap` | horizontal gap from the first fault's lower tip to the second fault | `faults` (2) |
  | `well_fault_distance` | horizontal distances from a well to faults at a level in a unit (`fraction` 0 top, 0.5 mid, 1 base) on the unfaulted stratigraphy | `well`, `faults`, `unit` |
  | `block_highest_at_fault` | a well's block: horizon depth where it meets a fault | `well`, `fault`, `horizon` |
  | `relief_to_end` | a well's block: horizon relief from that fault to the section `end` | `well`, `fault`, `horizon` |

  Thickness at a fault (`fault_overlap`, `throw_ratio`) is the vertical distance between the two
  horizons' cut-offs on the fault's footwall: the interval the fault cuts. Blocks are lettered A,
  B, ... in the order they hold the key's horizon along the section. A `basis` may use
  placeholders: `{x}`, `{y}`, `{s}` (location of the value, with `:N` decimals), `{block}`,
  `{fw_block}`, `{hw_block}`. The sheet is written in the four figure formats beside its summary
  table; text is set in DejaVu Sans (if the font is not installed the system substitutes another
  sans-serif font). Notes, if any, are drawn under the legend.

## 2.6 Curved faults and fault systems (2026.2)

- **Geometry.** `geometry.type` is `planar`, `listric` (`dip` at the reference depth flattening to
  `detachment_depth`), `ramp_flat` (flats and ramps from a `profile`) or `sticks` (a curved
  surface interpolated through fault sticks from `input`). `trace` gives a curved surface trace.
- **Slip.** `slip_mode` is `dip_slip`, `vertical_offset`, `oblique` (with `rake`) or
  `strike_slip`. `throw.sense` is `normal`, `reverse`, `left_lateral` or `right_lateral`.
- **Hanging-wall fields.** Beside `rigid`, `decay`, `linear_between_faults` and `layer_rigid`:
  - `vertical_shear` and `inclined_shear` (`shear_angle`, `shear_sense`), for rollover on
    listric faults;
  - `fault_parallel_flow`, for ramp-flat faults.
- **Fault systems.**
  - `abuts: {fault, side, transfer}` ends a fault against another; with `transfer` the throw
    passes across the branch line (splays).
  - `relay: {fault, breached}` hands throw to an overlapping fault across a relay ramp.
  - `drag: {type, fraction, width}` bends the beds next to the fault.
- **Checks.** C10 (displacement-length scaling, a warning) and C11 (throw continuous across
  branch lines and relays).

Examples: `listric_rollover`, `fault_system_relay_splay`, `thrust_repeated_section`,
`strike_slip_offset`, `oblique_fault`.

## 2.7 Folding mechanisms (2026.2)

- **Trishear fault-propagation folds** (`fold: {type: trishear, apical_angle,
  propagation_slip_ratio}` on a fault). The beds ahead of the tip fold in a triangular zone.
- **Fault-bend folds** come from `ramp_flat` faults with `fault_parallel_flow` hanging walls.
- **Detachment folds** (`type: detachment_fold`, with `shortening` above a `detachment` horizon).
  Their relief comes from area balance.
- **Fold primitives.**
  - `mechanics`: `parallel` (constant bed thickness) or `similar`.
  - `array: {count, step_along, step_across}`: en echelon fold trains.
  - `growth: {onset, end}`: folds that grew during deposition.
- **Check.** C12 (area conserved by the fold kinematics).

Examples: `trishear_fold`, `fault_bend_fold`, `detachment_folds`.

## 2.8 Stratigraphic architecture (2026.2)

- **Series** (`stratigraphy.series`) group horizons. Each series' `base` is `conformable`,
  `erosional` (truncates older series: unconformities) or `onlap` (younger beds end against it).
  `tilt: {dip, azimuth}` tilts a series before the younger ones are laid on it.
- **Clinoforms** (`clinoforms: {profile, height, width, direction, first_break, progradation}`;
  since 2026.4 also `profile: tangential` and `trajectory`, below)
  build prograding sigmoidal or oblique surfaces that downlap and toplap.
- **Isochores** (`isochore: {from, thickness, trend}`) place a horizon by its thickness from
  another, with a `linear` or `radial` trend; the thickness never goes below zero (pinch-out).
- **Growth strata.** `growth: {onset, end}` on a fault or a fold makes it move during the
  deposition of the horizons between them, so the units thicken into the hanging wall.
- **Outputs.**
  - Stratal termination lines (truncation, onlap, downlap, toplap, pinch-out) on maps.
  - Isopach maps of terminating units, and of those in `outputs.isopachs`.
- **Checks.** C13 (no negative thickness after the series relations) and C14 (truncation
  consistent with age order).

Examples: `unconformity_onlap_pinchout`, `growth_fault`, `prograding_wedge`.

## 2.9 Salt and burial (2026.2)

- **Salt** (`structure.salt: {top, base, source_thickness, bodies}`). Bodies are `pillow`,
  `diapir` or `wall`, each with a withdrawal rim that conserves the salt volume. Overburden over
  a body is raised; strata older than the salt top's rise are clipped against it. Bodies may
  grow (`growth`). Rule F24 refuses bodies the source layer cannot feed.
- **Compaction** (`stratigraphy.compaction`).
  - Laws per unit: `exponential` (`surface_porosity`, `decay`) or `table`.
  - `geometry: present` treats the stated depths as present-day.
  - `geometry: depositional` stacks the units from their solid thickness, so drape and
    differential compaction appear over highs.
- **Outputs.**
  - Burial-history curves at the wells.
  - Porosity grids per compacting unit (`exports/porosity/`).
- **Checks.** C15 (salt volume conserved) and C16 (solid thickness conserved through compaction).

Examples: `salt_pillow_withdrawal`, `salt_diapir_crestal_faults`, `compaction_drape`.

## 2.10 Deformation history and restoration (2026.2)

- **Events** (top-level `events`) put the model on a timeline in Ma. Types:
  - `deposition` (a horizon at an `age`);
  - `fault` (`from`, `to`, `throw`, `sense`);
  - `fold` (a primitive);
  - `salt` (a body).

  A fault may have several phases: they sum into one signed throw, so inversion (normal, then
  reverse) works. Rule F26 checks the timeline.
- **The model at any event.** `stratonix state model.yaml --age 45 -o out/45Ma` builds and
  packages the model as it was at that time.
- **Restoration report.** `restoration/` in every build package of a model with events holds:
  - a section and a map per event;
  - the sequence as one deck;
  - a flexural-slip unfolded section;
  - `balance.csv` and `states_checks.csv`.
- **Checks.** C17 (every state builds and the present state reproduces the model) and C18
  (line-length and area balance, a warning).

Example: `rift_inversion_history`.

## 2.11 Faults as rock, stress and curvature (2026.2)

A `fault_rock` block, with `clay_fraction` given on the units:

```yaml
fault_rock:
  shale_threshold: "0.5 fraction | assumed | clay fraction of a smearing shale"
  thickness_displacement_ratio: "0.01 ratio | literature_derived | Hull (1988)"
  seal_calibration: {type: sgr_power, a: "27 ratio | literature_derived | ...", c: "0.5 ratio | literature_derived | ..."}
```

- **Juxtaposition and fault rock.** With the block above, every fault gets:
  - an Allan diagram;
  - shale gouge ratio (SGR), clay smear potential (CSP) and fault-zone thickness maps;
  - with a calibration, seal capacity (bar). Without one, no seal capacity is reported.

  An attribute whose inputs are missing is reported as missing, with the reason; it is never
  defaulted.
- **Stress** (`stress`). Give either a `regime` (`sv_gradient`, `shmax_ratio`, `shmin_ratio`,
  `shmax_azimuth`) or a `tensor_gradient` (`[xx, yy, zz, xy, xz, yz]` in east, north, up), not
  both (rule F28). Add `pore_pressure_gradient`, `friction` and, optionally, `cohesion`. Units
  include MPa/km, kPa/m, bar/m and psi/ft. Every fault gets:
  - slip tendency;
  - dilation tendency;
  - the pore-pressure increase to reactivation (MPa).
- **Curvature** (`outputs.curvature: {horizons, fracture_threshold}`). Curvature maps with
  fracture-trend ticks, and nine curvature grids per horizon (`exports/curvature/`). Curvature is
  positive where the surface is convex upwards.
- **Summary table.** It lists every fault-rock and stress input with its evidence and source.
  Below them, each fault's attribute ranges and its self-juxtaposed and sand-on-sand shares.
- **Checks.** C19 (juxtaposition consistent with the restored stratigraphy) and C20 (tendencies
  within their bounds).

Examples: `fault_system_relay_splay`, `trishear_fold`, `salt_diapir_crestal_faults`,
`detachment_folds`.

## 2.12 Models conditioned to data (2026.3)

The implicit builder (2.4) now takes large data sets, series and faults from data, and the
parametric builder can fit its parameters to data.

- **Large data.** `implicit.solver: {type: auto | dense | partition_of_unity, points_per_patch,
  overlap}`. `auto` uses the dense solver up to 3000 interface data and a partition of unity
  above, which keeps exact interpolation: 100 000 points build in about 12 s.
- **Series.** With `stratigraphy.series`, one field per series, combined by their relations
  (erosional, onlap). A horizon that thins out within its series is listed in
  `implicit.isochore_horizons`: it is built from the horizon above minus a thickness from its own
  picks.
- **Faults from data.**
  - `throw: {estimate: {band, station_spacing}}` estimates throw, length, centre and vertical tips
    from horizon cut-offs (labelled calibrated).
  - `geometry: {type: sticks, fit: profile}` fits a listric or ramp-flat profile to the sticks.
- **Uncertainty of data.** Give `uncertainty` on an `implicit.data` entry, `sigma` (m) and
  `sigma_deg` columns in a CSV, or `wells[].top_uncertainty`. A point with an uncertainty is
  approached within it rather than honoured exactly. `implicit.data_influence: true` adds
  leave-one-out errors and the ten most influential data.
- **Hybrid conditioning.** A top-level `conditioning` block fits any scalar parameter of a
  parametric description to points, orientations and well tops with stated uncertainties:

```yaml
conditioning:
  data: [{input: picks, kind: points, horizon: H2, uncertainty: "3 m | interpreted | pick accuracy"}]
  parameters:
    - {path: faults.F1.throw.max, min: "60 m | assumed | ...", max: "200 m | assumed | ..."}
```

  Fitted values are labelled calibrated, with their uncertainty (`conditioning.csv`) and the data
  that control them (`influence.csv`). A fit that does not converge inside its range, or that no
  datum controls, is refused with the reason. The least-squares fit updates its Jacobian by
  Broyden's formula, so an iteration costs one model build (2026.4).
- **Bayesian conditioning (2026.4).** `conditioning.method: bayesian` samples the posterior of
  the same parameters: the prior is the parameter's distribution under `uncertainty.parameters`,
  or uniform on [min, max]; the likelihood is Gaussian in the stated uncertainties.

```yaml
conditioning:
  method: bayesian
  sampler: {walkers: "16 ordinal", steps: "2000 ordinal", burn_in: "500 ordinal", thin: "10 ordinal", seed: "1 ordinal"}
```

  The build writes `posterior_samples.csv` and `posterior_summary.csv` (prior, least-squares
  value with its Laplace sigma, mean, 95 % interval, R-hat, effective sample size), trace and
  corner figures, and runs check C28. `stratonix ensemble` on the same description builds one
  member per retained sample (or `--members N` of them), so statistics, sensitivity and
  probability maps describe the posterior.
- **Checks.** C21 (data misfit within the stated uncertainty) and C22 (the conditioning
  converged inside its ranges).

Examples: `implicit_unconformity_from_samples`, `implicit_relay_splay_from_samples`,
`implicit_detachment_folds_from_samples`, `implicit_faults_series_from_data`,
`hybrid_fault_fold_conditioning`.

## 2.13 Grids and properties (2026.3)

A top-level `grid` block grids any model into a corner-point grid that honours its faults and
horizons:

```yaml
grid:
  spacing: "[200, 200] m | user_instruction | areal cell size"
  margin: "50 m | user_instruction | pillars above the top and below the base horizon"
  zones:
    - {unit: Unit_1, layers: "4 ordinal | user_instruction | ...", layering: proportional}
    - {unit: Unit_2, layers: "3 ordinal | user_instruction | ...", layering: follow_top, layer_thickness: "70 m | user_instruction | ..."}
  facies: [{unit: Unit_1, code: "1 ordinal | user_instruction | ...", name: sand}]
  porosity: compaction                       # from the stage 2D compaction laws
  transforms:
    - {name: PERMX, unit: mD, from: PORO, law: exponential, a: "0.01 mD | assumed | ...", b: "30 ratio | assumed | ..."}
  fault_transmissibility: {calibration: manzocchi_1999}   # needs fault_rock and PERMX, PERMY
  tetrahedral: true
```

- **Pillars.** Pillars lie on the faults; cells take the zone and fault block they are in. Cells
  thinner than `min_thickness` are inactive, with the reason.
- **Properties.** Each carries a unit, an evidence class and a source. A value that cannot be
  computed (a unit without a compaction law, a face without a shale gouge ratio) is left missing,
  counted and explained in the summary table; it is never defaulted.
- **Transmissibility.** Multipliers (MULTX, MULTY) follow Manzocchi et al. (1999) from the
  shale gouge ratio and fault-zone thickness of 2.11. `grid/fault_faces.csv` lists every fault
  face with its inputs and multiplier.
- **Outputs.** `grid/cells.csv`, grid sections, layer maps and 3D views. Exports `grdecl`,
  `grid_vtk` and `resqml` (5).
- **Checks.** C23 (positive cell volumes), C24 (no inverted cells), C25 (fault faces on the
  fault surfaces).
- **Rule F32** checks the block.

Examples: `fault_system_relay_splay`, `rift_inversion_history`.

## 2.14 Uncertainty and ensembles (2026.3)

A top-level `uncertainty` block lets labelled parameters vary. A parameter without a
distribution stays at its value.

```yaml
uncertainty:
  members: "200 ordinal | user_instruction | ensemble size"
  seed: "20261115 ordinal | user_instruction | recorded in the manifest"
  parameters:
    - {path: faults.F1.throw.max, distribution: normal, mean: "120 m | assumed | ...", sd: "15 m | assumed | ..."}
    - {path: wells.W1.tops.H2, distribution: uniform, min: "1180 m | assumed | ...", max: "1220 m | assumed | ..."}
  correlations: [{a: faults.F1.throw.max, b: faults.F3.throw.max, rho: "0.6 ratio | assumed | ..."}]
  scenarios:
    - name: splay
      options:
        - {name: with_F2, weight: "0.7 fraction | assumed | ..."}
        - {name: without_F2, weight: "0.3 fraction | assumed | ...", remove: [faults.F2]}
  one_at_a_time: true
  probability_maps: [{kind: unit, unit: Unit_1, depth: "1150 m | ...", spacing: "200 m | ..."}]
```

- **Paths and distributions.** Paths name faults, primitives, horizons and wells by name.
  Distributions are `uniform` (min, max), `normal` (mean, sd, optional min and max), `triangular`
  (min, mode, max) and `discrete` (values, weights). Scenario options `set` description nodes or
  `remove` list items.
- **Running.** `stratonix ensemble model.yaml -o out/ens` builds every member with its checks,
  in parallel. The same seed gives the same files on every computer and for any `--threads`.
- **Outputs.**
  - `ensemble/members.csv`: every member, with seed, values, options, status, reason, failed
    checks and key numbers. A failed member is listed, never dropped.
  - `statistics.csv`: P10, P50 and P90 of each key number, with N and the seed.
  - `sensitivity_sobol.csv`: first-order indices.
  - `sensitivity_oat.csv`: one-at-a-time swings.
  - `maps/` probability maps, and `figures/` histograms, index charts, tornado charts and maps.
  - The base case's full build in `base/`.
- **Rule F33** checks the block.

Example: `fault_system_uncertainty`.

## 2.15 Geological time and topology (2026.4)

- **Geological time.** Give every horizon an age (`stratigraphy.horizons[].age`, Ma, young to
  old). Each point then has depositional coordinates (u, v, t): (u, v) from the model's own
  restoration and t interpolated between the horizon ages in the restored column. The build
  writes a Wheeler diagram for every section (`figures/*_wheeler_<section>`): distance against
  time, coloured by unit, with hiatuses blank and fault cut-outs grey.
- **Topology.** Every build writes `topology.json`: the units in each fault block and the
  surfaces between them, with sampled areas (lattice `outputs.voxel_spacing`, or the largest
  extent / 64). Check C26 refuses contacts between non-adjacent units that no unconformity or
  termination explains; check C27 measures how far horizon edges stop short of the fault,
  horizon or boundary that ends them. In an ensemble, `members.csv` states for every member
  whether its contacts match the base case's, and lists those that differ.

## 2.16 Gravity and magnetics (2026.4)

A top-level `geophysics` block gives the units' densities (kg/m3) and susceptibilities (SI), the
background values every other unit takes, and the inducing field for magnetics:

```yaml
geophysics:
  cell: "500 m | ..."                       # prism size; default the voxel spacing
  background_density: "2650 kg/m3 | ..."
  inducing_field: {inclination: "60 deg | ...", declination: "10 deg | ...", intensity: "50000 nT | ..."}
  units:
    - {unit: U1, density: "2350 kg/m3 | ...", susceptibility: "0.01 SI | ..."}
  gravity: {input: GRAV, format: csv_stations}
  magnetics: {input: MAG, format: csv_stations}
```

- **Forward model.** The model is divided into columns of `cell` x `cell`; in each column the unit
  boundaries are found at the column centre, and each run of one unit is a right rectangular prism
  (gravity after Nagy et al. 2000; magnetics from the prism's potential tensor, induced
  magnetisation only). Anomalies are contrasts against the background.
- **Stations.** `csv_stations`: columns X, Y, Z (elevation), value (mGal or nT) and optional
  sigma, read as delivered: state their corrections in the input's source; Stratonix applies none.
- **Outputs.** `geophysics_stations.csv` (observed, modelled, residual), anomaly maps and profiles
  along the sections, and check C29 (reduced chi-square, warning above 2).
- **Inversion.** As conditioning data (`kind: gravity` or `magnetics`) the stations fit the
  structural parameters listed under `conditioning.parameters`, alone or with picks and well tops,
  by least squares or by the Bayesian sampler. Densities and susceptibilities stay as stated.
- **Rule F37.** Listed units must exist; gravity needs the background density and one unit
  density; magnetics needs one susceptibility and the inducing field.

Example: `bayesian_gravity_anticline` (picks and gravity, Bayesian, with Wheeler diagrams).

## 2.17 Bodies: salt sheets, sills, dykes and injectites (2026.4)

A top-level `bodies` list places closed bodies that cut across or fold over the stratigraphy.
Each has a name, a `type` and an emplacement `age` on the faults' ordinal scale (larger =
younger):

| Type | Keys | Effect on the host |
| --- | --- | --- |
| `salt_sheet` | `footprint` polygon, `base` (`horizon` + `offset`, or `depth`), `thickness`, `taper` | strata above the sheet drape it (lifted by its local thickness) |
| `sill`, `saucer_sill` | `centre`, `depth` (mid-plane), `radius`, `thickness`; saucer: `outer_radius`, `climb` | the roof is lifted by `accommodation` x the local thickness, decaying upward by `uplift_decay` per km |
| `dyke` | `trace`, `width`, `top_depth`, `base_depth` | none |
| `injectite` | `parent` unit, `ceiling` horizon, `area`, `sheets`, `length`, `thickness`, `climb`, `net_to_gross`, `max_fraction`, `seed` | the injected volume is taken from the parent's sand (mass balance in the log) |

- Bodies are units of their own in sections, maps, grids, exports and checks.
- Restoration evaluates the host without them.
- A body crossed by a younger fault is refused with the reason: displacing bodies by younger
  faults is not built, so give the body an age after the fault's.
- Check C30 reports each body's volume and flags bodies thinner than half the grid spacing.
- Rule F38 covers names, types, ranges and parent/ceiling order.

Examples: `salt_sheet_drape`, `saucer_sill_forced_fold`, `injectite_complex`.

## 2.18 Fault systems, breached relays and growth on more hanging walls (2026.4)

A top-level `fault_systems` list writes a whole network once. Before the description is typed,
each system expands into ordinary faults named `<system>_01`, `<system>_02`, ... They are
appended to `faults`, carry the system's evidence class and a source naming the system, and are
listed in the resolved description. Every rule and check then applies to each of them as to a
fault written by hand (C10 scaling included). Rule F43 covers the systems.

| `type` | Keys | Expands into |
| --- | --- | --- |
| `polygonal` | `unit_top`, `unit_base` (horizons), `area` [[x0, y0], [x1, y1]], `reference_depth`, `spacing` (100 to 2000 m), `dip_mean` (30 to 70 degrees), `dip_sd`, `throw_ratio` (throw / length, at most 0.1), `seed` | planar normal faults on a seeded jittered lattice, tips on the two horizons, conjugate dip directions, throws scaled by length; at most 48 per system |
| `flower` | `master` (a planar fault dipping at least 80 degrees), `sense` (`positive` or `negative`), `splays` (2 to 9), `spread` (10 to 60 degrees), `vertical_component`, `branch_depth` | the non-vertical splays, reverse (positive) or normal (negative), each abutting the master on its own side and sharing the vertical component; the master keeps the horizontal slip |
| `duplex` | `start` (the hinterland ramp, midway between floor and roof), `strike` (ramps dip to its right), `floor_depth`, `roof_depth`, `horses` (at most 12), `horse_length`, `ramp_angle` (10 to 45 degrees), `slip` (per horse), `length` | ramp-flat thrusts (floor flat, ramp, roof flat) on fault-parallel-flow hanging walls, each younger than the one behind it; the style (hinterland-dipping, antiformal stack, foreland-dipping) follows from slip / horse length and is named in the source |

Every system takes a `name` and an `age`; generated faults take the age plus 0.001 per fault in
order.

**Breached relays.** `relay: {fault: G, breached: true}` replaces the soft-linked ramp with a
connecting fault `<F>_<G>_breach` across the middle of the overlap. It contains the relay faults'
slip vector, so neither carries its hanging wall across it. Both relay faults terminate against
it, each on its own side, so the throw passes from one to the other through the breach. The
connecting fault has no slip of its own (`linking: true`): its offset is the linked faults'
throw, so C4 and C10 skip it.

**Growth strata and phases** move each layer by its own throw. Since 2026.4 they work on planar
faults with `rigid`, `decay` and `linear_between_faults` hanging walls as well as `layer_rigid`.
The hanging wall's distance factor scales each layer's throw. Curved and fold-kinematic hanging
walls are not combined with growth yet. Faults that move level by level must be all the faults
of a model or none.

**Grids of dense networks.** `grid.pillars: vertical` keeps every pillar vertical, with the
faults stair-stepped across cells. Use it when many intersecting faults would tilt pillars into
each other. The default is `fault_conforming`.

Examples: `polygonal_faults`, `negative_flower`, `thrust_duplex`, `breached_relay`,
`growth_fault_decay`.

## 2.19 Clinoform trajectories and geobodies (2026.4)

**Clinoforms.** A series' `clinoforms` block also takes:
- `profile: tangential`, a foreset that leaves the break at the oblique profile's slope and
  flattens exponentially into the bottomset;
- `trajectory`, the angle of rise of the shelf-break trajectory (0 to 15 degrees, rule F23). Each
  younger clinoform's topset stands higher by `progradation` x tan(trajectory), so the topsets
  aggrade. A falling trajectory would cut the older topsets and is not built.

Example: `tangential_clinoforms`.

**Geobodies.** A top-level `geobodies` list places sets of channel belts or lobes inside a unit
between two horizons. Each set takes a `name`, a `type` (`channels` or `lobes`), the host `unit`,
a `count`, a `seed`, a `thickness`, an `orientation` and an optional `color`:

| `type` | Keys | Geometry |
| --- | --- | --- |
| `channels` | `width`, `wavelength`, `amplitude` (0 = straight), `shape` (`rectangular` or `parabolic`) | sine centrelines across the model, seeded offsets, levels and phases; sinuosity at most 3 |
| `lobes` | `radius` (semi-major axis), `aspect` (major / minor, at least 1) | ellipses in plan with a flat base and a lens-shaped top, seeded centres and levels |

- **Placement.** Sets are placed in geological-time coordinates: the restored position of a
  point and its fraction of the unit's restored thickness. They follow the unit's folding,
  faulting and thickness changes.
- **Overlap.** A later set takes a point before an earlier one; within a set, bodies overlap.
- **Model.** Each set becomes a unit of its own (listed under "Geobodies" in figure legends).
- **Checks.** C30 reports each set's net-to-gross in its unit, with the closed-form value beside
  it.
- **Rules.** Rule F39 covers names, types, ranges (width-to-thickness 5 to 500) and the host unit.
- **Model types.** Geobodies need a model with a restored stratigraphic column (parametric,
  implicit and conditioned models have one).

Example: `channels_and_lobes`.

## 2.20 Restoration, compaction across faults, backstripping and eroded section (2026.4)

**Volumetric restoration** (`restoration.volumetric`) restores the whole model to one horizon.
- **Keys:** `horizon`, its `datum` (depth), a `pin` [x, y], `young`, `poisson`, and optionally
  per-unit `units: [{unit, young, poisson}]`, `mesh_size` and `max_strain` (default 10 %).
- **Mesh:** the model below the horizon is meshed in tetrahedra on a lattice, cut by fault block.
  Faults slide (their normal gap is held).
- **Solve:** the horizon is brought flat to its datum by minimising elastic strain energy, in
  eight increments so rotations are not taken as strain.
- **Outputs:** `restoration_strain.csv` (the tetrahedra the horizon crosses: dilatation, maximum
  shear, bed stretch) and two strain maps.
- **C31:** fails when the largest Green-Lagrange strain exceeds `max_strain`, or the restored
  horizon's in-plane stretch exceeds `restoration.bed_length_limit` (default 5 %).
- **Limits:** the lattice does not conform to the faults, so strain peaks in the cells a fault
  crosses. A doubly curved horizon (a dome) cannot flatten without stretch.

**Bed length along sections.** In models with fault-parallel flow or a `restoration` block, C31
also reports, for every section, the largest change in bed length of any horizon (present
against restored).

**Compaction across faults.** `stratigraphy.compaction.burial: block` (the default since
2026.4) takes each unit's porosity at its maximum burial in its own fault block. That is the
deeper of its present and restored depths, so a hanging wall buried further by a fault compacts
more. `burial: restored` keeps the stage 2D column-wise behaviour.

**Backstripping** (`backstripping`) removes the sediment load from the burial history:
- **Keys:** `basement` horizon, `mantle_density`, `water_density`, `grain_density`, `te` (0 =
  Airy, up to 100 km; with `young` and `poisson` for flexure), and `steps` (paleo-water depth
  `water_depth` and `sea_level` per horizon at the surface, and optionally `present`).
- **Outputs:** tectonic subsidence (Airy, and flexural through the plate's response on a map grid
  when Te > 0) in `tectonic_subsidence.csv` at the wells, `tectonic_subsidence_maps.csv`, and
  dashed and dotted curves on the burial-history figures.

**Eroded section** (`eroded_section`) estimates the section removed at an unconformity:
- **Keys:** a `well`, the `unconformity` horizon, a porosity `log` input (depth below the datum,
  porosity) and the normal trend (`surface_porosity`, `decay`).
- **Method:** each sample's equivalent burial on the trend minus its depth below the
  unconformity, averaged (Magara 1976).
- **Output:** `eroded_section.csv` gives the residual; an optional `vitrinite_estimate` is
  reported beside the estimate. Underburial is reported as zero erosion with the reason.

Rule F40 covers all three. Example: `restoration_and_backstripping`.

## 2.21 Mechanics, fractures, heat and pressure (2026.4)

**Mechanical stratigraphy** (`mechanics`) gives units their own horizontal stresses. Each unit
takes `young`, `poisson` and optionally `ucs` and `biot`, with optional `tectonic_strain: {min,
max}` (along Shmin and SHmax, from -0.005 to 0.005). With a `stress.regime` (which gives the
vertical stress, the pore pressure and the SHmax azimuth), the horizontal stresses follow the
uniaxial-strain relation plus the tectonic strains (Thiercelin and Plumb 1994), unit by unit.
Slip and dilation tendencies on faults (C20) use the footwall's unit, and `stress_profile.csv`
lists the stresses along every well.

**Fracture intensity** (`fractures.units`). The law is P32 = `p32_0` + `a_k` |k_max|, clipped
to `p32_max`, where k_max is the maximum absolute principal curvature of the unit's top horizon.
It is drawn as a map per unit and tabled in `fracture_intensity.csv`; cells without curvature are
reported missing.

**Fracture networks** (`fractures.dfn`) place Baecher disks in a box: the `area` over the full
thickness of the `unit` at its centre.
- **Sets:** each takes a `pole` [trend, plunge] with a Fisher `kappa`, a truncated power-law
  radius (`exponent` > 2, `rmin`, `rmax`) and a `weight`.
- **Generation:** disks are drawn until the in-box area reaches the box's mean P32 from the
  intensity map. A network above `limit` disks (default 10^6) is refused with the reason.
- **Outputs:** `dfn_disks.csv`; C34 compares the realised P32 with the target.

**Heat** (`heat`) runs transient conduction along 1D columns through the burial history.
- **Inputs:** `surface_temperature`, `basal_heat_flow` (or a `history` of ages and flows),
  `water_conductivity`, and per unit the grain `conductivity` (mixed with water by the geometric
  mean of porosity) and `heat_production`. Every unit with a compaction law needs them.
- **Outputs:** temperature (and EASY%Ro with `vitrinite: true`) of every horizon at every step
  at the wells, present profiles, and maps of the present temperature and %Ro at each horizon.
- **C32** checks the heat balance of the present column in steady state.

**Pressure** (`pressure`) models overpressure from disequilibrium compaction through the burial
history, in Gibson's (1958) moving-boundary form.
- **Inputs:** compressibility from the unit's compaction law in effective stress, and
  permeability from its law (`exponential`: k0 exp(b (phi - phi0)); `kozeny_carman`), with the
  `viscosity` and densities.
- **Outputs:** overpressure per horizon and step at the wells, present profiles with their
  hydrostatic and lithostatic bounds, and overpressure maps.
- **C33** fails a pore pressure outside those bounds.

Rule F41 covers heat and pressure, F42 mechanics and fractures. Units: degC, W/m/K, mW/m2, uW/m3,
mPa.s (cP). Example: `basin_heat_pressure_fractures`.

# 3. Running

| Command | What it does |
|---|---|
| `stratonix validate model.yaml` | schema and semantic rules only |
| `stratonix build model.yaml -o out/run1` | build, run every check, write summary, checks report, summary table, figures and the exports listed in the description, and `manifest.json` |
| `stratonix section model.yaml -o out/s1 --line "x1,y1;x2,y2;..."` | one section along any polyline (or `--name` of a section in the description) |
| `stratonix export model.yaml -o out/e1 --formats petrel,vtk` | exports only |
| `stratonix package model.yaml -o out/delivery --config package.yaml --zip` | delivery package (section 6) |
| `stratonix state model.yaml --age 45 -o out/45Ma` | the model at an event of its timeline, with its package (2026.2) |
| `stratonix ensemble model.yaml -o out/ens [--members N] [--seed S]` | the seeded ensemble of the description's uncertainty block (2026.3) |
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
| C10 fault scaling | WARNING: a fault's displacement-length ratio is outside the stated range |
| C11 fault system | throw jumps across a branch line or a relay |
| C12 area | a fold mechanism changes a unit's area beyond the tolerance |
| C13 thickness restored | a series relation leaves a negative thickness |
| C14 truncation order | a truncation contradicts the age order |
| C15 salt volume | salt is created or lost by the bodies and their withdrawal |
| C16 solid thickness | compaction and decompaction change a unit's grain thickness |
| C17 restorability | a state of the timeline does not build, or the present state differs from the model |
| C18 balance | WARNING: line length or area differs between present and restored |
| C19 juxtaposition | the units on a fault's walls disagree with the restored stratigraphy |
| C20 tendencies | a slip or dilation tendency lies outside its bounds |
| C21 data misfit | fewer than 95 % of the data with an uncertainty lie within 3 sigma, or an exact datum is missed |
| C22 conditioning | a fitted parameter did not converge strictly inside its range (WARNING when the reduced chi2 exceeds 4) |
| C23 cell volumes | an active grid cell has a volume of zero or less |
| C24 inverted cells | an active grid cell is inverted (corner order or corner Jacobian) |
| C25 fault faces | a fault face of the grid lies off its fault surface by more than 1 % of the cell size |
| C26 topology | units in contact inside one fault block are not adjacent in the column, and no unconformity (change of series) or terminating unit explains it |
| C27 watertight | a horizon edge ends further than the tolerance from any fault, other horizon or the domain boundary (`checks.tolerances.watertight`, default half the surface spacing) |
| C30 bodies | a body has no volume inside the model (bodies thinner than half the grid spacing are reported) |
| C29 potential fields | WARNING: the stations' reduced chi-square exceeds 2 |
| C28 posterior | a parameter's split R-hat is 1.01 or more, or its effective sample size is below 100: run the sampler longer |

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

From 2026.2 the exports add, where the description asks for them:

- `porosity/` (compacting units), `isopachs/` and `lines/` (stratal terminations), and
  `curvature/` (nine grids per horizon), all IRAP classic;
- fault attributes as VTK point data on the fault surfaces;
- `csv/<fault>_fault_attributes.csv`, one row per surface node.

From 2026.3, with a `grid` block (2.13):

| Id | Files |
|---|---|
| `grdecl` | `grid/<model>_grid.grdecl`: SPECGRID, COORD, ZCORN, ACTNUM (depth positive down), ZONE, BLOCK, FACIES and the real properties with unit and evidence in comments, MULTX, MULTY and FAULTS; missing values written as -999.25 and explained in the header |
| `grid_vtk` | `grid/grid_cells.vtu` (hexahedra) and, with `tetrahedral: true`, `grid/grid_tetrahedra.vtu` |
| `resqml` | RESQML 2.0.1 EPC and HDF5 (IJK grid, split coordinate lines, properties). It needs a build with FESAPI (`STRATONIX_WITH_FESAPI`); other builds say so and write nothing |

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

# 7. Known limitations

- Hanging-wall fields: phases and growth need `layer_rigid` hanging walls; rigid, shear and
  fold-kinematic faults take one event, after all deposition.
- Restoration is geometric (kinematic), not mechanical. Area balance of growth strata is reported
  for the syn- and post-kinematic units, not tested.
- Compaction is not carried across faults that moved during burial.
- Fault rock needs a restored stratigraphic column: implicit models report SGR and CSP as
  missing.
- Stress and pore pressure are proportional to depth and laterally uniform; there is no
  mechanical stratigraphy.
- Dip-slip moves the hanging wall down the dip; vertical offset moves it straight down. They differ
  by up to a few metres on dipping horizons.
- A fault's throw must change by less than 1 m per metre of slip; elliptical dip profiles are
  refused for that reason.
- Where a reverse fault repeats a horizon, grids keep the upper occurrence (reported as overlap
  nodes).
- At most 64 faults per model.
- Grids follow faults oblique to the lattice in stair steps, and are of the present-day model
  only. RESQML export needs a build with FESAPI; validation against the RESQML schemas needs the
  Energistics schema files.
- A unit that thins out within a series must be declared in `implicit.isochore_horizons`.
- Sensitivity gives first-order indices only; the key numbers are those of the summary table.
