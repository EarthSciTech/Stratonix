---
title: Stratonix roadmap v2 to v5
doc_id: STX-RM
version: 1.0
date: 2026-11-13
---

# Stratonix roadmap v2 to v5

Updated at release 1.0 from what building v1 taught. Effort is not estimated here; each stage
starts with its own specification.

# 1. What v1 taught

| Lesson | Consequence for later stages |
|---|---|
| Restoration through several faults needs a well-posed inverse: the throw had to be read at the restored position, tips of younger faults had to be resolved on displaced stratigraphy, and a throw gradient of 1 m per metre of slip or more makes the inverse non-unique (now refused). A floating-point edge case at the end of the solver's bracket made a check fail at five positions. | Any new kinematics must state its inverse and its admissible range, and come with an analytical test case. Keep root finding bracketed. |
| The checks found real defects every time the model grew (C4 throw, C6 unoffset base, the elliptical profile). | Each new engine plugs into the same checks; add checks before features. |
| Optimised kernels gave 22 to 31 times on the full build, almost all in faulted grids and checks; parity tests caught nothing wrong but made the changes safe. | Keep the reference-and-parity rule; the implicit engine's solver is the next heavy kernel. |
| Golden hashes made byte identity testable; the C++ number formatting (`to_chars`) is the only formatting used in files. GCC and Clang wrote identical bytes. | Keep golden hashes for every new format; run them on Windows in CI. |
| The external cross-check showed that an interpolator without faults is unusable within about 500 m of the faults and that block-wise interpolation is good to 0.1 m away from faults and 0.5 m near them on the synthetic model. | v2 must be fault-aware from the start (per-block fields or fault-displacement functions), and its tests should use the same cross-check. |
| No network, no Windows host: VTK, PROJ, Eigen, FESAPI, xtgeo, LoopStructural and Petrel were unavailable, so their parts were replaced by own code or deferred. | Before v2: one networked Windows and Linux build with vcpkg pinned; RESQML through FESAPI; Petrel import validated. |

# 2. v2: implicit modelling

**Adds.** A native implicit interpolator: horizons from points and orientations (well tops,
interpreted points, dips), faults as interfaces from sticks or surfaces, with the fault network
reused from v1.

**Builds on.** The `Model` interface (a scalar per horizon, unit, block); the fault network object;
importers for points, grids, sticks, TSurf and tops; evidence classes on inputs; the model suite;
checks C1-C9; the cross-check script as an acceptance test.

**Interface changes.** None required for the model interface. Additions: an `ImplicitBuilder`
structure-builder plugin; a solver kernel (reference dense, optimised sparse or matrix-free) with
parity tests; `Param` for interpolation settings (range, nugget, regularisation) with evidence
`assumed` or `calibrated`; per-point uncertainty on inputs (an optional field on imported data).
Prompt 14 starts this stage.

# 3. v3: grids and properties

**Adds.** Corner-point and unstructured grids conforming to faults and horizons; properties;
upscaling hooks.

**Builds on.** The extractor plugin kind; the voxel extractor as the first grid extractor;
structure-of-arrays storage; fault surfaces and cut-offs from `repr`.

**Interface changes.** A grid representation type (`GridRepr`) and an `Extractor` plugin interface
returning it; a property container keyed by grid cell with units and evidence; exporters for
GRDECL and RESQML grids. The model interface stays as it is.

# 4. v4: uncertainty

**Adds.** Parameter sampling, ensembles, sensitivity of key numbers (closure, spill, throw at wells).

**Builds on.** Every parameter is a typed `Param` with value, unit, evidence and source; builds are
deterministic; a build is callable from the library and the C API without the command line;
key numbers are already a table.

**Interface changes.** A distribution on `Param` (optional; absent = fixed); an ensemble runner in
the package layer; seeds recorded in the manifest. The C API gains functions to set a parameter on
a loaded description before building (additive, minor C API version).

# 5. v5: graphical application, Petrel plug-in, bindings

**Adds.** A Qt application; a Petrel plug-in (Ocean); Python bindings.

**Builds on.** The thin command-line layer; the stable C API with opaque handles and status codes;
figures drawn on a display list (a GUI view is another backend).

**Interface changes.** Python bindings over the C API (no C++ ABI exposure); a canvas backend for
Qt; the Petrel plug-in calls the C API from .NET through P/Invoke. Event callbacks (progress,
cancel) are added to the C API as additive functions.
