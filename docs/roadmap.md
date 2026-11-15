---
title: Stratonix roadmap v3 to v5
doc_id: STX-RM
version: 2.0
date: 2026-11-15
---

# Stratonix roadmap v3 to v5

Updated at release 2026.2 from what building v1 and v2 taught. Effort is not estimated here; each
stage starts with its own specification and its tests, written before the code.

# 1. Releases so far

| Release | Content |
|---|---|
| 2026.1 | v1: parametric design models with planar faults, checks C1-C9, the Petrel set and neutral exports, figures, delivery packages; the implicit builder opened (1.1) |
| 2026.2 | v2: fault systems and curved faults (2A), fault-related folding (2B), stratigraphic architecture (2C), salt and burial (2D), deformation history and restoration (2E), faults as rock and the stress field (2F); checks C10-C20 |

# 2. What v1 taught (kept from version 1.1 of this roadmap)

| Lesson | Consequence |
|---|---|
| Restoration through several faults needs a well-posed inverse; a throw gradient of 1 m per metre of slip or more makes it non-unique (refused). | Every new mechanism states its inverse and admissible range, with an analytical test. Root finding stays bracketed. |
| The checks found real defects every time the model grew. | Add checks before features. |
| Golden hashes made byte identity testable. | Keep golden hashes for every format; new outputs are gated so older descriptions keep their bytes. |
| An interpolator without faults is unusable near faults. | Data conditioning must be fault-aware (stage 3A). |

# 3. What v2 taught

| Lesson | Consequence for later stages |
|---|---|
| Restoration-based evaluation scaled to every v2 mechanism: curved faults, shear and fault-parallel-flow hanging walls, trishear, series relations, salt, compaction and events. Each needed only an exact inverse and a restored-frame rule. | Grids (3B) and data conditioning (3A) read the model through the same two steps. A mechanism without an exact inverse is not admitted. |
| Engine-specific queries from upper layers broke the layering twice (stratigraphy in 2C, the restored column in 2F). Model-layer capability interfaces solved both. | New cross-layer needs get a capability interface in the model layer, offered or not by each engine; callers report "missing" when it is absent. |
| Timing rules unified growth strata, phases, inversion and backstripping: the part of an event after a horizon's deposition. | Ensembles (3C) can vary ages and phases like any parameter. |
| Area balance in the restored frame flags growth strata, which are thicker in the hanging wall by deposition. | Balance tests apply to pre-kinematic units; syn-kinematic units need a partial restoration (open item). |
| Several checks were first too strict at designed terminations (C1 near pinch-outs, C4 on trishear faults). | A check states the geometry it applies to and skips the rest explicitly, in the report. |
| Missing inputs must never be defaulted (fault rock, stress). | Applications (4A) report every quantity with the inputs it depends on and "missing" when one is absent. |
| Fault attributes are recomputed for each output (checks, figures, summary, exports). Builds of fault-rich models take tens of seconds. | Cache per-fault attribute surfaces in the export context before stage 3C multiplies builds. |

# 4. v3: data conditioning, grids and uncertainty

- **Stage 3A, data-conditioned modelling (Prompts 34 to 36).** A scalable implicit solver
  (partition of unity or sparse) with series and faults taken from data. Hybrid conditioning of
  parametric models to wells and points, with data uncertainty. It builds on the implicit builder
  (1.1), the fault network, the `column` capability and the restoration-based evaluation.
- **Stage 3B, grids and properties (Prompts 37 to 39).** Corner-point grids conforming to faults
  and horizons, with their export (GRDECL, RESQML). Properties, fault transmissibility from SGR
  and fault-zone thickness (stage 2F), and unstructured grids. It builds on the voxel extractor,
  fault surfaces and cut-offs, and the porosity from compaction (stage 2D).
- **Stage 3C, uncertainty (Prompts 40 and 41).** Parameter distributions on `Param`, ensembles
  with every member's checks recorded, sensitivity and probability maps. Prompt R then releases
  v3.

# 5. v4: applications (stage 4A, Prompts 42 to 47)

Shared tools (conditions at depth, traps, migration, volumes), then one application family per
prompt:

- petroleum, conventional and unconventional;
- CO2 and hydrogen storage;
- geothermal and groundwater.

The stage ends with synthetic seismic and the model report. It builds on traps and closure (v1),
fault seal and stress (2F), burial and porosity (2D) and grids (3B). Prompt R then releases v4.

# 6. v5: platform and reach (stage 5A, Prompts 48 to 50)

- **Browser viewer.**
- **Desktop application.**
- **Python bindings and the Petrel plug-in.**

These build on the thin command line, the stable C API with opaque handles and the display-list
figures. Prompt R then releases v5.
