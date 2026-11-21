---
title: Stratonix roadmap v4 to v6
doc_id: STX-RM
version: 3.1
date: 2026-11-21
---

# Stratonix roadmap v4 to v6

Updated at release 2026.3 from what building v1, v2 and v3 taught, and re-planned on 2026-11-21: a new v4
(geological engine and real data) comes before the applications, which move to v5, and the platform, which
moves to v6. Effort is not estimated here; each
stage starts with its own specification and its tests, written before the code.

# 1. Releases so far

| Release | Content |
|---|---|
| 2026.1 | v1: parametric design models with planar faults, checks C1-C9, the Petrel set and neutral exports, figures, delivery packages; the implicit builder opened (1.1) |
| 2026.2 | v2: fault systems and curved faults (2A), fault-related folding (2B), stratigraphic architecture (2C), salt and burial (2D), deformation history and restoration (2E), faults as rock and the stress field (2F); checks C10-C20 |
| 2026.3 | v3: data-conditioned modelling (3A), grids and properties (3B), uncertainty (3C); checks C21-C25; C API 1.2 |

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
| Missing inputs must never be defaulted (fault rock, stress). | Applications (5A) report every quantity with the inputs it depends on and "missing" when one is absent. |
| Fault attributes are recomputed for each output (checks, figures, summary, exports). Builds of fault-rich models take tens of seconds. | Cache per-fault attribute surfaces in the export context before stage 3C multiplies builds. |

# 4. What v3 taught

| Lesson | Consequence for later stages |
|---|---|
| One field per series was not enough: a unit that thins out within a series distorts a single field by up to 200 m, in Stratonix and in LoopStructural alike. Declared isochore horizons fixed it. | A potential-field formulation with free horizon values (open item) would remove the declaration. External cross-checks state how each code is set up for such units. |
| Fault rock read at a grid face had to be read on the fault, not at the face centre; corner-point stair steps put faces up to a cell away. | Applications (5A) that read fault properties at grid faces use the projection onto the fault. |
| Description documents edited by path (values, scenario options) and validated again made ensembles, C API parameter setting and hybrid conditioning one mechanism. | v4 parameters, v5 applications and v6 bindings set values through the same paths. |
| Byte-identical ensembles needed a generator and inverse CDFs written out in integer and IEEE arithmetic, and results written in member order, never in completion order. | Every parallel output keeps a fixed order; tests compare files from 1 and 4 threads. |
| Optional heavy dependencies (FESAPI, LoopStructural) are best built and exercised in their own CI jobs, with their interface printed in the log. | v4 to v6 dependencies (well-log and SEG-Y readers, PROJ, MODFLOW and simulator decks, viewers) follow the same pattern. |
| Ensembles multiply build cost: about 1 s per member of a three-fault model. | Cache per-fault attribute surfaces and checks' lattices in stage 4D, before v5 ensembles of application models. |

# 5. v4: geological engine and real data (stages 4A to 4D, Prompts 42 to 61)

- **4A, core engine:** one potential field per series with free horizon values, at field scale (closes the
  isochore declaration, open item I-037); folds and faults learnt from data (fold frames, fault displacement
  fields, fault networks); Bayesian conditioning with model gradients and a sampler (closes I-039);
  geological-time coordinates and Wheeler diagrams; topology and watertightness checks; gravity and
  magnetics with joint inversion.
- **4B, geological completeness:** salt overhangs and sheets (I-028); sills, dykes and sand injectites;
  polygonal faults, flower structures and thrust duplexes; breached relays (I-023), trishear on curved
  faults (I-025) and growth on every hanging wall (I-027, I-030); onlap, downlap and toplap, clinoform sets,
  channel and lobe geobodies; flexural-slip and 3D volumetric restoration with strain (I-026, I-020);
  compaction across moving faults (I-029); flexural backstripping and eroded section; fractures, mechanical
  stratigraphy and lateral stress (I-033), fault rock for implicit models (I-032); temperature and
  overpressure through time.
- **4C, data real projects start from:** well logs, surveys, tops and ages; seismic picks and SEG-Y;
  time-to-depth conversion calibrated to wells; real coordinates through PROJ; geological maps, outcrop
  data and elevation models; pick and velocity uncertainty carried into depth; RESQML schema validation
  (I-034) and RESQML in the release builds.
- **4D, properties, uncertainty and scale:** variograms, sequential Gaussian and indicator simulation and
  object-based facies in geological time; grids that follow fault traces (I-035); total-effect sensitivity
  (I-036); caching, grids larger than memory and parallel builds; Python bindings on PyPI.

Prompt R then releases v4.

# 6. v5: applications (stage 5A, Prompts 62 to 68)

Shared tools (conditions at depth from the 4B histories, traps, migration, volumes), then one application
family per prompt:

- petroleum, conventional and unconventional;
- CO2 and hydrogen storage;
- geothermal and groundwater.

Then a simulation link (decks for OPM Flow, run as an external program, and history matching by ensemble
smoother), synthetic seismic, the model report and a geological validation report for any model. Prompt R
then releases v5.

# 7. v6: platform and reach (stage 6A, Prompts 69 to 72)

- **Browser viewer.**
- **Desktop application that edits models:** every edit is written back as a labelled change to the
  description.
- **Full Python API and the Petrel plug-in,** including validating Petrel models from inside Petrel.
- **Model provenance:** comparing model versions, a decision record, and checks run on every model change.

These build on the thin command line, the stable C API with opaque handles and the display-list figures.
Prompt R then releases v6.

# 8. After v6

A benchmark against GemPy and LoopStructural on their public example models; a real-data case study if the
data owners agree; and a software paper.
