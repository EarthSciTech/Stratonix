# Stratonix

Stratonix is a general 2D/3D geomodelling toolkit written in C++. It builds structural models
(horizons and fault planes) from a parametric design, checks them, and exports them in formats
other software can load, Petrel first.

Release 2026.2 (version 2) adds:

- fault systems: listric, ramp-flat and curved faults, splays, relays, drag, and oblique and
  strike-slip movement;
- fault-related folding: trishear, fault-bend and detachment folds;
- stratigraphic architecture: unconformities, onlap, clinoforms, isochores and growth strata;
- salt bodies and compaction;
- a deformation history with sequential restoration;
- faults as rock: Allan diagrams, shale gouge ratio, clay smear, seal capacity, slip and
  dilation tendency, and horizon curvature with fracture trends.

Next come data conditioning, grids and uncertainty (v3), applications (v4) and a graphical
application with bindings (v5).

- Library: `libstratonix` (namespace `stratonix`), one static library per layer.
- Program: `stratonix` (`stratonix.exe` on Windows).
- C API: `stratonix_c.h` (shared library `stratonix_c`).

## Quick start

```
stratonix validate examples/synthetic_fold_two_faults.yaml
stratonix build    examples/synthetic_fold_two_faults.yaml -o out/synthetic_v1
stratonix section  examples/synthetic_fold_two_faults.yaml -o out/line_1 --line "500,500;4000,3000;7500,5500"
stratonix bench    --spacing 100,50 --mode both -o out/bench
stratonix plugins
stratonix version
```

`build` writes the model summary, the checks report, a summary table, a manifest with SHA-256
checksums and, when every check passes, the figure set (sections, maps per horizon, two 3D views,
a fence diagram) in SVG, PNG and TIFF at 300 dpi and PPTX. Nothing is drawn or exported when a
check fails.

A model description is YAML (or JSON). Every numeric parameter carries a value, a unit, an
evidence class and a source, in full or compact form:

```yaml
throw_max: {value: 60, unit: m, evidence: assumed, source: "design case 2"}
throw_max: "60 m | assumed | design case 2"
```

Evidence classes: measured, interpreted, calibrated, upscaled, literature_derived, assumed,
user_instruction. A description with an unlabelled parameter is refused. Design models are
labelled DESIGN in every exported file.

## Documentation

- `docs/user_guide.md`: installing, writing a description (worked examples, and every version 2
  feature in sections 2.6 to 2.11), running, reading the checks report, exporting and packaging,
  known limitations. `examples/` holds 21 worked descriptions.
- `docs/developer_guide.md`: architecture, model interface and suite, build, adding plugins.
- `docs/api_reference.md` (generated from the headers; `doxygen docs/Doxyfile` for HTML).
- `docs/roadmap.md`: v3 to v5, and what v1 and v2 taught. `docs/performance_and_threading.md`, `docs/package_config.md`.
- C API: `capi/include/stratonix_c.h`, example `capi/examples/stratonix_example.c`.

## Building

Linux (system packages):

```
sudo apt-get install cmake ninja-build g++ libyaml-dev libcairo2-dev zlib1g-dev
cmake --preset linux-gcc-release
cmake --build build/linux-gcc-release
ctest --test-dir build/linux-gcc-release
```

Windows (Visual Studio 2022 and vcpkg):

```
set VCPKG_ROOT=C:\vcpkg
cmake --preset windows-msvc-release
cmake --build build/windows-msvc-release --config Release
ctest --test-dir build/windows-msvc-release -C Release
cd build/windows-msvc-release && cpack -C Release
```

Options: `STRATONIX_WITH_CAIRO` (PNG/TIFF figures, default ON), `STRATONIX_WITH_PROJ`
(reprojection, default OFF), `STRATONIX_SANITIZE` (address/UB sanitisers).

## Layout

| Folder | Contents |
|---|---|
| `include/stratonix/` | public headers, one folder per layer |
| `src/` | one folder and one CMake target per layer |
| `apps/stratonix/` | command-line program |
| `capi/` | C API, header and C test |
| `tests/` | unit, interface-suite, regression, parity and round-trip tests |
| `benchmarks/` | timing benchmarks |
| `examples/` | synthetic descriptions only (`broken/` holds deliberately invalid ones) |
| `schema/` | JSON Schema of model descriptions |
| `tools/` | external test tools (Python), not part of the product |
| `docs/` | Doxygen configuration and developer notes |

See `CONTRIBUTING.md` for how to add a plugin or an optimised kernel, and the design
specification STX-TN-01 for the architecture.

## Container image

```
docker run --rm -v "$PWD:/work" ghcr.io/earthscitech/stratonix:latest validate model.yaml
docker run --rm -v "$PWD:/work" ghcr.io/earthscitech/stratonix:latest build model.yaml -o out
```

Windows installer, portable zip and Linux archive: see the Releases page.

## Licence

MIT, see `LICENSE`. Third-party components and their licences are listed in `THIRD_PARTY.md`.

## Developer

Amirsaman Rezaeyan

## Citation

How to cite:

Amirsaman Rezaeyan, 2026, Stratonix v.2026.1. https://github.com/EarthSciTech/Stratonix
