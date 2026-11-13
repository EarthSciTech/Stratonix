# Third-party components

| Component | Use | Licence | Linking |
|---|---|---|---|
| libyaml | YAML and JSON reading | MIT | static where available |
| zlib | PNG, zip and PPTX compression | zlib | static where available |
| Cairo (optional) | rasterising figures to PNG and TIFF | LGPL-2.1 or MPL-1.1 | **dynamic only** (DLLs shipped beside the program) |
| PROJ (optional, off by default) | coordinate reprojection | MIT | static |
| C++ standard library | everything | per toolchain (GCC runtime exception, MSVC redistributable terms) | static runtime |

Not used: CGAL (GPL parts), any GPL library. Planned when a networked build is available:
Eigen (MPL-2.0), VTK (BSD-3), FESAPI (Apache-2.0) with HDF5 (BSD-style), CLI11 (BSD-3),
Catch2 (BSL-1.0), Google Benchmark (Apache-2.0). See STX-TN-01 section 11.
