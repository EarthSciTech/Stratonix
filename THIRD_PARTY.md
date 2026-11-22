# Third-party components

| Component | Use | Licence | Linking |
|---|---|---|---|
| libyaml | YAML and JSON reading | MIT | static where available |
| zlib | PNG, zip and PPTX compression | zlib | static where available |
| Cairo (optional) | rasterising figures to PNG and TIFF | LGPL-2.1 or MPL-1.1 | **dynamic only** (DLLs shipped beside the program) |
| PROJ (optional; in the release packages from 2026.4) | coordinate transformations | MIT | dynamic (Windows: DLL shipped; Linux: system library) |
| SQLite (optional; in the release packages from 2026.4) | GeoPackage maps | public domain | dynamic (Windows: DLL shipped; Linux: system library) |
| FESAPI (optional; in the release packages from 2026.4) | RESQML 2.0.1 grid export | Apache-2.0 (licence text shipped as `FESAPI_LICENSE.txt`) | dynamic (library shipped beside the program) |
| HDF5, minizip (with FESAPI) | RESQML arrays and packages | BSD-style, zlib | dynamic (Windows: DLLs shipped; Linux: system libraries) |
| C++ standard library | everything | per toolchain (GCC runtime exception, MSVC redistributable terms) | static runtime (Windows builds with FESAPI: the MSVC DLL runtime) |

Not used: CGAL (GPL parts), any GPL library. Test tools only (never shipped): lasio, dlisio,
dliswriter, segyio, resqpy, LoopStructural, harmonica. The RESQML 2.0.1 XML schemas (Energistics
Product Licence) are never committed or shipped; CI uses them when the owner supplies them.
Planned when needed: Eigen (MPL-2.0), VTK (BSD-3), CLI11 (BSD-3), Catch2 (BSL-1.0), Google
Benchmark (Apache-2.0). See STX-TN-01 section 11.
