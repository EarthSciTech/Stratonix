# Contributing to Stratonix

## Ground rules

- **No fabrication.** Every parameter in a description carries value, unit, evidence class and
  source. Code must never invent a default that changes a model; built-in defaults are limited to
  check tolerances and are reported as such.
- **Generality.** No site, field or project name in code, defaults, examples, tests or docs.
  Real models live outside the source tree.
- **Layers.** Lower layers never depend on higher ones (`cmake/LayerCheck.cmake` enforces it).
  Representations, checks, exports and figures use only the `Model` interface.
- **Portability.** Standard C++20. No compiler extensions without a fallback. Text outputs are
  written in binary mode with `\n`, numbers with `fmt_fixed` / `fmt_shortest` (locale
  independent), so files are byte-identical on Windows and Linux.
- **Determinism.** No timestamps outside `manifest.json`. Parallel code uses `parallel_for`, whose
  chunks do not depend on the thread count, and writes each output element once.

## Adding a plugin

1. Derive from the plugin interface of its kind (`StructureBuilder`, `Exporter`, `Importer`,
   `Check`, ...), in the layer that owns that kind.
2. Register it in its `.cpp` file:
   ```cpp
   STRATONIX_REGISTER(my_exporter, registry::Kind::Exporter, "my_format", "1.0.0",
                      "One line on what it writes", [] { return std::make_unique<MyExporter>(); })
   ```
3. Add `STRATONIX_PLUGIN_ANCHOR(my_exporter);` to the layer's `register_*_plugins()` function so
   the linker keeps it.
4. Add tests: a unit test, and for exporters a round trip (write, read back, compare) plus an
   independent reader in `tools/` when one exists.
5. `stratonix plugins` must list it.

Run-time plugins (shared libraries) will use the C ABI in `include/stratonix/registry/plugin_abi.h`.

## Adding an optimised kernel

1. Keep (or first write) the **reference implementation**: plain, readable, single-threaded.
   It remains the specification of correct behaviour.
2. Write the optimised version behind the same function signature; select it with
   `optimised_kernels()` so tests and users can switch it off (`--reference-kernels`).
3. Add a **parity test**: synthetic models and randomised inputs, including points on a fault
   plane, at a tip and outside the domain; state the tolerance.
4. Add a **benchmark** reporting time and memory for both versions at two sizes.
5. Results must be identical for any thread count; document thread safety in the header.
6. The optimised version is enabled only when parity and benchmark both pass. Optimise what the
   profile shows, not what you expect.

## Tests

`ctest` runs unit tests, the model-interface suite, the layer check, the C API test, command-line
tests and (from later stages) regression, parity and round-trip tests. Sanitiser builds:
`cmake --preset linux-asan`.

## Release numbers

Releases are numbered `<year>.<major>`: the year the release is made, then the major version
(`2026.1`, then `2026.2` or `2027.2`). A fix to a published release adds a patch number
(`2026.1.1`). The number is set in one place, `project(... VERSION ...)` in `CMakeLists.txt`,
and the git tag is `v` followed by it.
