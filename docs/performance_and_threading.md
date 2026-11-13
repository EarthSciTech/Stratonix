# Performance, optimised kernels and thread safety

## Reference and optimised kernels

Every heavy kernel has a plain reference implementation, which is the specification of correct
behaviour, and an optimised implementation behind the same function. The process-wide switch
`set_optimised_kernels(bool)` selects between them; the command line exposes it as
`--reference-kernels` (`build`) and `--mode reference|optimised|both` (`bench`).

| Kernel | Reference | Optimised | Parity test (tests/unit/test_parity.cpp) |
|---|---|---|---|
| K1 point evaluation (structure + fault displacement) | per point, each horizon evaluated independently | shared structure terms across horizons; parallel over points (chunks of 1024, n >= 2048) | bitwise equal; equal for 1, 3 and 8 threads |
| K2 per-block horizon grids | bisection per node | regula falsi (Illinois) seeded from a coarse grid (stride 4) | within 1e-6 m |
| K3 faulted horizon grids | dense sampling of 96 levels per node, then bisection | block iteration with Illinois seeded from coarse grid | within 1e-6 m where both are defined; gap masks identical |
| K4 cut-off lines | bisection on the extrapolated scalar | Illinois on the same bracket | within 1e-4 m |
| K5 voxel grids | per voxel unit lookup | column-wise lookup from the horizon roots | identical |
| Tip resolution (build) | serial table of 513 samples | same samples, parallel, Illinois | covered by K1 (identical model) |

Edge cases in the parity tests: points exactly on a fault plane, at a tip line, on a horizon, and
outside the domain (all must agree, including NaN for undefined).

## Determinism

`parallel_for(n, chunk, fn)` splits `[0, n)` into fixed chunks whose boundaries do not depend on
the number of threads; each index is written by exactly one call and no reduction depends on
scheduling order. Results are therefore identical for any thread count (`--threads N`, `set_thread_count`). Tests run the K1 kernel with 1, 3 and 8 threads and compare bitwise.

## Thread safety

- **Models** are immutable after `build_model` returns. `Model::evaluate`, `scalar_in_block`,
  `evaluate_fault` and `restore` are `const` and re-entrant: they may be called from any number
  of threads on the same model.
- **Builders** (`StructureBuilder::build`, the fault network's tip resolution) are not re-entrant
  on the same output object; build one model per thread.
- **Registry**: plugins are registered once by `register_all_plugins()` (call it before starting
  threads). After that, lookups are read-only and safe from any thread.
- **Global switches** (`set_thread_count`, `set_optimised_kernels`, log level) are atomics but are
  process-wide; do not change them while kernels are running.
- **Logging** is serialised by a mutex.
- **Check contexts** (`CheckContext`) cache grids lazily and are not thread-safe; use one per
  thread.
- **Exceptions** thrown inside workers are captured and the first one is rethrown on the caller
  after all workers finish.

## Benchmarks

`stratonix bench --spacing 100,50 --mode both` builds the synthetic production model (40 x 40 km,
six horizons, ten faults) and reports, for each kernel and for the full build with checks, run
time and peak resident memory. Results for each release are kept outside the source tree with the
project outputs; the CHANGELOG records the headline numbers.
