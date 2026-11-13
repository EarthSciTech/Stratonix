# Prototype regression (Prompts 3 and 4)

The prototype section script's reference outputs are the regression target: Stratonix must
reproduce its horizon and fault-cut tables to 0.01 m (vertical-offset mode with linear block
displacement) before anything new is trusted.

**Status:** the reference files have not been supplied (register issue I-001). The test
`regression_prototype` is skipped until they are.

To run it, place in one folder:

- `prototype_description.yaml` - the prototype's section expressed as a Stratonix description
  (vertical_offset slip mode, `hanging_wall: {type: linear_between_faults}`, fixed block by the
  control well), with every parameter labelled;
- `horizons.csv` - columns `horizon,x,y,z` (z elevation, metres) sampled along the section;
- `fault_cuts.csv` - columns `fault,horizon,wall,x,y,z` (`wall` = hw or fw);

then run `ctest -R regression_prototype` with `STRATONIX_PROTOTYPE_REF=<folder>`. The column
names above are a proposal: confirm them against the prototype's actual tables when they arrive,
and adapt `tests/unit/test_regression.cpp` if they differ.
