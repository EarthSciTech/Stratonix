# Stratonix for Python

Python bindings over the stable C API of the Stratonix geomodelling toolkit (`ctypes`, no
compiler needed). The wheel carries the Stratonix library for its platform (Windows x64, Linux
x86-64), so every build from Python writes the same files as the `stratonix` program.

```python
import stratonix
stratonix.validate("model.yaml")
stratonix.build("model.yaml", "out/model")          # as `stratonix build`
with stratonix.Description("model.yaml") as d:
    d.set("faults.F1.throw.max", 120.0)             # labelled "set through the C API"
    with d.build() as m:
        print(m.horizons, m.evaluate([[500, 500, -900]])["unit"])
```

Licence: MIT. Documentation: the Stratonix user guide, section 2.26.
