"""Stratonix from Python (stage 4D, STX-TN-17 section 5, decision 25).

A thin binding over the stable C API (``stratonix_c.h``) through ``ctypes``: no compiler, no
third-party dependency (NumPy arrays are returned when NumPy is installed, lists otherwise).
Every call goes through the same library as the ``stratonix`` program, so a build from Python
writes the same files as ``stratonix build``.

    import stratonix
    stratonix.validate("model.yaml")
    stratonix.build("model.yaml", "out/model")            # as `stratonix build`
    with stratonix.Description("model.yaml") as d:
        d.set("faults.F1.throw.max", 120.0)               # labelled "set through the C API"
        with d.build() as m:
            units = m.evaluate([[500, 500, -900]])["unit"]
"""
import ctypes
import glob
import os
import sys

__all__ = ["StratonixError", "Model", "Description", "validate", "build", "export", "package", "ensemble",
           "library_version", "capi_version", "library_path"]

_REQUIRED_CAPI = (1, 2)  # major, minimum minor


class StratonixError(RuntimeError):
    """A C API call returned a status other than STX_OK; ``status`` holds its name."""

    def __init__(self, status, message):
        super().__init__("%s: %s" % (status, message))
        self.status = status
        self.message = message


def _candidates():
    env = os.environ.get("STRATONIX_LIB")
    if env:
        yield env
    here = os.path.dirname(os.path.abspath(__file__))
    names = {"win32": ["stratonix_c.dll"], "darwin": ["libstratonix_c.dylib"]}.get(sys.platform, ["libstratonix_c.so", "libstratonix_c.so.1"])
    for d in (os.path.join(here, "_lib"), here):
        for n in names:
            yield os.path.join(d, n)
        for p in sorted(glob.glob(os.path.join(d, "*stratonix_c*"))):
            yield p


def _load():
    tried = []
    for p in _candidates():
        if not os.path.isfile(p):
            tried.append(p)
            continue
        if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
            os.add_dll_directory(os.path.dirname(p))  # the library's own dependencies beside it
        return ctypes.CDLL(p), p
    raise ImportError("the Stratonix C library was not found (set STRATONIX_LIB); tried: " + ", ".join(tried))


_lib, library_path = _load()
_c = ctypes
_s = _c.c_char_p
_lib.stx_capi_version.restype = _c.c_int
_lib.stx_library_version.restype = _s
_lib.stx_last_error.restype = _s
_lib.stx_status_name.restype = _s
_lib.stx_status_name.argtypes = [_c.c_int]
for _name, _args in {
    "stx_validate": [_s],
    "stx_model_build": [_s, _c.POINTER(_c.c_void_p)],
    "stx_model_counts": [_c.c_void_p, _c.POINTER(_c.c_int), _c.POINTER(_c.c_int), _c.POINTER(_c.c_int)],
    "stx_model_evaluate": [_c.c_void_p, _c.c_size_t, _c.POINTER(_c.c_double), _c.POINTER(_c.c_int32), _c.POINTER(_c.c_uint64), _c.POINTER(_c.c_double)],
    "stx_model_horizon_grid": [_c.c_void_p, _c.c_int, _c.c_double, _c.c_double, _c.c_double, _c.c_double, _c.c_size_t, _c.c_size_t, _c.c_int64,
                               _c.POINTER(_c.c_double)],
    "stx_build_to_folder": [_s, _s],
    "stx_export": [_s, _s, _s],
    "stx_package": [_s, _s, _s],
    "stx_description_load": [_s, _c.POINTER(_c.c_void_p)],
    "stx_description_get_parameter": [_c.c_void_p, _s, _c.POINTER(_c.c_double), _c.c_char_p, _c.c_size_t],
    "stx_description_set_parameter": [_c.c_void_p, _s, _c.c_double, _s],
    "stx_model_build_from": [_c.c_void_p, _c.POINTER(_c.c_void_p)],
    "stx_ensemble_run": [_s, _s, _s],
}.items():
    getattr(_lib, _name).argtypes = _args
    getattr(_lib, _name).restype = _c.c_int
for _name in ("stx_model_horizon_name", "stx_model_fault_name"):
    getattr(_lib, _name).argtypes = [_c.c_void_p, _c.c_int]
    getattr(_lib, _name).restype = _s
_lib.stx_model_free.argtypes = [_c.c_void_p]
_lib.stx_model_free.restype = None
_lib.stx_description_free.argtypes = [_c.c_void_p]
_lib.stx_description_free.restype = None

_v = _lib.stx_capi_version()
if _v // 100 != _REQUIRED_CAPI[0] or _v % 100 < _REQUIRED_CAPI[1]:
    raise ImportError("Stratonix C API %d.%d found at %s; these bindings need %d.%d or a later minor version"
                      % (_v // 100, _v % 100, library_path, *_REQUIRED_CAPI))


def capi_version():
    """The C API version as (major, minor)."""
    return (_v // 100, _v % 100)


def library_version():
    """The library's version ("2026.4.1")."""
    return _lib.stx_library_version().decode()


def _b(x):
    return None if x is None else os.fspath(x).encode("utf-8")


def _check(status):
    if status != 0:
        raise StratonixError(_lib.stx_status_name(status).decode(), _lib.stx_last_error().decode("utf-8", "replace"))


def _array(values, kind):
    try:
        import numpy as np
        return np.array(values, dtype=kind)
    except ImportError:
        return list(values)


def validate(description):
    """Validate a description file; raises StratonixError with the rule broken."""
    _check(_lib.stx_validate(_b(description)))


def build(description, out_dir):
    """Build, check and write every output into a new folder (as `stratonix build`)."""
    _check(_lib.stx_build_to_folder(_b(description), _b(out_dir)))


def export(description, out_dir, formats=None):
    """Export the model into a new folder; formats as `stratonix export --formats`."""
    _check(_lib.stx_export(_b(description), _b(out_dir), _b(formats)))


def package(description, out_dir, options_json=None):
    """Write a delivery package (as `stratonix package`)."""
    _check(_lib.stx_package(_b(description), _b(out_dir), _b(options_json)))


def ensemble(description, out_dir, options_json=None):
    """Run the description's ensemble (as `stratonix ensemble`)."""
    _check(_lib.stx_ensemble_run(_b(description), _b(out_dir), _b(options_json)))


class Model:
    """A built model; free it with close() or a with-block."""

    def __init__(self, handle):
        self._h = handle

    def close(self):
        if self._h:
            _lib.stx_model_free(self._h)
            self._h = None

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def __del__(self):
        self.close()

    def counts(self):
        h, u, f = _c.c_int(), _c.c_int(), _c.c_int()
        _check(_lib.stx_model_counts(self._h, _c.byref(h), _c.byref(u), _c.byref(f)))
        return {"horizons": h.value, "units": u.value, "faults": f.value}

    @property
    def horizons(self):
        return [_lib.stx_model_horizon_name(self._h, i).decode() for i in range(self.counts()["horizons"])]

    @property
    def faults(self):
        return [_lib.stx_model_fault_name(self._h, i).decode() for i in range(self.counts()["faults"])]

    def evaluate(self, points):
        """Unit, fault block and horizon scalars at points [[x, y, z], ...] (z elevation, m)."""
        pts = [float(v) for p in points for v in p]
        n = len(pts) // 3
        nh = self.counts()["horizons"]
        xyz = (_c.c_double * len(pts))(*pts)
        unit = (_c.c_int32 * n)()
        block = (_c.c_uint64 * n)()
        scal = (_c.c_double * max(1, n * nh))()
        _check(_lib.stx_model_evaluate(self._h, n, xyz, unit, block, scal))
        return {"unit": _array(unit, "int32"), "block": _array(block, "uint64"), "scalars": _array(list(scal)[: n * nh], "float64")}

    def horizon_grid(self, horizon, x0, y0, dx, dy, nx, ny, block=-1):
        """Elevations of a horizon on a grid (x fastest); NaN in fault gaps."""
        out = (_c.c_double * (nx * ny))()
        hi = horizon if isinstance(horizon, int) else self.horizons.index(horizon)
        _check(_lib.stx_model_horizon_grid(self._h, hi, x0, y0, dx, dy, nx, ny, block, out))
        return _array(out, "float64")


class Description:
    """A loaded description whose parameters can be read and set before building."""

    def __init__(self, path):
        h = _c.c_void_p()
        _check(_lib.stx_description_load(_b(path), _c.byref(h)))
        self._h = h

    def close(self):
        if self._h:
            _lib.stx_description_free(self._h)
            self._h = None

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def __del__(self):
        self.close()

    def get(self, path):
        """(value, unit) of a scalar parameter."""
        v = _c.c_double()
        u = _c.create_string_buffer(64)
        _check(_lib.stx_description_get_parameter(self._h, _b(path), _c.byref(v), u, 64))
        return v.value, u.value.decode()

    def set(self, path, value, unit=None):
        """Set a scalar parameter (labelled user_instruction, "set through the C API")."""
        _check(_lib.stx_description_set_parameter(self._h, _b(path), float(value), _b(unit)))

    def build(self):
        h = _c.c_void_p()
        _check(_lib.stx_model_build_from(self._h, _c.byref(h)))
        return Model(h)


def build_model(description):
    """Build a model from a description file and return it (free with close() or a with-block)."""
    h = _c.c_void_p()
    _check(_lib.stx_model_build(_b(description), _c.byref(h)))
    return Model(h)


__all__.append("build_model")
__version__ = library_version()
