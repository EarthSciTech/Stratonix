"""Stage 4D (STX-TN-17 section 5): the Python bindings reproduce the command line's outputs and
expose the C API. Run with STRATONIX_EXE (the program) and STRATONIX_EXAMPLES (the examples)."""
import hashlib
import os
import subprocess

import pytest

import stratonix

EXE = os.environ.get("STRATONIX_EXE")
EX = os.environ.get("STRATONIX_EXAMPLES", os.path.join(os.path.dirname(__file__), "..", "..", "examples"))
DESC = os.path.join(EX, "one_fault.yaml")


def digest(folder):
    out = {}
    for root, _, files in os.walk(folder):
        for f in files:
            if f == "manifest.json":
                continue  # records the command and the time
            p = os.path.join(root, f)
            out[os.path.relpath(p, folder).replace(os.sep, "/")] = hashlib.sha256(open(p, "rb").read()).hexdigest()
    return out


def test_versions():
    assert stratonix.capi_version() >= (1, 2)
    assert stratonix.library_version()


def test_validate_and_errors(tmp_path):
    stratonix.validate(DESC)
    bad = tmp_path / "bad.yaml"
    bad.write_text("stratonix_schema: 1\n")
    with pytest.raises(stratonix.StratonixError) as e:
        stratonix.validate(str(bad))
    assert e.value.status == "STX_E_VALIDATION"


def test_model(tmp_path):
    with stratonix.build_model(DESC) as m:
        assert m.counts() == {"horizons": 3, "units": 4, "faults": 1}
        assert m.horizons == ["Top_A", "Top_B", "Top_C"] and m.faults == ["F1"]
        r = m.evaluate([[1000, 2000, -1000], [1000, 2000, -1100]])
        assert list(r["unit"]) == [1, 2]
        z = m.horizon_grid("Top_B", 0, 0, 500, 500, 3, 3)
        assert len(z) == 9


def test_description_parameters():
    with stratonix.Description(DESC) as d:
        v, u = d.get("faults.F1.throw.max")
        assert (v, u) == (30.0, "m")
        d.set("faults.F1.throw.max", 45.0)
        assert d.get("faults.F1.throw.max")[0] == 45.0
        with d.build() as m:
            assert m.counts()["faults"] == 1


@pytest.mark.skipif(not EXE, reason="STRATONIX_EXE not set")
def test_same_files_as_the_command_line(tmp_path):
    a, b = tmp_path / "python", tmp_path / "cli"
    stratonix.build(DESC, str(a))
    subprocess.run([EXE, "build", DESC, "-o", str(b)], check=True, capture_output=True)
    da, db = digest(a), digest(b)
    assert da.keys() == db.keys() and len(da) > 5
    assert da == db


def test_numpy_loaded_first():
    # a process that already has the system's C++ runtime (NumPy brings it) must not mix it with
    # the library's own: building a model in such a process crashed in 2026.4 (fixed in 2026.4.1)
    pytest.importorskip("numpy")
    code = "import numpy, stratonix; m = stratonix.build_model(%r); print(m.counts()['horizons'])" % DESC
    r = subprocess.run([os.sys.executable, "-c", code], capture_output=True, text=True, env=dict(os.environ))
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "3"
