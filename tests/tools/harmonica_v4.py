#!/usr/bin/env python3
"""Stage 4A external check (STX-TN-14 section 10): gravity and magnetics against Harmonica.

A flat three-horizon model whose unit U1 (1000 to 1200 m deep, the whole 10 x 10 km domain) has a
density contrast of +300 kg/m3 and a susceptibility of 0.02 SI against the background. Stratonix
(`stratonix build`, geophysics_stations.csv) and Harmonica (prism_gravity, prism_magnetic on the
same slab as one prism) give the anomaly at 121 stations on a 1 km lattice at the surface,
domain edges included. The largest difference is reported as a fraction of the largest anomaly;
TN-14 asks for agreement within 1 %.

Runs outside the product. Without Harmonica (--require-harmonica not given) an independent NumPy
implementation of the prism formulas stands in and the report says so.
Usage: harmonica_v4.py <stratonix> <report folder> [--require-harmonica]
"""
import pathlib, subprocess, sys, tempfile

import numpy as np

L = "| assumed | crosscheck"
YAML = f"""stratonix_schema: 1
model: {{name: harmonica_slab, version: "1.0.0", status: DESIGN}}
crs: {{name: local, vertical_datum: datum, length_unit: m}}
domain:
  origin: "[0, 0] m {L}"
  size: "[10000, 10000] m {L}"
  depth_range: "[0, 3000] m {L}"
  grid_spacing: "100 m {L}"
wells:
  - {{name: W1, role: control, head: "[5000, 5000] m {L}", tops: {{H1: "1000 m {L}", H2: "1200 m {L}", H3: "1400 m {L}"}}}}
stratigraphy:
  reference: {{well: W1}}
  horizons: [{{name: H1}}, {{name: H2}}, {{name: H3}}]
  units: [{{name: U1, top: H1, base: H2}}, {{name: U2, top: H2, base: H3}}]
structure:
  regional: {{type: flat}}
inputs:
  - {{id: G, path: stations.csv, role: gravity stations, sharing: shareable, source: crosscheck, evidence: interpreted}}
  - {{id: M, path: stations.csv, role: magnetic stations, sharing: shareable, source: crosscheck, evidence: interpreted}}
geophysics:
  cell: "250 m {L}"
  background_density: "2650 kg/m3 {L}"
  background_susceptibility: "0 SI {L}"
  inducing_field: {{inclination: "60 deg {L}", declination: "10 deg {L}", intensity: "50000 nT {L}"}}
  units:
    - {{unit: U1, density: "2950 kg/m3 {L}", susceptibility: "0.02 SI {L}"}}
  gravity: {{input: G, format: csv_stations}}
  magnetics: {{input: M, format: csv_stations}}
"""
G = 6.67430e-11


def corners(e, n, u, p):
    w, ea, s, no, b, t = p
    for i, x in enumerate((w, ea)):
        for j, y in enumerate((s, no)):
            for k, z in enumerate((b, t)):
                yield (-1) ** (i + j + k + 1), x - e, y - n, z - u


def log_ar(a, b, c, r):
    """ln(a + r), stable for a < 0; the ln(b^2 + c^2) part cancels between corners where b = c = 0."""
    q = b * b + c * c
    with np.errstate(divide="ignore", invalid="ignore"):
        neg = np.where(q > 0, np.log(np.where(q > 0, q, 1.0)), 0.0) - np.log(r - np.minimum(a, 0))
        return np.where(a >= 0, np.log(np.maximum(a, 0) + r), neg)


def np_gz(e, n, u, p, rho):
    out = 0
    for mu, x, y, z in corners(e, n, u, p):
        r = np.sqrt(x * x + y * y + z * z)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(x != 0, x * log_ar(y, x, z, r), 0) + np.where(y != 0, y * log_ar(x, y, z, r), 0) - np.where(z != 0, z * np.arctan(x * y / (z * r)), 0)
        out = out + mu * t
    return G * rho * out * 1e5


def np_dt(e, n, u, p, M, F):
    K = {k: 0 for k in ("xx", "yy", "zz", "xy", "xz", "yz")}
    for mu, x, y, z in corners(e, n, u, p):
        r = np.sqrt(x * x + y * y + z * z)
        with np.errstate(divide="ignore", invalid="ignore"):
            K["xx"] += mu * np.where(x != 0, np.arctan(y * z / (x * r)), 0)
            K["yy"] += mu * np.where(y != 0, np.arctan(x * z / (y * r)), 0)
            K["zz"] += mu * np.where(z != 0, np.arctan(x * y / (z * r)), 0)
            K["xy"] += mu * log_ar(z, x, y, r)
            K["xz"] += mu * log_ar(y, x, z, r)
            K["yz"] += mu * log_ar(x, y, z, r)
    B = [-K["xx"] * M[0] + K["xy"] * M[1] + K["xz"] * M[2], K["xy"] * M[0] - K["yy"] * M[1] + K["yz"] * M[2], K["xz"] * M[0] + K["yz"] * M[1] - K["zz"] * M[2]]
    return (B[0] * F[0] + B[1] * F[1] + B[2] * F[2]) * 1e-7 * 1e9


def main(exe, report, require):
    try:
        import harmonica as hm
        tool = "Harmonica " + hm.__version__
    except ImportError:
        if require:
            raise SystemExit("Harmonica is required (--require-harmonica) but not installed")
        hm, tool = None, "substitute: NumPy prism formulas (Harmonica not installed)"
    report = pathlib.Path(report)
    report.mkdir(parents=True, exist_ok=True)
    xs = np.arange(0, 10001, 1000.0)
    E, N = np.meshgrid(xs, xs, indexing="ij")
    E, N, U = E.ravel(), N.ravel(), np.zeros(E.size)
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        (td / "stations.csv").write_text("X,Y,Z,value,sigma\n" + "".join(f"{e:.1f},{n:.1f},0,0,1\n" for e, n in zip(E, N)))
        (td / "model.yaml").write_text(YAML)
        r = subprocess.run([exe, "build", str(td / "model.yaml"), "-o", str(td / "out")], capture_output=True, text=True)
        if not (td / "out" / "geophysics_stations.csv").exists():
            raise SystemExit("stratonix build failed:\n" + r.stdout + r.stderr)
        rows = [l.split(",") for l in (td / "out" / "geophysics_stations.csv").read_text().splitlines()[1:]]
    sx = {k: np.array([float(r[6]) for r in rows if r[0] == k]) for k in ("gravity", "magnetics")}
    slab = (0.0, 10000.0, 0.0, 10000.0, -1200.0, -1000.0)
    inc, dec, B0, k = np.radians(60), np.radians(10), 50000.0, 0.02
    F = np.array([np.cos(inc) * np.sin(dec), np.cos(inc) * np.cos(dec), -np.sin(inc)])
    Mv = F * (k * B0 * 1e-9 / (4e-7 * np.pi))
    if hm is not None:
        gz = hm.prism_gravity((E, N, U), [slab], [300.0], field="g_z")
        mag = (np.array([Mv[0]]), np.array([Mv[1]]), np.array([Mv[2]]))  # per-prism easting, northing, upward (A/m)
        be, bn, bu = hm.prism_magnetic((E, N, U), [slab], mag, field="b")
        dt = be * F[0] + bn * F[1] + bu * F[2]
    else:
        gz = np_gz(E, N, U, slab, 300.0)
        dt = np_dt(E, N, U, slab, Mv, F)
    lines = [f"Stage 4A geophysics cross-check against {tool}", f"121 stations, slab 1000 to 1200 m deep, contrast 300 kg/m3 and 0.02 SI", ""]
    ok = True
    for name, ref, unit in (("gravity", gz, "mGal"), ("magnetics", dt, "nT")):
        diff = np.abs(sx[name] - ref)
        frac = diff.max() / np.abs(ref).max()
        ok = ok and frac <= 0.01
        lines.append(f"{name}: largest {np.abs(ref).max():.4f} {unit}; largest difference {diff.max():.5f} {unit} = {100 * frac:.3f} % (limit 1 %)")
    lines.append("")
    lines.append("PASS" if ok else "FAIL")
    (report / "geophysics_crosscheck.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    if len(a) != 2:
        raise SystemExit(__doc__)
    sys.exit(main(a[0], a[1], "--require-harmonica" in sys.argv))
