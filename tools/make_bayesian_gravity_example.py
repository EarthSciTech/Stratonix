#!/usr/bin/env python3
"""Regenerate the gravity stations of the stage 4A example bayesian_gravity_anticline (STX-TN-14
sections 4 and 7).

Usage: make_bayesian_gravity_example.py <stratonix executable>

Builds examples/data/hybrid_fault_fold/source.yaml (fault throw 120 m, fold relief 250 m) with the
example's geophysics block, takes the modelled gravity anomaly at 49 stations on a 1.25 km
lattice at the surface, adds Gaussian noise of 0.02 mGal (seed 2026) and writes them with a sigma
column of 0.02 mGal to examples/data/bayesian_gravity/gravity.csv.
"""
import csv, pathlib, random, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "examples" / "data" / "hybrid_fault_fold" / "source.yaml"
OUT = ROOT / "examples" / "data" / "bayesian_gravity" / "gravity.csv"
L = "| user_instruction | worked example"
GEO = f"""inputs:
  - {{id: G, path: stations.csv, role: gravity stations, sharing: shareable, source: generator, evidence: interpreted}}
geophysics:
  cell: "500 m {L}"
  background_density: "2650 kg/m3 {L}"
  units:
    - {{unit: U1, density: "2350 kg/m3 {L}"}}
    - {{unit: U2, density: "2550 kg/m3 {L}"}}
  gravity: {{input: G, format: csv_stations}}
"""


def main(exe):
    xs = [1250.0 * i for i in range(1, 8)]
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        (td / "stations.csv").write_text("X,Y,Z,value\n" + "".join(f"{x:.1f},{y:.1f},0,0\n" for x in xs for y in xs))
        (td / "model.yaml").write_text(SRC.read_text().rstrip() + "\n" + GEO)
        subprocess.run([exe, "build", str(td / "model.yaml"), "-o", str(td / "out"), "--no-figures"], check=False, stdout=subprocess.DEVNULL)
        rows = list(csv.DictReader(open(td / "out" / "geophysics_stations.csv")))
    rng = random.Random(2026)
    with open(OUT, "w", newline="\n") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["X", "Y", "Z", "value", "sigma"])
        for r in rows:
            w.writerow([r["x"], r["y"], r["z"], f"{float(r['modelled']) + rng.gauss(0.0, 0.02):.4f}", "0.02"])
    print("wrote", OUT, len(rows), "stations")


if __name__ == "__main__":
    main(sys.argv[1])
