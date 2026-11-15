#!/usr/bin/env python3
"""Regenerate the data of the stage 3A example hybrid_fault_fold_conditioning (STX-TN-11 section 4).

Usage: make_hybrid_example.py <stratonix executable>

Exports the H1 and H2 picks of examples/data/hybrid_fault_fold/source.yaml (GemPy set, every 250 m),
adds Gaussian noise of 3 m (seed 2026) and writes them with a sigma column of 3 m to
examples/data/hybrid_fault_fold/picks.csv. The source's fault throw (120 m) and fold relief
(250 m) are what the example recovers.
"""
import csv, pathlib, random, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
DD = ROOT / "examples" / "data" / "hybrid_fault_fold"


def main(exe):
    with tempfile.TemporaryDirectory() as td:
        subprocess.run([exe, "export", str(DD / "source.yaml"), "-o", td + "/out", "--formats", "gempy"], check=True, stdout=subprocess.DEVNULL)
        rows = list(csv.DictReader(open(td + "/out/exports/gempy/surface_points.csv")))
    rng = random.Random(2026)
    with open(DD / "picks.csv", "w", newline="\n") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["X", "Y", "Z", "surface", "sigma"])
        for r in rows:
            if r["surface"] not in ("H1", "H2"):
                continue
            # every other pick, so the data are sparse like interpreted picks
            if (int(float(r["X"])) // 250 + int(float(r["Y"])) // 250) % 2:
                continue
            w.writerow([r["X"], r["Y"], f"{float(r['Z']) + rng.gauss(0.0, 3.0):.3f}", r["surface"], "3.0"])
    print("wrote", DD / "picks.csv")


if __name__ == "__main__":
    main(sys.argv[1])
