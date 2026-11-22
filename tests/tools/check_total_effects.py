#!/usr/bin/env python3
"""Stage 4D (STX-TN-17 section 3): the Saltelli design of examples/sensitivity_one_fault.yaml.
"Max throw F1" is the fault's maximum throw, so the throw's total effect is 1 within its interval
and the other parameters' indices are 0 (within 0.05); every row states the number of runs.
Usage: check_total_effects.py <ensemble output folder>
"""
import csv
import pathlib
import sys

rows = [r for r in csv.DictReader(l for l in open(pathlib.Path(sys.argv[1]) / "ensemble" / "sensitivity_total.csv") if not l.startswith("#"))]
bad = []
for r in rows:
    if r["runs"] != "320":
        bad.append("%s: %s runs, expected 320" % (r["parameter"], r["runs"]))
    s1, st = float(r["first_order"]), float(r["total_effect"])
    if r["parameter"] == "faults.F1.throw.max":
        if not float(r["total_effect_low"]) <= 1.0 <= float(r["total_effect_high"]):
            bad.append("throw: total effect interval [%s, %s] misses 1" % (r["total_effect_low"], r["total_effect_high"]))
    elif abs(s1) > 0.05 or abs(st) > 0.05:
        bad.append("%s: indices %.3f, %.3f, expected 0" % (r["parameter"], s1, st))
print("\n".join(bad) if bad else "total effects: %d parameters as expected" % len(rows))
sys.exit(1 if bad or len(rows) != 3 else 0)
