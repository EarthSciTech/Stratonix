#!/usr/bin/env python3
"""Stage 4C (STX-TN-16 section 3): the ensemble of examples/time_depth_uncertainty.yaml reports the
depth variance from the conversion (picks and velocity) separately from the structure, for every
horizon of the velocity model; every member is ok; the pick error on H1 gives a conversion spread
of the size the picks imply (6 ms two-way at about 2200 m/s: 6.6 m, within 30 % at 40 members).
Usage: check_variance_split.py <ensemble output folder>
"""
import csv
import pathlib
import sys

e = pathlib.Path(sys.argv[1]) / "ensemble"
bad = []
rows = {r["horizon"]: r for r in csv.DictReader(open(e / "depth_variance_split.csv"))}
if sorted(rows) != ["H1", "H2", "H3"]:
    bad.append("horizons %s, expected H1 H2 H3" % sorted(rows))
for h, r in rows.items():
    if int(r["members"]) != 40:
        bad.append("%s: %s members, expected 40" % (h, r["members"]))
    if not float(r["conversion_variance_m2"]) > 0:
        bad.append("%s: no conversion variance" % h)
    tot = float(r["conversion_variance_m2"]) + float(r["structural_variance_m2"]) + 2 * float(r["covariance_m2"])
    if abs(tot - float(r["total_variance_m2"])) > 1e-3 * max(1.0, tot):
        bad.append("%s: the parts do not add up to the total" % h)
if "H1" in rows and abs(float(rows["H1"]["conversion_sd_m"]) - 6.6) > 0.3 * 6.6:
    bad.append("H1 conversion sd %s m, expected about 6.6 m" % rows["H1"]["conversion_sd_m"])
for r in csv.DictReader(open(e / "members.csv")):
    if r["member"] != "base" and r["status"] != "ok":
        bad.append("member %s: %s %s" % (r["member"], r["status"], r["reason"]))
print("\n".join(bad) if bad else "variance split: H1 %s m, H2 %s m, H3 %s m from the conversion" %
      tuple(rows[h]["conversion_sd_m"] for h in ("H1", "H2", "H3")))
sys.exit(1 if bad else 0)
