#!/usr/bin/env python3
"""Stage 4D cross-check (STX-TN-17 section 1.7): Stratonix's experimental variogram against
GeostatsPy's gamv (MIT; the GSLIB program converted to Python), used only as a test tool, on
scattered points (no pair on a lag boundary), omnidirectional and along an azimuth.

GeostatsPy's gamv counts every pair in both orders (twice again when omnidirectional) and puts
lag k at (k - 1) x lag (its first lag is centred on 0); Stratonix counts each pair once and
centres lag k on k x lag. Semivariance and mean distance are compared lag by lag, and the pair
counts must be proportional.

Usage: geostat_crosscheck.py <stratonix executable> <work folder>
Exit 1 on any disagreement.
"""
import csv
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from geostatspy import geostats

EXE, WORK = sys.argv[1], sys.argv[2]
os.makedirs(WORK, exist_ok=True)
rng = np.random.default_rng(2026)
n = 300
x, y = rng.uniform(0, 1000, n), rng.uniform(0, 1000, n)
v = np.sin(x / 150.0) + np.cos(y / 220.0) + rng.normal(0, 0.3, n)
pts = os.path.join(WORK, "points.csv")
with open(pts, "w") as f:
    f.write("x,y,z,value\n")
    for a, b, c in zip(x, y, v):
        f.write("%.17g,%.17g,0,%.17g\n" % (a, b, c))
df = pd.DataFrame({"X": x, "Y": y, "V": v})
bad, lines = [], []
for name, azm, atol, band in [("omnidirectional", 0.0, 90.0, 1e9), ("azimuth 30, tolerance 22.5", 30.0, 22.5, 1e9),
                              ("azimuth 120, tolerance 15, bandwidth 80", 120.0, 15.0, 80.0)]:
    lag, nlag = 50.0, 10
    out = os.path.join(WORK, "stx_%d.csv" % int(azm * 10 + atol))
    args = [EXE, "variogram", pts, "--lag", str(lag), "--lags", str(nlag + 1), "--tolerance", str(lag / 2), "-o", out]
    if atol < 90:
        args += ["--azimuth", str(azm), "--angle-tolerance", str(atol), "--bandwidth", str(band)]
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode:
        bad.append("%s: Stratonix failed: %s" % (name, r.stderr.strip()))
        continue
    stx = list(csv.DictReader(open(out)))
    dis, gam, npp = geostats.gamv(df, "X", "Y", "V", -1e21, 1e21, lag, lag / 2, nlag + 1, azm, atol, band, 0)
    worst, ratios = 0.0, set()
    for k in range(1, nlag + 1):  # Stratonix lag k = GeostatsPy index k + 1
        s, g = stx[k - 1], k + 1
        if npp[g] == 0 and int(s["pairs"]) == 0:
            continue
        if int(s["pairs"]) == 0 or npp[g] == 0:
            bad.append("%s lag %d: pairs %s against %d" % (name, k, s["pairs"], npp[g]))
            continue
        ratios.add(float(round(npp[g] / int(s["pairs"]), 9)))
        worst = max(worst, abs(float(s["gamma"]) - gam[g]), abs(float(s["distance"]) - dis[g]) / lag)
    if worst > 1e-9:
        bad.append("%s: largest difference %.3g" % (name, worst))
    if len(ratios) != 1:
        bad.append("%s: pair counts not proportional (%s)" % (name, sorted(ratios)))
    lines.append("%s: %d lags, largest difference %.2e, pairs x%s in GeostatsPy" % (name, nlag, worst, sorted(ratios)))
report = "Stratonix experimental variogram against GeostatsPy gamv (stage 4D)\n" + "\n".join(lines) + "\n" + "\n".join(bad) + "\n%d failed\n" % len(bad)
open(os.path.join(WORK, "geostat_crosscheck.txt"), "w").write(report)
print(report)
sys.exit(1 if bad else 0)
