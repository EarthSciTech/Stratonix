#!/usr/bin/env python3
"""Write the well files of examples/wells_logs_ages.yaml (stage 4C, STX-TN-16 section 1).

Synthetic DESIGN data: logs follow the layering of the example (flat horizons at the depths below),
with seeded noise. W1 and W3 are LAS 2.0 (W3 wrapped), W2 is LAS 3.0 with an inclinometry section,
and W2_survey.csv is its deviation survey (minimum curvature, grid azimuths).

Usage: python3 tools/make_wells_example.py   (writes examples/data/wells/)
"""
import math
import os

import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "..", "examples", "data", "wells")
TOPS = [("H1", 800.0), ("H2", 1200.0), ("H3", 1500.0), ("H4", 1900.0)]
# per unit above H1, H1-H2 (shale), H2-H3 (sand), H3-H4 (shale), below H4 (carbonate)
GR = [55.0, 95.0, 35.0, 100.0, 20.0]
RHOB = [2.20, 2.40, 2.25, 2.50, 2.65]
DT = [120.0, 95.0, 85.0, 80.0, 60.0]


def unit_of(tvd):
    for k, (_, d) in enumerate(TOPS):
        if tvd < d:
            return k
    return len(TOPS)


def min_curvature(st):
    """Positions (north, east, tvd) by minimum curvature (Sawaryn and Thorogood 2005)."""
    out = [(0.0, 0.0, 0.0)]
    for (m1, i1, a1), (m2, i2, a2) in zip(st[:-1], st[1:]):
        i1, a1, i2, a2 = map(math.radians, (i1, a1, i2, a2))
        t1 = np.array([math.sin(i1) * math.cos(a1), math.sin(i1) * math.sin(a1), math.cos(i1)])
        t2 = np.array([math.sin(i2) * math.cos(a2), math.sin(i2) * math.sin(a2), math.cos(i2)])
        b = 2 * math.asin(min(1.0, np.linalg.norm(t2 - t1) / 2))
        rf = 1.0 if b < 1e-9 else 2 / b * math.tan(b / 2)
        d = 0.5 * (m2 - m1) * (t1 + t2) * rf
        n, e, v = out[-1]
        out.append((n + d[0], e + d[1], v + d[2]))
    return out


def tvd_at(st, pos, md):
    for k in range(1, len(st)):
        if md <= st[k][0]:
            t = (md - st[k - 1][0]) / (st[k][0] - st[k - 1][0])
            return pos[k - 1][2] + t * (pos[k][2] - pos[k - 1][2])
    return pos[-1][2] + (md - st[-1][0])


def logs(md, tvd, rng):
    u = [unit_of(t) for t in tvd]
    gr = [GR[k] + rng.normal(0, 4) for k in u]
    rhob = [RHOB[k] + rng.normal(0, 0.02) for k in u]
    dt = [DT[k] + rng.normal(0, 2) for k in u]
    return gr, rhob, dt


def las2(path, well, md, cols, wrap=False):
    names = [("DEPT", "M", "Measured depth"), ("GR", "GAPI", "Gamma ray"), ("RHOB", "G/CM3", "Bulk density"), ("DT", "US/FT", "Sonic")]
    step = md[1] - md[0]
    with open(path, "w") as f:
        f.write("~VERSION INFORMATION\n VERS.   2.0 : CWLS LOG ASCII STANDARD - VERSION 2.0\n")
        f.write(" WRAP.   %s : %s\n" % ("YES" if wrap else "NO", "Multiple lines per depth step" if wrap else "One line per depth step"))
        f.write("~WELL INFORMATION\n")
        f.write(" STRT.M %.4f : START DEPTH\n STOP.M %.4f : STOP DEPTH\n STEP.M %.4f : STEP\n NULL.  -999.25 : NULL VALUE\n" % (md[0], md[-1], step))
        f.write(" WELL.  %s : WELL\n COMP.  STRATONIX WORKED EXAMPLE : COMPANY\n FLD .  SYNTHETIC : FIELD\n" % well)
        f.write("~CURVE INFORMATION\n")
        for n, u, d in names:
            f.write(" %-4s.%-6s : %s\n" % (n, u, d))
        f.write("~PARAMETER INFORMATION\n BHT .DEGC   62.0 : Bottom hole temperature\n")
        f.write("~A  DEPT GR RHOB DT\n")
        for i in range(len(md)):
            vals = [cols[c][i] for c in range(3)]
            if wrap:
                f.write("%.2f\n %.3f %.4f\n %.2f\n" % (md[i], vals[0], vals[1], vals[2]))
            else:
                f.write("%.2f %.3f %.4f %.2f\n" % (md[i], *vals))


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(2026)
    # W1: vertical, 0 to 2300 m; a null patch in RHOB (washout)
    md = [600 + 0.5 * i for i in range(int((2300 - 600) / 0.5) + 1)]
    gr, rhob, dt = logs(md, md, rng)
    for i in range(len(md)):
        if 1340 <= md[i] <= 1346:
            rhob[i] = -999.25
    las2(os.path.join(OUT, "W1.las"), "W1", md, [gr, rhob, dt])
    # W2: kick-off at 500 m, build 3 deg/30 m to 30 deg, azimuth 060 grid
    st = [(0.0, 0.0, 60.0), (500.0, 0.0, 60.0)]
    m = 500.0
    while m < 800.0:
        m += 30.0
        st.append((m, min(30.0, (m - 500.0) / 30.0 * 3.0), 60.0))
    while m < 2400.0:
        m += 100.0
        st.append((m, 30.0, 60.0))
    with open(os.path.join(OUT, "W2_survey.csv"), "w") as f:
        f.write("# W2 deviation survey (synthetic): measured depth from the well datum, inclination, grid azimuth\n")
        f.write("md [m],inc [deg],azi [deg]\n")
        for s in st:
            f.write("%.1f,%.3f,%.3f\n" % s)
    pos = min_curvature(st)
    md2 = [600 + 1.0 * i for i in range(int((2350 - 600) / 1.0) + 1)]
    tvd2 = [tvd_at(st, pos, x) for x in md2]
    gr2, rhob2, _ = logs(md2, tvd2, rng)
    with open(os.path.join(OUT, "W2.las"), "w") as f:
        f.write("~Version\n VERS. 3.0 : CWLS LOG ASCII STANDARD - VERSION 3.0\n WRAP. NO : One line per depth step\n DLM . COMMA : Delimiting character\n")
        f.write("~Well\n STRT.M %.1f : START\n STOP.M %.1f : STOP\n STEP.M 1.0 : STEP\n NULL. -999.25 : NULL VALUE\n WELL. W2 : WELL\n" % (md2[0], md2[-1]))
        f.write("~Log_Definition\n DEPT.M : Measured depth {F}\n GR.gAPI : Gamma ray {F}\n RHOB.g/cm3 : Bulk density {F}\n")
        f.write("~Log_Data | Log_Definition\n")
        for i in range(len(md2)):
            f.write("%.1f,%.3f,%.4f\n" % (md2[i], gr2[i], rhob2[i]))
        f.write("~Inclinometry_Definition\n MD.m : Measured depth {F}\n INC.deg : Inclination {F}\n AZI.deg : Grid azimuth {F}\n")
        f.write("~Inclinometry_Data | Inclinometry_Definition\n")
        for s in st:
            f.write("%.1f,%.3f,%.3f\n" % s)
    # W3: vertical check well, wrapped LAS 2.0
    md3 = [700 + 1.0 * i for i in range(int((2200 - 700) / 1.0) + 1)]
    gr3, rhob3, dt3 = logs(md3, md3, rng)
    las2(os.path.join(OUT, "W3.las"), "W3", md3, [gr3, rhob3, dt3], wrap=True)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
