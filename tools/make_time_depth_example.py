#!/usr/bin/env python3
"""Write the inputs of examples/time_depth_conversion.yaml (stage 4C, STX-TN-16 section 2.3).

Synthetic DESIGN data from a known depth model and a known layer-cake velocity model
v(z) = v0 + k z: two-way-time grids of three horizons (IRAP classic ASCII, ms), the wells' tops
in depth and time, and a small SEG-Y revision 2.0 volume in time (IEEE float, big-endian) with the
reflections of the three horizons (a 25 Hz Ricker wavelet).

Usage: python3 tools/make_time_depth_example.py   (writes examples/data/time_depth/)
"""
import math
import os
import struct

import numpy as np

OUT = os.path.join(os.path.dirname(__file__), "..", "examples", "data", "time_depth")
SIZE = 4000.0


def depths(x, y):
    r2 = (x - 2000.0) ** 2 + (y - 2000.0) ** 2
    z1 = 900.0 - 150.0 * math.exp(-r2 / 1500.0 ** 2)
    z2 = z1 + 350.0 + 50.0 * x / SIZE
    z3 = z2 + 450.0
    return [z1, z2, z3]


def velocity(x, y):
    # true model: v0 of layer 2 varies laterally; k as stated in the description
    return [(1750.0, 0.6), (2300.0 + 150.0 * x / SIZE, 0.4), (2900.0, 0.25)]


def layer_time(v0, k, zt, zb):
    vt = v0 + k * zt
    return math.log1p(k * (zb - zt) / vt) / k


def twt(x, y):
    zs, vs = depths(x, y), velocity(x, y)
    t, zt, out = 0.0, 0.0, []
    for z, (v0, k) in zip(zs, vs):
        t += layer_time(v0, k, zt, z)
        zt = z
        out.append(2.0 * t)
    return out


def irap(path, x0, y0, dx, nx, ny, vals):
    with open(path, "w") as f:
        f.write("-996 %d %.6f %.6f\n" % (ny, dx, dx))
        f.write("%.6f %.6f %.6f %.6f\n" % (x0, x0 + dx * (nx - 1), y0, y0 + dx * (ny - 1)))
        f.write("%d 0.000000 %.6f %.6f\n0 0 0 0 0 0 0\n" % (nx, x0, y0))
        k = 0
        for j in range(ny):
            for i in range(nx):
                f.write("%.4f" % vals[j][i])
                k += 1
                f.write("\n" if k % 6 == 0 else " ")
        f.write("\n")


def ricker(t, f=25.0):
    a = (math.pi * f * t) ** 2
    return (1 - 2 * a) * math.exp(-a)


def segy(path, spacing=200.0, dt=0.004, tmax=1.6):
    n = int(SIZE / spacing) + 1
    ns = int(round(tmax / dt)) + 1
    text = ("C 1 STRATONIX SYNTHETIC TIME VOLUME (DESIGN), tools/make_time_depth_example.py").ljust(80)
    text += ("C 2 REFLECTIONS OF H1 H2 H3, 25 HZ RICKER, IEEE FLOAT, SEG-Y REV 2.0").ljust(80)
    text = text.ljust(3200)
    b = bytearray(text.encode("ascii"))
    bh = bytearray(400)
    struct.pack_into(">H", bh, 16, int(dt * 1e6))
    struct.pack_into(">H", bh, 20, ns)
    struct.pack_into(">H", bh, 24, 5)
    struct.pack_into(">I", bh, 96, 0x01020304)
    struct.pack_into(">H", bh, 300, 0x0200)
    b += bh
    # reflection coefficients of the three horizons (impedance contrasts, DESIGN values)
    rc = [0.12, -0.08, 0.15]
    for i in range(n):
        for j in range(n):
            x, y = j * spacing, i * spacing
            th = twt(x, y)
            tr = bytearray(240)
            struct.pack_into(">h", tr, 70, -100)
            struct.pack_into(">i", tr, 180, int(round(x * 100)))
            struct.pack_into(">i", tr, 184, int(round(y * 100)))
            struct.pack_into(">H", tr, 114, ns)
            struct.pack_into(">H", tr, 116, int(dt * 1e6))
            struct.pack_into(">i", tr, 188, 100 + i)
            struct.pack_into(">i", tr, 192, 200 + j)
            s = np.zeros(ns, dtype=">f4")
            for t0, r in zip(th, rc):
                for k in range(ns):
                    s[k] += r * ricker(k * dt - t0)
            b += tr + s.tobytes()
    with open(path, "wb") as f:
        f.write(b)


def main():
    os.makedirs(OUT, exist_ok=True)
    dx = 100.0
    n = int(SIZE / dx) + 1
    grids = [[[0.0] * n for _ in range(n)] for _ in range(3)]
    for j in range(n):
        for i in range(n):
            t = twt(i * dx, j * dx)
            for h in range(3):
                grids[h][j][i] = 1000.0 * t[h]
    for h in range(3):
        irap(os.path.join(OUT, "H%d_twt_ms.irap" % (h + 1)), 0.0, 0.0, dx, n, n, grids[h])
    with open(os.path.join(OUT, "wells_tops.txt"), "w") as f:
        f.write("# well x y: depth (m) and two-way time (ms) of H1 H2 H3 from the synthetic model\n")
        for name, (x, y) in [("W1", (2000, 2000)), ("W2", (800, 900)), ("W3", (3300, 1000)), ("W4", (1200, 3400)), ("W5", (3000, 3100))]:
            z, t = depths(x, y), twt(x, y)
            f.write("%s %d %d: %s | %s\n" % (name, x, y, " ".join("%.2f" % v for v in z), " ".join("%.2f" % (1000 * v) for v in t)))
    segy(os.path.join(OUT, "volume_twt.sgy"))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
