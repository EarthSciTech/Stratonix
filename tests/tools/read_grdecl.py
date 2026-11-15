#!/usr/bin/env python3
"""Independent GRDECL reader (test tool only, STX-TN-12 section 6).

Reads a Stratonix build folder: exports/grid/*.grdecl and grid/cells.csv. Parses SPECGRID,
COORD, ZCORN, ACTNUM and the integer and real properties written by Stratonix, rebuilds every
cell's eight corners from its pillars, computes each cell's volume as the exact volume of the
trilinear hexahedron (Gauss quadrature, written here without Stratonix code), and compares:
  - ACTNUM with the 'active' column of cells.csv (equal);
  - the volume of every active cell with 'volume_m3' (within 1e-4 relative plus the rounding of
    the written depths, 3 decimals, over the cell's area);
  - ZONE and BLOCK with cells.csv (ZONE as written, BLOCK equal);
  - the FAULTS keyword: every face lies between cells of different BLOCK values.
Exit 0 on agreement, 1 on a mismatch, 77 (skip) when the folder has no grid.
Pure Python; no third-party packages.
"""
import csv
import glob
import math
import os
import sys


def tokens(path):
    with open(path) as f:
        for line in f:
            line = line.split("--", 1)[0]
            for t in line.split():
                yield t


def read_grdecl(path):
    kw, cur, faults, fault_line = {}, None, [], []
    for t in tokens(path):
        if cur is None:
            if t[0].isalpha():
                cur = t
                kw.setdefault(cur, [])
            continue
        if cur == "FAULTS":
            if t == "/":
                if fault_line:
                    faults.append(fault_line)
                    fault_line = []
                else:
                    cur = None
                continue
            fault_line.append(t.strip("'"))
            continue
        end = t.endswith("/")
        t = t.rstrip("/")
        if t:
            if "*" in t:
                n, v = t.split("*")
                kw[cur].extend([float(v)] * int(n))
            elif cur not in ("MAPUNITS", "GRIDUNIT"):
                try:
                    kw[cur].append(float(t))
                except ValueError:
                    kw[cur].append(t)
        if end or t == "":
            if end:
                cur = None
    return kw, faults


def corners(kw, nx, ny, i, j, k):
    coord, zc = kw["COORD"], kw["ZCORN"]
    pts = []
    for c in range(8):
        ib, jb, bot = c & 1, (c >> 1) & 1, (c >> 2) & 1
        p = (i + ib) + (nx + 1) * (j + jb)
        xt, yt, dt, xb, yb, db = coord[6 * p:6 * p + 6]
        d = zc[((2 * k + bot) * 2 * ny + (2 * j + jb)) * 2 * nx + 2 * i + ib]
        t = (d - dt) / (db - dt) if abs(db - dt) > 1e-12 else 0.0
        pts.append((xt + t * (xb - xt), yt + t * (yb - yt), -d))  # elevation
    return pts


def volume(X):
    g = 0.5 / math.sqrt(3.0)
    q = (0.5 - g, 0.5 + g)
    v = 0.0
    for a in q:
        for b in q:
            for cz in q:
                du, dv, dw = [0.0] * 3, [0.0] * 3, [0.0] * 3
                for c in range(8):
                    ib, jb, kb = c & 1, (c >> 1) & 1, (c >> 2) & 1
                    Na, Nb, Nc = (a if ib else 1 - a), (b if jb else 1 - b), (cz if kb else 1 - cz)
                    sa, sb, sc = (1 if ib else -1), (1 if jb else -1), (1 if kb else -1)
                    for d in range(3):
                        du[d] += X[c][d] * sa * Nb * Nc
                        dv[d] += X[c][d] * sb * Na * Nc
                        dw[d] += X[c][d] * sc * Na * Nb
                cr = (du[1] * dv[2] - du[2] * dv[1], du[2] * dv[0] - du[0] * dv[2], du[0] * dv[1] - du[1] * dv[0])
                v -= (cr[0] * dw[0] + cr[1] * dw[1] + cr[2] * dw[2]) * 0.125
    return v


def main(folder):
    files = glob.glob(os.path.join(folder, "exports", "grid", "*.grdecl"))
    table = os.path.join(folder, "grid", "cells.csv")
    if not files or not os.path.exists(table):
        print("skipped: no exports/grid/*.grdecl and grid/cells.csv in", folder)
        return 77
    kw, faults = read_grdecl(files[0])
    nx, ny, nz = (int(v) for v in kw["SPECGRID"][:3])
    n = nx * ny * nz
    bad = []
    for key, size in (("COORD", 6 * (nx + 1) * (ny + 1)), ("ZCORN", 8 * n), ("ACTNUM", n)):
        if len(kw.get(key, [])) != size:
            bad.append("%s has %d values, expected %d" % (key, len(kw.get(key, [])), size))
    if bad:
        print("\n".join(bad))
        return 1
    rows = list(csv.DictReader(open(table)))
    if len(rows) != n:
        print("cells.csv has %d rows, expected %d" % (len(rows), n))
        return 1
    worst = 0.0
    nact = 0
    for r in rows:
        i, j, k = int(r["i"]) - 1, int(r["j"]) - 1, int(r["k"]) - 1
        c = i + nx * (j + ny * k)
        act = int(kw["ACTNUM"][c])
        if act != int(r["active"]):
            bad.append("ACTNUM differs at cell %d %d %d" % (i + 1, j + 1, k + 1))
        if "BLOCK" in kw and int(kw["BLOCK"][c]) != int(r["block"]):
            bad.append("BLOCK differs at cell %d %d %d" % (i + 1, j + 1, k + 1))
        if not act:
            continue
        nact += 1
        X = corners(kw, nx, ny, i, j, k)
        v = volume(X)
        ref = float(r["volume_m3"])
        err = abs(v - ref)
        worst = max(worst, err / max(abs(ref), 1.0))
        # rounding of the written depths (3 decimals: +-0.0005 m on each corner) over the cell's area
        area = abs((X[1][0] - X[0][0]) * (X[2][1] - X[0][1]) - (X[2][0] - X[0][0]) * (X[1][1] - X[0][1]))
        if err > 1e-4 * abs(ref) + 0.0015 * area + 1.0:
            bad.append("volume of cell %d %d %d: %.3f (GRDECL) vs %.3f (cells.csv)" % (i + 1, j + 1, k + 1, v, ref))
    nf = 0
    for f in faults:
        name, i1, i2, j1, j2, k1, k2, face = f[0], *map(int, f[1:7]), f[7]
        for k in range(k1, k2 + 1):
            a = (i1 - 1) + nx * ((j1 - 1) + ny * (k - 1))
            b = a + 1 if face.startswith("X") else a + nx
            nf += 1
            if "BLOCK" in kw and kw["BLOCK"][a] == kw["BLOCK"][b]:
                bad.append("FAULTS face of %s at %d %d %d %s has the same block on both sides" % (name, i1, j1, k, face))
    print("GRDECL %s: %d x %d x %d cells, %d active; worst relative volume difference %.2e; %d fault faces"
          % (os.path.basename(files[0]), nx, ny, nz, nact, worst, nf))
    if bad:
        print("\n".join(bad[:20]))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
