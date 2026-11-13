#!/usr/bin/env python3
"""Independent check of a Stratonix Petrel export folder (own parsers, numpy only; no Stratonix code).

xtgeo is not available in the build environment, so the formats are parsed here from their
published layouts and cross-checked against each other:
  - IRAP classic ASCII and ZMAP+ grids of each horizon hold the same nodes and the same gaps;
  - TSurf horizon vertices are IRAP nodes; fault TSurf vertices are coplanar for planar faults;
  - fault sticks (IRAP lines and Petrel points) lie on their fault's TSurf plane;
  - fault polygons are closed and lie inside the fault gaps of their horizon grid;
  - well tops lie on the IRAP grid of their horizon at the well (where the node cell is defined);
  - every file is ASCII with LF line endings.
Exit 0 when all pass, 1 otherwise, 77 when numpy is missing.
"""
import sys, pathlib, re

try:
    import numpy as np
except ImportError:
    print("SKIP: numpy missing"); sys.exit(77)

UNDEF_IRAP, NULL_ZMAP, SEP = 9999900.0, -99999.0, 999.0
fails = []

def check(cond, msg):
    if not cond:
        fails.append(msg)
    return cond

def read_irap(p):
    t = p.read_text().split()
    assert t[0] == "-996", p
    ny, dx, dy = int(t[1]), float(t[2]), float(t[3])
    xmin, xmax, ymin, ymax = map(float, t[4:8])
    nx = int(t[8])
    v = np.array(list(map(float, t[19:19 + nx * ny]))).reshape(ny, nx)  # x fastest from (xmin, ymin)
    v[v >= UNDEF_IRAP - 1] = np.nan
    return dict(nx=nx, ny=ny, x0=xmin, y0=ymin, dx=dx, dy=dy, z=v)

def read_zmap(p):
    lines = [l for l in p.read_text().splitlines() if l.strip() and not l.startswith("!")]
    ats = [i for i, l in enumerate(lines) if l.startswith("@")]
    h1 = [s.strip() for s in lines[ats[0] + 1].split(",")]
    h2 = [s.strip() for s in lines[ats[0] + 2].split(",")]
    null = float(h1[1]); ny, nx = int(h2[0]), int(h2[1])
    vals = np.array(" ".join(lines[ats[1] + 1:]).split(), dtype=float)
    z = vals[:nx * ny].reshape(nx, ny)[:, ::-1].T  # columns from the top down -> [row from ymin, col]
    z[np.abs(z - null) < 1e-6] = np.nan
    return dict(nx=nx, ny=ny, x0=float(h2[2]), y0=float(h2[4]), z=z)

def read_lines(p):
    out, cur = [], []
    for l in p.read_text().splitlines():
        f = l.split()
        if len(f) < 3: continue
        x, y, z = map(float, f[:3])
        if x == SEP and y == SEP and z == SEP:
            if cur: out.append(np.array(cur)); cur = []
        else:
            cur.append((x, y, z))
    if cur: out.append(np.array(cur))
    return out

def read_tsurf(p):
    vid, v, t = {}, [], []
    for l in p.read_text().splitlines():
        f = l.split()
        if not f: continue
        if f[0] in ("VRTX", "PVRTX"):
            vid[int(f[1])] = len(v); v.append(tuple(map(float, f[2:5])))
        elif f[0] == "TRGL":
            t.append(tuple(vid[int(k)] for k in f[1:4]))
    return np.array(v), np.array(t, dtype=int)

def quoted_fields(l):
    return [a or b for a, b in re.findall(r'"([^"]*)"|(\S+)', l)]

def read_table(p):
    rows, data = [], False
    for l in p.read_text().splitlines():
        if not l.strip() or l.startswith("#"): continue
        if l.strip() == "END HEADER": data = True; continue
        if data: rows.append(quoted_fields(l))
    return rows

def plane_fit(pts):
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c)
    return c, vt[2]

def point_in_polygon(X, Y, px, py):
    inside = np.zeros(X.shape, bool)
    n = len(px)
    for k in range(n - 1):
        x1, y1, x2, y2 = px[k], py[k], px[k + 1], py[k + 1]
        cond = (y1 > Y) != (y2 > Y)
        with np.errstate(divide="ignore", invalid="ignore"):
            xc = x1 + (Y - y1) * (x2 - x1) / (y2 - y1)
        inside ^= cond & (X < xc)
    return inside

def distance_to_polyline(X, Y, px, py):
    d = np.full(X.shape, np.inf)
    for k in range(len(px) - 1):
        ax, ay, bx, by = px[k], py[k], px[k + 1], py[k + 1]
        vx, vy = bx - ax, by - ay
        L = vx * vx + vy * vy
        t = np.clip(((X - ax) * vx + (Y - ay) * vy) / L, 0, 1) if L > 0 else 0
        d = np.minimum(d, np.hypot(X - (ax + t * vx), Y - (ay + t * vy)))
    return d

def bilinear(g, x, y):
    fx, fy = (x - g["x0"]) / g["dx"], (y - g["y0"]) / g["dy"]
    i, j = int(np.floor(fx)), int(np.floor(fy))
    i, j = min(max(i, 0), g["nx"] - 2), min(max(j, 0), g["ny"] - 2)
    u, w = fx - i, fy - j
    z = g["z"]
    c = [z[j, i], z[j, i + 1], z[j + 1, i], z[j + 1, i + 1]]
    if any(np.isnan(c)): return np.nan
    return (1 - u) * (1 - w) * c[0] + u * (1 - w) * c[1] + (1 - u) * w * c[2] + u * w * c[3]

def main(folder):
    root = pathlib.Path(folder)
    for f in root.rglob("*"):
        if f.is_file():
            b = f.read_bytes()
            check(b"\r" not in b, f"{f.name}: CR line endings")
            check(all(c < 128 for c in b), f"{f.name}: non-ASCII bytes")
    grids = {}
    for p in sorted((root / "horizons").glob("*_faulted.irap")):
        name = p.name[:-len("_faulted.irap")]
        a, z = read_irap(p), read_zmap(p.with_name(name + "_faulted.zmap"))
        grids[name] = a
        check((a["nx"], a["ny"]) == (z["nx"], z["ny"]), f"{name}: IRAP/ZMAP size")
        same_gap = np.array_equal(np.isnan(a["z"]), np.isnan(z["z"]))
        check(same_gap, f"{name}: IRAP/ZMAP gaps differ")
        if same_gap:
            d = np.nanmax(np.abs(a["z"] - z["z"]))
            check(d < 1e-6, f"{name}: IRAP/ZMAP values differ by {d}")
        vts, tri = read_tsurf(root / "horizons" / (name + "_faulted.ts"))
        i = np.rint((vts[:, 0] - a["x0"]) / a["dx"]).astype(int); j = np.rint((vts[:, 1] - a["y0"]) / a["dy"]).astype(int)
        check(np.allclose(a["z"][j, i], vts[:, 2], atol=1e-6), f"{name}: TSurf vertices are not IRAP nodes")
        check(np.nanmax(a["z"]) < 0, f"{name}: elevations should be negative (z is elevation)")
        print(f"ok   {name}: {a['nx']}x{a['ny']} nodes, {int(np.isnan(a['z']).sum())} gap nodes, {len(tri)} TSurf triangles")
    planes = {}
    for p in sorted((root / "faults").glob("*.ts")):
        v, t = read_tsurf(p)
        c, n = plane_fit(v)
        planes[p.stem] = (c, n)
        dev = np.abs((v - c) @ n).max()
        check(dev < 0.05, f"{p.stem}: TSurf vertices not coplanar ({dev:.3f} m)")
        sticks = read_lines(root / "faults" / (p.stem + "_sticks.dat"))
        sd = max(np.abs((s - c) @ n).max() for s in sticks)
        check(sd < 0.05, f"{p.stem}: sticks off the TSurf plane ({sd:.3f} m)")
        print(f"ok   {p.stem}: {len(t)} triangles, {len(sticks)} sticks, plane deviation {dev:.4f} m, sticks {sd:.4f} m")
    pts = read_table(root / "faults" / "fault_sticks_points.txt")
    for r in pts:
        c, n = planes[r[3]]
        check(abs((np.array(list(map(float, r[:3]))) - c) @ n) < 0.05, f"points: {r} off fault {r[3]}")
    print(f"ok   fault_sticks_points: {len(pts)} points on their faults")
    for p in sorted((root / "horizons" / "fault_polygons").glob("*.dat")):
        name = p.name[:-len("_fault_polygons.dat")]
        g = grids[name]
        polys = read_lines(p)
        for poly in polys:
            check(np.allclose(poly[0], poly[-1]), f"{name}: polygon not closed")
        # Grid nodes strictly inside a fault polygon must be gaps (the property Petrel relies on), and
        # gap nodes must lie inside or next to (within one cell of) a polygon.
        X, Y = np.meshgrid(g["x0"] + g["dx"] * np.arange(g["nx"]), g["y0"] + g["dy"] * np.arange(g["ny"]))
        inside = np.zeros(X.shape, bool)
        near = np.zeros(X.shape, bool)
        for poly in polys:
            inside |= point_in_polygon(X, Y, poly[:, 0], poly[:, 1])
            near |= distance_to_polyline(X, Y, poly[:, 0], poly[:, 1]) <= 1.0 * max(g["dx"], g["dy"])
        valued_inside = int((inside & ~np.isnan(g["z"])).sum())
        check(valued_inside == 0, f"{name}: {valued_inside} grid nodes inside fault polygons carry values")
        gaps = np.isnan(g["z"])
        stray = int((gaps & ~(inside | near)).sum())
        check(stray == 0, f"{name}: {stray} gap nodes away from every fault polygon")
        print(f"ok   {name}: {len(polys)} fault polygons; {int((inside & gaps).sum())} gap nodes inside, none with values; no stray gaps")
    tops = read_table(root / "wells" / "well_tops.txt")
    worst = 0.0
    for r in tops:
        x, y, z = map(float, r[:3])
        zg = bilinear(grids[r[5]], x, y)
        if not np.isnan(zg):
            worst = max(worst, abs(zg - z))
    tol = max(1.0, 0.02 * max(g["dx"] for g in grids.values()))  # bilinear interpolation of a curved surface
    check(worst < tol, f"well tops differ from grids by {worst:.3f} m (tolerance {tol:.2f} m)")
    print(f"ok   well tops: {len(tops)} tops, max difference to the grids {worst:.3f} m (grid interpolation)")
    heads = read_table(root / "wells" / "well_heads.txt")
    print(f"ok   well heads: {len(heads)} wells")
    for f in fails: print("FAIL", f)
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
