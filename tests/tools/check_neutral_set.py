#!/usr/bin/env python3
"""Independent check of Stratonix neutral exports (VTK XML, OBJ, STL, CSV/XYZ, GeoJSON, GemPy and
LoopStructural sets), with lxml, json, csv and numpy only; no Stratonix code.

VTK and ParaView are not installable in the build environment, so the VTK XML files are checked
against the VTK XML file-format rules (element structure, array sizes, offsets, cell types) and
cross-checked against the IRAP grids of the same export folder.
Exit 0 when all pass, 1 otherwise, 77 when a package is missing.
"""
import sys, pathlib, json, csv, io

try:
    import numpy as np
    from lxml import etree
except ImportError as e:
    print("SKIP:", e); sys.exit(77)

fails = []
def check(c, m):
    if not c: fails.append(m)
    return c

def arr(el):
    return np.array(el.text.split(), dtype=float) if el is not None and el.text else np.array([])

def irap(p):
    t = p.read_text().split()
    ny, nx = int(t[1]), int(t[8])
    dx, dy, x0, y0 = float(t[2]), float(t[3]), float(t[4]), float(t[6])
    z = np.array(t[19:19 + nx * ny], dtype=float).reshape(ny, nx)
    z[z >= 9999899] = np.nan
    return x0, y0, dx, dy, z

def check_vtu(p, ref=None):
    root = etree.parse(str(p)).getroot()
    check(root.tag == "VTKFile" and root.get("type") == "UnstructuredGrid", f"{p.name}: not a VTU file")
    piece = root.find("UnstructuredGrid/Piece")
    npts, ncell = int(piece.get("NumberOfPoints")), int(piece.get("NumberOfCells"))
    pts = arr(piece.find("Points/DataArray")).reshape(-1, 3)
    check(len(pts) == npts, f"{p.name}: point count")
    cells = {d.get("Name"): arr(d) for d in piece.find("Cells")}
    off, conn, types = cells["offsets"].astype(int), cells["connectivity"].astype(int), cells["types"].astype(int)
    check(len(off) == ncell and len(types) == ncell, f"{p.name}: cell arrays length")
    check(np.all(np.diff(off) > 0) and off[-1] == len(conn), f"{p.name}: offsets")
    check(conn.min() >= 0 and conn.max() < npts, f"{p.name}: connectivity out of range")
    check(set(types.tolist()) <= {3, 4, 5}, f"{p.name}: unexpected cell types {set(types.tolist())}")
    for d in piece.find("PointData"):
        check(len(arr(d)) == npts, f"{p.name}: point array {d.get('Name')} length")
    for d in piece.find("CellData"):
        check(len(arr(d)) == ncell, f"{p.name}: cell array {d.get('Name')} length")
    if ref is not None:
        x0, y0, dx, dy, z = ref
        i = np.rint((pts[:, 0] - x0) / dx).astype(int); j = np.rint((pts[:, 1] - y0) / dy).astype(int)
        check(np.allclose(z[j, i], pts[:, 2], atol=1e-6), f"{p.name}: points are not the IRAP nodes")
        dep = arr(piece.find("PointData/DataArray[@Name='depth_m']"))
        check(np.allclose(dep, -pts[:, 2], atol=1e-6), f"{p.name}: depth_m is not -z")
    return pts, conn, types

def main(folder):
    root = pathlib.Path(folder)
    grids = {p.name[:-len("_faulted.irap")]: irap(p) for p in (root / "horizons").glob("*_faulted.irap")}
    vtk = root / "vtk"
    for p in sorted(vtk.glob("*.vtu")):
        name = p.stem
        pts, conn, types = check_vtu(p, grids.get(name))
        print(f"ok   vtk/{p.name}: {len(pts)} points, {len(types)} cells")
    v = etree.parse(str(vtk / "voxels.vti")).getroot()
    img = v.find("ImageData")
    e = list(map(int, img.get("WholeExtent").split()))
    ncell = e[1] * e[3] * e[5]
    unit = arr(img.find("Piece/CellData/DataArray[@Name='unit']"))
    check(len(unit) == ncell, "voxels.vti: unit array length")
    # top layer of cells (k = nz-1) should be unit 0 or above in the stack; bottom the deepest units
    u3 = unit.reshape(e[5], e[3], e[1])
    check(u3[-1].mean() <= u3[0].mean(), "voxels.vti: units should increase downwards")
    print(f"ok   vtk/voxels.vti: {e[1]} x {e[3]} x {e[5]} cells, units {int(unit.min())}..{int(unit.max())}")
    pvd = etree.parse(str(vtk / "model.pvd")).getroot()
    for d in pvd.iter("DataSet"):
        check((vtk / d.get("file")).exists(), f"model.pvd: missing {d.get('file')}")
    compile((vtk / "paraview_state.py").read_text(), "paraview_state.py", "exec")
    print("ok   vtk/model.pvd and paraview_state.py (compiles)")
    for p in sorted((root / "obj").glob("*.obj")):
        vs = [l for l in p.read_text().splitlines() if l.startswith("v ")]
        fs_ = [l for l in p.read_text().splitlines() if l.startswith("f ")]
        idx = np.array([list(map(int, l.split()[1:4])) for l in fs_])
        check(idx.min() >= 1 and idx.max() <= len(vs), f"{p.name}: face index out of range")
        vtu = vtk / (p.stem + ".vtu")
        if vtu.exists():
            pts = arr(etree.parse(str(vtu)).getroot().find("UnstructuredGrid/Piece/Points/DataArray")).reshape(-1, 3)
            ov = np.array([list(map(float, l.split()[1:4])) for l in vs])
            check(np.allclose(ov, pts, atol=1e-9), f"{p.name}: vertices differ from the VTU")
        stl = (root / "stl" / (p.stem + ".stl")).read_text()
        check(stl.count("facet normal") == len(fs_), f"{p.stem}.stl: facet count differs from OBJ")
    print(f"ok   obj/stl: {len(list((root / 'obj').glob('*.obj')))} surfaces, OBJ = VTU = STL")
    for name, (x0, y0, dx, dy, z) in grids.items():
        rows = [r for r in csv.reader(io.StringIO("\n".join(l for l in (root / "csv" / (name + ".csv")).read_text().splitlines() if not l.startswith("#"))))]
        check(rows[0] == ["x", "y", "elevation_m", "depth_m", "block"], f"csv/{name}.csv header")
        a = np.array([[float(c) for c in r[:4]] for r in rows[1:]])
        check(len(a) == int((~np.isnan(z)).sum()), f"csv/{name}.csv row count")
        i = np.rint((a[:, 0] - x0) / dx).astype(int); j = np.rint((a[:, 1] - y0) / dy).astype(int)
        check(np.allclose(z[j, i], a[:, 2], atol=1e-6) and np.allclose(a[:, 3], -a[:, 2]), f"csv/{name}.csv values")
        xyz = np.loadtxt(root / "xyz" / (name + ".xyz"))
        check(np.allclose(xyz, a[:, :3]), f"xyz/{name}.xyz differs from the CSV")
    print(f"ok   csv/xyz: {len(grids)} horizons equal to the IRAP grids")
    gj = json.loads((root / "geojson" / "outlines.geojson").read_text())
    check(gj["type"] == "FeatureCollection", "geojson type")
    for f in gj["features"]:
        if f["geometry"]["type"] == "Polygon":
            ring = f["geometry"]["coordinates"][0]
            check(ring[0] == ring[-1] and len(ring) >= 4, "geojson: polygon ring not closed")
    print(f"ok   geojson: {len(gj['features'])} features, rings closed")
    for fn in ("gempy/surface_points.csv", "gempy/orientations.csv", "loopstructural/data.csv"):
        rows = list(csv.reader(io.StringIO((root / fn).read_text())))
        check(all(len(r) == len(rows[0]) for r in rows), f"{fn}: ragged rows")
        print(f"ok   {fn}: {len(rows) - 1} rows, columns {rows[0]}")
    ls = list(csv.DictReader(io.StringIO((root / "loopstructural/data.csv").read_text())))
    setup = json.loads((root / "loopstructural/model_setup.json").read_text())
    for r in ls:
        if r["feature_name"] == "strati" and r["val"]:
            check(abs(float(r["val"]) - setup["strati_values"][r["stx_horizon"]]) < 1e-3, "loopstructural: val differs from strati_values")
    for f in fails: print("FAIL", f)
    return 1 if fails else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
