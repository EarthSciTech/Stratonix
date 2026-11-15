#!/usr/bin/env python3
"""RESQML round trip of a Stratonix grid export (test tool only, STX-TN-12 section 6).

Usage: resqml_roundtrip.py <export folder containing grid/*.epc and grid/*.grdecl>

Reads the EPC package written by Stratonix through FESAPI with resqpy (BP, MIT), an
independent reader, and compares it with the GRDECL of the same export (read by
read_grdecl.py, also independent of Stratonix):
  - the grid extent;
  - the depth of every cell's top corners (within 0.002 m: GRDECL carries 3 decimals) and of
    the last layer's base corners;
  - the active cells;
  - the ZONE and BLOCK properties.
Then, when the RESQML 2.0.1 XML schemas are available (folder in STX_RESQML_XSD_DIR, from the
Energistics download package), every XML part of the EPC is validated against them with lxml.
Exit 0 on agreement, 1 on a mismatch, 2 when resqpy is missing.
"""
import os
import pathlib
import sys
import zipfile

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import read_grdecl  # noqa: E402


def main(folder):
    folder = pathlib.Path(folder)
    epcs = sorted((folder / "grid").glob("*.epc"))
    grdecls = sorted((folder / "grid").glob("*.grdecl"))
    if not epcs or not grdecls:
        print("no grid/*.epc and grid/*.grdecl in", folder)
        return 1
    try:
        import numpy as np
        import resqpy.model as rq
    except ImportError as e:
        print("resqpy is required:", e)
        return 2
    kw, _ = read_grdecl.read_grdecl(str(grdecls[0]))
    nx, ny, nz = (int(v) for v in kw["SPECGRID"][:3])
    model = rq.Model(str(epcs[0]))
    grid = model.grid()
    bad = []
    if tuple(grid.extent_kji) != (nz, ny, nx):
        print("extent", tuple(grid.extent_kji), "vs GRDECL", (nz, ny, nx))
        return 1
    cp = grid.corner_points(cache_cp_array=True)  # [k, j, i, kp, jp, ip, xyz]
    worst = 0.0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                X = read_grdecl.corners(kw, nx, ny, i, j, k)
                for c in range(8):
                    ib, jb, bot = c & 1, (c >> 1) & 1, (c >> 2) & 1
                    if bot and k + 1 < nz:
                        continue  # bases of inner layers: the next layer's top in RESQML
                    p = cp[k, j, i, bot, jb, ib]
                    d = max(abs(p[0] - X[c][0]), abs(p[1] - X[c][1]), abs(p[2] - (-X[c][2])))
                    worst = max(worst, d)
                    if d > 0.002 and len(bad) < 20:
                        bad.append("corner %d of cell %d %d %d differs by %.4f m" % (c, i + 1, j + 1, k + 1, d))
    pc = grid.property_collection
    for name in ("ACTIVE", "ZONE", "BLOCK"):
        a = pc.single_array_ref(citation_title=name)
        if a is None:
            bad.append("property %s missing from the EPC" % name)
            continue
        ref = np.asarray(kw["ACTNUM" if name == "ACTIVE" else name]).reshape(nz, ny, nx)
        if not np.array_equal(np.asarray(a).astype(np.int64), ref.astype(np.int64)):
            bad.append("property %s differs" % name)
    print("RESQML %s read by resqpy: %d x %d x %d cells, worst corner difference %.2e m, %d properties"
          % (epcs[0].name, nx, ny, nz, worst, len(pc.parts())))
    xsd = os.environ.get("STX_RESQML_XSD_DIR")
    if xsd and pathlib.Path(xsd).is_dir():
        from lxml import etree
        schemas = {}
        nparts = 0
        with zipfile.ZipFile(epcs[0]) as z:
            for name in z.namelist():
                if not name.endswith(".xml") or name.startswith("_rels") or "docProps" in name or name == "[Content_Types].xml":
                    continue
                doc = etree.fromstring(z.read(name))
                ns = etree.QName(doc).namespace
                main_xsd = "ResqmlAllObjects.xsd" if "resqml" in ns.lower() else "EmlAllObjects.xsd"
                path = next(pathlib.Path(xsd).rglob(main_xsd), None)
                if path is None:
                    bad.append("no %s under %s" % (main_xsd, xsd))
                    break
                if path not in schemas:
                    schemas[path] = etree.XMLSchema(etree.parse(str(path)))
                if not schemas[path].validate(doc):
                    bad.append("%s: %s" % (name, schemas[path].error_log.last_error))
                nparts += 1
        print("schema validation: %d XML parts against the RESQML 2.0.1 schemas in %s" % (nparts, xsd))
    else:
        print("schema validation not run: set STX_RESQML_XSD_DIR to the RESQML 2.0.1 schemas (Energistics download package)")
    if bad:
        print("\n".join(bad))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
