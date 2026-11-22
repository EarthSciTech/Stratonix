#!/usr/bin/env python3
"""Stage 4C cross-check (STX-TN-16 section 7): Stratonix's SEG-Y reader against segyio (LGPL 3.0),
used only as a test tool, on files segyio writes here (IEEE and IBM floats, 2- and 4-byte
integers, both byte orders) and on the example volume.

Usage: segy_crosscheck.py <stratonix executable> <example .sgy> <work folder>
Exit 1 on any disagreement.
"""
import csv
import os
import subprocess
import sys

import numpy as np
import segyio

EXE, EXAMPLE, WORK = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(WORK, exist_ok=True)
failures, lines = 0, []


def stx(path):
    out = os.path.join(WORK, os.path.basename(path) + ".csv")
    r = subprocess.run([EXE, "read-seismic", path, "-o", out], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError((r.stdout + r.stderr).strip())
    with open(out, newline="") as fh:
        rows = list(csv.reader(fh))[1:]
    il = np.array([int(r[0]) for r in rows])
    xl = np.array([int(r[1]) for r in rows])
    xy = np.array([[float(r[2]), float(r[3])] for r in rows])
    s = np.array([[float(v) for v in r[4:]] for r in rows])
    return il, xl, xy, s


def ref(path, endian="big"):
    with segyio.open(path, ignore_geometry=True, endian=endian) as f:
        il = np.array(f.attributes(segyio.TraceField.INLINE_3D)[:])
        xl = np.array(f.attributes(segyio.TraceField.CROSSLINE_3D)[:])
        sc = np.array(f.attributes(segyio.TraceField.SourceGroupScalar)[:]).astype(float)
        cx = np.array(f.attributes(segyio.TraceField.CDP_X)[:]).astype(float)
        cy = np.array(f.attributes(segyio.TraceField.CDP_Y)[:]).astype(float)
        sx = np.array(f.attributes(segyio.TraceField.SourceX)[:]).astype(float)
        sy = np.array(f.attributes(segyio.TraceField.SourceY)[:]).astype(float)
        mul = np.where(sc > 0, sc, np.where(sc < 0, -1.0 / np.where(sc == 0, 1, sc), 1.0))
        zero = (cx == 0) & (cy == 0)
        x = np.where(zero, sx, cx) * mul
        y = np.where(zero, sy, cy) * mul
        s = np.array(f.trace.raw[:], dtype=float)
    return il, xl, np.column_stack([x, y]), s


def compare(name, path, endian="big"):
    global failures
    try:
        a = stx(path)
    except RuntimeError as e:
        failures += 1
        lines.append(f"FAIL {name}: Stratonix did not read it: {e}")
        return
    b = ref(path, endian)
    ok = (np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1]) and np.allclose(a[2], b[2], atol=1e-6)
          and a[3].shape == b[3].shape and np.allclose(a[3], b[3], rtol=1e-6, atol=1e-6))
    diff = float(np.max(np.abs(a[3] - b[3]))) if a[3].shape == b[3].shape else float("nan")
    lines.append(f"{'PASS' if ok else 'FAIL'} {name}: {len(a[0])} traces x {a[3].shape[1] if a[3].ndim == 2 else 0} samples, largest sample difference {diff:g}")
    failures += 0 if ok else 1


def write(path, fmt, endian):
    spec = segyio.spec()
    spec.format, spec.samples, spec.ilines, spec.xlines, spec.endian = fmt, list(range(50)), [10, 11, 12], [100, 101, 102, 103], endian
    rng = np.random.default_rng(fmt)
    with segyio.create(path, spec) as f:
        k = 0
        for i in spec.ilines:
            for j in spec.xlines:
                tr = rng.normal(0, 100 if fmt in (2, 3) else 1, 50)
                f.header[k] = {segyio.TraceField.INLINE_3D: i, segyio.TraceField.CROSSLINE_3D: j,
                               segyio.TraceField.SourceGroupScalar: -100,
                               segyio.TraceField.CDP_X: 50000000 + 2500 * j, segyio.TraceField.CDP_Y: 700000000 + 2500 * i}
                f.trace[k] = np.round(tr) if fmt in (2, 3) else tr.astype(np.float32)
                k += 1
        f.bin.update(tsort=segyio.TraceSortingFormat.INLINE_SORTING, hdt=4000)


for fmt, label in [(1, "IBM float"), (5, "IEEE float"), (2, "4-byte integer"), (3, "2-byte integer")]:
    for endian in ("big", "little"):
        p = os.path.join(WORK, f"segyio_f{fmt}_{endian}.sgy")
        write(p, fmt, endian)
        compare(f"segyio-written {label}, {endian}-endian", p, endian)
compare("example volume", EXAMPLE)

report = "Stratonix SEG-Y reader against segyio (stage 4C)\n" + "\n".join(lines) + f"\n{failures} failed\n"
open(os.path.join(WORK, "segy_crosscheck.txt"), "w").write(report)
print(report)
sys.exit(1 if failures else 0)
