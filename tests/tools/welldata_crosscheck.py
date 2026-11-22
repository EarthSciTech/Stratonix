#!/usr/bin/env python3
"""Stage 4C cross-check (STX-TN-16 section 7): Stratonix's LAS and DLIS readers against independent
readers used only as test tools - lasio (MIT) on its public example files, and dlisio (LGPL 3.0) on
DLIS files written here by dliswriter (MIT) from known logs.

Usage: welldata_crosscheck.py <stratonix executable> <lasio examples folder> <work folder>
Exit 1 on any disagreement on a file both readers read.
"""
import csv
import glob
import math
import os
import subprocess
import sys

import numpy as np

EXE, LASDIR, WORK = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(WORK, exist_ok=True)
# files the CWLS standard does not define the way lasio reads them (reason recorded)
EXCLUDE = {
    "null_policy_ind.las": "'-1.#IND0000' is the Microsoft C runtime's text for an undefined number: "
                           "Stratonix keeps it as text (missing), lasio reads -1",
    "issue79.las": "its data follow a section title without a name ('~   DEPTH ...'): lasio reads no data, "
                   "Stratonix reads them as the log data and reports the title",
    "null_policy_runon.las": "run-on numbers ('7.020101130.188'): lasio splits them by a heuristic; "
                             "Stratonix reports the lines and does not guess where one number ends",
}


def stx_read(path, out):
    r = subprocess.run([EXE, "read-log", path, "-o", out], capture_output=True, text=True)
    if r.returncode != 0:
        return None, (r.stdout + r.stderr).strip()
    outs = sorted(glob.glob(out[:-4] + "_*.csv")) if not os.path.exists(out) else [out]
    sets = []
    for f in outs:
        with open(f, newline="") as fh:
            rows = list(csv.reader(fh))
        head = [h.split(" [")[0] for h in rows[0]]
        cols = list(zip(*rows[1:])) if len(rows) > 1 else [[] for _ in head]

        def num(x):
            try:
                return float(x)
            except ValueError:
                return math.nan
        sets.append([(h, np.array([num(x) for x in c])) for h, c in zip(head, cols)])  # duplicates kept
    return sets, ""


def same(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return a.shape == b.shape and np.allclose(a, b, rtol=1e-9, atol=1e-9, equal_nan=True)


bad = []
compared = 0
import lasio  # noqa: E402

for path in sorted(glob.glob(os.path.join(LASDIR, "*.las"))):
    name = os.path.basename(path)
    if name in EXCLUDE:
        print("skip %s: %s" % (name, EXCLUDE[name]))
        continue
    try:
        las = lasio.read(path)
    except Exception as e:  # lasio refuses: not a comparison
        print("lasio refuses %s (%s)" % (name, type(e).__name__))
        continue
    mn = [c.mnemonic for c in las.curves]
    if len(set(mn)) != len(mn) or len(mn) < 1 or any(c.data.dtype.kind not in "fiu" for c in las.curves):
        print("skip %s: duplicate or text curves" % name)
        continue
    sets, err = stx_read(path, os.path.join(WORK, name + ".csv"))
    if sets is None:
        bad.append("%s: Stratonix refused a file lasio reads: %s" % (name, err))
        continue
    # curves compared by position (lasio renames duplicate, empty and differently cased
    # mnemonics; both keep the file's order); names that differ are listed, not failed
    s = sets[0]
    for k, c in enumerate(las.curves):
        if k >= len(s):
            bad.append("%s: curve %d (%s) missing" % (name, k + 1, c.mnemonic))
            continue
        if s[k][0] != c.mnemonic:
            print("%s: curve %d named %s by Stratonix, %s by lasio" % (name, k + 1, s[k][0], c.mnemonic))
        if not same(s[k][1], c.data):
            bad.append("%s: curve %d (%s) differs" % (name, k + 1, c.mnemonic))
    compared += 1
print("LAS: %d files compared with lasio" % compared)

# DLIS written by dliswriter, read by Stratonix and dlisio
from dliswriter import DLISFile, enums  # noqa: E402
from dlisio import dlis  # noqa: E402

n = 400
depth = np.linspace(1000.0, 1199.5, n)
gr = 40 + 30 * np.sin(np.linspace(0, 12, n))
rhob = (2.2 + 0.3 * np.cos(np.linspace(0, 7, n))).astype(np.float32)
df = DLISFile()
lf = df.add_logical_file()
try:
    lf.add_origin("ORIGIN", well_name="XCHECK-1")
except TypeError:
    lf.add_origin("ORIGIN")
c1 = lf.add_channel("DEPT", data=depth, units=enums.Unit.METER)
c2 = lf.add_channel("GR", data=gr)
c3 = lf.add_channel("RHOB", data=rhob)
lf.add_frame("MAIN", channels=(c1, c2, c3), index_type=enums.FrameIndexType.BOREHOLE_DEPTH)
dpath = os.path.join(WORK, "xcheck.dlis")
df.write(dpath)
sets, err = stx_read(dpath, os.path.join(WORK, "xcheck_dlis.csv"))
if sets is None:
    bad.append("DLIS: Stratonix refused a dliswriter file: " + err)
else:
    with dlis.load(dpath) as files:
        ref = files[0].frames[0].curves()
    s = dict(sets[0])
    for nm in ("DEPT", "GR", "RHOB"):
        if nm not in s or not same(s[nm], ref[nm]):
            bad.append("DLIS: channel %s differs from dlisio" % nm)
    print("DLIS: %d frames compared with dlisio" % len(s.get("DEPT", [])))

for b in bad:
    print("MISMATCH", b)
sys.exit(1 if bad else 0)
