#!/usr/bin/env python3
"""External cross-check of the stage 3A implicit examples against LoopStructural (STX-TN-11 section 11).

Usage: crosscheck_v3.py <stratonix executable> <examples folder> <report folder> [--require-loopstructural]

For each implicit example, the script:
  1. exports it with `stratonix export <example> --formats petrel,csv,loopstructural`;
  2. builds a LoopStructural GeologicalModel from the exported input set alone (data.csv and
     model_setup.json): one fault per exported fault with the exported displacement, then one
     foliation per series (youngest first, an unconformity on each erosional base), or one
     foliation "strati" without series;
  3. compares each horizon with Stratonix's grid (exports/horizons/<h>_faulted.irap) where the
     horizon has data within 300 m: away from faults (more than 500 m from a fault gap) and near
     them;
  4. compares fault throws, measured on both models' horizons at stations along each fault, with
     the same rule for both.

Tolerances (STX-TN-11 section 11):
  - away from faults: mean |dz| <= 10 m and 95th percentile <= 25 m;
  - near faults: 95th percentile <= 50 m;
  - per fault: throws within 20 % of each other.

Without LoopStructural, a SciPy substitute is used (thin-plate RBF per horizon and fault block). It
is labelled as such and does not count as the acceptance; with --require-loopstructural the script
then fails. Runs outside the product: it reads only exported files. Also kept in 02_CODE/tools.
"""
import csv, io, json, pathlib, subprocess, sys, tempfile

import numpy as np

NEAR = 500.0       # m: near a fault = within this distance of a fault gap
STRIDE = 4         # compare on every 4th grid node (200 m on the 50 m examples)
FOOTPRINT = 300.0  # m: compare a horizon only within this distance of its data
EXAMPLES = ["implicit_unconformity_from_samples", "implicit_relay_splay_from_samples",
            "implicit_detachment_folds_from_samples", "implicit_faults_series_from_data"]
TOL = dict(away_mean=10.0, away_p95=25.0, near_p95=50.0, throw_rel=0.20)


def read_irap(p):
    t = p.read_text().split()
    ny, nx = int(t[1]), int(t[8])
    dx, dy, x0, y0 = float(t[2]), float(t[3]), float(t[4]), float(t[6])
    z = np.array(t[19:19 + nx * ny], dtype=float).reshape(ny, nx)
    z[z >= 9999899] = np.nan
    X, Y = np.meshgrid(x0 + dx * np.arange(nx), y0 + dy * np.arange(ny))
    return X, Y, z, dx


def distance_to_gaps(z, dx):
    from scipy.ndimage import distance_transform_edt
    gaps = np.isnan(z)
    return np.full(z.shape, np.inf) if not gaps.any() else distance_transform_edt(~gaps) * dx


def near_data(X, Y, px, py, r):
    from scipy.spatial import cKDTree
    if len(px) == 0:
        return np.zeros(X.shape, dtype=bool)
    d, _ = cKDTree(np.c_[px, py]).query(np.c_[X.ravel(), Y.ravel()])
    return (d <= r).reshape(X.shape)


def stats(d):
    d = np.abs(d[np.isfinite(d)])
    if len(d) == 0:
        return dict(n=0, mean=float("nan"), p95=float("nan"), max=float("nan"))
    return dict(n=int(len(d)), mean=float(d.mean()), p95=float(np.percentile(d, 95)), max=float(d.max()))


def first_crossing(values, zs, target):
    """Elevation of the first downward crossing of `target` per column (values: columns x levels)."""
    out = np.full(values.shape[0], np.nan)
    for k in range(len(zs) - 1):
        a, b = values[:, k] - target, values[:, k + 1] - target
        m = np.isnan(out) & np.isfinite(a) & np.isfinite(b) & (a >= 0) & (b < 0)
        t = a[m] / (a[m] - b[m])
        out[m] = zs[k] + t * (zs[k + 1] - zs[k])
    return out


# LoopStructural set-ups tried by --diagnose (the first is the acceptance set-up):
#   order: "young_first" creates the series youngest first, each erosional base added as an
#          unconformity right after its series; "old_first" creates the older series first and
#          adds the unconformity after the younger one; "none" adds no unconformity.
#   abut:  apply the exported fault abutting relations (add_abutting_fault).
VARIANTS = [dict(order="young_first", abut=True), dict(order="old_first", abut=True), dict(order="none", abut=True),
            dict(order="young_first", abut=False)]


def loopstructural_surfaces(root, setup, rows, horizons, grids, variant=None):  # pragma: no cover - needs LoopStructural
    import pandas as pd
    from LoopStructural import GeologicalModel
    variant = variant or VARIANTS[0]
    df = pd.read_csv(root / "loopstructural/data.csv")
    df = df.rename(columns={"gx": "nx", "gy": "ny", "gz": "nz"})
    df = df.drop(columns=[c for c in df.columns if c.startswith("stx_")])
    model = GeologicalModel(setup["origin"], setup["maximum"])
    model.data = df
    faults = {}
    for f in setup["faults_oldest_first"]:
        faults[f["name"]] = model.create_and_add_fault(f["name"], f["displacement_max_m"], nelements=20000, interpolatortype="FDI")
        if variant["abut"] and f.get("abuts") in faults:
            try:
                faults[f["name"]].add_abutting_fault(faults[f["abuts"]])
            except Exception as e:  # noqa: BLE001 - report and go on without the relation
                print("  abutting %s on %s not applied: %s" % (f["name"], f["abuts"], e))
    series = setup.get("series_youngest_first") or [{"name": "strati", "base": "conformable", "horizons": horizons}]
    if variant["order"] == "old_first":
        feats = {}
        for s in reversed(series):
            feats[s["name"]] = model.create_and_add_foliation(s["name"], interpolatortype="FDI", nelements=100000)
            if s["base"] == "erosional":
                model.add_unconformity(feats[s["name"]], setup["strati_values"][s["horizons"][-1]])
    else:
        for s in series:
            feat = model.create_and_add_foliation(s["name"], interpolatortype="FDI", nelements=100000)
            if s["base"] == "erosional" and variant["order"] == "young_first":
                base = s["horizons"][-1]
                model.add_unconformity(feat, setup["strati_values"][base])
    model.update()
    feature_of = {h: s["name"] for s in series for h in s["horizons"]}
    out = {}
    for h in horizons:
        X, Y, z, dx = grids[h]
        zs = np.linspace(setup["maximum"][2], setup["origin"][2], 251)
        pts = np.c_[np.repeat(X.ravel(), len(zs)), np.repeat(Y.ravel(), len(zs)), np.tile(zs, X.size)]
        v = model.evaluate_feature_value(feature_of[h], pts).reshape(X.size, len(zs))
        out[h] = first_crossing(v, zs, setup["strati_values"][h]).reshape(X.shape)
    return out


def substitute_surfaces(rows, horizons, grids, root):
    """SciPy thin-plate RBF per horizon and per Stratonix fault block (not LoopStructural)."""
    from scipy.interpolate import RBFInterpolator
    out = {}
    for h in horizons:
        X, Y, z, dx = grids[h]
        pts = np.array([[float(r["X"]), float(r["Y"]), float(r["Z"]), int(r["stx_block"] or -1)] for r in rows
                        if r["stx_horizon"] == h and r["val"]])
        blk = np.full(X.shape, -1, dtype=np.int64)
        for r in csv.DictReader(io.StringIO("\n".join(l for l in (root / "csv" / f"{h}.csv").read_text().splitlines() if not l.startswith("#")))):
            i, j = int(round((float(r["x"]) - X[0, 0]) / dx)), int(round((float(r["y"]) - Y[0, 0]) / dx))
            blk[j, i] = int(r["block"])
        s = np.full(X.shape, np.nan)
        for b in np.unique(pts[:, 3]).astype(int):
            sel = pts[:, 3] == b
            if sel.sum() < 4:
                continue
            nb = None if sel.sum() <= 5000 else 400  # global fit where it fits in memory
            m = blk == b
            try:
                f = RBFInterpolator(pts[sel, :2], pts[sel, 2], kernel="thin_plate_spline", neighbors=nb)
                s[m] = f(np.c_[X[m], Y[m]])
            except np.linalg.LinAlgError:  # picks of a small block in a line: nearest pick
                from scipy.interpolate import NearestNDInterpolator
                s[m] = NearestNDInterpolator(pts[sel, :2], pts[sel, 2])(np.c_[X[m], Y[m]])
        out[h] = s
    return out


def fault_throws(rows, horizons, surfaces, fault, feet):
    """Median throw at stations along a fault, measured on the given horizon surfaces where the
    horizon has data (feet): the vertical separation of each horizon 400 m either side of the
    fault trace at that horizon's level."""
    fp = np.array([[float(r["X"]), float(r["Y"]), float(r["Z"])] for r in rows if r["stx_horizon"] == "" and r.get("feature_name") == fault])
    if len(fp) < 3:
        return float("nan")
    c = fp.mean(axis=0)
    _, _, vt = np.linalg.svd(fp - c)
    n = vt[2] if vt[2][2] > 0 else -vt[2]
    dipdir = np.array([n[0], n[1]]) / max(1e-9, np.hypot(n[0], n[1]))
    sdir = np.array([-dipdir[1], dipdir[0]])
    tan_dip = np.hypot(n[0], n[1]) / max(1e-9, n[2])
    s_lo, s_hi = np.percentile((fp[:, :2] - c[:2]) @ sdir, [10, 90])
    throws = []
    for h in horizons:
        X, Y, z = surfaces[h]
        x0, y0, dx = X[0, 0], Y[0, 0], X[0, 1] - X[0, 0]

        foot = feet[h]

        def at(x, y):
            i, j = int(round((x - x0) / dx)), int(round((y - y0) / dx))
            return z[j, i] if 0 <= i < z.shape[1] and 0 <= j < z.shape[0] and foot[j, i] else np.nan

        for s in np.linspace(s_lo, s_hi, 15):
            p = c[:2] + sdir * s
            zl = np.nanmedian([at(*(p - dipdir * 400)), at(*(p + dipdir * 400))])
            if not np.isfinite(zl):
                continue
            q = p + dipdir * (c[2] - zl) / max(tan_dip, 1e-9)  # the trace at the horizon's level
            fw, hw = at(*(q - dipdir * 400)), at(*(q + dipdir * 400))
            if np.isfinite(fw) and np.isfinite(hw):
                throws.append(fw - hw)
    return float(np.median(throws)) if throws else float("nan")


def diagnose(exe, exdir):  # pragma: no cover - needs LoopStructural
    """Every LoopStructural set-up in VARIANTS on every example: mean |dz| away from faults per
    horizon and throws, printed for choosing the set-up."""
    for ex in EXAMPLES:
        with tempfile.TemporaryDirectory() as td:
            out = pathlib.Path(td) / "out"
            subprocess.run([exe, "export", str(exdir / f"{ex}.yaml"), "-o", str(out), "--formats", "petrel,csv,loopstructural"], check=True,
                           stdout=subprocess.DEVNULL)
            root = out / "exports"
            setup = json.loads((root / "loopstructural/model_setup.json").read_text())
            rows = list(csv.DictReader(io.StringIO((root / "loopstructural/data.csv").read_text())))
            horizons = [h for h in setup["strati_values"]]
            grids, feet = {}, {}
            for h in horizons:
                X, Y, z, dx = read_irap(root / "horizons" / f"{h}_faulted.irap")
                grids[h] = (X[::STRIDE, ::STRIDE], Y[::STRIDE, ::STRIDE], z[::STRIDE, ::STRIDE], dx * STRIDE)
                px = [float(r["X"]) for r in rows if r["stx_horizon"] == h and r["val"]]
                py = [float(r["Y"]) for r in rows if r["stx_horizon"] == h and r["val"]]
                feet[h] = near_data(grids[h][0], grids[h][1], px, py, FOOTPRINT)
            for v in VARIANTS:
                if not setup.get("series_youngest_first") and v["order"] != "young_first":
                    continue
                if not any(f.get("abuts") for f in setup["faults_oldest_first"]) and not v["abut"]:
                    continue
                try:
                    other = loopstructural_surfaces(root, setup, rows, horizons, grids, v)
                except Exception as e:  # noqa: BLE001
                    print("DIAG %s %s: failed: %s" % (ex, v, e))
                    continue
                parts = []
                for h in horizons:
                    X, Y, z, dx = grids[h]
                    d = other[h] - z
                    m = feet[h] & np.isfinite(d) & (distance_to_gaps(z, dx) > NEAR)
                    parts.append("%s %.1f/%.1f" % (h, stats(d[m])["mean"], stats(d[m])["p95"]))
                th = []
                for f in setup["faults_oldest_first"]:
                    t_s = fault_throws(rows, horizons, {h: grids[h][:3] for h in horizons}, f["name"], feet)
                    t_o = fault_throws(rows, horizons, {h: (grids[h][0], grids[h][1], other[h]) for h in horizons}, f["name"], feet)
                    th.append("%s %.1f/%.1f" % (f["name"], t_s, t_o))
                print("DIAG %s order=%s abut=%s: away mean/p95 %s; throws stx/ls %s" % (ex, v["order"], v["abut"], ", ".join(parts), ", ".join(th)))
    return 0


def main(argv):
    exe, exdir, report = argv[1], pathlib.Path(argv[2]), pathlib.Path(argv[3])
    require = "--require-loopstructural" in argv
    report.mkdir(parents=True, exist_ok=True)
    try:
        import LoopStructural
        engine = "LoopStructural " + LoopStructural.__version__
        have_ls = True
    except ImportError:
        engine = "SciPy RBFInterpolator (thin-plate spline, per fault block) - substitute: LoopStructural not installed"
        have_ls = False
        if require:
            print("LoopStructural is required (--require-loopstructural) but cannot be imported")
            return 2
    if "--diagnose" in argv and have_ls:
        return diagnose(exe, exdir)
    table, ftable, verdicts = [], [], []
    for ex in EXAMPLES:
        with tempfile.TemporaryDirectory() as td:
            out = pathlib.Path(td) / "out"
            subprocess.run([exe, "export", str(exdir / f"{ex}.yaml"), "-o", str(out), "--formats", "petrel,csv,loopstructural"], check=True,
                           stdout=subprocess.DEVNULL)
            root = out / "exports"
            setup = json.loads((root / "loopstructural/model_setup.json").read_text())
            rows = list(csv.DictReader(io.StringIO((root / "loopstructural/data.csv").read_text())))
            horizons = [h for h in setup["strati_values"]]
            grids = {}
            for h in horizons:
                X, Y, z, dx = read_irap(root / "horizons" / f"{h}_faulted.irap")
                grids[h] = (X[::STRIDE, ::STRIDE], Y[::STRIDE, ::STRIDE], z[::STRIDE, ::STRIDE], dx * STRIDE)
            other = loopstructural_surfaces(root, setup, rows, horizons, grids) if have_ls else substitute_surfaces(rows, horizons, grids, root)
            ok = True
            feet = {}
            for h in horizons:
                X, Y, z, dx = grids[h]
                px = [float(r["X"]) for r in rows if r["stx_horizon"] == h and r["val"]]
                py = [float(r["Y"]) for r in rows if r["stx_horizon"] == h and r["val"]]
                foot = near_data(X, Y, px, py, FOOTPRINT)
                feet[h] = foot
                d = other[h] - z
                dist = distance_to_gaps(z, dx)
                near, away = foot & np.isfinite(d) & (dist <= NEAR), foot & np.isfinite(d) & (dist > NEAR)
                sa, sn = stats(d[away]), stats(d[near])
                passed = (sa["n"] == 0 or (sa["mean"] <= TOL["away_mean"] and sa["p95"] <= TOL["away_p95"])) and (sn["n"] == 0 or sn["p95"] <= TOL["near_p95"])
                ok &= passed
                table.append([ex, h, sa["n"], sa["mean"], sa["p95"], sa["max"], sn["n"], sn["mean"], sn["p95"], sn["max"], "pass" if passed else "FAIL"])
            for f in setup["faults_oldest_first"]:
                fn = f["name"]
                t_stx = fault_throws(rows, horizons, {h: grids[h][:3] for h in horizons}, fn, feet)
                t_oth = fault_throws(rows, horizons, {h: (grids[h][0], grids[h][1], other[h]) for h in horizons}, fn, feet)
                rel = abs(t_oth - t_stx) / max(abs(t_stx), 1e-9) if np.isfinite(t_stx) and np.isfinite(t_oth) else float("nan")
                passed = bool(np.isfinite(rel) and rel <= TOL["throw_rel"])
                ok &= passed
                ftable.append([ex, fn, t_stx, t_oth, rel, "pass" if passed else "FAIL"])
            verdicts.append((ex, ok))
    with open(report / "crosscheck_horizons.csv", "w", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["example", "horizon", "away_n", "away_mean_abs_m", "away_p95_m", "away_max_m", "near_n", "near_mean_abs_m", "near_p95_m", "near_max_m", "verdict"])
        for r in table:
            w.writerow([r[0], r[1], r[2]] + [f"{v:.3f}" for v in r[3:6]] + [r[6]] + [f"{v:.3f}" for v in r[7:10]] + [r[10]])
    with open(report / "crosscheck_faults.csv", "w", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["example", "fault", "throw_stratonix_m", "throw_other_m", "relative_difference", "verdict"])
        for r in ftable:
            w.writerow([r[0], r[1], f"{r[2]:.2f}", f"{r[3]:.2f}", f"{r[4]:.3f}", r[5]])
    allok = all(ok for _, ok in verdicts)
    lines = [f"Stratonix stage 3A external cross-check (STX-TN-11 section 11)", f"Interpolator: {engine}",
             "Acceptance: " + ("yes (LoopStructural)" if have_ls else "NO - substitute interpolator, not the acceptance"),
             f"Tolerances: away mean <= {TOL['away_mean']} m and p95 <= {TOL['away_p95']} m; near p95 <= {TOL['near_p95']} m; throws within {TOL['throw_rel']:.0%}", ""]
    lines += [f"{ex}: {'within tolerance' if ok else 'OUTSIDE tolerance'}" for ex, ok in verdicts]
    (report / "crosscheck_report.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    # the tables too, so the log carries the numbers
    print("\n" + (report / "crosscheck_horizons.csv").read_text())
    print((report / "crosscheck_faults.csv").read_text())
    if have_ls:
        return 0 if allok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
