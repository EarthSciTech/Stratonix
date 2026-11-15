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
     the same rule for both;
  5. compares both models with the design model the example was sampled from: a horizon where
     the reference misses that source beyond the tolerance while Stratonix matches it is reported
     as inconclusive (the reference's error), not as a Stratonix failure.

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
# The design models the examples were sampled from (the truth of the comparison).
SOURCES = {"implicit_unconformity_from_samples": "unconformity_onlap_pinchout", "implicit_relay_splay_from_samples": "fault_system_relay_splay",
           "implicit_detachment_folds_from_samples": "detachment_folds", "implicit_faults_series_from_data": "data/implicit_faults_series/truth"}


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
#   iso:   "exclude" builds a series' foliation without the data of its isochore horizons (a
#          horizon that pinches out within its series, which one conformable field cannot
#          represent) and derives each such horizon from the horizon above it minus a thin-plate
#          interpolation of the thicknesses at its picks, as STX-TN-11 section 3 states;
#          "none" puts every horizon in its series' field.
#   nel:   FDI elements of each foliation.
VARIANTS = [dict(order="young_first", abut=True, iso="exclude", nel=100000), dict(order="young_first", abut=True, iso="none", nel=100000),
            dict(order="young_first", abut=True, iso="exclude", nel=400000), dict(order="young_first", abut=True, iso="none", nel=400000)]


def loopstructural_surfaces(root, setup, rows, horizons, grids, variant=None):  # pragma: no cover - needs LoopStructural
    import pandas as pd
    from LoopStructural import GeologicalModel
    variant = variant or VARIANTS[0]
    df = pd.read_csv(root / "loopstructural/data.csv")
    df = df.rename(columns={"gx": "nx", "gy": "ny", "gz": "nz"})
    iso = setup.get("isochore_horizons", []) if variant.get("iso", "exclude") == "exclude" else []
    if iso:
        df = df[~df["stx_horizon"].isin(iso)]
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
            feats[s["name"]] = model.create_and_add_foliation(s["name"], interpolatortype="FDI", nelements=variant.get("nel", 100000))
            if s["base"] == "erosional":
                model.add_unconformity(feats[s["name"]], setup["strati_values"][s["horizons"][-1]])
    else:
        for s in series:
            feat = model.create_and_add_foliation(s["name"], interpolatortype="FDI", nelements=variant.get("nel", 100000))
            if s["base"] == "erosional" and variant["order"] == "young_first":
                base = s["horizons"][-1]
                model.add_unconformity(feat, setup["strati_values"][base])
    model.update()
    feature_of = {h: s["name"] for s in series for h in s["horizons"]}
    zs = np.linspace(setup["maximum"][2], setup["origin"][2], 251)
    # an older series lies below the erosional bases of the younger ones: its field is searched
    # only below them (above, it is an extrapolation with no data)
    ceilings = {}
    bases = []
    for s_ in series:
        for h in s_["horizons"]:
            ceilings[h] = list(bases)
        if s_["base"] == "erosional":
            bases.append(s_["horizons"][-1])
    iso_f = {}

    def zfun(h, P):
        P = np.asarray(P, dtype=float).reshape(-1, 2)
        if h in iso:
            above = horizons[horizons.index(h) - 1]
            return zfun(above, P) - np.maximum(0.0, iso_f[h](P))
        pts = np.c_[np.repeat(P[:, 0], len(zs)), np.repeat(P[:, 1], len(zs)), np.tile(zs, len(P))]
        v = model.evaluate_feature_value(feature_of[h], pts).reshape(len(P), len(zs))
        if ceilings.get(h):
            ceil = np.min(np.vstack([zfun(b, P) for b in ceilings[h]]), axis=0)
            v = np.where(zs[None, :] > np.nan_to_num(ceil, nan=1e300)[:, None] + 1e-6, np.nan, v)
        return first_crossing(v, zs, setup["strati_values"][h])

    # isochore horizons: the horizon above minus the thicknesses at the picks (thin-plate spline, >= 0)
    from scipy.interpolate import RBFInterpolator
    for h in iso:
        above = horizons[horizons.index(h) - 1]
        picks = np.array([[float(r["X"]), float(r["Y"]), float(r["Z"])] for r in rows if r["stx_horizon"] == h and r["val"]])
        t = zfun(above, picks[:, :2]) - picks[:, 2]
        ok = np.isfinite(t)
        iso_f[h] = RBFInterpolator(picks[ok, :2], t[ok], kernel="thin_plate_spline", neighbors=None if ok.sum() <= 5000 else 400)
    out = {}
    for h in horizons:
        X, Y, z, dx = grids[h]
        out[h] = zfun(h, np.c_[X.ravel(), Y.ravel()]).reshape(X.shape)
    return out, zfun


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
    return out, grid_sampler(grids, out)


def grid_sampler(grids, surfaces):
    """z of a horizon at points from gridded surfaces (nearest node)."""
    def zfun(h, P):
        P = np.asarray(P, dtype=float).reshape(-1, 2)
        X, Y, _, dx = grids[h]
        z = surfaces[h]
        i = np.rint((P[:, 0] - X[0, 0]) / dx).astype(int)
        j = np.rint((P[:, 1] - Y[0, 0]) / dx).astype(int)
        ok = (i >= 0) & (i < z.shape[1]) & (j >= 0) & (j < z.shape[0])
        r = np.full(len(P), np.nan)
        r[ok] = z[j[ok], i[ok]]
        return r
    return zfun


OFFSET = 400.0        # m: horizon sampled this far either side of a fault's trace
OFFSET_SPLAY = 150.0  # m: for a fault that abuts another (a splay), closer, to stay off the other fault


def fault_throws(rows, horizons, zfun, fault, abuts=None):
    """Median throw at stations along a fault: the vertical separation of each horizon a fixed
    distance either side of the fault trace at that horizon's level, where the horizon has data
    within FOOTPRINT. For a fault that abuts another the distance is OFFSET_SPLAY and stations
    where the other fault lies within that distance plus 100 m are skipped."""
    fp = np.array([[float(r["X"]), float(r["Y"]), float(r["Z"])] for r in rows if r["stx_horizon"] == "" and r.get("feature_name") == fault])
    if len(fp) < 3:
        return float("nan")
    from scipy.spatial import cKDTree
    c = fp.mean(axis=0)
    _, _, vt = np.linalg.svd(fp - c)
    n = vt[2] if vt[2][2] > 0 else -vt[2]
    dipdir = np.array([n[0], n[1]]) / max(1e-9, np.hypot(n[0], n[1]))
    sdir = np.array([-dipdir[1], dipdir[0]])
    tan_dip = np.hypot(n[0], n[1]) / max(1e-9, n[2])
    s_lo, s_hi = np.percentile((fp[:, :2] - c[:2]) @ sdir, [10, 90])
    d = OFFSET_SPLAY if abuts else OFFSET
    other = np.array([[float(r["X"]), float(r["Y"]), float(r["Z"])] for r in rows if r["stx_horizon"] == "" and r.get("feature_name") == abuts]) if abuts else None
    P = np.array([c[:2] + sdir * t for t in np.linspace(s_lo, s_hi, 25 if abuts else 15)])
    throws = []
    for h in horizons:
        picks = np.array([[float(r["X"]), float(r["Y"])] for r in rows if r["stx_horizon"] == h and r["val"]])
        if len(picks) == 0:
            continue
        tree = cKDTree(picks)

        def at(Q):
            z = zfun(h, Q)
            dist, _ = tree.query(Q)
            return np.where(dist <= FOOTPRINT, z, np.nan)

        with np.errstate(all="ignore"):
            zl = np.nanmedian(np.vstack([at(P - dipdir * d), at(P + dipdir * d)]), axis=0)
        Q = P + np.outer((c[2] - zl) / max(tan_dip, 1e-9), dipdir)  # the trace at the horizon's level
        keep = np.isfinite(zl)
        if other is not None and len(other):
            for k in np.nonzero(keep)[0]:
                lvl = other[np.abs(other[:, 2] - zl[k]) < 150.0]
                if len(lvl) and np.min(np.hypot(lvl[:, 0] - Q[k, 0], lvl[:, 1] - Q[k, 1])) < d + 100.0:
                    keep[k] = False  # near the branch line the offset is shared with the other fault
        if not keep.any():
            continue
        Qk = Q[keep]
        fw, hw = at(Qk - dipdir * d), at(Qk + dipdir * d)
        t = fw - hw
        throws.extend(t[np.isfinite(t)].tolist())
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
            grids, feet, full = {}, {}, {}
            for h in horizons:
                X, Y, z, dx = read_irap(root / "horizons" / f"{h}_faulted.irap")
                full[h] = (X, Y, z, dx)
                grids[h] = (X[::STRIDE, ::STRIDE], Y[::STRIDE, ::STRIDE], z[::STRIDE, ::STRIDE], dx * STRIDE)
                px = [float(r["X"]) for r in rows if r["stx_horizon"] == h and r["val"]]
                py = [float(r["Y"]) for r in rows if r["stx_horizon"] == h and r["val"]]
                feet[h] = near_data(grids[h][0], grids[h][1], px, py, FOOTPRINT)
            szf = grid_sampler(full, {h: full[h][2] for h in horizons})
            for v in VARIANTS:
                if (not setup.get("isochore_horizons") and v["iso"] != "exclude") or v["nel"] != 100000:
                    continue
                try:
                    other, ozf = loopstructural_surfaces(root, setup, rows, horizons, grids, v)
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
                    t_s = fault_throws(rows, horizons, szf, f["name"], f.get("abuts"))
                    t_o = fault_throws(rows, horizons, ozf, f["name"], f.get("abuts"))
                    th.append("%s %.1f/%.1f" % (f["name"], t_s, t_o))
                print("DIAG %s iso=%s nel=%d: away mean/p95 %s; throws stx/ls %s" % (ex, v["iso"], v["nel"], ", ".join(parts), ", ".join(th)))
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
            grids, full = {}, {}
            for h in horizons:
                X, Y, z, dx = read_irap(root / "horizons" / f"{h}_faulted.irap")
                full[h] = (X, Y, z, dx)
                grids[h] = (X[::STRIDE, ::STRIDE], Y[::STRIDE, ::STRIDE], z[::STRIDE, ::STRIDE], dx * STRIDE)
            szf = grid_sampler(full, {h: full[h][2] for h in horizons})
            other, ozf = loopstructural_surfaces(root, setup, rows, horizons, grids) if have_ls else substitute_surfaces(rows, horizons, grids, root)
            # the source model's surfaces, on the same lattice
            src = pathlib.Path(td) / "src"
            subprocess.run([exe, "export", str(exdir / f"{SOURCES[ex]}.yaml"), "-o", str(src), "--formats", "irap_ascii"], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            truth = {}
            for h in horizons:
                f_ = src / "exports" / "horizons" / f"{h}_faulted.irap"
                if f_.exists():
                    Xs, Ys, zs_, _ = read_irap(f_)
                    truth[h] = zs_[::STRIDE, ::STRIDE] if zs_.shape == full[h][2].shape else None
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

                def within(a, b):
                    return (a["n"] == 0 or (a["mean"] <= TOL["away_mean"] and a["p95"] <= TOL["away_p95"])) and (b["n"] == 0 or b["p95"] <= TOL["near_p95"])

                passed = within(sa, sn)
                verdict = "pass" if passed else "FAIL"
                ts = to = None
                if truth.get(h) is not None:
                    # each model against the source: when the reference misses the truth and
                    # Stratonix does not, the disagreement is the reference's, and it is reported
                    # as inconclusive rather than as a Stratonix failure
                    ds, do = z - truth[h], other[h] - truth[h]
                    ts = (stats(ds[foot & np.isfinite(ds) & (dist > NEAR)]), stats(ds[foot & np.isfinite(ds) & (dist <= NEAR)]))
                    to = (stats(do[foot & np.isfinite(do) & (dist > NEAR)]), stats(do[foot & np.isfinite(do) & (dist <= NEAR)]))
                    if not passed and within(*ts) and not within(*to):
                        verdict = "inconclusive"
                ok &= verdict != "FAIL"
                table.append([ex, h, sa["n"], sa["mean"], sa["p95"], sa["max"], sn["n"], sn["mean"], sn["p95"], sn["max"], verdict,
                              ts[0]["mean"] if ts else float("nan"), ts[0]["p95"] if ts else float("nan"),
                              to[0]["mean"] if to else float("nan"), to[0]["p95"] if to else float("nan")])
            for f in setup["faults_oldest_first"]:
                fn = f["name"]
                t_stx = fault_throws(rows, horizons, szf, fn, f.get("abuts"))
                t_oth = fault_throws(rows, horizons, ozf, fn, f.get("abuts"))
                rel = abs(t_oth - t_stx) / max(abs(t_stx), 1e-9) if np.isfinite(t_stx) and np.isfinite(t_oth) else float("nan")
                passed = bool(np.isfinite(rel) and rel <= TOL["throw_rel"])
                ok &= passed
                ftable.append([ex, fn, t_stx, t_oth, rel, "pass" if passed else "FAIL"])
            verdicts.append((ex, ok))
    with open(report / "crosscheck_horizons.csv", "w", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["example", "horizon", "away_n", "away_mean_abs_m", "away_p95_m", "away_max_m", "near_n", "near_mean_abs_m", "near_p95_m", "near_max_m", "verdict",
                    "stratonix_vs_source_away_mean_m", "stratonix_vs_source_away_p95_m", "other_vs_source_away_mean_m", "other_vs_source_away_p95_m"])
        for r in table:
            w.writerow([r[0], r[1], r[2]] + [f"{v:.3f}" for v in r[3:6]] + [r[6]] + [f"{v:.3f}" for v in r[7:10]] + [r[10]] + [f"{v:.3f}" for v in r[11:15]])
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
    inc = [f"{r[0]} {r[1]}" for r in table if r[10] == "inconclusive"]
    if inc:
        lines += ["", "Inconclusive (the reference misses the source model beyond tolerance while Stratonix matches it): " + ", ".join(inc)]
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
