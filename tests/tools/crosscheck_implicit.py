#!/usr/bin/env python3
"""External cross-check: rebuild the horizons of a Stratonix model with an implicit interpolator from
the exported LoopStructural input set, and compare depths with Stratonix's own grids.

Usage: crosscheck_implicit.py <export folder> <report folder>

The export folder is one written by `stratonix export ... --formats petrel,csv,loopstructural`.

Interpolator:
  - LoopStructural, when it can be imported (GeologicalModel with one fault per exported fault and a
    'strati' foliation). Not installable in the environment this script was written in, so that
    path has not been run.
  - Otherwise a substitute implicit interpolator from SciPy: radial basis functions (thin-plate
    spline) fitted to each horizon's interface points, in two variants:
      A "fault-blind": all points of the horizon in one surface (no fault knowledge);
      B "per-block":   one surface per Stratonix fault block (block of each point from the
                       stx_block column, block of each grid node from the CSV export).
The report labels which interpolator was used. Differences are reported per horizon, for nodes
near faults (within `near` metres of a fault gap) and away from them.
"""
import sys, pathlib, json, csv, io, datetime

import numpy as np

NEAR = 500.0  # m: "near a fault" = within this distance of a fault-gap node


def read_irap(p):
    t = p.read_text().split()
    ny, nx = int(t[1]), int(t[8])
    dx, dy, x0, y0 = float(t[2]), float(t[3]), float(t[4]), float(t[6])
    z = np.array(t[19:19 + nx * ny], dtype=float).reshape(ny, nx)
    z[z >= 9999899] = np.nan
    X, Y = np.meshgrid(x0 + dx * np.arange(nx), y0 + dy * np.arange(ny))
    return X, Y, z, dx


def read_block_grid(p, X, Y):
    blk = np.full(X.shape, -1, dtype=np.int64)
    x0, y0, dx, dy = X[0, 0], Y[0, 0], X[0, 1] - X[0, 0], Y[1, 0] - Y[0, 0]
    for r in csv.DictReader(io.StringIO("\n".join(l for l in p.read_text().splitlines() if not l.startswith("#")))):
        i, j = int(round((float(r["x"]) - x0) / dx)), int(round((float(r["y"]) - y0) / dy))
        blk[j, i] = int(r["block"])
    return blk


def distance_to_gaps(z, dx):
    """Distance (m) from every node to the nearest undefined node (inf if none)."""
    from scipy.ndimage import distance_transform_edt
    gaps = np.isnan(z)
    if not gaps.any():
        return np.full(z.shape, np.inf)
    return distance_transform_edt(~gaps) * dx


def stats(d):
    d = np.abs(d[np.isfinite(d)])
    if len(d) == 0:
        return dict(n=0, mean=np.nan, p95=np.nan, max=np.nan)
    return dict(n=int(len(d)), mean=float(d.mean()), p95=float(np.percentile(d, 95)), max=float(d.max()))


def rbf_surface(px, py, pz, X, Y):
    from scipy.interpolate import RBFInterpolator
    f = RBFInterpolator(np.c_[px, py], pz, kernel="thin_plate_spline", smoothing=0.0)
    return f(np.c_[X.ravel(), Y.ravel()]).reshape(X.shape)


def loopstructural_surfaces(root, horizons, grids):  # pragma: no cover - LoopStructural not installed here
    import pandas as pd
    from LoopStructural import GeologicalModel
    setup = json.loads((root / "loopstructural/model_setup.json").read_text())
    df = pd.read_csv(root / "loopstructural/data.csv")
    df = df.rename(columns={"gx": "nx", "gy": "ny", "gz": "nz"})
    model = GeologicalModel(setup["origin"], setup["maximum"])
    model.data = df.drop(columns=[c for c in df.columns if c.startswith("stx_")])
    for f in setup["faults_oldest_first"]:
        model.create_and_add_fault(f["name"], f["displacement_max_m"])
    model.create_and_add_foliation("strati")
    model.update()
    out = {}
    for h in horizons:
        X, Y, z, dx = grids[h]
        zs = np.linspace(setup["origin"][2], setup["maximum"][2], 400)
        target = setup["strati_values"][h]
        best = np.full(X.shape, np.nan)
        prev = None
        for zk in zs:
            v = model.evaluate_feature_value("strati", np.c_[X.ravel(), Y.ravel(), np.full(X.size, zk)]).reshape(X.shape)
            if prev is not None:
                cross = np.isnan(best) & ((prev[0] - target) * (v - target) <= 0)
                t = (target - prev[0]) / np.where(v != prev[0], v - prev[0], 1)
                best[cross] = (prev[1] + t * (zk - prev[1]))[cross]
            prev = (v, zk)
        out[h] = best
    return out


def main(folder, report):
    root, rep = pathlib.Path(folder), pathlib.Path(report)
    rep.mkdir(parents=True, exist_ok=False)
    rows = list(csv.DictReader(io.StringIO((root / "loopstructural/data.csv").read_text())))
    horizons = sorted({r["stx_horizon"] for r in rows if r["feature_name"] == "strati" and r["val"]},
                      key=lambda h: -float(next(r["val"] for r in rows if r["stx_horizon"] == h and r["val"])))
    grids = {h: read_irap(root / "horizons" / f"{h}_faulted.irap") for h in horizons}
    try:
        import LoopStructural  # noqa: F401
        engine = "LoopStructural " + LoopStructural.__version__
        surfaces = {"LoopStructural": loopstructural_surfaces(root, horizons, grids)}
    except ImportError:
        engine = "SciPy RBFInterpolator (thin-plate spline) - substitute: LoopStructural not installable here"
        surfaces = {"A fault-blind": {}, "B per-block": {}}
        for h in horizons:
            X, Y, z, dx = grids[h]
            pts = np.array([[float(r["X"]), float(r["Y"]), float(r["Z"]), int(r["stx_block"])]
                            for r in rows if r["feature_name"] == "strati" and r["val"] and r["stx_horizon"] == h])
            surfaces["A fault-blind"][h] = rbf_surface(pts[:, 0], pts[:, 1], pts[:, 2], X, Y)
            blk = read_block_grid(root / "csv" / f"{h}.csv", X, Y)
            sb = np.full(X.shape, np.nan)
            for b in np.unique(pts[:, 3]).astype(int):
                sel = pts[:, 3] == b
                if sel.sum() < 4:
                    continue
                s = rbf_surface(pts[sel, 0], pts[sel, 1], pts[sel, 2], X, Y)
                sb[blk == b] = s[blk == b]
            surfaces["B per-block"][h] = sb
    lines = []
    table = []
    for variant, surf in surfaces.items():
        for h in horizons:
            X, Y, z, dx = grids[h]
            d = surf[h] - z
            dist = distance_to_gaps(z, dx)
            near, away = np.isfinite(d) & (dist <= NEAR), np.isfinite(d) & (dist > NEAR)
            sn, sa = stats(d[near]), stats(d[away])
            table.append([variant, h, sa["n"], sa["mean"], sa["p95"], sa["max"], sn["n"], sn["mean"], sn["p95"], sn["max"]])
    with open(rep / "crosscheck_table.csv", "w", newline="\n") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["variant", "horizon", "away_n", "away_mean_abs_m", "away_p95_m", "away_max_m", "near_n", "near_mean_abs_m", "near_p95_m", "near_max_m"])
        for r in table:
            w.writerow([r[0], r[1], r[2]] + [f"{v:.3f}" if isinstance(v, float) and np.isfinite(v) else ("-" if isinstance(v, float) else v) for v in r[3:6]] +
                       [r[6]] + [f"{v:.3f}" if np.isfinite(v) else "-" for v in r[7:10]])
    lines.append("# External cross-check of horizon depths\n")
    lines.append(f"Date: {datetime.date.today().isoformat()}  ")
    lines.append(f"Export folder: {root.name}  ")
    lines.append(f"Interpolator: {engine}  ")
    lines.append(f"Near a fault: within {NEAR:.0f} m of a fault-gap node of the Stratonix grid.\n")
    lines.append("| Variant | Horizon | Away: nodes | mean abs (m) | p95 (m) | max (m) | Near: nodes | mean abs (m) | p95 (m) | max (m) |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in table:
        f = lambda v: "-" if not np.isfinite(v) else f"{v:.2f}"
        lines.append(f"| {r[0]} | {r[1]} | {r[2]} | {f(r[3])} | {f(r[4])} | {f(r[5])} | {r[6]} | {f(r[7])} | {f(r[8])} | {f(r[9])} |")
    (rep / "crosscheck_report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
