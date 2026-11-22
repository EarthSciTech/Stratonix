#!/usr/bin/env python3
"""Stage 4A external check (STX-TN-14 section 10): LoopStructural's single-fold example.

The Noddy single-fold data set shipped with LoopStructural is shuffled with its own `random`
column; the first N rows (as in LoopStructural's "adding folds" example) and the axial-foliation
rows (s1) are the training data. Both codes then predict the bedding (s0) normal at every s0
orientation not used for training:

  LoopStructural  fold frame on s1 and a folded foliation on s0 (create_and_add_fold_frame,
                  create_and_add_folded_foliation), gradient of s0 at the held-out points
  Stratonix       `stratonix fold-frame` on the s0 orientations of the training rows, with the s1
                  orientations as axial foliation

The angle between each prediction and the measured normal (sign ignored) is reported as p50 and
p95. Runs outside the product; needs LoopStructural. Usage: ls_fold_v4.py <stratonix> <report>
"""
import json, pathlib, subprocess, sys, tempfile

import numpy as np


def normals(df):
    cols = set(df.columns)
    for c in (("nx", "ny", "nz"), ("gx", "gy", "gz")):
        if set(c) <= cols:
            return df[list(c)].to_numpy(dtype=float)
    if {"strike", "dip"} <= cols:
        s, d = np.radians(df["strike"].to_numpy(float)), np.radians(df["dip"].to_numpy(float))
        dd = s + np.pi / 2
        return np.stack([np.sin(dd) * np.sin(d), np.cos(dd) * np.sin(d), np.cos(d)], axis=1)
    raise SystemExit("no orientation columns in " + ", ".join(df.columns))


def angles(a, b):
    a = a / np.linalg.norm(a, axis=1, keepdims=True)
    b = b / np.linalg.norm(b, axis=1, keepdims=True)
    return np.degrees(np.arccos(np.clip(np.abs(np.sum(a * b, axis=1)), 0, 1)))


def main(exe, report, npoints=100):
    from LoopStructural import GeologicalModel, __version__ as lsv
    from LoopStructural.datasets import load_noddy_single_fold
    report = pathlib.Path(report)
    report.mkdir(parents=True, exist_ok=True)
    data, bb = load_noddy_single_fold()
    print("LoopStructural", lsv, "columns:", list(data.columns))
    if "random" in data.columns:
        data = data.sort_values("random")
    s0 = data[data["feature_name"] == "s0"]
    s1 = data[data["feature_name"] == "s1"]
    train = data[:npoints]
    n_all = normals(s0)
    has_n = np.all(np.isfinite(n_all), axis=1)
    train_idx = set(train.index)
    hold = s0[has_n & ~s0.index.isin(train_idx)]
    tr0 = s0[has_n & s0.index.isin(train_idx)]
    n_s1 = normals(s1)
    s1n = s1[np.all(np.isfinite(n_s1), axis=1)]
    print(f"training: {len(tr0)} s0 orientations, {len(s1n)} s1 orientations; held out: {len(hold)}")
    truth = normals(hold)
    # LoopStructural
    import pandas as pd
    model = GeologicalModel(bb[0, :], bb[1, :])
    model.set_model_data(pd.concat([train, s1]))
    frame = model.create_and_add_fold_frame("s1", interpolatortype="PLI", nelements=10000, buffer=0.5, solver="pyamg", damp=True)
    model.create_and_add_folded_foliation("s0", frame, nelements=10000, buffer=0.5)
    xyz = hold[["X", "Y", "Z"]].to_numpy(float)
    g_ls = model.evaluate_feature_gradient("s0", xyz)
    e_ls = angles(g_ls, truth)
    # Stratonix
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)

        def write(path, df, nn):
            with open(path, "w") as f:
                f.write("X,Y,Z,G_x,G_y,G_z,surface\n")
                for (x, y, z), n in zip(df[["X", "Y", "Z"]].to_numpy(float), nn):
                    f.write(f"{x},{y},{z},{n[0]},{n[1]},{n[2]},s\n")
        write(td / "bed.csv", tr0, normals(tr0))
        write(td / "fol.csv", s1n, normals(s1n))
        with open(td / "at.csv", "w") as f:
            f.write("X,Y,Z,surface\n")
            for x, y, z in xyz:
                f.write(f"{x},{y},{z},s\n")
        args = [exe, "fold-frame", str(td / "bed.csv"), "--at", str(td / "at.csv"), "-o", str(td / "frame.json")]
        if len(s1n) >= 3:
            args[3:3] = ["--foliation", str(td / "fol.csv")]
        print(subprocess.run(args, capture_output=True, text=True).stdout)
        fr = json.loads((td / "frame.json").read_text())
    if not fr.get("ok"):
        print("Stratonix fold frame not fitted:", fr.get("reason"))
        e_st = np.full(len(truth), np.nan)
    else:
        e_st = angles(np.array(fr["predicted_normals"], float), truth)
    res = {"loopstructural": lsv, "held_out": int(len(truth)), "training_s0": int(len(tr0)), "training_s1": int(len(s1n)),
           "ls_p50_deg": float(np.nanpercentile(e_ls, 50)), "ls_p95_deg": float(np.nanpercentile(e_ls, 95)),
           "stx_p50_deg": float(np.nanpercentile(e_st, 50)) if np.isfinite(e_st).any() else None,
           "stx_p95_deg": float(np.nanpercentile(e_st, 95)) if np.isfinite(e_st).any() else None, "frame": {k: fr.get(k) for k in
           ("trend_deg", "plunge_deg", "interlimb_deg", "wavelength_m", "order", "cylindricity", "rms_deg", "reason")}}
    ok = res["stx_p95_deg"] is not None and res["stx_p95_deg"] <= max(15.0, res["ls_p95_deg"] + 5.0)
    res["verdict"] = "within tolerance" if ok else "OUTSIDE tolerance"
    (report / "ls_fold_v4.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
