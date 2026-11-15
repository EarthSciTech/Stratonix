#!/usr/bin/env python3
"""Regenerate the stage 3A implicit examples (STX-TN-11 section 10) from their parametric sources.

Usage: make_implicit_examples.py <stratonix executable>

For each example, the parametric source description is exported with
`stratonix export <source> --formats gempy,fault_sticks_irap`. The surface points, orientations
and fault sticks are copied to examples/data/<data folder>/. The implicit description is written
from the source by:
  - removing the design structure (structure, isochores, series tilts, clinoforms, events);
  - turning each fault with a band into sticks, with its throw estimated from the data; other
    faults are kept as stated in the source;
  - keeping horizons defined by isochores as isochore horizons, their thickness from the data;
  - adding the inputs and an implicit block.
Run from anywhere; paths are relative to this script's repository.
"""
import pathlib, re, shutil, subprocess, sys, tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
EX = ROOT / "examples"

# name: (source description, data folder, {fault: throw-estimate band in m} (others keep their stated throws), title)
# A fault within a band's width of another (the F2 splay, 250 m from F1) keeps its stated throw:
# its band would take in the other fault's offset.
EXAMPLES = {
    "implicit_unconformity_from_samples": ("unconformity_onlap_pinchout.yaml", "implicit_unconformity", None,
                                           "Angular unconformity and onlap rebuilt from sampled data (implicit builder, one field per series)"),
    "implicit_relay_splay_from_samples": ("fault_system_relay_splay.yaml", "implicit_relay_splay", {"F1": 600, "F3": 600},
                                          "Relay and splay fault system rebuilt from sampled data, throws estimated from the data"),
    "implicit_detachment_folds_from_samples": ("detachment_folds.yaml", "implicit_detachment_folds", None,
                                               "Detachment folds rebuilt from sampled data (implicit builder)"),
    "implicit_faults_series_from_data": ("data/implicit_faults_series/truth.yaml", "implicit_faults_series", {"F1": 400},
                                         "Unconformity, onlap and a fault rebuilt from data: one field per series, throw from the data"),
}


def blocks(text):
    """Top-level YAML blocks: {key: text}, in order (comments before the first key kept apart)."""
    out, cur, head = {}, None, []
    for line in text.splitlines(keepends=True):
        m = re.match(r"^([a-z_]+):", line)
        if m:
            cur = m.group(1)
            out[cur] = line
        elif cur is None:
            head.append(line)
        else:
            out[cur] += line
    return "".join(head), out


def fault_entries(faults_block):
    lines = faults_block.splitlines(keepends=True)[1:]
    entries, cur = [], None
    for l in lines:
        if l.startswith("  - name:"):
            cur = [l]
            entries.append(cur)
        elif cur is not None:
            cur.append(l)
    return ["".join(e) for e in entries]


def implicit_fault(entry, bands):
    name = re.search(r"- name: (\S+)", entry).group(1)
    band = (bands or {}).get(name)
    age = re.search(r"\n    age: (.+)\n", entry).group(1)
    pdip = re.search(r"profile_dip: (\{.*\})", entry)
    pstr = re.search(r"profile_strike: (\{.*\})", entry)
    keep = []
    for key in ("relay", "abuts", "drag"):
        m = re.search(r"\n    (" + key + r": .+)\n", entry)
        if m:
            keep.append("    " + m.group(1) + "\n")
    if not band:
        return entry  # kept as interpreted: geometry and throw as stated in the source
    planar = "type: planar" in entry and "slip_mode: dip_slip" in entry
    s = f"  - name: {name}\n    age: {age}\n"
    s += f"    geometry: {{type: sticks, input: STICKS_{name}, input_format: irap_lines, fit: plane}}\n"
    if band and planar:
        s += "    throw:\n"
        s += f"      estimate: {{band: \"{band} m | assumed | data within {band} m of the fault\", station_spacing: \"250 m | assumed | the data spacing\"}}\n"
    else:
        ext = re.search(r"\n    extent:\n((?:      .+\n)+)", entry).group(0)
        s += ext.lstrip("\n")
        thr = re.search(r"\n    throw:\n((?:      .+\n)+)", entry).group(0)
        s += thr.lstrip("\n")
        return s + re.search(r"\n    slip_mode: .+\n    hanging_wall: .+\n", entry).group(0).lstrip("\n") + "".join(keep)
    if pdip:
        s += f"      profile_dip: {pdip.group(1)}\n"
    if pstr:
        s += f"      profile_strike: {pstr.group(1)}\n"
    s += "    slip_mode: dip_slip\n    hanging_wall: {type: rigid}\n" + "".join(keep)
    return s


def main(exe):
    for name, (src, folder, band, title) in EXAMPLES.items():
        srcp = EX / src
        text = srcp.read_text()
        with tempfile.TemporaryDirectory() as td:
            subprocess.run([exe, "export", str(srcp), "-o", td + "/out", "--formats", "gempy,fault_sticks_irap"], check=True,
                           stdout=subprocess.DEVNULL)
            dd = EX / "data" / folder
            dd.mkdir(parents=True, exist_ok=True)
            shutil.copy(td + "/out/exports/gempy/surface_points.csv", dd / "surface_points.csv")
            shutil.copy(td + "/out/exports/gempy/orientations.csv", dd / "orientations.csv")
            sticks = sorted(pathlib.Path(td + "/out/exports/faults").glob("*_sticks.dat")) if pathlib.Path(td + "/out/exports/faults").exists() else []
            for st in sticks:
                shutil.copy(st, dd / st.name)
        head, b = blocks(text)
        model_name = re.search(r"name: (\S+?),", b["model"]).group(1)
        src_version = re.search(r'version: "([^"]+)"', b["model"]).group(1)
        out = (f"# Stage 3A example (STX-TN-11): {title}.\n"
               f"# The data were sampled from the parametric example {src} (model {model_name} v{src_version}) with\n"
               f"#   stratonix export {src} -o out --formats gempy,fault_sticks_irap\n"
               f"# by tools/make_implicit_examples.py, and are rebuilt here from data alone.\n")
        out += "stratonix_schema: 1\n"
        out += f'model: {{name: {name}, version: "1.0.0", status: DESIGN,\n        title: "{title}"}}\n'
        out += b["crs"] + b["domain"]
        src_note = f"sampled from {src} ({model_name} v{src_version}) every 250 m"
        out += "inputs:\n"
        out += f'  - {{id: POINTS, path: data/{folder}/surface_points.csv, role: surface points, sharing: shareable, source: "{src_note}", evidence: user_instruction}}\n'
        out += f'  - {{id: ORIENT, path: data/{folder}/orientations.csv, role: orientations, sharing: shareable, source: "{src_note}", evidence: user_instruction}}\n'
        faults = fault_entries(b["faults"]) if "faults" in b else []
        for fe in faults:
            fn = re.search(r"- name: (\S+)", fe).group(1)
            if not (band or {}).get(fn):
                continue
            out += f'  - {{id: STICKS_{fn}, path: data/{folder}/{fn}_sticks.dat, role: fault sticks, sharing: shareable, source: "sticks of {fn} from {src}", evidence: user_instruction}}\n'
        out += b["wells"]
        strat = b["stratigraphy"]
        iso = re.findall(r"\{name: (\w+)[^\n]*isochore:", strat)                                 # isochore horizons stay isochores,
        strat = re.sub(r", isochore: \{.*\}\}(\})", r"\1", strat)                               # their thickness from the data
        strat = re.sub(r", tilt: \{[^}]*\}\}", "}", strat)                                     # series tilts: design only
        strat = re.sub(r"\n  compaction:\n(    .+\n)+", "\n", strat)
        out += strat
        out += ("implicit:\n  data:\n    - {input: POINTS, kind: points, format: csv_points}\n"
                "    - {input: ORIENT, kind: orientations, format: csv_orientations}\n  use_well_tops: true\n  kernel: cubic\n"
                '  nugget: "0 m | assumed | exact interpolation of the data"\n'
                '  surface_spacing: "100 m | assumed | twice the grid spacing"\n'
                '  iterations: "10 ordinal | assumed | maximum restoration passes"\n'
                '  vertical_stretch: "5 ratio | assumed | anisotropy: horizons vary slowly laterally compared with their spacing"\n')
        if iso:
            out += "  isochore_horizons: [" + ", ".join(iso) + "]\n"
        if faults:
            out += "faults:\n" + "".join(implicit_fault(fe, band) for fe in faults)
        for key in ("frame", "checks"):
            if key in b:
                out += b[key]
        if "outputs" in b:
            o = b["outputs"]
            o = re.sub(r"\n  isopachs: .+", "", o)
            o = re.sub(r"\n  curvature:\n(    .+\n)+", "\n", o)
            out += o
        (EX / f"{name}.yaml").write_text(out)
        print("wrote", name)


if __name__ == "__main__":
    main(sys.argv[1])
