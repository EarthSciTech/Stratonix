#!/usr/bin/env python3
"""Build the Stratonix documents (.docx) from the Markdown sources in docs/.

Usage: build_docs.py <repository root> <output folder> <version tag, e.g. v1_0>
Writes: user guide, developer guide, roadmap and API reference (regenerated from the headers).
"""
import subprocess, sys, pathlib

DOCS = [("user_guide.md", "STX-UG_user_guide"), ("developer_guide.md", "STX-DG_developer_guide"),
        ("roadmap.md", "STX-RM_roadmap_v2-v5"), ("api_reference.md", "STX-API_reference")]

def main(root, out, tag):
    root, out = pathlib.Path(root), pathlib.Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tools = root / "tools"
    subprocess.run([sys.executable, str(tools / "gen_api_reference.py"), str(root), str(root / "docs/api_reference.md")], check=True)
    for src, name in DOCS:
        dst = out / f"{name}_{tag}.docx"
        if dst.exists():
            raise SystemExit(f"refusing to overwrite an issued document: {dst}")
        subprocess.run([sys.executable, str(tools / "md2docx.py"), str(root / "docs" / src), str(dst)], check=True)

if __name__ == "__main__":
    main(*sys.argv[1:4])
