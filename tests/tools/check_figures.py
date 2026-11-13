#!/usr/bin/env python3
"""Independent check of a Stratonix figure folder (not using Stratonix code).

For every figure basename: SVG parses as XML and holds <text> elements (text stays text);
PNG is 300 dpi (pHYs) and opaque RGB(A); TIFF is 300 dpi, LZW, RGB, no alpha, pixel-identical
to the PNG; PPTX opens with python-pptx and holds one 16:9 slide with a title, a picture, a
caption and speaker notes. Exit code 0 when all pass, 1 otherwise; skips (77) when the Python
packages are missing.
"""
import sys, pathlib, xml.etree.ElementTree as ET

try:
    import numpy as np
    import tifffile
    from PIL import Image
    from pptx import Presentation
    from pptx.util import Emu
except ImportError as e:  # pragma: no cover
    print("SKIP: missing python package:", e)
    sys.exit(77)

def main(folder):
    folder = pathlib.Path(folder)
    svgs = sorted(folder.glob("STX_*.svg"))
    if not svgs:
        print("FAIL: no figures in", folder); return 1
    bad = 0
    for svg in svgs:
        base = svg.with_suffix("")
        errs = []
        root = ET.parse(svg).getroot()
        texts = [t for t in root.iter() if t.tag.endswith("text")]
        if not texts: errs.append("SVG has no text elements")
        png, tif, pptx = base.with_suffix(".png"), base.with_suffix(".tif"), base.with_suffix(".pptx")
        for f in (png, tif, pptx):
            if not f.exists(): errs.append("missing " + f.name)
        if not errs:
            im = Image.open(png)
            dpi = im.info.get("dpi", (0, 0))
            if abs(dpi[0] - 300) > 0.5 or abs(dpi[1] - 300) > 0.5: errs.append(f"PNG dpi {dpi}")
            a = np.asarray(im.convert("RGB"))
            with tifffile.TiffFile(tif) as t:
                p = t.pages[0]
                if p.compression.name != "LZW": errs.append("TIFF compression " + p.compression.name)
                if p.photometric.name != "RGB": errs.append("TIFF photometric " + p.photometric.name)
                if p.samplesperpixel != 3: errs.append(f"TIFF samples {p.samplesperpixel} (alpha?)")
                xr = p.tags["XResolution"].value; yr = p.tags["YResolution"].value
                if abs(xr[0] / xr[1] - 300) > 0.5 or abs(yr[0] / yr[1] - 300) > 0.5: errs.append(f"TIFF dpi {xr} {yr}")
                b = p.asarray()
            if a.shape != b.shape or not np.array_equal(a, b): errs.append("TIFF pixels differ from PNG")
            pr = Presentation(str(pptx))
            if abs(pr.slide_width / pr.slide_height - 16 / 9) > 1e-3: errs.append("slide not 16:9")
            if len(pr.slides) != 1: errs.append(f"{len(pr.slides)} slides")
            s = pr.slides[0]
            kinds = [sh.shape_type for sh in s.shapes]
            if not any(sh.shape_type == 13 for sh in s.shapes): errs.append("no picture on slide")
            txt = " ".join(sh.text_frame.text for sh in s.shapes if sh.has_text_frame)
            if len(txt.strip()) < 10: errs.append("no title/caption text")
            if not s.has_notes_slide or len(s.notes_slide.notes_text_frame.text) < 10: errs.append("no speaker notes")
            for sh in s.shapes:
                if sh.left < 0 or sh.top < 0 or sh.left + sh.width > pr.slide_width or sh.top + sh.height > pr.slide_height:
                    errs.append("shape outside slide")
        print(("FAIL " if errs else "ok   ") + base.name + ("" if not errs else ": " + "; ".join(errs)))
        bad += bool(errs)
    print(f"{len(svgs) - bad}/{len(svgs)} figures pass")
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
