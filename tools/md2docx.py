#!/usr/bin/env python3
"""Build a technical note (.docx) from a Markdown source.

Supports the subset used by the Stratonix notes: # headings (1-3), paragraphs,
'-' and '1.' lists, pipe tables, fenced code blocks, **bold**, `code`, and a
front-matter block (title, subtitle, doc_id, version, date, status).

Usage: python3 md2docx.py <source.md> <output.docx>
"""
import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt, RGBColor, Cm

INLINE = re.compile(r"(\*\*[^*]+\*\*|`[^`]+`|\*[^*\s][^*]*\*)")


def add_runs(par, text):
    for part in INLINE.split(text):
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            par.add_run(part[2:-2]).bold = True
        elif len(part) > 2 and part.startswith("*") and part.endswith("*"):
            par.add_run(part[1:-1]).italic = True
        elif part.startswith("`") and part.endswith("`"):
            r = par.add_run(part[1:-1])
            r.font.name = "Consolas"
            r.font.size = Pt(9)
        else:
            par.add_run(part)


def shade(cell, hex_color):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hex_color)
    tcPr.append(shd)


def parse_front(lines):
    meta = {}
    if lines and lines[0].strip() == "---":
        end = lines.index("---", 1) if "---" in lines[1:] else -1
        for l in lines[1:end]:
            if ":" in l:
                k, v = l.split(":", 1)
                meta[k.strip()] = v.strip()
        lines = lines[end + 1:]
    return meta, lines


def continued(lines, i, text):
    """Join indented continuation lines of a list item (not a new item, not blank)."""
    while i < len(lines) and lines[i].startswith("  ") and lines[i].strip() and not re.match(r"^\s*([-*]|\d+\.) ", lines[i]):
        text += " " + lines[i].strip()
        i += 1
    return text, i


def build(src: Path, out: Path):
    lines = src.read_text(encoding="utf-8").splitlines()
    lines = [l.rstrip("\n") for l in lines]
    meta, lines = parse_front([l.rstrip() for l in lines])
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Calibri"
    st.font.size = Pt(10.5)
    for s in doc.sections:
        s.left_margin = s.right_margin = Cm(2.2)
    if meta:
        p = doc.add_paragraph()
        r = p.add_run(meta.get("title", ""))
        r.bold = True
        r.font.size = Pt(20)
        r.font.color.rgb = RGBColor(0x1F, 0x3A, 0x5F)
        if meta.get("subtitle"):
            p = doc.add_paragraph()
            r = p.add_run(meta["subtitle"])
            r.italic = True
            r.font.size = Pt(12)
        t = doc.add_table(rows=0, cols=2)
        t.style = "Light Grid Accent 1"
        for k in ("doc_id", "version", "date", "status", "author"):
            if k in meta:
                row = t.add_row().cells
                row[0].text = k.replace("_", " ").capitalize()
                row[1].text = meta[k]
        doc.add_paragraph()
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.startswith("```"):
            i += 1
            buf = []
            while i < n and not lines[i].startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Cm(0.4)
            r = p.add_run("\n".join(buf))
            r.font.name = "Consolas"
            r.font.size = Pt(8)
            pPr = p._p.get_or_add_pPr()
            shd = OxmlElement("w:shd")
            shd.set(qn("w:val"), "clear")
            shd.set(qn("w:fill"), "F2F4F7")
            pPr.append(shd)
            continue
        m = re.match(r"^(#{1,3}) (.*)", line)
        if m:
            doc.add_heading(m.group(2).strip(), level=len(m.group(1)))
            i += 1
            continue
        if line.startswith("|"):
            rows = []
            while i < n and lines[i].startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells):
                    rows.append(cells)
                i += 1
            ncol = max(len(r) for r in rows)
            t = doc.add_table(rows=len(rows), cols=ncol)
            t.style = "Table Grid"
            t.alignment = WD_TABLE_ALIGNMENT.CENTER
            for ri, r in enumerate(rows):
                for ci in range(ncol):
                    cell = t.rows[ri].cells[ci]
                    cell.text = ""
                    par = cell.paragraphs[0]
                    add_runs(par, r[ci] if ci < len(r) else "")
                    for run in par.runs:
                        run.font.size = Pt(8.5)
                        if ri == 0:
                            run.bold = True
                    if ri == 0:
                        shade(cell, "D9E2F3")
            doc.add_paragraph()
            continue
        m = re.match(r"^(\s*)[-*] (.*)", line)
        if m:
            p = doc.add_paragraph(style="List Bullet" if len(m.group(1)) < 2 else "List Bullet 2")
            text, i = continued(lines, i + 1, m.group(2))
            add_runs(p, text)
            continue
        m = re.match(r"^(\s*)\d+\. (.*)", line)
        if m:
            p = doc.add_paragraph(style="List Number")
            text, i = continued(lines, i + 1, m.group(2))
            add_runs(p, text)
            continue
        if not line.strip():
            i += 1
            continue
        buf = [line.strip()]
        i += 1
        while i < n and lines[i].strip() and not re.match(r"^(#|\||```|\s*[-*] |\s*\d+\. )", lines[i]):
            buf.append(lines[i].strip())
            i += 1
        p = doc.add_paragraph()
        add_runs(p, " ".join(buf))
    # footer with doc id
    if meta.get("doc_id"):
        for s in doc.sections:
            fp = s.footer.paragraphs[0]
            fp.text = f"{meta['doc_id']} v{meta.get('version', '')} | {meta.get('status', '')}"
            fp.alignment = WD_ALIGN_PARAGRAPH.RIGHT
            for r in fp.runs:
                r.font.size = Pt(8)
    out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    build(Path(sys.argv[1]), Path(sys.argv[2]))
