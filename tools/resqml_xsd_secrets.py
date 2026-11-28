#!/usr/bin/env python3
"""RESQML 2.0.1 schemas for CI, carried in repository secrets (I-034).

The Energistics schemas are licensed: they are never committed to the repository and never hosted
anywhere. The owner packs them once into a few repository secrets; the `resqml` CI job rebuilds the
schema folder from those secrets and validates every RESQML export against it.

A GitHub secret holds at most 48 KB, and the schema package is larger even when compressed, so:

  * the documentation (xs:annotation) is removed from every .xsd file; validation ignores it;
  * the .xsd files, with their folder layout (the schemas import each other by relative path), go
    into one tar.xz archive, written byte for byte the same for the same input;
  * the archive is base64-encoded and split into parts of at most 48 000 characters;
  * part 1 starts with a header line "STXXSD1 <parts> <sha256 of the archive>", so CI uses exactly
    the parts of the latest packing and checks the archive before unpacking it.

Owner, once (and again only if the schemas change):

    python3 tools/resqml_xsd_secrets.py pack <downloaded package .zip or folder> -o xsd_parts

(on Windows: py tools\resqml_xsd_secrets.py pack ...), then set each part as a secret,
STX_RESQML_XSD_1, STX_RESQML_XSD_2, ... (up to 8), either with the commands the script prints
(GitHub CLI, PowerShell or bash) or by pasting each file's text in Settings > Secrets and
variables > Actions. Delete xsd_parts afterwards.

CI:

    python3 tools/resqml_xsd_secrets.py unpack <folder>     # reads STX_RESQML_XSD_1..8

prints nothing of the schemas themselves, only the number of files and the archive checksum.
"""
import argparse
import base64
import hashlib
import io
import lzma
import os
import pathlib
import re
import sys
import tarfile
import zipfile
import xml.etree.ElementTree as ET

MAGIC = "STXXSD1"
PART_CHARS = 48000          # GitHub's limit is 48 KB per secret
MAX_PARTS = 8               # the CI job reads STX_RESQML_XSD_1 .. STX_RESQML_XSD_8
MAIN_SCHEMAS = ("ResqmlAllObjects.xsd", "EmlAllObjects.xsd")
XSD_NS = "http://www.w3.org/2001/XMLSchema"


def read_package(src):
    """All .xsd files of a package (zip, nested zips or folder) as {relative path: bytes}."""
    files = {}

    def from_zip(data, prefix=""):
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            for name in z.namelist():
                if name.endswith("/"):
                    continue
                if name.lower().endswith(".xsd"):
                    files[prefix + name] = z.read(name)
                elif name.lower().endswith(".zip"):
                    from_zip(z.read(name), prefix + name[:-4] + "/")

    p = pathlib.Path(src)
    if p.is_dir():
        for f in sorted(p.rglob("*")):
            if f.is_file() and f.suffix.lower() == ".xsd":
                files[f.relative_to(p).as_posix()] = f.read_bytes()
            elif f.is_file() and f.suffix.lower() == ".zip":
                from_zip(f.read_bytes(), f.relative_to(p).as_posix()[:-4] + "/")
    elif p.is_file() and zipfile.is_zipfile(p):
        from_zip(p.read_bytes())
    else:
        sys.exit("%s is neither a folder nor a zip file" % src)
    return files


def strip_annotations(text):
    """Remove xs:annotation elements by text, keeping every namespace prefix as written: the schemas
    refer to types by prefixed name (type="resqml2:...") and an XML rewrite would rename prefixes."""
    prefixes = re.findall(r'xmlns:([A-Za-z_][\w.-]*)\s*=\s*["\']' + re.escape(XSD_NS) + r'["\']', text)
    for p in prefixes:
        text = re.sub(r"<%s:annotation\b[^>]*/>" % re.escape(p), "", text)
        text = re.sub(r"<%s:annotation\b.*?</%s:annotation\s*>" % (re.escape(p), re.escape(p)), "", text, flags=re.S)
    if re.search(r'xmlns\s*=\s*["\']' + re.escape(XSD_NS) + r'["\']', text):  # schema namespace as the default
        text = re.sub(r"<annotation\b[^>]*/>", "", text)
        text = re.sub(r"<annotation\b.*?</annotation\s*>", "", text, flags=re.S)
    return re.sub(r"\n\s*\n+", "\n", text)


def common_root(paths):
    """Drop the leading folders every file shares, so the archive starts at the package's top folder."""
    parts = [p.split("/")[:-1] for p in paths]
    root = []
    for level in zip(*parts):
        if len(set(level)) != 1:
            break
        root.append(level[0])
    return "/".join(root) + "/" if root else ""


def make_archive(files):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.USTAR_FORMAT) as tar:
        for name in sorted(files):
            info = tarfile.TarInfo(name)
            info.size = len(files[name])
            info.mtime = 0
            info.mode = 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tar.addfile(info, io.BytesIO(files[name]))
    return lzma.compress(buf.getvalue(), preset=9 | lzma.PRESET_EXTREME)


def check_schemas(folder):
    """Compile the main schemas when lxml is available (the CI job's validator)."""
    try:
        from lxml import etree
    except ImportError:
        return "lxml not installed: schemas not compiled here"
    out = []
    for main in MAIN_SCHEMAS:
        path = next(pathlib.Path(folder).rglob(main), None)
        if path is None:
            out.append("%s missing" % main)
            continue
        etree.XMLSchema(etree.parse(str(path)))  # raises on a broken schema set
        out.append("%s compiles" % main)
    return ", ".join(out)


def pack(args):
    files = read_package(args.package)
    if not files:
        sys.exit("no .xsd files in %s" % args.package)
    names = {pathlib.PurePosixPath(n).name for n in files}
    missing = [m for m in MAIN_SCHEMAS if m not in names]
    if missing:
        sys.exit("the package has no %s: is it the RESQML 2.0.1 schema package?" % ", ".join(missing))
    root = common_root(list(files))
    stripped, raw = {}, 0
    for name, data in files.items():
        raw += len(data)
        text = strip_annotations(data.decode("utf-8-sig"))
        ET.fromstring(text.encode("utf-8"))  # still well-formed
        stripped[name[len(root):]] = text.encode("utf-8")
    archive = make_archive(stripped)
    sha = hashlib.sha256(archive).hexdigest()
    b64 = base64.b64encode(archive).decode("ascii")
    body_first = PART_CHARS - 100  # room for the header line
    chunks = [b64[:body_first]] + [b64[i:i + PART_CHARS] for i in range(body_first, len(b64), PART_CHARS)]
    chunks = [c for c in chunks if c] or [""]
    if len(chunks) > MAX_PARTS:
        sys.exit("the packed schemas need %d parts, more than the %d the CI job reads" % (len(chunks), MAX_PARTS))
    chunks[0] = "%s %d %s\n%s" % (MAGIC, len(chunks), sha, chunks[0])
    out = pathlib.Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    for old in out.glob("STX_RESQML_XSD_*.txt"):
        old.unlink()
    for i, c in enumerate(chunks, 1):
        (out / ("STX_RESQML_XSD_%d.txt" % i)).write_text(c)
    # the parts must rebuild the same schema set, and that set must still compile
    check_dir = out / "_check"
    unpacked = unpack_parts(chunks, check_dir)
    status = check_schemas(check_dir)
    print("%d .xsd files, %d KB with documentation, %d KB archive (sha256 %s), %d secret part(s); %s"
          % (unpacked, raw // 1024, len(archive) // 1024, sha, len(chunks), status))
    print("\nSet the secrets with the GitHub CLI, from the repository folder. PowerShell:")
    for i in range(1, len(chunks) + 1):
        print("  Get-Content -Raw %s | gh secret set STX_RESQML_XSD_%d" % (out / ("STX_RESQML_XSD_%d.txt" % i), i))
    print("bash:")
    for i in range(1, len(chunks) + 1):
        print("  gh secret set STX_RESQML_XSD_%d < %s" % (i, (out / ("STX_RESQML_XSD_%d.txt" % i)).as_posix()))
    print("or paste each file's text as a secret of the same name in Settings > Secrets and variables > Actions.")
    print("Parts above %d from an earlier packing are ignored (the header gives the count). Delete %s afterwards."
          % (len(chunks), out.as_posix()))
    for f in sorted(check_dir.rglob("*"), reverse=True):
        f.unlink() if f.is_file() else f.rmdir()
    check_dir.rmdir()
    return 0


def unpack_parts(parts, folder):
    head, _, first = parts[0].partition("\n")
    fields = head.split()
    if len(fields) != 3 or fields[0] != MAGIC:
        sys.exit("STX_RESQML_XSD_1 does not start with the packing header: pack the schemas again")
    n, sha = int(fields[1]), fields[2]
    if len(parts) < n:
        sys.exit("the packing has %d parts but only %d secrets are set" % (n, len(parts)))
    b64 = "".join([first] + [p.strip() for p in parts[1:n]])
    archive = base64.b64decode("".join(b64.split()))
    if hashlib.sha256(archive).hexdigest() != sha:
        sys.exit("the parts do not rebuild the packed archive (checksum differs): set all %d parts again" % n)
    folder = pathlib.Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    count = 0
    with tarfile.open(fileobj=io.BytesIO(lzma.decompress(archive)), mode="r") as tar:
        for m in tar.getmembers():
            target = (folder / m.name).resolve()
            if not m.isfile() or not str(target).startswith(str(folder.resolve()) + os.sep) or not m.name.lower().endswith(".xsd"):
                sys.exit("unexpected member %s in the schema archive" % m.name)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(tar.extractfile(m).read())
            count += 1
    return count


def unpack(args):
    parts = []
    for i in range(1, MAX_PARTS + 1):
        v = os.environ.get("STX_RESQML_XSD_%d" % i, "")
        if not v.strip():
            break
        parts.append(v)
    if not parts:
        print("no schema secrets set (STX_RESQML_XSD_1 is empty)")
        return 3
    count = unpack_parts(parts, args.folder)
    sha = parts[0].split("\n", 1)[0].split()[2]
    print("RESQML schemas rebuilt from secrets: %d .xsd files (archive sha256 %s); %s"
          % (count, sha, check_schemas(args.folder)))
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pack", help="pack a downloaded schema package into secret parts")
    p.add_argument("package", help="the Energistics RESQML 2.0.1 package (.zip, or its unpacked folder)")
    p.add_argument("-o", "--output", default="xsd_parts", help="folder for the part files (default xsd_parts)")
    u = sub.add_parser("unpack", help="rebuild the schema folder from STX_RESQML_XSD_1..8")
    u.add_argument("folder")
    args = ap.parse_args()
    return pack(args) if args.cmd == "pack" else unpack(args)


if __name__ == "__main__":
    sys.exit(main())
