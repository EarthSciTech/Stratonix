#!/usr/bin/env python3
"""Self-test of tools/resqml_xsd_secrets.py (I-034) on a made-up schema package laid out like the
Energistics one: nested zip, schemas importing each other by relative path, types referred to by
prefixed name, long documentation (several secret parts), the schema namespace under two prefixes.
Checks that pack + unpack rebuild a schema set that accepts and rejects the same documents as the
original, that a stale extra part is ignored, and that a damaged part is refused."""
import io
import os
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "resqml_xsd_secrets.py"


def doc(n):
    rnd = random.Random(n)
    words = ["horizon", "fault", "grid", "property", "crs", "interpretation", "feature", "unit"]
    return " ".join(rnd.choice(words) + str(rnd.randrange(10**6)) for _ in range(4000))


EML = """<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:eml="http://www.energistics.org/energyml/data/commonv2"
  targetNamespace="http://www.energistics.org/energyml/data/commonv2" elementFormDefault="qualified">
  <xs:annotation><xs:documentation>%s</xs:documentation></xs:annotation>
  <xs:simpleType name="String64"><xs:annotation><xs:documentation>%s</xs:documentation></xs:annotation>
    <xs:restriction base="xs:string"><xs:maxLength value="64"/></xs:restriction></xs:simpleType>
  <xs:complexType name="Citation"><xs:sequence>
    <xs:element name="Title" type="eml:String64"><xs:annotation><xs:appinfo>x</xs:appinfo></xs:annotation></xs:element>
  </xs:sequence></xs:complexType>
  <xs:element name="Citation" type="eml:Citation"/>
%s
</xs:schema>
""" % (doc(1), doc(2), "\n".join('  <xs:element name="E%016x" type="eml:String64"/>' % random.Random(i).getrandbits(64)
                               for i in range(12000)))

RESQML = """<?xml version="1.0" encoding="UTF-8"?>
<xsd:schema xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:resqml2="http://www.energistics.org/energyml/data/resqmlv2"
  xmlns:eml="http://www.energistics.org/energyml/data/commonv2"
  targetNamespace="http://www.energistics.org/energyml/data/resqmlv2" elementFormDefault="qualified">
  <xsd:import namespace="http://www.energistics.org/energyml/data/commonv2" schemaLocation="../../../common/v2.0/xsd_schemas/EmlAllObjects.xsd"/>
  <xsd:annotation>
    <xsd:documentation>%s</xsd:documentation>
  </xsd:annotation>
  <xsd:complexType name="Grid"><xsd:annotation><xsd:documentation>%s</xsd:documentation></xsd:annotation>
    <xsd:sequence><xsd:element ref="eml:Citation"/><xsd:element name="Ni" type="xsd:positiveInteger"/></xsd:sequence>
  </xsd:complexType>
  <xsd:element name="IjkGridRepresentation" type="resqml2:Grid"/>
</xsd:schema>
""" % (doc(3), doc(4))

GOOD = b"""<IjkGridRepresentation xmlns="http://www.energistics.org/energyml/data/resqmlv2"
  xmlns:eml="http://www.energistics.org/energyml/data/commonv2"><eml:Citation><eml:Title>g</eml:Title></eml:Citation><Ni>3</Ni></IjkGridRepresentation>"""
BAD = GOOD.replace(b"<Ni>3</Ni>", b"<Ni>-3</Ni>")


def verdicts(folder):
    from lxml import etree
    s = etree.XMLSchema(etree.parse(str(next(pathlib.Path(folder).rglob("ResqmlAllObjects.xsd")))))
    return s.validate(etree.fromstring(GOOD)), s.validate(etree.fromstring(BAD))


def main():
    try:
        import lxml  # noqa: F401
    except ImportError:
        print("lxml not installed: skipped")
        return 77
    tmp = pathlib.Path(tempfile.mkdtemp(prefix="stx_xsd_"))
    try:
        inner = io.BytesIO()
        with zipfile.ZipFile(inner, "w") as z:
            z.writestr("energyml/data/resqmlv2/v2.0.1/xsd_schemas/ResqmlAllObjects.xsd", RESQML)
            z.writestr("energyml/data/common/v2.0/xsd_schemas/EmlAllObjects.xsd", EML)
            z.writestr("energyml/data/resqmlv2/v2.0.1/README.txt", "not a schema")
        with zipfile.ZipFile(tmp / "package.zip", "w") as z:
            z.writestr("RESQML_v2.0.1/schemas.zip", inner.getvalue())
        ref = tmp / "ref"
        with zipfile.ZipFile(io.BytesIO(inner.getvalue())) as z:
            z.extractall(ref)
        out = tmp / "parts"
        r = subprocess.run([sys.executable, str(TOOL), "pack", str(tmp / "package.zip"), "-o", str(out)],
                           capture_output=True, text=True)
        print(r.stdout.strip())
        if r.returncode:
            print(r.stderr)
            return 1
        parts = sorted(out.glob("STX_RESQML_XSD_*.txt"), key=lambda p: int(p.stem.rsplit("_", 1)[1]))
        bad = []
        if len(parts) < 2:
            bad.append("expected several parts, got %d" % len(parts))
        if any(len(p.read_text()) > 48000 for p in parts):
            bad.append("a part exceeds 48 000 characters")
        env = {k: v for k, v in os.environ.items() if not k.startswith("STX_RESQML_XSD_")}
        for i, p in enumerate(parts, 1):
            env["STX_RESQML_XSD_%d" % i] = p.read_text()
        env["STX_RESQML_XSD_%d" % (len(parts) + 1)] = "stale part of an older packing"
        r = subprocess.run([sys.executable, str(TOOL), "unpack", str(tmp / "rebuilt")], capture_output=True, text=True, env=env)
        print(r.stdout.strip())
        if r.returncode:
            bad.append("unpack failed: " + r.stderr.strip())
        else:
            want, got = verdicts(ref), verdicts(tmp / "rebuilt")
            if want != (True, False) or got != want:
                bad.append("verdicts differ: original %s, rebuilt %s" % (want, got))
            text = next((tmp / "rebuilt").rglob("ResqmlAllObjects.xsd")).read_text()
            if "documentation" in text or "xsd:import" not in text or "resqml2:Grid" not in text:
                bad.append("the rebuilt schema lost its imports or prefixes, or kept its documentation")
            if any((tmp / "rebuilt").rglob("*.txt")):
                bad.append("a file other than a schema was packed")
        env["STX_RESQML_XSD_2"] = env["STX_RESQML_XSD_2"][:-10] + "AAAAAAAAAA"
        r = subprocess.run([sys.executable, str(TOOL), "unpack", str(tmp / "damaged")], capture_output=True, text=True, env=env)
        if r.returncode == 0 or "checksum" not in (r.stderr + r.stdout):
            bad.append("a damaged part was not refused")
        print("\n".join(bad) if bad else "PASS: the secret parts rebuild an equivalent schema set")
        return 1 if bad else 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
