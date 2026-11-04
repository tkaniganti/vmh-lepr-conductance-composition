#!/usr/bin/env python3
"""docx2md.py - extract a .docx to markdown using only the standard library.

Written because python-docx is not installed on the cluster and the checks need
a text version of the manuscript. A .docx is a zip holding word/document.xml;
paragraphs are <w:p>, text runs are <w:t>, and the paragraph style sits in
<w:pStyle w:val="..."/>. Heading styles become # marks so 39_undeclared_numbers
can find the References boundary and 54_completeness can read the structure.

Usage: python docx2md.py in.docx out.md
"""
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def main(src, dst):
    xml = zipfile.ZipFile(src).read("word/document.xml")
    root = ET.fromstring(xml)
    out = []
    for p in root.iter(f"{W}p"):
        txt = "".join(t.text or "" for t in p.iter(f"{W}t"))
        if not txt.strip():
            continue
        style = ""
        ps = p.find(f"{W}pPr/{W}pStyle")
        if ps is not None:
            style = ps.get(f"{W}val", "")
        m = re.match(r"[Hh]eading\s*(\d)", style)
        if m:
            out.append("#" * int(m.group(1)) + " " + txt.strip())
        else:
            out.append(txt.strip())
    open(dst, "w").write("\n\n".join(out) + "\n")
    print(f"  {len(out)} paragraphs, "
          f"{sum(1 for l in out if l.startswith('#'))} headings -> {dst}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
