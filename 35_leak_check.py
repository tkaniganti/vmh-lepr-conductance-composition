#!/usr/bin/env python3
"""
35_leak_check.py - did the annotation layer leak into the manuscript?

Drafts are assembled from a working file that carries editorial notes in brackets. The
assembly strips those notes by dropping lines that begin with "[". Markdown is hard-wrapped,
so every note that ran onto a second line kept everything from that line onward: the opening
bracket went, the text and the closing bracket stayed. Ten such fragments reached a built
docx and were found by accident, through a sweep looking for something else.

The failure generalizes past this filter. Any transform that drops lines by a leading marker
keeps the continuation of every wrapped line, and the residue reads as broken prose rather
than as a visible artifact. A comment stripper, a blockquote filter or a diff cleaner over
hard-wrapped text fails the same way.

Three checks, none of which is a general "does this sentence belong" judgment:

  brackets    a closing bracket with no opener within the preceding 600 characters
  vocabulary  words that belong to the annotation layer and never to manuscript prose
  structure   a paragraph ending in a bare "]", which is the signature of a halved note

Reports the section each hit falls in, because a fragment in a figure legend or in Methods
reads as a methodological statement and survives a read-through, while the same fragment in
a Results paragraph is obviously broken.

Exits nonzero on any hit.

Usage: python3 35_leak_check.py path/to/manuscript.md
"""
import re
import sys

VOCAB = [
    "DECISION", "RESOLVED", "CHANGED", "REMOVED", "CHECK:",
    "declared in the manifest", "the manifest", "reviewer notes",
    "previous draft", "the rewrite", "Pre-arguing", "Recorded so",
    "gene_scan.json", "trap.json", "discreteness.json", "class_nulls.json",
    "detection.json", "results/", "34_manifest", "33_reconcile",
    "verification run", "pending note",
]


def sections(text):
    return [(m.group(1).strip(), m.start())
            for m in re.finditer(r"^#{1,4}\s+(.+)$", text, re.M)]


# A scoped whitelist for "manifest" inside Data and code availability was added and
# then removed on 3 Nov 2026. It was written to silence a hit in Methods, but that
# hit was in Methods rather than in the availability statement, so the whitelist
# could never have suppressed it; the Methods sentence was reworded instead. The
# statement's own wording is "carries a manifest", which this list's "the manifest"
# entry does not match, so the whitelist had no case to act on in either place.
# Dead code with a comment describing behaviour it does not produce is worse than
# no code: it reads as a tested exception. If a future draft writes "the manifest"
# in that section, add the scope rule back WITH a fixture case that proves it still
# fires elsewhere.


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: 35_leak_check.py manuscript.md")
    path = sys.argv[1]
    s = open(path).read()
    secs = sections(s)

    def where(i):
        last = "(before the first heading)"
        for name, pos in secs:
            if pos <= i:
                last = name
            else:
                break
        return last[:50]

    hits = 0

    print("=" * 78)
    print(f"Leak check: {path}")
    print("=" * 78)

    # ---- 1. closing brackets with no opener nearby
    print("\n1. Unbalanced closing brackets")
    found = 0
    for m in re.finditer(r"\]", s):
        seg = s[max(0, m.start() - 600):m.start()]
        if seg.count("[") <= seg.count("]"):
            found += 1
            ctx = ' '.join(s[max(0, m.start() - 120):m.start() + 2].split())
            print(f"   [{where(m.start())}]")
            print(f"      ...{ctx}")
    print(f"   {found} found")
    hits += found

    # ---- 2. annotation vocabulary
    print("\n2. Annotation vocabulary in the body")
    body = s[:s.index("## References")] if "## References" in s else s
    found = 0
    for v in VOCAB:
        for m in re.finditer(re.escape(v), body):
            found += 1
            ctx = ' '.join(body[max(0, m.start() - 70):m.start() + 70].split())
            print(f"   [{where(m.start())}] {v}")
            print(f"      ...{ctx}")
    print(f"   {found} found")
    hits += found

    # ---- 3. paragraphs ending in a bare bracket
    print("\n3. Paragraphs ending in a bare ']'")
    found = 0
    pos = 0
    for para in s.split("\n\n"):
        if para.rstrip().endswith("]") and not para.lstrip().startswith("|"):
            found += 1
            print(f"   [{where(pos)}] ...{' '.join(para.split())[-100:]}")
        pos += len(para) + 2
    print(f"   {found} found")
    hits += found

    print("\n" + "=" * 78)
    if hits:
        print(f"{hits} possible leaks. A fragment in a figure legend or in Methods reads")
        print("as a methodological statement and will survive a read-through.")
        return 1
    print("No leaks found.")
    print("\nNOTE: this checks that the annotation layer did not leak. It does not check")
    print("that the prose is correct, only that it was written on purpose.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
