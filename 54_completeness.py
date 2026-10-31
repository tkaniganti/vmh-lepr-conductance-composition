#!/usr/bin/env python3
"""
54_completeness.py - what is defined and never used, and used and never defined?

THE GAP THIS FILLS. Every existing check runs from cited things to their sources:
34_manifest.py takes declarations and verifies them against results/, --scope takes
declarations and verifies their producers run, 44_manuscript_check.py takes declarations and
finds them in the text. A check that loops over what is cited cannot find what nobody cites.

Three things were found by hand that this would have caught:

  S5    produced by 25_supplementary_tables.py, defined in the Methods list, cited nowhere
        in the body for many versions
  S3    same shape, surfaced only when a Discussion sentence that cited it was removed
  Fig 2 the partition panel was deleted from 19_manuscript_figures.py and the text kept
        citing Figure 2B; meanwhile the contrasts panel became B and the text still said C.
        So the text cited a panel that did not exist AND mislabelled one that did

FOUR DIRECTIONS, because the Figure 2 problem ran both ways:

  1. supplementary tables the builder produces that the text never cites
  2. figure panels the figure script draws that the text never cites, and panels the text
     cites that the script does not draw
  3. panels named in a figure legend with no matching panel in the script, and panels the
     script draws with no legend entry
  4. results keys the pipeline emits that nothing in the manuscript or the manifest refers
     to. This is the direction no existing check covers, and it is where a quantity computed
     and forgotten would sit

WHAT THIS CANNOT DO, and it matters because having it invites the opposite belief. It checks
the document against the code. It cannot check a statement ABOUT the document made from
memory of what an edit was meant to do. Several errors in this project were of that kind:
a count inferred while writing a sentence about precision, a report that a sentence was
untouched when the edit had removed it, a sweep described as complete that had landed on one
paragraph. The file was internally consistent every time; the claim about it was not. The
only control for that is reading the live document before asserting about it.

Exit status is 1 if anything is defined-and-uncited or cited-and-undefined, 0 otherwise.

Usage: python 54_completeness.py <manuscript.md>
"""
import os
import re
import sys

FIGSCRIPT = "19_manuscript_figures.py"
# The release keeps this under supplementary/ while the working tree has it at
# the top level, so look in both. From a clean extract the top-level-only
# version found no produced tables and reported all nine as cited-and-not-
# produced: nine false failures for anyone running the checks from the archive.
TABSCRIPT = next((p_ for p_ in ("25_supplementary_tables.py",
                                "supplementary/25_supplementary_tables.py")
                  if os.path.exists(p_)), "25_supplementary_tables.py")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)
    path = sys.argv[1]
    if not os.path.exists(path):
        print(f"not found: {path}")
        sys.exit(2)
    text = open(path).read()

    # legends are definitions, not citations: split them off so a panel named only in
    # its own legend does not count as cited
    legend_spans = [m.span() for m in
                    re.finditer(r"\*\*Figure \d[^*]*\*\*.*?(?=\n\n|\Z)", text, re.S)]
    body = text
    for a, b in sorted(legend_spans, reverse=True):
        body = body[:a] + " " * (b - a) + body[b:]
    legends = {}
    for m in re.finditer(r"\*\*Figure (\d)[^*]*\*\*(.*?)(?=\n\n|\Z)", text, re.S):
        legends[m.group(1)] = m.group(2)

    problems = []

    # ---- 1. supplementary tables
    produced = set(re.findall(r'sheet\(wb,\s*"(S\d+)', open(TABSCRIPT).read())) \
        if os.path.exists(TABSCRIPT) else set()
    cited = set(re.findall(r"\bS(\d+)\b", body))
    cited = {f"S{n}" for n in cited}
    print(f"{'=' * 70}\n1. SUPPLEMENTARY TABLES\n{'=' * 70}")
    print(f"  produced by {TABSCRIPT}: {', '.join(sorted(produced)) or 'none found'}")
    for t in sorted(produced - cited):
        print(f"  DEFINED, NEVER CITED IN THE BODY: {t}")
        problems.append(f"table {t} produced and uncited")
    for t in sorted(cited - produced):
        print(f"  CITED, NOT PRODUCED: {t}")
        problems.append(f"table {t} cited and not produced")
    if not (produced - cited) and not (cited - produced):
        print("  every produced table is cited and every cited table is produced")

    # ---- 2 and 3. figure panels
    print(f"\n{'=' * 70}\n2. FIGURE PANELS\n{'=' * 70}")
    drawn = {}
    if os.path.exists(FIGSCRIPT):
        src = open(FIGSCRIPT).read()
        for fm in re.finditer(r'def (fig(\d))\(\):(.*?)(?=\ndef |\Z)', src, re.S):
            num = fm.group(2)
            drawn[num] = set(re.findall(r'set_title\(\s*(?:f?")([A-Z])\s', fm.group(3)))
            drawn[num] |= set(re.findall(r'set_title\(\s*"([A-Z])\s+', fm.group(3)))
    for num in sorted(set(drawn) | set(legends)):
        d = drawn.get(num, set())
        in_legend = set(re.findall(r"\(([A-Z])\)", legends.get(num, "")))
        in_body = set(re.findall(rf"Figure\s*{num}([A-Z])\b", body))
        print(f"  Figure {num}: script draws {sorted(d) or '-'}, "
              f"legend defines {sorted(in_legend) or '-'}, "
              f"text cites {sorted(in_body) or '-'}")
        for p in sorted(in_body - d):
            print(f"    TEXT CITES {num}{p}, SCRIPT DOES NOT DRAW IT")
            problems.append(f"Figure {num}{p} cited, not drawn")
        for p in sorted(in_body - in_legend):
            print(f"    TEXT CITES {num}{p}, LEGEND DOES NOT DEFINE IT")
            problems.append(f"Figure {num}{p} cited, not in legend")
        for p in sorted(in_legend - d):
            print(f"    LEGEND DEFINES {num}{p}, SCRIPT DOES NOT DRAW IT")
            problems.append(f"Figure {num}{p} in legend, not drawn")
        for p in sorted(d - in_legend):
            print(f"    SCRIPT DRAWS {num}{p}, LEGEND DOES NOT DEFINE IT")
            problems.append(f"Figure {num}{p} drawn, not in legend")

    # ---- 4. results keys nothing refers to
    print(f"\n{'=' * 70}\n4. RESULTS FILES NOTHING REFERS TO\n{'=' * 70}")
    # Search EVERY consumer, not three. The first version missed the figure script and
    # so listed trpc5_depth.json and trpc5_depth_null.csv as orphans although Figure 4C
    # reads both. A list that over-reports trains the reader to skim it, which is worse
    # than not having it.
    # A FILE'S OWN WRITER DOES NOT COUNT AS A CONSUMER. The first version concatenated
    # every .py in the project, including the script that writes each file, and every
    # script names its own output in its docstring and its write call. So almost
    # nothing could be an orphan and the list was structurally near-empty: it reported
    # three files, which were the anomaly rather than the signal. Checked against a
    # planted case (28_esr1_1_markers.py writes esr1_1_markers.csv and nothing reads
    # it) the first version said it was consumed.
    #
    # So: read each script separately, and when testing whether file F is consumed,
    # skip any script that writes F.
    sources = {}
    for fn in sorted(f for f in os.listdir(".") if f.endswith(".py")):
        try:
            sources[fn] = open(fn).read()
        except OSError:
            pass

    def writers_of(fname):
        """Scripts that WRITE this file. The test must be line-local: an earlier
        version asked whether the script mentioned the file anywhere AND contained
        open( anywhere, which nearly every script does, so every consumer was marked
        a writer and almost everything looked orphaned. Checked against a case known
        to be consumed (trpc5_depth.json, read by Figure 4C) that version listed it."""
        out = set()
        W = ('"w"', "'w'", "to_csv", "dump(", "savefig", "writer")
        for fn, src in sources.items():
            for line in src.splitlines():
                if fname in line and any(w in line for w in W):
                    out.add(fn)
                    break
        return out

    def consumed_by_another(fname):
        w = writers_of(fname)
        for fn, src in sources.items():
            if fn in w:
                continue
            if f'"{fname}"' in src or f"'{fname}'" in src:
                return fn
        return None

    manifest_src = sources.get("34_manifest.py", "")
    recon_src = sources.get("33_reconcile.py", "")

    # this project's results directory is shared with other work; those files are not
    # orphans of this manuscript and are reported separately rather than mixed in
    OTHER_WORK = ("sigma1_", "ap_subunits", "ap_share")

    if os.path.isdir("results"):
        orphan, other = [], []
        for f in sorted(os.listdir("results")):
            if not f.endswith((".json", ".csv")):
                continue
            if f in text or f in manifest_src or f in recon_src:
                continue
            if consumed_by_another(f):
                continue
            (other if f.startswith(OTHER_WORK) else orphan).append(f)
        if other:
            print(f"  {len(other)} file(s) belong to other work sharing this results")
            print("  directory and are not listed.")
        if orphan:
            print("  produced, but named in neither the manifest, the reconcile script")
            print("  nor the manuscript:")
            for f in orphan:
                print(f"    {f}")
            print("  These are not errors. A file can be an intermediate or a sensitivity")
            print("  that is deliberately unquoted. The list is here because a quantity")
            print("  computed and forgotten looks exactly like one of those.")
        else:
            print("  every results file is named somewhere")

    print(f"\n{'=' * 70}")
    if problems:
        print(f"  {len(problems)} completeness problem(s):")
        for p in problems:
            print(f"    {p}")
    else:
        print("  No completeness problems found.")
    print()
    print("  This checks the document against the code. It cannot check a claim ABOUT")
    print("  the document made from memory rather than from reading it, which is a")
    print("  different error class and has no automated control.")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
