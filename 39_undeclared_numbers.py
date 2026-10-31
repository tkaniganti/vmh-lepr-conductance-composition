#!/usr/bin/env python3
"""
39_undeclared_numbers.py - which numbers in the manuscript has nobody declared?

34_manifest.py answers "does each declared number match its source?" and --scope answers
"does the verification reach each declared producer?". Both loop over the 72 declarations,
so a quoted figure that was never declared is absent from both. The scope check says so in
its own closing note, and that note sat there without anyone acting on it.

This runs the comparison the other way: pull every numeral out of the manuscript, match it
against the declaration list, and print what nothing claims. It is how the claim that the
estimator "recovers about 1% where 2.9% is present at six levels" should have been caught,
since no script in the pipeline produced either figure.

Expect noise. Years, accessions, sample sizes, figure callouts, p values and section numbers
will all appear, and most are fine. The output is a list to skim once, not a pass/fail gate.
What matters is whether anything in it is a RESULT nobody can regenerate.

Usage:
    python 39_undeclared_numbers.py v50.md
    python 39_undeclared_numbers.py v50.md --context    # show the sentence around each

Exit status is always 0. This reports; it does not judge.
"""
import re
import sys

# numbers that are obviously not results
SKIP_EXACT = {
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "20", "100",
    "0.05", "0.1", "95", "50",
}
SKIP_PATTERNS = [
    re.compile(r"^(19|20)\d{2}$"),              # years
    re.compile(r"^1\.1\.\d+$"),                 # release tags
    # The section-number pattern that used to sit here, ^\d{1,2}\.\d{1,2}$, dropped
    # every token of that shape. It was discarding about 90 distinct values, among them
    # 3.17 (the depth coefficient 2.6 rests on), 0.74 (the injection recovery bounding
    # 2.3's central null), 0.77 / 0.92 / 0.99 (the cluster-mean reliabilities) and
    # 1.77 / 1.78 (the axis extremes). Shape cannot separate a section reference from a
    # result: see is_section_reference below, which uses position instead.
]
# lines that are citations or accessions rather than claims
SKIP_LINE = re.compile(
    r"GSE\d|doi:|https?://|^\s*\d+\.\s+[A-Z][a-z]+ [A-Z]{1,3}[,.]|Zenodo|ORCID|"
    r"grants? #|PMC\d|ISSN", re.I)

NUM = re.compile(r"(?<![\w.])[−-]?\d[\d,]*\.?\d*(?:\s*×\s*10[⁻⁰¹²³⁴⁵⁶⁷⁸⁹-]+)?(?![\w])")


# ---------------------------------------------------------------- section references
# A token shaped like a section number is a REFERENCE if its position says so, and a
# RESULT otherwise. Verified by reading every colliding token in the live manuscript:
# 44 references, all genuine, and 3 results, all genuine.
CUE = re.compile(r"\b(sections?|see)\b|[#\u00a7]", re.I)
CLAUSE_END = re.compile(r"[,.;:)]")
# Bare "in" only at distance one. A fixed word window was tried: at four words it
# classified this draft correctly, at five it called "the Nfib class with 2.6-fold" a
# reference, at three it lost the last member of "in sections 2.4 and 2.5". One word of
# margin either side, set by whichever sentences happen to be in the draft.
IN_ADJACENT = re.compile(r"\bin\s*$", re.I)
JOINER = re.compile(r"^[\s]*(and|,|or)[\s]*$", re.I)
# A TRAILING UNIT OVERRIDES THE POSITION RULES. "detected in 2.4% of nuclei" has a cue
# at distance one and would otherwise be hidden, which is the under-firing direction.
# Narrow and deliberate: no section reference is ever written "section 2.4%".
UNIT_SUFFIX = re.compile(r"^(%|-[A-Za-z]|\s*\u00d7)")
SECTION_SHAPE = re.compile(r"^\d{1,2}\.\d{1,2}$")
HEADING_LINE = re.compile(r"^\s*#{1,6}\s")


def heading_numbers(text):
    """Section numbers AS STRINGS, read from the live document rather than hard-coded.

    Strings because "2.80" appears in the range on the 3.57% component and equals
    section 2.8 as a float. From the document because the next renumbering would
    otherwise reclassify every cross-reference as a result.
    """
    out = set()
    for line in text.split("\n"):
        m = re.match(r"^\s*#{1,6}\s+(\d{1,2}\.\d{1,2})\b", line)
        if m:
            out.add(m.group(1))
    return out


def is_section_reference(before, after, tok, headings, line_refs, in_heading):
    if not SECTION_SHAPE.match(tok):
        return False
    if UNIT_SUFFIX.match(after):
        return False
    if tok not in headings:
        return False
    if in_heading and not before.strip(" #"):
        return True
    if IN_ADJACENT.search(before):
        return True
    for m in reversed(list(CUE.finditer(before))):
        if not CLAUSE_END.search(before[m.end():]):
            return True
        break
    for pos, ptok in reversed(line_refs):
        gap = before[pos:]
        if JOINER.match(gap) or JOINER.match(gap.replace(ptok, "", 1)):
            return True
        break
    return False


def normalize(tok):
    t = tok.replace(",", "").replace("−", "-").strip()
    t = re.sub(r"\s*×\s*10.*$", "", t)
    try:
        return float(t)
    except ValueError:
        return None


def roundings(v):
    """Every shorter rounding of a declared value, since prose quotes 0.136 for a
    stored 0.1364. 44_manuscript_check.py accepts these in the other direction; this
    script matched exactly, so it reported declared values as unclaimed."""
    # Python rounds halves to even: round(-0.0925, 3) is -0.092, while prose writes
    # -0.093, rounding half away from zero. A declared value sitting exactly on a
    # midpoint was therefore reported as unclaimed. Generate BOTH conventions.
    from decimal import Decimal, ROUND_HALF_UP
    a = abs(v)
    out = {round(a, 4)}
    for dp in range(0, 5):
        out.add(round(round(a, dp), 4))
        q = Decimal(1).scaleb(-dp) if dp else Decimal(1)
        out.add(float(Decimal(repr(a)).quantize(q, rounding=ROUND_HALF_UP)))
    return {round(x, 6) for x in out}


def declared_values():
    """Every value the manifest declares, as floats. 34_manifest.py cannot be imported
    under that name, so its DECLARATIONS tuples are read from the source instead."""
    src = open("34_manifest.py").read()
    vals = set()
    for mm in re.finditer(r'\(\s*"(?:[^"\\]|\\.)*"\s*,\s*([−-]?[\d.]+)\s*,', src):
        v = normalize(mm.group(1))
        if v is not None:
            vals |= roundings(v)
    return vals


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)
    path = sys.argv[1]
    want_context = "--context" in sys.argv

    declared = declared_values()
    print(f">> {len(declared)} distinct values declared in 34_manifest.py")

    text = open(path).read()
    # drop the references section and anything after it
    cut = re.search(r"^##\s*References\s*$", text, re.M)
    if cut:
        text = text[:cut.start()]
        print(">> references section excluded")

    headings = heading_numbers(text)
    seen, undeclared = set(), []
    for line in text.split("\n"):
        if SKIP_LINE.search(line):
            continue
        in_heading = bool(HEADING_LINE.match(line))
        line_refs = []
        # finditer, not findall: the position rules need the text either side
        for m in NUM.finditer(line):
            raw = m.group(0).strip()
            before, after = line[:m.start()], line[m.end():m.end() + 4]
            if is_section_reference(before, after, raw, headings,
                                    line_refs, in_heading):
                line_refs.append((m.end(), raw))
                continue
            if raw in SKIP_EXACT or any(p.match(raw) for p in SKIP_PATTERNS):
                continue
            v = normalize(raw)
            if v is None or v == 0:
                continue
            key = round(abs(v), 4)
            if key in declared or key in seen:
                continue
            seen.add(key)
            undeclared.append((raw, line.strip()))

    # ---- position-weighted triage.
    # A list of 100 is long enough that nobody reads it, which is how the claim a
    # section was TITLED for ("no transcript accounts for more than about 2% of the
    # variance") sat undeclared for 200 revisions. Driving the count to zero is not the
    # goal. A number in a section title, the abstract, or the first sentence of a
    # paragraph is load-bearing and should fail the run; one in a parenthetical or a
    # sensitivity range can carry a backlog indefinitely.
    def position_of(line, raw):
        st = line.strip()
        if st.startswith("#"):
            return "HEADING"
        before = st[:st.find(raw)] if raw in st else ""
        # PARENTHESIS DEPTH, not presence. The first version asked whether "(" appeared
        # in `before` with no ")" after it, but `before` spans the whole preceding
        # paragraph, so any earlier closed parenthesis defeated the test. The POMC and
        # AgRP sentences both put their numbers inside a parenthesis that opens
        # mid-sentence, and both produced false load-bearing entries: eight of the
        # fourteen reported on v242, which is the state that makes a list unread.
        if before.count("(") > before.count(")"):
            return "parenthetical"
        # first sentence of the paragraph: no sentence-ending punctuation before it
        if not any(p_ in before for p_ in (". ", "? ", "! ")):
            return "topic sentence"
        return "body"

    abstract_end = text.find("\n## ")
    load_bearing, ordinary = [], []
    for raw, line in undeclared:
        pos = position_of(line, raw)
        in_abstract = abstract_end > 0 and text.find(line) < abstract_end
        if pos == "HEADING" or in_abstract or pos == "topic sentence":
            load_bearing.append((raw, line, "abstract" if in_abstract else pos))
        else:
            ordinary.append((raw, line))

    print(f">> {len(undeclared)} distinct numbers in the body that no declaration claims")
    print(f">> of which {len(load_bearing)} are load-bearing by position "
          f"(heading, abstract, or topic sentence)\n")
    if load_bearing:
        print("  LOAD-BEARING AND UNDECLARED:")
        for raw, line, pos in load_bearing:
            i = line.find(raw)
            a, b = max(0, i - 50), min(len(line), i + len(raw) + 50)
            print(f"    [{pos}] {raw}")
            print(f"            ...{line[a:b]}...")
        print()

    for raw, line in ordinary:
        if want_context:
            ctx = line if len(line) < 300 else line[:300] + "..."
            print(f"  {raw:>12}   {ctx}")
            print()
        else:
            # the few words either side, enough to recognise it
            i = line.find(raw)
            a, b = max(0, i - 55), min(len(line), i + len(raw) + 55)
            print(f"  {raw:>12}   ...{line[a:b]}...")

    print("\n  The list above is a backlog and can stay one. The load-bearing set is not:")
    print("  a number in a title, the abstract or a topic sentence is carrying a claim.")
    print("\n  Skim for anything that is a RESULT. Years, counts, accessions, p values,")
    print("  figure callouts and sample sizes are expected here and are usually fine.")
    print("  A result in this list is a number the archive cannot regenerate.")
    if load_bearing:
        sys.exit(1)


if __name__ == "__main__":
    main()
