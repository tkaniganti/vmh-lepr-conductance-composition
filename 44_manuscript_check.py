#!/usr/bin/env python3
"""
44_manuscript_check.py - does the manuscript say what the manifest declares?

THE GAP THIS CLOSES. 34_manifest.py compares a value typed into its own DECLARATIONS list
against the value stored in results/. It never opens the manuscript. Its closing line,
"Every quoted figure matches the object the manifest declares for it", is therefore not
what it sounds like: it means the manifest agrees with the pipeline, not that the paper
does. The step from manifest to manuscript is where six of the eight numerical errors in
the v49 review happened, and nothing checks it.

Today made that concrete. Seven cluster-level values were declared at four decimals
(-0.3215, 0.4477) and the manuscript quotes them at three (-0.322, 0.448). The manifest
reported 78 of 78 agreeing while 39_undeclared_numbers.py listed all seven as unclaimed,
because one compares against results/ and the other matches literal digits in the text.

WHAT THIS DOES. For each declaration, look for its value in the manuscript at the precision
the manuscript uses. A declared value is REACHED if the text contains it at full precision
or at any shorter rounding of it. Anything not reached is either absent from the paper or
written differently from what the manifest claims, and both need a human.

Rounding is accepted DELIBERATELY: a manuscript that writes 0.448 for a stored 0.4477 is
correct, and demanding four decimals in prose would be worse writing. What is not accepted
is a value that rounds to something else, which is the transcription error this is for.

Run it with the same manuscript path you give 35_leak_check.py and
39_undeclared_numbers.py. Those three answer different questions:

    34_manifest.py            declaration  vs  results/      (is the manifest right?)
    44_manuscript_check.py    declaration  vs  manuscript    (does the paper say it?)
    39_undeclared_numbers.py  manuscript   vs  declaration   (what is unclaimed?)

Exit status is 1 if any declaration is unreached, so it can gate a submission.
"""
import importlib.util
import re
import sys


def load_declarations():
    """DECLARATIONS from 34_manifest.py, which cannot be imported by that name."""
    sp = importlib.util.spec_from_file_location("m34", "34_manifest.py")
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    for attr in ("DECLARATIONS", "MANIFEST"):
        if hasattr(m, attr):
            return getattr(m, attr)
    raise SystemExit("no DECLARATIONS or MANIFEST in 34_manifest.py")


def renderings(v):
    """How this value might legitimately appear in prose."""
    out = set()
    if isinstance(v, bool):
        return out
    if isinstance(v, int):
        out.add(str(v))
        out.add(f"{v:,}")
        return out
    s = f"{v}"
    out.add(s)
    out.add(s.lstrip("-"))
    from decimal import Decimal, ROUND_HALF_UP
    for dp in range(0, 7):
        # f-strings round halves to even; prose rounds half away from zero. Generate
        # both so a declared value on a midpoint is not reported as missing.
        q = Decimal(1).scaleb(-dp) if dp else Decimal(1)
        hu = Decimal(repr(v)).quantize(q, rounding=ROUND_HALF_UP)
        for r in (f"{v:.{dp}f}", f"{hu:.{dp}f}"):
            if float(r) != 0 or v == 0:
                out.add(r)
                out.add(r.lstrip("-"))
                if "." in r:
                    out.add(r.rstrip("0").rstrip("."))
    # percentages quoted without the leading zero, and comma forms
    extra = set()
    for t in out:
        if t.startswith("0."):
            extra.add(t[1:])
        try:
            f = float(t)
            if abs(f) >= 1000:
                extra.add(f"{f:,.0f}")
        except ValueError:
            pass
    return {t for t in (out | extra) if t not in ("", "-", "0")}


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)
    path = sys.argv[1]
    text = open(path).read()
    # normalise the minus signs and thousands separators the doc exports
    norm = text.replace("\u2212", "-").replace("\u2013", "-")

    decls = load_declarations()
    print(f">> {len(decls)} declarations, checked against {path}\n")

    unreached, not_quoted = [], []
    for e in decls:
        label, value = e[0], e[1]
        hit = None
        for cand in sorted(renderings(value), key=len, reverse=True):
            # prose writes 243rd and 98.0th; a bare \w boundary rejects both
            if re.search(rf'(?<![\w.]){re.escape(cand)}(?:st|nd|rd|th)?(?![\w])', norm):
                hit = cand
                break
        proc = e[4] if len(e) > 4 else ""
        if hit is None:
            # NOT QUOTED marks a value the paper deliberately never states: an
            # intermediate, or a count at a different stage of filtering. It still has
            # to match results/, which 34_manifest.py checks.
            if "NOT QUOTED" in proc.upper():
                not_quoted.append((label, value))
            else:
                unreached.append((label, value, e[2] if len(e) > 2 else "", proc))

    if not_quoted:
        print(f"  {len(not_quoted)} marked NOT QUOTED, absence expected:")
        for lb, v in not_quoted:
            print(f"    {lb:<42} {v}")
        print()

    if not unreached:
        print(f"  All {len(decls) - len(not_quoted)} quoted declarations appear in "
              f"the manuscript,")
        print("  at full precision or at a shorter rounding of it.")
        sys.exit(0)

    print(f"  {len(unreached)} of {len(decls)} declared values do NOT appear:\n")
    for label, value, src, proc in unreached:
        print(f"    {label}")
        print(f"      declared {value}   from {src}")
        if proc:
            print(f"      {proc[:110]}")
        print()

    print("  Each is one of three things:")
    print("    - a number the manuscript no longer quotes, so the declaration is stale")
    print("    - a number written differently from what the manifest claims, which is")
    print("      the transcription error this check exists to find")
    print("    - a rounding that does not round to the declared value, same thing")
    print("\n  34_manifest.py cannot see any of these: it never opens the manuscript.")
    sys.exit(1)


if __name__ == "__main__":
    main()
