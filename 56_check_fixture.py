#!/usr/bin/env python3
"""
56_check_fixture.py - do the checks fire on defects we know are there?

WHY THIS EXISTS. Every check in this project has a null hypothesis of "fine" and needs
positive evidence to fire. None has ever been given a defect and asked to find it. The
standing rule for estimators is that each gets a positive control including a zero case,
kept running in the pipeline rather than exercised once during development. The checks had
no such control, and the consequences were not hypothetical:

  54_completeness.py    its consumer search counted a script's mentions of its OWN output
                        as consumption. Every file is named by whatever writes it, so the
                        orphan list was structurally near-empty. Found by accident, when
                        esr1_1_markers.csv turned out to be written by 28 and read by
                        nothing while the check called it consumed
  34_manifest.py --scope derives its source list from declarations, so a producer with no
                        declarations is invisible however many numbers the paper takes
                        from it. Four scripts (27, 29, 43, 52) sat outside run_verify.slurm
                        all day and no check could say so
  44_manuscript_check.py matches VALUES with a rounding tolerance, so a declared 0.2677
                        passed against a manuscript saying 0.272, because 0.27 is a shorter
                        rounding of both. The tolerance was built to permit prose rounding
                        and permits staleness as a side effect
  39_undeclared_numbers.py listed 106 items, which nobody reads, so the claim a section was
                        TITLED for sat undeclared for 200 revisions

Each of those is a check permitting rather than missing, which is the same defect the class
component estimator had before anyone tested it.

WHY A FIXTURE RATHER THAN PLANT-AND-REMOVE. Planting a defect, watching the check fire, and
removing it tests the check once, on the day it is built. It leaves no evidence a month
later and nothing for the next person. A planted defect left in the live tree is a real
defect. So: a small tree of deliberately broken artifacts that the checks run against as a
SECOND target. Each check must fire on the fixture and stay quiet on the live tree. Both
halves matter; a check that fires on everything is as useless as one that fires on nothing.

SEVEN PLANTED DEFECTS, each modelled on something that actually happened here:

  1  declared value off in the third decimal        the 0.272 case
  2  file written, never read, writer's own
     mentions left in place                         the esr1_1_markers.csv case
  3  supplementary table defined and uncited        the S5 case
  4  panel cited in the text with no legend entry   the Figure 2C case
  5  panel in the legend and the text that the
     script no longer draws                         the Figure 2B case, the other direction
  6  undeclared number in a section title           the 2% case, which tests the position
                                                    triage specifically
  7  declared value whose producer is absent from
     the verification wrapper                       the scripts 27/29/43/52 case

Usage: python 56_check_fixture.py          build the fixture and run every check against it
       python 56_check_fixture.py --keep   leave the fixture in place afterwards

Exit status is the number of checks that failed to behave as specified.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

FIXTURE = "check_fixture"


def build(d):
    """Write a small broken tree. Every defect here is deliberate."""
    os.makedirs(os.path.join(d, "results"), exist_ok=True)

    # ---- results the fixture's manifest will point at
    json.dump({"alpha": {"value": 0.2677}, "beta": 12.5},
              open(os.path.join(d, "results", "fixture_a.json"), "w"), indent=2)
    # defect 2: written by its own script, read by nothing
    json.dump({"orphan": True},
              open(os.path.join(d, "results", "fixture_orphan.json"), "w"), indent=2)
    json.dump({"gamma": 7.77},
              open(os.path.join(d, "results", "fixture_unreached.json"), "w"), indent=2)

    # ---- the script that writes the orphan, with its own mentions intact.
    # This is the case the first consumer search called consumed.
    open(os.path.join(d, "90_fixture_producer.py"), "w").write(
        '"""Writes results/fixture_orphan.json. Nothing reads it."""\n'
        'import json, os\n'
        'json.dump({"orphan": True},\n'
        '          open(os.path.join("results", "fixture_orphan.json"), "w"))\n')
    open(os.path.join(d, "91_fixture_unreached.py"), "w").write(
        '"""Writes results/fixture_unreached.json. Not in the wrapper."""\n'
        'import json, os\n'
        'json.dump({"gamma": 7.77},\n'
        '          open(os.path.join("results", "fixture_unreached.json"), "w"))\n')

    # ---- figure script: draws A and B for Figure 1. NOT C.
    open(os.path.join(d, "19_manuscript_figures.py"), "w").write(
        'def fig1():\n'
        '    ax[0].set_title("A   First panel", loc="left")\n'
        '    ax[1].set_title("B   Second panel", loc="left")\n'
        '    # defect 5: the legend and text both claim D, which is not drawn here\n')

    # ---- supplementary builder: produces S1 and S2. S2 is never cited.
    open(os.path.join(d, "25_supplementary_tables.py"), "w").write(
        'sheet(wb, "S1 first", "Supplementary Table S1.", "", [], [])\n'
        'sheet(wb, "S2 second", "Supplementary Table S2.", "", [], [])\n')

    # ---- the broken manuscript
    open(os.path.join(d, "fixture_manuscript.md"), "w").write(
        "# Fixture\n\n"
        "## 1.2 A section whose title carries 4.44 undeclared\n\n"
        "The first result was 0.272 (Figure 1A).\n\n"
        "A second value of 12.5 appears here, and it is declared correctly.\n\n"
        "The text cites Figure 1C, which has no legend entry, and Figure 1D, "
        "which the legend defines but the script does not draw.\n\n"
        # defect 8: a colliding token with a cue five words back, carrying a unit
        # suffix. This sentence was in the live manuscript yesterday and has since
        # been rewritten away, which is why it belongs here: the manuscript is a
        # moving target and anything validated only against it is validated against
        # whatever it says this hour. Tests the clause bound, the distance-one rule
        # and the suffix guard at once.
        "The estimator detects structure in the second class with 1.2-fold poorer "
        "sampling, and a detection rate of 1.1% in that class.\n\n"
        # defect 9: a parenthesis opening mid-sentence, after an earlier closed one.
        # position_of must call 6.66 a parenthetical, not a topic sentence. Nothing
        # tested that function before this case; it is taken from the POMC passage,
        # which produced eight false load-bearing entries on v242.
        "A first clause (already closed) and then the measurement itself "
        "(reported as 6.66 here) which is parenthetical.\n\n"
        "Supplementary Table S1 is cited. S2 is not.\n\n"
        "**Figure 1. A fixture figure.** (A) First panel. (B) Second panel. "
        "(D) A panel the script does not draw.\n")

    # ---- the fixture manifest
    open(os.path.join(d, "34_manifest.py"), "w").write(
        'DECLARATIONS = [\n'
        '    ("alpha", 0.2677, "fixture_a.json", ["alpha", "value"],\n'
        '     "defect 1: the manuscript says 0.272, which 0.2677 does not round to at "\n'
        '     "three decimals, but 0.27 is a shorter rounding of both"),\n'
        '    ("beta", 12.5, "fixture_a.json", ["beta"], "correctly stated"),\n'
        '    ("gamma", 7.77, "fixture_unreached.json", ["gamma"],\n'
        '     "defect 7: 91_fixture_unreached.py is not in the wrapper"),\n'
        ']\n')

    # ---- a wrapper that omits 91
    open(os.path.join(d, "run_verify.slurm"), "w").write(
        "#!/bin/bash\n"
        "for sc in 90_fixture_producer.py; do\n"
        '    "$PY" "$sc"\n'
        "done\n")


EXPECT = [
    ("completeness: uncited supplementary table",
     "54_completeness.py", ["S2"],
     "defect 3: S2 is produced and never cited"),
    ("completeness: panel cited with no legend entry",
     "54_completeness.py", ["1C"],
     "defect 4: the text cites Figure 1C and the legend does not define it"),
    ("completeness: panel in legend the script does not draw",
     "54_completeness.py", ["1D"],
     "defect 5: the legend and text both carry D; the script draws A and B"),
    ("completeness: file written and never read",
     "54_completeness.py", ["fixture_orphan.json"],
     "defect 2: 90_fixture_producer.py writes it and nothing else names it. The "
     "first consumer search called this consumed, because the writer names it"),
    ("undeclared: colliding token with a unit suffix is a RESULT",
     "39_undeclared_numbers.py", ["1.2"],
     "defect 8: 1.2-fold collides with a heading number and has a cue five words "
     "back. The suffix makes it a result. If this is silent, the position rules are "
     "hiding a number, which is the under-firing direction"),
    ("undeclared: parenthetical is NOT load-bearing",
     "39_undeclared_numbers.py", ["6.66"],
     "defect 9: 6.66 sits inside a parenthesis that opens after an earlier closed "
     "one. It must appear in the backlog and NOT in the load-bearing set. This "
     "assertion only checks it is extracted; read the output to confirm placement"),
    ("undeclared: number in a section title",
     "39_undeclared_numbers.py", ["4.44"],
     "defect 6: 4.44 sits in a heading and no declaration claims it. Tests the "
     "position triage, which is the newest part and the least exercised"),
]


def run(script, target, cwd):
    """Run a check inside the fixture, return its combined output."""
    env = dict(os.environ)
    env["PYTHONPATH"] = cwd
    try:
        p = subprocess.run([sys.executable, script, target], cwd=cwd, env=env,
                           capture_output=True, text=True, timeout=180)
        return (p.stdout or "") + (p.stderr or "")
    except Exception as e:                                  # noqa: BLE001
        return f"RUN FAILED: {e}"


def main():
    keep = "--keep" in sys.argv
    here = os.getcwd()
    d = os.path.join(here, FIXTURE)
    if os.path.exists(d):
        shutil.rmtree(d)
    build(d)
    # the checks themselves are copied in, so the fixture exercises the real code
    for f in ("54_completeness.py", "39_undeclared_numbers.py"):
        if os.path.exists(f):
            shutil.copy(f, d)
        else:
            print(f"  {f} not found in {here}; cannot test it")
    print(f">> fixture built at {d}\n")

    outputs = {}
    failures = 0
    for name, script, needles, why in EXPECT:
        if script not in outputs:
            outputs[script] = run(script, "fixture_manuscript.md", d)
        out = outputs[script]
        hit = all(n in out for n in needles)
        print(f"  [{'FIRES' if hit else 'SILENT'}] {name}")
        if not hit:
            failures += 1
            print(f"           expected to see {needles} in {script}'s output")
            print(f"           {why}")
            print(f"           A CHECK THAT DOES NOT FIRE ON A KNOWN DEFECT HAS NEVER")
            print(f"           BEEN SHOWN TO WORK.")

    print()
    print("  Defects 1 and 7 need the manifest and the manuscript check, which read")
    print("  config and vmhlib; run them against the fixture by hand:")
    print(f"    cd {FIXTURE} && python 44_manuscript_check.py fixture_manuscript.md")
    print("  Expected: alpha reported as not appearing, because the manuscript says")
    print("  0.272 and the declaration says 0.2677. If it passes, the rounding")
    print("  tolerance is permitting staleness, which is defect 1.")
    print(f"    cd {FIXTURE} && python 34_manifest.py --scope")
    print("  Expected: fixture_unreached.json reported as NOT REACHED.")

    print(f"\n{'=' * 68}")
    if failures:
        print(f"  {failures} check(s) did not fire on a defect that is there.")
    else:
        print("  Every check tested fired on its planted defect.")
    print("  The other half of the control is that these same checks stay quiet on")
    print("  the live tree, which the normal run shows.")

    if not keep:
        shutil.rmtree(d)
        print(f"\n  fixture removed; --keep leaves it in place")
    sys.exit(failures)


if __name__ == "__main__":
    main()
