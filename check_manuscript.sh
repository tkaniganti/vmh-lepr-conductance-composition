#!/bin/bash
# check_manuscript.sh - run every manuscript-facing check against one export.
#
# Five checks, four directions. They are separate scripts because they answer separate
# questions, and running them by hand means running three of five and assuming the rest.
#
#   33_reconcile.py           prints the stored values behind each section, for eyeballing
#   34_manifest.py            declaration -> results/     is the manifest right?
#   34_manifest.py --scope    declaration -> producer      does verification run it?
#   44_manuscript_check.py    declaration -> manuscript    does the paper say it?
#   39_undeclared_numbers.py  manuscript  -> declaration   what is unclaimed?
#   54_completeness.py        document    <-> code         what is defined and uncited?
#   35_leak_check.py          annotation layer             did the drafting leak in?
#
# Usage: bash check_manuscript.sh "VMH leptin preprint v219.md"
#
# Exit status is the number of checks that failed, so it can gate a deposit.

set -u
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh

M="${1:-}"
if [ -z "$M" ]; then
    echo "usage: bash check_manuscript.sh <manuscript.md>"
    echo
    echo "exports present in the project root:"
    ls -1 *.md 2>/dev/null | grep -v '^README' | sed 's/^/  /'
    exit 2
fi
if [ ! -f "$M" ]; then
    echo "not found: $M"
    exit 2
fi

# More than one export in the root is the hazard these checks cannot see: every one of
# them takes a path, so whichever is named is what gets checked.
N_MD=$(ls -1 *.md 2>/dev/null | grep -vc '^README')
if [ "$N_MD" -gt 1 ]; then
    echo "WARNING: $N_MD manuscript files in the project root. Checking only:"
    echo "  $M"
    ls -1 *.md | grep -v '^README' | grep -vF "$M" | sed 's/^/  stale? /'
    echo
fi

FAILED=0
run () {
    local name="$1"; shift
    echo
    echo "################################################################"
    echo "## $name"
    echo "################################################################"
    "$@"
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo ">> $name EXITED $rc"
        FAILED=$((FAILED + 1))
    fi
    return 0
}

echo "Checking: $M"
echo "Date:     $(date -u '+%Y-%m-%d %H:%M UTC')"

run "manifest: declarations against results/"      "$PY" 34_manifest.py
run "scope: does verification reach each producer" "$PY" 34_manifest.py --scope
run "manuscript: does the paper state each declared value" \
                                                   "$PY" 44_manuscript_check.py "$M"
run "completeness: defined and uncited, cited and undefined" \
                                                   "$PY" 54_completeness.py "$M"
run "leak: did the annotation layer reach the body" "$PY" 35_leak_check.py "$M"
run "undeclared: numbers the manuscript states that nothing claims" \
                                                   "$PY" 39_undeclared_numbers.py "$M"

echo
echo "################################################################"
if [ $FAILED -eq 0 ]; then
    echo "## All checks passed."
else
    echo "## $FAILED check(s) failed. See above."
fi
cat <<'NOTE'
##
## What none of these can do:
##   - tell a number read from a file from one computed while writing the sentence.
##     Both look identical on the page, and several errors in this project were of
##     that kind. The control is reading the live document before asserting about it.
##   - catch a declaration that passes on a coincidence. 44_manuscript_check matches
##     values, so a declared 0.9288 absent from the text passes if any other number
##     rounds to it. Distinctive values are checked; round ones are not.
##   - judge whether a sentence means what it says.
################################################################
NOTE
exit $FAILED
