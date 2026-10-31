#!/bin/bash
#SBATCH -J vmh-pde3b2
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 8
#SBATCH -t 03:00:00
#SBATCH -o logs/vmh-pde3b2-%j.out
#SBATCH -e logs/vmh-pde3b2-%j.err

# Pde3b sensitivity, using the real estimators.
#
# A reimplementation of the nested partition and the per-class excess disagreed with
# scripts 22 and 23 twice, in the classes that sit at zero, which are exactly the ones
# the central claim rests on. Any further divergence would be noise about that code
# rather than about Pde3b.
#
# This instead runs the actual scripts twice, editing config.py between the two passes
# and restoring it at the end. config.py is backed up first and restored even if a run
# fails, so the repository is left as it was found.

set -e
set -o pipefail

cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh

BK="../scratch/config.py.pde3b-backup"
OUT="results/pde3b_variants"
mkdir -p "$OUT" ../scratch
cp config.py "$BK"

# restore config.py whatever happens
restore() {
    cp "$BK" config.py
    # pass 2 leaves the without-Pde3b values in results/; the manuscript's numbers
    # are the with-Pde3b ones, so put them back or the next figure build, manifest
    # check or supplementary table silently uses the wrong variant.
    if [ -f "$OUT/nested_with_pde3b.json" ]; then
        cp "$OUT/nested_with_pde3b.json"      results/nested_partition.json
        cp "$OUT/class_nulls_with_pde3b.json" results/class_nulls.json
        echo "results/ restored to the with-Pde3b variant"
        fi
    echo; echo "config.py restored from $BK"
    "$PY" -c "
import importlib, config as C
importlib.reload(C)
print('  HYPER now:', len(C.PANELS['HYPER']), 'genes;', 'Pde3b' in C.PANELS['HYPER'])"
}
trap restore EXIT

run_pass () {
    local tag="$1"
    echo; echo "=================================================================="
    echo "pass: $tag"
    echo "=================================================================="
    "$PY" -c "
import config as C
h = C.PANELS['HYPER']
print(f'  HYPER {len(h)} genes, Pde3b present: {\"Pde3b\" in h}')"
    for sc in 22_nested_partition.py 23_class_nulls.py; do
        echo "--- $sc"
        "$PY" "$sc" 2>&1 | tail -8
    done
    cp results/nested_partition.json "$OUT/nested_${tag}.json"
    cp results/class_nulls.json      "$OUT/class_nulls_${tag}.json"
    echo "  saved to $OUT/*_${tag}.json"
}

# ---- pass 1: as the manuscript stands
run_pass with_pde3b

# ---- pass 2: Pde3b removed
"$PY" - << 'EOF'
p = 'config.py'
s = open(p).read()
old = '"Kcnj11", "Abcc8", "Pde3b",'
assert old in s, "panel line not matched; config.py untouched"
open(p, 'w').write(s.replace(old, '"Kcnj11", "Abcc8",', 1))
print("Pde3b removed from config.py for pass 2")
EOF
run_pass without_pde3b

# ---- compare, and check pass 1 against the published values
echo; echo "=================================================================="
echo "COMPARISON"
echo "=================================================================="
"$PY" - << 'EOF'
import json
D = "results/pde3b_variants"
def load(kind, tag):
    return json.load(open(f"{D}/{kind}_{tag}.json"))

nw, nx = load("nested", "with_pde3b"), load("nested", "without_pde3b")
cw, cx = load("class_nulls", "with_pde3b"), load("class_nulls", "without_pde3b")

PUBLISHED = {"class": 2.6884, "cluster": 3.5668, "Dlk1": 2.9352, "Esr1": 15.6489,
             "Fezf1": 2.6859, "Foxp2": 0.1326, "Lepr": -0.0275, "Nfib": 6.2888}

def gv(d, *p):
    for k in p:
        d = d[k]
    return d

print("  FIRST: does pass 1 reproduce the manuscript?")
checks = [("class", gv(nw, "nested_all_vmh", "pct_class", "mean")),
          ("cluster", gv(nw, "nested_all_vmh", "pct_cluster_in_class", "mean"))]
for c in sorted(cw["by_class"]):
    checks.append((c, cw["by_class"][c]["excess_over_null"]))
ok = True
for name, got in checks:
    exp = PUBLISHED.get(name)
    match = exp is not None and abs(got - exp) < 0.001
    ok &= match
    print(f"    {'OK  ' if match else 'FAIL'} {name:<9} {got:>10} vs published {exp}")
print(f"\n  pass 1 reproduces the manuscript: {ok}")
if not ok:
    print("  STOP. The comparison below means nothing until this matches.")

print("\n  WHAT Pde3b CHANGES")
print(f"    {'quantity':<12}{'with':>10}{'without':>10}{'delta':>10}")
rows = [("class", gv(nw,"nested_all_vmh","pct_class","mean"),
                  gv(nx,"nested_all_vmh","pct_class","mean")),
        ("cluster", gv(nw,"nested_all_vmh","pct_cluster_in_class","mean"),
                    gv(nx,"nested_all_vmh","pct_cluster_in_class","mean"))]
for c in sorted(cw["by_class"]):
    rows.append((c, cw["by_class"][c]["excess_over_null"],
                    cx["by_class"][c]["excess_over_null"]))
for name, a, b in rows:
    print(f"    {name:<12}{a:>10}{b:>10}{b-a:>+10.4f}")

print("\n  draws landing at zero (the quantity the reimplementation got wrong)")
for c in sorted(cw["by_class"]):
    zw = cw["by_class"][c]["observed_pct"].get("n_draws_at_zero")
    zx = cx["by_class"][c]["observed_pct"].get("n_draws_at_zero")
    print(f"    {c:<9} with {zw:>3} of 20   without {zx:>3} of 20")

order_w = sorted(cw["by_class"], key=lambda c: -cw["by_class"][c]["excess_over_null"])
order_x = sorted(cx["by_class"], key=lambda c: -cx["by_class"][c]["excess_over_null"])
lep_w = cw["by_class"]["Lepr"]
lep_x = cx["by_class"]["Lepr"]
print(f"\n  class ordering unchanged: {order_w == order_x}")
print(f"    with:    {' > '.join(order_w)}")
print(f"    without: {' > '.join(order_x)}")
print(f"  Lepr inside its own null, with:    "
      f"{lep_w['observed_pct']['mean']} <= {lep_w['null_pct']['p95']} "
      f"{lep_w['observed_pct']['mean'] <= lep_w['null_pct']['p95']}")
print(f"  Lepr inside its own null, without: "
      f"{lep_x['observed_pct']['mean']} <= {lep_x['null_pct']['p95']} "
      f"{lep_x['observed_pct']['mean'] <= lep_x['null_pct']['p95']}")
EOF

echo; echo "Both variants are in $OUT. config.py is restored on exit."
