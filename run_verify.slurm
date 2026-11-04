#!/usr/bin/env bash
#
# To run this elsewhere, change four things:
#   --account, to your own Slurm account
#   --output and --error paths, to a directory you can write
#   BASE below, to the directory holding both this checkout and the cached downloads
# BASE must contain a `vmh-leptin-project/data` directory with the accessions already
# downloaded by 01_download.sh; the job links to it rather than fetching again.
#
#SBATCH --job-name=vmh-verify
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=/ocean/projects/bio260082p/kanigant/verify-%j.out
#SBATCH --error=/ocean/projects/bio260082p/kanigant/verify-%j.err

# Clean-checkout verification of the regeneration claim in Data availability.
# Clones the archived release into a fresh directory, points it at the cached downloads,
# reruns scripts 22 to 25, and compares the numbers against what the manuscript reports.
#
# Set FULL=1 to re-download the accessions as well, which also tests 01_download.sh.

set -euo pipefail
BASE=/ocean/projects/bio260082p/kanigant
VERIFY=$BASE/verify-run
FULL=${FULL:-0}

rm -rf "$VERIFY"
mkdir -p "$VERIFY"
cd "$VERIFY"

echo "=== clean checkout of v1.1.5 ==="
git clone --depth 1 --branch v1.1.5 \
    https://github.com/tkaniganti/vmh-lepr-conductance-composition.git repo
cd repo
git log --oneline -1

echo; echo "=== point at this checkout, reuse the existing conda env ==="
export PROJECT="$VERIFY/repo"
source ./env.sh
mkdir -p results figures logs

if [ "$FULL" = "1" ]; then
    echo; echo "=== re-downloading accessions (login node normally; may fail here) ==="
    bash 01_download.sh
else
    echo; echo "=== linking cached downloads ==="
    rm -rf data
    ln -s "$BASE/vmh-leptin-project/data" data
    ls data | head || true
fi

echo; echo "=== core: population definition and per-cell scores ==="
echo "--- 02_snrna_core.py"
"$PY" 02_snrna_core.py 2>&1 | tail -12

echo; echo "=== the genome-wide passes and the cluster axis (expensive block) ==="
# 41 does its own matrix pass over all 6,048 VMH nuclei across 20 draws, about
# 21 min, and nothing downstream reads its output. It is here because its seven
# declarations were previously unverified, which is how the cluster-level values
# came to be quoted from the uncorrected pseudobulk file.
for sc in 08_gene_scan.py 18_doublet_sensitivity.py 41_cluster_axis_corrected.py; do
    echo "--- $sc"
    "$PY" "$sc" 2>&1 | tail -12
done

echo; echo "=== conserved classes, discreteness, translatome ==="
# 51 belongs with 04: it replaces the power grid 04 used to write, reads the stored
# per-nucleus scores rather than the matrix, and is cheap. 04 now produces the
# observed-data results only; the grid behind section 2.4 comes from 51.
# 52 reads percell_scores.npz and cluster_axis_corrected.json, so it runs after
# 02 and 41. It emits the reliability figures 2.2 states and the ceiling 2.5
# argues from, including every derived value, computed from its own inputs.
# 27 holds the random-panel controls 2.2 and 2.3 quote. It was never in this
# wrapper, and the scope check could not say so: nothing declared from its
# output, and that check derives its source list from declarations.
# 29 produces the Esr1_1 identity figures section 2.9 quotes. Like 27 it was
# outside this wrapper, and for the same reason the scope check could not say
# so: nothing was declared from its output.
for sc in 04_discreteness.py 27_correctness_checks.py 29_esr1_1_identity.py \
          06_crosspop.py 03_trpc5_depth.py 43_axis_distribution.py 57_contrast_specificity.py 58_class_excess_specificity.py \
          51_power_grid.py \
          52_axis_reliability.py \
          53_effector_covariation.py 05_trap.py 10_conserved_classes.py; do
    echo "--- $sc"
    "$PY" "$sc" 2>&1 | tail -12
done

echo; echo "=== rerunning 22 to 25 ==="
for sc in 22_nested_partition.py 23_class_nulls.py 24_mincluster_sensitivity.py \
          25_supplementary_tables.py; do
    echo "--- $sc"
    "$PY" "$sc" 2>&1 | tail -12
done

echo; echo "=== checking the claims the manuscript makes, not point estimates ==="
"$PY" - << 'EOF'
import json, sys
R = "results"
np_ = json.load(open(f"{R}/nested_partition.json"))
cn  = json.load(open(f"{R}/class_nulls.json"))
mc  = json.load(open(f"{R}/mincluster_sensitivity.json"))
bc  = cn["by_class"]
res = []

def check(name, ok, detail):
    res.append((name, ok, detail))

# A. Draw-averaged magnitudes. These are means over 20 control-gene draws, so a tight
#    tolerance is meaningful; drift here would indicate the code differs.
for name, got, want, tol in (
        ("nested class component",        np_["nested_all_vmh"]["pct_class"]["mean"], 2.69, .25),
        ("nested cluster-in-class",       np_["nested_all_vmh"]["pct_cluster_in_class"]["mean"], 3.57, .25),
        ("Esr1 excess",                   bc["Esr1"]["excess_over_null"], 15.65, .80),
        ("Nfib excess",                   bc["Nfib"]["excess_over_null"], 6.29, .50),
        ("Nfib excess at threshold 20",   mc["by_class"]["Nfib"]["20"]["excess"], 6.89, .55)):
    check(f"[magnitude] {name}", abs(got - want) <= tol, f"{got:.3f} vs {want} (+/-{tol})")

# B. The claims themselves. Seed-insensitive: these fail only if something real changed.
check("[claim] Lepr excess inside its own null",
      bc["Lepr"]["observed_pct"]["mean"] <= bc["Lepr"]["null_pct"]["p95"],
      f'obs {bc["Lepr"]["observed_pct"]["mean"]:.3f} <= null p95 {bc["Lepr"]["null_pct"]["p95"]:.3f}')
check("[claim] Foxp2 excess inside its own null",
      bc["Foxp2"]["observed_pct"]["mean"] <= bc["Foxp2"]["null_pct"]["p95"],
      f'obs {bc["Foxp2"]["observed_pct"]["mean"]:.3f} <= null p95 {bc["Foxp2"]["null_pct"]["p95"]:.3f}')
check("[claim] Lepr and Foxp2 reach zero on multiple draws",
      bc["Lepr"]["observed_pct"]["n_draws_at_zero"] >= 5
      and bc["Foxp2"]["observed_pct"]["n_draws_at_zero"] >= 5,
      f'Lepr {bc["Lepr"]["observed_pct"]["n_draws_at_zero"]}/20, '
      f'Foxp2 {bc["Foxp2"]["observed_pct"]["n_draws_at_zero"]}/20')
check("[claim] Nfib exceeds its null by a clear margin",
      bc["Nfib"]["observed_pct"]["mean"] > bc["Nfib"]["null_pct"]["p95"],
      f'obs {bc["Nfib"]["observed_pct"]["mean"]:.3f} vs null p95 {bc["Nfib"]["null_pct"]["p95"]:.3f}')

# C. Ordering and structure. What the argument actually rests on.
order = sorted(bc, key=lambda c: -bc[c]["excess_over_null"])
check("[order] Lepr lowest of six classes", order[-1] == "Lepr", " > ".join(order))
check("[order] Nfib highest of the distributed classes",
      bc["Nfib"]["excess_over_null"] > max(bc[c]["excess_over_null"]
                                           for c in ("Dlk1", "Fezf1", "Foxp2")),
      f'Nfib {bc["Nfib"]["excess_over_null"]:.2f}')
loo = bc["Esr1"]["leave_one_cluster_out"]
check("[structure] Esr1 collapses only on Esr1_1 removal",
      loo["Esr1_1"] < 1.0 and min(v for k, v in loo.items() if k != "Esr1_1") > 10.0,
      f'minus Esr1_1 {loo["Esr1_1"]:.2f}, others >= '
      f'{min(v for k, v in loo.items() if k != "Esr1_1"):.2f}')
check("[structure] Lepr returns zero on every single-cluster removal",
      all(v == 0.0 for v in bc["Lepr"]["leave_one_cluster_out"].values()),
      f'{len(bc["Lepr"]["leave_one_cluster_out"])} removals')
check("[structure] Nfib never collapses on a single removal",
      bc["Nfib"]["loo_min"] > 2.0,
      f'range {bc["Nfib"]["loo_min"]:.2f} to {bc["Nfib"]["loo_max"]:.2f}')

bad = 0
for name, ok, detail in res:
    bad += not ok
    print(f"   {'OK  ' if ok else 'FAIL'} {name:52s} {detail}")
print()
if bad:
    print(f"{bad} of {len(res)} checks failed.")
    print("Magnitude failures alone can reflect the stochastic control-gene draws.")
    print("Claim, order or structure failures mean the archived code differs from what")
    print("produced the manuscript, and should be resolved before posting.")
else:
    print(f"All {len(res)} checks pass: magnitudes, claims, ordering and structure.")
sys.exit(1 if bad else 0)
EOF

echo; echo "=== does the verification reach every producer of a quoted number? ==="
"$PY" 34_manifest.py --scope
scope_status=$?
if [ "$scope_status" -ne 0 ]; then
    echo
    echo "A declared number has no producer in this verification. Either extend this"
    echo "file or remove the declaration; do not ship past this."
    exit "$scope_status"
fi
