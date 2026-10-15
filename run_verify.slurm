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

echo "=== clean checkout of v1.1.3 ==="
git clone --depth 1 --branch v1.1.3 \
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
    ls data | head
fi

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
