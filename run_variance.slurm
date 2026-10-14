#!/usr/bin/env bash
#SBATCH --job-name=vmh-vc
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=logs/vmh-vc-%j.out
#SBATCH --error=logs/vmh-vc-%j.err

# Applies the model respecification to 20_variance_components.py, then runs it, then runs
# the BH correction on the TRAP enrichments. Everything that was being lost to terminated
# interactive sessions, in one unattended job.
#
# Model change: class becomes a FIXED effect and cluster-within-class the only random term,
# because six class levels are too few to estimate a variance component (the earlier
# specification gave non-positive-definite Hessians on three of twenty seeds). Adds a
# per-class decomposition of the cluster term and a model restricted to the conserved
# Lepr class, which is the comparison the manuscript actually claims.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./env.sh

echo "host=$(hostname)  started=$(date)"
cp -n 20_variance_components.py 20_variance_components.py.bak || true

# ---------------------------------------------------------------- patch (idempotent)
"$PY" - << 'EOF'
p = '20_variance_components.py'
s = open(p).read()
changed = []

old_fit = '''    r = smf.mixedlm("score ~ logdepth", df, groups=df["cls"], re_formula="1",
                    vc_formula=vcf).fit(reml=True, method="lbfgs")
    var_class = float(np.asarray(r.cov_re)[0, 0])
    var_cluster = float(r.vcomp[0]) if len(r.vcomp) else 0.0
    var_resid = float(r.scale)
    return r, var_class, var_cluster, var_resid'''
new_fit = '''    r = smf.mixedlm("score ~ logdepth + C(cls)", df,
                    groups=df["cluster_in_class"]).fit(reml=True, method="lbfgs")
    var_cluster = float(np.asarray(r.cov_re)[0, 0])
    var_resid = float(r.scale)
    return r, 0.0, var_cluster, var_resid'''
if old_fit in s:
    s = s.replace(old_fit, new_fit); changed.append("fit respecified")
elif 'C(cls)' in s:
    changed.append("fit already respecified")
else:
    raise SystemExit("FATAL: fit block not found and not already patched")

if "per_class" not in s:
    anchor = "    # Permutation null for the nested term"
    block = '''    # Cluster-level variance computed separately within each conserved class
    sc0 = V.module_scores(expr, cols, colsum, *seed_sets[0])["HminusD"]
    per_class = {}
    for cl in sorted(set(cls.tolist())):
        m = cls == cl
        if len(set(cluster[m])) < 2:
            continue
        d = pd.DataFrame({"score": sc0[m], "grp": cluster[m], "logdepth": logdepth[m]})
        try:
            mod = smf.mixedlm("score ~ logdepth", d, groups=d["grp"]).fit(reml=True,
                                                                         method="lbfgs")
            vg = float(np.asarray(mod.cov_re)[0, 0]); vr = float(mod.scale)
            per_class[cl] = round(100 * vg / (vg + vr), 3)
        except Exception:
            per_class[cl] = None
    print(">> cluster variance (%) within each conserved class:", per_class, flush=True)

    # Model restricted to the conserved Lepr class, the comparison the manuscript claims
    lep = cls == "Lepr"
    print(f">> restricted model: {int(lep.sum())} nuclei, "
          f"{len(set(cluster[lep]))} clusters in the conserved Lepr class", flush=True)
    lep_recs = []
    for pr_, cr_ in seed_sets:
        sc = V.module_scores(expr, cols, colsum, pr_, cr_)["HminusD"][lep]
        d = pd.DataFrame({"score": sc, "grp": cluster[lep], "logdepth": logdepth[lep]})
        try:
            mod = smf.mixedlm("score ~ logdepth", d, groups=d["grp"]).fit(reml=True,
                                                                         method="lbfgs")
            vg = float(np.asarray(mod.cov_re)[0, 0]); vr = float(mod.scale)
            lep_recs.append(100 * vg / (vg + vr))
        except Exception:
            pass
    if lep_recs:
        print(">> cluster variance within the Lepr class: "
              f"{sum(lep_recs)/len(lep_recs):.3f}% "
              f"(min {min(lep_recs):.3f}, max {max(lep_recs):.3f}, n={len(lep_recs)})",
              flush=True)

'''
    if anchor not in s:
        raise SystemExit("FATAL: permutation anchor not found")
    s = s.replace(anchor, block + anchor, 1); changed.append("per-class and restricted added")

if '"cluster_pct_by_class"' not in s:
    s = s.replace('"n_class_by_cluster": len(set(zip(cls.tolist(), cluster.tolist()))),',
                  '"n_class_by_cluster": len(set(zip(cls.tolist(), cluster.tolist()))),\n'
                  '           "cluster_pct_by_class": per_class,\n'
                  '           "restricted_lepr_class_cluster_pct": {\n'
                  '               "mean": round(float(np.mean(lep_recs)), 4) if lep_recs else None,\n'
                  '               "min": round(float(np.min(lep_recs)), 4) if lep_recs else None,\n'
                  '               "max": round(float(np.max(lep_recs)), 4) if lep_recs else None,\n'
                  '               "n_seeds": len(lep_recs),\n'
                  '               "note": "nine Lepr clusters within the conserved Lepr class only"},',
                  1)
    changed.append("json fields added")

s = s.replace('"model": "score ~ logdepth + (1|class) + (1|class:cluster), REML",',
              '"model": "score ~ logdepth + C(class) fixed; (1|class:cluster) random; REML",')
s = s.replace('"nested_cluster_term_LR": summ("LR_cluster_term"),\n', '')
s = s.replace('("variance_explained_pct", "nested_cluster_term_LR",\n'
              '                              "nested_cluster_term_p_max")',
              '("variance_explained_pct", "cluster_pct_by_class",\n'
              '                              "restricted_lepr_class_cluster_pct",\n'
              '                              "nested_term_permutation")')
open(p, 'w').write(s)
print("PATCH:", "; ".join(changed))
EOF

"$PY" -c "import ast; ast.parse(open('20_variance_components.py').read()); print('parses ok')"
grep -n "C(cls)\|per_class = {}\|restricted_lepr_class_cluster_pct\|cluster_pct_by_class" \
     20_variance_components.py || true

# ---------------------------------------------------------------- run the model
echo; echo "=== variance components ==="
"$PY" 20_variance_components.py 2>/dev/null

# ---------------------------------------------------------------- BH on the TRAP panel
echo; echo "=== BH correction, family = effector panel ==="
"$PY" - << 'EOF'
import csv, numpy as np, config as cfg
from scipy import stats
panel = set(cfg.PANELS['HYPER']) | set(cfg.PANELS['DEPOL'])
rows = [r for r in csv.DictReader(open('results/trap_enrichment.csv'))
        if r['gene'] in panel and r['SF1_welch_p_unadjusted'] not in ('', 'None')]
print('m =', len(rows), 'of', len(panel), 'panel genes with a computable p')
print('no computable p:', sorted(panel - {r['gene'] for r in rows}))
p = np.array([float(r['SF1_welch_p_unadjusted']) for r in rows])
q = stats.false_discovery_control(p, method='bh')
for r, pi, qi in sorted(zip(rows, p, q), key=lambda t: t[1]):
    print('%-9s log2=%7s p=%.4f q=%.4f %s' % (r['gene'], r['SF1_log2_IP_over_Sup'],
                                              pi, qi, 'PASS' if qi < 0.05 else ''))
EOF

echo; echo "finished=$(date)"
echo "results/variance_components.json:"
cat results/variance_components.json
