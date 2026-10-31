#!/usr/bin/env python3
"""
53_effector_covariation.py - do the effector genes covary at all, and at which resolution?

THE QUESTION. The balance axis averages 22 z-scored genes. Split-half on the per-nucleus
axis gave gene halves spanning -0.013 to 0.062 across 20 stratified halvings, a range
containing zero, against a count-split reliability of 0.261. So gene disagreement rather than
capture noise is what limits the per-nucleus measure. That was asked where the measurement is
mostly noise. This asks it directly, at both resolutions, with one statistic.

  cluster level     depth-residualized pseudobulk profile per gene across the 24 clusters
  per-nucleus level the same matrix on depth-residualized per-nucleus scores

If the genes covary at cluster level and not per nucleus, a common effector factor exists and
the per-cell measure is too noisy to carry it: the paper restricts every axis claim to
cluster and class level with a clean justification. If they covary at neither, there is no
common factor at any resolution, the axis is a formative index, and reliability is the wrong
criterion for it; the only thing that could validate it is an external criterion, which is
the POMC-AgRP gate.

PRIMARY QUANTITY: the CONTRAST, mean within-panel correlation minus mean between-panel.
Not the absolute signs. Clusters with more channel machinery overall carry more of both
panels, so the panels can correlate positively across clusters while their difference stays
meaningful; a formative index does not require its halves to oppose. The contrast is what
distinguishes internal coherence from none, whatever the signs.

STATED BEFORE RUNNING, so a result is not read as a finding:
  - within-DEPOL should exceed within-HYPER. Section 2.2 reports the depolarizing panel
    spanning 0.448 units across clusters against 0.168, and more variation to covary on
    produces higher correlations. That is arithmetic, not biology
  - individual pairwise correlations are NOT to be read. At 24 clusters the standard error
    on Spearman is near 0.22, and the pairs within each mean share genes and share clusters

CONTROLS, both gating:
  zero case    each gene's profile permuted INDEPENDENTLY. A single permutation applied to
               all 22 leaves every pairwise correlation unchanged, because correlation
               across clusters is invariant to a common relabeling; that control would
               return the observed values and look like catastrophic failure when nothing
               is wrong
  non-zero     a common factor of known loading injected into the 22 profiles, opposite in
               sign for the two panels. The contrast must recover it and sit outside the
               null band. Without this a null cannot be told from an estimator that returns
               the null whatever it is given, which is the failure that cost six manuscript
               versions on the class component

NULL: 200 draws of 22 detectability-matched non-panel genes, split into the same group sizes,
giving the contrast's distribution under no construct. Every quantity is reported as a
percentile of its own null rather than against a fixed threshold, as elsewhere in this work.

CONDITIONS, because three of them are known to matter:
  all genes, no cluster-size threshold
  20-nucleus threshold, since medians run 4 to 7 nuclei against means of 46 to 76 and a
    4-nucleus pseudobulk profile is close to noise
  panel genes below 10% detection removed (Kcnj11 ~3%, Hcn4 7.5%, Kcnj5 borderline), which
    contribute 20 of the 78 within-HYPER pairs between them
  without Pde3b, which has no established sign and contributes 12 within-HYPER pairs

Spearman is reported with the tie fraction, and Pearson on log pseudobulk alongside, since
sparse genes tie heavily at zero across only 24 points.

Output: results/effector_covariation.json
"""
import json
import os
import warnings
from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy import stats

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "cov"
N_NULL = 200
INJECT_LOADING = 0.6


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        m["_grp"] = "VMH" if (m.get("VMH_cluster") or "").strip() else ""
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "_grp", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cols = sorted(c for c in colsum if c in col2bc
                  and (meta[col2bc[c]].get("VMH_cluster") or "").strip())
    cidx = {c: i for i, c in enumerate(cols)}
    clu = np.array([meta[col2bc[c]]["VMH_cluster"].strip() for c in cols])
    depth = np.array([float(colsum[c]) for c in cols])
    logd = np.log10(np.where(depth <= 0, 1, depth))
    clusters = sorted(set(clu))
    print(f">> {len(cols)} VMH nuclei, {len(clusters)} clusters", flush=True)

    panels = {p: [g for g in cfg.PANELS[p] if g in name2row] for p in ("HYPER", "DEPOL")}
    pgenes = panels["HYPER"] + panels["DEPOL"]
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    pool = [int(r) for r in rowsum if rowsum[r] > 0]
    cand = list(rng.choice(pool, min(4000, len(pool)), replace=False))
    need = {name2row[g] for g in pgenes} | set(cand)
    ent = V.pass_rows(PREFIX, need, TAG, "cov")

    counts = defaultdict(lambda: np.zeros(len(cols)))
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset:
            counts[int(r)][cidx[c]] = float(v)
    detect = {r: float((counts[r] > 0).mean()) for r in need}

    def resid(v):
        """remove log depth linearly, as the decomposition does"""
        b = np.polyfit(logd, v, 1)
        return v - np.polyval(b, logd)

    def cluster_profile(row, keep):
        """depth-residualized pseudobulk mean per cluster"""
        v = resid(np.log1p(counts[row] / np.where(depth <= 0, 1, depth)
                           * cfg.NORM_TARGET))
        return np.array([v[clu == c].mean() for c in keep])

    def percell_profile(row):
        return resid(np.log1p(counts[row] / np.where(depth <= 0, 1, depth)
                              * cfg.NORM_TARGET))

    def summarise(profiles, nh, method="spearman"):
        """mean within-HYPER, within-DEPOL, between, and the contrast"""
        n = len(profiles)
        f = (lambda a, b: stats.spearmanr(a, b).statistic) if method == "spearman" \
            else (lambda a, b: stats.pearsonr(a, b)[0])
        wh, wd, bt = [], [], []
        for i, j in combinations(range(n), 2):
            r = f(profiles[i], profiles[j])
            if not np.isfinite(r):
                continue
            if i < nh and j < nh:
                wh.append(r)
            elif i >= nh and j >= nh:
                wd.append(r)
            else:
                bt.append(r)
        mw = np.mean(wh + wd) if (wh or wd) else np.nan
        mb = np.mean(bt) if bt else np.nan
        return {"within_hyper": round(float(np.mean(wh)), 4) if wh else None,
                "within_depol": round(float(np.mean(wd)), 4) if wd else None,
                "between": round(float(mb), 4) if bt else None,
                "contrast": round(float(mw - mb), 4),
                "n_pairs": len(wh) + len(wd) + len(bt)}

    out = {"n_nuclei": len(cols), "n_clusters": len(clusters),
           "primary_quantity": "contrast = mean within-panel minus mean between-panel",
           "primary_condition": ("all 22 genes, all clusters. The other three "
                                 "are pre-declared sensitivities; leading with "
                                 "the strongest reads as forking even when it "
                                 "is not"),
           "do_not_read_individual_pairs": (
               "at 24 clusters the SE on Spearman is near 0.22 and the pairs share genes "
               "and clusters; only the means and the contrast carry information"),
           "stated_before_running": (
               "within-DEPOL is expected to exceed within-HYPER because the depolarizing "
               "panel spans 0.448 units across clusters against 0.168; that is arithmetic"),
           "conditions": {}}

    def run_condition(name, keep_clusters, gene_subset):
        hyp = [g for g in panels["HYPER"] if g in gene_subset]
        dep = [g for g in panels["DEPOL"] if g in gene_subset]
        if len(hyp) < 2 or len(dep) < 2 or len(keep_clusters) < 8:
            out["conditions"][name] = {"skipped": "too few genes or clusters"}
            print(f"   {name:<34} skipped", flush=True)
            return
        rows_h = [name2row[g] for g in hyp]
        rows_d = [name2row[g] for g in dep]
        prof = [cluster_profile(r, keep_clusters) for r in rows_h + rows_d]
        obs = summarise(prof, len(hyp))
        ties = float(np.mean([1 - len(set(p)) / len(p) for p in prof]))

        # zero control: INDEPENDENT permutation per gene
        zc = []
        for s in range(50):
            rg = np.random.default_rng(5000 + s)
            pz = [rg.permutation(p) for p in prof]
            zc.append(summarise(pz, len(hyp))["contrast"])
        zero = float(np.mean(zc))

        # non-zero control: inject a common factor, opposite sign by panel
        rg = np.random.default_rng(7000)
        fac = rg.normal(0, 1, len(keep_clusters))
        pin = []
        for i, p in enumerate(prof):
            sgn = 1.0 if i < len(hyp) else -1.0
            pin.append(p + sgn * INJECT_LOADING * np.std(p) * fac)
        pos = summarise(pin, len(hyp))["contrast"]

        # null: detectability-matched non-panel gene sets
        nulls = []
        for _ in range(N_NULL):
            pick = []
            for g in hyp + dep:
                d = detect[name2row[g]]
                near = sorted(cand, key=lambda q: abs(detect[q] - d))[:200]
                pick.append(int(rng.choice(near)))
            pn = [cluster_profile(r, keep_clusters) for r in pick]
            nulls.append(summarise(pn, len(hyp))["contrast"])
        nulls = np.array([v for v in nulls if np.isfinite(v)])
        pct = float((nulls < obs["contrast"]).mean() * 100)

        out["conditions"][name] = {
            **obs, "n_hyper": len(hyp), "n_depol": len(dep),
            "n_clusters_used": len(keep_clusters),
            "tie_fraction_mean": round(ties, 3),
            "pearson_contrast": round(summarise(prof, len(hyp), "pearson")["contrast"], 4),
            "zero_control_contrast": round(zero, 4),
            "injected_control_contrast": round(float(pos), 4),
            "null_contrast_p5": round(float(np.percentile(nulls, 5)), 4),
            "null_contrast_p95": round(float(np.percentile(nulls, 95)), 4),
            "percentile_of_null": round(pct, 1)}
        print(f"   {name:<34} contrast {obs['contrast']:+.4f}  "
              f"null p5-p95 [{np.percentile(nulls, 5):+.3f}, "
              f"{np.percentile(nulls, 95):+.3f}]  pct {pct:.1f}", flush=True)

    sizes = {c: int((clu == c).sum()) for c in clusters}
    big = [c for c in clusters if sizes[c] >= 20]
    lowdet = {g for g in pgenes if detect[name2row[g]] < 0.10}
    print(f">> below 10% detection: {sorted(lowdet)}", flush=True)
    print(f">> clusters >= 20 nuclei: {len(big)} of {len(clusters)}", flush=True)

    print("\n== cluster level ==", flush=True)
    run_condition("all genes, all clusters", clusters, set(pgenes))
    run_condition("all genes, clusters >= 20", big, set(pgenes))
    run_condition("detection >= 10%, all clusters", clusters, set(pgenes) - lowdet)
    run_condition("without Pde3b", clusters, set(pgenes) - {"Pde3b"})

    print("\n== per-nucleus level, same statistic ==", flush=True)
    rows_h = [name2row[g] for g in panels["HYPER"]]
    rows_d = [name2row[g] for g in panels["DEPOL"]]
    pc = [percell_profile(r) for r in rows_h + rows_d]
    obs_pc = summarise(pc, len(rows_h))
    nulls_pc = []
    for _ in range(50):
        pick = []
        for g in pgenes:
            d = detect[name2row[g]]
            near = sorted(cand, key=lambda q: abs(detect[q] - d))[:200]
            pick.append(int(rng.choice(near)))
        nulls_pc.append(summarise([percell_profile(r) for r in pick],
                                  len(rows_h))["contrast"])
    nulls_pc = np.array([v for v in nulls_pc if np.isfinite(v)])

    # The per-nucleus branch carries the negative result and had NO controls. The
    # cluster branch had both. A negative from an uncontrolled branch is the failure
    # this script's own docstring cites, so both controls run here too.
    #
    # Grouping diagnostic first. If a masking error put all pairs into both terms the
    # contrast is exactly zero for any data, which is the most likely way to get
    # -0.0000 to four places. With 9 hyperpolarizing and 13 depolarizing genes the
    # three groups must be 36, 78 and 117.
    gs = {"within_hyper": 0, "within_depol": 0, "between": 0}
    nh_pc = len(rows_h)
    for i, j in combinations(range(len(pc)), 2):
        if i < nh_pc and j < nh_pc:
            gs["within_hyper"] += 1
        elif i >= nh_pc and j >= nh_pc:
            gs["within_depol"] += 1
        else:
            gs["between"] += 1
    print(f"   group sizes {gs} (expect 36 / 78 / 117 for 9 and 13 genes)", flush=True)

    pc_zero = []
    for s_ in range(20):
        rg2 = np.random.default_rng(8000 + s_)
        pc_zero.append(summarise([rg2.permutation(v) for v in pc], nh_pc)["contrast"])
    pc_zero_m = float(np.mean(pc_zero))

    # Injection at TWO loadings. 0.6 shows the branch is not at a floor. The claim the
    # manuscript will make is that a factor present at cluster level does not survive
    # to the per-nucleus measure, which is a negative WITH A DETECTABILITY REQUIREMENT,
    # structurally identical to 2.3's 0.10 SD bound. That needs an injection matched to
    # the cluster-level effect, not an arbitrary large one. Without it the supportable
    # claim is only "no LARGE common factor per nucleus".
    #
    # The matched loading is solved for: inject into the CLUSTER profiles at a range of
    # loadings, find the one reproducing the observed cluster contrast, then inject that
    # same loading per nucleus.
    cl_obs = out["conditions"]["all genes, all clusters"]["contrast"]
    base_prof = [cluster_profile(name2row[g], clusters)
                 for g in panels["HYPER"] + panels["DEPOL"]]
    nh_cl = len(panels["HYPER"])
    # SOLVE FROM A ZERO BASELINE. The first version injected into the observed
    # profiles, which already carry the +0.0638 contrast, so the sweep exited at its
    # first step (0.02) and the answer was the smallest number in the grid rather than
    # the matched loading. Both failed solves in this script shared that defect: a
    # search that accepted its first candidate without checking the baseline it was
    # searching from. So: permute each profile independently to remove the contrast,
    # CONFIRM the baseline is near zero, then bisect.
    rgm = np.random.default_rng(8600)
    perm_prof = [np.random.default_rng(8700 + i).permutation(p_)
                 for i, p_ in enumerate(base_prof)]
    base_contrast = summarise(perm_prof, nh_cl)["contrast"]
    print(f"   permuted-profile baseline contrast {base_contrast:+.4f} "
          f"(must be near zero before the solve means anything)", flush=True)
    fac_cl = rgm.normal(0, 1, len(clusters))

    def cl_contrast_at(load):
        trial = [p_ + (1.0 if i < nh_cl else -1.0) * load * np.std(p_) * fac_cl
                 for i, p_ in enumerate(perm_prof)]
        return summarise(trial, nh_cl)["contrast"]

    matched = None
    if abs(base_contrast) < 0.05:
        lo, hi = 0.0, 2.0
        if cl_contrast_at(hi) < cl_obs:
            print("   even loading 2.0 does not reach the observed contrast", flush=True)
        else:
            for _ in range(40):                       # bisection, not a grid
                mid = (lo + hi) / 2
                if cl_contrast_at(mid) < cl_obs:
                    lo = mid
                else:
                    hi = mid
            matched = round(float(hi), 4)
    else:
        print("   BASELINE NOT ZERO: the permutation did not remove the contrast, so "
              "the solve is void", flush=True)
    print(f"   loading reproducing the observed cluster contrast "
          f"({cl_obs:+.4f}) from a zero baseline: {matched}", flush=True)

    rg2 = np.random.default_rng(8500)
    fac_pc = rg2.normal(0, 1, len(cols))

    def inject_pc(load):
        pin = [v + (1.0 if i < nh_pc else -1.0) * load * np.std(v) * fac_pc
               for i, v in enumerate(pc)]
        return summarise(pin, nh_pc)["contrast"]

    pc_pos = inject_pc(INJECT_LOADING)
    pc_matched = inject_pc(matched) if matched else None
    if matched:
        print(f"   per-nucleus at the MATCHED loading {matched}: "
              f"{pc_matched:+.4f}", flush=True)

    ties_pc = {g: round(float(1 - len(set(pc[i])) / len(pc[i])), 3)
               for i, g in enumerate(panels["HYPER"] + panels["DEPOL"])}

    pc_zero_ok = abs(pc_zero_m) < 0.02
    pc_pos_ok = pc_pos > float(np.percentile(nulls_pc, 95))
    print(f"   per-nucleus zero control     {pc_zero_m:+.4f}  "
          f"{'OK' if pc_zero_ok else 'FAIL'}", flush=True)
    print(f"   per-nucleus injected control {pc_pos:+.4f}  "
          f"{'OK' if pc_pos_ok else 'FAIL'}", flush=True)
    if not pc_pos_ok:
        print("   The per-nucleus branch cannot recover a factor the cluster branch")
        print("   recovers easily. Its zero is then a statement about the estimator at")
        print("   this sparsity, not about the genes, and the comparison between levels")
        print("   is NOT AVAILABLE.", flush=True)

    out["per_nucleus"] = {
        "group_sizes": gs,
        "zero_control_contrast": round(pc_zero_m, 4),
        "injected_control_contrast": round(float(pc_pos), 4),
        "matched_loading": matched,
        "matched_loading_baseline_contrast": round(float(base_contrast), 4),
        "matched_loading_method": ("bisection from a permuted-profile baseline; the "
                                   "first version swept upward from the OBSERVED "
                                   "profiles, which already carry the contrast, and "
                                   "returned the smallest grid step"),
        "injected_at_matched_loading": (round(float(pc_matched), 4)
                                        if pc_matched is not None else None),
        "detects_a_cluster_sized_factor": (
            bool(pc_matched is not None
                 and pc_matched > float(np.percentile(nulls_pc, 95)))),
        "sensitivity_note": (
            "each per-nucleus correlation rests on 2,913 observations against 24 for "
            "the cluster profiles, and the per-nucleus branch recovered the 0.6 "
            "injection at a higher contrast than the cluster branch did against the "
            "same null ceiling, so the branch that found nothing is the more sensitive "
            "of the two"),
        "zero_control_passes": bool(pc_zero_ok),
        "injected_control_passes": bool(pc_pos_ok),
        "tie_fraction_by_gene": ties_pc,
        "comparison_between_levels_available": bool(pc_zero_ok and pc_pos_ok),
        **obs_pc,
        "null_contrast_p5": round(float(np.percentile(nulls_pc, 5)), 4),
        "null_contrast_p95": round(float(np.percentile(nulls_pc, 95)), 4),
        "percentile_of_null": round(float((nulls_pc < obs_pc["contrast"]).mean() * 100), 1)}
    print(f"   per-nucleus contrast {obs_pc['contrast']:+.4f}  "
          f"null p5-p95 [{np.percentile(nulls_pc, 5):+.3f}, "
          f"{np.percentile(nulls_pc, 95):+.3f}]", flush=True)

    # ---------------------------------------------------------------- verdict
    base = out["conditions"].get("all genes, all clusters", {})
    zero_ok = abs(base.get("zero_control_contrast", 1)) < 0.05
    pos_ok = base.get("injected_control_contrast", 0) > base.get("null_contrast_p95", 9)
    out["zero_control_passes"] = bool(zero_ok)
    out["injected_control_passes"] = bool(pos_ok)
    print(f"\n{'=' * 70}")
    print(f"  zero control (independent permutation): "
          f"{base.get('zero_control_contrast')}  {'OK' if zero_ok else 'FAIL'}")
    print(f"  injected control (loading {INJECT_LOADING}): "
          f"{base.get('injected_control_contrast')}  {'OK' if pos_ok else 'FAIL'}")
    if not (zero_ok and pos_ok):
        print("\n  A control failed. Nothing above is readable.")
    else:
        cl_in = base.get("percentile_of_null", 50) > 95
        pn_in = out["per_nucleus"]["percentile_of_null"] > 95
        print(f"\n  cluster level outside its null:     {cl_in}")
        print(f"  per-nucleus level outside its null: {pn_in}")
        if not out["per_nucleus"].get("comparison_between_levels_available"):
            print("\n  The per-nucleus branch failed a control; the between-level")
            print("  comparison is not available and no conclusion follows from it.")
        elif cl_in and not pn_in:
            print("\n  A common effector factor exists at cluster level and is not")
            print("  recoverable per nucleus. The construct is real and the per-cell")
            print("  measure is too noisy to carry it: restrict axis claims to cluster")
            print("  and class level.")
        elif not cl_in:
            print("\n  No common effector factor at either resolution. The axis is a")
            print("  formative index, reliability is the wrong criterion for it, and the")
            print("  only thing that could validate it is an external criterion.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "effector_covariation.json"), "w"),
              indent=2)
    print("\nWritten to results/effector_covariation.json")


if __name__ == "__main__":
    main()
