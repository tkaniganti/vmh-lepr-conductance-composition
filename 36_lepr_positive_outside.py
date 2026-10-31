#!/usr/bin/env python3
"""
36_lepr_positive_outside.py - do Lepr-detected nuclei outside the Lepr clusters resemble
the Lepr clusters on the effector axis, or their own host class?

Lepr is detected in 23.5% of nuclei inside the nine Lepr clusters and 5.9% of VMH neurons
outside them. The outside group was never scored separately: those nuclei sit in the Esr1,
Nfib, Dlk1, Fezf1 and Foxp2 classes and contribute to those classes' scores.

The question this answers: if receptor expression rather than class identity tracked
effector composition, Lepr-positive nuclei outside the Lepr clusters should sit closer to
the Lepr clusters on the balance axis than their classmates do. If class identity is what
organizes composition, they should look like the class they are in.

READ THIS BEFORE USING THE RESULT. The grouping is a per-cell detection call, which section
2.1 of the manuscript argues cannot define a population: a nucleus scoring Lepr-positive
outside the Lepr clusters is some mixture of genuine expression, ambient contamination and
chance. At 5.9% detection the group is small, roughly 185 nuclei across five classes, and
is biased toward deeply sequenced nuclei, which the depth comparison below reports. This is
a check on what the gate selects, not a result about a population.

The comparison is made WITHIN each host class, so it is not confounded by between-class
differences in composition, and a depth-matched comparison is run alongside because
Lepr-positive nuclei are deeper by construction.

Output: results/lepr_positive_outside.json
"""
import os, json, warnings
from collections import defaultdict
import numpy as np
from scipy import stats
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "lpo"
MIN_PER_GROUP = 10          # a host class below this in either group is reported, not tested
N_BOOT = 2000


def boot_diff(a, b, seed=0):
    """Bootstrap CI on the difference of means, a minus b."""
    rng = np.random.default_rng(seed)
    d = np.array([rng.choice(a, a.size, True).mean() - rng.choice(b, b.size, True).mean()
                  for _ in range(N_BOOT)])
    return float(a.mean() - b.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def depth_match(a_idx, b_idx, logdepth, rng):
    """Pair each a with its nearest unused b by log depth. Returns paired index arrays."""
    pa, pb, used = [], [], set()
    for i in sorted(a_idx, key=lambda i: logdepth[i]):
        best, bd = None, np.inf
        for j in b_idx:
            if j in used:
                continue
            d = abs(logdepth[i] - logdepth[j])
            if d < bd:
                best, bd = j, d
        if best is not None:
            used.add(best)
            pa.append(i); pb.append(best)
    return np.array(pa, int), np.array(pb, int)


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

    def cls(c):
        raw = (meta[col2bc[c]].get("Conserved_cluster") or "").strip()
        return "" if (not raw or raw.startswith("Doublets")) else raw.split("_")[0]

    def vmhclu(c):
        return meta[col2bc[c]]["VMH_cluster"].strip()

    in_lepr = np.array([vmhclu(c).startswith("Lepr") for c in cols])
    host = np.array([cls(c) for c in cols])
    print(f">> {len(cols)} VMH nuclei; {int(in_lepr.sum())} in Lepr clusters", flush=True)

    # Lepr detection per nucleus
    lrow = name2row["Lepr"]
    ent = V.pass_rows(PREFIX, {lrow}, TAG, "lepr")
    pos = set()
    for l in open(ent):
        r, c, v = l.split()
        if float(v) > 0:
            pos.add(int(c))
    lepr_pos = np.array([c in pos for c in cols])
    print(f">> Lepr detected: {int((lepr_pos & in_lepr).sum())} inside the Lepr clusters "
          f"({100*(lepr_pos & in_lepr).sum()/max(in_lepr.sum(),1):.1f}%), "
          f"{int((lepr_pos & ~in_lepr).sum())} outside "
          f"({100*(lepr_pos & ~in_lepr).sum()/max((~in_lepr).sum(),1):.1f}%)", flush=True)

    # the effector axis, pipeline's own first draw
    panel_rows, ctrl_rows = V.choose_controls(name2row, rowsum, len(cols), cfg.RANDOM_SEED)
    need = {r for rs_ in list(panel_rows.values()) + list(ctrl_rows.values()) for r in rs_}
    ent2 = V.pass_rows(PREFIX, need, TAG, "union")
    expr = defaultdict(dict)
    colset = set(cols)
    for l in open(ent2):
        r, c, v = l.split()
        c = int(c)
        if c in colset:
            expr[c][int(r)] = float(v)
    axis = V.module_scores(expr, cols, colsum, panel_rows, ctrl_rows)["HminusD"]
    depth = np.array([colsum[c] for c in cols], float)

    logdepth = np.log10(depth)

    # Does the axis track capture depth at all? If it does, every difference below is
    # suspect, because Lepr-positive nuclei are deeper by selection.
    rho_all, p_all = stats.spearmanr(axis, logdepth)
    rho_out, p_out = stats.spearmanr(axis[~in_lepr], logdepth[~in_lepr])
    print(f"\n>> axis vs log10 depth: rho {rho_all:+.3f} (p {p_all:.2g}) across all VMH, "
          f"{rho_out:+.3f} outside the Lepr clusters", flush=True)
    if abs(rho_all) > 0.1:
        print("   the axis tracks depth; read every difference below as partly depth",
              flush=True)

    lepr_axis = float(axis[in_lepr].mean())
    out = {"n_vmh": len(cols),
           "n_in_lepr_clusters": int(in_lepr.sum()),
           "n_lepr_positive_outside": int((lepr_pos & ~in_lepr).sum()),
           "mean_axis_lepr_clusters": round(lepr_axis, 4),
           "seed": cfg.RANDOM_SEED,
           "caveat": ("the grouping is a per-cell detection call, which section 2.1 argues "
                      "cannot define a population; treat as a check on what the gate "
                      "selects, not as a result about a population"),
           "axis_vs_log_depth_spearman_all_vmh": None,
           "axis_vs_log_depth_spearman_outside_lepr": None,
           "by_host_class": {}}

    out["axis_vs_log_depth_spearman_all_vmh"] = round(float(rho_all), 4)
    out["axis_vs_log_depth_spearman_outside_lepr"] = round(float(rho_out), 4)
    print(f"\n>> Lepr-cluster mean on the axis: {lepr_axis:+.4f}")
    print(f"\n{'host class':<9}{'n Lepr+':>8}{'n other':>9}{'Lepr+ mean':>12}"
          f"{'other mean':>12}{'difference':>12}{'95% CI':>22}{'depth ratio':>13}")

    for cl in sorted({c for c in host if c and c != "Lepr"}):
        sel = (host == cl) & ~in_lepr
        a = axis[sel & lepr_pos]
        b = axis[sel & ~lepr_pos]
        da, db = depth[sel & lepr_pos], depth[sel & ~lepr_pos]
        rec = {"n_lepr_positive": int(a.size), "n_other": int(b.size),
               "mean_lepr_positive": round(float(a.mean()), 4) if a.size else None,
               "mean_other": round(float(b.mean()), 4) if b.size else None,
               "median_depth_lepr_positive": round(float(np.median(da))) if da.size else None,
               "median_depth_other": round(float(np.median(db))) if db.size else None}
        if a.size >= MIN_PER_GROUP and b.size >= MIN_PER_GROUP:
            d, lo, hi = boot_diff(a, b, seed=cfg.RANDOM_SEED)
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            rec.update({"difference": round(d, 4), "ci_low": round(lo, 4),
                        "ci_high": round(hi, 4), "mannwhitney_p": float(p),
                        "depth_ratio": round(float(np.median(da) / np.median(db)), 2),
                        "ci_excludes_zero": bool(lo > 0 or hi < 0)})
            # the same comparison with depth equalised by nearest-neighbour pairing
            ia = np.where(sel & lepr_pos)[0]
            ib = np.where(sel & ~lepr_pos)[0]
            pa, pb = depth_match(ia, ib, logdepth, None)
            if pa.size >= MIN_PER_GROUP:
                dm, lom, him = boot_diff(axis[pa], axis[pb], seed=cfg.RANDOM_SEED)
                rec.update({
                    "n_pairs": int(pa.size),
                    "difference_depth_matched": round(dm, 4),
                    "ci_low_depth_matched": round(lom, 4),
                    "ci_high_depth_matched": round(him, 4),
                    "depth_ratio_after_matching": round(
                        float(np.median(depth[pa]) / np.median(depth[pb])), 2),
                    "ci_excludes_zero_depth_matched": bool(lom > 0 or him < 0)})
            print(f"{cl:<9}{a.size:>8}{b.size:>9}{a.mean():>+12.4f}{b.mean():>+12.4f}"
                  f"{d:>+12.4f}  [{lo:+.4f}, {hi:+.4f}]{rec['depth_ratio']:>13}")
            if "difference_depth_matched" in rec:
                print(f"{'  depth-matched':<9}{rec['n_pairs']:>8}{rec['n_pairs']:>9}"
                      f"{'':>12}{'':>12}{rec['difference_depth_matched']:>+12.4f}"
                      f"  [{rec['ci_low_depth_matched']:+.4f}, "
                      f"{rec['ci_high_depth_matched']:+.4f}]"
                      f"{rec['depth_ratio_after_matching']:>13}")
        else:
            rec["skipped"] = f"fewer than {MIN_PER_GROUP} nuclei in one group"
            print(f"{cl:<9}{a.size:>8}{b.size:>9}   skipped, too few")
        out["by_host_class"][cl] = rec

    # The first version counted host classes shifting "toward the Lepr-cluster mean".
    # That mean is near zero and mid-range among the classes, so the direction it measured
    # was an artifact of where each class sat. The sign of each difference, before and
    # after depth matching, is the thing to read.
    tested = [v for v in out["by_host_class"].values() if "difference" in v]
    matched = [v for v in tested if "difference_depth_matched" in v]
    if tested:
        out["summary"] = {
            "n_classes_tested": len(tested),
            "n_positive_raw": sum(1 for v in tested if v["difference"] > 0),
            "n_ci_excludes_zero_raw": sum(1 for v in tested if v["ci_excludes_zero"]),
            "n_positive_depth_matched": sum(
                1 for v in matched if v["difference_depth_matched"] > 0),
            "n_ci_excludes_zero_depth_matched": sum(
                1 for v in matched if v["ci_excludes_zero_depth_matched"]),
            "reading": ("Lepr-positive nuclei are deeper by selection, so a consistent "
                        "positive difference before matching and none after it means the "
                        "gate selects deep nuclei rather than a compositionally distinct "
                        "set; differences surviving matching are the interpretable ones")}
        sm = out["summary"]
        print(f"\n>> raw: {sm['n_positive_raw']} of {len(tested)} differences positive, "
              f"{sm['n_ci_excludes_zero_raw']} with a CI excluding zero")
        print(f">> depth-matched: {sm['n_positive_depth_matched']} of {len(matched)} "
              f"positive, {sm['n_ci_excludes_zero_depth_matched']} with a CI excluding zero")

    json.dump(out, open(os.path.join(cfg.RESULTS, "lepr_positive_outside.json"), "w"),
              indent=2)
    print("\nWritten to results/lepr_positive_outside.json")
    print("Caveat:", out["caveat"])


if __name__ == "__main__":
    main()
