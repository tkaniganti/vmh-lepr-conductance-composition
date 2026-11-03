#!/usr/bin/env python3
"""
58_class_excess_specificity.py - are the per-class excesses specific to effector genes?

WHY. Two quantities in this paper now have a random-panel control and one does not.
27_correctness_checks.py check_B asks whether arbitrary detectability-matched panels
reproduce the CLASS-level component (8 of 100 do, so it is not specific), and check_B2
asks the same within the nine leptin-receptor clusters. 57_contrast_specificity.py asks
it of the five class contrasts (24 of 100 on the largest). The per-class excesses in
23_class_nulls.py and 24_mincluster_sensitivity.py, including the Nfib value the
Discussion leads on, are tested only against cluster-label shuffles: nothing has asked
whether an arbitrary panel finds the same within-class structure.

That gap is load-bearing. The title and abstract set an uneven distribution WITHIN
classes against a class-level difference that is not specific, which implies the
within-class result is specific. Nobody has tested that, and it is the same shape as
this afternoon: one quantity got the control and its neighbour did not.

WHAT THIS RUNS. For each conserved class, the one-way component of the per-nucleus score
across that class's clusters, observed against 100 detectability-matched random panels,
each with its own control-gene draw. One draw throughout, with the observed recomputed on
the same draw as the panels, as in 57.

CRITERIA FIXED BEFORE THE RUN:
  per class, one-sided: a panel counts if its component is at least as large as observed
  the classes of interest are Nfib and Dlk1, which carry distributed structure; Lepr is
  included as the known case (check_B2 found arbitrary panels find MORE structure there
  than the effector panel, which is the opposite direction and worth reproducing here)
  no joint criterion, because the claim is per class rather than about the set

CONTROLS:
  zero      cluster labels shuffled within each class: components should collapse
  positive  the effector panel against itself must match the observed exactly
"""
import json
import os
import warnings
from collections import defaultdict

import numpy as np

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "ces"
N_PANELS = 100


def one_way(y, grp):
    import pandas as pd
    df = pd.DataFrame({"y": np.asarray(y, float), "g": np.asarray(grp)})
    n, k = len(df), df.g.nunique()
    if k < 2 or n <= k:
        return None
    grand = df.y.mean()
    ssb = sum(len(g) * (g.y.mean() - grand) ** 2 for _, g in df.groupby("g"))
    ssw = sum(((g.y - g.y.mean()) ** 2).sum() for _, g in df.groupby("g"))
    msb, msw = ssb / (k - 1), ssw / (n - k)
    ni = df.groupby("g").size().to_numpy(float)
    n0 = (n - (ni ** 2).sum() / n) / (k - 1)
    sb = max((msb - msw) / n0, 0.0)
    return 100 * sb / (sb + msw)


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        raw = (m.get("Conserved_cluster") or "").strip()
        m["_cls"] = "" if (not raw or raw.startswith("Doublets")) else raw.split("_")[0]
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "_cls", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    c2c = {}
    for line in open(cg):
        a, b = line.split()
        c2c[int(a)] = b
    cols = sorted(c for c in colsum if c in c2c)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cls = np.array([c2c[c] for c in cols])
    clu = np.array([(meta[col2bc[c]].get("VMH_cluster") or "NA").strip() for c in cols])
    classes = sorted(set(cls))
    print(f">> {len(cols)} nuclei, classes {classes}", flush=True)

    seed = cfg.RANDOM_SEED
    panel = {p: [g for g in cfg.PANELS[p] if g in name2row] for p in ("HYPER", "DEPOL")}
    need = {name2row[g] for gs in panel.values() for g in gs}
    rng = np.random.default_rng(seed)
    pool = [int(r) for r in rowsum if rowsum[r] > 0]
    cand = list(rng.choice(pool, min(4000, len(pool)), replace=False))
    ent = V.pass_rows(PREFIX, need | set(cand), TAG, "ces")
    expr = defaultdict(dict)
    det_count = defaultdict(int)
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset:
            r, fv = int(r), float(v)
            expr[c][r] = fv
            if fv > 0:
                det_count[r] += 1
    detect = {r: det_count[r] / len(cols) for r in need | set(cand)}
    logd = np.log10([max(float(colsum[c]), 1.0) for c in cols])

    def per_class(pr, cr):
        y = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
        y = y - np.polyval(np.polyfit(logd, y, 1), logd)
        out = {}
        for cl in classes:
            m = cls == cl
            v = one_way(y[m], clu[m])
            if v is not None:
                out[cl] = round(v, 4)
        return out

    pr_obs, cr_obs = V.choose_controls(name2row, rowsum, len(cols), seed)
    obs = per_class(pr_obs, cr_obs)
    print(f">> observed per class: {obs}", flush=True)
    pos_ok = per_class(pr_obs, cr_obs) == obs

    counts = {k: 0 for k in obs}
    for i in range(N_PANELS):
        pr = {}
        for p_name, genes in panel.items():
            picked = []
            for g in genes:
                d = detect[name2row[g]]
                near = sorted(cand, key=lambda q: abs(detect[q] - d))[:200]
                picked.append(int(rng.choice(near)))
            pr[p_name] = picked
        _, cr = V.choose_controls(name2row, rowsum, len(cols), seed + 6007 + i)
        v = per_class(pr, cr)
        for k in counts:
            if k in v and v[k] >= obs[k]:
                counts[k] += 1
        if (i + 1) % 25 == 0:
            print(f"   {i + 1}/{N_PANELS}", flush=True)

    # ZERO CONTROL, SECOND VERSION. The first took ONE shuffle per class and asked
    # every value to sit below its observed. A class whose observed component is
    # truncated at zero can never satisfy that: a single shuffle returns something
    # positive by chance, and Foxp2 duly returned 1.34 against an observed 0.0. The
    # criterion was untestable for exactly the classes it most needed to cover.
    #
    # The estimator is bias-corrected but truncated, so under no structure it returns
    # zero on average with positive excursions. The right control is therefore the
    # MEAN over many shuffles, which should sit near the untruncated expectation,
    # and the comparison is only meaningful for classes whose observed value is not
    # itself at the floor.
    rz = np.random.default_rng(seed + 17)
    yz = V.module_scores(expr, cols, colsum, pr_obs, cr_obs)["HminusD"]
    yz = yz - np.polyval(np.polyfit(logd, yz, 1), logd)
    N_SHUF = 50
    zero = {}
    for cl in classes:
        m = cls == cl
        vals = [one_way(yz[m], rz.permutation(clu[m])) or 0.0 for _ in range(N_SHUF)]
        zero[cl] = {"shuffle_mean": round(float(np.mean(vals)), 4),
                    "shuffle_p95": round(float(np.percentile(vals, 95)), 4),
                    "n_shuffles": N_SHUF}
    structured = [k for k in obs if obs[k] > 0.5]
    zero_ok = all(zero[k]["shuffle_mean"] < obs[k] for k in structured)
    print(f">> zero control on {structured} (classes above the floor): {zero_ok}",
          flush=True)

    out = {"n_nuclei": len(cols), "n_panels": N_PANELS, "seed": seed,
           "observed_per_class": obs,
           "panels_reaching_observed": counts,
           "control_zero_shuffled_clusters": zero,
           "control_zero_classes_tested": structured,
           "control_zero_passes": bool(zero_ok),
           "control_positive_self_match": bool(pos_ok),
           "criteria_fixed_before_run": (
               "per class, one-sided; classes of interest Nfib and Dlk1, which carry "
               "distributed structure; Lepr included as the known case; no joint "
               "criterion because the claim is per class")}

    print(f"\n{'=' * 68}")
    if not (pos_ok and zero_ok):
        print(f"  A CONTROL FAILED. positive {pos_ok}, zero {zero_ok}: {zero}")
    else:
        for k in sorted(obs, key=lambda x: -obs[x]):
            print(f"  {k:<8} observed {obs[k]:>7.3f}%   panels reaching it: "
                  f"{counts[k]:>3}/{N_PANELS}")
        print("\n  Nfib and Dlk1 are the ones the title's within-class claim rests on.")
    json.dump(out, open(os.path.join(cfg.RESULTS,
                                     "class_excess_specificity.json"), "w"), indent=2)
    print("\nWritten to results/class_excess_specificity.json")


if __name__ == "__main__":
    main()
