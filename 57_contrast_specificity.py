#!/usr/bin/env python3
"""
57_contrast_specificity.py - are the five class contrasts specific to effector genes?

WHY. Section 2.2 rests on two things: a class-level variance component, and five direct
contrasts between the conserved Lepr class and the other five. On 2 Nov 2026 the component
was recomputed with corrected expected-mean-square coefficients and then failed a
specificity test: 8 of 100 detectability-matched random panels reached it
(27_correctness_checks.py, check_B_random_panels, pass=False). The contrasts have a
permutation null over class labels (10_conserved_classes.py, perm_p_seed1, 1e-4 to 0.018)
but have never been asked the specificity question. Leading 2.2 on them without this test
would put the weight on the one quantity whose specificity nobody checked, in the same
paragraph that reports the other failing it.

CRITERIA FIXED BEFORE THE RUN, so the result cannot be read whichever way suits the section:

  Single draw. The observed contrasts are recomputed on the SAME draw and control-gene
  seed as the panels, not taken from the 20-draw means in conserved_classes.json.
  Comparing a 20-draw mean against 100 single-draw panels would let the panels spread
  wider than the observed from draw noise alone.

  Per contrast, SIGNED. A panel counts as reaching the observed only if its difference is
  at least as large IN THE SAME DIRECTION. Lepr_vs_Esr1 is observed positive, the other
  four negative or small positive; an absolute-value criterion would let a panel pass by
  being extreme the wrong way.

  Joint criterion. The number of panels whose LARGEST of the five contrasts reaches the
  observed largest. Without it, one panel passing on one contrast reads as a panel
  passing on all five, and five per-contrast counts invite that reading.

CONTROLS, both gating:
  zero case     class labels permuted among nuclei: observed contrasts should collapse
                toward zero, and the per-contrast counts should rise toward 100
  positive case the effector panel scored against itself must reach the observed on every
                contrast by construction; anything else means the scoring path differs
                between the observed and panel arms

Output: results/contrast_specificity.json
"""
import json
import os
import warnings
from collections import defaultdict

import numpy as np

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "csp"
N_PANELS = 100
FOCUS = "Lepr"


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
    cls = np.array([c2c[c] for c in cols])
    classes = sorted(set(cls))
    others = [c for c in classes if c != FOCUS]
    print(f">> {len(cols)} nuclei, {len(classes)} conserved classes: {classes}",
          flush=True)

    seed = cfg.RANDOM_SEED
    panel = {p: [g for g in cfg.PANELS[p] if g in name2row] for p in ("HYPER", "DEPOL")}
    n_h, n_d = len(panel["HYPER"]), len(panel["DEPOL"])
    print(f">> panel sizes: HYPER {n_h}, DEPOL {n_d}; single draw at seed {seed}",
          flush=True)

    logd = np.log10([max(float(colsum[c]), 1.0) for c in cols])

    def contrasts_from(pr, cr, residualize=True):
        """Class-mean differences, optionally with log depth removed first.

        Script 10 computes these on raw HminusD: resid_on_depth is defined there
        but applied only in the variance-partition block, never to the contrasts.
        If capture depth differs by class, an arbitrary detectability-matched
        panel separates the classes through depth alone, which would explain why
        23 of 100 such panels reproduced the observed contrasts. Residualizing
        both arms asks whether anything survives that.
        """
        y = V.module_scores(expr_cache, cols, colsum, pr, cr)["HminusD"]
        if residualize:
            y = y - np.polyval(np.polyfit(logd, y, 1), logd)
        mu = {c: float(y[cls == c].mean()) for c in classes}
        return {f"{FOCUS}_vs_{o}": round(mu[FOCUS] - mu[o], 4) for o in others}

    # detectability pool for the random panels, matched to the effector panel
    need = {name2row[g] for gs in panel.values() for g in gs}
    pool_rows = [int(r) for r in rowsum if rowsum[r] > 0]
    rng = np.random.default_rng(seed)
    cand = list(rng.choice(pool_rows, min(4000, len(pool_rows)), replace=False))
    ent = V.pass_rows(PREFIX, need | set(cand), TAG, "csp")
    # module_scores indexes expr[column][row], nested dicts. The row-keyed array
    # layout that 53_effector_covariation.py uses is for its own scoring code and
    # is not interchangeable; passing it here raised
    # AttributeError: 'numpy.ndarray' object has no attribute 'get'.
    expr_cache = defaultdict(dict)
    det_count = defaultdict(int)
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset:
            r = int(r)
            fv = float(v)
            expr_cache[c][r] = fv
            if fv > 0:
                det_count[r] += 1
    n_cols = len(cols)
    detect = {r: det_count[r] / n_cols for r in need | set(cand)}

    # ---- observed, on this draw and this seed
    pr_obs, cr_obs = V.choose_controls(name2row, rowsum, len(cols), seed)
    obs_raw = contrasts_from(pr_obs, cr_obs, residualize=False)
    obs = contrasts_from(pr_obs, cr_obs)
    print(f">> observed, raw          : {obs_raw}", flush=True)
    print(f">> observed, depth-resid  : {obs}", flush=True)

    # ---- positive control: the effector panel against itself
    pos = contrasts_from(pr_obs, cr_obs)
    pos_ok = all(abs(pos[k] - obs[k]) < 1e-9 for k in obs)
    print(f">> positive control, effector panel vs itself: {pos_ok}", flush=True)

    # ---- random panels, matched to each effector gene's detection rate
    keys = list(obs)
    counts = {k: 0 for k in keys}
    panel_vals = []
    for i in range(N_PANELS):
        pr = {}
        for p_name, genes in panel.items():
            picked = []
            for g in genes:
                d = detect[name2row[g]]
                near = sorted(cand, key=lambda q: abs(detect[q] - d))[:200]
                picked.append(int(rng.choice(near)))
            pr[p_name] = picked
        _, cr = V.choose_controls(name2row, rowsum, len(cols), seed + 7919 + i)
        v = contrasts_from(pr, cr)
        panel_vals.append(v)
        for k in keys:
            if (obs[k] >= 0 and v[k] >= obs[k]) or (obs[k] < 0 and v[k] <= obs[k]):
                counts[k] += 1
        if (i + 1) % 25 == 0:
            print(f"   {i + 1}/{N_PANELS} panels", flush=True)

    # ---- joint criterion, fixed before the run
    obs_max = max(abs(obs[k]) for k in keys)
    joint = sum(1 for v in panel_vals if max(abs(v[k]) for k in keys) >= obs_max)

    # ---- zero control: class labels permuted
    rngz = np.random.default_rng(seed + 31)
    cls_perm = rngz.permutation(cls)
    yz = V.module_scores(expr_cache, cols, colsum, pr_obs, cr_obs)["HminusD"]
    muz = {c: float(yz[cls_perm == c].mean()) for c in classes}
    zero = {f"{FOCUS}_vs_{o}": round(muz[FOCUS] - muz[o], 4) for o in others}
    zero_ok = all(abs(zero[k]) < abs(obs[k]) for k in keys)

    out = {"n_nuclei": len(cols), "n_panels": N_PANELS, "seed": seed,
           "design": ("single control-gene draw; observed contrasts recomputed on the "
                      "same draw and seed as the panels, not the 20-draw means"),
           "criteria_fixed_before_run": (
               "per contrast, signed: a panel counts only if its difference is at least "
               "as large in the SAME direction as observed. Joint: panels whose largest "
               "of the five reaches the observed largest"),
           "observed_single_draw": obs,
           "observed_single_draw_raw": obs_raw,
           "panels_reaching_observed": counts,
           "panels_reaching_observed_pct": {k: 100.0 * counts[k] / N_PANELS
                                            for k in keys},
           "joint_largest_of_five": joint,
           "joint_pct": 100.0 * joint / N_PANELS,
           "observed_largest_abs": round(obs_max, 4),
           "control_zero_permuted_labels": zero,
           "control_zero_passes": bool(zero_ok),
           "control_positive_self_match": bool(pos_ok)}

    print(f"\n{'=' * 70}")
    if not (pos_ok and zero_ok):
        print("  A CONTROL FAILED. Nothing below is readable.")
        print(f"    positive (panel vs itself): {pos_ok}")
        print(f"    zero (permuted labels)    : {zero_ok}  {zero}")
    else:
        for k in keys:
            print(f"  {k:<22} observed {obs[k]:+.4f}   panels reaching it: "
                  f"{counts[k]}/{N_PANELS}")
        print(f"\n  joint, largest of five:   {joint}/{N_PANELS}")
        print(f"  (the component itself:    8/100, which it failed)")
        if joint <= 5 and max(counts.values()) <= 10:
            print("\n  The contrasts are specific to the effector panel in a way the")
            print("  class component is not. 2.2 can lead on them.")
        else:
            print("\n  Arbitrary panels reproduce the contrasts at a similar rate to the")
            print("  component. 2.2 cannot lead on them either, and the abstract's")
            print("  ordering needs revisiting.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "contrast_specificity.json"), "w"),
              indent=2)
    print("\nWritten to results/contrast_specificity.json")


if __name__ == "__main__":
    main()
