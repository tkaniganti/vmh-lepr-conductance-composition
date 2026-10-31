#!/usr/bin/env python3
"""
43_axis_distribution.py - where does the per-nucleus effector axis sit relative to zero?

WHY. The Discussion argues that a graded distribution of effector balance centered near the
point where opposing conductances balance would produce mixed response signs without
requiring discrete molecular identities. That claim is currently qualitative. The per-nucleus
scores exist, so the distribution's position, spread and shape can be stated.

WHAT ZERO IS AND IS NOT. Zero on this axis means the two panels sit equally far above their
own matched control genes. It is a property of the control-gene draw, not a membrane
potential, and not a point where opposing conductances physiologically balance. A nucleus at
-0.01 is not closer to depolarizing than one at +0.01 in any electrical sense. So the
fraction below zero is reported here as a description of the distribution, NOT as an
assignment of cells to response classes.

THE REASON THAT CAUTION IS NOT HYPOTHETICAL. Section 2.8 reports that arcuate POMC and AgRP,
whose leptin responses are known and opposite, did not separate on this axis: the contrast
was not consistent in sign across control-gene draws, and discrimination came from Trpc5
alone. An axis that fails to distinguish two populations with known opposite responses cannot
be used to assign response sign to individual cells.

WHAT IS REPORTED. Mean, median, SD, the fraction below zero, skew, and how far the mean sits
from zero in units of the distribution's own SD. The last is the one that matters for the
Discussion claim: a distribution centered near zero relative to its own spread is what the
balance model predicts, and a distribution displaced by a large multiple of its SD is not.

Reads data/percell_scores.npz, which holds one control-gene draw for the 2,913 Lepr-cluster
nuclei. Single draw, so the figures carry no across-draw range; that is stated in the output
rather than hidden.

Output: results/axis_distribution.json
"""
import json
import os
import warnings

import numpy as np
from scipy import stats

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")


def main():
    z = np.load(V.dpath("percell_scores.npz"), allow_pickle=True)
    y = np.asarray(z["HminusD"], float)
    lab = np.asarray(z["labels"])
    print(f">> {y.size} Lepr-cluster nuclei, {len(set(lab))} clusters, "
          f"one control-gene draw", flush=True)

    sd = float(y.std())
    out = {
        "n": int(y.size),
        "n_clusters": int(len(set(lab))),
        "draws": 1,
        "caveat": ("zero is where the two panels sit equally far above their own matched "
                   "control genes; it is not a membrane potential and not a physiological "
                   "balance point, so the fraction below zero describes the distribution "
                   "and does not assign cells to response classes"),
        "mean": round(float(y.mean()), 4),
        "median": round(float(np.median(y)), 4),
        "sd": round(sd, 4),
        "min": round(float(y.min()), 4),
        "max": round(float(y.max()), 4),
        "frac_below_zero": round(float((y < 0).mean()), 4),
        # The manuscript quotes this as a PERCENTAGE, 47.0%, and
        # 44_manuscript_check.py compares values, so a declaration against the
        # fraction can never match the text. Same shape as the 85.4% and the
        # Trpc5 range: the quoted form needs its own key.
        "pct_below_zero": round(100 * float((y < 0).mean()), 1),
        "skew": round(float(stats.skew(y)), 4),
        "kurtosis_excess": round(float(stats.kurtosis(y)), 4),
        "mean_in_sd_units_from_zero": round(float(abs(y.mean()) / sd), 4),
    }
    for q in (1, 5, 25, 50, 75, 95, 99):
        out[f"p{q}"] = round(float(np.percentile(y, q)), 4)

    print(f"\n  mean    {out['mean']:>9.4f}")
    print(f"  median  {out['median']:>9.4f}")
    print(f"  SD      {out['sd']:>9.4f}")
    print(f"  range   {out['min']:>9.4f} to {out['max']:.4f}")
    print(f"  below zero                {100 * out['frac_below_zero']:.1f}% of nuclei")
    print(f"  mean distance from zero   {out['mean_in_sd_units_from_zero']:.3f} SD")
    print(f"  skew {out['skew']:+.3f}, excess kurtosis {out['kurtosis_excess']:+.3f}")
    print("\n  percentiles: " + ", ".join(f"p{q}={out[f'p{q}']:+.3f}"
                                          for q in (1, 5, 25, 50, 75, 95, 99)))

    # is the distribution centered near zero relative to its own spread?
    near = out["mean_in_sd_units_from_zero"] < 0.25
    out["centered_near_zero"] = bool(near)
    print(f"\n  centered near zero relative to its own spread: {near}")
    if near:
        print("  This supports the Discussion's graded-distribution claim in a way the")
        print("  current text states only qualitatively. It does NOT mean the cells split,")
        print("  and 2.4 reports the same distribution as unimodal.")
    else:
        print("  The distribution is displaced from zero by an appreciable fraction of its")
        print("  own spread. The Discussion's wording should be checked against that.")

    # per-cluster, for reference only
    out["by_cluster_frac_below_zero"] = {
        str(c): round(float((y[lab == c] < 0).mean()), 4) for c in sorted(set(lab))}
    print("\n  fraction below zero, by cluster (reference only, not an assignment):")
    for c in sorted(set(lab)):
        print(f"    {c:<10} {100 * out['by_cluster_frac_below_zero'][str(c)]:>5.1f}%  "
              f"n={int((lab == c).sum())}")

    json.dump(out, open(os.path.join(cfg.RESULTS, "axis_distribution.json"), "w"),
              indent=2)
    print("\nWritten to results/axis_distribution.json")


if __name__ == "__main__":
    main()
