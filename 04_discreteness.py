#!/usr/bin/env python3
"""
04_discreteness.py - test "graded, not discrete" directly, and check that the test has power.

The Kruskal-Wallis across the nine deposited LepR clusters asks whether effector balance
tracks labels that were built without effector genes. A small epsilon-squared is the
expected result whether or not discrete effector subtypes exist, so it cannot support the
graded claim on its own. Two things are needed:

  1. An unsupervised test on the effector axis itself: Gaussian mixture with BIC over
     k = 1..4, and Hartigan's dip test for unimodality, on the per-cell H-D score.
  2. A power check: inject a synthetic two-subtype split of known size into the observed
     counts, re-score under the same dropout, and record what epsilon-squared and what
     BIC-selected k are recovered. This says how large a real split would have to be
     before this pipeline could see it.

Output: results/discreteness.json
"""
import os, json
import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import kruskal
import config as C
import vmhlib as V

try:
    import diptest
    HAVE_DIP = True
except Exception:
    HAVE_DIP = False

EFFECT_SIZES = [0.25, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0]   # injected split sizes, in units of the observed SD (see line ~87)
N_SIM = 25


def fit_mixture(x, kmax=4, seed=0):
    x = np.asarray(x, float).reshape(-1, 1)
    out = []
    for k in range(1, kmax + 1):
        gm = GaussianMixture(n_components=k, covariance_type="full",
                             random_state=seed, n_init=5).fit(x)
        out.append({"k": k, "bic": float(gm.bic(x)), "aic": float(gm.aic(x)),
                    "weights": [round(float(w), 3) for w in gm.weights_],
                    "means": [round(float(m[0]), 3) for m in gm.means_]})
    best = min(out, key=lambda r: r["bic"])
    return out, best


def main():
    z = np.load(V.dpath("percell_scores.npz"), allow_pickle=True)
    hd = np.asarray(z["HminusD"], float)
    labels = np.asarray(z["labels"])
    print(f">> {hd.size} LepR nuclei")

    # ------------------------------------------------------------ observed structure
    models, best = fit_mixture(hd)
    res = {"n": int(hd.size),
           "observed": {"mixture": models, "best_k_by_bic": best["k"],
                        "bic_gap_k1_minus_kbest": round(models[0]["bic"] - best["bic"], 1)}}
    if HAVE_DIP:
        d, p = diptest.diptest(hd)
        res["observed"]["dip"] = {"D": float(d), "p": float(p),
                                  "unimodal_not_rejected": bool(p > 0.05)}
    else:
        res["observed"]["dip"] = "diptest not installed; mixture/BIC only"

    groups = [hd[labels == k] for k in sorted(set(labels.tolist()))]
    Hst, p = kruskal(*groups)
    res["observed"]["kruskal_over_deposited_clusters"] = {
        "H": round(float(Hst), 2), "p": float(p),
        "eps2": round(V.eps2_from_kruskal(float(Hst), len(groups), hd.size), 4)}

    # ------------------------------------------------------------ power simulation
    # Split the observed cells at random into two halves and shift one half's score by a
    # known amount in units of the observed within-cluster SD. The observed score already
    # carries the real dropout, so the recovered statistics are the ones this pipeline
    # would report if a true two-subtype split of that size existed.
    sd = float(np.std(hd))
    rng = np.random.default_rng(C.RANDOM_SEED)
    power = []
    # ---------------------------------------------------------------- power grid
    # RETIRED 29 Oct 2026. This loop simulated an exactly balanced two-component split
    # (idx[:hd.size // 2]) and scored recovery by BIC, and section 2.4 quoted its 2.0 SD
    # figure as the bound on discreteness. Two problems, both in the manuscript's
    # reviewer pass:
    #
    #   The configuration was wrong. Sohn et al. report 14 of 80 depolarized and 14
    #   hyperpolarized, which is two minority groups either side of a dominant centre,
    #   not a balanced split. That configuration is not recovered at ANY separation,
    #   because opposite displacements broaden the marginal density without moving its
    #   centre and mixture selection fits a broadened density as one component.
    #
    #   The stored dip column was 0.0 at every level and was read as low power. It is
    #   not: a 50/50 equal-variance mixture is unimodal below 2 sigma, so the dip test
    #   was answering its own question correctly.
    #
    # 51_power_grid.py replaces this with six configurations and both criteria. Leaving
    # the old grid here would put two objects in results/ for one quantity, and a
    # verification run would reproduce and certify a number the manuscript no longer
    # makes. The observed-data results above are unaffected and still current.
    res["interpretation"] = (
        "Observed-data results only. The power grid that used to live here is retired; "
        "see 51_power_grid.py and results/power_grid.json. What BIC and the dip test "
        "report on the observed score does not on its own bound the subpopulation "
        "proposal, because the configuration that proposal describes is not recoverable "
        "by either criterion at any separation (section 2.4).")

    json.dump(res, open(os.path.join(C.RESULTS, "discreteness.json"), "w"), indent=2)
    print("\nWritten to", os.path.join(C.RESULTS, "discreteness.json"))


if __name__ == "__main__":
    main()
