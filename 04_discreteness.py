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

EFFECT_SIZES = [0.25, 0.5, 1.0, 2.0]   # synthetic log-fold shifts in the effector panels
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
    for eff in EFFECT_SIZES:
        eps_rec, k_rec, dip_rec = [], [], []
        for _ in range(N_SIM):
            idx = rng.permutation(hd.size)
            half = idx[: hd.size // 2]
            sim = hd.copy()
            sim[half] += eff * sd
            fake = np.array(["A"] * hd.size, dtype=object)
            fake[half] = "B"
            Hs, _ = kruskal(sim[fake == "A"], sim[fake == "B"])
            eps_rec.append(V.eps2_from_kruskal(float(Hs), 2, sim.size))
            _, b = fit_mixture(sim, kmax=3, seed=int(rng.integers(1e6)))
            k_rec.append(b["k"])
            if HAVE_DIP:
                dip_rec.append(diptest.diptest(sim)[1])
        power.append({
            "effect_sd_units": eff,
            "eps2_median": round(float(np.median(eps_rec)), 4),
            "pct_sims_bic_picks_k_ge_2": round(100 * float(np.mean(np.array(k_rec) >= 2)), 1),
            "pct_sims_dip_rejects_unimodal": (
                round(100 * float(np.mean(np.array(dip_rec) < 0.05)), 1) if HAVE_DIP else None),
        })
        print(f"   shift {eff} SD -> eps2 {power[-1]['eps2_median']:.3f}, "
              f"BIC k>=2 in {power[-1]['pct_sims_bic_picks_k_ge_2']}% of sims")
    res["power_simulation"] = {"n_sim_per_level": N_SIM, "levels": power}

    res["interpretation"] = (
        "If BIC selects k = 1 and the dip test does not reject unimodality on the observed "
        "H-D score, while the simulation shows that a split of around X SD would have been "
        "detected, the graded claim is supported at the resolution of this measurement and "
        "should be stated with that bound. If BIC selects k >= 2, the graded framing must be "
        "revised rather than defended.")

    json.dump(res, open(os.path.join(C.RESULTS, "discreteness.json"), "w"), indent=2)
    print("\nWritten to", os.path.join(C.RESULTS, "discreteness.json"))


if __name__ == "__main__":
    main()
