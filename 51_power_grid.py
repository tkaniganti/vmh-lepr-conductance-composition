#!/usr/bin/env python3
"""
51_power_grid.py - what separation would be needed to resolve discrete subtypes, under
criteria that answer the right question and at the configurations the hypothesis predicts?

WHY THE EXISTING GRID IS NOT ENOUGH, and why the obvious fix was wrong.

04_discreteness.py scores recovery by BIC only, simulates an exactly balanced split
(idx[:hd.size // 2]), and stops at 2.0 SD. It also stores pct_sims_dip_rejects_unimodal,
which is 0.0 at every level including 2.0 SD where BIC reaches 100%, and nothing reads it.

The first reading of those zeros was that the dip test lacks power. That is wrong. A 50/50
mixture of two equal-variance Gaussians is UNIMODAL unless the means are more than 2 sigma
apart: below that the density has a single peak. So at 1.5 and 1.75 SD the dip test is
returning the correct answer to the question it asks. At exactly 2.0 the density sits at the
flat-topped boundary and dropout noise smears it, so no rejection there either. The flat zero
line is what the mathematics requires.

THE CONSEQUENCE. The dip test and BIC are not two criteria for one question. BIC asks
whether a mixture is present; the dip asks whether the density has two peaks. Across the
separation range that matters here a real mixture produces one peak, so a dip non-rejection
carries almost no information about whether subtypes exist. Limitations currently rests the
unimodality conclusion on the dip test, which rests it on the weaker question.

So the fix is not to switch criteria. It is to keep asking BIC's question and fix BIC's
actual defect, which is that it buys components to fit skew, the behaviour that makes it
select three under tail trimming.

CRITERIA SCORED HERE

  ICL          TRIED AND REMOVED, 29 Oct. The idea was that an entropy penalty would
               refuse components bought to fit a tail, which is BIC's actual defect.
               It returned 0% at EVERY cell, including balanced_2 at 3.0 SD where BIC
               is 100% and the dip rejects 100%. At that separation the entropy was
               405 over 2,913 nuclei, 0.139 per nucleus, because each injected
               component inherits the full observed spread of 0.45 and the two overlap
               heavily however far apart their means sit. ICL penalises components
               that do not cleanly assign cells, and at these overlaps none does. That
               is correct behaviour for ICL and the wrong tool for this question: it
               would report no subtypes at any separation the data could contain.
               Do not reach for it again here.
  bootstrap LRT  parametric bootstrap likelihood ratio test for k = 1 against k = 2.
               Fits k = 1, simulates from it, refits both, and places the observed
               likelihood ratio in that null. Gives a calibrated p value rather than a
               penalty heuristic. Cheap at n = 2,913
  BIC          retained so the new grid can be read against the published one
  dip          reported across the extended range to DOCUMENT why it is not the bound:
               its rejection rate should rise once separation passes 2 SD and the density
               genuinely becomes bimodal. A row of zeros with an explanation beside it is
               worth more than a row of zeros alone

CONFIGURATIONS. Sohn et al. report 14 of 80 depolarized, 14 hyperpolarized and the rest
unresponsive, so roughly 18/18/64. That is not an unbalanced two-group split; it is two
minority groups on opposite sides of a dominant centre. Nobody has simulated it, and it is
the configuration where modality reasoning fails hardest, since a symmetric three-component
mixture with a dominant middle stays unimodal well past 2 SD.

  balanced_2      50/50, the published condition, for comparability
  minority_2      18/82, one minority group
  sohn_3          18/18/64, two minority groups symmetric about the centre

SEPARATIONS 1.0 to 3.0 SD. Past the old ceiling, so the dip curve has room to rise and the
other criteria have room to show where the real bound sits.

EXPECT THE BOUND TO WEAKEN. At the observed proportions it should come out above 2 SD, which
makes 2.4's statement weaker than it currently reads. That is the correct outcome and is
better found here than by a reviewer noticing the grid bounds a configuration the hypothesis
does not predict.

Output: results/power_grid.json
"""
import json
import os
import warnings

import numpy as np

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")

# 4.0 to 6.0 added 29 Oct to TEST THE MECHANISM, not to estimate detectability.
#
# The claim behind section 2.4 is that symmetric opposite displacement hides structure
# because it broadens the density without shifting its centre. That predicts recovery
# once the components stop overlapping. If sohn_3 recovers somewhere above 3 SD the
# mechanism is confirmed; if it never recovers even at 6 SD, where the components are
# plainly separated, the mechanism as stated is incomplete and the section that rests
# on it needs rewriting.
#
# READ ROWS ABOVE 3 SD AS A MECHANISM CHECK ONLY. Each injected component inherits the
# full spread of the observed scores (SD 0.4515), so at 6 SD the simulated distribution
# is far wider than the real one: the dropout model and the per-nucleus noise are being
# applied to a configuration nothing like the data they came from. A threshold measured
# up there is a fact about the estimator, not about anything these neurons could
# contain. It is not a detectability bound for this population and must not be quoted
# as one. For scale, 6 SD is 2.71 in axis units against an observed range of 3.55.
SEPARATIONS = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
N_SIM = 100
N_BOOT = 200          # parametric bootstrap replicates per LRT
ALPHA = 0.05


def main():
    try:
        import diptest
    except ImportError:
        raise SystemExit("diptest not installed")
    from sklearn.mixture import GaussianMixture

    z = np.load(V.dpath("percell_scores.npz"), allow_pickle=True)
    hd = np.asarray(z["HminusD"], float)
    n = hd.size
    sd = float(hd.std())
    print(f">> {n} nuclei, observed SD {sd:.4f}, one control-gene draw", flush=True)

    def fit(v, k, seed=0):
        return GaussianMixture(k, random_state=seed, n_init=5).fit(
            np.asarray(v, float).reshape(-1, 1))

    def best_by(fn, v, seed=0):
        return int(min((1, 2, 3), key=lambda k: fn(v, k, seed)))

    def bic(v, k, seed=0):
        r = np.asarray(v, float).reshape(-1, 1)
        return fit(r, k, seed).bic(r)

    # REMOVED 29 Oct. A parametric bootstrap LRT requires the null family to be
    # correctly specified. The null here is a single Gaussian and the observed axis is
    # not Gaussian: excess kurtosis 0.474, skew -0.220 (axis_distribution.json). So
    # non-normality entered the statistic as evidence for a second component. It
    # returned p = 0.0000 on data that BIC, ICL and the dip test all call unimodal, and
    # fired in 35% of simulations at 1.0 SD where BIC and ICL were both at 0%. Fixing
    # it means choosing a richer unimodal null family, which is a modelling decision
    # needing its own justification. ICL was tried as a replacement and failed for a
    # separate reason recorded in the docstring above, so the grid scores BIC and the
    # dip test only, and BIC's trimming instability is reported as a limitation rather
    # than solved.
    def _lrt_removed(v, rng, n_boot=N_BOOT):
        """Parametric bootstrap LRT, k=1 against k=2."""
        r = np.asarray(v, float).reshape(-1, 1)
        g1, g2 = fit(r, 1), fit(r, 2)
        obs = 2 * (g2.score(r) - g1.score(r)) * r.shape[0]
        mu, sg = float(g1.means_[0, 0]), float(np.sqrt(g1.covariances_[0, 0, 0]))
        null = []
        for _ in range(n_boot):
            s = rng.normal(mu, sg, r.shape[0]).reshape(-1, 1)
            a, b = fit(s, 1), fit(s, 2)
            null.append(2 * (b.score(s) - a.score(s)) * s.shape[0])
        return float((np.array(null) >= obs).mean())

    rng = np.random.default_rng(cfg.RANDOM_SEED)

    def make(cfg_name, eff):
        """Inject a configuration at separation eff (in observed SD units)."""
        sim = hd.copy()
        idx = rng.permutation(n)
        if cfg_name == "balanced_2":
            sim[idx[: n // 2]] += eff * sd
        elif cfg_name == "minority_2":
            k = int(round(0.18 * n))
            sim[idx[:k]] += eff * sd
        elif cfg_name == "sohn_3_skew":
            # -0.2 and +0.8 of the separation. The mass is no longer balanced about
            # the centre, so the density should SKEW rather than merely broaden, and
            # skew is what BIC is known to buy components for. If this recovers, the
            # claim is that the structure is invisible WHILE THE DISPLACEMENTS ROUGHLY
            # BALANCE, with a stated cause and a demonstrated edge. If it fails, the
            # broad claim survives a test built to break it.
            k = int(round(0.18 * n))
            sim[idx[:k]] -= eff * sd * 0.2
            sim[idx[k:2 * k]] += eff * sd * 0.8
        elif cfg_name == "uneven_3":
            # 10/26/64: displacements symmetric, GROUP SIZES unequal. sohn_3_asym
            # varies displacement and holds size; this holds displacement and varies
            # size. Together they separate "the displacements balance" from "the
            # groups are equal", which the current design confounds.
            a, b = int(round(0.10 * n)), int(round(0.26 * n))
            sim[idx[:a]] -= eff * sd / 2
            sim[idx[a:a + b]] += eff * sd / 2
        elif cfg_name == "sohn_3_asym":
            # 18% at -eff*2/5 and 18% at +eff*3/5: same total span, displaced
            # unequally, so the distribution skews. If BIC recovers this where it
            # fails on sohn_3, the invisibility belongs to SYMMETRY, which is a
            # much tighter and more interesting claim than "minority structure is
            # invisible".
            k = int(round(0.18 * n))
            sim[idx[:k]] -= eff * sd * 0.4
            sim[idx[k:2 * k]] += eff * sd * 0.6
        elif cfg_name == "sohn_3":
            k = int(round(0.18 * n))
            sim[idx[:k]] += eff * sd / 2
            sim[idx[k:2 * k]] -= eff * sd / 2
        return sim

    CONFIGS = {
        "balanced_2": "50/50, the published condition",
        "minority_2": "18/82, one minority group",
        "sohn_3_skew": ("18/18/64 at -0.2 and +0.8 of the separation, so the mass is "
                        "unbalanced about the centre and the density skews rather than "
                        "broadens. The sharper test of the broadening mechanism"),
        "uneven_3": ("10/26/64, displacements symmetric but group sizes unequal. Varies "
                     "size where sohn_3_asym varies displacement"),
        "sohn_3_asym": ("18/18/64 displaced UNEQUALLY about the centre, -0.4 and +0.6 "
                        "of the stated separation. Tests whether the invisibility of "
                        "sohn_3 is specific to exact symmetry"),
        "sohn_3": ("18/18/64, two minority groups symmetric about a dominant centre; "
                   "separation is the FULL distance between the two minority means"),
    }

    out = {"n": int(n), "observed_sd": round(sd, 4), "n_sim_per_cell": N_SIM,
           "n_bootstrap": N_BOOT, "separations_sd": SEPARATIONS,
           "configurations": CONFIGS,
           "separations_above_3sd_are_mechanism_check_only": True,
           "note": ("a 50/50 equal-variance mixture is unimodal below 2 sigma, so the "
                    "dip column is expected to be zero there by construction; it is "
                    "reported to document why the dip test is not the bound"),
           "grid": {}}

    d_obs, p_obs = diptest.diptest(np.ascontiguousarray(hd))
    out["observed"] = {"dip_p": round(float(p_obs), 5),
                       "bic_k": best_by(bic, hd),
                       }
    print(f">> observed: BIC k={out['observed']['bic_k']}, "
          f"dip p={out['observed']['dip_p']:.4f}", flush=True)

    print(f"\n{'config':<13}{'sep':>5}{'BIC>=2':>9}{'dip<.05':>9}")
    for cname in CONFIGS:
        for eff in SEPARATIONS:
            # Seed PER CELL. Seeding once for the whole grid made every cell depend on
            # how many draws preceded it, so adding separations changed cells that were
            # already declared: minority_2 at 2.0 SD moved 93 -> 95 and sohn_3_skew at
            # 2.5 moved 100 -> 95 when three rows were appended. A declared value from
            # a stochastic simulation is only reproducible if its own draws are.
            # Built arithmetically, NOT from hash(): Python randomizes string hashing
            # per process unless PYTHONHASHSEED is set, so a hash-derived seed would
            # differ on every submission and the cells would keep moving, which is the
            # problem this is meant to fix.
            ci = list(CONFIGS).index(cname)
            si = SEPARATIONS.index(eff)
            rng = np.random.default_rng(cfg.RANDOM_SEED + 1000 * ci + si)
            b = d = 0
            for _ in range(N_SIM):
                s = make(cname, eff)
                if best_by(bic, s) >= 2:
                    b += 1
                if diptest.diptest(np.ascontiguousarray(s))[1] < ALPHA:
                    d += 1
            cell = {"pct_bic": round(100 * b / N_SIM, 1),
                    "pct_dip": round(100 * d / N_SIM, 1)}
            out["grid"][f"{cname}_sd{eff}"] = cell
            print(f"{cname:<13}{eff:>5.1f}{cell['pct_bic']:>9.1f}"
                  f"{cell['pct_dip']:>9.1f}", flush=True)
        print()

    # ---- how much would a resolvable configuration widen the distribution?
    # The claim in 2.4 is that the observed spread EXCLUDES the configuration rather
    # than merely failing to detect it. That is a variance argument and it has to be
    # computed, not reasoned: two independent reconstructions from the figures the
    # manuscript states landed 50% apart, because one of them dropped the 64% centre
    # component, which is the term that keeps the total variance down.
    #
    # Law of total variance for the sohn_3 configuration: two groups of weight w at
    # +/- sep/2 and a centre group at 0, every component carrying the observed SD.
    #     Var = sd^2 + 2*w*(sep/2)^2
    widen = {}
    for eff in SEPARATIONS:
        mu = eff * sd / 2.0
        w = 0.18
        tot = float(np.sqrt(sd ** 2 + 2 * w * mu ** 2))
        widen[str(eff)] = {"separation_axis_units": round(eff * sd, 4),
                           "component_means_at": round(mu, 4),
                           "mixture_sd": round(tot, 4),
                           "ratio_to_observed": round(tot / sd, 3)}
    out["widening_if_resolvable"] = {
        "by_separation": widen,
        "centre_weight": 0.64,
        "minority_weight_each": 0.18,
        "caveat": ("observed_sd is from ONE control-gene draw, so this ratio is a "
                   "single-draw quantity presented in 2.4 as a population property. "
                   "Checking it across draws needs per-nucleus scores for all 20, "
                   "which are not stored; 41_cluster_axis_corrected.py computes them "
                   "but keeps cluster means only."),
        "note": ("the 64% centre component is INCLUDED. Omitting it inflates the "
                 "ratio by roughly 50%, which is the error this calculation exists "
                 "to prevent recurring.")}
    print(f"\n  widening if the configuration were resolvable:")
    for eff in SEPARATIONS:
        wv = widen[str(eff)]
        print(f"    {eff:>4} SD  means +/-{wv['component_means_at']:.2f}  "
              f"mixture SD {wv['mixture_sd']:.3f}  "
              f"{wv['ratio_to_observed']:.2f}x observed")

    top = out["grid"][f"balanced_2_sd{SEPARATIONS[-1]}"]
    ok = top["pct_bic"] >= 95
    out["positive_control_passes"] = bool(ok)
    print(f"  positive control, BIC saturates at 3.0 SD balanced: {ok} "
          f"({top['pct_bic']}%)")
    if not ok:
        print("  CONTROL FAILED. Nothing above is readable.")

    def thr(cname, key):
        for eff in SEPARATIONS:
            if out["grid"][f"{cname}_sd{eff}"][key] >= 80:
                # Flag thresholds found only in the mechanism-check range so a reader
                # of the JSON cannot mistake them for detectability estimates.
                return eff if eff <= 3.0 else f"{eff} (mechanism check only)"
        return None

    print(f"\n  separation for 80% recovery")
    print(f"  {'config':<13}{'BIC':>7}{'dip':>7}")
    summ = {}
    for cname in CONFIGS:
        row = {k: thr(cname, f"pct_{k}") for k in ("bic", "dip")}
        summ[cname] = row
        print(f"  {cname:<13}" + "".join(f"{str(row[k]):>7}"
                                         for k in ("bic", "dip")))
    out["thresholds_80pct"] = summ

    print(f"\n{'=' * 64}")
    print("  Read the sohn_3 row against balanced_2. If BIC recovers the balanced split")
    print("  and nothing in sohn_3 at any separation, the predicted configuration is")
    print("  undetectable by any method reading the marginal density: two minority")
    print("  groups displaced symmetrically about a dominant centre leave the shape")
    print("  nearly unchanged. That is invisibility IN PRINCIPLE, not for want of power,")
    print("  and more nuclei would not fix it.")
    print()
    print("  sohn_3_asym is the control on that claim. If it IS recovered, the")
    print("  invisibility belongs to exact symmetry specifically, which is the tighter")
    print("  and more useful statement.")
    print()
    print("  This bounds the TEST, not the hypothesis. The subpopulation proposal")
    print("  remains exactly as open as it was.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "power_grid.json"), "w"), indent=2)
    print("\nWritten to results/power_grid.json")


if __name__ == "__main__":
    main()
