#!/usr/bin/env python3
"""
38_upper_level_recovery.py - does the nested estimator report class-level variance that is
not there?

ORIGIN. Section 2.2 claims simulation shows the estimator is conservative at the upper
level, recovering about 1% where 2.9% is present at six levels, and concludes the class
figure is a floor. That claim had no producer: script 22 injects CLUSTER-level effects into
the focus set, which is the Figure 3C bound, and nothing ever injected a class-level effect.
A first attempt at this simulation failed its own zero control, returning 0.888% class
variance when none was injected, so it is rewritten here to call script 22's own nested()
instead of a reimplementation and to vary both levels independently.

WHY THE ZERO CONTROL MIGHT FAIL FOR A REAL REASON. sigma2_A = (MSA - MSB(A)) / n0_A is the
standard unbalanced estimator, and it cancels cleanly only when the design is balanced. Here
the six classes hold 9 to 16 clusters each and 461 to 1,858 nuclei, so the coefficient on the
cluster variance differs between MSA and MSB(A) and some cluster variance can land in the
class term. If that is what the first run saw, the direction is the opposite of the
manuscript's claim: the class component would be inflated, not conservative, and 2.69% and
every per-class excess would be affected.

THE DECISIVE ROWS are the two with no class effect injected. With no cluster effect either,
the estimator must return zero. With a cluster effect and still no class effect, anything
appreciable is leakage.

The permutation control already in 2.2 cannot see this, because shuffling both label sets at
once destroys both levels together and never creates the case of real cluster structure with
no class structure.

Output: results/upper_level_recovery.json
"""
import importlib.util
import json
import os
import warnings

import numpy as np

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX = "GSE172204"
N_SIM = 200
CLASS_GRID = [0.0, 1.0, 2.9, 5.0]      # true class share, percent
CLUSTER_GRID = [0.0, 1.8, 3.6, 7.2]    # true cluster-within-class share, percent


def load_nested():
    """Script 22's own nested(), so this measures that estimator and not a copy of it."""
    sp = importlib.util.spec_from_file_location("m22", "22_nested_partition.py")
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m.nested


def main():
    nested = load_nested()
    print(">> using nested() imported from 22_nested_partition.py", flush=True)

    meta = V.load_metadata(PREFIX)
    cls, clu = [], []
    for m in meta.values():
        vc = (m.get("VMH_cluster") or "").strip()
        cc = (m.get("Conserved_cluster") or "").strip()
        if not vc or not cc or cc.startswith("Doublets"):
            continue
        cls.append(cc.split("_")[0])
        clu.append(vc)
    cls, clu = np.array(cls), np.array(clu)
    classes = sorted(set(cls))
    groups = sorted({f"{a}|{b}" for a, b in zip(cls, clu)})
    print(f">> {cls.size} nuclei, {len(classes)} classes, "
          f"{len(groups)} class-by-cluster groups", flush=True)
    for c in classes:
        print(f"   {c:<7} {int((cls == c).sum()):>5} nuclei, "
              f"{len({g for g in groups if g.startswith(c + '|')}):>3} groups", flush=True)

    rng = np.random.default_rng(cfg.RANDOM_SEED)
    out = {"n_nuclei": int(cls.size), "n_classes": len(classes),
           "n_class_by_cluster": len(groups), "n_sim_per_level": N_SIM,
           "estimator": "nested() imported from 22_nested_partition.py",
           "design": ("synthetic scores on the real class-by-cluster structure; class and "
                      "cluster effects of known size injected independently"),
           "grid": {}}

    print(f"\n{'true class':>11}{'true clus':>11}{'class out':>12}{'clus out':>11}"
          f"{'class 95% range':>22}")
    for tc in CLASS_GRID:
        for tb in CLUSTER_GRID:
            ve = 1.0 - tc / 100 - tb / 100
            if ve <= 0:
                continue
            ra, rb = [], []
            for _ in range(N_SIM):
                a = ({c: rng.normal(0, np.sqrt(tc / 100)) for c in classes}
                     if tc > 0 else {c: 0.0 for c in classes})
                b = ({g: rng.normal(0, np.sqrt(tb / 100)) for g in groups}
                     if tb > 0 else {g: 0.0 for g in groups})
                y = (np.array([a[c] for c in cls])
                     + np.array([b[f"{c}|{s}"] for c, s in zip(cls, clu)])
                     + rng.normal(0, np.sqrt(ve), cls.size))
                r = nested(y, cls, clu)
                ra.append(r["pct_class"])
                rb.append(r["pct_cluster_in_class"])
            lo, hi = np.percentile(ra, [2.5, 97.5])
            out["grid"][f"class{tc}_cluster{tb}"] = {
                "true_class_pct": tc, "true_cluster_pct": tb,
                "recovered_class_mean": round(float(np.mean(ra)), 4),
                "recovered_class_p2_5": round(float(lo), 4),
                "recovered_class_p97_5": round(float(hi), 4),
                "recovered_cluster_mean": round(float(np.mean(rb)), 4),
                "n_class_at_zero": int(sum(1 for v in ra if v == 0.0))}
            print(f"{tc:>11.1f}{tb:>11.1f}{np.mean(ra):>12.3f}{np.mean(rb):>11.3f}"
                  f"  [{lo:>6.3f}, {hi:>6.3f}]")

    g = out["grid"]
    null_null = g["class0.0_cluster0.0"]["recovered_class_mean"]
    null_clus = g["class0.0_cluster3.6"]["recovered_class_mean"]
    out["control_no_effects"] = null_null
    out["leak_with_cluster_only"] = null_clus
    out["control_passes"] = bool(null_null < 0.1)
    out["leak_detected"] = bool(null_clus > 0.3)

    print(f"\n  no effects injected:              class returns {null_null:.3f}%  "
          f"{'OK' if out['control_passes'] else 'FAIL, the simulation has a bug'}")
    print(f"  cluster 3.6%, no class injected:  class returns {null_clus:.3f}%")
    if out["control_passes"] and out["leak_detected"]:
        print("\n  LEAKAGE. With no class effect present the estimator still reports one,")
        print("  because the unbalanced design leaves cluster variance in the class term.")
        print("  The class component is inflated, not conservative. The manuscript's")
        print("  'floor' reading is backwards and 2.69% needs a correction, as do the")
        print("  per-class excesses in 2.3.")
    elif out["control_passes"]:
        print("\n  No leakage. Compare the class0.0 rows against the class2.9 rows to see")
        print("  whether recovery is conservative, and rewrite 2.2 to whatever holds.")
    else:
        print("\n  Fix the simulation before reading anything else here.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "upper_level_recovery.json"), "w"),
              indent=2)
    print("\nWritten to results/upper_level_recovery.json")


if __name__ == "__main__":
    main()
