#!/usr/bin/env python3
# SUPERSEDED IN PART, RETAINED FOR PROVENANCE. The 3.57% and 2.69% in the
# docstring below are pre-3-Nov-2026 values from the uncorrected estimator; the
# corrected figures are 4.73% and 1.63%. The diagnosis that part C draws from them
# is superseded: the apparent class-level leak was the expected-mean-square
# coefficients, not the design's imbalance. See 22_nested_partition.py.
"""
40_leak_diagnosis.py - three questions the leak result leaves open.

38_upper_level_recovery.py established that nested() reports class-level variance when none
is present, scaling with the cluster component: 0.400% at a true cluster share of 1.8%,
0.873% at 3.6%, 1.433% at 7.2%, with the zero-zero control clean at 0.025%. The observed
cluster component is 3.57%, so roughly 0.87 of the observed 2.69% class component is
spurious.

That is enough to know the manuscript is wrong. It is not enough to know what to write
instead. This answers three things:

  A. IS IT THE IMBALANCE? Rerun the same injection on a balanced design with the same totals
     (6 classes, equal nuclei, equal groups per class). If the leak vanishes, the cause is
     confirmed as unbalancedness and a correction is principled. If it survives, the formula
     is being applied wrongly and that is a different and larger problem.

  B. DOES IT REACH THE PER-CLASS EXCESSES IN 2.3? Those are cluster-level estimates computed
     inside one class, so the upward leak into the class term should not touch them. But
     each is reported as an excess over a null built by shuffling cluster labels within the
     class, which preserves group count and sizes, so the null may absorb the same artifact
     and the excess may be close to right. Checked here by injecting known cluster effects
     within a single class and seeing whether observed minus null recovers them.

  C. WHAT IS THE CORRECTED CLASS COMPONENT? Invert the leak curve at the observed cluster
     component to get a leak-subtracted estimate, with an interval, so 2.2 can state a
     number instead of withdrawing one.

Output: results/leak_diagnosis.json
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
N_SIM = 300
OBSERVED_CLUSTER_PCT = 3.5668
OBSERVED_CLASS_PCT = 2.6884


def load_nested():
    sp = importlib.util.spec_from_file_location("m22", "22_nested_partition.py")
    m = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(m)
    return m.nested


def simulate(nested, cls, clu, groups, classes, tc, tb, rng, n_sim=N_SIM):
    """Inject class share tc and cluster share tb, return recovered class percentages."""
    ve = 1.0 - tc / 100 - tb / 100
    out = []
    for _ in range(n_sim):
        a = ({c: rng.normal(0, np.sqrt(tc / 100)) for c in classes}
             if tc > 0 else {c: 0.0 for c in classes})
        b = ({g: rng.normal(0, np.sqrt(tb / 100)) for g in groups}
             if tb > 0 else {g: 0.0 for g in groups})
        y = (np.array([a[c] for c in cls])
             + np.array([b[f"{c}|{s}"] for c, s in zip(cls, clu)])
             + rng.normal(0, np.sqrt(ve), cls.size))
        out.append(nested(y, cls, clu)["pct_class"])
    return np.array(out)


def main():
    nested = load_nested()
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
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    out = {"n_sim": N_SIM, "observed_class_pct": OBSERVED_CLASS_PCT,
           "observed_cluster_pct": OBSERVED_CLUSTER_PCT}

    # ---------------------------------------------------------------- A: imbalance?
    print("=" * 66)
    print("A. is the leak caused by the imbalance?")
    print("=" * 66)
    sizes = [int((cls == c).sum()) for c in classes]
    gps = [len({g for g in groups if g.startswith(c + "|")}) for c in classes]
    print(f"   real design: nuclei per class {min(sizes)} to {max(sizes)}, "
          f"groups per class {min(gps)} to {max(gps)}")

    # balanced design with the same totals
    n_bal = (cls.size // (len(classes) * 11)) * (len(classes) * 11)
    per_class = n_bal // len(classes)
    per_group = per_class // 11
    bcls, bclu = [], []
    for c in classes:
        for g in range(11):
            bcls += [c] * per_group
            bclu += [f"g{g}"] * per_group
    bcls, bclu = np.array(bcls), np.array(bclu)
    bgroups = sorted({f"{a}|{b}" for a, b in zip(bcls, bclu)})
    print(f"   balanced design: {bcls.size} nuclei, {len(classes)} classes, "
          f"{per_group} nuclei in each of 11 groups per class")

    out["imbalance_test"] = {}
    print(f"\n   {'true clus %':>12}{'unbalanced':>13}{'balanced':>11}")
    for tb in (0.0, 3.6, 7.2):
        u = simulate(nested, cls, clu, groups, classes, 0.0, tb, rng).mean()
        b = simulate(nested, bcls, bclu, bgroups, classes, 0.0, tb, rng).mean()
        out["imbalance_test"][str(tb)] = {"unbalanced": round(float(u), 4),
                                          "balanced": round(float(b), 4)}
        print(f"   {tb:>12.1f}{u:>13.3f}{b:>11.3f}")
    ub = out["imbalance_test"]["3.6"]
    out["leak_is_imbalance"] = bool(ub["balanced"] < 0.2 < ub["unbalanced"])
    print(f"\n   leak explained by imbalance: {out['leak_is_imbalance']}")
    if not out["leak_is_imbalance"]:
        print("   NOTE: the leak survives balancing, so unbalancedness is not the cause.")
        print("   The estimator or its application needs re-deriving before any rewrite.")

    # ---------------------------------------------------------------- B: 2.3 excesses
    print()
    print("=" * 66)
    print("B. does the leak reach the per-class excesses in 2.3?")
    print("=" * 66)
    print("   Within one class the class term does not exist, so the upward leak cannot")
    print("   apply. What matters is whether observed minus shuffled-null recovers a known")
    print("   cluster effect. Tested on the leptin-receptor class at its own n and k.\n")

    def within_excess(y, g, rng_, n_null=200):
        def comp(lab):
            gt, subs = y.mean(), sorted(set(lab))
            ssb = sse = 0.0
            sizes_ = []
            for s in subs:
                ys = y[lab == s]
                ssb += ys.size * (ys.mean() - gt) ** 2
                sse += ((ys - ys.mean()) ** 2).sum()
                sizes_.append(ys.size)
            dfb, dfe = len(subs) - 1, y.size - len(subs)
            if dfb <= 0 or dfe <= 0:
                return 0.0
            msb, mse = ssb / dfb, sse / dfe
            n0 = (y.size - sum(v * v for v in sizes_) / y.size) / dfb
            s2 = max((msb - mse) / n0, 0.0)
            return 100 * s2 / (s2 + mse)
        obs = comp(g)
        lab = g.copy()
        nulls = []
        for _ in range(n_null):
            rng_.shuffle(lab)
            nulls.append(comp(lab))
        return obs, float(np.mean(nulls))

    m = cls == "Lepr"
    gl = clu[m]
    subs = sorted(set(gl))
    out["within_class_test"] = {}
    print(f"   {'injected %':>12}{'observed':>11}{'null':>9}{'excess':>10}{'recovers':>11}")
    for true_b in (0.0, 1.0, 3.0, 6.0):
        exc = []
        for _ in range(40):
            eff = ({s: rng.normal(0, np.sqrt(true_b / 100)) for s in subs}
                   if true_b > 0 else {s: 0.0 for s in subs})
            y = (np.array([eff[s] for s in gl])
                 + rng.normal(0, np.sqrt(1 - true_b / 100), gl.size))
            o, nl = within_excess(y, gl, rng)
            exc.append((o, nl, o - nl))
        o_m = float(np.mean([e[0] for e in exc]))
        n_m = float(np.mean([e[1] for e in exc]))
        e_m = float(np.mean([e[2] for e in exc]))
        ok = abs(e_m - true_b) < max(0.5, 0.25 * true_b)
        out["within_class_test"][str(true_b)] = {
            "observed": round(o_m, 4), "null": round(n_m, 4),
            "excess": round(e_m, 4), "recovers_truth": bool(ok)}
        print(f"   {true_b:>12.1f}{o_m:>11.3f}{n_m:>9.3f}{e_m:>10.3f}{str(ok):>11}")
    out["excesses_sound"] = bool(all(v["recovers_truth"]
                                     for v in out["within_class_test"].values()))
    print(f"\n   per-class excesses recover injected truth: {out['excesses_sound']}")
    print("   If true, 2.3 stands as written and only 2.2's class figure changes.")

    # ---------------------------------------------------------------- C: corrected value
    print()
    print("=" * 66)
    print("C. what should 2.2 report?")
    print("=" * 66)
    print("   Leak measured at the observed cluster component, then subtracted.\n")
    leak = simulate(nested, cls, clu, groups, classes, 0.0, OBSERVED_CLUSTER_PCT, rng)
    lm = float(leak.mean())
    llo, lhi = (float(v) for v in np.percentile(leak, [2.5, 97.5]))
    corrected = OBSERVED_CLASS_PCT - lm
    out["leak_at_observed_cluster"] = {"mean": round(lm, 4),
                                       "p2_5": round(llo, 4), "p97_5": round(lhi, 4)}
    out["corrected_class_pct"] = round(corrected, 4)
    out["corrected_range"] = [round(OBSERVED_CLASS_PCT - lhi, 4),
                              round(OBSERVED_CLASS_PCT - llo, 4)]
    print(f"   observed class component      {OBSERVED_CLASS_PCT:>8.4f}%")
    print(f"   leak at cluster {OBSERVED_CLUSTER_PCT:.2f}%          {lm:>8.4f}%  "
          f"(95% of simulations {llo:.3f} to {lhi:.3f})")
    print(f"   leak-subtracted estimate      {corrected:>8.4f}%  "
          f"(range {out['corrected_range'][0]:.3f} to {out['corrected_range'][1]:.3f})")
    print(f"\n   roughly {100 * lm / OBSERVED_CLASS_PCT:.0f}% of the reported class "
          f"component is leakage.")
    print("   The class level still carries structure, but less than stated, and the")
    print("   'floor' reading must be removed: the estimator inflates this term.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "leak_diagnosis.json"), "w"), indent=2)
    print("\nWritten to results/leak_diagnosis.json")


if __name__ == "__main__":
    main()
