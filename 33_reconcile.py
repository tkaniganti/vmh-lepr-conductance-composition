#!/usr/bin/env python3
"""
33_reconcile.py - print every number the Results quote, next to the key it came from.

A review of v49 found seven numerical items that do not reconcile: a percentage whose
denominator is unclear, one quantity reported with two values, a threshold series starting
from a figure different from the one first reported, leave-one-out values that all exceed
the full estimate, two classes sharing a suspiciously identical number, a correlation given
as 0.137 in one place and 0.143 in another, and a sentence whose stated maximum p
contradicts its own claim.

Each of those is checkable against the deposited results. This prints the stored value and
its full key path for every quantity the Results section quotes, so the manuscript can be
corrected against the files rather than against memory. It computes nothing and is not part
of the analysis: it reads results/ and reports.

Where a quantity exists in more than one form, draw-mean and single-draw for instance, both
are printed with their paths, because the point is to make the distinction visible rather
than to pick one.

Usage: run after the pipeline. Output is printed, not written.
"""
import os, json, sys
import config as cfg

R = cfg.RESULTS


def load(fn):
    p = os.path.join(R, fn)
    if not os.path.exists(p):
        print(f"  !! {fn} missing")
        return None
    return json.load(open(p))


def show(label, obj, path, note=""):
    """Walk a key path and print what is there."""
    cur = obj
    try:
        for k in path:
            cur = cur[k]
    except (KeyError, TypeError, IndexError):
        print(f"  {label:<46} KEY NOT FOUND  {'.'.join(map(str, path))}")
        return None
    if isinstance(cur, dict):
        print(f"  {label:<46} {'.'.join(map(str, path))}")
        for k, v in list(cur.items())[:8]:
            if not isinstance(v, (dict, list)):
                print(f"      {k:<34} {v}")
    else:
        print(f"  {label:<46} {cur}   <- {'.'.join(map(str, path))}"
              + (f"   [{note}]" if note else ""))
    return cur


def main():
    print("=" * 92)
    print("Every Results-quoted number against its source key")
    print("=" * 92)

    # ---------------------------------------------------------------- 2.1
    print("\n--- 2.1 detection and gating")
    d = load("detection.json")
    if d:
        show("VMH nuclei analysed", d, ["n_vmh_with_counts"])
        show("VMH labelled", d, ["n_vmh_labelled"])
        show("neurons", d, ["n_neurons"])
        show("Lepr-cluster nuclei", d, ["n_lepr_labelled"])
        show("detection, all VMH", d, ["vmh_detect"])
        show("detection, Lepr clusters", d, ["lepr_detect"])
        show("naive gate", d, ["naive_gate"])
    db = load("doublet_sensitivity.json")
    if db:
        show("doublet-flagged among Lepr", db, ["n_lepr_doublet_flagged"])
        show("gate, full set", db, ["with_doublets", "naive_gate_lepr"])
        show("gate, doublet-free", db, ["without_doublets", "naive_gate_lepr"])
        show("doublet profile", db, ["doublet_profile"])

    # ---------------------------------------------------------------- 2.2
    print("\n--- 2.2 effector axis and per-cell effect sizes")
    pc = load("percell_stats.json")
    if pc:
        for mod in ("HYPER", "DEPOL", "HminusD"):
            print(f"  {mod}:")
            for k in ("eps2_mean", "eps2_min", "eps2_max", "p_max", "label"):
                show(f"    {k}", pc, [mod, k],
                     note=("NOTE: a MAXIMUM p below 0.05 means every draw was "
                           "significant") if k == "p_max" else "")

    print("\n--- 2.2 conserved classes (the 75.5% / 1,857 question)")
    cc = load("conserved_classes.json")
    if cc:
        print("  top-level keys:", list(cc.keys()))
        for k in cc:
            if isinstance(cc[k], (int, float, str)):
                print(f"    {k:<42} {cc[k]}")
        for k in ("class_sizes", "membership", "lepr_class"):
            if k in cc:
                show(f"  {k}", cc, [k])

    # ---------------------------------------------------------------- nested
    print("\n--- 2.2 nested partition (the 2.69 vs 2.636 question)")
    np_ = load("nested_partition.json")
    if np_:
        for sub in ("nested_all_vmh",):
            if sub in np_:
                for comp in ("pct_class", "pct_cluster_in_class"):
                    show(f"{comp}, draw mean", np_, [sub, comp, "mean"],
                         note="this is the 20-draw mean")
                    show(f"{comp}, range", np_, [sub, comp, "min"])
                    show(f"{comp}, range", np_, [sub, comp, "max"])
        for k in np_:
            if k != "nested_all_vmh":
                show(k, np_, [k])

    # ---------------------------------------------------------------- per class
    print("\n--- 2.2 per-class excesses (6.29 vs 6.19, Esr1 LOO, Fezf1 0.38)")
    cn = load("class_nulls.json")
    if cn and "by_class" in cn:
        for cl, v in cn["by_class"].items():
            print(f"  {cl}:")
            for k in ("observed_pct", "null_pct", "excess_over_null",
                      "n_per_cluster", "k_clusters"):
                if k in v:
                    val = v[k]
                    if isinstance(val, dict):
                        print(f"      {k:<30} {dict(list(val.items())[:4])}")
                    else:
                        print(f"      {k:<30} {val}")
            loo = v.get("leave_one_cluster_out")
            if isinstance(loo, dict):
                lo = min(loo.values()); hi = max(loo.values())
                print(f"      leave-one-out                  {lo} to {hi}"
                      f"   (full excess {v.get('excess_over_null')})")
                if lo > v.get("excess_over_null", lo):
                    print("        !! every removal EXCEEDS the full value; "
                          "check whether this is real")
                for name, val in sorted(loo.items(), key=lambda t: t[1])[:3]:
                    print(f"        lowest: {name:<12} {val}")

    print("\n--- 2.2 threshold series (does it start at the reported excess?)")
    mc = load("mincluster_sensitivity.json")
    if mc:
        print("  top-level keys:", list(mc.keys())[:8])
        for cl in ("Nfib", "Esr1", "Lepr"):
            for k in mc:
                if cl.lower() in str(k).lower():
                    show(f"  {k}", mc, [k])

    # ---------------------------------------------------------------- 2.4, 2.5
    print("\n--- 2.4 discreteness and power")
    disc = load("discreteness.json")
    if disc:
        # These sat at the top level and the values live under "observed", so all
        # three printed KEY NOT FOUND. The power block below them printed fine, so
        # nothing drew attention to it. Found 29 Oct when that block was retired.
        show("BIC by k", disc, ["observed", "mixture"])
        show("best k", disc, ["observed", "best_k_by_bic"])
        show("dip", disc, ["observed", "dip"])
        if "power_simulation" in disc:
            print("  NOTE: this discreteness.json still carries the retired power grid.")
            print("  Rerun 04_discreteness.py to drop it. The current grid lives in")
            print("  results/power_grid.json; see 51_power_grid.py.")
        pg = load("power_grid.json")
        if pg:
            print("  power grid (51_power_grid.py), separation for 80% recovery:")
            for cname, row in pg.get("thresholds_80pct", {}).items():
                print(f"      {cname:<14} BIC {str(row.get('bic')):>6}  "
                      f"dip {str(row.get('dip')):>6}")

    print("\n--- 2.5 within-population scan (0.137 vs 0.143)")
    for fn, keys in (("gene_scan_axis.csv", None), ("outlier_sensitivity.json",
                                                    ["scan", "as_reported"])):
        p = os.path.join(R, fn)
        if not os.path.exists(p):
            print(f"  !! {fn} missing")
            continue
        if fn.endswith(".csv"):
            import csv
            rows = list(csv.DictReader(open(p)))
            key = next((k for k in rows[0] if "rho" in k.lower()
                        or "corr" in k.lower()), None)
            if key:
                rows.sort(key=lambda r: -abs(float(r[key] or 0)))
                print(f"  top 3 in {fn} by |{key}|:")
                for r in rows[:3]:
                    g = r.get("gene") or list(r.values())[0]
                    print(f"      {g:<14} {key} {r[key]}   q {r.get('q', '')}")
        else:
            o = load(fn)
            if o:
                show("outlier-check scan", o, keys)

    print("\n" + "=" * 92)
    print("Compare each value above against the manuscript. Where two values exist for")
    print("one quantity, the manuscript must say which it is quoting.")
    print("=" * 92)


if __name__ == "__main__":
    main()
