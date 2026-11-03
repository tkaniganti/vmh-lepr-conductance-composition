#!/usr/bin/env python3
"""
52_axis_reliability.py - how much signal does the per-nucleus effector axis carry?

WHY THIS IS THE GATING TEST. Sections 2.3, 2.4 and 2.5 all report absences measured on the
per-nucleus balance axis: no within-class variance component above chance, a single mixture
component, no transcript above about 2% of the variance. Every one of those is a statement
about a quantity whose reliability has never been estimated. If a nucleus's axis value is
mostly sampling noise, all three absences are guaranteed by the measurement and none of them
is about biology.

The existing positive control cannot detect this. It injects a shift INTO THE MODULE SCORE,
in the score's own units, and recovers it. That shows the estimator responds to a shift in
the score. It says nothing about whether a real difference in channel composition produces a
shift of that size, or whether a single nucleus's score is stable at all.

TWO SPLITS, measuring different things.

  PANEL SPLIT   Divide each panel into halves by alternating detection rank, build two
                independent axes, correlate them across nuclei. This asks whether the genes
                agree with each other. Spearman-Brown then estimates the reliability of the
                full panel from the half-panel correlation.

  COUNT SPLIT   Split each nucleus's UMI counts binomially into two halves, build an axis on
                each, correlate. This asks how much of a nucleus's score survives resampling
                its own transcripts, and isolates capture noise from gene disagreement. It is
                the stricter test and the one that bears on sparsity.

CONTROLS

  upper bound   the same axis under two different control-gene draws. Control-gene choice is
                a nuisance the paper already averages over, so this correlation is the
                ceiling any split can reach, and a split near it means the split is not the
                limiting factor
  shuffled      counts permuted across nuclei within each gene, destroying per-nucleus
                structure while preserving every gene's marginal. The split-half correlation
                there is what this procedure returns on noise

WHAT THE NUMBER MEANS FOR THE PAPER. An observed correlation between the axis and anything
else is attenuated by the square root of the axis's reliability. If reliability is 0.5, a
true correlation of 0.2 reads as 0.14, and section 2.5's ceiling of about 2% of variance
corresponds to something larger in truth. The script prints that correction so the bound can
be stated in disattenuated terms rather than left implicit.

This does NOT rescue a low value. Low reliability means the absences are uninformative, and
the honest response is to say so rather than to disattenuate a null into a result.

Output: results/axis_reliability.json
"""
import json
import os
import warnings
from collections import defaultdict

import numpy as np
from scipy import stats

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "rel"
N_SPLITS = 20


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
    cols = sorted(c for c in colsum
                  if c in col2bc
                  and (meta[col2bc[c]].get("VMH_cluster") or "").strip().startswith("Lepr"))
    cidx = {c: i for i, c in enumerate(cols)}
    depth = np.array([float(colsum[c]) for c in cols])
    print(f">> {len(cols)} Lepr-cluster nuclei", flush=True)

    panels = {p: [g for g in cfg.PANELS[p] if g in name2row] for p in ("HYPER", "DEPOL")}
    allrows = {name2row[g] for gs in panels.values() for g in gs}
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    pool = [int(r) for r in rowsum if rowsum[r] > 0]
    ctrl_pool = list(rng.choice(pool, min(3000, len(pool)), replace=False))
    ent = V.pass_rows(PREFIX, allrows | set(ctrl_pool), TAG, "rel")

    counts = defaultdict(lambda: np.zeros(len(cols)))
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset:
            counts[int(r)][cidx[c]] = float(v)

    detect = {r: float((counts[r] > 0).mean()) for r in allrows}

    def z(vec, dep):
        v = np.log1p(vec / np.where(dep <= 0, 1, dep) * cfg.NORM_TARGET)
        sd = v.std()
        return (v - v.mean()) / (sd if sd > 0 else 1.0)

    def axis_from(gene_rows, cnt, dep, rng_):
        """Module-score axis from a given gene set and count matrix."""
        out = {}
        for p in ("HYPER", "DEPOL"):
            rows_p = [r for r in gene_rows[p]]
            if not rows_p:
                return None
            pan = np.mean([z(cnt[r], dep) for r in rows_p], axis=0)
            ctl_rows = []
            for r in rows_p:
                d = detect.get(r, 0.0)
                near = sorted(ctrl_pool,
                              key=lambda q: abs((cnt[q] > 0).mean() - d))[:40]
                ctl_rows += list(rng_.choice(near, min(10, len(near)), replace=False))
            ctl = np.mean([z(cnt[r], dep) for r in ctl_rows], axis=0)
            out[p] = pan - ctl
        return out["HYPER"] - out["DEPOL"]

    full_rows = {p: [name2row[g] for g in panels[p]] for p in panels}
    out = {"n_nuclei": len(cols), "n_splits": N_SPLITS,
           "panel_sizes": {p: len(panels[p]) for p in panels},
           "why": ("2.3, 2.4 and 2.5 all report absences on this axis; none of them is "
                   "interpretable without knowing how much signal a single nucleus's "
                   "value carries")}

    # ---------------------------------------------------------------- upper bound
    a1 = axis_from(full_rows, counts, depth, np.random.default_rng(1))
    a2 = axis_from(full_rows, counts, depth, np.random.default_rng(2))
    ceil = float(stats.spearmanr(a1, a2).statistic)
    out["control_gene_draw_ceiling"] = round(ceil, 4)
    print(f">> ceiling (two control-gene draws, same genes): rho {ceil:.4f}", flush=True)

    # ---------------------------------------------------------------- panel split
    print(f">> panel split, {N_SPLITS} alternating partitions", flush=True)
    pr = []
    for s in range(N_SPLITS):
        rg = np.random.default_rng(100 + s)
        ha, hb = {}, {}
        ok = True
        for p in panels:
            rows_p = sorted((name2row[g] for g in panels[p]),
                            key=lambda r: detect.get(r, 0.0))
            rg.shuffle(rows_p)
            half = len(rows_p) // 2
            ha[p], hb[p] = rows_p[:half], rows_p[half:2 * half]
            if half < 2:
                ok = False
        if not ok:
            continue
        xa = axis_from(ha, counts, depth, rg)
        xb = axis_from(hb, counts, depth, rg)
        if xa is not None and xb is not None:
            pr.append(float(stats.spearmanr(xa, xb).statistic))
    r_half = float(np.mean(pr)) if pr else float("nan")
    r_full = 2 * r_half / (1 + r_half) if r_half == r_half else float("nan")
    out["panel_split"] = {"half_rho_mean": round(r_half, 4),
                          "half_rho_min": round(float(np.min(pr)), 4) if pr else None,
                          "half_rho_max": round(float(np.max(pr)), 4) if pr else None,
                          "spearman_brown_full_panel": round(r_full, 4),
                          "n_splits_scored": len(pr)}
    print(f"   half-panel rho {r_half:.4f}, Spearman-Brown full-panel "
          f"{r_full:.4f}", flush=True)

    # ---------------------------------------------------------------- count split
    print(f">> count split, {N_SPLITS} binomial partitions", flush=True)
    cr = []
    for s in range(N_SPLITS):
        rg = np.random.default_rng(200 + s)
        ca, cb = {}, {}
        for r in allrows | set(ctrl_pool):
            v = counts[r]
            first = rg.binomial(v.astype(int), 0.5).astype(float)
            ca[r], cb[r] = first, v - first
        # Each half gets ITS OWN realised depth. Using depth/2 for both assumed the
        # thinning split every nucleus exactly evenly, which binomial thinning does
        # not do, and wrote a deterministic quantity into a stochastic one.
        # Only the panel and control genes are thinned here, so the unthinned
        # remainder of the transcriptome is halved alongside them to keep the
        # normalisation on the same scale.
        thinned = set(ca)
        rest = depth - np.array([sum(counts[r][i] for r in thinned)
                                 for i in range(len(cols))])
        da = (np.array([sum(ca[r][i] for r in ca) for i in range(len(cols))])
              + rest / 2.0)
        db_ = (np.array([sum(cb[r][i] for r in cb) for i in range(len(cols))])
               + rest / 2.0)
        xa = axis_from(full_rows, ca, da, rg)
        xb = axis_from(full_rows, cb, db_, rg)
        if xa is not None and xb is not None:
            cr.append(float(stats.spearmanr(xa, xb).statistic))
    c_half = float(np.mean(cr)) if cr else float("nan")
    c_full = 2 * c_half / (1 + c_half) if c_half == c_half else float("nan")
    out["count_split"] = {"half_rho_mean": round(c_half, 4),
                          "half_rho_min": round(float(np.min(cr)), 4) if cr else None,
                          "half_rho_max": round(float(np.max(cr)), 4) if cr else None,
                          "spearman_brown_full_depth": round(c_full, 4),
                          "n_splits_scored": len(cr)}
    print(f"   half-count rho {c_half:.4f}, Spearman-Brown full-depth "
          f"{c_full:.4f}", flush=True)

    # ---------------------------------------------------------------- noise floor
    # TWO INDEPENDENT permutations. The first version permuted once and scored that
    # one matrix twice with different control-gene draws, which is exactly what the
    # ceiling measures, and it returned 0.9216 against a ceiling of 0.9288. That was
    # the correct answer to the wrong question and it invalidated the run.
    print(">> shuffled control, independent permutations", flush=True)
    sh_a = {r: np.random.default_rng(900 + i).permutation(counts[r])
            for i, r in enumerate(sorted(allrows | set(ctrl_pool)))}
    sh_b = {r: np.random.default_rng(1900 + i).permutation(counts[r])
            for i, r in enumerate(sorted(allrows | set(ctrl_pool)))}
    xa = axis_from(full_rows, sh_a, depth, np.random.default_rng(3))
    xb = axis_from(full_rows, sh_b, depth, np.random.default_rng(4))
    floor = float(stats.spearmanr(xa, xb).statistic)
    out["shuffled_floor"] = round(floor, 4)
    print(f"   rho on shuffled counts: {floor:.4f}", flush=True)

    # ---------------------------------------------------------------- what it means
    rel = c_full if c_full == c_full else r_full
    out["reliability_used"] = round(rel, 4)
    if rel > 0:
        out["attenuation_factor"] = round(float(np.sqrt(rel)), 4)
        out["disattenuation_note"] = (
            "an observed correlation between the axis and any other variable is "
            "attenuated by sqrt(reliability); divide by this factor for the "
            "disattenuated value")
    out["floor_ok"] = bool(abs(floor) < 0.15)
    # ---------------------------------------------------------------- derived keys
    # Every number the reliability paragraphs state is emitted here, COMPUTED FROM the
    # emitted inputs rather than written in. A key holding 0.511 is checkable only
    # against itself; a key holding sqrt of the emitted 0.261 makes the manifest check
    # the chain. The prose arithmetic in 2.2 and 2.5 is otherwise exactly where an
    # inferred number can sit and look like a read one.
    import json as _j
    import vmhlib as _V

    # between-cluster and within-cluster variance, across the 24 VMH clusters
    # THREE THINGS WERE WRONG HERE until 3 Nov 2026, all inflating between_cluster:
    #   (1) cl_means came from cluster_axis_corrected.json, over all 24 VMH clusters
    #       and 6,048 nuclei, while tot_var came from a1, over the 2,913 Lepr-cluster
    #       nuclei this script runs on. The ratio divided one population's
    #       between-cluster variance by another's total. The 24-cluster value is
    #       0.01076 against 0.00465 for the nine Lepr clusters.
    #   (2) the two also came from DIFFERENT AXES: a1 is built here, the stored means
    #       from script 41's own control-gene draw.
    #   (3) n_min/n_med/n_max were taken over the 24 clusters, so cluster_rel was
    #       reported at averaging sizes that do not occur in this population.
    # Cluster means are now computed from a1 itself, over this script's own nuclei.
    clu = np.array([(meta[col2bc[c]].get("VMH_cluster") or "").strip() for c in cols])
    a1a = np.asarray(a1, float)
    assert len(clu) == len(a1a), "cluster labels and axis disagree in length"
    cl_keys = sorted(set(clu))
    cl_means = np.array([a1a[clu == k].mean() for k in cl_keys])
    cl_ns = np.array([(clu == k).sum() for k in cl_keys], float)
    tot_var = float(a1a.var(ddof=1))
    between_raw = float(np.average((cl_means - np.average(cl_means, weights=cl_ns)) ** 2,
                                   weights=cl_ns))

    # The raw weighted variance of cluster means carries the sampling noise of those
    # means: each is an average of n_c noisy scores. With W = tot_var - true,
    # E[raw] = true + k*W/N, so true = (raw - k*tot_var/N) / (1 - k/N). Unlike the
    # raw statistic this CAN return a negative value under no structure, which is
    # the property that makes a zero control meaningful.
    k_cl, n_tot = len(cl_keys), float(len(a1a))
    between = (between_raw - k_cl * tot_var / n_tot) / (1.0 - k_cl / n_tot)
    noise_v = (1.0 - rel) * tot_var
    within_v = tot_var - noise_v - between

    def _between_control(reps=60):
        """Inject a KNOWN between-cluster variance and check recovery, zero included."""
        rngc = np.random.default_rng(99)
        W = max(tot_var - max(between, 0.0), 1e-9)
        res = []
        for true_b in (0.0, 0.002, 0.005):
            got, realized = [], []
            for _ in range(reps):
                eff = {kk: rngc.normal(0, np.sqrt(true_b)) for kk in cl_keys}
                ev = np.array([eff[kk] for kk in cl_keys])
                realized.append(float(np.average(
                    (ev - np.average(ev, weights=cl_ns)) ** 2, weights=cl_ns)))
                y = (np.array([eff[x] for x in clu])
                     + rngc.normal(0, np.sqrt(W), len(clu)))
                m = np.array([y[clu == kk].mean() for kk in cl_keys])
                raw = float(np.average((m - np.average(m, weights=cl_ns)) ** 2,
                                       weights=cl_ns))
                tv = float(y.var(ddof=1))
                got.append((raw - k_cl * tv / n_tot) / (1.0 - k_cl / n_tot))
            res.append({"injected": true_b,
                        # The nominal parameter is not what nine clusters actually
                        # realize. Recovery must be judged against the variance the
                        # injected effects really have, not against true_b.
                        "realized_mean": round(float(np.mean(realized)), 6),
                        "recovered_mean": round(float(np.mean(got)), 6),
                        "recovered_se": round(
                            float(np.std(got, ddof=1) / np.sqrt(len(got))), 6)})
        return res

    between_ctrl = _between_control()
    print(f">> between-cluster control: {between_ctrl}", flush=True)

    def cluster_rel(n):
        return between / (between + (within_v + noise_v) / n)

    n_min, n_med, n_max = int(cl_ns.min()), int(np.median(cl_ns)), int(cl_ns.max())
    class_n = 1857          # conserved Lepr class, the largest averaging unit quoted

    out["derived"] = {
        "attenuation_ceiling": round(float(np.sqrt(rel)), 4),
        "attenuation_ceiling_from": ("sqrt of count_split.spearman_brown_full_depth; "
                                     "the largest correlation this axis can show with "
                                     "a perfectly measured variable"),
        "variance_split_pct": {
            "capture_noise": round(100 * noise_v / tot_var, 1),
            "between_cluster": round(100 * between / tot_var, 1),
            "between_cluster_raw_pct": round(100 * between_raw / tot_var, 1),
            "between_raw": round(between_raw, 6),
            "between_noise_corrected": round(between, 6),
            "between_control": between_ctrl,
            "n_clusters_used": k_cl,
            "within_cluster": round(100 * within_v / tot_var, 1),
            "note": ("noise is (1 - reliability) of the total; between-cluster is the "
                     "size-weighted variance of the nine Lepr cluster means, "
                     "corrected for the sampling noise in those means; within-cluster "
                     "is the remainder")},
        "cluster_mean_reliability": {
            "smallest_cluster": {"n": n_min, "reliability": round(cluster_rel(n_min), 4)},
            "median_cluster": {"n": n_med, "reliability": round(cluster_rel(n_med), 4)},
            "largest_cluster": {"n": n_max, "reliability": round(cluster_rel(n_max), 4)},
            "at_200": {"n": 200, "reliability": round(cluster_rel(200), 4)},
            "conserved_class": {"n": class_n,
                                "reliability": round(cluster_rel(class_n), 4)},
            "note": ("estimates the cluster's value RELATIVE TO OTHER CLUSTERS, so "
                     "within-cluster variation counts against it. Sizes are the nine "
                     "Lepr clusters' own, computed here rather than read from "
                     "cluster_axis_corrected.json, so the floor follows the data rather "
                     "than a literal in the text")},
        "detectability_bound_conversion": {
            "observed_sd_units": [0.10, 0.20],
            "underlying_sd_units": [round(0.10 / np.sqrt(rel), 3),
                                    round(0.20 / np.sqrt(rel), 3)],
            "note": "injected level divided by sqrt(reliability)"},
        "controls": {
            "count_split_floor": out["shuffled_floor"],
            "control_gene_draw_ceiling": out["control_gene_draw_ceiling"],
            "note": ("the floor is two independently permuted matrices and must be near "
                     "zero; the ceiling is the same genes and counts under two "
                     "control-gene draws and is the most any split can reach")}}

    d = out["derived"]
    print(f"\n  derived quantities, for the manifest:")
    print(f"    attenuation ceiling            {d['attenuation_ceiling']}")
    print(f"    variance split (noise/btw/wtn) "
          f"{d['variance_split_pct']['capture_noise']}% / "
          f"{d['variance_split_pct']['between_cluster']}% / "
          f"{d['variance_split_pct']['within_cluster']}%")
    for k in ("smallest_cluster", "at_200", "conserved_class"):
        e = d["cluster_mean_reliability"][k]
        print(f"    reliability, {k:<18} n={e['n']:<6} {e['reliability']}")
    print(f"    bound 0.10 observed ->         "
          f"{d['detectability_bound_conversion']['underlying_sd_units'][0]} underlying")
    print(f"    bound 0.20 observed ->         "
          f"{d['detectability_bound_conversion']['underlying_sd_units'][1]} underlying")

    print(f"\n{'=' * 66}")
    if not out["floor_ok"]:
        print(f"  CONTROL FAILED: shuffled counts give rho {floor:.4f}, which should be")
        print("  near zero. The procedure is finding agreement that is not in the data.")
        print("  Nothing below is readable.")
        json.dump(out, open(os.path.join(cfg.RESULTS, "axis_reliability.json"), "w"),
                  indent=2)
        return
    print(f"  reliability (count split, Spearman-Brown): {rel:.4f}")
    print(f"  observed correlations are attenuated by sqrt(r) = "
          f"{np.sqrt(max(rel, 0)):.4f}")
    if rel < 0.3:
        print("\n  LOW. A single nucleus's axis value is mostly noise, and the absences")
        print("  in 2.3, 2.4 and 2.5 are then guaranteed by the measurement rather than")
        print("  informative about biology. Say so; do not disattenuate a null into a")
        print("  result.")
    elif rel < 0.6:
        print("\n  MODERATE. The absences stand but their bounds should be stated in")
        print("  disattenuated terms, since the measured ceiling understates the true one.")
    else:
        print("\n  The axis carries substantial per-nucleus signal and the absences in")
        print("  2.3, 2.4 and 2.5 are about the quantity they claim to be about.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "axis_reliability.json"), "w"),
              indent=2)
    print("\nWritten to results/axis_reliability.json")


if __name__ == "__main__":
    main()
