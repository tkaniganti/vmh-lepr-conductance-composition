#!/usr/bin/env python3
"""
45_patchseq_within_family.py - does a transcript-derived conductance axis find structure
that measured electrophysiology says is there, inside a single family?

WHY THE FIRST ATTEMPT FAILED. 42_patchseq_calibration.py pooled Pvalb and Sst and asked
whether the axis split them. It did, but so did 99 of 100 arbitrary ion-channel panels:
two families that differ in hundreds of channel transcripts separate on almost anything.
The test was too easy to be informative, and the result said nothing about conductance
composition.

WHAT CHANGES HERE. Two things.

Everything happens INSIDE one family, so the between-family difference that swamped the
first attempt is gone.

And the ground truth comes from the measured electrophysiology rather than from a claim in
the paper. The axis is still built from transcripts alone, exactly as in the manuscript;
the ephys is used only to decide whether a real split exists in those same cells. That is
the comparison the manuscript cannot make for VMH, because no hypothalamic dataset has both.

THE TEST. Within the family:

  TRUTH:     cluster the cells on their measured ephys features. Does a two-component
             solution beat one? This says whether the family really is electrically split.
  CANDIDATE: build the conductance axis from transcripts only and run the manuscript's own
             tests on it, dip and Gaussian mixture.

  Four outcomes:
    ephys splits, axis splits      the method detects real structure. 2.4's null gains
                                   weight, because the test works where truth exists
    ephys splits, axis flat        the method misses real structure. 2.4's bound is
                                   optimistic and the stated two-SD figure overstates it
    ephys flat, axis flat          consistent, and consistent with Scala's report of
                                   continuous within-family variation. Weak but not nothing
    ephys flat, axis splits        the axis manufactures structure. Worst case

CONTROL, AGAIN. Arbitrary ion-channel panels are run through the same pipeline. If they
split the family as often as the chosen panel does, the panel is not what is doing the work
and no outcome above can be read. This is the control the first attempt failed, and it is
the first thing to look at in the output.

PANEL. Kv3 (Kcnc1-4) against HCN, chosen from general knowledge of cortical interneuron
physiology, not from anything measured here. Stated plainly because the control is what
makes it interpretable.

Data: berenslab/mini-atlas. Output: results/patchseq_within_family.json
"""
import argparse
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

KV3 = ["Kcnc1", "Kcnc2", "Kcnc3", "Kcnc4"]
HCN = ["Hcn1", "Hcn2", "Hcn3", "Hcn4"]
N_RANDOM = 100


def fail(msg):
    print(msg)
    sys.exit(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--ttypes", required=True)
    ap.add_argument("--ephys", required=True)
    ap.add_argument("--family", default="Pvalb")
    ap.add_argument("--cell-col", default="Cell",
                    help="cell identifier in the ttypes table")
    ap.add_argument("--ephys-cell-col", default="cell id",
                    help="cell identifier in the ephys table, which uses a different "
                         "column name from the metadata table")
    ap.add_argument("--family-col", default="RNA family")
    ap.add_argument("--ion-channels", default=None)
    a = ap.parse_args()

    import pandas as pd
    from sklearn.mixture import GaussianMixture
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    try:
        import diptest
    except ImportError:
        fail("diptest not installed")

    cnt = pd.read_csv(a.counts, sep=None, engine="python", index_col=0)
    tt = pd.read_csv(a.ttypes, sep=None, engine="python")
    eph = pd.read_csv(a.ephys, sep=None, engine="python")
    if a.cell_col not in tt.columns:
        fail(f"--cell-col '{a.cell_col}' not in ttypes: {list(tt.columns)[:12]}")
    if a.ephys_cell_col not in eph.columns:
        fail(f"--ephys-cell-col '{a.ephys_cell_col}' not in ephys: "
             f"{list(eph.columns)[:12]}")
    tt = tt.set_index(a.cell_col)
    eph = eph.set_index(a.ephys_cell_col)

    if KV3[0] not in cnt.columns:
        if KV3[0] in cnt.index:
            cnt = cnt.T
        else:
            fail(f"{KV3[0]} in neither axis of the counts table")

    fam = tt[a.family_col].astype(str).str.strip()
    cells = cnt.index.intersection(tt.index).intersection(eph.index)
    cells = [c for c in cells if fam.get(c) == a.family]
    if len(cells) < 60:
        fail(f"only {len(cells)} {a.family} cells with counts and ephys; too few")
    print(f">> {len(cells)} {a.family} cells with both modalities", flush=True)

    # ---------------------------------------------------------------- TRUTH: ephys
    E = eph.loc[cells].select_dtypes(include=[np.number])
    E = E.loc[:, E.notna().mean() > 0.9].fillna(E.median())
    Z = StandardScaler().fit_transform(E.to_numpy(float))
    npc = min(5, Z.shape[1])
    P = PCA(n_components=npc, random_state=0).fit_transform(Z)
    bic_e = {k: GaussianMixture(k, random_state=0, n_init=10).fit(P).bic(P)
             for k in (1, 2, 3)}
    best_e = int(min(bic_e, key=bic_e.get))
    d_e, p_e = diptest.diptest(np.ascontiguousarray(P[:, 0]))
    print(f">> ephys: {E.shape[1]} features, {npc} PCs, best k = {best_e}, "
          f"dip p on PC1 = {p_e:.4f}", flush=True)

    # POSITIVE CONTROL ON THE GROUND TRUTH ARM. Everything downstream is read against
    # this call, and BIC on PCs of 29 correlated features is where BIC adds components
    # to fit shape. Shuffle each feature independently, destroying joint structure and
    # preserving every marginal, then recluster.
    rng_e = np.random.default_rng(0)
    null_k = []
    for _ in range(100):
        Zs = np.column_stack([rng_e.permutation(Z[:, j]) for j in range(Z.shape[1])])
        Ps = PCA(n_components=npc, random_state=0).fit_transform(Zs)
        b = {k: GaussianMixture(k, random_state=0, n_init=5).fit(Ps).bic(Ps)
             for k in (1, 2, 3)}
        null_k.append(int(min(b, key=b.get)))
    ephys_null_rate = sum(1 for k in null_k if k > 1) / len(null_k)
    print(f">> ephys control: BIC picks k > 1 on {ephys_null_rate:.0%} of shuffled "
          f"feature sets", flush=True)
    if ephys_null_rate > 0.2:
        print("   GROUND TRUTH ARM UNRELIABLE. Read the dip test on PC1 instead and")
        print("   treat any outcome resting on 'ephys says split' as unsupported.")

    # ---------------------------------------------------------------- the axis
    X = cnt.loc[cells].to_numpy(float)
    dep = X.sum(1, keepdims=True); dep[dep == 0] = 1
    L = np.log1p(X / dep * 1e4)
    Zg = (L - L.mean(0)) / np.where(L.std(0) > 0, L.std(0), 1)
    genes = list(cnt.columns); gi = {g: i for i, g in enumerate(genes)}
    detect = (X > 0).mean(0)
    kv3 = [g for g in KV3 if g in gi]
    hcn = [g for g in HCN if g in gi]
    if not kv3 or not hcn:
        fail("a panel is empty after gene matching")

    def ctrls(panel, rng, n_per=25):
        out = []
        for g in panel:
            order = np.argsort(np.abs(detect - detect[gi[g]]))
            pool = [j for j in order[:400] if genes[j] not in panel]
            out += list(rng.choice(pool, min(n_per, len(pool)), replace=False))
        return out

    def axis(pa, pb, rng):
        return (Zg[:, [gi[g] for g in pa]].mean(1) - Zg[:, ctrls(pa, rng)].mean(1)) - \
               (Zg[:, [gi[g] for g in pb]].mean(1) - Zg[:, ctrls(pb, rng)].mean(1))

    def test1d(v):
        v = np.asarray(v, float)
        r = v.reshape(-1, 1)
        b = {k: GaussianMixture(k, random_state=0, n_init=10).fit(r).bic(r)
             for k in (1, 2, 3)}
        D, p = diptest.diptest(np.ascontiguousarray(v))
        return int(min(b, key=b.get)), float(D), float(p)

    rng = np.random.default_rng(0)
    k_ax, D_ax, p_ax = test1d(axis(kv3, hcn, rng))
    print(f">> axis:  best k = {k_ax}, dip D = {D_ax:.4f}, p = {p_ax:.4f}", flush=True)

    # ---------------------------------------------------------------- control
    pool = genes
    desc = "all genes"
    if a.ion_channels:
        ic = pd.read_csv(a.ion_channels)
        col = next((c for c in ic.columns if "symbol" in c.lower()), None)
        if col is None:
            fail(f"no symbol column in {list(ic.columns)}")
        want = {str(x).strip().title() for x in ic[col]}
        pool = [g for g in genes if g in want]
        desc = f"{len(pool)} ion-channel genes"
        if len(pool) < 40:
            fail("too few ion-channel genes matched")
    print(f">> control pool: {desc}", flush=True)

    hits = 0
    for _ in range(N_RANDOM):
        pick = list(rng.choice(pool, len(kv3) + len(hcn), replace=False))
        try:
            k, _, _ = test1d(axis(pick[:len(kv3)], pick[len(kv3):], rng))
            hits += int(k > 1)
        except Exception:
            pass
    rate = hits / N_RANDOM

    out = {"family": a.family, "n_cells": len(cells),
           "ephys_n_features": int(E.shape[1]), "ephys_best_k": best_e,
           "ephys_bic": {str(k): round(float(v), 2) for k, v in bic_e.items()},
           "ephys_dip_p_pc1": round(float(p_e), 5),
           "ephys_null_split_rate": ephys_null_rate,
           "ephys_truth_reliable": bool(ephys_null_rate <= 0.2),
           "axis_best_k": k_ax, "axis_dip_D": round(D_ax, 5),
           "axis_dip_p": round(p_ax, 5),
           "random_panel_split_rate": rate, "control_pool": desc,
           "panel_note": "Kv3 against HCN, chosen from general knowledge, not measured here"}

    print(f"\n{'=' * 64}\n  control: arbitrary panels splitting this family "
          f"{hits}/{N_RANDOM}")
    if rate > 0.2:
        out["readable"] = False
        print("  TOO PERMISSIVE. Arbitrary panels split this family as readily as the")
        print("  chosen one, so nothing below can be attributed to the panel. This is")
        print("  the same failure as the pooled design and the result is not usable.")
    else:
        out["readable"] = True
        ephys_splits = best_e > 1
        axis_splits = k_ax > 1
        out["ephys_splits"] = bool(ephys_splits)
        out["axis_splits"] = bool(axis_splits)
        print(f"  ephys says split: {ephys_splits}")
        print(f"  axis says split:  {axis_splits}")
        if ephys_splits and axis_splits:
            print("\n  The method detects structure an independent modality shows exists.")
            print("  2.4's null gains weight: the test works where truth is present.")
        elif ephys_splits and not axis_splits:
            print("\n  The method MISSES structure the ephys shows is there. 2.4's bound")
            print("  is optimistic and the two-SD figure overstates what it can detect.")
            print("  This is the outcome that would change what the manuscript claims.")
        elif not ephys_splits and not axis_splits:
            print("\n  Neither splits. Consistent with Scala's continuous within-family")
            print("  variation, but a null against a null bounds nothing on its own.")
        else:
            print("\n  The axis splits where the ephys does not. The axis may be")
            print("  manufacturing structure; nothing here supports using it.")

    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/patchseq_within_family.json", "w"), indent=2)
    print("\nWritten to results/patchseq_within_family.json")


if __name__ == "__main__":
    main()
