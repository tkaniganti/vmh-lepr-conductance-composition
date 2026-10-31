#!/usr/bin/env python3
"""
42_patchseq_calibration.py - calibrate the discreteness test against a dataset where the
answer is known from measured phenotype.

WHY. Section 2.4 reports that unsupervised testing on the effector-balance axis selects a
single component and does not reject unimodality, with a detectability bound from simulation.
The only calibration in the paper is arcuate POMC, which is a worked negative case: the test
fails to resolve subdivision known to exist, at a sample size where its bound is weak. It
bounds the test and does not show the test working.

Scala et al. (2021, Nature 598:144-150) recorded electrophysiology, morphology and
transcriptome from the same 1,300 mouse motor cortex neurons, and report that transcriptomic
families have largely distinct morpho-electric phenotypes while variation within most
families is continuous. That is a known answer in the same cells, from an independent
modality. Running this test on their transcriptomes alone, with their electrophysiology set
aside, asks whether a transcript-derived conductance axis reproduces it.

TWO TESTS, OPPOSITE EXPECTED ANSWERS.

  POSITIVE: Pvalb and Sst pooled. These families differ unambiguously in measured
  electrophysiology, fast-spiking against adapting. A test that works should find two
  components. If it does not, the null in 2.4 carries much less weight than stated.

  NEGATIVE: within Pvalb alone. Scala et al. report continuous variation within families, so
  a single component here agrees with the measured answer.

Passing both would let 2.4 say the method detects discreteness where an independent modality
shows it and reports continuity where that modality reports continuity. Neither of those can
be said now.

THE PANEL IS A JUDGMENT AND IS TREATED AS ONE. The hyperpolarizing and depolarizing panels
in the manuscript were chosen for the leptin response and are not the axis that separates
cortical interneuron families. This uses Kv3 channels (Kcnc1-4), the canonical determinant
of fast spiking in Pvalb interneurons, against HCN. I am reasoning from general knowledge
of cortical interneuron physiology in choosing that, not reporting a value I traced, so the
random-panel control below matters more than usual: if arbitrary detectability-matched panels
separate these families just as well, the test is not measuring conductance composition and
the positive result would mean nothing.

DATA. Counts and t-type labels from https://github.com/berenslab/mini-atlas. Download before
running; the script reports what it cannot find rather than guessing at file names.

Output: results/patchseq_calibration.json
"""
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

KV3 = ["Kcnc1", "Kcnc2", "Kcnc3", "Kcnc4"]
HCN = ["Hcn1", "Hcn2", "Hcn3", "Hcn4"]
N_RANDOM_PANELS = 100
DATA_HINT = """
Expected a counts matrix and a per-cell t-type table from berenslab/mini-atlas, placed in
data/patchseq/. The repository layout is not assumed here. Point the two paths below at
whatever the download actually produced:

    --counts  <genes x cells or cells x genes table, gene symbols somewhere>
    --ttypes  <per-cell table with a t-type or cluster column>

Both are read with pandas, so csv, tsv and gzipped variants all work.
"""


def fail(msg):
    print(msg)
    print(DATA_HINT)
    sys.exit(1)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--ttypes", required=True)
    ap.add_argument("--ttype-col", default="RNA family",
                    help="column holding the family label")
    ap.add_argument("--cell-col", default="Cell",
                    help="column in the ttypes table holding the cell identifier that "
                         "matches the counts table; the meta table is keyed on Number, "
                         "which does not match, so this has to be set explicitly")
    ap.add_argument("--ion-channels", default=None,
                    help="ion-channel-genes-group-177.csv from mini-atlas. Random control "
                         "panels are drawn from this list instead of from all genes, which "
                         "asks whether Kv3 against HCN separates the families or whether "
                         "any channel panel does. The list is human HGNC symbols and is "
                         "matched to mouse by title case, which is approximate.")
    args = ap.parse_args()

    import pandas as pd
    from scipy import stats
    from sklearn.mixture import GaussianMixture
    try:
        import diptest
    except ImportError:
        fail("diptest not installed in this environment")

    for p in (args.counts, args.ttypes):
        if not os.path.exists(p):
            fail(f"not found: {p}")

    cnt = pd.read_csv(args.counts, sep=None, engine="python", index_col=0)
    tt = pd.read_csv(args.ttypes, sep=None, engine="python")
    if args.cell_col not in tt.columns:
        fail(f"--cell-col '{args.cell_col}' not in {list(tt.columns)[:15]}")
    tt = tt.set_index(args.cell_col)
    print(f">> counts {cnt.shape}, ttypes {tt.shape}", flush=True)

    # orient so rows are cells
    if not set(KV3[:1]).issubset(cnt.columns):
        if set(KV3[:1]).issubset(cnt.index):
            cnt = cnt.T
            print(">> transposed counts so rows are cells", flush=True)
        else:
            fail(f"gene {KV3[0]} not found in either axis of the counts table")

    col = args.ttype_col
    if col not in tt.columns:
        fail(f"--ttype-col '{col}' not in {list(tt.columns)[:15]}")
    print(f">> using family column '{col}'", flush=True)

    shared = cnt.index.intersection(tt.index)
    if len(shared) < 100:
        fail(f"only {len(shared)} cells shared between the two tables; check the indices")
    cnt, tt = cnt.loc[shared], tt.loc[shared]
    fam = tt[col].astype(str).str.strip()
    print(f">> {len(shared)} cells; families: "
          f"{dict(fam.value_counts().head(10))}", flush=True)

    present = {g: g for g in KV3 + HCN if g in cnt.columns}
    missing = [g for g in KV3 + HCN if g not in cnt.columns]
    if missing:
        print(f">> not in the matrix, dropped: {missing}", flush=True)
    kv3 = [g for g in KV3 if g in present]
    hcn = [g for g in HCN if g in present]
    if not kv3 or not hcn:
        fail("one of the two panels is empty after matching gene symbols")

    # CP10k, log1p, z-score per gene, then panel minus matched controls
    X = cnt.to_numpy(float)
    depth = X.sum(axis=1, keepdims=True)
    depth[depth == 0] = 1
    L = np.log1p(X / depth * 1e4)
    Z = (L - L.mean(0)) / np.where(L.std(0) > 0, L.std(0), 1)
    genes = list(cnt.columns)
    gi = {g: i for i, g in enumerate(genes)}
    detect = (X > 0).mean(0)

    def matched_controls(panel, rng, n_per=25):
        out = []
        for g in panel:
            d = detect[gi[g]]
            order = np.argsort(np.abs(detect - d))
            pool = [j for j in order[:400] if genes[j] not in panel]
            out += list(rng.choice(pool, size=min(n_per, len(pool)), replace=False))
        return out

    def axis_for(panel_a, panel_b, rng):
        ca, cb = matched_controls(panel_a, rng), matched_controls(panel_b, rng)
        a = Z[:, [gi[g] for g in panel_a]].mean(1) - Z[:, ca].mean(1)
        b = Z[:, [gi[g] for g in panel_b]].mean(1) - Z[:, cb].mean(1)
        return a - b

    def test(v):
        v = np.asarray(v, float).reshape(-1, 1)
        bic = {k: GaussianMixture(k, random_state=0, n_init=5).fit(v).bic(v)
               for k in (1, 2, 3)}
        D, p = diptest.diptest(v.ravel())
        return {"best_k": int(min(bic, key=bic.get)),
                "bic": {str(k): round(float(x), 2) for k, x in bic.items()},
                "dip_D": round(float(D), 5), "dip_p": float(p), "n": int(v.size)}

    # control pool: the paper's own ion-channel list if given, else every gene
    pool = genes
    pool_desc = "all genes"
    if args.ion_channels:
        ic = pd.read_csv(args.ion_channels)
        symcol = next((c for c in ic.columns if "symbol" in c.lower()), None)
        if symcol is None:
            fail(f"no symbol column in {list(ic.columns)}")
        mouse = {str(x).strip().title() for x in ic[symcol]}
        pool = [g for g in genes if g in mouse]
        pool_desc = (f"{len(pool)} of {len(mouse)} ion-channel genes present "
                     f"(human symbols title-cased to mouse)")
        print(f">> control pool: {pool_desc}", flush=True)
        if len(pool) < 40:
            fail("too few ion-channel genes matched; check the symbol conversion")

    rng = np.random.default_rng(0)
    out = {"n_cells": int(len(shared)), "kv3_panel": kv3, "hcn_panel": hcn,
           "control_pool": pool_desc,
           "source": "Scala et al. 2021, mini-atlas",
           "panel_choice_note": ("Kv3 against HCN, chosen from general knowledge of "
                                 "cortical interneuron physiology and not from any "
                                 "measurement in this dataset; the random-panel control "
                                 "is what makes the positive result interpretable")}

    m_pos = fam.isin(["Pvalb", "Sst"]).to_numpy()
    m_neg = (fam == "Pvalb").to_numpy()
    print(f"\n>> positive test on {int(m_pos.sum())} Pvalb+Sst cells")
    print(f">> negative test on {int(m_neg.sum())} Pvalb cells")
    if m_pos.sum() < 50 or m_neg.sum() < 30:
        fail("too few cells in those families; check how the t-type column is formatted")

    ax = axis_for(kv3, hcn, rng)
    out["positive_pvalb_sst"] = test(ax[m_pos])
    out["negative_within_pvalb"] = test(ax[m_neg])

    print(f"\n   {'test':<22}{'best k':>8}{'dip p':>10}{'n':>7}")
    for k, lab in (("positive_pvalb_sst", "Pvalb + Sst"),
                   ("negative_within_pvalb", "within Pvalb")):
        r = out[k]
        print(f"   {lab:<22}{r['best_k']:>8}{r['dip_p']:>10.4f}{r['n']:>7}")

    # how often do arbitrary panels separate these families?
    print(f"\n>> random-panel control, {N_RANDOM_PANELS} panels", flush=True)
    hits = 0
    for _ in range(N_RANDOM_PANELS):
        pick = list(rng.choice(pool, len(kv3) + len(hcn), replace=False))
        pa, pb = pick[:len(kv3)], pick[len(kv3):]
        try:
            r = test(axis_for(pa, pb, rng)[m_pos])
            hits += int(r["best_k"] > 1)
        except Exception:
            pass
    out["random_panel_two_component_rate"] = hits / N_RANDOM_PANELS
    print(f"   arbitrary panels finding two components on Pvalb+Sst: "
          f"{hits}/{N_RANDOM_PANELS}  (pool: {pool_desc})")

    pos_ok = out["positive_pvalb_sst"]["best_k"] > 1
    neg_ok = out["negative_within_pvalb"]["best_k"] == 1
    rare = out["random_panel_two_component_rate"] < 0.2
    out["calibration_passes"] = bool(pos_ok and neg_ok and rare)
    print(f"\n   separates families:            {pos_ok}")
    print(f"   continuous within Pvalb:       {neg_ok}")
    print(f"   arbitrary panels rarely do:    {rare}")
    if out["calibration_passes"]:
        print("\n   The test detects discreteness an independent modality shows exists,")
        print("   and reports continuity where that modality reports continuity. 2.4 can")
        print("   say so, which is more than the POMC case supports.")
    else:
        print("\n   Calibration does not pass. Read which of the three failed before")
        print("   drawing any conclusion: a failed positive means the test is weaker")
        print("   than 2.4 implies, a failed random control means the panel is not what")
        print("   is doing the work.")

    os.makedirs("results", exist_ok=True)
    json.dump(out, open("results/patchseq_calibration.json", "w"), indent=2)
    print("\nWritten to results/patchseq_calibration.json")


if __name__ == "__main__":
    main()
